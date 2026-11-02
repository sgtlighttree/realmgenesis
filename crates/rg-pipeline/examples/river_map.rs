//! Top-down hillshade + macro river network of a box, as a PNG (to judge drainage patterns
//! by eye): `cargo run -p rg-pipeline --release --example river_map -- <seed> <lat> <lon> <half_km> <out.png>`
//! Same `RG_*` env overrides and cache as `probe_realism`.
#[path = "common/mod.rs"]
mod common;
use glam::DVec3;

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let seed: u64 = a[0].parse().unwrap();
    let (lat0, lon0, half) = (a[1].parse::<f64>().unwrap().to_radians(), a[2].parse::<f64>().unwrap().to_radians(), a[3].parse::<f64>().unwrap() * 1000.0);
    let out = &a[4];
    let w = common::world(seed);
    let r = w.params.tectonics.radius_m;
    let hm = &w.tectonics.height_m;
    let hy = w.hydrology.as_ref().expect("erosion on");
    let px = 1000usize;
    let m_per_px = 2.0 * half / px as f64;
    let at = |x: f64, y: f64| {
        let (la, lo) = (lat0 + y / r, lon0 + x / (r * lat0.cos()));
        DVec3::new(la.cos() * lo.cos(), la.cos() * lo.sin(), la.sin())
    };
    let to_px = |d: DVec3| {
        let (la, lo) = (d.z.asin(), d.y.atan2(d.x));
        let x = (lo - lon0) * r * lat0.cos();
        let y = (la - lat0) * r;
        ((x + half) / m_per_px, (half - y) / m_per_px)
    };
    // Hillshade of the macro surface (drained interpolant), or with RG_MAP_FULL=1 the final
    // surface (macro + detail) and no river overlay.
    let full = std::env::var("RG_MAP_FULL").is_ok();
    let sampler = rg_pipeline::terrain_sampler(&w);
    // Full mode with a network: evaluate through 8×8 sub-patches (far cheaper than per-point
    // cached patches).
    let cut = m_per_px.min(200.0);
    let sub = 8usize;
    let patches: Vec<Option<rg_terrain::detail_net::Patch>> = (0..sub * sub)
        .map(|k| {
            let (bx, by) = (k % sub, k / sub);
            let bs = 2.0 * half / sub as f64;
            let (x, y) = (-half + (bx as f64 + 0.5) * bs, half - (by as f64 + 0.5) * bs);
            if full { sampler.network_patch(at(x, y), bs * 0.75 + 4.0 * m_per_px, cut) } else { None }
        })
        .collect();
    let hgt = |d: DVec3| {
        if !full {
            return hm.sample_drained(d);
        }
        let (px_x, px_y) = to_px(d);
        let (bx, by) = (((px_x / px as f64) * sub as f64).clamp(0.0, sub as f64 - 1.0) as usize, ((px_y / px as f64) * sub as f64).clamp(0.0, sub as f64 - 1.0) as usize);
        match &patches[by * sub + bx] {
            Some(p) => p.heights(d)[0],
            None => sampler.height(d, cut).0,
        }
    };
    let mut img = vec![0u8; px * px * 3];
    for j in 0..px {
        for i in 0..px {
            let (x, y) = (-half + (i as f64 + 0.5) * m_per_px, half - (j as f64 + 0.5) * m_per_px);
            let h = hgt(at(x, y));
            let e = if full { m_per_px } else { 2.0 * m_per_px };
            let gx = (hgt(at(x + e, y)) - hgt(at(x - e, y))) / (2.0 * e);
            let gy = (hgt(at(x, y + e)) - hgt(at(x, y - e))) / (2.0 * e);
            let ex = if full { 1.5 } else { 8.0 };
            let n = DVec3::new(-gx * ex, -gy * ex, 1.0).normalize();
            let shade = n.dot(DVec3::new(-0.5, 0.5, 0.7).normalize()).max(0.0);
            let k = (j * px + i) * 3;
            let c = if h <= 0.0 { [20.0, 40.0, 80.0] } else { let t = (h / 4000.0).min(1.0); [90.0 + 120.0 * t, 110.0 + 80.0 * t, 70.0 + 110.0 * t] };
            for ch in 0..3 {
                img[k + ch] = (c[ch] * (0.35 + 0.85 * shade)).min(255.0) as u8;
            }
        }
    }
    // Rivers: every node above the threshold, drawn to its downstream neighbour.
    let ang = hy.drainage_km2.spacing_rad();
    let v = |f: [f32; 3]| DVec3::new(f[0] as f64, f[1] as f64, f[2] as f64);
    let a_min = if full { f64::MAX } else { 300.0 };
    for (k, _, _, _, d) in hy.drainage_km2.iter_samples() {
        let area = hy.drainage_km2.data[k] as f64;
        if area < a_min || d.angle_between(at(0.0, 0.0)) * r > half * 1.5 {
            continue;
        }
        let f = v(hy.flow_dir.data[k]);
        if f.length() < 0.5 {
            continue;
        }
        // Snap the downstream point to the receiver sample.
        let to = (d + f * ang * 1.2).normalize();
        let recv = nearest_sample(&hy.drainage_km2, to);
        let (x0, y0) = to_px(d);
        let (x1, y1) = to_px(recv);
        let width = (0.6 * (area / a_min).log10() + 0.5).min(4.0);
        let steps = ((x1 - x0).hypot(y1 - y0) * 2.0).ceil().max(1.0) as usize;
        for s in 0..=steps {
            let t = s as f64 / steps as f64;
            let (x, y) = (x0 + (x1 - x0) * t, y0 + (y1 - y0) * t);
            let rad = width.ceil() as i64;
            for dy in -rad..=rad {
                for dx in -rad..=rad {
                    if ((dx * dx + dy * dy) as f64).sqrt() > width {
                        continue;
                    }
                    let (xi, yi) = (x as i64 + dx, y as i64 + dy);
                    if xi >= 0 && yi >= 0 && (xi as usize) < px && (yi as usize) < px {
                        let k = (yi as usize * px + xi as usize) * 3;
                        img[k..k + 3].copy_from_slice(&[30, 90, 200]);
                    }
                }
            }
        }
    }
    let file = std::fs::File::create(out).unwrap();
    let mut enc = png::Encoder::new(std::io::BufWriter::new(file), px as u32, px as u32);
    enc.set_color(png::ColorType::Rgb);
    enc.write_header().unwrap().write_image_data(&img).unwrap();
    eprintln!("wrote {out}");
}

/// The cube-field sample nearest a direction (its unit direction).
fn nearest_sample(f: &rg_geo::CubeField<f32>, d: DVec3) -> DVec3 {
    let (face, u, v) = rg_geo::CubeFace::from_dir(d);
    let r = f.res as f64;
    let i = ((u + 1.0) * 0.5 * r).round().clamp(0.0, r) as u32;
    let j = ((v + 1.0) * 0.5 * r).round().clamp(0.0, r) as u32;
    f.sample_dir(face, i, j)
}
