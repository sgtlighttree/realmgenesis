//! `cargo run -p rg-erosion --release --example erode_preview -- <seed> <res> <out_dir> [key=val ...]`
//!
//! Runs tectonics then erosion, prints timings + geomorphic stats, and writes equirect PNGs:
//! `before.png` / `after.png` (hillshaded hypsometric tint, rivers on `after`).
//! Keys: duration, steps, k, m, umax, maxres, uexp, ubase, belt, elev, hyps, kvar, rough, rho8, seed.
use glam::DVec3;
use rg_erosion::{erode, ErosionParams};
use rg_geo::CubeField;
use std::f64::consts::PI;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let seed: u64 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(3);
    let res: u32 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(320);
    let out = args.get(3).cloned().unwrap_or_else(|| ".".into());
    let mut p = ErosionParams::default();
    for kv in args.iter().skip(4) {
        let (k, v) = kv.split_once('=').expect("key=val");
        let v: f64 = v.parse().expect("number");
        match k {
            "duration" => p.duration_myr = v,
            "steps" => p.steps = v as u32,
            "k" => p.k = v,
            "m" => p.m = v,
            "umax" => p.uplift_max_mm_yr = v,
            "maxres" => p.max_res = v as u32,
            "uexp" => p.uplift_exponent = v,
            "ubase" => p.uplift_base_mm_yr = v,
            "belt" => p.belt_width_km = v,
            "elev" => p.uplift_source = if v > 0.0 { rg_erosion::UpliftSource::Elevation } else { rg_erosion::UpliftSource::Tectonic },
            "hyps" => p.hypsometry_scale = v,
            "kvar" => p.k_variability = v,
            "rough" => p.initial_roughness_m = v,
            "rho8" => p.rho8 = v > 0.0,
            "seed" => p.seed = v as u64,
            _ => panic!("unknown key {k}"),
        }
    }
    std::fs::create_dir_all(&out).unwrap();
    let tp = rg_tectonics::TectonicsParams { seed, res, ..Default::default() };
    let t0 = std::time::Instant::now();
    let snap = rg_tectonics::generate(&tp);
    println!("tectonics {:.2} s (res {res})", t0.elapsed().as_secs_f64());
    let t1 = std::time::Instant::now();
    let r = erode(&snap, tp.radius_m, &p);
    println!("erosion   {:.2} s  {:?}", t1.elapsed().as_secs_f64(), p);
    println!("{:#?}", r.stats);
    println!("rivers: {} segments", r.rivers.len());
    // Area-weighted land hypsometry of the output.
    {
        let g = rg_tectonics::grid::Grid::new(r.height_m.res);
        let mut land: Vec<(f32, f64)> = (0..g.n)
            .map(|i| (r.height_m.data[g.node_sample[i] as usize], g.weight[i] as f64))
            .filter(|(h, _)| *h > 0.0)
            .collect();
        land.sort_by(|a, b| a.0.total_cmp(&b.0));
        let tot: f64 = land.iter().map(|x| x.1).sum();
        let mean = land.iter().map(|x| x.0 as f64 * x.1).sum::<f64>() / tot;
        let above = |t: f32| land.iter().filter(|x| x.0 > t).map(|x| x.1).sum::<f64>() / tot * 100.0;
        let (mut c, mut median) = (0.0, 0.0);
        for x in &land {
            c += x.1;
            if c >= 0.5 * tot {
                median = x.0;
                break;
            }
        }
        println!(
            "hypsometry: mean {mean:.0} m, median {median:.0} m, >1km {:.1}%, >2km {:.1}%, >4km {:.1}%, max {:.0} m",
            above(1000.0), above(2000.0), above(4000.0), land.last().map(|x| x.0).unwrap_or(0.0)
        );
    }

    let w = 1600u32;
    save(&format!("{out}/before.png"), w, &snap.height_m, tp.radius_m, None);
    save(&format!("{out}/after.png"), w, &r.height_m, tp.radius_m, Some(&r));
}

fn save(path: &str, w: u32, f: &CubeField<f32>, radius: f64, er: Option<&rg_erosion::ErosionResult>) {
    let h = w / 2;
    let px = |x: u32, y: u32| -> DVec3 {
        let lat = PI / 2.0 - (y as f64 + 0.5) / h as f64 * PI;
        let lon = -PI + (x as f64 + 0.5) / w as f64 * 2.0 * PI;
        DVec3::new(lat.cos() * lon.cos(), lat.cos() * lon.sin(), lat.sin())
    };
    let hs: Vec<f32> = (0..w * h).map(|k| f.sample(px(k % w, k / w))).collect();
    let mut img = vec![0u8; (w * h * 3) as usize];
    let m_per_px = 2.0 * PI * radius / w as f64;
    for y in 0..h {
        for x in 0..w {
            let k = (y * w + x) as usize;
            let v = hs[k];
            let c: [f32; 3] = if v <= 0.0 {
                let t = (-v / 6000.0).clamp(0.0, 1.0);
                [0.15 - 0.1 * t, 0.35 - 0.25 * t, 0.55 - 0.3 * t]
            } else {
                let t = (v / 4500.0).clamp(0.0, 1.0);
                let lo = [0.33, 0.48, 0.25];
                let mid = [0.62, 0.55, 0.36];
                let hi = [0.92, 0.92, 0.92];
                if t < 0.5 {
                    lerp3(lo, mid, t * 2.0)
                } else {
                    lerp3(mid, hi, (t - 0.5) * 2.0)
                }
            };
            // Hillshade (NW light), 8× exaggeration for legibility at this scale.
            let xl = (x + w - 1) % w;
            let xr = (x + 1) % w;
            let yu = y.saturating_sub(1);
            let yd = (y + 1).min(h - 1);
            let lat = PI / 2.0 - (y as f64 + 0.5) / h as f64 * PI;
            let dzdx = (hs[(y * w + xr) as usize] - hs[(y * w + xl) as usize]) as f64 / (2.0 * m_per_px * lat.cos().max(0.05));
            let dzdy = (hs[(yu * w + x) as usize] - hs[(yd * w + x) as usize]) as f64 / (2.0 * m_per_px);
            let ex = 8.0;
            let n = DVec3::new(-dzdx * ex, -dzdy * ex, 1.0).normalize();
            let l = DVec3::new(-0.6, 0.6, 0.7).normalize();
            let shade = if v > 0.0 { (n.dot(l) / l.z).clamp(0.3, 1.6) as f32 } else { 1.0 };
            for ch in 0..3 {
                img[k * 3 + ch] = ((c[ch] * shade).clamp(0.0, 1.0) * 255.0) as u8;
            }
        }
    }
    if let Some(er) = er {
        for s in &er.rivers {
            let width = ((s.area_km2 as f64 / 2.0e4).sqrt() * 0.6).clamp(0.6, 3.0);
            draw_line(&mut img, w, h, s.from, s.to, width, [40, 90, 200]);
        }
    }
    let file = std::fs::File::create(path).unwrap();
    let mut enc = png::Encoder::new(std::io::BufWriter::new(file), w, h);
    enc.set_color(png::ColorType::Rgb);
    enc.set_depth(png::BitDepth::Eight);
    enc.write_header().unwrap().write_image_data(&img).unwrap();
}

fn lerp3(a: [f32; 3], b: [f32; 3], t: f32) -> [f32; 3] {
    [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t]
}

fn to_px(d: DVec3, w: u32, h: u32) -> (f64, f64) {
    let lat = d.z.clamp(-1.0, 1.0).asin();
    let lon = d.y.atan2(d.x);
    ((lon + PI) / (2.0 * PI) * w as f64, (PI / 2.0 - lat) / PI * h as f64)
}

fn draw_line(img: &mut [u8], w: u32, h: u32, a: DVec3, b: DVec3, width: f64, c: [u8; 3]) {
    let (x0, y0) = to_px(a, w, h);
    let (mut x1, y1) = to_px(b, w, h);
    if (x1 - x0).abs() > w as f64 / 2.0 {
        x1 += if x1 < x0 { w as f64 } else { -(w as f64) };
    }
    let steps = ((x1 - x0).abs().max((y1 - y0).abs()) * 2.0).ceil().max(1.0) as i32;
    let r = width * 0.5;
    for s in 0..=steps {
        let t = s as f64 / steps as f64;
        let (x, y) = (x0 + (x1 - x0) * t, y0 + (y1 - y0) * t);
        let ri = r.ceil() as i32;
        for dy in -ri..=ri {
            for dx in -ri..=ri {
                let (fx, fy) = (dx as f64, dy as f64);
                if fx * fx + fy * fy > r * r + 0.25 {
                    continue;
                }
                let px = ((x + fx).rem_euclid(w as f64)) as u32;
                let py = (y + fy).clamp(0.0, h as f64 - 1.0) as u32;
                let k = ((py * w + px) * 3) as usize;
                img[k..k + 3].copy_from_slice(&c);
            }
        }
    }
}
