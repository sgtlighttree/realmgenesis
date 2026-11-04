//! The fields behind one spot, side by side, to find which one carries a feature seen in the
//! app (s03: a stair-stepped escarpment). Writes one PNG of 3 × 2 panels:
//!
//! | hillshade (displayed surface) | uplift (mm/yr)    | relative runoff (erosion's) |
//! | lake depth (m)                | plate id + bounds | crust type                  |
//!
//! `cargo run -p rg-pipeline --release --example field_maps -- <seed> <lat> <lon> <half_km> <out.png>`
//! Same `RG_*` overrides and cache as the other probes (`RG_RES=1024` = the app's Ultra).
//! Every panel samples its field at the nearest node, so cell edges show as they are.
#[path = "common/mod.rs"]
mod common;
use glam::DVec3;
use rg_geo::CubeField;

const PANEL: usize = 600;

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let seed: u64 = a[0].parse().unwrap();
    let lat0 = a[1].parse::<f64>().unwrap().to_radians();
    let lon0 = a[2].parse::<f64>().unwrap().to_radians();
    let half = a[3].parse::<f64>().unwrap() * 1000.0;
    let out = &a[4];
    let w = common::world(seed);
    let r = w.params.tectonics.radius_m;
    let t = &w.tectonics;
    let hy = w.hydrology.as_ref().expect("erosion on");
    let sampler = rg_pipeline::terrain_sampler(&w);
    eprintln!("computing runoff (reruns tier A climate on the tectonic surface)");
    let runoff = rg_pipeline::relative_runoff_for_probe(&w.params.planet, t);

    let m_per_px = 2.0 * half / PANEL as f64;
    let at = |i: usize, j: usize| {
        let (x, y) = (-half + (i as f64 + 0.5) * m_per_px, half - (j as f64 + 0.5) * m_per_px);
        let (la, lo) = (lat0 + y / r, lon0 + x / (r * lat0.cos()));
        DVec3::new(la.cos() * lo.cos(), la.cos() * lo.sin(), la.sin())
    };

    let mut img = vec![0u8; 3 * PANEL * 2 * PANEL * 3];
    let mut put = |panel: usize, i: usize, j: usize, c: [u8; 3]| {
        let (px, py) = (panel % 3 * PANEL + i, panel / 3 * PANEL + j);
        let k = (py * 3 * PANEL + px) * 3;
        img[k..k + 3].copy_from_slice(&c);
    };

    // Value ranges inside the box, for the colour ramps.
    let mut ranges = [(f64::MAX, f64::MIN); 3];
    for j in (0..PANEL).step_by(4) {
        for i in (0..PANEL).step_by(4) {
            let d = at(i, j);
            let vals = [
                hy.uplift_mm_yr.as_ref().map_or(0.0, |u| nearest(u, d) as f64),
                nearest(&runoff, d) as f64,
                nearest(&hy.lake_depth_m, d) as f64,
            ];
            for (rg, v) in ranges.iter_mut().zip(vals) {
                *rg = (rg.0.min(v), rg.1.max(v));
            }
        }
    }
    eprintln!("uplift {:?} mm/yr · runoff {:?} (× land mean) · lake depth {:?} m", ranges[0], ranges[1], ranges[2]);

    for j in 0..PANEL {
        for i in 0..PANEL {
            let d = at(i, j);
            // 0: hillshade of the displayed surface, coloured by height.
            let h = sampler.macro_height(d);
            let e = 2.0 * m_per_px;
            let hx = |dx: f64, dy: f64| {
                let (la, lo) = (d.z.asin() + dy / r, d.y.atan2(d.x) + dx / (r * lat0.cos()));
                sampler.macro_height(DVec3::new(la.cos() * lo.cos(), la.cos() * lo.sin(), la.sin()))
            };
            let gx = (hx(e, 0.0) - hx(-e, 0.0)) / (2.0 * e);
            let gy = (hx(0.0, e) - hx(0.0, -e)) / (2.0 * e);
            let n = DVec3::new(-gx * 8.0, -gy * 8.0, 1.0).normalize();
            let shade = n.dot(DVec3::new(-0.5, 0.5, 0.7).normalize()).max(0.0);
            let c = if h <= 0.0 {
                [20.0, 40.0, 80.0]
            } else {
                let t = (h / 4000.0).min(1.0);
                [90.0 + 120.0 * t, 110.0 + 80.0 * t, 70.0 + 110.0 * t]
            };
            put(0, i, j, c.map(|v| (v * (0.35 + 0.85 * shade)).min(255.0) as u8));
            // 1–3: scalar fields on a ramp.
            let up = hy.uplift_mm_yr.as_ref().map_or(0.0, |u| nearest(u, d) as f64);
            put(1, i, j, ramp(up, ranges[0]));
            put(2, i, j, ramp(nearest(&runoff, d) as f64, ranges[1]));
            put(3, i, j, ramp(nearest(&hy.lake_depth_m, d) as f64, ranges[2]));
            // 4: plate id as a hashed colour, boundaries dark.
            let pid = nearest(&t.plate_id, d) as u64;
            let mut c = hash_color(pid);
            if nearest(&t.boundary, d) != 0 {
                c = c.map(|v| v / 3);
            }
            put(4, i, j, c);
            // 5: crust type as a hashed colour.
            put(5, i, j, hash_color(1000 + nearest(&t.crust, d) as u64));
        }
    }
    let file = std::fs::File::create(out).unwrap();
    let mut enc = png::Encoder::new(std::io::BufWriter::new(file), 3 * PANEL as u32, 2 * PANEL as u32);
    enc.set_color(png::ColorType::Rgb);
    enc.write_header().unwrap().write_image_data(&img).unwrap();
    eprintln!("wrote {out}");
}

/// The value of the node nearest `d` (no interpolation, so cell edges stay visible).
fn nearest<T: Copy>(f: &CubeField<T>, d: DVec3) -> T {
    let (face, u, v) = rg_geo::CubeFace::from_dir(d);
    let r = f.res as f64;
    let i = ((u + 1.0) * 0.5 * r).round().clamp(0.0, r) as u32;
    let j = ((v + 1.0) * 0.5 * r).round().clamp(0.0, r) as u32;
    f.data[f.index(face, i, j)]
}

/// Dark blue → yellow over `range`.
fn ramp(v: f64, range: (f64, f64)) -> [u8; 3] {
    let t = ((v - range.0) / (range.1 - range.0).max(1e-9)).clamp(0.0, 1.0);
    [(30.0 + 220.0 * t) as u8, (30.0 + 200.0 * t) as u8, (120.0 - 90.0 * t) as u8]
}

fn hash_color(k: u64) -> [u8; 3] {
    let h = k.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ (k >> 7);
    [(h >> 8) as u8 / 2 + 80, (h >> 24) as u8 / 2 + 80, (h >> 40) as u8 / 2 + 80]
}
