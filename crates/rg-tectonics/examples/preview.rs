//! Writes equirectangular preview PNGs of a tectonic snapshot.
//!
//! Usage: `cargo run -p rg-tectonics --example preview --release -- <seed> <out_dir> [key=value ...]`
//! Optional keys: `res`, `plates`, `land`, `continents`, `mountain`, `rough`, `w` (image width).
//! Writes `height.png`, `plates.png`, `crust.png` into `<out_dir>`.

use glam::DVec3;
use rg_tectonics::{TectonicsParams, generate};
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let seed: u64 = args.first().and_then(|s| s.parse().ok()).unwrap_or(1);
    let out_dir = args.get(1).cloned().unwrap_or_else(|| ".".to_string());
    let mut p = TectonicsParams { seed, ..Default::default() };
    let mut width = 1024usize;
    for kv in args.iter().skip(2) {
        let Some((k, v)) = kv.split_once('=') else { continue };
        match k {
            "res" => p.res = v.parse().unwrap(),
            "plates" => p.plate_count = v.parse().unwrap(),
            "land" => p.land_fraction = v.parse().unwrap(),
            "continents" => p.continent_count_hint = v.parse().unwrap(),
            "mountain" => p.mountain_scale = v.parse().unwrap(),
            "rough" => p.roughness = v.parse().unwrap(),
            "w" => width = v.parse().unwrap(),
            _ => eprintln!("unknown key {k}"),
        }
    }
    let height = width / 2;
    std::fs::create_dir_all(&out_dir).unwrap();

    let t0 = Instant::now();
    let snap = generate(&p);
    let dt = t0.elapsed();
    let (lo, hi) = snap.height_m.min_max();
    println!(
        "seed {seed}: generate {:.0} ms, land {:.2}%, offset {:.0} m, height {lo:.0}..{hi:.0} m, {} plates",
        dt.as_secs_f64() * 1000.0,
        snap.land_fraction * 100.0,
        snap.sea_level_offset_m,
        snap.plates.len()
    );

    let mut land: Vec<f32> = snap.height_m.data.iter().copied().filter(|&v| v > 0.0).collect();
    land.sort_by(|a, b| a.total_cmp(b));
    let pct = |q: f64| land[((land.len() - 1) as f64 * q) as usize];
    let mut sea: Vec<f32> = snap.height_m.data.iter().copied().filter(|&v| v < 0.0).collect();
    sea.sort_by(|a, b| a.total_cmp(b));
    let spct = |q: f64| sea[((sea.len() - 1) as f64 * q) as usize];
    println!(
        "  land m: p25 {:.0} p50 {:.0} p75 {:.0} p90 {:.0} p99 {:.0} | sea m: p1 {:.0} p10 {:.0} p50 {:.0} p90 {:.0}",
        pct(0.25), pct(0.5), pct(0.75), pct(0.9), pct(0.99), spct(0.01), spct(0.1), spct(0.5), spct(0.9)
    );

    // Sample the fields on an equirectangular grid.
    let dir_at = |x: usize, y: usize| -> DVec3 {
        let lon = ((x as f64 + 0.5) / width as f64 - 0.5) * std::f64::consts::TAU;
        let lat = (0.5 - (y as f64 + 0.5) / height as f64) * std::f64::consts::PI;
        DVec3::new(lat.cos() * lon.cos(), lat.cos() * lon.sin(), lat.sin())
    };
    let mut h = vec![0f32; width * height];
    for y in 0..height {
        for x in 0..width {
            h[y * width + x] = snap.height_m.sample(dir_at(x, y));
        }
    }

    // --- height.png: hypsometric tint + hillshade
    let mut img = vec![0u8; width * height * 3];
    let m_per_px = p.radius_m * std::f64::consts::TAU / width as f64;
    for y in 0..height {
        let lat = (0.5 - (y as f64 + 0.5) / height as f64) * std::f64::consts::PI;
        let dx = m_per_px * lat.cos().max(0.05);
        let dy = m_per_px;
        for x in 0..width {
            let at = |xx: isize, yy: isize| {
                let xx = xx.rem_euclid(width as isize) as usize;
                let yy = yy.clamp(0, height as isize - 1) as usize;
                h[yy * width + xx].max(-200.0) as f64 // flatten the sea for shading
            };
            let (xi, yi) = (x as isize, y as isize);
            let dzdx = (at(xi + 1, yi) - at(xi - 1, yi)) / (2.0 * dx);
            let dzdy = (at(xi, yi - 1) - at(xi, yi + 1)) / (2.0 * dy);
            let ex = 40.0;
            // Light from the north-west.
            let (nx, ny, nz) = (-dzdx * ex, -dzdy * ex, 1.0);
            let nl = (nx * nx + ny * ny + nz * nz).sqrt();
            let light = ((-0.55 * nx + 0.55 * ny + 0.63 * nz) / nl).clamp(0.0, 1.0);
            let shade = 0.55 + 0.6 * (light - 0.63);
            let hv = h[y * width + x] as f64;
            let mut c = hypsometric(hv);
            let s = if hv > 0.0 { shade } else { 1.0 };
            for v in &mut c {
                *v = (*v * s).clamp(0.0, 255.0);
            }
            let o = (y * width + x) * 3;
            img[o] = c[0] as u8;
            img[o + 1] = c[1] as u8;
            img[o + 2] = c[2] as u8;
        }
    }
    write_png(&Path::new(&out_dir).join("height.png"), width, height, &img);
    let height_img = img.clone();

    // --- plates.png
    let mut img = vec![0u8; width * height * 3];
    for y in 0..height {
        for x in 0..width {
            let d = dir_at(x, y);
            let id = snap.plate_id.nearest(d) as u32;
            let hash = id.wrapping_mul(2654435761).wrapping_add(12345);
            let mut c = [
                80.0 + ((hash >> 8) & 0x7f) as f64,
                80.0 + ((hash >> 16) & 0x7f) as f64,
                80.0 + ((hash >> 24) & 0x7f) as f64,
            ];
            // Darken the ocean slightly so coasts are readable.
            if h[y * width + x] < 0.0 {
                for v in &mut c {
                    *v *= 0.75;
                }
            }
            let b = snap.boundary.nearest(d);
            if b != 0 {
                c = match b {
                    1 => [230.0, 30.0, 30.0],
                    2 => [30.0, 90.0, 255.0],
                    _ => [30.0, 220.0, 60.0],
                };
            }
            let o = (y * width + x) * 3;
            img[o] = c[0] as u8;
            img[o + 1] = c[1] as u8;
            img[o + 2] = c[2] as u8;
        }
    }
    write_png(&Path::new(&out_dir).join("plates.png"), width, height, &img);
    let mut both = height_img;
    both.extend_from_slice(&img);
    write_png(&Path::new(&out_dir).join("height_plates.png"), width, height * 2, &both);

    // --- crust.png
    let mut img = vec![0u8; width * height * 3];
    for y in 0..height {
        for x in 0..width {
            let cont = snap.crust.nearest(dir_at(x, y)) != 0;
            let c = if cont { [196u8, 170, 120] } else { [30u8, 60, 110] };
            let o = (y * width + x) * 3;
            img[o..o + 3].copy_from_slice(&c);
        }
    }
    write_png(&Path::new(&out_dir).join("crust.png"), width, height, &img);
}

fn mix(a: [f64; 3], b: [f64; 3], t: f64) -> [f64; 3] {
    [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t, a[2] + (b[2] - a[2]) * t]
}

/// Piecewise-linear colour ramp over (height, rgb) stops.
fn ramp(stops: &[(f64, [f64; 3])], h: f64) -> [f64; 3] {
    if h <= stops[0].0 {
        return stops[0].1;
    }
    for w in stops.windows(2) {
        if h <= w[1].0 {
            return mix(w[0].1, w[1].1, (h - w[0].0) / (w[1].0 - w[0].0));
        }
    }
    stops[stops.len() - 1].1
}

fn hypsometric(h: f64) -> [f64; 3] {
    const SEA: [(f64, [f64; 3]); 8] = [
        (-10000.0, [4.0, 8.0, 40.0]),
        (-7000.0, [8.0, 20.0, 80.0]),
        (-5500.0, [14.0, 40.0, 120.0]),
        (-4000.0, [28.0, 75.0, 165.0]),
        (-2500.0, [55.0, 125.0, 205.0]),
        (-1000.0, [85.0, 155.0, 220.0]),
        (-200.0, [115.0, 185.0, 228.0]),
        (0.0, [140.0, 205.0, 235.0]),
    ];
    const LAND: [(f64, [f64; 3]); 8] = [
        (0.0, [70.0, 145.0, 75.0]),
        (400.0, [110.0, 170.0, 80.0]),
        (1000.0, [185.0, 185.0, 95.0]),
        (2000.0, [175.0, 140.0, 85.0]),
        (3200.0, [140.0, 105.0, 75.0]),
        (4500.0, [120.0, 95.0, 90.0]),
        (6000.0, [215.0, 210.0, 210.0]),
        (8000.0, [255.0, 255.0, 255.0]),
    ];
    if h < 0.0 { ramp(&SEA, h) } else { ramp(&LAND, h) }
}

fn write_png(path: &Path, w: usize, h: usize, rgb: &[u8]) {
    let file = BufWriter::new(File::create(path).unwrap());
    let mut enc = png::Encoder::new(file, w as u32, h as u32);
    enc.set_color(png::ColorType::Rgb);
    enc.set_depth(png::BitDepth::Eight);
    let mut writer = enc.write_header().unwrap();
    writer.write_image_data(rgb).unwrap();
    println!("wrote {}", path.display());
}
