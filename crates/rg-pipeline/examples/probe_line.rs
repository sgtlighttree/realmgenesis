//! Probe terrain continuity: `cargo run -p rg-pipeline --release --example probe_line -- <seed> <lat> <lon>`
use rg_pipeline::{generate, WorldParams};
use std::sync::Arc;

fn main() {
    let a: Vec<f64> = std::env::args().skip(1).map(|s| s.parse().unwrap()).collect();
    let (seed, lat, lon) = (a[0] as u64, a[1].to_radians(), a[2].to_radians());
    let mut p = WorldParams::default();
    p.tectonics.seed = seed;
    let w = generate(&p);
    let s = rg_terrain::TerrainSampler::new(Arc::new(w.tectonics.height_m.clone()), p.tectonics.radius_m, seed);
    let r = p.tectonics.radius_m;
    let dir = |la: f64, lo: f64| glam::DVec3::new(la.cos() * lo.cos(), la.cos() * lo.sin(), la.sin());
    for min_wl in [20.0, 200.0, 2000.0] {
        let step = 5.0;
        let n = 8000; // 40 km line, eastward
        let mut prev = None;
        let mut worst = (0.0f64, 0usize);
        for i in 0..n {
            let lo = lon + (i as f64 * step) / (r * lat.cos());
            let h = s.levels(dir(lat, lo), min_wl).h;
            if let Some(ph) = prev {
                let d: f64 = h - ph;
                if d.abs() > worst.0 {
                    worst = (d.abs(), i);
                }
            }
            prev = Some(h);
        }
        // Slope stats over 25 m baselines.
        let hs: Vec<f64> = (0..n).map(|i| s.levels(dir(lat, lon + (i as f64 * step) / (r * lat.cos())), min_wl).h).collect();
        let mut sl: Vec<f64> = hs.windows(6).map(|w| (w[5] - w[0]).abs() / 25.0).collect();
        sl.sort_by(f64::total_cmp);
        let q = |f: f64| sl[((sl.len() - 1) as f64 * f) as usize];
        println!("  slope (25 m) p50 {:.2} p90 {:.2} p99 {:.2} max {:.2}  (45° = 1.0)", q(0.5), q(0.9), q(0.99), q(1.0));
        println!("min_wl {min_wl:6.0} m: worst 5 m jump {:.1} m at {:.2} km", worst.0, worst.1 as f64 * step / 1000.0);
        if min_wl == 20.0 {
            let x0 = worst.1 as f64 * step - step;
            for sub in [1.0, 0.1, 0.01, 0.001] {
                let mut wj = (0.0f64, 0.0);
                let mut prev = None;
                for t in 0..((2.0 * step / sub) as usize) {
                    let x = x0 + t as f64 * sub;
                    let h = s.levels(dir(lat, lon + x / (r * lat.cos())), 20.0).h;
                    if let Some(p) = prev { let d: f64 = h - p; if d.abs() > wj.0 { wj = (d.abs(), x); } }
                    prev = Some(h);
                }
                println!("   step {sub} m: worst jump {:.4} m at x {:.3} m", wj.0, wj.1);
                if sub == 0.001 {
                    let mut bj = 0.0f64; let mut mj = 0.0f64;
                    let mut pb: Option<f64> = None; let mut pm: Option<f64> = None;
                    for t in 0..((2.0 * step / sub) as usize) {
                        let x = x0 + t as f64 * sub;
                        let d = dir(lat, lon + x / (r * lat.cos()));
                        let (b, m) = (s.base_surface(d), s.macro_height(d));
                        if let Some(p) = pb { bj = bj.max((b - p).abs()); }
                        if let Some(p) = pm { mj = mj.max((m - p).abs()); }
                        pb = Some(b); pm = Some(m);
                    }
                    println!("   base max jump {bj:.5} m, macro max jump {mj:.5} m");

                }
            }
        }
    }
}
