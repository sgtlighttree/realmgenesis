//! Ocean hypsometry of the committed surface against Earth:
//! `cargo run -p rg-pipeline --release --example probe_bathymetry -- <seed>`
//! (same `RG_*` overrides and cache as `probe_realism`).
//!
//! Earth reference, share of ocean floor per depth band (Costello et al. 2010, ETOPO-class
//! data; mean depth 3688 m from ETOPO1, NOAA NGDC): 0–200 m 7.3%, 200–1000 m 4.4%,
//! 1–4 km 35.2%, 4–6 km 52.1%, 6–7 km 0.9%, > 7 km 0.13%. Ocean = 71% of the surface.
#[path = "common/mod.rs"]
mod common;
use glam::DVec3;

const EARTH: [(f64, f64, f64); 6] =
    [(0.0, 200.0, 7.30), (200.0, 1000.0, 4.38), (1000.0, 4000.0, 35.22), (4000.0, 6000.0, 52.08), (6000.0, 7000.0, 0.89), (7000.0, 1.0e9, 0.13)];

fn main() {
    let seed: u64 = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(3);
    let w = common::world(seed);
    let h = &w.tectonics.height_m;
    // Fibonacci sphere: equal-area samples.
    let n = 2_000_000usize;
    let golden = std::f64::consts::PI * (3.0 - 5f64.sqrt());
    let mut depths: Vec<f64> = Vec::new();
    for i in 0..n {
        let z = 1.0 - 2.0 * (i as f64 + 0.5) / n as f64;
        let r = (1.0 - z * z).sqrt();
        let t = golden * i as f64;
        let v = h.sample_drained(DVec3::new(r * t.cos(), r * t.sin(), z));
        if v <= 0.0 {
            depths.push(-v);
        }
    }
    let ocean = depths.len() as f64 / n as f64;
    depths.sort_by(f64::total_cmp);
    let mean = depths.iter().sum::<f64>() / depths.len() as f64;
    let q = |f: f64| depths[((depths.len() - 1) as f64 * f) as usize];
    println!("ocean {:.1}% of surface (Earth 71%) · mean depth {mean:.0} m (Earth 3688) · median {:.0} m · p99 {:.0} m · max {:.0} m", ocean * 100.0, q(0.5), q(0.99), q(1.0));
    println!("band            ours    Earth");
    for (lo, hi, e) in EARTH {
        let c = depths.iter().filter(|&&d| d > lo && d <= hi).count() as f64 / depths.len() as f64 * 100.0;
        let label = if hi > 1e8 { format!("> {:.0} km", lo / 1000.0) } else { format!("{lo:.0}–{hi:.0} m") };
        println!("{label:<14} {c:5.1}%  {e:5.1}%");
    }
    // Finer histogram (500 m bins) to see the shape: Earth's mode is 4–5 km.
    let mut line = String::from("500 m bins:");
    for b in 0..16 {
        let (lo, hi) = (b as f64 * 500.0, (b + 1) as f64 * 500.0);
        let c = depths.iter().filter(|&&d| d > lo && d <= hi).count() as f64 / depths.len() as f64 * 100.0;
        line += &format!(" {c:.0}");
    }
    println!("{line}");
}
