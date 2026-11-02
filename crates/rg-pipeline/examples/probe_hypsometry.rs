//! Global land hypsometry of the final terrain vs. Earth:
//! `cargo run -p rg-pipeline --release --example probe_hypsometry -- <seed>`
use rg_pipeline::{generate, WorldParams};
use std::sync::Arc;

fn main() {
    let seed: u64 = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(3);
    let mut p = WorldParams::default();
    p.tectonics.seed = seed;
    let w = generate(&p);
    let s = rg_terrain::TerrainSampler::new(Arc::new(w.tectonics.height_m.clone()), p.tectonics.radius_m, seed);
    let mut rng = 99u64;
    let mut next = || { rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); (rng >> 11) as f64 / (1u64 << 53) as f64 };
    let (mut macro_h, mut full) = (Vec::new(), Vec::new());
    while full.len() < 6000 {
        // Uniform on the sphere.
        let z = next() * 2.0 - 1.0;
        let t = next() * std::f64::consts::TAU;
        let r = (1.0 - z * z).sqrt();
        let d = glam::DVec3::new(r * t.cos(), r * t.sin(), z);
        let m = s.macro_height(d);
        if m > 0.0 {
            macro_h.push(m);
            full.push(s.height(d, 100.0).0);
        }
    }
    for (name, v) in [("macro", &mut macro_h), ("final", &mut full)] {
        v.sort_by(f64::total_cmp);
        let q = |f: f64| v[((v.len() - 1) as f64 * f) as usize];
        let mean = v.iter().sum::<f64>() / v.len() as f64;
        let above = |h: f64| v.iter().filter(|&&x| x > h).count() as f64 / v.len() as f64 * 100.0;
        println!("{name:6} land: mean {mean:5.0} m  p50 {:5.0}  p90 {:5.0}  p99 {:5.0}  max {:5.0}  | >1 km {:4.1}%  >2 km {:4.1}%  >4 km {:4.1}%", q(0.5), q(0.9), q(0.99), q(1.0), above(1000.0), above(2000.0), above(4000.0));
    }
    println!("Earth  land: mean   840 m  p50 ~440  (ETOPO; ~25% > 1 km, ~9% > 2 km, ~2% > 4 km, max 8848)");
}
