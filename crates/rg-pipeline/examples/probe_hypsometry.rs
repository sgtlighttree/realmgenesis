//! Global land hypsometry of the final terrain vs. the macro surface (the hypsometric remap
//! targets a smoothed, ice-free Earth curve: mean 656 m, median 382 m):
//! `cargo run -p rg-pipeline --release --example probe_hypsometry -- <seed>`
//! Clustered sampling: 150 random land centres × 40 points in a shared 4 km patch each, so
//! the network detail builds 150 patches instead of 6000.
#[path = "common/mod.rs"]
mod common;

fn main() {
    let seed: u64 = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(3);
    let w = common::world(seed);
    let s = rg_pipeline::terrain_sampler(&w);
    let r = w.params.tectonics.radius_m;
    let mut rng = 99u64;
    let mut next = || { rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); (rng >> 11) as f64 / (1u64 << 53) as f64 };
    let (mut macro_h, mut full) = (Vec::new(), Vec::new());
    while full.len() < 6000 {
        // A cluster centre, uniform on the sphere.
        let z = next() * 2.0 - 1.0;
        let t = next() * std::f64::consts::TAU;
        let rr = (1.0 - z * z).sqrt();
        let c = glam::DVec3::new(rr * t.cos(), rr * t.sin(), z);
        if s.macro_height(c) <= 0.0 {
            continue;
        }
        let east = glam::DVec3::Z.cross(c).normalize_or(glam::DVec3::X);
        let north = c.cross(east);
        let patch = s.network_patch(c, 2900.0, 100.0);
        for _ in 0..40 {
            let (x, y) = ((next() - 0.5) * 4000.0, (next() - 0.5) * 4000.0);
            let d = (c + east * (x / r) + north * (y / r)).normalize();
            let m = s.macro_height(d);
            if m > 0.0 {
                macro_h.push(m);
                full.push(match &patch {
                    Some(p) => s.levels_in_patch(p, d).h,
                    None => s.height(d, 100.0).0,
                });
            }
        }
    }
    // Paired: detail's mean effect (low variance, unlike comparing two noisy means).
    let diff = full.iter().zip(&macro_h).map(|(f, m)| f - m).sum::<f64>() / full.len() as f64;
    println!("detail − macro: mean {diff:+.0} m (should be ≈ 0: the macro already carries the calibrated hypsometry)");
    for (name, v) in [("macro", &mut macro_h), ("final", &mut full)] {
        v.sort_by(f64::total_cmp);
        let q = |f: f64| v[((v.len() - 1) as f64 * f) as usize];
        let mean = v.iter().sum::<f64>() / v.len() as f64;
        let above = |h: f64| v.iter().filter(|&&x| x > h).count() as f64 / v.len() as f64 * 100.0;
        println!("{name:6} land: mean {mean:5.0} m  p50 {:5.0}  p90 {:5.0}  p99 {:5.0}  max {:5.0}  | >1 km {:4.1}%  >2 km {:4.1}%  >4 km {:4.1}%", q(0.5), q(0.9), q(0.99), q(1.0), above(1000.0), above(2000.0), above(4000.0));
    }
    println!("target land: mean   656 m  p50  382  (smoothed, ice-free Earth curve used by the remap)");
}
