//! Closed-depression (pit) density of the detail terrain in a box:
//! `cargo run -p rg-pipeline --release --example probe_pits -- <seed> <lat> <lon> <half_km> <spacing_m>`
//! A pit = a land cell lower than all 8 neighbours. Natural fluvial terrain has very few
//! at 0.5–1 km spacing; raw noise terrain has many.
use rg_pipeline::{generate, WorldParams};

fn main() {
    let a: Vec<f64> = std::env::args().skip(1).map(|s| s.parse().unwrap()).collect();
    let (seed, lat0, lon0, half, sp) = (a[0] as u64, a[1].to_radians(), a[2].to_radians(), a[3] * 1000.0, a[4]);
    let mut p = WorldParams::default();
    p.tectonics.seed = seed;
    let w = generate(&p);
    let s = rg_pipeline::terrain_sampler(&w);
    let r = p.tectonics.radius_m;
    let n = (2.0 * half / sp) as usize;
    let mut h = vec![0.0; n * n];
    for j in 0..n {
        for i in 0..n {
            let (x, y) = (-half + i as f64 * sp, -half + j as f64 * sp);
            let (la, lo) = (lat0 + y / r, lon0 + x / (r * lat0.cos()));
            h[j * n + i] = s.height(glam::DVec3::new(la.cos() * lo.cos(), la.cos() * lo.sin(), la.sin()), sp).0;
        }
    }
    let (mut pits, mut land) = (0usize, 0usize);
    for j in 1..n - 1 {
        for i in 1..n - 1 {
            let c = h[j * n + i];
            if c <= 0.0 {
                continue;
            }
            land += 1;
            let lowest = (-1i64..=1).flat_map(|dj| (-1i64..=1).map(move |di| (di, dj))).filter(|&(di, dj)| di != 0 || dj != 0).all(|(di, dj)| h[(j as i64 + dj) as usize * n + (i as i64 + di) as usize] > c);
            if lowest {
                pits += 1;
            }
        }
    }
    println!("pits {pits} / {land} land cells = {:.2}% at {sp} m spacing", pits as f64 / land.max(1) as f64 * 100.0);
}
