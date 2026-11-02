//! Find sharp features near a point: `cargo run -p rg-pipeline --release --example probe_grid -- <seed> <lat> <lon> <half_km>`
//! Scans a grid at 200 m spacing and reports the steepest spots per component.
use rg_pipeline::{generate, WorldParams};
use std::sync::Arc;

fn main() {
    let a: Vec<f64> = std::env::args().skip(1).map(|s| s.parse().unwrap()).collect();
    let (seed, lat0, lon0, half) = (a[0] as u64, a[1].to_radians(), a[2].to_radians(), a[3] * 1000.0);
    let mut p = WorldParams::default();
    p.tectonics.seed = seed;
    let w = generate(&p);
    let s = rg_terrain::TerrainSampler::new(Arc::new(w.tectonics.height_m.clone()), p.tectonics.radius_m, seed);
    let r = p.tectonics.radius_m;
    let dir = |x: f64, y: f64| {
        let (la, lo) = (lat0 + y / r, lon0 + x / (r * lat0.cos()));
        glam::DVec3::new(la.cos() * lo.cos(), la.cos() * lo.sin(), la.sin())
    };
    let step = 200.0;
    let n = (2.0 * half / step) as i64;
    let mut worst = [(0.0f64, 0.0, 0.0); 3];
    for j in 0..n {
        for i in 0..n {
            let (x, y) = (-half + i as f64 * step, -half + j as f64 * step);
            let comp = |d: glam::DVec3| (s.macro_height(d), s.base_surface(d), s.height(d, 50.0).0);
            let c0 = comp(dir(x, y));
            let c1 = comp(dir(x + step, y));
            for (k, (a, b)) in [(c0.0, c1.0), (c0.1, c1.1), (c0.2, c1.2)].into_iter().enumerate() {
                let sl = (b - a).abs() / step;
                if sl > worst[k].0 {
                    worst[k] = (sl, x, y);
                }
            }
        }
    }
    for (k, name) in ["macro", "macro+hills", "final"].iter().enumerate() {
        println!("{name:12}: max slope {:.3} at ({:.1} km, {:.1} km)", worst[k].0, worst[k].1 / 1000.0, worst[k].2 / 1000.0);
    }
}
