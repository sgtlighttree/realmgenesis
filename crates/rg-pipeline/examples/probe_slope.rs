//! Slope distribution in a box, compared with Earth: `cargo run -p rg-pipeline --release --example probe_slope -- <seed> <lat> <lon> <half_km>`
//! Samples 2000 random points, central differences at 90 m (SRTM-like), reports degrees.
use rg_pipeline::{generate, WorldParams};

fn main() {
    let a: Vec<f64> = std::env::args().skip(1).map(|s| s.parse().unwrap()).collect();
    let (seed, lat0, lon0, half) = (a[0] as u64, a[1].to_radians(), a[2].to_radians(), a[3] * 1000.0);
    let mut p = WorldParams::default();
    p.tectonics.seed = seed;
    if let Some(r) = std::env::var("RG_RES").ok().and_then(|v| v.parse().ok()) {
        p.tectonics.res = r;
    }
    if let Some(r) = std::env::var("RG_EROSION_RES").ok().and_then(|v| v.parse().ok()) {
        p.erosion.max_res = r;
    }
    let t0 = std::time::Instant::now();
    let w = generate(&p);
    eprintln!("generate {:.1} s (tectonics {:.0} ms, erosion {:.0} ms)", t0.elapsed().as_secs_f64(), w.timings.tectonics_ms, w.timings.erosion_ms);
    let s = rg_pipeline::terrain_sampler(&w);
    let r = p.tectonics.radius_m;
    let at = |x: f64, y: f64| {
        let (la, lo) = (lat0 + y / r, lon0 + x / (r * lat0.cos()));
        s.height(glam::DVec3::new(la.cos() * lo.cos(), la.cos() * lo.sin(), la.sin()), 40.0).0
    };
    let detail = |x: f64, y: f64| {
        let (la, lo) = (lat0 + y / r, lon0 + x / (r * lat0.cos()));
        let d = glam::DVec3::new(la.cos() * lo.cos(), la.cos() * lo.sin(), la.sin());
        s.height(d, 40.0).0 - s.macro_height(d)
    };
    let mut rng = 12345u64;
    let mut next = || { rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); (rng >> 11) as f64 / (1u64 << 53) as f64 };
    let d = 45.0;
    let mut deg: Vec<f64> = (0..2000)
        .map(|_| {
            let (x, y) = ((next() * 2.0 - 1.0) * half, (next() * 2.0 - 1.0) * half);
            let gx = (at(x + d, y) - at(x - d, y)) / (2.0 * d);
            let gy = (at(x, y + d) - at(x, y - d)) / (2.0 * d);
            (gx * gx + gy * gy).sqrt().atan().to_degrees()
        })
        .collect();
    deg.sort_by(f64::total_cmp);
    let mut dd: Vec<f64> = (0..400).map(|i| detail((i % 20) as f64 / 20.0 * 2.0 * half - half, (i / 20) as f64 / 20.0 * 2.0 * half - half)).collect();
    dd.sort_by(f64::total_cmp);
    println!("detail above macro: p1 {:.0} m  p50 {:.0} m  p99 {:.0} m", dd[4], dd[200], dd[395]);
    let q = |f: f64| deg[((deg.len() - 1) as f64 * f) as usize];
    let mean = deg.iter().sum::<f64>() / deg.len() as f64;
    println!("slope° mean {mean:.1}  p50 {:.1}  p90 {:.1}  p99 {:.1}  max {:.1}", q(0.5), q(0.9), q(0.99), q(1.0));
}
