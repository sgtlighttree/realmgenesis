//! Hydrological consistency of the final surface in a box:
//! `cargo run -p rg-pipeline --release --example probe_drainage -- <seed> <lat> <lon> <half_km> <spacing_m>`
//! - Priority-Flood fill (box edges = outlets): share of area in closed depressions, depth.
//! - On macro river cells (drainage > 3000 km²): share where the final surface's steepest
//!   descent (D8) points within 60° of the macro flow direction.
use rg_pipeline::{generate, WorldParams};
use std::cmp::Ordering;
use std::collections::BinaryHeap;

fn main() {
    let a: Vec<f64> = std::env::args().skip(1).map(|s| s.parse().unwrap()).collect();
    let (seed, lat0, lon0, half, sp) = (a[0] as u64, a[1].to_radians(), a[2].to_radians(), a[3] * 1000.0, a[4]);
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
    let hy = w.hydrology.as_ref().expect("erosion on");
    let r = p.tectonics.radius_m;
    let n = (2.0 * half / sp) as usize;
    let dir_at = |i: usize, j: usize| {
        let (x, y) = (-half + i as f64 * sp, -half + j as f64 * sp);
        let (la, lo) = (lat0 + y / r, lon0 + x / (r * lat0.cos()));
        glam::DVec3::new(la.cos() * lo.cos(), la.cos() * lo.sin(), la.sin())
    };
    let h: Vec<f64> = (0..n * n).map(|k| s.height(dir_at(k % n, k / n), sp).0).collect();

    // Priority-Flood fill with the box edge as outlet.
    #[derive(PartialEq)]
    struct It(f64, usize);
    impl Eq for It {}
    impl Ord for It {
        fn cmp(&self, o: &Self) -> Ordering { o.0.total_cmp(&self.0).then(o.1.cmp(&self.1)) }
    }
    impl PartialOrd for It {
        fn partial_cmp(&self, o: &Self) -> Option<Ordering> { Some(self.cmp(o)) }
    }
    let mut filled = h.clone();
    let mut done = vec![false; n * n];
    let mut heap = BinaryHeap::new();
    for k in 0..n * n {
        let (i, j) = (k % n, k / n);
        if i == 0 || j == 0 || i == n - 1 || j == n - 1 || h[k] <= 0.0 {
            done[k] = true;
            heap.push(It(h[k], k));
        }
    }
    while let Some(It(hc, c)) = heap.pop() {
        let (ci, cj) = ((c % n) as i64, (c / n) as i64);
        for (di, dj) in [(-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1)] {
            let (ni, nj) = (ci + di, cj + dj);
            if ni < 0 || nj < 0 || ni >= n as i64 || nj >= n as i64 { continue; }
            let k = nj as usize * n + ni as usize;
            if done[k] { continue; }
            done[k] = true;
            filled[k] = filled[k].max(hc);
            heap.push(It(filled[k], k));
        }
    }
    let land = h.iter().filter(|&&x| x > 0.0).count().max(1);
    let depths: Vec<f64> = (0..n * n).map(|k| filled[k] - h[k]).filter(|&d| d > 0.5).collect();
    let mean_d = depths.iter().sum::<f64>() / depths.len().max(1) as f64;
    let max_d = depths.iter().cloned().fold(0.0, f64::max);
    println!("depressions: {:.2}% of land, mean depth {:.0} m, max {:.0} m (at {sp} m)", depths.len() as f64 / land as f64 * 100.0, mean_d, max_d);

    // Flow agreement on macro river cells.
    let (mut agree, mut total) = (0usize, 0usize);
    for j in 1..n - 1 {
        for i in 1..n - 1 {
            let k = j * n + i;
            let d = dir_at(i, j);
            if h[k] <= 0.0 || hy.drainage_km2.sample(d) < 3000.0 { continue; }
            let f = hy.flow_dir.nearest(d);
            let macro_dir = glam::DVec3::new(f[0] as f64, f[1] as f64, f[2] as f64);
            if macro_dir.length() < 0.5 { continue; }
            // Steepest descent on the filled surface (so flats in depressions still route).
            let mut best = (0.0, glam::DVec3::ZERO);
            for (di, dj) in [(-1i64, -1i64), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1)] {
                let (ni, nj) = ((i as i64 + di) as usize, (j as i64 + dj) as usize);
                let kk = nj * n + ni;
                let dist = sp * ((di * di + dj * dj) as f64).sqrt();
                let s = (filled[k] - filled[kk]) / dist;
                if s > best.0 { best = (s, (dir_at(ni, nj) - d).normalize()); }
            }
            if best.0 > 0.0 {
                total += 1;
                if best.1.dot(macro_dir) > 0.5 { agree += 1; }
            }
        }
    }
    println!("macro-river flow agreement: {:.0}% of {} river cells within 60°", agree as f64 / total.max(1) as f64 * 100.0, total);
}
