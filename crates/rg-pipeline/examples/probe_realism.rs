//! Realism scorecard for the **final** surface (macro + streamed detail), in several boxes:
//! `cargo run -p rg-pipeline --release --example probe_realism -- <seed> [half_km]`
//!
//! Env: `RG_RES` (tectonics res), `RG_EROSION_RES`, `RG_EROSION_STEPS`, plus the detail
//! overrides `TerrainSampler` reads (`RG_HILL_SLOPE`, `RG_GULLY_STRENGTH`, ...).
//! The generated world is cached in `target/probe-cache/`, so re-runs that only change detail
//! settings skip generation.
//!
//! Per box (three of the most rugged macro spots, one ordinary lowland):
//! - **depressions:** Priority-Flood fill on a 250 m grid, box edges = outlets. Share of land
//!   in closed depressions and their mean depth. Earth's fluvial landscapes: ~1% at this scale.
//! - **coherence / aspect:** structure-tensor coherence of the gradient in 1.5 km windows,
//!   and gradient-weighted aspect concentration, on a 10 km sub-box at 60 m. Fishbone,
//!   herringbone and grid-aligned patterns score high; isotropic terrain low.
//! - **slope:** 90 m central differences (SRTM-like), degrees.
#[path = "common/mod.rs"]
mod common;
#[path = "common/metrics.rs"]
mod metrics;
use common::world;
use rg_pipeline::World;
use std::cmp::Ordering;
use std::collections::BinaryHeap;

fn unit(lat: f64, lon: f64) -> glam::DVec3 {
    glam::DVec3::new(lat.cos() * lon.cos(), lat.cos() * lon.sin(), lat.sin())
}

/// Box centres: the three most rugged macro spots at least 1500 km apart, then a lowland.
fn pick_boxes(w: &World) -> Vec<(&'static str, f64, f64)> {
    let f = &w.tectonics.height_m;
    let res = f.res;
    let mut rugged: Vec<(f32, glam::DVec3)> = Vec::new();
    let mut low: Vec<glam::DVec3> = Vec::new();
    for (k, face, i, j, d) in f.iter_samples() {
        if i == 0 || j == 0 || i == res || j == res || f.data[k] <= 0.0 {
            continue;
        }
        let nb = [f.index(face, i + 1, j), f.index(face, i - 1, j), f.index(face, i, j + 1), f.index(face, i, j - 1)];
        let relief = nb.iter().map(|&n| (f.data[n] - f.data[k]).abs()).fold(0.0f32, f32::max);
        if f.data[k] > 500.0 {
            rugged.push((relief, d));
        } else if f.data[k] > 50.0 && f.data[k] < 300.0 && nb.iter().all(|&n| f.data[n] > 30.0) {
            low.push(d);
        }
    }
    rugged.sort_by(|a, b| b.0.total_cmp(&a.0));
    let r = w.params.tectonics.radius_m;
    let mut picked: Vec<glam::DVec3> = Vec::new();
    for (_, d) in rugged {
        if picked.iter().all(|p| p.angle_between(d) * r > 1.5e6) {
            picked.push(d);
        }
        if picked.len() == 3 {
            break;
        }
    }
    let ll = |d: glam::DVec3| (d.z.asin(), d.y.atan2(d.x));
    let names = ["mountain A", "mountain B", "mountain C"];
    let mut out: Vec<_> = picked.iter().zip(names).map(|(d, n)| (n, ll(*d).0, ll(*d).1)).collect();
    if let Some(d) = low.get(low.len() / 2) {
        out.push(("lowland", ll(*d).0, ll(*d).1));
    }
    out
}

struct Score {
    depress_pct: f64,
    depth_m: f64,
    slope: [f64; 4], // mean, p50, p90, p99
    macro_slope_mean: f64,
    /// Mean structure-tensor coherence of the gradient in 1.5 km windows (0 = isotropic,
    /// 1 = one orientation): fishbone/herringbone and grid-aligned patterns score high.
    coherence: f64,
    /// Gradient-weighted aspect concentration |Σ|g|e^{2iθ}| / Σ|g| over the 10 km sub-box.
    aspect: f64,
    /// 1D power-spectrum exponent and slope–area concavity (see `common/metrics.rs`).
    beta: f64,
    theta: f64,
}

fn score(w: &World, s: &rg_terrain::TerrainSampler, lat0: f64, lon0: f64, half: f64) -> Score {
    let r = w.params.tectonics.radius_m;
    let at = |x: f64, y: f64| unit(lat0 + y / r, lon0 + x / (r * lat0.cos()));
    let sp: f64 = std::env::var("RG_PROBE_SP").ok().and_then(|v| v.parse().ok()).unwrap_or(250.0);
    let cut: f64 = std::env::var("RG_PROBE_CUT").ok().and_then(|v| v.parse().ok()).unwrap_or(sp);
    let n = (2.0 * half / sp) as usize;
    let dir_at = |i: usize, j: usize| at(-half + i as f64 * sp, -half + j as f64 * sp);
    // RG_PROBE_SURFACE=macro scores the bare macro surface (bilinear committed heights).
    let mode = std::env::var("RG_PROBE_SURFACE").unwrap_or_default();
    let macro_only = mode == "macro";
    // RG_PROBE_SURFACE=gully: the old stateless detail (hills + gully filter), for comparison.
    let gully = {
        let mut g = s.clone();
        g.net = None;
        g
    };
    let s = if mode == "gully" { &gully } else { s };
    // RG_PROBE_SURFACE=net: the drained detail network, through cached 5 km sub-patches
    // (keyed by band-limit and block).
    let dnet = (mode == "net").then(|| {
        let hy = w.hydrology.as_ref().expect("erosion on");
        let mn = rg_terrain::network::MacroNet::new(
            std::sync::Arc::new(w.tectonics.height_m.clone()),
            std::sync::Arc::new(hy.flow_dir.clone()),
            std::sync::Arc::new(hy.drainage_km2.clone()),
        )
        .expect("hydrology at macro res");
        let mn = match &hy.uplift_mm_yr {
            Some(u) => mn.with_uplift(std::sync::Arc::new(u.clone())),
            None => mn,
        };
        rg_terrain::detail_net::DetailNet::new(mn, w.params.tectonics.radius_m, w.params.tectonics.seed).with_erodibility(w.params.erosion.k / hy.height_scale.unwrap_or(1.0))
    });
    let patches: std::cell::RefCell<std::collections::HashMap<(i64, i64, u64), rg_terrain::detail_net::Patch>> = Default::default();
    let rr = w.params.tectonics.radius_m;
    let surf = |d: glam::DVec3, wl: f64| {
        if macro_only {
            return s.macro_height(d);
        }
        let Some(dn) = &dnet else { return s.height(d, wl).0 };
        // Block of the local box frame (lat/lon offsets in m).
        let (la, lo) = (d.z.asin(), d.y.atan2(d.x));
        let (x, y) = ((lo - lon0) * rr * lat0.cos(), (la - lat0) * rr);
        let bs = 5000.0;
        let key = ((x / bs).floor() as i64, (y / bs).floor() as i64, wl.to_bits());
        let mut ps = patches.borrow_mut();
        // Bounded: fine patches are several MB each, and random slope samples touch hundreds.
        if ps.len() >= 24 && !ps.contains_key(&key) {
            ps.clear();
        }
        let patch = ps.entry(key).or_insert_with(|| {
            let (cx, cy) = ((key.0 as f64 + 0.5) * bs, (key.1 as f64 + 0.5) * bs);
            let (cla, clo) = (lat0 + cy / rr, lon0 + cx / (rr * lat0.cos()));
            let c = glam::DVec3::new(cla.cos() * clo.cos(), cla.cos() * clo.sin(), cla.sin());
            dn.patch(c, bs * 0.75, wl)
        });
        patch.heights(d)[0]
    };
    let h: Vec<f64> = (0..n * n).map(|k| surf(dir_at(k % n, k / n), cut)).collect();

    // Priority-Flood fill with the box edge (and the sea) as outlet.
    #[derive(PartialEq)]
    struct It(f64, usize);
    impl Eq for It {}
    impl Ord for It {
        fn cmp(&self, o: &Self) -> Ordering {
            o.0.total_cmp(&self.0).then(o.1.cmp(&self.1))
        }
    }
    impl PartialOrd for It {
        fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
            Some(self.cmp(o))
        }
    }
    const NB: [(i64, i64); 8] = [(-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1)];
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
        for (di, dj) in NB {
            let (ni, nj) = (ci + di, cj + dj);
            if ni < 0 || nj < 0 || ni >= n as i64 || nj >= n as i64 {
                continue;
            }
            let k = nj as usize * n + ni as usize;
            if done[k] {
                continue;
            }
            done[k] = true;
            filled[k] = filled[k].max(hc);
            heap.push(It(filled[k], k));
        }
    }
    let land = h.iter().filter(|&&x| x > 0.0).count().max(1);
    let depths: Vec<f64> = (0..n * n).map(|k| filled[k] - h[k]).filter(|&d| d > 0.5).collect();
    let depth_m = depths.iter().sum::<f64>() / depths.len().max(1) as f64;

    // 90 m slopes at 2000 pseudo-random points (only on land).
    let mut rng = 12345u64;
    let mut next = || {
        rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (rng >> 11) as f64 / (1u64 << 53) as f64
    };
    let d = 45.0;
    let (mut deg, mut mdeg) = (Vec::new(), Vec::new());
    for _ in 0..2000 {
        let (x, y) = ((next() * 2.0 - 1.0) * half, (next() * 2.0 - 1.0) * half);
        let hh = |x: f64, y: f64| surf(at(x, y), 40.0);
        if hh(x, y) <= 0.0 {
            continue;
        }
        let gx = (hh(x + d, y) - hh(x - d, y)) / (2.0 * d);
        let gy = (hh(x, y + d) - hh(x, y - d)) / (2.0 * d);
        deg.push((gx * gx + gy * gy).sqrt().atan().to_degrees());
        let mh = |x: f64, y: f64| s.macro_height(at(x, y));
        let gx = (mh(x + d, y) - mh(x - d, y)) / (2.0 * d);
        let gy = (mh(x, y + d) - mh(x, y - d)) / (2.0 * d);
        mdeg.push((gx * gx + gy * gy).sqrt().atan().to_degrees());
    }
    deg.sort_by(f64::total_cmp);
    let q = |f: f64| deg.get(((deg.len().max(1) - 1) as f64 * f) as usize).copied().unwrap_or(0.0);
    let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len().max(1) as f64;
    // Orientation: a 10 km sub-box at 60 m, structure tensor in 25×25-sample windows.
    let (m, sp2) = (167usize, 60.0);
    let g: Vec<f64> = (0..m * m).map(|k| surf(at(-5000.0 + (k % m) as f64 * sp2, -5000.0 + (k / m) as f64 * sp2), sp2)).collect();
    // Shared with earth_reference (identical code on real DEMs).
    let (coherence, aspect) = metrics::orientation(&g, m, sp2);
    let beta = metrics::spectrum_beta(&g, m, sp2);
    let (theta, _) = metrics::slope_area_theta(&g, m, sp2);
    Score {
        coherence,
        aspect,
        beta,
        theta,
        depress_pct: depths.len() as f64 / land as f64 * 100.0,
        depth_m,
        slope: [mean(&deg), q(0.5), q(0.9), q(0.99)],
        macro_slope_mean: mean(&mdeg),
    }
}

/// Channel straightness on the macro network: for channel samples (> 1000 km²), the share of
/// downstream steps that keep exactly the same flow direction, and the longest straight run.
/// Plain D8 on a smooth tilted surface gives long parallel straight runs.
fn straightness(w: &World) -> (f64, usize, usize) {
    let hm = &w.tectonics.height_m;
    let hy = w.hydrology.as_ref().expect("erosion on");
    let (a, fd) = (&hy.drainage_km2, &hy.flow_dir);
    let ang = a.spacing_rad();
    let v = |f: [f32; 3]| glam::DVec3::new(f[0] as f64, f[1] as f64, f[2] as f64);
    let next = |d: glam::DVec3| -> Option<glam::DVec3> {
        let f = v(fd.nearest(d));
        if f.length() < 0.5 || a.nearest(d) < 1000.0 {
            return None;
        }
        // Step to the downstream sample: diagonal steps are ~1.41 cells.
        let step = (d + f * ang * 1.2).normalize();
        Some(step)
    };
    let (mut same, mut total, mut longest) = (0usize, 0usize, 0usize);
    for (k, _, _, _, d) in a.iter_samples() {
        if a.data[k] < 1000.0 || v(fd.data[k]).length() < 0.5 {
            continue;
        }
        let f0 = v(fd.data[k]);
        if let Some(n) = next(d) {
            let f1 = v(fd.nearest(n));
            if f1.length() > 0.5 {
                total += 1;
                let grad = (hm.sample_drained(d) - hm.sample_drained(n)) / (ang * 1.2 * w.params.tectonics.radius_m);
                let flat = grad < 2e-4;
                FLAT.with(|c| { let mut c = c.borrow_mut(); if flat { c.0 += 1 } else { c.2 += 1 } });
                if f0.dot(f1) > 0.9995 {
                    FLAT.with(|c| { let mut c = c.borrow_mut(); if flat { c.1 += 1 } else { c.3 += 1 } });
                    same += 1;
                    // Walk the run (only from its upstream end would be cheaper; fine here).
                    let (mut run, mut p, mut fp) = (1usize, n, f1);
                    while run < 200 {
                        let Some(q) = next(p) else { break };
                        let fq = v(fd.nearest(q));
                        if fq.length() < 0.5 || fp.dot(fq) <= 0.9995 {
                            break;
                        }
                        run += 1;
                        p = q;
                        fp = fq;
                    }
                    longest = longest.max(run);
                }
            }
        }
    }
    (same as f64 / total.max(1) as f64 * 100.0, longest, total)
}

thread_local! {
    /// (flat steps, flat straight, sloped steps, sloped straight)
    static FLAT: std::cell::RefCell<(usize, usize, usize, usize)> = const { std::cell::RefCell::new((0, 0, 0, 0)) };
}

fn main() {
    let seed: u64 = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(3);
    let half_km: f64 = std::env::args().nth(2).and_then(|s| s.parse().ok()).unwrap_or(40.0);
    let w = world(seed);
    let s = rg_pipeline::terrain_sampler(&w);
    let (pct, longest, steps) = straightness(&w);
    println!("channels: {pct:.0}% of {steps} steps continue straight, longest straight run {longest} cells");
    let (f, fs, sl, ss) = FLAT.with(|c| *c.borrow());
    println!("  flats (< 0.2 m/km): {:.0}% of steps, {:.0}% straight; slopes: {:.0}% straight", f as f64 / steps.max(1) as f64 * 100.0, fs as f64 / f.max(1) as f64 * 100.0, ss as f64 / sl.max(1) as f64 * 100.0);
    let boxes = pick_boxes(&w);
    let scores: Vec<Score> = std::thread::scope(|sc| {
        let (w, s) = (&w, &s);
        let hs: Vec<_> = boxes.iter().map(|&(_, la, lo)| sc.spawn(move || score(w, s, la, lo, half_km * 1000.0))).collect();
        hs.into_iter().map(|h| h.join().unwrap()).collect()
    });
    println!("box          lat     lon    depress  depth   coh   asp    β     θ    slope° mean/p50/p90/p99   macro°");
    for ((name, la, lo), sc) in boxes.iter().zip(&scores) {
        println!(
            "{name:<11} {:6.1} {:7.1}   {:5.2}%  {:4.0} m  {:.2}  {:.2}  {:4.2}  {:4.2}   {:4.1} {:4.1} {:4.1} {:4.1}      {:4.1}",
            la.to_degrees(),
            lo.to_degrees(),
            sc.depress_pct,
            sc.depth_m,
            sc.coherence,
            sc.aspect,
            sc.beta,
            sc.theta,
            sc.slope[0],
            sc.slope[1],
            sc.slope[2],
            sc.slope[3],
            sc.macro_slope_mean
        );
    }
}
