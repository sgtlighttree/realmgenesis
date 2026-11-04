//! Regional detail spike (s02, plan B milestone M0; `docs/design/regional-detail.md`): one
//! window of the globe re-eroded coarse-to-fine at high resolution, written as hillshade
//! PNGs plus stats, so the look can be judged before any schema or UI work.
//!
//! `cargo run -p rg-pipeline --release --example regional_spike -- <seed> <lat> <lon> <half_km> <out_prefix>`
//!
//! Writes `<out_prefix>.png` (the regional result, with rivers) and `<out_prefix>_macro.png`
//! (the same window from the smooth macro surface, for comparison).
//!
//! Env (all optional):
//! - `RG_REG_LEVELS` (3): cascade levels. The finest has 2^(10 + levels − 2) + 1 nodes a
//!   side, so 3 → 2049² (2^11 + 1) and 2 → 1025².
//! - `RG_REG_STEPS` ("16,8,6"): implicit steps per level, coarse to fine.
//! - `RG_REG_DT_MYR` ("1.0,0.5,0.25"): step length per level, Myr.
//! - `RG_REG_REPOSE_DEG` (33): threshold hillslope angle; 0 disables the talus cap.
//! - `RG_REG_NOISE` (0.05): new-octave noise amplitude per level, as a fraction of the node
//!   spacing (land only, keyed to position).
//! - `RG_REG_SHADE_EX` (2.0): hillshade vertical exaggeration.
//! - `RG_REG_UPLIFT` (1.0): multiplier on the macro's uplift. 0 = detail erosion only
//!   (carve into the upsampled base, as Gaea/WM/Terragen's detail stages do).
//! - `RG_REG_ANCHOR` (0): 1 = after each level, replace everything coarser than a macro
//!   cell with the macro's own shape (h += blur(macro − h)), keeping only the window's
//!   fine detail, so the macro stays the authority (no drift, no edge step).
//! - `RG_REG_EDGE_BAND` (0.03125): with anchoring, fade the detail to the macro across this
//!   fraction of the window at each edge (smoothstep), so the window has no seam.
//! - the world's usual `RG_RES` etc. (cached in `target/probe-cache/`).
#[path = "common/mod.rs"]
mod common;
#[path = "common/metrics.rs"]
mod metrics;
use glam::DVec3;
use rg_geo::{CubeFace, CubeField};
use rg_tectonics::grid::Grid;
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::time::Instant;

/// A square window of one cube face, in that face's (u, v).
struct Window {
    face: CubeFace,
    u0: f64,
    v0: f64,
    /// Half-width in u and v.
    h_uv: f64,
}

impl Window {
    /// Unit direction of node (i, j) of an n×n raster spanning the window, corners included.
    fn dir(&self, n: usize, i: usize, j: usize) -> DVec3 {
        let t = |k: usize| -1.0 + 2.0 * k as f64 / (n - 1) as f64;
        self.face.uv_to_dir(self.u0 + self.h_uv * t(i), self.v0 + self.h_uv * t(j))
    }

    /// Position-keyed lattice coordinates of node (i, j) at this level: the face's (u, v)
    /// quantised at the node spacing, so the same place hashes the same way in any window.
    fn lattice(&self, n: usize, i: usize, j: usize) -> (i64, i64) {
        let du = 2.0 * self.h_uv / (n - 1) as f64;
        let u = self.u0 - self.h_uv + i as f64 * du;
        let v = self.v0 - self.h_uv + j as f64 * du;
        ((u / du).round() as i64, (v / du).round() as i64)
    }

    /// (i, j) of the node nearest to `d` on an n×n raster, if `d` is inside the window.
    fn node_of(&self, n: usize, d: DVec3) -> Option<(usize, usize)> {
        let (face, u, v) = CubeFace::from_dir(d);
        if face != self.face || (u - self.u0).abs() > self.h_uv || (v - self.v0).abs() > self.h_uv {
            return None;
        }
        let s = (n - 1) as f64 / (2.0 * self.h_uv);
        let i = ((u - self.u0 + self.h_uv) * s).round().clamp(0.0, (n - 1) as f64) as usize;
        let j = ((v - self.v0 + self.h_uv) * s).round().clamp(0.0, (n - 1) as f64) as usize;
        Some((i, j))
    }
}

/// An 8-connected raster `Grid` over the window. Neighbour lists are sorted (row-major
/// order), which `diagonal_mask` relies on. Only `dir`, `n`, `weight` and the adjacency
/// are used by the solver; the cube-sample maps stay empty.
fn regional_grid(w: &Window, n: usize) -> Grid {
    let dir: Vec<DVec3> = (0..n * n).map(|k| w.dir(n, k % n, k / n)).collect();
    let weight: Vec<f32> = (0..n * n)
        .map(|k| {
            let (i, j) = (k % n, k / n);
            let (i0, i1) = (i.saturating_sub(1), (i + 1).min(n - 1));
            let (j0, j1) = (j.saturating_sub(1), (j + 1).min(n - 1));
            let du = (dir[j * n + i1] - dir[j * n + i0]) / (i1 - i0) as f64;
            let dv = (dir[j1 * n + i] - dir[j0 * n + i]) / (j1 - j0) as f64;
            (du.cross(dv).length() / (4.0 * std::f64::consts::PI)) as f32
        })
        .collect();
    let mut adj_start = Vec::with_capacity(n * n + 1);
    let mut adj = Vec::with_capacity(n * n * 8);
    let (mut min_edge, mut max_edge) = (f32::MAX, 0.0f32);
    for k in 0..n * n {
        adj_start.push(adj.len() as u32);
        let (i, j) = ((k % n) as i64, (k / n) as i64);
        for dj in -1..=1i64 {
            for di in -1..=1i64 {
                let (x, y) = (i + di, j + dj);
                if (di, dj) == (0, 0) || x < 0 || y < 0 || x >= n as i64 || y >= n as i64 {
                    continue;
                }
                let b = y as usize * n + x as usize;
                adj.push(b as u32);
                let e = dir[k].angle_between(dir[b]) as f32;
                min_edge = min_edge.min(e);
                max_edge = max_edge.max(e);
            }
        }
    }
    adj_start.push(adj.len() as u32);
    Grid { res: (n - 1) as u32, n: n * n, dir, weight, adj_start, adj, min_edge, max_edge, sample_node: vec![], node_sample: vec![] }
}

/// splitmix64 of (seed, key) to [0, 1).
fn hash01(seed: u64, key: u64) -> f64 {
    let mut z = seed ^ key.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^= z >> 31;
    (z >> 11) as f64 / (1u64 << 53) as f64
}

/// Min-heap entry by height.
struct Low(f64, usize);
impl PartialEq for Low {
    fn eq(&self, o: &Self) -> bool {
        self.cmp(o) == Ordering::Equal
    }
}
impl Eq for Low {}
impl PartialOrd for Low {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Low {
    fn cmp(&self, o: &Self) -> Ordering {
        o.0.total_cmp(&self.0).then(o.1.cmp(&self.1))
    }
}

/// Threshold hillslopes: no node stands more than `s_max · distance` above a neighbour
/// (landslides; the material is assumed carried off by rivers, Burbank et al. 1996 style).
/// Processed from the lowest node up (Dijkstra), so one pass is exact. Fixed nodes don't
/// move; sea neighbours count as 0 m so the coast never moves inland or out.
fn talus(grid: &Grid, h: &mut [f64], fixed: &[bool], radius_m: f64, s_max: f64) -> usize {
    let mut heap: BinaryHeap<Low> = (0..grid.n).map(|i| Low(h[i], i)).collect();
    let mut lowered = 0;
    while let Some(Low(hi, i)) = heap.pop() {
        if hi != h[i] {
            continue; // stale entry
        }
        let base = hi.max(0.0);
        for &b in grid.neighbors(i) {
            let b = b as usize;
            if fixed[b] {
                continue;
            }
            let cap = base + s_max * grid.dir[i].angle_between(grid.dir[b]) * radius_m;
            if h[b] > cap {
                h[b] = cap;
                heap.push(Low(cap, b));
                lowered += 1;
            }
        }
    }
    lowered
}

fn list(var: &str, default: &str) -> Vec<f64> {
    std::env::var(var).unwrap_or_else(|_| default.into()).split(',').map(|s| s.trim().parse().expect(var)).collect()
}

fn main() {
    let t_all = Instant::now();
    let a: Vec<String> = std::env::args().skip(1).collect();
    let seed: u64 = a[0].parse().unwrap();
    let (lat, lon) = (a[1].parse::<f64>().unwrap().to_radians(), a[2].parse::<f64>().unwrap().to_radians());
    let half_m = a[3].parse::<f64>().unwrap() * 1000.0;
    let prefix = &a[4];
    let levels: usize = common::env("RG_REG_LEVELS").unwrap_or(3);
    let steps = list("RG_REG_STEPS", "16,8,6");
    let dts = list("RG_REG_DT_MYR", "1.0,0.5,0.25");
    let repose: f64 = common::env("RG_REG_REPOSE_DEG").unwrap_or(33.0);
    let noise: f64 = common::env("RG_REG_NOISE").unwrap_or(0.05);
    let shade_ex: f64 = common::env("RG_REG_SHADE_EX").unwrap_or(2.0);
    let uplift_mult: f64 = common::env("RG_REG_UPLIFT").unwrap_or(1.0);
    let anchor = common::env::<u32>("RG_REG_ANCHOR").unwrap_or(0) == 1;
    let edge_band: f64 = common::env("RG_REG_EDGE_BAND").unwrap_or(0.03125);
    assert!(steps.len() >= levels && dts.len() >= levels, "RG_REG_STEPS / RG_REG_DT_MYR need one entry per level");

    let w = common::world(seed);
    let r = w.params.tectonics.radius_m;
    let hm = &w.tectonics.height_m;
    let hy = w.hydrology.as_ref().expect("the world needs erosion (hydrology)");
    let k_eff = w.params.erosion.k / hy.height_scale.unwrap_or(1.0);

    // The window: centred on (lat, lon), sized so its half-width is ~half_km on the ground.
    let centre = DVec3::new(lat.cos() * lon.cos(), lat.cos() * lon.sin(), lat.sin());
    let (face, u0, v0) = CubeFace::from_dir(centre);
    let mut h_uv = half_m / r;
    for _ in 0..3 {
        let span = face.uv_to_dir(u0 - h_uv, v0).angle_between(face.uv_to_dir(u0 + h_uv, v0)) * r / 2.0;
        h_uv *= half_m / span;
    }
    assert!(u0.abs() + h_uv < 1.0 && v0.abs() + h_uv < 1.0, "window crosses a cube-face edge; move it or shrink it (v1 limit)");
    let win = Window { face, u0, v0, h_uv };

    // Relative runoff (land mean = 1), as the macro erosion used, sampled per node.
    let wb = rg_pipeline::water_balance(&w.params.planet, &w.climate, hm);
    let mean_ro = wb.mean_land_runoff_mm.max(1.0);

    let n_final = (1usize << (levels + 8)) + 1;
    eprintln!(
        "window: face {:?}, {:.0} km wide, levels {levels} → {n_final}² nodes ({:.0} m spacing); K_eff {k_eff:.3e}",
        face,
        2.0 * half_m / 1000.0,
        2.0 * half_m / (n_final - 1) as f64
    );

    let mut h_prev: Vec<f64> = Vec::new();
    let mut n_prev = 0usize;
    let mut last_area: Vec<f64> = Vec::new();
    let mut grid_final: Option<Grid> = None;
    for level in 0..levels {
        let t_level = Instant::now();
        let n = (1usize << (level + 9)) + 1;
        let grid = regional_grid(&win, n);
        let sp = 2.0 * half_m / (n - 1) as f64;
        let macro_h: Vec<f64> = grid.dir.iter().map(|&d| hm.sample_drained(d)).collect();
        // Sea is fixed from the macro, and so is the edge ring (base level, outlets).
        let edge = |k: usize| {
            let (i, j) = (k % n, k / n);
            i == 0 || j == 0 || i == n - 1 || j == n - 1
        };
        let fixed: Vec<bool> = (0..n * n).map(|k| macro_h[k] <= 0.0 || edge(k)).collect();

        // Initial surface: the macro (level 0) or the previous level upsampled ×2 with the
        // drained interpolant (a non-drained start makes thousands of pits).
        let mut h: Vec<f64> = if level == 0 {
            macro_h.clone()
        } else {
            let c = |i: usize, j: usize| h_prev[j.min(n_prev - 1) * n_prev + i.min(n_prev - 1)];
            (0..n * n)
                .map(|k| {
                    let (i, j) = (k % n, k / n);
                    let (ci, cj) = (i / 2, j / 2);
                    let (fx, fy) = ((i % 2) as f64 * 0.5, (j % 2) as f64 * 0.5);
                    rg_geo::field::drained_interp(c(ci, cj), c(ci + 1, cj), c(ci, cj + 1), c(ci + 1, cj + 1), fx, fy)
                })
                .collect()
        };
        // New octave: white noise keyed to position, land only; land stays land.
        for k in 0..n * n {
            if fixed[k] {
                h[k] = macro_h[k];
                continue;
            }
            let (li, lj) = win.lattice(n, k % n, k / n);
            let key = (face as u64) << 60 ^ ((li as u64) << 30) ^ (lj as u64) ^ ((level as u64) << 56);
            h[k] += (hash01(seed ^ 0x5EED_0FDE_7A11, key) - 0.5) * 2.0 * noise * sp;
            h[k] = h[k].max(1.0);
        }

        let uplift: Vec<f64> = grid
            .dir
            .iter()
            .map(|&d| hy.uplift_mm_yr.as_ref().map_or(0.0, |u| u.sample(d) as f64) * 1e-3 * uplift_mult)
            .collect();
        let k_field = vec![k_eff; n * n];
        let cell_area: Vec<f64> = grid.weight.iter().map(|&wt| wt as f64 * 4.0 * std::f64::consts::PI * r * r).collect();
        let mut runoff: Vec<f64> = grid.dir.iter().map(|&d| (wb.runoff_mm.sample(d) as f64 / mean_ro).clamp(0.02, 6.0)).collect();
        // Inflow: where a macro river crosses into the window, inject its drainage area at
        // the first interior node (once per macro edge, not per regional edge node).
        let ang = hy.drainage_km2.spacing_rad();
        let mut injected = 0.0;
        let mut inflows = 0;
        for (k, _, _, _, d) in hy.drainage_km2.iter_samples() {
            let area_km2 = hy.drainage_km2.data[k] as f64;
            if area_km2 < 10.0 || win.node_of(n, d).is_some() || d.angle_between(centre) * r > 2.0 * half_m {
                continue;
            }
            let f = hy.flow_dir.data[k];
            let f = DVec3::new(f[0] as f64, f[1] as f64, f[2] as f64);
            if f.length() < 0.5 {
                continue;
            }
            let to = (d + f * ang * 1.2).normalize();
            if let Some((i, j)) = win.node_of(n, to) {
                let (i, j) = (i.clamp(1, n - 2), j.clamp(1, n - 2));
                let node = j * n + i;
                if !fixed[node] {
                    runoff[node] += area_km2 * 1e6 / cell_area[node];
                    injected += area_km2;
                    inflows += 1;
                }
            }
        }

        let mut ep = w.params.erosion.clone();
        ep.steps = 1;
        ep.duration_myr = dts[level];
        let mut lowered = 0;
        for s in 0..steps[level] as usize {
            let solved = rg_erosion::solve_full_with_progress(
                &grid, r, &mut h, &fixed, &uplift, &k_field, Some(&runoff), None, None, &ep, &|_, _| true,
            )
            .expect("no cancel callback");
            if repose > 0.0 {
                lowered += talus(&grid, &mut h, &fixed, r, repose.to_radians().tan());
            }
            if s + 1 == steps[level] as usize {
                last_area = solved.net.area;
            }
        }
        if anchor {
            // Low-pass anchoring: blur the difference to the macro at about one macro cell
            // (three box passes ≈ Gaussian) and add it back, on the free nodes only.
            let radius = ((hm.spacing_rad() * r / sp) * 0.5).round().max(1.0) as usize;
            let mut diff: Vec<f64> = (0..n * n).map(|k| if fixed[k] { 0.0 } else { macro_h[k] - h[k] }).collect();
            for _ in 0..3 {
                diff = box_blur(&diff, n, radius);
            }
            // Edge taper: detail fades to the macro near the window's edges (no seam).
            let band = (edge_band * (n - 1) as f64).max(1.0);
            for k in 0..n * n {
                if !fixed[k] {
                    let (i, j) = (k % n, k / n);
                    let e = i.min(j).min(n - 1 - i).min(n - 1 - j) as f64 / band;
                    let t = if e >= 1.0 { 1.0 } else { e * e * (3.0 - 2.0 * e) };
                    let anchored = h[k] + diff[k];
                    h[k] = (macro_h[k] + (anchored - macro_h[k]) * t).max(1.0);
                }
            }
            if repose > 0.0 {
                lowered += talus(&grid, &mut h, &fixed, r, repose.to_radians().tan());
            }
        }
        let land: Vec<usize> = (0..n * n).filter(|&k| macro_h[k] > 0.0).collect();
        let drift = land.iter().map(|&k| h[k] - macro_h[k]).sum::<f64>() / land.len().max(1) as f64;
        eprintln!(
            "level {level}: {n}² ({sp:.0} m), {} steps × {} Myr, inflows {inflows} ({injected:.0} km²), talus lowered {lowered}, mean drift vs macro {drift:+.0} m, {:.1} s",
            steps[level],
            dts[level],
            t_level.elapsed().as_secs_f64()
        );
        h_prev = h;
        n_prev = n;
        if level + 1 == levels {
            grid_final = Some(grid);
        }
    }

    // ---- outputs and stats on the final level
    let n = n_prev;
    let sp = 2.0 * half_m / (n - 1) as f64;
    let grid = grid_final.unwrap();
    let macro_h: Vec<f64> = grid.dir.iter().map(|&d| hm.sample_drained(d)).collect();
    let smooth_h: Vec<f64> = grid.dir.iter().map(|&d| smooth_macro(hm, d)).collect();
    let h = &h_prev;
    // Slopes over land, at the raster spacing (degrees).
    let mut deg: Vec<f64> = Vec::new();
    for j in 1..n - 1 {
        for i in 1..n - 1 {
            let k = j * n + i;
            if h[k] <= 0.0 {
                continue;
            }
            let gx = (h[k + 1] - h[k - 1]) / (2.0 * sp);
            let gy = (h[k + n] - h[k - n]) / (2.0 * sp);
            deg.push((gx * gx + gy * gy).sqrt().atan().to_degrees());
        }
    }
    deg.sort_by(f64::total_cmp);
    let q = |f: f64| deg.get(((deg.len().max(1) - 1) as f64 * f) as usize).copied().unwrap_or(0.0);
    let (dep_pct, dep_mean) = metrics::depressions(h, n);
    // Edge band: nodes 2–8 in from the ring, against the macro.
    let band: Vec<f64> = (0..n * n)
        .filter(|&k| {
            let (i, j) = (k % n, k / n);
            let e = i.min(j).min(n - 1 - i).min(n - 1 - j);
            (2..=8).contains(&e) && macro_h[k] > 0.0
        })
        .map(|k| (h[k] - macro_h[k]).abs())
        .collect();
    let max_h = h.iter().cloned().fold(f64::MIN, f64::max);
    eprintln!(
        "final: land slopes mean {:.1}°, p50 {:.1}, p90 {:.1}, p99 {:.1} (at {sp:.0} m); depressions {dep_pct:.2}% (mean {dep_mean:.1} m); edge band |Δ| mean {:.0} m; max {max_h:.0} m (macro max {:.0}); total {:.1} s",
        deg.iter().sum::<f64>() / deg.len().max(1) as f64,
        q(0.5),
        q(0.9),
        q(0.99),
        band.iter().sum::<f64>() / band.len().max(1) as f64,
        macro_h.iter().cloned().fold(f64::MIN, f64::max),
        t_all.elapsed().as_secs_f64()
    );

    // Rivers: drainage above ~50 km² (wet area), drawn as blue pixels.
    let river_min_m2 = 50.0e6;
    write_png(&format!("{prefix}.png"), n, &hillshade(h, n, sp, shade_ex, Some((&last_area, river_min_m2))));
    write_png(&format!("{prefix}_macro.png"), n, &hillshade(&smooth_h, n, sp, shade_ex, None));
    eprintln!("wrote {prefix}.png and {prefix}_macro.png");
}

/// Separable box blur of an n×n raster with half-width `rad` (edges clamped).
fn box_blur(a: &[f64], n: usize, rad: usize) -> Vec<f64> {
    let pass = |src: &[f64], horizontal: bool| -> Vec<f64> {
        let mut out = vec![0.0; n * n];
        for line in 0..n {
            let at = |t: usize| if horizontal { src[line * n + t] } else { src[t * n + line] };
            // Running sum over the window [t − rad, t + rad], clamped at the ends.
            let mut sum: f64 = (0..=rad).map(|t| at(t.min(n - 1))).sum::<f64>() + at(0) * rad as f64;
            for t in 0..n {
                let v = sum / (2 * rad + 1) as f64;
                if horizontal { out[line * n + t] = v } else { out[t * n + line] = v }
                sum += at((t + rad + 1).min(n - 1)) - at(t.saturating_sub(rad));
            }
        }
        out
    };
    pass(&pass(a, true), false)
}

/// The macro-only display surface (bicubic, coast kept on the drained sign), as the app
/// draws it with streamed detail off.
fn smooth_macro(hm: &CubeField<f32>, d: DVec3) -> f64 {
    let (drained, smooth) = (hm.sample_drained(d), hm.sample_smooth(d));
    if (smooth > 0.0) == (drained > 0.0) { smooth } else { drained }
}

/// RGB hillshade of an n×n raster (row 0 = south, so rows are flipped to put north up),
/// coloured by height like `river_map`, optionally with rivers.
fn hillshade(h: &[f64], n: usize, sp: f64, ex: f64, rivers: Option<(&[f64], f64)>) -> Vec<u8> {
    let mut img = vec![0u8; n * n * 3];
    let light = DVec3::new(-0.5, 0.5, 0.7).normalize();
    for j in 0..n {
        for i in 0..n {
            let k = j * n + i;
            let (il, ir) = (i.saturating_sub(1), (i + 1).min(n - 1));
            let (jd, ju) = (j.saturating_sub(1), (j + 1).min(n - 1));
            let gx = (h[j * n + ir] - h[j * n + il]) / ((ir - il) as f64 * sp);
            let gy = (h[ju * n + i] - h[jd * n + i]) / ((ju - jd) as f64 * sp);
            let nrm = DVec3::new(-gx * ex, -gy * ex, 1.0).normalize();
            let shade = nrm.dot(light).max(0.0);
            let c = if h[k] <= 0.0 {
                [20.0, 40.0, 80.0]
            } else {
                let t = (h[k] / 4000.0).min(1.0);
                [90.0 + 120.0 * t, 110.0 + 80.0 * t, 70.0 + 110.0 * t]
            };
            let row = n - 1 - j;
            let o = (row * n + i) * 3;
            let river = rivers.is_some_and(|(a, min)| h[k] > 0.0 && a[k] > min);
            for ch in 0..3 {
                img[o + ch] = if river { [40, 100, 210][ch] } else { (c[ch] * (0.35 + 0.85 * shade)).min(255.0) as u8 };
            }
        }
    }
    img
}

fn write_png(path: &str, n: usize, img: &[u8]) {
    let file = std::fs::File::create(path).unwrap();
    let mut enc = png::Encoder::new(std::io::BufWriter::new(file), n as u32, n as u32);
    enc.set_color(png::ColorType::Rgb);
    enc.write_header().unwrap().write_image_data(img).unwrap();
}
