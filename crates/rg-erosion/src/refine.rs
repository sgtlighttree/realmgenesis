//! Refinement: detail a fine grid from a coarser, already-eroded **canonical** world
//! (s04, `docs/design/canonical-refine.md`; Matt: "Ultra basically upscales Normal").
//!
//! Erosion on its own is resolution-dependent: every cell is a channel, tiny hollows held
//! as lakes multiply on finer grids, and the global height calibration spreads both to
//! every region. So the large-scale shape comes from one canonical solve (Normal), and a
//! finer world only adds what Normal can't hold:
//!
//! 1. Start: `h = up(H_canon) + s·D` on land, where `up` is the smooth (Catmull-Rom)
//!    upsampling the display uses, `D = T_fine − up(T_canon)` is the sub-cell tectonic
//!    detail, and `s` is the canonical height factor (so `D` is in committed metres).
//!    Ocean nodes keep the fine tectonic sea floor.
//! 2. Lakes: canonical lake levels are upsampled and held as base levels. The lake model
//!    doesn't run here (on a fine grid it multiplies tiny lakes). Every other pit is filled.
//! 3. Incision: a few implicit stream-power steps with **no uplift**. Then everything added
//!    to the upsampled canonical surface (detail + incision) is made **zero-mean at the
//!    canonical cell scale**: its average over each canonical cell (smoothly upsampled) is
//!    subtracted. Without uplift, a few Myr of stream power also lowers everything by
//!    hundreds of metres (s04: world median 272 → 12 m), and `D` alone doesn't average to
//!    zero over a cell; the correction keeps only sub-cell relief, so the canonical shape
//!    can't move. There is no
//!    recalibration. With n = 1 and U = 0 the incision is linear in height, so running on
//!    calibrated (committed) heights equals running in solver units and scaling.
//!
//! Deterministic like the solver (no hash maps, no clock; threads only change speed).

use crate::solver::{implicit_step_based, lithology, priority_flood_fill_to, Column, ErosionParams, Network, RiverSegment};
use rg_geo::CubeField;
use rg_tectonics::grid::Grid;
use rg_tectonics::util::par_map;

/// The canonical world a refinement starts from (all at the canonical resolution, except
/// `fine_t`).
pub struct RefineInput<'a> {
    /// Eroded, calibrated canonical heights (committed metres).
    pub canon_h: &'a CubeField<f32>,
    /// Canonical tectonic heights, before erosion.
    pub canon_t: &'a CubeField<f32>,
    /// Canonical standing-water depth, m.
    pub canon_lake_m: &'a CubeField<f32>,
    /// 1 where the canonical lake is closed (terminal), else 0.
    pub canon_closed: &'a CubeField<u8>,
    /// Fine tectonic heights: their resolution is the output resolution.
    pub fine_t: &'a CubeField<f32>,
    /// The canonical solved → committed height factor (1 if none).
    pub height_scale: f64,
    /// Relative runoff (1 = land mean), any resolution, as the canonical solve used it.
    pub runoff: Option<&'a CubeField<f32>>,
}

/// The refined fields, all at the fine resolution.
pub struct Refined {
    pub height_m: CubeField<f32>,
    pub drainage_km2: CubeField<f32>,
    pub flow_dir: CubeField<[f32; 3]>,
    pub lake_depth_m: CubeField<f32>,
    pub rivers: Vec<RiverSegment>,
}

/// Below this canonical water depth (m) a node doesn't count as lake (interpolation fringe).
const LAKE_MIN_DEPTH_M: f64 = 0.5;

/// Lakes smaller than this (km²) are not drawn (`rg-pipeline` layers: one- to three-cell
/// solver lakes read as speckle, Matt s02) and are refined as **playas**: the bed is raised
/// flat to the water level, which keeps the basin closed without leaving a square pit.
pub const MIN_VISIBLE_LAKE_KM2: f64 = 3_900.0;

pub fn refine(inp: &RefineInput, radius_m: f64, p: &ErosionParams) -> Refined {
    let res = inp.fine_t.res;
    let grid = Grid::new(res);
    let n = grid.n;
    let area_m2: Vec<f64> = grid.weight.iter().map(|w| *w as f64 * 4.0 * std::f64::consts::PI * radius_m * radius_m).collect();
    let s = inp.height_scale;

    // 1. Start surface.
    let t_fine: Vec<f64> = (0..n).map(|i| inp.fine_t.data[grid.node_sample[i] as usize] as f64).collect();
    let ocean: Vec<bool> = t_fine.iter().map(|&t| t <= 0.0).collect();
    let canon_up: Vec<f64> = par_map(n, |i| inp.canon_h.sample_smooth(grid.dir[i]));
    let mut h: Vec<f64> = par_map(n, |i| {
        if ocean[i] {
            return t_fine[i];
        }
        let detail = t_fine[i] - inp.canon_t.sample_smooth(grid.dir[i]);
        (canon_up[i] + s * detail).max(1.0)
    });

    // 2. Canonical lakes: each keeps its flat water level and floods the fine terrain below
    //    it, inside the lake's canonical cells, so shores follow the fine relief (s04:
    //    an upsampled depth field drew square lakes).
    let cg = Grid::new(inp.canon_h.res);
    let canon_level: Vec<f64> = (0..cg.n)
        .map(|c| {
            let k = cg.node_sample[c] as usize;
            let depth = inp.canon_lake_m.data[k] as f64;
            if depth > LAKE_MIN_DEPTH_M { inp.canon_h.data[k] as f64 + depth } else { f64::NAN }
        })
        .collect();
    // Only inside the lake's own canonical cells: spreading into neighbouring cells grew
    // lakes past their canonical size (s04: shown at Ultra, hidden at Normal).
    let near_lake = |c: usize| canon_level[c].is_finite();
    let canon_of: Vec<u32> = par_map(n, |i| cg.nearest_node(grid.dir[i]) as u32);
    let mut level = vec![f64::NAN; n];
    let mut closed = vec![false; n];
    let mut queue: std::collections::VecDeque<usize> = std::collections::VecDeque::new();
    for i in 0..n {
        let c = canon_of[i] as usize;
        if !ocean[i] && h[i] < canon_level[c] {
            level[i] = canon_level[c];
            closed[i] = inp.canon_closed.data[cg.node_sample[c] as usize] == 1;
            queue.push_back(i);
        }
    }
    while let Some(i) = queue.pop_front() {
        for &j in grid.neighbors(i) {
            let j = j as usize;
            if ocean[j] || level[j].is_finite() || h[j] >= level[i] || !near_lake(canon_of[j] as usize) {
                continue;
            }
            level[j] = level[i];
            closed[j] = closed[i];
            queue.push_back(j);
        }
    }
    // Canonical lake sizes: lakes below MIN_VISIBLE_LAKE_KM2 become playas.
    let small_lake: Vec<bool> = {
        let mut comp_small = vec![false; cg.n];
        let mut seen = vec![false; cg.n];
        let cell_km2 = |c: usize| cg.weight[c] as f64 * 4.0 * std::f64::consts::PI * (radius_m * 1e-3).powi(2);
        for start in 0..cg.n {
            if seen[start] || !canon_level[start].is_finite() {
                continue;
            }
            seen[start] = true;
            let (mut members, mut stack, mut km2) = (vec![start], vec![start], 0.0);
            while let Some(c) = stack.pop() {
                km2 += cell_km2(c);
                for &o in cg.neighbors(c) {
                    let o = o as usize;
                    if !seen[o] && canon_level[o].is_finite() {
                        seen[o] = true;
                        members.push(o);
                        stack.push(o);
                    }
                }
            }
            if km2 < MIN_VISIBLE_LAKE_KM2 {
                for c in members {
                    comp_small[c] = true;
                }
            }
        }
        (0..n).map(|i| level[i].is_finite() && comp_small[canon_of[i] as usize]).collect()
    };
    let flooded: Vec<bool> = (0..n).map(|i| level[i] > h[i]).collect();
    for i in (0..n).filter(|&i| flooded[i] && small_lake[i]) {
        h[i] = level[i];
    }
    let lake: Vec<bool> = (0..n).map(|i| flooded[i] && !small_lake[i]).collect();
    // Closed lakes are outlets (their water ends there); open ones pass flow on.
    let sink: Vec<bool> = (0..n).map(|i| ocean[i] || (flooded[i] && closed[i])).collect();
    let base: Vec<bool> = (0..n).map(|i| ocean[i] || flooded[i]).collect();
    let base_h: Vec<f64> = (0..n).map(|i| if flooded[i] { level[i] } else { 0.0 }).collect();

    let k = lithology(&grid, radius_m, p);
    let wet: Vec<f64> = match inp.runoff {
        Some(r) => (0..n).map(|i| area_m2[i] * r.sample(grid.dir[i]).clamp(0.02, 6.0) as f64).collect(),
        None => area_m2.clone(),
    };
    let mut net = Network::new(n);
    net.rho8_seed = p.rho8.then_some(p.seed);
    let uplift = vec![0.0; n];
    let mut length = vec![0.0f64; n];

    // The routing surface: lakes at their level, every other pit filled to its spill.
    let routed = |h: &[f64]| -> Vec<f64> {
        let mut f: Vec<f64> = (0..n).map(|i| if flooded[i] { level[i] } else { h[i] }).collect();
        priority_flood_fill_to(&grid, &mut f, &sink, Some(&base_h));
        f
    };

    // 3. Incision without uplift.
    let steps = p.refine_steps.max(1);
    let dt = p.refine_myr.max(0.0) * 1.0e6 / steps as f64;
    if dt > 0.0 {
        for _ in 0..steps {
            let f = routed(&h);
            // Filled pits are sediment: the rock rises to the routing surface (not lakes).
            for i in (0..n).filter(|&i| !base[i]) {
                h[i] = f[i];
            }
            net.route(&grid, &f, &sink);
            net.build_stack();
            net.accumulate(&wet, &sink);
            for (i, l) in length.iter_mut().enumerate() {
                *l = (grid.dir[i] - grid.dir[net.receiver[i] as usize]).length() * radius_m;
            }
            let col = Column {
                receiver: &net.receiver,
                stack: &net.stack,
                drainage_m2: &net.area,
                length_m: &length,
                uplift_m_yr: &uplift,
                base: &base,
                cell_area_m2: &area_m2,
            };
            implicit_step_based(&col, &mut h, |i| k[i], p.m, dt, &base_h);
        }
    }

    // Keep only the sub-cell part of everything added to the canonical surface (tectonic
    // detail and incision): subtract its canonical-cell average. `D` alone doesn't average
    // to zero over a cell (s04 test: 919 m rms at 64 → 128 without this).
    {
        let (mut sum, mut wsum) = (vec![0.0f64; cg.n], vec![0.0f64; cg.n]);
        for i in (0..n).filter(|&i| !base[i]) {
            let c = canon_of[i] as usize;
            sum[c] += (h[i] - canon_up[i]) * area_m2[i];
            wsum[c] += area_m2[i];
        }
        let mut mean = CubeField::<f32>::new(inp.canon_h.res);
        for (k, v) in mean.data.iter_mut().enumerate() {
            let c = cg.sample_node[k] as usize;
            *v = if wsum[c] > 0.0 { (sum[c] / wsum[c]) as f32 } else { 0.0 };
        }
        let lift: Vec<f64> = par_map(n, |i| if base[i] { 0.0 } else { mean.sample_smooth(grid.dir[i]) });
        for i in (0..n).filter(|&i| !base[i]) {
            h[i] = (h[i] - lift[i]).max(1.0);
        }
    }

    // Final surface: drained everywhere except under lakes (the bed stays the bed).
    let f = routed(&h);
    for i in (0..n).filter(|&i| !base[i]) {
        h[i] = f[i];
    }
    net.route(&grid, &f, &sink);
    net.build_stack();
    net.accumulate(&wet, &sink);

    let to_field = |vals: &dyn Fn(usize) -> f32| -> CubeField<f32> {
        let mut out = CubeField::<f32>::new(res);
        for (k, v) in out.data.iter_mut().enumerate() {
            *v = vals(grid.sample_node[k] as usize);
        }
        out
    };
    let height_m = to_field(&|i| h[i] as f32);
    let drainage_km2 = to_field(&|i| if ocean[i] { 0.0 } else { (net.area[i] * 1.0e-6) as f32 });
    let lake_depth_m = to_field(&|i| if lake[i] { (level[i] - h[i]) as f32 } else { 0.0 });
    let mut flow_dir = CubeField::<[f32; 3]>::new(res);
    for (kk, v) in flow_dir.data.iter_mut().enumerate() {
        let i = grid.sample_node[kk] as usize;
        let r = net.receiver[i] as usize;
        if !ocean[i] && r != i {
            let d = grid.dir[r] - grid.dir[i];
            *v = (d - grid.dir[i] * d.dot(grid.dir[i])).normalize_or_zero().as_vec3().to_array();
        }
    }
    let min_area = p.river_min_area_km2 * 1.0e6;
    let rivers = (0..n)
        .filter(|&i| !ocean[i] && net.area[i] >= min_area)
        .map(|i| RiverSegment { from: grid.dir[i], to: grid.dir[net.receiver[i] as usize], area_km2: (net.area[i] * 1.0e-6) as f32 })
        .collect();
    Refined { height_m, drainage_km2, flow_dir, lake_depth_m, rivers }
}
