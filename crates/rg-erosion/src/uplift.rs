//! Rock uplift rate U(x), derived from **tectonic process**, not from elevation.
//!
//! The first version took `U ∝ (h_tectonic / h95)^0.8`. That is circular: a high plateau far
//! from any plate boundary got the *maximum* uplift just for being high, and came out as the
//! steepest terrain. Real uplift is driven by plate convergence, so:
//!
//! ```text
//! U = U_base + U_belt · exp(−d / L_belt) · w_crust · c
//! ```
//!
//! - `d`: great-circle distance from the node to the nearest **convergent boundary** sample
//!   (`Grid::nearest_source`, which measures true arc distance, not distance along grid lines);
//! - `L_belt`: e-folding width of an orogenic belt (Himalaya, Andes, Alps: 300-500 km);
//! - `w_crust`: 1 on continental crust, less on oceanic crust (island arcs, volcanic islands);
//! - `c`: how fast the two plates close, relative to the planet's mean convergence rate.
//!   For a boundary sample, `c ∝ max(0, (v_a − v_b) · n)`, where `v = ω (e × p)` is a plate's
//!   rigid-rotation velocity at the point, and `n` is the unit direction from plate `a` across
//!   the border into plate `b`. Slow or oblique convergence therefore builds lower mountains.
//!
//! A plateau far from convergent boundaries gets only `U_base`.

use crate::solver::{ErosionParams, UpliftSource};
use glam::DVec3;
use rg_geo::{CubeFace, CubeField};
use rg_tectonics::boundary_class;
use rg_tectonics::grid::Grid;
use rg_tectonics::TectonicSnapshot;

/// Uplift rate in m/yr at every node (0 on ocean nodes), and the source actually used (the
/// tectonic one falls back to the elevation one when the snapshot has no convergent boundary).
pub(crate) fn compute(
    grid: &Grid,
    snap: &TectonicSnapshot,
    radius_m: f64,
    h0: &[f64],
    ocean: &[bool],
    p: &ErosionParams,
) -> (Vec<f64>, UpliftSource) {
    let (mut u, src) = if p.uplift_source == UpliftSource::Tectonic
        && let Some(u) = from_boundaries(grid, snap, radius_m, ocean, p)
    {
        (u, UpliftSource::Tectonic)
    } else {
        (from_elevation(h0, ocean, p), UpliftSource::Elevation)
    };
    // Continental dynamic topography as a rate: sags subside and swells rise over the run,
    // so arid sags stay closed basins against erosion (s01: as a static offset, 30 Myr of
    // erosion captured most closed basins, 14.7% → 2.3% endorheic).
    if let Some(dt) = &snap.dynamic_topo_m {
        let years = (p.duration_myr * 1.0e6).max(1.0);
        for i in 0..u.len() {
            if !ocean[i] {
                u[i] += dt.sample(grid.dir[i]) as f64 / years;
            }
        }
    }
    (u, src)
}

/// Legacy / fallback: uplift ∝ (h / h95)^exponent. Circular (see module docs), but needs no
/// boundary data.
pub(crate) fn from_elevation(h0: &[f64], ocean: &[bool], p: &ErosionParams) -> Vec<f64> {
    let mut land: Vec<f64> = (0..h0.len()).filter(|&i| !ocean[i]).map(|i| h0[i]).collect();
    land.sort_by(f64::total_cmp);
    let p95 = land.get(((land.len().max(1) - 1) as f64 * 0.95).round() as usize).copied().unwrap_or(1.0).max(1.0);
    let mm = 1.0e-3;
    (0..h0.len())
        .map(|i| {
            if ocean[i] {
                return 0.0;
            }
            let w = (h0[i] / p95).clamp(0.0, 1.0);
            (p.uplift_base_mm_yr + (p.uplift_max_mm_yr - p.uplift_base_mm_yr).max(0.0) * w.powf(p.uplift_exponent)) * mm
        })
        .collect()
}

fn from_boundaries(grid: &Grid, snap: &TectonicSnapshot, radius_m: f64, ocean: &[bool], p: &ErosionParams) -> Option<Vec<f64>> {
    let n = grid.n;
    // The solve grid may be coarser than the snapshot: map each convergent snapshot sample
    // to its nearest solve-grid node, so thin boundary lines are never missed.
    let probe = CubeField::<u8>::new(grid.res);
    let node_of = |d: DVec3| -> usize {
        let (face, u, v) = CubeFace::from_dir(d);
        let r = grid.res as f64;
        let i = ((u + 1.0) * 0.5 * r).round().clamp(0.0, r) as u32;
        let j = ((v + 1.0) * 0.5 * r).round().clamp(0.0, r) as u32;
        grid.sample_node[probe.index(face, i, j)] as usize
    };

    let b = &snap.boundary;
    let (mut rate_sum, mut rate_cnt) = (vec![0.0f64; n], vec![0u32; n]);
    let mut is_seed = vec![false; n];
    let res = b.res;
    for (idx, face, i, j, dir) in b.iter_samples() {
        if b.data[idx] != boundary_class::CONVERGENT {
            continue;
        }
        let node = node_of(dir);
        is_seed[node] = true;
        // Convergence speed from the plates' rigid rotations, using an in-face neighbour that
        // belongs to the other plate to get the direction across the border.
        let pa = snap.plate_id.data[idx] as usize;
        for (di, dj) in [(1i32, 0i32), (-1, 0), (0, 1), (0, -1)] {
            let (ii, jj) = (i as i32 + di, j as i32 + dj);
            if ii < 0 || jj < 0 || ii > res as i32 || jj > res as i32 {
                continue;
            }
            let k = snap.plate_id.index(face, ii as u32, jj as u32);
            let pb = snap.plate_id.data[k] as usize;
            if pb == pa || pa >= snap.plates.len() || pb >= snap.plates.len() {
                continue;
            }
            let dn = b.sample_dir(face, ii as u32, jj as u32);
            let t = dn - dir * dn.dot(dir);
            if t.length_squared() < 1e-24 {
                continue;
            }
            let nrm = t.normalize();
            let (a, c) = (&snap.plates[pa], &snap.plates[pb]);
            let va = a.euler_pole.cross(dir) * (a.omega * radius_m);
            let vb = c.euler_pole.cross(dir) * (c.omega * radius_m);
            rate_sum[node] += (va - vb).dot(nrm).max(0.0);
            rate_cnt[node] += 1;
            break;
        }
    }
    let seeds: Vec<u32> = (0..n as u32).filter(|&i| is_seed[i as usize]).collect();
    if seeds.is_empty() {
        return None;
    }

    // Convergence-rate factor per seed: rate relative to the mean over seeds that have one.
    let rates: Vec<f64> = seeds
        .iter()
        .filter(|&&s| rate_cnt[s as usize] > 0)
        .map(|&s| rate_sum[s as usize] / rate_cnt[s as usize] as f64)
        .collect();
    let mean_rate = rates.iter().sum::<f64>() / rates.len().max(1) as f64;
    let factor = |s: usize| -> f64 {
        if !p.scale_by_convergence_rate || rate_cnt[s] == 0 || mean_rate <= 0.0 {
            1.0
        } else {
            (rate_sum[s] / rate_cnt[s] as f64 / mean_rate).clamp(0.25, 2.5)
        }
    };

    let l_belt = (p.belt_width_km * 1.0e3).max(1.0);
    let max_rad = ((6.0 * l_belt / radius_m).min(std::f64::consts::PI)) as f32;
    let near = grid.nearest_source(&seeds, max_rad);
    let mm = 1.0e-3;
    Some(
        (0..n)
            .map(|i| {
                if ocean[i] {
                    return 0.0;
                }
                let d = near.dist_rad[i] as f64 * radius_m;
                let belt = if near.src[i] == u32::MAX { 0.0 } else { (-d / l_belt).exp() * factor(near.src[i] as usize) };
                let w_crust = if snap.crust.nearest(grid.dir[i]) == 1 { 1.0 } else { p.oceanic_crust_uplift };
                (p.uplift_base_mm_yr + p.uplift_max_mm_yr * belt * w_crust) * mm
            })
            .collect(),
    )
}
