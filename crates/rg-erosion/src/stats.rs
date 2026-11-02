//! Geomorphic sanity checks: do the generated landscapes obey the laws real ones do?
//!
//! - **Hack's law** `L = c · A^h`: the length of the longest stream above a point grows
//!   with drainage area at h ≈ 0.45–0.7 (natural basins; Schott et al. 2023 report
//!   0.54–0.58 for their generated terrains).
//! - **Concavity** θ: channel slope falls with drainage area, `S ∝ A^−θ`, θ ≈ 0.45–0.5.
//! - **Hypsometric integral** HI = (mean − min) / (max − min) of land: 0.3–0.6 for
//!   real landscapes (youthful high, old low).
//! - **Every land node drains to the sea** (after filling): a hard invariant.

use crate::solver::{ErosionParams, Network};
use rg_tectonics::grid::Grid;

#[derive(Clone, Debug, Default, PartialEq)]
pub struct ErosionStats {
    /// Hack exponent h (log–log fit of mainstem length vs. area on channel nodes).
    pub hack_h: f64,
    /// Hack coefficient c, with L in km and A in km².
    pub hack_c: f64,
    /// Concavity θ (negative slope of log S vs. log A on channel nodes).
    pub concavity: f64,
    /// Uplift-normalised concavity: slope of log(S·K/U) vs. log A. At steady state with
    /// n = 1 this must equal m exactly, so it checks the solver, not the scenery.
    pub concavity_normalised: f64,
    /// Hypsometric integral of all land together (Earth's land ≈ 0.095).
    pub hypsometric_integral: f64,
    /// Median hypsometric integral over drainage basins of ≥ 200 nodes (real basins:
    /// 0.3–0.6; old, worn-down basins lower, young ones higher).
    pub basin_hypsometric_median: f64,
    /// Number of basins in that median.
    pub basins_measured: usize,
    /// Mean land elevation after rescale, m (Earth ≈ 840 m).
    pub mean_land_m: f64,
    /// Share of land nodes whose flow path reaches the sea (must be 1.0).
    pub drained_fraction: f64,
    /// Largest drainage basin, km².
    pub largest_basin_km2: f64,
    /// How many steps needed a Priority-Flood fill.
    pub fill_steps: u32,
    /// The final amplitude rescale applied to land.
    pub amplitude_scale: f64,
    pub grid_res: u32,
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn measure(
    grid: &Grid,
    h: &[f64],
    ocean: &[bool],
    net: &Network,
    cell_area: &[f64],
    uplift: &[f64],
    radius_m: f64,
    p: &ErosionParams,
    fill_steps: u32,
    scale: f64,
) -> ErosionStats {
    let n = grid.n;
    // Longest upstream flow length, by walking the stack from the top down.
    let mut lup = vec![0.0f64; n];
    for &i in net.stack.iter().rev() {
        let i = i as usize;
        let r = net.receiver[i] as usize;
        if r != i && !ocean[r] {
            let l = (grid.dir[i] - grid.dir[r]).length() * radius_m;
            lup[r] = lup[r].max(lup[i] + l);
        }
    }
    // Channel nodes: well above a few cells of drainage area.
    let typical_cell = cell_area.iter().zip(ocean).filter(|(_, o)| !**o).map(|(a, _)| *a).sum::<f64>()
        / ocean.iter().filter(|o| !**o).count().max(1) as f64;
    let channel = |i: usize| !ocean[i] && net.area[i] > 20.0 * typical_cell;

    let (mut hx, mut hy, mut sx, mut sy) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    let (mut nx, mut ny) = (Vec::new(), Vec::new());
    for i in 0..n {
        if !channel(i) {
            continue;
        }
        let a_km2 = net.area[i] * 1.0e-6;
        if lup[i] > 0.0 {
            hx.push(a_km2.ln());
            hy.push((lup[i] * 1.0e-3).ln());
        }
        let r = net.receiver[i] as usize;
        let hr = if ocean[r] { 0.0 } else { h[r] };
        let l = (grid.dir[i] - grid.dir[r]).length() * radius_m;
        let s = (h[i] - hr) / l.max(1.0);
        if s > 1.0e-6 {
            sx.push(a_km2.ln());
            sy.push(s.ln());
            if uplift[i] > 0.0 {
                nx.push(a_km2.ln());
                ny.push((s / uplift[i]).ln());
            }
        }
    }
    let (hack_h, hack_ln_c) = fit(&hx, &hy);
    let (slope_sa, _) = fit(&sx, &sy);
    let (slope_norm, _) = fit(&nx, &ny);

    let land: Vec<f64> = (0..n).filter(|&i| !ocean[i]).map(|i| h[i]).collect();
    let (lo, hi) = land.iter().fold((f64::MAX, f64::MIN), |(a, b), &v| (a.min(v), b.max(v)));
    let mean = land.iter().sum::<f64>() / land.len().max(1) as f64;
    let hi_int = if hi > lo { (mean - lo) / (hi - lo) } else { 0.0 };

    // Drainage check: follow receivers to the sea (bounded walk).
    let mut drained = 0usize;
    for i in (0..n).filter(|&i| !ocean[i]) {
        let mut c = i;
        for _ in 0..n {
            let r = net.receiver[c] as usize;
            if ocean[r] {
                drained += 1;
                break;
            }
            if r == c {
                break;
            }
            c = r;
        }
    }
    let largest = net.area.iter().cloned().fold(0.0, f64::max) * 1.0e-6;

    // Per-basin hypsometry: label each land node with its outlet (the last land node before
    // the sea), walking the stack from the coast inland.
    let mut basin = vec![u32::MAX; n];
    for &i in &net.stack {
        let i = i as usize;
        if ocean[i] {
            continue;
        }
        let r = net.receiver[i] as usize;
        basin[i] = if ocean[r] || r == i { i as u32 } else { basin[r] };
    }
    let mut per: std::collections::BTreeMap<u32, (f64, f64, f64, usize)> = Default::default();
    for i in 0..n {
        if basin[i] == u32::MAX {
            continue;
        }
        let e = per.entry(basin[i]).or_insert((f64::MAX, f64::MIN, 0.0, 0));
        e.0 = e.0.min(h[i]);
        e.1 = e.1.max(h[i]);
        e.2 += h[i];
        e.3 += 1;
    }
    let mut his: Vec<f64> = per
        .values()
        .filter(|v| v.3 >= 200 && v.1 > v.0)
        .map(|v| (v.2 / v.3 as f64 - v.0) / (v.1 - v.0))
        .collect();
    his.sort_by(f64::total_cmp);
    let basin_hi = if his.is_empty() { 0.0 } else { his[his.len() / 2] };
    let _ = p;
    ErosionStats {
        hack_h,
        hack_c: hack_ln_c.exp(),
        concavity: -slope_sa,
        concavity_normalised: -slope_norm,
        hypsometric_integral: hi_int,
        basin_hypsometric_median: basin_hi,
        basins_measured: his.len(),
        mean_land_m: mean,
        drained_fraction: drained as f64 / land.len().max(1) as f64,
        largest_basin_km2: largest,
        fill_steps,
        amplitude_scale: scale,
        grid_res: grid.res,
    }
}

/// Least-squares line y = a·x + b. Returns (a, b).
fn fit(x: &[f64], y: &[f64]) -> (f64, f64) {
    let n = x.len() as f64;
    if n < 3.0 {
        return (0.0, 0.0);
    }
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let sxy: f64 = x.iter().zip(y).map(|(a, b)| (a - mx) * (b - my)).sum();
    let sxx: f64 = x.iter().map(|a| (a - mx) * (a - mx)).sum();
    let a = if sxx > 0.0 { sxy / sxx } else { 0.0 };
    (a, my - a * mx)
}
