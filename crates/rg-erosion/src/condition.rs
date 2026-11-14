//! Drainage conditioning for imported heightmaps (s09).
//!
//! A drawn heightmap has no river valleys: its slopes are smooth ramps, its plains are
//! exactly flat (8-bit grey steps) and it is full of small pits. Routed as it is, rivers
//! run ruler-straight and parallel down the ramps and every pit becomes a lake.
//!
//! This builds the surface the **final routing and the lakes** see on pinned (heightmap)
//! land. The committed heights are not touched.
//!
//! 1. **Shallow pits are drawing noise:** a depression on pinned land no deeper than
//!    `pit_m` is filled, so only basins the map really draws hold a lake.
//! 2. **Wander:** position-keyed fBm, with an amplitude that follows the local relief and a
//!    floor for exactly flat plains. It is *not* added to the surface (that would dig new
//!    pits): the router adds it when it picks a receiver among the neighbours that are
//!    truly lower (`Network::wander`). So water always runs downhill on the real surface,
//!    and the wander only chooses which way down.
//!
//! Position-keyed, so every resolution bends the same rivers the same way.
use crate::jitter::roughness_at;
use crate::solver::priority_flood_fill_eps;
use crate::uplift::smooth_on_sphere;
use rg_tectonics::grid::Grid;

#[derive(Clone, Copy, Debug)]
pub struct Conditioning {
    /// Longest wavelength of the wander, km.
    pub wavelength_km: f64,
    /// Wander amplitude as a multiple of the smoothed drop per cell.
    pub strength: f64,
    /// Smallest wander amplitude, m (flat plains).
    pub floor_m: f64,
    /// Depressions on pinned land up to this deep (m) are filled.
    pub pit_m: f64,
}

impl Default for Conditioning {
    fn default() -> Self {
        let env = |k: &str, d: f64| std::env::var(k).ok().and_then(|v| v.parse().ok()).unwrap_or(d);
        Self {
            wavelength_km: env("RG_COND_WL", 400.0),
            strength: env("RG_COND_AMP", 2.0),
            floor_m: env("RG_COND_FLOOR", 3.0),
            pit_m: env("RG_COND_PIT", 40.0),
        }
    }
}

/// The routing surface for `h` (m) with heightmap pin weights `pin` (0..1 per node).
pub fn condition(grid: &Grid, radius_m: f64, h: &[f64], ocean: &[bool], pin: &[f64], seed: u64, c: &Conditioning) -> Conditioned {
    let n = grid.n;
    let cell_m = std::f64::consts::FRAC_PI_2 * radius_m / grid.res as f64;
    // Local relief: the biggest drop to a neighbour, smoothed over ~3 canonical cells.
    let drop: Vec<f64> = (0..n)
        .map(|i| if ocean[i] { 0.0 } else { grid.neighbors(i).iter().map(|&j| (h[i] - h[j as usize].max(0.0)).abs()).fold(0.0, f64::max) })
        .collect();
    let land: Vec<f64> = ocean.iter().map(|&o| if o { 0.0 } else { 1.0 }).collect();
    let sigma = 100_000.0;
    let (sd, sl) = (smooth_on_sphere(grid, &drop, sigma, radius_m), smooth_on_sphere(grid, &land, sigma, radius_m));
    // Basins the map really draws: depressions of the raw import deeper than `pit_m`.
    // Shallower ones are drawing noise and are filled, so only real basins hold a lake.
    let mut raw_fill = h.to_vec();
    priority_flood_fill_eps(grid, &mut raw_fill, ocean, None, 0.0, None);
    let mut basin = vec![false; n];
    let mut seen = vec![false; n];
    for s in 0..n {
        if seen[s] || ocean[s] || raw_fill[s] - h[s] <= 1e-9 {
            continue;
        }
        seen[s] = true;
        let (mut members, mut stack, mut deepest) = (vec![s], vec![s], 0.0f64);
        while let Some(i) = stack.pop() {
            deepest = deepest.max(raw_fill[i] - h[i]);
            for &j in grid.neighbors(i) {
                let j = j as usize;
                if !seen[j] && !ocean[j] && raw_fill[j] - h[j] > 1e-9 {
                    seen[j] = true;
                    members.push(j);
                    stack.push(j);
                }
            }
        }
        if deepest > c.pit_m {
            for i in members {
                basin[i] = true;
            }
        }
    }
    let salt = (seed as u32) ^ 0xC0ED_17E5;
    let wl_m = c.wavelength_km * 1e3;
    let wander: Vec<f64> = (0..n)
        .map(|i| {
            if ocean[i] || pin[i] <= 0.0 {
                return 0.0;
            }
            // Drop per cell → per metre, so the amplitude is the same at every resolution.
            let slope = sd[i] / sl[i].max(1e-9) / cell_m;
            let amp = (c.strength * slope * wl_m * 0.25).max(c.floor_m) * pin[i];
            roughness_at(grid.dir[i], radius_m, wl_m, 2.0 * cell_m, salt) * amp
        })
        .collect();
    let mut surface = h.to_vec();
    for i in 0..n {
        if !ocean[i] && pin[i] >= 0.5 && !basin[i] {
            surface[i] = raw_fill[i];
        }
    }
    // Flats (terraced maps are mostly flats): a 25× spread of flood costs, position-keyed.
    let flat_wl = env_or("RG_COND_FLAT_WL", 250.0) * 1e3;
    let spread = env_or("RG_COND_FLAT_SPREAD", 5.0f64).ln();
    let flat_cost: Vec<f64> = (0..n)
        .map(|i| if ocean[i] || pin[i] <= 0.0 { 1.0 } else { (spread * pin[i] * roughness_at(grid.dir[i], radius_m, flat_wl, 2.0 * cell_m, salt ^ 0x00F1_A700)).exp() })
        .collect();
    Conditioned { surface, wander, flat_cost }
}

/// What [`condition`] returns.
pub struct Conditioned {
    /// The surface lakes are worked out on: the import with its shallow pits filled.
    pub surface: Vec<f64>,
    /// Per node, metres added to the surface when choosing a receiver (0 off the import).
    pub wander: Vec<f64>,
    /// Per node, the multiplier on the flat-routing gradient (1 off the import).
    pub flat_cost: Vec<f64>,
}

fn env_or(k: &str, d: f64) -> f64 {
    std::env::var(k).ok().and_then(|v| v.parse().ok()).unwrap_or(d)
}
