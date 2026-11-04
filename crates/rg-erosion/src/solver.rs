use crate::stats::{self, ErosionStats};
use crate::{hypsometry, uplift};
use glam::DVec3;
use rg_geo::CubeField;
use rg_tectonics::grid::Grid;
use rg_tectonics::util::par_map;
use rg_tectonics::TectonicSnapshot;
use std::cmp::Ordering;
use std::collections::{BinaryHeap, VecDeque};

/// Where the rock uplift rate comes from.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum UpliftSource {
    /// Uplift decays with distance from convergent plate boundaries, scaled by crust type and
    /// the plates' closing speed (see `uplift.rs`). Falls back to `Elevation` if the
    /// snapshot has no convergent boundary.
    #[default]
    Tectonic,
    /// Legacy: uplift ∝ (tectonic height / h95)^`uplift_exponent`. Circular (plateaus get the
    /// maximum uplift), kept only as a fallback for data with no boundary information.
    Elevation,
}

#[derive(Clone, Debug, PartialEq)]
pub struct ErosionParams {
    /// Run erosion at all.
    pub enabled: bool,
    /// The solve runs on a grid of at most this many cells per face edge; finer macro
    /// fields get the erosion *change* upsampled onto them (Schott 2023: blend changes,
    /// never two eroded surfaces).
    pub max_res: u32,
    /// Erosion grid resolution (cells per face edge). 0 = the tectonics resolution capped at
    /// `max_res`. A value above the tectonics resolution makes the committed macro surface
    /// finer than tectonics: heights are bilinearly sampled onto the finer grid, eroded
    /// there, and returned at that resolution (the coast is the tectonic 0 m contour).
    pub res: u32,
    /// Simulated time, million years. Longer = deeper, more mature valleys.
    pub duration_myr: f64,
    /// Implicit time steps (the scheme is stable at any step; more steps let the drainage
    /// network reorganise as the surface changes).
    pub steps: u32,
    /// Erodibility K, in m^(1−2m)/yr. Typical bedrock 1e-7…1e-5.
    pub k: f64,
    /// Area exponent m (with n = 1 the concavity θ = m).
    pub m: f64,
    /// Peak uplift rate **above the base rate** at a convergent boundary (`Tectonic`), or the
    /// total peak rate (`Elevation`), mm/yr. Himalaya ~5–10, Andes ~1–3.
    pub uplift_max_mm_yr: f64,
    /// Background uplift on stable continent, mm/yr.
    pub uplift_base_mm_yr: f64,
    /// `Elevation` source only: U ∝ (h / h95)^exponent.
    pub uplift_exponent: f64,
    /// Where uplift comes from (see [`UpliftSource`]).
    pub uplift_source: UpliftSource,
    /// `Tectonic` source: e-folding distance of uplift away from a convergent boundary, km.
    pub belt_width_km: f64,
    /// `Tectonic` source: uplift weight on oceanic crust (continental = 1).
    pub oceanic_crust_uplift: f64,
    /// `Tectonic` source: scale belt uplift by the plates' closing speed (mean = 1).
    pub scale_by_convergence_rate: bool,
    /// `Tectonic` source: how far (km, Gaussian σ) the closing-speed factor is blended
    /// between neighbouring boundary stretches. 0 = each place takes the factor of its single
    /// nearest boundary point, which steps the uplift along straight Voronoi lines and erodes
    /// into flat "floors" with escarpments between them (s03; Matt: "interesting, not
    /// natural"). Larger values ramp smoothly instead. Old world files load as 0.
    pub convergence_smoothing_km: f64,
    /// Lithology: standard deviation of ln K across the planet (0 = uniform rock). Default 0.3:
    /// at 13–30 km cells rock strength is averaged over many units; 0.6 pushed basin
    /// hypsometry and Hack's exponent to the edge of Earth's range (s01 log). Rock
    /// erodibility varies ~10× between lithologies (e.g. Stock & Montgomery 1999); rivers
    /// route around resistant rock, which breaks up straight, parallel drainage.
    pub k_variability: f64,
    /// Correlation length of the lithology field, km.
    pub lithology_scale_km: f64,
    /// Roughness of the initial land surface, m (white, per node). Braun & Willett (2013)
    /// start from small random noise; on a perfectly smooth tilted surface the drainage
    /// network has nothing to organise around.
    pub initial_roughness_m: f64,
    /// R1 spike (s03): wavelength (km) of the initial roughness as position-keyed fBm, down
    /// to two cells. 0 = the original white noise per node, which is different at every
    /// resolution (one cause of "terrain differs between resolutions", s02).
    pub roughness_wavelength_km: f64,
    /// R1 spike (s03): move every erosion node up to ±this many cells (0..0.4) so rivers
    /// aren't bound to the lattice's 8 directions; results are resampled back to the raster
    /// (`jitter.rs`). Replaces Rho8 when on. 0 = the regular grid.
    pub node_jitter: f64,
    /// Route with Rho8 (Fairfield & Leymarie 1991) instead of plain D8: each node's diagonal
    /// neighbours get a random effective distance in [1, 2] axis spacings instead of √2.
    /// Plain D8 sends every node of a smooth tilted surface the same grid direction, which
    /// gives long, straight, parallel rivers; Rho8 follows the true aspect on average.
    pub rho8: bool,
    /// Seed for Rho8 (the pipeline sets it from the world seed).
    pub seed: u64,
    /// Smallest drainage area drawn as a river, km².
    pub river_min_area_km2: f64,
    /// Hypsometric calibration: land heights are quantile-mapped onto the land hypsometric
    /// curve (`hypsometry.rs`) scaled by this factor (1.0 = Earth-like; 2.0 = everything
    /// twice as high). 0 disables it (land is then linearly rescaled to the tectonic 99th
    /// percentile).
    pub hypsometry_scale: f64,
    /// Physical calibration instead of the quantile remap: if > 0, land heights from the
    /// solve are scaled by ONE factor so the area-weighted mean land height equals this (m,
    /// × `hypsometry_scale`). At stream-power steady state h ≈ (U/K)·χ is linear in 1/K, so
    /// this is fitting K; the hypsometric curve's *shape* then comes from the physics and is
    /// only checked against Earth (s01 review: "heights fitted, not modelled").
    pub linear_mean_m: f64,
    /// Elevation (m, solve units) at which rock uplift has tapered to zero:
    /// U_eff = U·(1 − h/limit). Lithospheric strength caps plateau elevation (Earth's highest
    /// plateaus ~5 km); 0 disables.
    pub elevation_limit_m: f64,
    /// Lakes and closed (endorheic) basins from the water balance (needs climate inputs).
    pub lakes: bool,
}

impl Default for ErosionParams {
    fn default() -> Self {
        Self {
            enabled: true,
            max_res: 2048,
            res: 0,
            duration_myr: 30.0,
            steps: 20,
            k: 1.0e-5,
            m: 0.45,
            uplift_max_mm_yr: 2.5,
            uplift_base_mm_yr: 0.15,
            uplift_exponent: 0.8,
            uplift_source: UpliftSource::Tectonic,
            belt_width_km: 400.0,
            oceanic_crust_uplift: 0.5,
            scale_by_convergence_rate: true,
            convergence_smoothing_km: 200.0,
            k_variability: 0.3,
            lithology_scale_km: 250.0,
            initial_roughness_m: 20.0,
            roughness_wavelength_km: 0.0,
            node_jitter: 0.0,
            rho8: true,
            seed: 0,
            river_min_area_km2: 20_000.0,
            hypsometry_scale: 1.0,
            linear_mean_m: 656.0,
            elevation_limit_m: 7000.0,
            lakes: true,
        }
    }
}

/// A river reach from a node to the node it drains into (unit directions).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RiverSegment {
    pub from: DVec3,
    pub to: DVec3,
    /// Drainage area at the upstream end, km².
    pub area_km2: f32,
}

pub struct ErosionResult {
    /// Eroded heights (m above sea level), at the erosion grid's resolution when that is at
    /// least the input's, otherwise at the input's resolution.
    pub height_m: CubeField<f32>,
    /// Upstream drainage area, km², at the erosion grid's resolution (0 at sea).
    pub drainage_km2: CubeField<f32>,
    /// Downstream unit flow direction (tangent), xyz, at the erosion grid's resolution.
    /// Zero at sea and on flats.
    pub flow_dir: CubeField<[f32; 3]>,
    /// Depth of standing water in closed basins at the end (lakes), m.
    pub lake_depth_m: CubeField<f32>,
    /// Rock uplift rate, mm/yr, at the erosion grid's resolution. At stream-power steady
    /// state it is also the erosion rate, and `k_sn = U/K`: streamed detail reads it.
    pub uplift_mm_yr: CubeField<f32>,
    /// With linear calibration (`linear_mean_m` > 0): the single factor from solved to
    /// committed heights, so slopes scale by it and the committed surface's steepness is
    /// k_sn = factor·U/K. `None` with the (non-linear) quantile remap.
    pub height_scale: Option<f64>,
    pub rivers: Vec<RiverSegment>,
    pub stats: ErosionStats,
}

// ---------------------------------------------------------------------------------------
// The numerical core: arrays in, arrays out. `erode()` and the tests both call these.
// ---------------------------------------------------------------------------------------

/// Per-node erodibility: `K · exp(σ · z)`, with z a smooth noise field of unit variance and
/// correlation length `lithology_scale_km` (log-normal, as rock strength data are).
fn lithology(grid: &Grid, radius_m: f64, p: &ErosionParams) -> Vec<f64> {
    if p.k_variability <= 0.0 {
        return vec![p.k; grid.n];
    }
    // Value noise on a jittered 3D lattice is enough here; smoothness matters, not isotropy.
    let cell = p.lithology_scale_km * 1000.0 / radius_m; // lattice spacing on the unit sphere
    let seed = p.seed ^ 0x0001_1740_106E;
    let lattice = |x: i64, y: i64, z: i64| {
        let key = (x as u64).wrapping_mul(73_856_093) ^ (y as u64).wrapping_mul(19_349_663) ^ (z as u64).wrapping_mul(83_492_791);
        unit_hash(seed, key) * 2.0 - 1.0
    };
    let smooth = |t: f64| t * t * (3.0 - 2.0 * t);
    let value = |q: DVec3| {
        let b = q.floor();
        let f = q - b;
        let (x, y, z) = (b.x as i64, b.y as i64, b.z as i64);
        let (u, v, w) = (smooth(f.x), smooth(f.y), smooth(f.z));
        let l = |dx, dy, dz| lattice(x + dx, y + dy, z + dz);
        let lerp = |a: f64, b: f64, t: f64| a + (b - a) * t;
        lerp(
            lerp(lerp(l(0, 0, 0), l(1, 0, 0), u), lerp(l(0, 1, 0), l(1, 1, 0), u), v),
            lerp(lerp(l(0, 0, 1), l(1, 0, 1), u), lerp(l(0, 1, 1), l(1, 1, 1), u), v),
            w,
        )
    };
    // Two octaves: the correlation length and half of it.
    let z: Vec<f64> = par_map(grid.n, |i| {
        let q = grid.dir[i] / cell;
        value(q) + 0.5 * value(q * 2.0 + DVec3::splat(17.3))
    });
    let mean = z.iter().sum::<f64>() / z.len() as f64;
    let sd = (z.iter().map(|v| (v - mean) * (v - mean)).sum::<f64>() / z.len() as f64).sqrt().max(1e-9);
    z.iter().map(|v| p.k * (p.k_variability * (v - mean) / sd).exp()).collect()
}

/// Which neighbours are diagonal, by topology rather than distance (the tangent warp makes
/// cells up to 1:√2 rectangles, so a length test misclassifies near face edges). In an
/// 8-connected quad grid, axis neighbours share 4 common neighbours and diagonal ones share
/// 2. Bit k of entry `a` is set when slot k of `grid.neighbors(a)` is diagonal.
pub fn diagonal_mask(grid: &Grid) -> Vec<u16> {
    par_map(grid.n, |a| {
        let na = grid.neighbors(a);
        let mut mask = 0u16;
        for (slot, &b) in na.iter().enumerate() {
            let nb = grid.neighbors(b as usize);
            // Both lists are sorted: count common entries with a merge.
            let (mut x, mut y, mut common) = (0, 0, 0);
            while x < na.len() && y < nb.len() {
                match na[x].cmp(&nb[y]) {
                    Ordering::Less => x += 1,
                    Ordering::Greater => y += 1,
                    Ordering::Equal => {
                        common += 1;
                        x += 1;
                        y += 1;
                    }
                }
            }
            if common <= 2 {
                mask |= 1 << slot;
            }
        }
        mask
    })
}

/// Deterministic hash of (seed, i) to [0, 1) (splitmix64 finaliser).
fn unit_hash(seed: u64, i: u64) -> f64 {
    let mut z = seed ^ i.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^= z >> 31;
    (z >> 11) as f64 / (1u64 << 53) as f64
}

/// A routed drainage network and the per-node data one implicit step needs.
pub struct Column<'a> {
    /// Receiver of each node (itself for base-level nodes).
    pub receiver: &'a [u32],
    /// Topological order: every node comes after its receiver.
    pub stack: &'a [u32],
    /// Upstream drainage area including the node's own cell, m².
    pub drainage_m2: &'a [f64],
    /// Distance to the receiver, m.
    pub length_m: &'a [f64],
    /// Rock uplift rate, m/yr.
    pub uplift_m_yr: &'a [f64],
    /// Base-level nodes (ocean): held fixed at height 0 for their donors.
    pub base: &'a [bool],
    /// Cell area, m² (only used for the volume bookkeeping).
    pub cell_area_m2: &'a [f64],
}

/// Volumes moved in one step (or summed over a run), m³, over non-base nodes.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MassBudget {
    /// Rock raised by uplift: Σ U·Δt·area.
    pub uplift_m3: f64,
    /// Rock removed by rivers: Σ K·A^m·S·Δt·area, with S the *new* slope.
    pub eroded_m3: f64,
    /// Material added by Priority-Flood filling of closed basins (sediment fill).
    pub filled_m3: f64,
}

impl MassBudget {
    /// The net change in land volume the books say happened.
    pub fn net_m3(&self) -> f64 {
        self.uplift_m3 - self.eroded_m3 + self.filled_m3
    }
    fn add(&mut self, o: MassBudget) {
        self.uplift_m3 += o.uplift_m3;
        self.eroded_m3 += o.eroded_m3;
        self.filled_m3 += o.filled_m3;
    }
}

/// One implicit stream-power step, `n = 1` (Braun & Willett 2013):
///
/// `h_i ← (h_i + U_i·Δt + F_i·h_r) / (1 + F_i)`, `F_i = K·Δt·A_i^m / L_i`,
///
/// visiting nodes in stack order so the receiver `h_r` is already updated. Rearranged,
/// `Δh_i = U_i·Δt − F_i·(h_i' − h_r')`, and the last term is exactly the river's
/// incision `K·A^m·S'·Δt` with the *new* slope `S'`. That is what we sum as eroded volume.
pub fn implicit_step(c: &Column, h: &mut [f64], k: f64, m: f64, dt: f64) -> MassBudget {
    implicit_step_varying(c, h, |_| k, m, dt)
}

/// As [`implicit_step`], with erodibility `k(node)` varying in space (lithology).
pub fn implicit_step_varying(c: &Column, h: &mut [f64], k: impl Fn(usize) -> f64, m: f64, dt: f64) -> MassBudget {
    let mut b = MassBudget::default();
    for &i in c.stack {
        let i = i as usize;
        if c.base[i] {
            continue;
        }
        let r = c.receiver[i] as usize;
        let hr = if c.base[r] { 0.0 } else { h[r] };
        let f = k(i) * dt * c.drainage_m2[i].powf(m) / c.length_m[i].max(1.0);
        let new = (h[i] + c.uplift_m_yr[i] * dt + f * hr) / (1.0 + f);
        b.uplift_m3 += c.uplift_m_yr[i] * dt * c.cell_area_m2[i];
        b.eroded_m3 += f * (new - hr) * c.cell_area_m2[i];
        h[i] = new;
    }
    b
}

/// As [`implicit_step_varying`], with base-level nodes held at `base_h` (0 for the ocean,
/// the water level for lakes) instead of 0, and base nodes uplifted but never eroded.
fn implicit_step_based(c: &Column, h: &mut [f64], k: impl Fn(usize) -> f64, m: f64, dt: f64, base_h: &[f64]) -> MassBudget {
    let mut b = MassBudget::default();
    for &i in c.stack {
        let i = i as usize;
        if c.base[i] {
            continue;
        }
        let r = c.receiver[i] as usize;
        let hr = if c.base[r] { base_h[r] } else { h[r] };
        let f = k(i) * dt * c.drainage_m2[i].powf(m) / c.length_m[i].max(1.0);
        let new = (h[i] + c.uplift_m_yr[i] * dt + f * hr) / (1.0 + f);
        b.uplift_m3 += c.uplift_m_yr[i] * dt * c.cell_area_m2[i];
        b.eroded_m3 += f * (new - hr) * c.cell_area_m2[i];
        h[i] = new;
    }
    b
}

/// Water inputs for lakes inside the solve: per node runoff and open-water net evaporation
/// (E − P), m/yr.
pub struct LakeInputs<'a> {
    pub runoff_m_yr: &'a [f64],
    pub net_evap_m_yr: &'a [f64],
}

/// What [`solve`] leaves behind.
pub struct Solved {
    /// The drainage network of the last step (receivers, stack, areas).
    pub net: Network,
    /// Volumes over the whole run, so `Σ_land area·(h_end − h_start) = budget.net_m3()`.
    pub budget: MassBudget,
    /// How many steps needed a Priority-Flood fill.
    pub filled_steps: u32,
    /// Lakes at the end of the run (with lake inputs).
    pub lakes: Option<crate::lakes::Lakes>,
}

/// Run `p.steps` implicit steps of duration `p.duration_myr / p.steps` on `h` (metres, with
/// `ocean` nodes as base level), re-routing the network each step.
pub fn solve(grid: &Grid, radius_m: f64, h: &mut [f64], ocean: &[bool], uplift_m_yr: &[f64], p: &ErosionParams) -> Solved {
    solve_with_k(grid, radius_m, h, ocean, uplift_m_yr, &vec![p.k; grid.n], p)
}

/// As [`solve`], with a per-node erodibility `k` (m^(1−2m)/yr) instead of the uniform `p.k`.
pub fn solve_with_k(
    grid: &Grid,
    radius_m: f64,
    h: &mut [f64],
    ocean: &[bool],
    uplift_m_yr: &[f64],
    k: &[f64],
    p: &ErosionParams,
) -> Solved {
    solve_full(grid, radius_m, h, ocean, uplift_m_yr, k, None, None, None, p)
}

/// As [`solve_with_k`], with per-node relative runoff (1 = the land mean): the stream-power
/// "area" becomes runoff-weighted drainage area, a discharge proxy (Q ∝ Σ runoff·area), so
/// wet slopes incise faster than rain shadows and deserts.
#[allow(clippy::too_many_arguments)]
pub fn solve_full(
    grid: &Grid,
    radius_m: f64,
    h: &mut [f64],
    ocean: &[bool],
    uplift_m_yr: &[f64],
    k: &[f64],
    runoff: Option<&[f64]>,
    elevation_limit_m: Option<f64>,
    lake_in: Option<&LakeInputs>,
    p: &ErosionParams,
) -> Solved {
    solve_full_with_progress(grid, radius_m, h, ocean, uplift_m_yr, k, runoff, elevation_limit_m, lake_in, p, &|_, _| true)
        .expect("a solve without a cancelling callback always finishes")
}

/// As [`solve_full`], calling `on_step(k, n)` as implicit step k of n (1-based) starts.
/// If it returns false the solve stops there and returns `None` (cancelled; `h` is left
/// part-way). Called on this thread only.
/// Called on this thread only (progress for the app's console, s02).
#[allow(clippy::too_many_arguments)]
pub fn solve_full_with_progress(
    grid: &Grid,
    radius_m: f64,
    h: &mut [f64],
    ocean: &[bool],
    uplift_m_yr: &[f64],
    k: &[f64],
    runoff: Option<&[f64]>,
    elevation_limit_m: Option<f64>,
    lake_in: Option<&LakeInputs>,
    p: &ErosionParams,
    on_step: &dyn Fn(u32, u32) -> bool,
) -> Option<Solved> {
    let n = grid.n;
    let cell_area: Vec<f64> = grid.weight.iter().map(|w| *w as f64 * 4.0 * std::f64::consts::PI * radius_m * radius_m).collect();
    let wet_area: Vec<f64> = match runoff {
        Some(r) => cell_area.iter().zip(r).map(|(a, r)| a * r).collect(),
        None => cell_area.clone(),
    };
    let mut net = Network::new(n);
    // Jitter already breaks the 8-direction bias (and Rho8's diagonal test assumes a lattice).
    net.rho8_seed = (p.rho8 && p.node_jitter <= 0.0).then_some(p.seed);
    let mut budget = MassBudget::default();
    let mut filled_steps = 0u32;
    let steps = p.steps.max(1);
    let dt = p.duration_myr * 1.0e6 / steps as f64;
    let mut length = vec![0.0f64; n];
    let mut lakes_now: Option<crate::lakes::Lakes> = None;
    for step in 0..steps {
        if !on_step(step + 1, steps) {
            return None;
        }
        // Base level: the ocean, plus lake water surfaces when lakes are on. With lakes the
        // rock is never filled (no "sediment from nowhere"); routing uses a filled copy with
        // closed lakes as outlets, so all water reaches the sea or a closed lake.
        let (sink, base, base_h, route_h): (Vec<bool>, Vec<bool>, Vec<f64>, Option<Vec<f64>>) = match lake_in {
            Some(li) => {
                let lk = crate::lakes::compute(grid, h, ocean, li.runoff_m_yr, li.net_evap_m_yr, &cell_area);
                let sub: Vec<bool> = (0..n).map(|i| !ocean[i] && lk.level[i] > h[i]).collect();
                let sink: Vec<bool> = (0..n).map(|i| ocean[i] || lk.closed_lake[i]).collect();
                let base: Vec<bool> = (0..n).map(|i| ocean[i] || sub[i]).collect();
                let base_h: Vec<f64> = (0..n).map(|i| if sub[i] { lk.level[i] } else { 0.0 }).collect();
                let mut f: Vec<f64> = (0..n).map(|i| if sub[i] { lk.level[i] } else { h[i] }).collect();
                priority_flood_fill_to(grid, &mut f, &sink, Some(&base_h));
                lakes_now = Some(lk);
                (sink, base, base_h, Some(f))
            }
            None => (ocean.to_vec(), ocean.to_vec(), vec![0.0; n], None),
        };
        // 1. Receivers (without lakes: fill pits first if there are any, and always on step 0).
        match &route_h {
            Some(f) => {
                net.route(grid, f, &sink);
            }
            None => {
                if step == 0 || net.route(grid, h, ocean) > 0 {
                    let before = h.to_vec();
                    priority_flood_fill(grid, h, ocean);
                    budget.filled_m3 += (0..n).filter(|&i| !ocean[i]).map(|i| (h[i] - before[i]) * cell_area[i]).sum::<f64>();
                    filled_steps += 1;
                    let pits = net.route(grid, h, ocean);
                    debug_assert_eq!(pits, 0, "filled surface must drain everywhere");
                }
            }
        }
        // 2–3. Stack and drainage area (runoff-weighted when climate is given). Closed
        // lakes absorb their inflow; open lakes pass it on to their outlet.
        net.build_stack();
        net.accumulate(&wet_area, &sink);
        // 4. Implicit stream-power incision with uplift.
        for (i, l) in length.iter_mut().enumerate() {
            *l = (grid.dir[i] - grid.dir[net.receiver[i] as usize]).length() * radius_m;
        }
        // Elevation limit: crust can't hold relief above what its strength supports; uplift
        // tapers to zero there and further shortening would widen the plateau instead
        // (thin viscous sheet; s01: without it arid, barely eroding highs grew to 25 km).
        let u_eff: Vec<f64>;
        let uplift_now: &[f64] = match elevation_limit_m {
            Some(lim) if lim > 0.0 => {
                u_eff = (0..n).map(|i| uplift_m_yr[i] * (1.0 - h[i].max(0.0) / lim).max(0.0)).collect();
                &u_eff
            }
            _ => uplift_m_yr,
        };
        let col = Column {
            receiver: &net.receiver,
            stack: &net.stack,
            drainage_m2: &net.area,
            length_m: &length,
            uplift_m_yr: uplift_now,
            base: &base,
            cell_area_m2: &cell_area,
        };
        if lake_in.is_some() {
            budget.add(implicit_step_based(&col, h, |i| k[i], p.m, dt, &base_h));
            // Lake beds still rise with the rock (but don't erode).
            for i in (0..n).filter(|&i| base[i] && !ocean[i]) {
                h[i] += uplift_now[i] * dt;
                budget.uplift_m3 += uplift_now[i] * dt * cell_area[i];
            }
        } else {
            budget.add(implicit_step_varying(&col, h, |i| k[i], p.m, dt));
        }
    }
    Some(Solved { net, budget, filled_steps, lakes: lakes_now })
}

/// Erode the tectonic snapshot. Deterministic: same input, same output, on any thread count.
///
/// Heights come from `tectonic.height_m`; uplift comes from the snapshot's boundaries, crust
/// and plate motions (see [`UpliftSource`]).
pub fn erode(tectonic: &TectonicSnapshot, radius_m: f64, p: &ErosionParams) -> ErosionResult {
    erode_with_runoff(tectonic, radius_m, p, None)
}

/// As [`erode`], with relative runoff (1 = land mean; any resolution, sampled per node).
/// Drainage areas (and the stored `drainage_km2`) are then runoff-weighted.
pub fn erode_with_runoff(tectonic: &TectonicSnapshot, radius_m: f64, p: &ErosionParams, runoff: Option<&CubeField<f32>>) -> ErosionResult {
    erode_climate(tectonic, radius_m, p, runoff, None)
}

/// As [`erode_with_runoff`], plus lakes and closed basins from the water balance:
/// `lakes` = (runoff, open-water net evaporation E − P), both mm/yr. Lakes are base level
/// in the solve (see `solve_full`); the committed height is the rock/lake bed and
/// `lake_depth_m` is the water depth (level − bed).
pub fn erode_climate(
    tectonic: &TectonicSnapshot,
    radius_m: f64,
    p: &ErosionParams,
    runoff: Option<&CubeField<f32>>,
    lakes: Option<(&CubeField<f32>, &CubeField<f32>)>,
) -> ErosionResult {
    erode_climate_with_progress(tectonic, radius_m, p, runoff, lakes, &|_, _| true)
        .expect("erosion without a cancelling callback always finishes")
}

/// As [`erode_climate`], reporting solver steps through `on_step(k, n)`; `None` if it
/// returned false (cancelled; see [`solve_full_with_progress`]).
pub fn erode_climate_with_progress(
    tectonic: &TectonicSnapshot,
    radius_m: f64,
    p: &ErosionParams,
    runoff: Option<&CubeField<f32>>,
    lakes: Option<(&CubeField<f32>, &CubeField<f32>)>,
    on_step: &dyn Fn(u32, u32) -> bool,
) -> Option<ErosionResult> {
    let tectonic_h = &tectonic.height_m;
    let res = if p.res > 0 { p.res.max(16) } else { tectonic_h.res.min(p.max_res.max(16)) };
    let mut grid = Grid::new(res);
    let jitter_offsets = (p.node_jitter > 0.0).then(|| crate::jitter::jitter_grid(&mut grid, p.node_jitter.min(0.4), p.seed));
    let grid = grid;
    let n = grid.n;
    let area_m2: Vec<f64> = grid.weight.iter().map(|w| *w as f64 * 4.0 * std::f64::consts::PI * radius_m * radius_m).collect();
    let weight: Vec<f64> = grid.weight.iter().map(|&w| w as f64).collect();

    // Input heights on the solve grid.
    let mut h0: Vec<f64> = par_map(n, |i| tectonic_h.sample(grid.dir[i]) as f64);
    let ocean: Vec<bool> = h0.iter().map(|&h| h <= 0.0).collect();
    // Initial roughness (land only; the coast stays where tectonics put it).
    let tectonic_relief = h0.clone();
    if p.initial_roughness_m > 0.0 {
        for i in (0..n).filter(|&i| !ocean[i]) {
            // Scaled by local relief: on a plain flatter than the noise, the step-0 fill would
            // turn it into ε-flats whose Priority-Flood routing is straight (s01 review).
            let relief = grid.neighbors(i).iter().map(|&j| (tectonic_relief[i] - tectonic_relief[j as usize]).abs()).fold(0.0, f64::max);
            let amp = p.initial_roughness_m.min(0.25 * relief);
            let r = if p.roughness_wavelength_km > 0.0 {
                // Position-keyed, so every resolution starts from the same large-scale bumps.
                let cell_m = std::f64::consts::FRAC_PI_2 * radius_m / res as f64;
                let salt = (p.seed as u32) ^ 0x5EED_FA11;
                crate::jitter::roughness_at(grid.dir[i], radius_m, p.roughness_wavelength_km * 1e3, 2.0 * cell_m, salt) * amp
            } else {
                (unit_hash(p.seed ^ 0x0005_EED0_FA11, i as u64) - 0.5) * 2.0 * amp
            };
            h0[i] = (h0[i] + r).max(h0[i].min(1.0));
        }
    }
    let k_field = lithology(&grid, radius_m, p);
    let runoff_w: Option<Vec<f64>> = runoff.map(|f| par_map(n, |i| f.sample(grid.dir[i]).clamp(0.02, 6.0) as f64));
    let lake_fields: Option<(Vec<f64>, Vec<f64>)> = lakes.map(|(ro, ev)| {
        (par_map(n, |i| ro.sample(grid.dir[i]) as f64 * 1e-3), par_map(n, |i| ev.sample(grid.dir[i]) as f64 * 1e-3))
    });
    let lake_in = lake_fields.as_ref().map(|(r, e)| LakeInputs { runoff_m_yr: r, net_evap_m_yr: e });

    let mut land_h: Vec<f64> = (0..n).filter(|&i| !ocean[i]).map(|i| h0[i]).collect();
    land_h.sort_by(f64::total_cmp);
    let p99 = percentile(&land_h, 0.99).max(1.0);

    // Uplift from tectonic process (or the elevation fallback).
    let (uplift, uplift_used) = uplift::compute(&grid, tectonic, radius_m, &h0, &ocean, p);

    let mut h = h0.clone();
    let mut net = Network::new(n);
    // Jitter already breaks the 8-direction bias (and Rho8's diagonal test assumes a lattice).
    net.rho8_seed = (p.rho8 && p.node_jitter <= 0.0).then_some(p.seed);
    let mut filled_steps = 0u32;
    let mut budget = MassBudget::default();
    if p.enabled && !land_h.is_empty() {
        let lim = (p.elevation_limit_m > 0.0).then_some(p.elevation_limit_m);
        let solved = solve_full_with_progress(&grid, radius_m, &mut h, &ocean, &uplift, &k_field, runoff_w.as_deref(), lim, lake_in.as_ref(), p, on_step)?;
        net = solved.net;
        budget = solved.budget;
        filled_steps = solved.filled_steps;
    }
    let start_land_m3: f64 = (0..n).filter(|&i| !ocean[i]).map(|i| h0[i] * area_m2[i]).sum();
    let end_land_m3: f64 = (0..n).filter(|&i| !ocean[i]).map(|i| h[i] * area_m2[i]).sum();
    let mass_balance_error = (end_land_m3 - start_land_m3 - budget.net_m3()).abs() / start_land_m3.abs().max(1.0);

    // Final routing on the result. Without lakes: measure depressions, then fill them for
    // drainage. With lakes: keep the bed, route on a filled copy to the sea or closed lakes.
    let wet: Vec<f64> = match &runoff_w {
        Some(r) => area_m2.iter().zip(r).map(|(a, r)| a * r).collect(),
        None => area_m2.clone(),
    };
    let mut endorheic_frac = 0.0;
    let mut lake_frac = 0.0;
    let mut lake: Vec<f64> = match &lake_in {
        Some(li) => {
            let lk = crate::lakes::compute(&grid, &h, &ocean, li.runoff_m_yr, li.net_evap_m_yr, &area_m2);
            let sub: Vec<bool> = (0..n).map(|i| !ocean[i] && lk.level[i] > h[i]).collect();
            let sink: Vec<bool> = (0..n).map(|i| ocean[i] || lk.closed_lake[i]).collect();
            let base_h: Vec<f64> = (0..n).map(|i| if sub[i] { lk.level[i] } else { 0.0 }).collect();
            let mut f: Vec<f64> = (0..n).map(|i| if sub[i] { lk.level[i] } else { h[i] }).collect();
            priority_flood_fill_to(&grid, &mut f, &sink, Some(&base_h));
            net.route(&grid, &f, &sink);
            net.build_stack();
            net.accumulate(&wet, &sink);
            let land: f64 = (0..n).filter(|&i| !ocean[i]).map(|i| area_m2[i]).sum::<f64>().max(1.0);
            endorheic_frac = (0..n).filter(|&i| !ocean[i] && lk.endorheic[i]).map(|i| area_m2[i]).sum::<f64>() / land;
            lake_frac = (0..n).filter(|&i| sub[i]).map(|i| area_m2[i]).sum::<f64>() / land;
            (0..n).map(|i| if sub[i] { lk.level[i] - h[i] } else { 0.0 }).collect()
        }
        None => {
            let pre_fill = h.clone();
            priority_flood_fill(&grid, &mut h, &ocean);
            net.route(&grid, &h, &ocean);
            net.build_stack();
            net.accumulate(&wet, &ocean);
            (0..n).map(|i| if ocean[i] { 0.0 } else { (h[i] - pre_fill[i]).max(0.0) }).collect()
        }
    };

    let uplift_eff: Vec<f64> = (0..n)
        .map(|i| if p.elevation_limit_m > 0.0 { uplift[i] * (1.0 - h[i].max(0.0) / p.elevation_limit_m).max(0.0) } else { uplift[i] })
        .collect();

    // Rescale land so the 99th percentile matches tectonics (see crate docs: calibration).
    let mut eroded_land: Vec<f64> = (0..n).filter(|&i| !ocean[i]).map(|i| h[i]).collect();
    eroded_land.sort_by(f64::total_cmp);
    let scale = if eroded_land.is_empty() { 1.0 } else { p99 / percentile(&eroded_land, 0.99).max(1.0) };
    for i in 0..n {
        if !ocean[i] {
            h[i] *= scale;
            lake[i] *= scale;
        }
    }

    // Geomorphic stats on the physics output (the remap below is monotonic, so drainage
    // is untouched, but slope magnitudes change).
    let mut stats = stats::measure(&grid, &h, &ocean, &net, &area_m2, &uplift, &k_field, radius_m, p, filled_steps, scale);
    stats.uplift_from_tectonics = uplift_used == UpliftSource::Tectonic;
    stats.mass_balance_error = mass_balance_error;
    stats.fill_to_erosion = budget.filled_m3 / budget.eroded_m3.max(1.0);
    stats.endorheic_fraction = endorheic_frac;
    stats.lake_fraction = lake_frac;

    let mut height_scale = None;
    if p.linear_mean_m > 0.0 {
        // One factor: undo the p99 rescale's arbitrariness and hit the target mean.
        let (mut sum, mut wsum) = (0.0, 0.0);
        for i in (0..n).filter(|&i| !ocean[i]) {
            sum += h[i] * weight[i];
            wsum += weight[i];
        }
        let mean = sum / wsum.max(1e-12);
        let target = p.linear_mean_m * p.hypsometry_scale.max(1e-6);
        let c = target / mean.max(1e-6);
        for i in (0..n).filter(|&i| !ocean[i]) {
            h[i] = (h[i] * c).max(1.0e-3);
            lake[i] *= c;
        }
        stats.mean_land_m = target;
        height_scale = Some(scale * c);
    } else if p.hypsometry_scale > 0.0 {
        hypsometry::remap(&mut h, &ocean, &weight, p.hypsometry_scale);
        let (mut sum, mut wsum) = (0.0, 0.0);
        for i in (0..n).filter(|&i| !ocean[i]) {
            sum += h[i] * weight[i];
            wsum += weight[i];
        }
        stats.mean_land_m = sum / wsum.max(1e-12);
    }

    // Back to CubeFields. If the solve ran coarser than the input, upsample the *change*.
    let to_field = |vals: &[f64]| -> CubeField<f32> {
        if let Some(offsets) = &jitter_offsets {
            return crate::jitter::resample(&grid, offsets, vals);
        }
        let mut f = CubeField::<f32>::new(res);
        for (k, v) in f.data.iter_mut().enumerate() {
            *v = vals[grid.sample_node[k] as usize] as f32;
        }
        f
    };
    let delta: Vec<f64> = (0..n).map(|i| h[i] - h0[i]).collect();
    let delta_f = to_field(&delta);
    let height_m = if res >= tectonic_h.res {
        let mut f = to_field(&h);
        if jitter_offsets.is_some() {
            // Interpolated between jittered nodes, the raster can hold small pits: fill them
            // on the regular grid so the committed field stays drained.
            let reg = Grid::new(res);
            let mut hr: Vec<f64> = (0..reg.n).map(|i| f.data[reg.node_sample[i] as usize] as f64).collect();
            let reg_ocean: Vec<bool> = hr.iter().map(|&v| v <= 0.0).collect();
            priority_flood_fill(&reg, &mut hr, &reg_ocean);
            for (k, v) in f.data.iter_mut().enumerate() {
                *v = hr[reg.sample_node[k] as usize] as f32;
            }
        }
        f
    } else {
        let mut out = tectonic_h.clone();
        for (k, d) in out.iter_samples().map(|(k, _, _, _, d)| (k, d)).collect::<Vec<_>>() {
            let t = tectonic_h.data[k];
            if t > 0.0 {
                // Keep land land: the coast is the committed tectonic contour.
                out.data[k] = (t + delta_f.sample(d)).max(0.5 * t.min(1.0)).max(0.1);
            }
        }
        // The upsampled change is not drained on the finer grid: fill it there, so the
        // committed field is pit-free at every resolution (s01 review).
        let fine = Grid::new(out.res);
        let mut hf: Vec<f64> = (0..fine.n).map(|i| out.data[fine.node_sample[i] as usize] as f64).collect();
        let fine_ocean: Vec<bool> = hf.iter().map(|&v| v <= 0.0).collect();
        priority_flood_fill(&fine, &mut hf, &fine_ocean);
        for (k, v) in out.data.iter_mut().enumerate() {
            *v = hf[fine.sample_node[k] as usize] as f32;
        }
        out
    };
    let drainage_km2 = to_field(&net.area.iter().map(|a| a * 1.0e-6).collect::<Vec<_>>());
    let lake_depth_m = to_field(&lake);
    // Committed uplift is the *effective* rate at the final surface, U·(1 − h/limit): at
    // steady state it equals the erosion rate, which sets the detail's hillslope regime
    // (plateaus at the strength limit erode slowly and stay gentle).
    let uplift_mm_yr = to_field(&uplift_eff.iter().map(|u| u * 1000.0).collect::<Vec<_>>());
    // Flow direction per node (unit tangent towards the receiver), then onto the raster.
    let node_flow: Vec<glam::DVec3> = (0..n)
        .map(|i| {
            let r = net.receiver[i] as usize;
            if ocean[i] || r == i {
                return glam::DVec3::ZERO;
            }
            let d = grid.dir[r] - grid.dir[i];
            (d - grid.dir[i] * d.dot(grid.dir[i])).normalize_or_zero()
        })
        .collect();
    let mut flow_dir = CubeField::<[f32; 3]>::new(res);
    if jitter_offsets.is_some() {
        let comp = |c: usize| to_field(&node_flow.iter().map(|v| v.to_array()[c]).collect::<Vec<_>>());
        let (fx, fy, fz) = (comp(0), comp(1), comp(2));
        for (k, v) in flow_dir.data.iter_mut().enumerate() {
            let t = glam::Vec3::new(fx.data[k], fy.data[k], fz.data[k]).normalize_or_zero();
            *v = t.to_array();
        }
    } else {
        for (k, v) in flow_dir.data.iter_mut().enumerate() {
            *v = node_flow[grid.sample_node[k] as usize].as_vec3().to_array();
        }
    }

    let min_area = p.river_min_area_km2 * 1.0e6;
    let rivers = (0..n)
        .filter(|&i| !ocean[i] && net.area[i] >= min_area)
        .map(|i| RiverSegment {
            from: grid.dir[i],
            to: grid.dir[net.receiver[i] as usize],
            area_km2: (net.area[i] * 1.0e-6) as f32,
        })
        .collect();

    Some(ErosionResult { height_m, drainage_km2, flow_dir, lake_depth_m, uplift_mm_yr, height_scale, rivers, stats })
}

/// The drainage network: receivers, the topological stack and drainage areas.
pub struct Network {
    /// Receiver of each node (itself for ocean nodes).
    pub receiver: Vec<u32>,
    /// Topological order, coast first.
    pub stack: Vec<u32>,
    /// Drainage area including the node's own cell, m² (0 at sea).
    pub area: Vec<f64>,
    /// `Some(seed)` routes with Rho8 (see [`ErosionParams::rho8`]); `None` is plain D8.
    pub rho8_seed: Option<u64>,
    /// Per node, bit k set if adjacency slot k is a diagonal neighbour (see
    /// [`diagonal_mask`]). Filled on the first Rho8 route.
    diag: Vec<u16>,
    donors_start: Vec<u32>,
    donors: Vec<u32>,
}

impl Network {
    pub fn new(n: usize) -> Self {
        Self {
            receiver: (0..n as u32).collect(),
            stack: Vec::with_capacity(n),
            area: vec![0.0; n],
            rho8_seed: None,
            diag: Vec::new(),
            donors_start: vec![0; n + 1],
            donors: vec![0; n],
        }
    }

    /// Steepest-descent receivers. Returns the number of land pits (no lower neighbour).
    /// Ties go to the lower node index, so routing is deterministic.
    pub fn route(&mut self, grid: &Grid, h: &[f64], ocean: &[bool]) -> usize {
        let rho8 = self.rho8_seed;
        if rho8.is_some() && self.diag.len() != grid.n {
            self.diag = diagonal_mask(grid);
        }
        let diag = &self.diag;
        let rec: Vec<u32> = par_map(grid.n, |i| {
            if ocean[i] {
                return i as u32;
            }
            // Rho8: one random number per node, fixed for the whole run.
            let diag_len = rho8.map(|seed| 2.0 - unit_hash(seed, i as u64));
            let mut best = i as u32;
            let mut best_s = 0.0;
            for (slot, &j) in grid.neighbors(i).iter().enumerate() {
                let j = j as usize;
                let hj = if ocean[j] { 0.0f64.min(h[j]) } else { h[j] };
                let mut d = (grid.dir[i] - grid.dir[j]).length();
                // Diagonal neighbours are ~√2 axis spacings away; Rho8 replaces √2 with a
                // random value in [1, 2].
                if let Some(f) = diag_len
                    && diag[i] >> slot & 1 == 1
                {
                    d *= f / std::f64::consts::SQRT_2;
                }
                let s = (h[i] - hj) / d;
                if s > best_s || (s == best_s && s > 0.0 && (j as u32) < best) {
                    best_s = s;
                    best = j as u32;
                }
            }
            best
        });
        self.receiver = rec;
        (0..grid.n).filter(|&i| !ocean[i] && self.receiver[i] as usize == i).count()
    }

    /// Topological order: base-level (ocean) nodes first, then every node after its receiver.
    pub fn build_stack(&mut self) {
        let n = self.receiver.len();
        self.donors_start.iter_mut().for_each(|c| *c = 0);
        for i in 0..n {
            let r = self.receiver[i] as usize;
            if r != i {
                self.donors_start[r + 1] += 1;
            }
        }
        for i in 0..n {
            self.donors_start[i + 1] += self.donors_start[i];
        }
        let mut fill = self.donors_start.clone();
        for i in 0..n {
            let r = self.receiver[i] as usize;
            if r != i {
                self.donors[fill[r] as usize] = i as u32;
                fill[r] += 1;
            }
        }
        self.stack.clear();
        let mut queue: VecDeque<u32> = (0..n as u32).filter(|&i| self.receiver[i as usize] == i).collect();
        while let Some(i) = queue.pop_front() {
            self.stack.push(i);
            let (a, b) = (self.donors_start[i as usize] as usize, self.donors_start[i as usize + 1] as usize);
            queue.extend(&self.donors[a..b]);
        }
    }

    /// Drainage area: each land node's own area plus everything upstream of it.
    pub fn accumulate(&mut self, cell_area: &[f64], ocean: &[bool]) {
        for (i, a) in self.area.iter_mut().enumerate() {
            *a = if ocean[i] { 0.0 } else { cell_area[i] };
        }
        for &i in self.stack.iter().rev() {
            let i = i as usize;
            let r = self.receiver[i] as usize;
            if r != i && !ocean[r] {
                self.area[r] += self.area[i];
            }
        }
    }
}

/// Priority-Flood+ε (Barnes et al. 2014): raise every closed depression on land to its
/// spill height plus a tiny gradient, so each land node has a strictly lower path to the sea.
pub fn priority_flood_fill(grid: &Grid, h: &mut [f64], ocean: &[bool]) {
    priority_flood_fill_to(grid, h, ocean, None)
}

/// As [`priority_flood_fill`], with outlets seeded at their own `base` height instead of 0
/// (closed lakes are outlets at their water level).
pub fn priority_flood_fill_to(grid: &Grid, h: &mut [f64], ocean: &[bool], base: Option<&[f64]>) {
    #[derive(PartialEq)]
    struct Item(f64, u32);
    impl Eq for Item {}
    impl Ord for Item {
        fn cmp(&self, o: &Self) -> Ordering {
            // Min-heap on height, then on index (determinism).
            o.0.total_cmp(&self.0).then_with(|| o.1.cmp(&self.1))
        }
    }
    impl PartialOrd for Item {
        fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
            Some(self.cmp(o))
        }
    }
    let n = grid.n;
    let mut done = ocean.to_vec();
    let mut heap = BinaryHeap::new();
    // Seeds: ocean nodes touching land, at base level 0.
    for i in 0..n {
        if ocean[i] && grid.neighbors(i).iter().any(|&j| !ocean[j as usize]) {
            heap.push(Item(base.map_or(0.0, |b| b[i]), i as u32));
        }
    }
    // A world with no ocean: drain to the lowest node.
    if heap.is_empty() {
        let lowest = (0..n).min_by(|&a, &b| h[a].total_cmp(&h[b])).unwrap_or(0);
        done[lowest] = true;
        heap.push(Item(h[lowest], lowest as u32));
    }
    // ε per metre of arc: ~1 mm per km, too small to see, enough to route flats.
    let eps_per_rad = 1.0e-3 * 6_371.0;
    while let Some(Item(hc, c)) = heap.pop() {
        for &j in grid.neighbors(c as usize) {
            let j = j as usize;
            if done[j] {
                continue;
            }
            done[j] = true;
            let d = (grid.dir[c as usize] - grid.dir[j]).length();
            let min_h = hc + eps_per_rad * d;
            if h[j] < min_h {
                h[j] = min_h;
            }
            heap.push(Item(h[j], j as u32));
        }
    }
}


fn percentile(sorted: &[f64], q: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    sorted[((sorted.len() - 1) as f64 * q).round() as usize]
}
