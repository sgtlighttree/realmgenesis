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
    /// Smallest drainage area drawn as a river, km².
    pub river_min_area_km2: f64,
    /// Hypsometric calibration: land heights are quantile-mapped onto the land hypsometric
    /// curve (`hypsometry.rs`) scaled by this factor (1.0 = Earth-like; 2.0 = everything
    /// twice as high). 0 disables it (land is then linearly rescaled to the tectonic 99th
    /// percentile).
    pub hypsometry_scale: f64,
}

impl Default for ErosionParams {
    fn default() -> Self {
        Self {
            enabled: true,
            max_res: 512,
            duration_myr: 30.0,
            steps: 50,
            k: 1.0e-5,
            m: 0.45,
            uplift_max_mm_yr: 2.5,
            uplift_base_mm_yr: 0.15,
            uplift_exponent: 0.8,
            uplift_source: UpliftSource::Tectonic,
            belt_width_km: 400.0,
            oceanic_crust_uplift: 0.5,
            scale_by_convergence_rate: true,
            river_min_area_km2: 20_000.0,
            hypsometry_scale: 1.0,
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
    /// Eroded heights at the input field's resolution (m above sea level).
    pub height_m: CubeField<f32>,
    /// Upstream drainage area, km², at the erosion grid's resolution (0 at sea).
    pub drainage_km2: CubeField<f32>,
    /// Downstream unit flow direction (tangent), xyz, at the erosion grid's resolution.
    /// Zero at sea and on flats.
    pub flow_dir: CubeField<[f32; 3]>,
    /// Depth of standing water in closed basins at the end (lakes), m.
    pub lake_depth_m: CubeField<f32>,
    pub rivers: Vec<RiverSegment>,
    pub stats: ErosionStats,
}

// ---------------------------------------------------------------------------------------
// The numerical core: arrays in, arrays out. `erode()` and the tests both call these.
// ---------------------------------------------------------------------------------------

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
    let mut b = MassBudget::default();
    for &i in c.stack {
        let i = i as usize;
        if c.base[i] {
            continue;
        }
        let r = c.receiver[i] as usize;
        let hr = if c.base[r] { 0.0 } else { h[r] };
        let f = k * dt * c.drainage_m2[i].powf(m) / c.length_m[i].max(1.0);
        let new = (h[i] + c.uplift_m_yr[i] * dt + f * hr) / (1.0 + f);
        b.uplift_m3 += c.uplift_m_yr[i] * dt * c.cell_area_m2[i];
        b.eroded_m3 += f * (new - hr) * c.cell_area_m2[i];
        h[i] = new;
    }
    b
}

/// What [`solve`] leaves behind.
pub struct Solved {
    /// The drainage network of the last step (receivers, stack, areas).
    pub net: Network,
    /// Volumes over the whole run, so `Σ_land area·(h_end − h_start) = budget.net_m3()`.
    pub budget: MassBudget,
    /// How many steps needed a Priority-Flood fill.
    pub filled_steps: u32,
}

/// Run `p.steps` implicit steps of duration `p.duration_myr / p.steps` on `h` (metres, with
/// `ocean` nodes as base level), re-routing the network each step.
pub fn solve(grid: &Grid, radius_m: f64, h: &mut [f64], ocean: &[bool], uplift_m_yr: &[f64], p: &ErosionParams) -> Solved {
    let n = grid.n;
    let cell_area: Vec<f64> = grid.weight.iter().map(|w| *w as f64 * 4.0 * std::f64::consts::PI * radius_m * radius_m).collect();
    let mut net = Network::new(n);
    let mut budget = MassBudget::default();
    let mut filled_steps = 0u32;
    let steps = p.steps.max(1);
    let dt = p.duration_myr * 1.0e6 / steps as f64;
    let mut length = vec![0.0f64; n];
    for step in 0..steps {
        // 1. Receivers (fill pits first if there are any, and always on step 0).
        if step == 0 || net.route(grid, h, ocean) > 0 {
            let before = h.to_vec();
            priority_flood_fill(grid, h, ocean);
            budget.filled_m3 += (0..n).filter(|&i| !ocean[i]).map(|i| (h[i] - before[i]) * cell_area[i]).sum::<f64>();
            filled_steps += 1;
            let pits = net.route(grid, h, ocean);
            debug_assert_eq!(pits, 0, "filled surface must drain everywhere");
        }
        // 2–3. Stack and drainage area.
        net.build_stack();
        net.accumulate(&cell_area, ocean);
        // 4. Implicit stream-power incision with uplift.
        for (i, l) in length.iter_mut().enumerate() {
            *l = (grid.dir[i] - grid.dir[net.receiver[i] as usize]).length() * radius_m;
        }
        let col = Column {
            receiver: &net.receiver,
            stack: &net.stack,
            drainage_m2: &net.area,
            length_m: &length,
            uplift_m_yr,
            base: ocean,
            cell_area_m2: &cell_area,
        };
        budget.add(implicit_step(&col, h, p.k, p.m, dt));
    }
    Solved { net, budget, filled_steps }
}

/// Erode the tectonic snapshot. Deterministic: same input, same output, on any thread count.
///
/// Heights come from `tectonic.height_m`; uplift comes from the snapshot's boundaries, crust
/// and plate motions (see [`UpliftSource`]).
pub fn erode(tectonic: &TectonicSnapshot, radius_m: f64, p: &ErosionParams) -> ErosionResult {
    let tectonic_h = &tectonic.height_m;
    let res = tectonic_h.res.min(p.max_res.max(16));
    let grid = Grid::new(res);
    let n = grid.n;
    let area_m2: Vec<f64> = grid.weight.iter().map(|w| *w as f64 * 4.0 * std::f64::consts::PI * radius_m * radius_m).collect();
    let weight: Vec<f64> = grid.weight.iter().map(|&w| w as f64).collect();

    // Input heights on the solve grid.
    let h0: Vec<f64> = par_map(n, |i| tectonic_h.sample(grid.dir[i]) as f64);
    let ocean: Vec<bool> = h0.iter().map(|&h| h <= 0.0).collect();

    let mut land_h: Vec<f64> = (0..n).filter(|&i| !ocean[i]).map(|i| h0[i]).collect();
    land_h.sort_by(f64::total_cmp);
    let p99 = percentile(&land_h, 0.99).max(1.0);

    // Uplift from tectonic process (or the elevation fallback).
    let (uplift, uplift_used) = uplift::compute(&grid, tectonic, radius_m, &h0, &ocean, p);

    let mut h = h0.clone();
    let mut net = Network::new(n);
    let mut filled_steps = 0u32;
    let mut budget = MassBudget::default();
    if p.enabled && !land_h.is_empty() {
        let solved = solve(&grid, radius_m, &mut h, &ocean, &uplift, p);
        net = solved.net;
        budget = solved.budget;
        filled_steps = solved.filled_steps;
    }
    let start_land_m3: f64 = (0..n).filter(|&i| !ocean[i]).map(|i| h0[i] * area_m2[i]).sum();
    let end_land_m3: f64 = (0..n).filter(|&i| !ocean[i]).map(|i| h[i] * area_m2[i]).sum();
    let mass_balance_error = (end_land_m3 - start_land_m3 - budget.net_m3()).abs() / start_land_m3.abs().max(1.0);

    // Final routing on the result, measuring lakes before filling them for drainage.
    let pre_fill = h.clone();
    priority_flood_fill(&grid, &mut h, &ocean);
    net.route(&grid, &h, &ocean);
    net.build_stack();
    net.accumulate(&area_m2, &ocean);
    let lake: Vec<f64> = (0..n).map(|i| if ocean[i] { 0.0 } else { (h[i] - pre_fill[i]).max(0.0) }).collect();

    // Rescale land so the 99th percentile matches tectonics (see crate docs: calibration).
    let mut eroded_land: Vec<f64> = (0..n).filter(|&i| !ocean[i]).map(|i| h[i]).collect();
    eroded_land.sort_by(f64::total_cmp);
    let scale = if eroded_land.is_empty() { 1.0 } else { p99 / percentile(&eroded_land, 0.99).max(1.0) };
    for i in 0..n {
        if !ocean[i] {
            h[i] *= scale;
        }
    }

    // Geomorphic stats on the physics output (the remap below is monotonic, so drainage
    // is untouched, but slope magnitudes change).
    let mut stats = stats::measure(&grid, &h, &ocean, &net, &area_m2, &uplift, radius_m, p, filled_steps, scale);
    stats.uplift_from_tectonics = uplift_used == UpliftSource::Tectonic;
    stats.mass_balance_error = mass_balance_error;

    if p.hypsometry_scale > 0.0 {
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
        let mut f = CubeField::<f32>::new(res);
        for (k, v) in f.data.iter_mut().enumerate() {
            *v = vals[grid.sample_node[k] as usize] as f32;
        }
        f
    };
    let delta: Vec<f64> = (0..n).map(|i| h[i] - h0[i]).collect();
    let delta_f = to_field(&delta);
    let height_m = if res == tectonic_h.res {
        to_field(&h)
    } else {
        let mut out = tectonic_h.clone();
        for (k, d) in out.iter_samples().map(|(k, _, _, _, d)| (k, d)).collect::<Vec<_>>() {
            let t = tectonic_h.data[k];
            if t > 0.0 {
                // Keep land land: the coast is the committed tectonic contour.
                out.data[k] = (t + delta_f.sample(d)).max(0.5 * t.min(1.0)).max(0.1);
            }
        }
        out
    };
    let drainage_km2 = to_field(&net.area.iter().map(|a| a * 1.0e-6).collect::<Vec<_>>());
    let lake_depth_m = to_field(&lake);
    let mut flow_dir = CubeField::<[f32; 3]>::new(res);
    for (k, v) in flow_dir.data.iter_mut().enumerate() {
        let i = grid.sample_node[k] as usize;
        let r = net.receiver[i] as usize;
        if !ocean[i] && r != i {
            let d = grid.dir[r] - grid.dir[i];
            let t = (d - grid.dir[i] * d.dot(grid.dir[i])).normalize_or_zero();
            *v = [t.x as f32, t.y as f32, t.z as f32];
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

    ErosionResult { height_m, drainage_km2, flow_dir, lake_depth_m, rivers, stats }
}

/// The drainage network: receivers, the topological stack and drainage areas.
pub struct Network {
    /// Receiver of each node (itself for ocean nodes).
    pub receiver: Vec<u32>,
    /// Topological order, coast first.
    pub stack: Vec<u32>,
    /// Drainage area including the node's own cell, m² (0 at sea).
    pub area: Vec<f64>,
    donors_start: Vec<u32>,
    donors: Vec<u32>,
}

impl Network {
    pub fn new(n: usize) -> Self {
        Self {
            receiver: (0..n as u32).collect(),
            stack: Vec::with_capacity(n),
            area: vec![0.0; n],
            donors_start: vec![0; n + 1],
            donors: vec![0; n],
        }
    }

    /// Steepest-descent receivers. Returns the number of land pits (no lower neighbour).
    /// Ties go to the lower node index, so routing is deterministic.
    pub fn route(&mut self, grid: &Grid, h: &[f64], ocean: &[bool]) -> usize {
        let rec: Vec<u32> = par_map(grid.n, |i| {
            if ocean[i] {
                return i as u32;
            }
            let mut best = i as u32;
            let mut best_s = 0.0;
            for &j in grid.neighbors(i) {
                let j = j as usize;
                let hj = if ocean[j] { 0.0f64.min(h[j]) } else { h[j] };
                let d = (grid.dir[i] - grid.dir[j]).length();
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
            heap.push(Item(0.0, i as u32));
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
