//! Solver tests that check the numerics against things the solver does not assume itself:
//! conservation laws, a closed-form 1-D solution, and convergence in the time step.
//!
//! (The "uplift-normalised concavity == m" statistic in `stats.rs` is nearly an algebraic
//! identity at the discrete steady state, so it is *not* relied on here.)

use rg_erosion::{implicit_step, solve, Column, ErosionParams, MassBudget, Network};
use rg_tectonics::grid::Grid;

const R: f64 = 6_371_000.0;

/// A small real world: tectonic heights on a Grid, with a deterministic, spatially varying
/// uplift field (so uplift is not trivially uniform).
struct World {
    grid: Grid,
    h0: Vec<f64>,
    ocean: Vec<bool>,
    uplift: Vec<f64>,
    cell_area: Vec<f64>,
}

fn world(seed: u64, res: u32) -> World {
    let tp = rg_tectonics::TectonicsParams { seed, res, ..Default::default() };
    let snap = rg_tectonics::generate(&tp);
    let grid = Grid::new(res);
    let h0: Vec<f64> = (0..grid.n).map(|i| snap.height_m.sample(grid.dir[i]) as f64).collect();
    let ocean: Vec<bool> = h0.iter().map(|&h| h <= 0.0).collect();
    let hmax = h0.iter().cloned().fold(1.0, f64::max);
    let uplift: Vec<f64> = (0..grid.n)
        .map(|i| if ocean[i] { 0.0 } else { (0.1 + 1.5 * (h0[i] / hmax).powf(0.7) + 0.2 * ((i % 7) as f64 / 7.0)) * 1e-3 })
        .collect();
    let cell_area = grid.weight.iter().map(|&w| w as f64 * 4.0 * std::f64::consts::PI * R * R).collect();
    World { grid, h0, ocean, uplift, cell_area }
}

fn params(steps: u32, duration_myr: f64) -> ErosionParams {
    ErosionParams { steps, duration_myr, ..Default::default() }
}

fn land_volume_change(w: &World, h_end: &[f64]) -> f64 {
    (0..w.grid.n).filter(|&i| !w.ocean[i]).map(|i| (h_end[i] - w.h0[i]) * w.cell_area[i]).sum()
}

// ---------------------------------------------------------------------------------------
// (a) Area conservation
// ---------------------------------------------------------------------------------------

fn check_area_conservation(w: &World, net: &Network, what: &str) {
    let total_land: f64 = (0..w.grid.n).filter(|&i| !w.ocean[i]).map(|i| w.cell_area[i]).sum();
    // Outlets: land nodes whose receiver is the ocean.
    let outlets: Vec<usize> = (0..w.grid.n).filter(|&i| !w.ocean[i] && w.ocean[net.receiver[i] as usize]).collect();
    assert!(!outlets.is_empty());
    let at_outlets: f64 = outlets.iter().map(|&i| net.area[i]).sum();
    assert!(
        ((at_outlets - total_land) / total_land).abs() < 1e-9,
        "{what}: outlet area {at_outlets} vs land {total_land}"
    );
    // Independent of the accumulation order: walk every land node's receiver chain to its
    // outlet and add its own cell there; each outlet's area must match.
    let mut by_outlet = vec![0.0f64; w.grid.n];
    for i in (0..w.grid.n).filter(|&i| !w.ocean[i]) {
        let mut c = i;
        let mut hops = 0;
        while !w.ocean[net.receiver[c] as usize] {
            assert_ne!(net.receiver[c] as usize, c, "{what}: node {c} is a pit");
            c = net.receiver[c] as usize;
            hops += 1;
            assert!(hops <= w.grid.n, "{what}: receiver cycle");
        }
        by_outlet[c] += w.cell_area[i];
    }
    for &o in &outlets {
        assert!(
            ((net.area[o] - by_outlet[o]) / by_outlet[o]).abs() < 1e-9,
            "{what}: outlet {o}: accumulated {} vs walked {}",
            net.area[o],
            by_outlet[o]
        );
    }
}

#[test]
fn drainage_area_is_conserved_at_the_outlets() {
    for seed in [3, 11] {
        let w = world(seed, 48);
        // Fresh network on the raw (filled) surface.
        let mut h = w.h0.clone();
        rg_erosion::priority_flood_fill(&w.grid, &mut h, &w.ocean);
        let mut net = Network::new(w.grid.n);
        assert_eq!(net.route(&w.grid, &h, &w.ocean), 0);
        net.build_stack();
        net.accumulate(&w.cell_area, &w.ocean);
        check_area_conservation(&w, &net, &format!("seed {seed}, initial"));
        // And the network the solver itself ends on, after 30 Myr of re-routing.
        let mut h = w.h0.clone();
        let solved = solve(&w.grid, R, &mut h, &w.ocean, &w.uplift, &params(30, 30.0));
        check_area_conservation(&w, &solved.net, &format!("seed {seed}, after solve"));
    }
}

// ---------------------------------------------------------------------------------------
// (b) Mass balance
// ---------------------------------------------------------------------------------------

#[test]
fn mass_balance_closes_over_a_run() {
    for seed in [3, 11] {
        let w = world(seed, 48);
        let p = params(40, 20.0);
        let mut h = w.h0.clone();
        let s = solve(&w.grid, R, &mut h, &w.ocean, &w.uplift, &p);
        let dv = land_volume_change(&w, &h);
        let b = s.budget;
        let scale = b.uplift_m3.max(b.eroded_m3);
        // Σ area·(h_end − h_start) = Σ U·T·area − eroded + filled.
        assert!(((dv - b.net_m3()) / scale).abs() < 1e-9, "seed {seed}: ΔV {dv} vs books {} ({b:?})", b.net_m3());
        // The uplift term is independently computable.
        let t = p.duration_myr * 1e6;
        let u_total: f64 = (0..w.grid.n).filter(|&i| !w.ocean[i]).map(|i| w.uplift[i] * t * w.cell_area[i]).sum();
        assert!(((b.uplift_m3 - u_total) / u_total).abs() < 1e-9, "uplift volume {} vs {u_total}", b.uplift_m3);
        // Rivers do remove rock; filling only ever adds.
        assert!(b.eroded_m3 > 0.0 && b.filled_m3 >= 0.0);
    }
}

#[test]
fn eroded_volume_equals_uplift_volume_at_steady_state() {
    // Run a coarse world long enough to reach a dynamic steady state, then measure one
    // more short stretch: rock out must equal rock in (a flux balance, independent of the
    // per-node algebra of the update).
    let w = world(3, 20);
    let mut h = w.h0.clone();
    solve(&w.grid, R, &mut h, &w.ocean, &w.uplift, &params(300, 1500.0));
    let probe = ErosionParams { steps: 1, duration_myr: 0.05, ..ErosionParams::default() };
    let s = solve(&w.grid, R, &mut h, &w.ocean, &w.uplift, &probe);
    let ratio = s.budget.eroded_m3 / s.budget.uplift_m3;
    assert!((ratio - 1.0).abs() < 0.03, "erosion/uplift = {ratio} at supposed steady state");
}

// ---------------------------------------------------------------------------------------
// (c) 1-D analytic checks on a single channel
// ---------------------------------------------------------------------------------------

/// A chain of nodes 0..n: node 0 is base level, node i drains to i−1, x_i = i·dx. The
/// drainage area follows a Hack-type power law, A(x) = a0·((x + x0)/1 km)^p, with a 40 km
/// headwater offset x0 so the integrand 1/A^m has no singularity at the divide.
const CHAIN_A0: f64 = 1.0e6;
const CHAIN_X0: f64 = 40_000.0;
const CHAIN_P: f64 = 1.67;

struct Chain {
    n: usize,
    dx: f64,
    receiver: Vec<u32>,
    stack: Vec<u32>,
    area: Vec<f64>,
    length: Vec<f64>,
    base: Vec<bool>,
    cell: Vec<f64>,
    k: f64,
    m: f64,
}

impl Chain {
    fn new(n: usize, dx: f64, k: f64, m: f64) -> Self {
        let area = (0..n).map(|i| CHAIN_A0 * ((i as f64 * dx + CHAIN_X0) / 1000.0).powf(CHAIN_P)).collect();
        Chain {
            n,
            dx,
            receiver: (0..n).map(|i| i.saturating_sub(1) as u32).collect(),
            stack: (0..n as u32).collect(),
            area,
            length: vec![dx; n],
            base: (0..n).map(|i| i == 0).collect(),
            cell: vec![dx * 1000.0; n],
            k,
            m,
        }
    }
    fn step(&self, h: &mut [f64], uplift: &[f64], dt: f64) -> MassBudget {
        let c = Column {
            receiver: &self.receiver,
            stack: &self.stack,
            drainage_m2: &self.area,
            length_m: &self.length,
            uplift_m_yr: uplift,
            base: &self.base,
            cell_area_m2: &self.cell,
        };
        implicit_step(&c, h, self.k, self.m, dt)
    }
    /// Travel time from the outlet up to node i of a signal moving at K·A^m: Σ dx/(K·A^m).
    fn travel_time(&self, i: usize) -> f64 {
        (1..=i).map(|j| self.dx / (self.k * self.area[j].powf(self.m))).sum()
    }
}

#[test]
fn steady_state_profile_matches_the_analytic_solution() {
    let ch = Chain::new(300, 2000.0, 1e-5, 0.45);
    let u = 1.0e-3;
    let uplift = vec![u; ch.n];
    let mut h = vec![0.0; ch.n];
    // One enormous step goes straight to the steady state of the implicit scheme.
    for _ in 0..3 {
        ch.step(&mut h, &uplift, 1.0e11);
    }
    // Discrete analytic profile: h_i = Σ_{j≤i} U·L/(K·A_j^m).
    let mut exact = 0.0;
    let mut worst: f64 = 0.0;
    for (hi, a) in h.iter().zip(&ch.area).skip(1) {
        exact += u * ch.dx / (ch.k * a.powf(ch.m));
        worst = worst.max(((hi - exact) / exact).abs());
    }
    assert!(worst < 1e-8, "steady profile deviates from Σ U·L/(K·A^m) by {worst}");
    // And against the continuous integral h(x) = ∫ U/(K·A^m) dx', A = a0·((x'+x0)/1 km)^p:
    //   = U/(K·a0^m) · (1 km)^(p·m) · [(x+x0)^(1−pm) − x0^(1−pm)] / (1−pm).
    // The scheme's Riemann sum differs from it by O(dx) (a first-order scheme).
    let q = 1.0 - CHAIN_P * ch.m;
    let x = (ch.n - 1) as f64 * ch.dx;
    let cont = u / (ch.k * CHAIN_A0.powf(ch.m)) * 1000.0f64.powf(CHAIN_P * ch.m) * ((x + CHAIN_X0).powf(q) - CHAIN_X0.powf(q)) / q;
    let rel = (h[ch.n - 1] - cont).abs() / cont;
    assert!(rel < 0.05, "discrete {} vs continuous {cont} (rel {rel})", h[ch.n - 1]);
}

#[test]
fn knickpoint_migrates_upstream_at_k_a_m() {
    let ch = Chain::new(300, 2000.0, 1e-5, 0.45);
    let u1 = 0.5e-3;
    let mut h = vec![0.0; ch.n];
    for _ in 0..3 {
        ch.step(&mut h, &vec![u1; ch.n], 1.0e11);
    }
    let h_old = h.clone();
    // Double the uplift rate everywhere. Characteristics of ∂h/∂t = U − K·A^m·∂h/∂x travel
    // upstream at c = K·A^m, so the profile at time t is h_old + ΔU·min(τ(x), t), where
    // τ(x) = ∫dx/c is the travel time from the outlet.
    let du = u1;
    let uplift = vec![2.0 * u1; ch.n];
    let tau_end = ch.travel_time(ch.n - 1);
    let dt = tau_end / 4000.0;
    let mut t = 0.0;
    let mut checked = 0;
    for step in 1..=2400 {
        ch.step(&mut h, &uplift, dt);
        t += dt;
        if step % 600 == 0 {
            // The knickpoint is where the slope jumps from the new steady value to the old one.
            // Excess slope over the old profile, relative to its fully-adjusted value
            // du·dx/(K·A^m) = du·Δτ: 1 downstream of the front, 0 upstream. A smeared front
            // is symmetric about its true position, so take where the ratio crosses 0.5.
            let ratio = |i: usize| (h[i] - h[i - 1] - (h_old[i] - h_old[i - 1])) / (du * (ch.travel_time(i) - ch.travel_time(i - 1)));
            let reached = (2..ch.n).find(|&i| ratio(i) < 0.5).expect("knickpoint left the domain");
            // Analytically: the first node whose travel time τ exceeds t. Allow one node of slack.
            let expected = (1..ch.n).find(|&i| ch.travel_time(i) > t).expect("expected knickpoint left the domain");
            assert!(
                (reached as i64 - expected as i64).abs() <= 1,
                "t = {t:.3e} yr: knickpoint at node {reached}, K·A^m celerity predicts node {expected}"
            );
            // Whole-profile check against h_old + ΔU·min(τ, t), as a share of the full response.
            let mut sq = 0.0;
            for i in 1..ch.n {
                let want = du * ch.travel_time(i).min(t);
                sq += (h[i] - h_old[i] - want).powi(2);
            }
            let rms = (sq / ch.n as f64).sqrt() / (du * t);
            assert!(rms < 0.02, "t = {t:.3e}: response RMS error {rms}");
            checked += 1;
        }
    }
    assert_eq!(checked, 4);
}

#[test]
fn knickpoint_error_shrinks_with_the_grid_spacing() {
    // Same physical channel (600 km, same A(x), K, U), sampled with dx = 4, 2, 1 km. The
    // front position is first-order accurate in dx, so the response error must fall.
    let err = |n: usize, dx: f64| {
        let ch = Chain::new(n, dx, 1e-5, 0.45);
        let u1 = 0.5e-3;
        let mut h = vec![0.0; n];
        ch.step(&mut h, &vec![u1; n], 1.0e11);
        ch.step(&mut h, &vec![u1; n], 1.0e11);
        let base = h.clone();
        let t_end = 0.5 * ch.travel_time(n - 1);
        let steps = 2000;
        for _ in 0..steps {
            ch.step(&mut h, &vec![2.0 * u1; n], t_end / steps as f64);
        }
        let sq: f64 = (1..n).map(|i| (h[i] - base[i] - u1 * ch.travel_time(i).min(t_end)).powi(2)).sum();
        (sq / n as f64).sqrt() / (u1 * t_end)
    };
    let (e1, e2, e3) = (err(150, 4000.0), err(300, 2000.0), err(600, 1000.0));
    assert!(e2 < 0.8 * e1 && e3 < 0.8 * e2, "errors {e1} {e2} {e3}");
}

#[test]
fn large_time_steps_stay_stable_and_bounded() {
    // Unconditional stability: with a Courant number far above 1 the response must still
    // lie between "no response yet" and the far-field rise ΔU·t (no overshoot, no ringing).
    let ch = Chain::new(300, 2000.0, 1e-5, 0.45);
    let u1 = 0.5e-3;
    let mut h = vec![0.0; ch.n];
    ch.step(&mut h, &vec![u1; ch.n], 1.0e11);
    ch.step(&mut h, &vec![u1; ch.n], 1.0e11);
    let base = h.clone();
    let dt = ch.travel_time(ch.n - 1) / 8.0;
    for k in 1..=8 {
        ch.step(&mut h, &vec![2.0 * u1; ch.n], dt);
        let t = dt * k as f64;
        for i in 1..ch.n {
            let rise = h[i] - base[i];
            assert!(rise >= -1e-9 && rise <= u1 * t * (1.0 + 1e-9), "step {k} node {i}: rise {rise} vs max {}", u1 * t);
        }
    }
}

// ---------------------------------------------------------------------------------------
// (d) Time-step convergence on a real world
// ---------------------------------------------------------------------------------------

#[test]
fn halving_dt_converges() {
    // A short, strongly transient run (1 Myr) on a real world. First-order scheme: the
    // land-mean height and eroded volume must converge as dt halves. Node-by-node heights
    // also converge but with a noise floor: D8 receivers flip and Priority-Flood fills pits
    // discontinuously, so the pointwise RMS is only required to fall overall.
    let w = world(3, 64);
    let steps = [20u32, 40, 80, 160, 320];
    let runs: Vec<(Vec<f64>, MassBudget)> = steps
        .iter()
        .map(|&s| {
            let mut h = w.h0.clone();
            let sol = solve(&w.grid, R, &mut h, &w.ocean, &w.uplift, &params(s, 1.0));
            (h, sol.budget)
        })
        .collect();
    let land: Vec<usize> = (0..w.grid.n).filter(|&i| !w.ocean[i]).collect();
    let mean = |h: &[f64]| land.iter().map(|&i| h[i] * w.cell_area[i]).sum::<f64>() / land.iter().map(|&i| w.cell_area[i]).sum::<f64>();
    let rms = |a: &[f64], b: &[f64]| (land.iter().map(|&i| (a[i] - b[i]).powi(2)).sum::<f64>() / land.len() as f64).sqrt();
    let means: Vec<f64> = runs.iter().map(|r| mean(&r.0)).collect();
    let dm: Vec<f64> = means.windows(2).map(|m| (m[1] - m[0]).abs()).collect();
    let de: Vec<f64> = runs.windows(2).map(|r| (r[1].1.eroded_m3 - r[0].1.eroded_m3).abs() / r[1].1.eroded_m3).collect();
    let dr: Vec<f64> = runs.windows(2).map(|r| rms(&r[0].0, &r[1].0)).collect();
    println!("mean land height {means:?}\nΔmean {dm:?}\nΔeroded (rel) {de:?}\nnodal RMS {dr:?}");
    // The mean converges geometrically (roughly halving) from the second halving on.
    assert!(dm[1] < dm[0] && dm[2] < dm[1] && dm[3] < dm[2], "mean not converging: {dm:?}");
    assert!(dm[3] < 0.5 * dm[0], "mean converges too slowly: {dm:?}");
    assert!(dm[3] < 0.005 * means[4], "final halving moves the mean by {} of {}", dm[3], means[4]);
    assert!(de[3] < de[0], "eroded volume not converging: {de:?}");
    // Pointwise: small relative to the surface, and smaller at the end than at the start.
    assert!(dr[3] < 0.06 * means[4], "final halving moves nodes by {} m RMS (mean {})", dr[3], means[4]);
    assert!(dr[3] < 0.8 * dr[0], "nodal RMS not falling overall: {dr:?}");
}

/// Diagonals are found by topology: every interior node of a face has exactly 4 diagonal
/// and 4 axis neighbours, however warped its cell is.
#[test]
#[allow(clippy::needless_range_loop)] // `a` is a node id used with several arrays
fn diagonal_mask_counts_four_diagonals() {
    let grid = rg_tectonics::grid::Grid::new(24);
    let mask = rg_erosion::diagonal_mask(&grid);
    let mut checked = 0;
    let per = (grid.res + 1) as usize;
    let mut seam_counts = std::collections::BTreeMap::new();
    for a in 0..grid.n {
        let s = grid.node_sample[a] as usize % (per * per);
        let (i, j) = (s % per, s / per);
        let on_edge = i == 0 || j == 0 || i == per - 1 || j == per - 1;
        if on_edge {
            *seam_counts.entry((grid.neighbors(a).len(), mask[a].count_ones())).or_insert(0) += 1;
            continue;
        }
        if grid.neighbors(a).len() == 8 {
            assert_eq!(mask[a].count_ones(), 4, "node {a}");
            // Diagonals are the longer pair in every cell (sanity check on geometry).
            let len = |b: u32| (grid.dir[a] - grid.dir[b as usize]).length();
            let (mut dmin, mut amax) = (f64::MAX, 0f64);
            for (k, &b) in grid.neighbors(a).iter().enumerate() {
                if mask[a] >> k & 1 == 1 { dmin = dmin.min(len(b)) } else { amax = amax.max(len(b)) }
            }
            assert!(dmin > amax * 0.99, "node {a}: diagonal {dmin} vs axis {amax}");
            checked += 1;
        }
    }
    assert!(checked > grid.n / 2);
    // Seam nodes: (neighbours, diagonals) → count. Printed for the record.
    eprintln!("seam nodes: {seam_counts:?}");
}

/// s07 (heightmap layers, advisor's first gate): uplift inverted from a surface makes it a
/// steady state. 20 full solver steps (30 Myr) leave it within 1e-6 m, without lakes and
/// with them (closed and open lakes recomputed every step from the same surface).
#[test]
fn inverted_uplift_holds_the_surface() {
    let w = world(3, 48);
    let n = w.grid.n;
    let p = params(20, 30.0);
    let k = vec![p.k; n];
    let runoff = vec![0.3; n];
    let evap = vec![1.2; n];
    let lakes = rg_erosion::LakeInputs { runoff_m_yr: &runoff, net_evap_m_yr: &evap };
    for lake_in in [None, Some(&lakes)] {
        let (u, hold) = rg_erosion::invert_uplift(&w.grid, R, &w.h0, &w.ocean, &k, None, lake_in, &p);
        let mut h = w.h0.clone();
        rg_erosion::solve_full(&w.grid, R, &mut h, &w.ocean, &u, &k, None, None, lake_in, &p);
        let worst = (0..n).filter(|&i| !w.ocean[i]).map(|i| (h[i] - hold[i]).abs()).fold(0.0, f64::max);
        assert!(worst < 1e-6, "lakes {}: moved {worst} m", lake_in.is_some());
        assert!(u.iter().any(|&v| v > 1e-5), "real uplift somewhere");
    }
}
