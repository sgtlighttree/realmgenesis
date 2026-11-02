//! Lakes and closed (endorheic) basins from a steady-state water balance, on the depression
//! hierarchy of the unfilled surface (fill–spill–merge by *fluxes*, not volumes).
//!
//! A depression's lake grows until its net evaporation balances its inflow:
//! `Q = (E − P)·A_lake`, i.e. `A* = Q / (E − P)` (Matsubara & Howard 2009). If the lake's area
//! at the depression's lowest pass is smaller than `A*`, it fills to the pass and spills the
//! surplus `Q − (E − P)·A(pass)` into the basin beyond; if both sides are full they merge and
//! rise together; a chain that reaches the ocean is exorheic (open lakes). Otherwise the lake
//! stops below its pass: a **closed, endorheic** basin (Earth: ~20% of land drains this way,
//! mostly in dry climates). Humid basins (E ≤ P over the lake) always fill and spill.
//!
//! The depression hierarchy is the basin graph merged in ascending pass height (as in Barnes
//! et al. 2020, "Computing water flow through complex landscapes"); a basin that later gains
//! inflow is re-evaluated, so cascades resolve. Deterministic: ties go by index.

use rg_tectonics::grid::Grid;
use std::cmp::Reverse;
use std::collections::BinaryHeap;

/// Per-node lake result.
pub struct Lakes {
    /// Water-surface level where the node is under a lake (bed below the level), else NaN.
    pub level: Vec<f64>,
    /// The node's water ends in a closed lake (endorheic land), not the ocean.
    pub endorheic: Vec<bool>,
    /// The node is under a closed (terminal) lake.
    pub closed_lake: Vec<bool>,
}

#[derive(Clone)]
struct Root {
    /// Lake bed heights and areas, sorted by height, with prefix area sums.
    beds: Vec<(f64, f64)>,
    prefix: Vec<f64>,
    q: f64,         // inflow, m³/yr (own runoff + received overflow)
    ea: f64,        // Σ (E − P)·area over nodes (for the mean net evaporation)
    area: f64,      // Σ area
    passes: BinaryHeap<Reverse<(OrdF, u32)>>, // (pass height, other original basin)
    state: State,
}

#[derive(Clone, Copy, PartialEq)]
enum State {
    Open,                  // not decided yet
    Closed(f64),           // terminal lake at this level
    Spilled { level: f64, into: u32, kept: f64 }, // full to `level`, overflow into basin `into`
    Exorheic(f64),         // full to `level`, drains to the ocean
}

/// An f64 ordered by `total_cmp`, for the pass heap.
#[derive(Clone, Copy, PartialEq)]
struct OrdF(f64);
impl Eq for OrdF {}
impl Ord for OrdF {
    fn cmp(&self, o: &Self) -> std::cmp::Ordering {
        self.0.total_cmp(&o.0)
    }
}
impl PartialOrd for OrdF {
    fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(o))
    }
}

impl Root {
    fn net_evap(&self) -> f64 {
        self.ea / self.area.max(1e-9)
    }
    /// Lake area if the water stood at `z`.
    fn area_below(&self, z: f64) -> f64 {
        let k = self.beds.partition_point(|b| b.0 < z);
        if k == 0 { 0.0 } else { self.prefix[k - 1] }
    }
    /// Would the lake reach `z` with the current inflow?
    fn reaches(&self, z: f64) -> bool {
        let e = self.net_evap();
        e <= 0.0 || self.q >= e * self.area_below(z)
    }
    /// Level of a lake that stops short: where its area equals Q/(E − P).
    fn closed_level(&self) -> f64 {
        let target = self.q / self.net_evap().max(1e-12);
        let k = self.prefix.partition_point(|&a| a < target);
        let k = k.min(self.beds.len() - 1);
        // Between bed k and k+1 the area is prefix[k]; the level sits just above bed k.
        let next = self.beds.get(k + 1).map_or(self.beds[k].0 + 1.0, |b| b.0);
        0.5 * (self.beds[k].0 + next)
    }
    fn merge(&mut self, o: Root) {
        let mut beds = std::mem::take(&mut self.beds);
        beds.extend(o.beds);
        beds.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut acc = 0.0;
        self.prefix = beds.iter().map(|b| {
            acc += b.1;
            acc
        }).collect();
        self.beds = beds;
        self.q += o.q;
        self.ea += o.ea;
        self.area += o.area;
        for p in o.passes {
            self.passes.push(p);
        }
        self.state = State::Open;
    }
}

fn find(parent: &mut [u32], mut x: u32) -> u32 {
    while parent[x as usize] != x {
        parent[x as usize] = parent[parent[x as usize] as usize];
        x = parent[x as usize];
    }
    x
}

/// Lakes for surface `h` (m), with `ocean` nodes as the ocean, per-node `runoff` (m/yr of
/// water reaching the drainage), `net_evap` (open-water E − P, m/yr) and cell `area` (m²).
pub fn compute(grid: &Grid, h: &[f64], ocean: &[bool], runoff: &[f64], net_evap: &[f64], area: &[f64]) -> Lakes {
    let n = grid.n;
    // 1. Steepest-descent receivers on the unfilled surface; pits drain to themselves.
    let rec: Vec<u32> = (0..n)
        .map(|i| {
            if ocean[i] {
                return i as u32;
            }
            let mut best = (h[i], i as u32);
            for &j in grid.neighbors(i) {
                let hj = if ocean[j as usize] { h[j as usize].min(0.0) } else { h[j as usize] };
                if hj < best.0 || (hj == best.0 && j < best.1 && hj < h[i]) {
                    best = (hj, j);
                }
            }
            best.1
        })
        .collect();
    // 2. Basin label: the terminal each node drains to (ocean nodes: the ocean basin 0).
    let mut order: Vec<u32> = (0..n as u32).collect();
    order.sort_by(|&a, &b| h[a as usize].total_cmp(&h[b as usize]).then(a.cmp(&b)));
    const OCEAN: u32 = 0;
    let mut label = vec![u32::MAX; n];
    let mut n_basins = 1u32; // 0 = ocean
    for &i in &order {
        let i = i as usize;
        let r = rec[i] as usize;
        label[i] = if ocean[i] {
            OCEAN
        } else if r == i {
            n_basins += 1;
            n_basins - 1
        } else {
            label[r]
        };
    }
    let nb = n_basins as usize;
    // 3. Basin data and passes (lowest saddle to each neighbouring basin).
    let mut roots: Vec<Root> = (0..nb)
        .map(|_| Root { beds: Vec::new(), prefix: Vec::new(), q: 0.0, ea: 0.0, area: 0.0, passes: BinaryHeap::new(), state: State::Open })
        .collect();
    for i in 0..n {
        let b = label[i] as usize;
        if b == OCEAN as usize {
            continue;
        }
        roots[b].beds.push((h[i], area[i]));
        roots[b].q += runoff[i] * area[i];
        roots[b].ea += net_evap[i] * area[i];
        roots[b].area += area[i];
    }
    let mut best_pass: std::collections::HashMap<(u32, u32), f64> = std::collections::HashMap::new();
    for i in 0..n {
        for &j in grid.neighbors(i) {
            let (a, b) = (label[i], label[j as usize]);
            if a == b || (a == OCEAN && b == OCEAN) {
                continue;
            }
            let key = (a.min(b), a.max(b));
            let z = h[i].max(h[j as usize]).max(0.0);
            let e = best_pass.entry(key).or_insert(f64::INFINITY);
            if z < *e {
                *e = z;
            }
        }
    }
    let mut keys: Vec<_> = best_pass.into_iter().collect();
    keys.sort_by_key(|x| x.0);
    for ((a, b), z) in keys {
        if a != OCEAN {
            roots[a as usize].passes.push(Reverse((OrdF(z), b)));
        }
        if b != OCEAN {
            roots[b as usize].passes.push(Reverse((OrdF(z), a)));
        }
    }
    for r in roots.iter_mut().skip(1) {
        r.beds.sort_by(|a, b| a.0.total_cmp(&b.0));
        let mut acc = 0.0;
        r.prefix = r.beds.iter().map(|b| {
            acc += b.1;
            acc
        }).collect();
    }
    // 4. Fill–spill–merge.
    let mut parent: Vec<u32> = (0..nb as u32).collect();
    let mut work: Vec<u32> = (1..nb as u32).rev().collect();
    let mut guard = 0usize;
    while let Some(b) = work.pop() {
        guard += 1;
        if guard > 50 * nb + 1000 {
            break; // safety: never loop forever on pathological ties
        }
        let b = find(&mut parent, b);
        if b == OCEAN || !matches!(roots[b as usize].state, State::Open) {
            continue;
        }
        // The lowest pass to another root.
        let pass = loop {
            let Some(Reverse((OrdF(z), other))) = roots[b as usize].passes.pop() else { break None };
            let o = find(&mut parent, other);
            if o != b {
                break Some((z, o));
            }
        };
        let Some((z, o)) = pass else {
            let lvl = roots[b as usize].closed_level();
            roots[b as usize].state = State::Closed(lvl);
            continue;
        };
        if !roots[b as usize].reaches(z) {
            // Stops short: closed. (Put the pass back: more inflow may arrive later.)
            roots[b as usize].passes.push(Reverse((OrdF(z), o)));
            let lvl = roots[b as usize].closed_level();
            roots[b as usize].state = State::Closed(lvl);
            continue;
        }
        if o == OCEAN || matches!(roots[o as usize].state, State::Exorheic(_)) {
            roots[b as usize].state = State::Exorheic(z);
            continue;
        }
        let other_full = matches!(roots[o as usize].state, State::Spilled { .. }) || roots[o as usize].reaches(z);
        if other_full {
            // Merge and rise together. A spilled side gives back what it kept (its own
            // evaporation), since its overflow already sits in the merged inflow.
            let mut taken = roots[o as usize].clone();
            if let State::Spilled { kept, .. } = taken.state {
                taken.q += kept;
            }
            taken.passes.push(Reverse((OrdF(z), b)));
            roots[b as usize].merge(taken);
            parent[o as usize] = b;
            work.push(b);
        } else {
            // Spill the surplus into the other basin, which re-evaluates.
            let kept = roots[b as usize].net_evap() * roots[b as usize].area_below(z);
            let surplus = (roots[b as usize].q - kept).max(0.0);
            roots[b as usize].state = State::Spilled { level: z, into: o, kept };
            roots[o as usize].q += surplus;
            if let State::Closed(_) = roots[o as usize].state {
                roots[o as usize].state = State::Open;
            }
            work.push(o);
        }
    }
    // 5. Per node: level and fate.
    let mut level = vec![f64::NAN; n];
    let mut endorheic = vec![false; n];
    let mut closed_lake = vec![false; n];
    // Terminal fate of each root: follow spills.
    let fate = |mut b: u32, parent: &mut Vec<u32>, roots: &Vec<Root>| -> bool {
        for _ in 0..nb {
            b = find(parent, b);
            match roots[b as usize].state {
                State::Spilled { into, .. } => b = into,
                State::Closed(_) => return true,
                _ => return false,
            }
        }
        false
    };
    for i in 0..n {
        if ocean[i] {
            continue;
        }
        let b = find(&mut parent, label[i]);
        let lvl = match roots[b as usize].state {
            State::Closed(l) => {
                if h[i] < l {
                    closed_lake[i] = true;
                }
                l
            }
            State::Spilled { level, .. } | State::Exorheic(level) => level,
            State::Open => f64::NAN,
        };
        if h[i] < lvl {
            level[i] = lvl;
        }
        endorheic[i] = fate(b, &mut parent, &roots);
    }
    Lakes { level, endorheic, closed_lake }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A bowl on land, far from the sea: dry → a small closed lake; wet → it fills and spills
    /// to the sea (no endorheic land).
    #[test]
    fn dry_bowl_is_closed_wet_bowl_spills() {
        let grid = Grid::new(24);
        let n = grid.n;
        let centre = glam::DVec3::new(0.0, 0.0, 1.0);
        // Land everywhere except a band of ocean around the equator; a bowl around +Z.
        let h: Vec<f64> = (0..n)
            .map(|i| {
                let d = grid.dir[i];
                if d.z.abs() < 0.15 {
                    -1000.0
                } else {
                    let a = d.angle_between(centre);
                    if a < 0.5 { 500.0 + 1500.0 * (a / 0.5).powi(2) } else { 2000.0 - 1000.0 * ((a - 0.5) / 0.9).min(1.0) }
                }
            })
            .collect();
        let ocean: Vec<bool> = h.iter().map(|&v| v <= 0.0).collect();
        let area = vec![1.0e8; n];
        let run = |runoff: f64, evap: f64| {
            compute(&grid, &h, &ocean, &vec![runoff; n], &vec![evap; n], &area)
        };
        let dry = run(0.01, 1.5);
        let in_bowl: Vec<usize> = (0..n).filter(|&i| grid.dir[i].angle_between(centre) < 0.3).collect();
        assert!(in_bowl.iter().all(|&i| dry.endorheic[i]), "dry bowl should drain internally");
        assert!(in_bowl.iter().any(|&i| dry.closed_lake[i]), "dry bowl should hold a lake");
        let wet = run(1.0, 0.5);
        assert!(in_bowl.iter().all(|&i| !wet.endorheic[i]), "wet bowl should spill to the sea");
        assert!(in_bowl.iter().any(|&i| !wet.level[i].is_nan()), "wet bowl is full (open lake)");
    }
}
