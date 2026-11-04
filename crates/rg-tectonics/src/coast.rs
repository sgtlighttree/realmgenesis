//! R3 step 1: coast classification. Seaways are classified first, and coasts inherit the class.
//!
//! Classifying each stretch of coast with its own noise gives short active/passive/active
//! runs, which Earth doesn't have. Instead (design: `docs/design/coastline-conditioning.md` §1):
//! 1. **Landmasses**: connected continental-crust components above a size floor.
//! 2. **Regions**: every ocean node belongs to the landmass that "reaches" it first (a noisy
//!    multi-source growth, so the borders aren't geometric).
//! 3. **Seaways**: the region borders form the medial axis between landmasses. One seaway is
//!    one connected piece of that axis for one landmass pair. Africa and South America face
//!    each other twice (Atlantic and Pacific); those are two seaways.
//! 4. **Class per seaway**: opening (both coasts passive, a ridge on the axis), closing
//!    (both rims active, oceanic plates and a ridge in between) or one-sided (one rim active).
//!    The `active_share` knob decides how much coast ends up active.
//!
//! **Everything here runs on a fixed res-128 grid**, whatever the generation resolution
//! (advisor, s05). Normal and Ultra generate tectonics separately, and Ultra's detail is
//! `T_fine − up(T_canon)`. A seaway class that flipped between the two runs would put a whole
//! trench into that detail. On one fixed grid every discrete decision is identical at every
//! resolution; the generation grid only samples the results.

use crate::Params;
use crate::crust::{self, Crust};
use crate::grid::Grid;
use crate::mask::StrokeOp;
use crate::plates::Kind;
use crate::rng::Rng;
use crate::util::*;
use glam::DVec3;
use rg_noise::salt;
use std::f64::consts::PI;

/// Resolution of the classification grid (spacing ~78 km on an Earth-size planet).
pub const CLASS_RES: u32 = 128;
pub const NONE: u16 = u16::MAX;
/// Landmasses smaller than this are islands (Britain 209,000 km², Iceland 103,000).
pub const LANDMASS_MIN_KM2: f64 = 200_000.0;
/// Landmasses larger than this are split into 2–3 plates along sutures (Australia 7.7 M km²).
pub const SUTURE_AREA_KM2: f64 = 8_000_000.0;
/// Seaways narrower than this (half-width) are one-sided when active (Mediterranean).
pub const ONE_SIDED_HALF_WIDTH_M: f64 = 300_000.0;
/// Coast with no other landmass within twice this faces open ocean.
pub const OPEN_OCEAN_M: f64 = 3_000_000.0;
/// Typical distance between segment cuts along a coast.
pub const SEGMENT_SPACING_M: f64 = 1_500_000.0;
/// Chance that a segment of a closing seaway's rim is passive (a rifted microcontinent).
pub const PASSIVE_SEGMENT_P: f64 = 0.15;

/// Margin classes, per node.
pub mod margin {
    pub const NONE: u8 = 0;
    pub const PASSIVE: u8 = 1;
    pub const ACTIVE: u8 = 2;
}

const STREAM_COAST: u64 = 60;
const STREAM_REGION_COST: u32 = 61;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SeawayClass {
    Opening,
    Closing,
    OneSided,
}

#[derive(Clone, Debug)]
pub struct Seaway {
    /// The two landmasses; `b == NONE` for an open-ocean pseudo-seaway (coast that faces no
    /// other landmass).
    pub a: u16,
    pub b: u16,
    pub class: SeawayClass,
    /// One-sided seaways: the landmass whose coast is active.
    pub active_side: u16,
    /// Median distance from the facing coasts to the axis, metres.
    pub half_width_m: f64,
    /// Length of coast facing this seaway, metres, per side (`a`, `b`).
    pub coast_len_m: [f64; 2],
}

#[derive(Clone, Debug)]
pub struct Landmass {
    pub area_km2: f64,
    /// Seed directions of the sub-plates (empty: the landmass is one plate).
    pub split_seeds: Vec<DVec3>,
    /// Kind of the internal (suture) boundaries.
    pub suture: Kind,
    /// Width of ocean a passive open-ocean coast keeps on its own plate, metres.
    pub apron_m: f64,
}

/// The classification, on the fixed grid `grid`.
pub struct Classification {
    pub grid: Grid,
    pub crust: Crust,
    pub landmass_of: Vec<u16>,
    pub landmasses: Vec<Landmass>,
    pub seaways: Vec<Seaway>,
    /// Per node: the seaway faced (index into `seaways`), or `u32::MAX`.
    pub facing: Vec<u32>,
    /// Per node: the margin class of the nearest coast node.
    pub coast_margin: Vec<u8>,
    /// Share of coast length that ended up active.
    pub active_share: f64,
}

/// Per-node noisy cost for the landmass regions (resolution-independent: a coarse field).
pub fn region_cost(grid: &Grid, seed: u64) -> Vec<f32> {
    let s = salt(seed, STREAM_REGION_COST);
    grid.coarse_field(64, |d| (0.8 + 0.4 * fbm01(d, s, 2.5, 3)) as f32)
}

/// Grow landmass regions over the whole grid from the labelled landmass nodes.
pub fn regions(grid: &Grid, landmass_of: &[u16], n_land: usize, cost: &[f32]) -> Vec<u16> {
    if n_land == 0 {
        return vec![NONE; grid.n];
    }
    let seeds: Vec<u32> = (0..grid.n as u32).filter(|&a| landmass_of[a as usize] != NONE).collect();
    let seed_region: Vec<u16> = seeds.iter().map(|&a| landmass_of[a as usize]).collect();
    let speed = vec![1.0f32; n_land];
    let axis = vec![DVec3::Z; n_land];
    let stretch = vec![0.0f32; n_land];
    grid.grow_regions(&seeds, &seed_region, &speed, &axis, &stretch, cost)
}

/// Landmass coast nodes: labelled nodes next to non-continental crust.
pub fn coast_nodes(grid: &Grid, cont: &[bool], landmass_of: &[u16]) -> Vec<u32> {
    (0..grid.n)
        .filter(|&a| landmass_of[a] != NONE && grid.neighbors(a).iter().any(|&b| !cont[b as usize]))
        .map(|a| a as u32)
        .collect()
}

/// Length of coast one node stands for, metres (the side of its cell).
fn cell_len_m(grid: &Grid, a: usize, radius_m: f64) -> f64 {
    (grid.weight[a] as f64 * 4.0 * PI).sqrt() * radius_m
}

impl Classification {
    pub fn build(p: &Params) -> Classification {
        let grid = Grid::new(CLASS_RES);
        let crust = crust::build(&grid, p);
        // What the landmasses are: the continental crust (procedural; the coast emerges
        // later) or the land mask (strict, R3.1b).
        let src: Vec<bool> = if p.conditioning.strict {
            let value = crust::continentalness(&grid, p);
            crate::mask::land_cover(&grid, p, &value).iter().map(|&c| c > 0.5).collect()
        } else {
            crust.cont.clone()
        };
        let r = p.radius_m;
        let n = grid.n;
        let sphere_km2 = 4.0 * PI * (r / 1000.0).powi(2);

        // --- 1. landmasses: crust components above the size floor (scan order is stable:
        // this grid never changes).
        let mut comp = vec![u32::MAX; n];
        let mut comps: Vec<(f64, Vec<u32>)> = Vec::new();
        for start in 0..n {
            if !src[start] || comp[start] != u32::MAX {
                continue;
            }
            let id = comps.len() as u32;
            comp[start] = id;
            let mut stack = vec![start as u32];
            let mut members = Vec::new();
            let mut area = 0.0;
            while let Some(a) = stack.pop() {
                members.push(a);
                area += grid.weight[a as usize] as f64 * sphere_km2;
                for &b in grid.neighbors(a as usize) {
                    if src[b as usize] && comp[b as usize] == u32::MAX {
                        comp[b as usize] = id;
                        stack.push(b);
                    }
                }
            }
            comps.push((area, members));
        }
        let mut landmass_of = vec![NONE; n];
        let mut areas = Vec::new();
        for (area, members) in &comps {
            if *area >= LANDMASS_MIN_KM2 {
                let l = areas.len() as u16;
                for &a in members {
                    landmass_of[a as usize] = l;
                }
                areas.push(*area);
            }
        }
        let n_land = areas.len();

        // --- 2. regions, coasts and the medial axis
        let cost = region_cost(&grid, p.seed);
        let region = regions(&grid, &landmass_of, n_land, &cost);
        let coast = coast_nodes(&grid, &src, &landmass_of);
        let is_border = |a: usize| {
            landmass_of[a] == NONE && region[a] != NONE && grid.neighbors(a).iter().any(|&b| region[b as usize] != region[a])
        };
        let pair_of = |a: usize| -> (u16, u16) {
            let other = grid
                .neighbors(a)
                .iter()
                .map(|&b| region[b as usize])
                .filter(|&q| q != region[a])
                .min()
                .unwrap_or(NONE);
            (region[a].min(other), region[a].max(other))
        };
        let border: Vec<bool> = (0..n).map(is_border).collect();
        let pair: Vec<(u16, u16)> = (0..n).map(|a| if border[a] { pair_of(a) } else { (NONE, NONE) }).collect();

        // Seaways: connected pieces of the axis, per landmass pair.
        let mut seaway_of = vec![u32::MAX; n];
        let mut seaways: Vec<Seaway> = Vec::new();
        for start in 0..n {
            if !border[start] || seaway_of[start] != u32::MAX {
                continue;
            }
            let id = seaways.len() as u32;
            let key = pair[start];
            seaway_of[start] = id;
            let mut stack = vec![start];
            while let Some(a) = stack.pop() {
                for &b in grid.neighbors(a) {
                    let b = b as usize;
                    if border[b] && seaway_of[b] == u32::MAX && pair[b] == key {
                        seaway_of[b] = id;
                        stack.push(b);
                    }
                }
            }
            seaways.push(Seaway {
                a: key.0,
                b: key.1,
                class: SeawayClass::Opening,
                active_side: NONE,
                half_width_m: 0.0,
                coast_len_m: [0.0; 2],
            });
        }
        let border_nodes: Vec<u32> = (0..n as u32).filter(|&a| border[a as usize]).collect();
        let near_border = grid.nearest_source(&border_nodes, PI as f32);

        // --- 3. which seaway each node faces; open-ocean pseudo-seaways per landmass.
        let mut open_of = vec![u32::MAX; n_land];
        let mut facing = vec![u32::MAX; n];
        for a in 0..n {
            let src = near_border.src[a];
            let far = src == u32::MAX || near_border.dist_rad[a] as f64 * r > OPEN_OCEAN_M;
            let own = if landmass_of[a] != NONE { landmass_of[a] } else { region[a] };
            if far {
                if own == NONE {
                    continue;
                }
                if open_of[own as usize] == u32::MAX {
                    open_of[own as usize] = seaways.len() as u32;
                    seaways.push(Seaway {
                        a: own,
                        b: NONE,
                        class: SeawayClass::Opening,
                        active_side: NONE,
                        half_width_m: f64::INFINITY,
                        coast_len_m: [0.0; 2],
                    });
                }
                facing[a] = open_of[own as usize];
            } else {
                facing[a] = seaway_of[src as usize];
            }
        }

        // --- 4. half-widths and coast lengths
        let mut dists: Vec<Vec<f64>> = vec![Vec::new(); seaways.len()];
        let mut total_coast = 0.0;
        for &c in &coast {
            let c = c as usize;
            let len = cell_len_m(&grid, c, r);
            total_coast += len;
            let s = facing[c];
            if s == u32::MAX {
                continue;
            }
            let sw = &mut seaways[s as usize];
            let side = if landmass_of[c] == sw.a { 0 } else { 1 };
            sw.coast_len_m[side] += len;
            dists[s as usize].push(near_border.dist_rad[c] as f64 * r);
        }
        for (sw, d) in seaways.iter_mut().zip(dists.iter_mut()) {
            if sw.b != NONE && !d.is_empty() {
                d.sort_by(f64::total_cmp);
                sw.half_width_m = d[d.len() / 2];
            }
        }

        // --- 5. class rolls. Each seaway has an "active option" (one-sided when narrow,
        // closing otherwise). Options are taken in order of `u / w(D)` (wide seaways first,
        // shuffled by `u`) while the active coast length stays within the target share.
        let mut rng = Rng::new(p.seed, STREAM_COAST);
        let rolls: Vec<(f64, f64)> = seaways.iter().map(|_| (rng.f64(), rng.f64())).collect();
        let mut order: Vec<(f64, usize, f64)> = Vec::new(); // (score, seaway, active length)
        for (s, sw) in seaways.iter_mut().enumerate() {
            let (u, side) = rolls[s];
            let total = sw.coast_len_m[0] + sw.coast_len_m[1];
            let active_len = if sw.b != NONE && sw.half_width_m < ONE_SIDED_HALF_WIDTH_M {
                sw.active_side = if side < 0.5 { sw.a } else { sw.b };
                sw.coast_len_m[if sw.active_side == sw.a { 0 } else { 1 }]
            } else {
                total * (1.0 - PASSIVE_SEGMENT_P)
            };
            let w = 0.5 + sw.half_width_m.min(4_000_000.0) / 2_000_000.0;
            if active_len > 0.0 {
                order.push((u / w, s, active_len));
            }
        }
        order.sort_by(|x, y| x.0.total_cmp(&y.0).then(x.1.cmp(&y.1)));
        let target = p.active_share.clamp(0.0, 1.0) as f64 * total_coast;
        let mut cum = 0.0;
        for &(_, s, len) in &order {
            if cum + 0.5 * len <= target {
                let sw = &mut seaways[s];
                sw.class = if sw.active_side != NONE { SeawayClass::OneSided } else { SeawayClass::Closing };
                cum += len;
            }
        }
        for sw in seaways.iter_mut() {
            if sw.class != SeawayClass::OneSided {
                sw.active_side = NONE;
            }
        }

        // --- 6. segments: coast nodes are grouped by the nearest of a seed-keyed point set
        // (~1,500 km apart, fixed on the sphere); a closing rim's segment may roll passive.
        let points = segment_points(p.seed, SEGMENT_SPACING_M / r);
        let mut node_margin = vec![margin::NONE; n];
        let mut active_len = 0.0;
        for &c in &coast {
            let c = c as usize;
            let s = facing[c];
            let m = if s == u32::MAX {
                margin::PASSIVE
            } else {
                let sw = &seaways[s as usize];
                match sw.class {
                    SeawayClass::Opening => margin::PASSIVE,
                    SeawayClass::OneSided => {
                        if landmass_of[c] == sw.active_side { margin::ACTIVE } else { margin::PASSIVE }
                    }
                    SeawayClass::Closing => {
                        let j = nearest_point(&points, grid.dir[c]);
                        let roll = Rng::new(p.seed ^ ((s as u64) << 32 | j as u64), STREAM_COAST + 1).f64();
                        if roll < PASSIVE_SEGMENT_P { margin::PASSIVE } else { margin::ACTIVE }
                    }
                }
            };
            // Painted overrides win (R3.3); the last stroke covering the node decides.
            let m = p
                .conditioning
                .strokes
                .iter()
                .rev()
                .find(|st| matches!(st.op, StrokeOp::Active | StrokeOp::Passive) && st.weight(grid.dir[c], r) > 0.5)
                .map_or(m, |st| if st.op == StrokeOp::Active { margin::ACTIVE } else { margin::PASSIVE });
            node_margin[c] = m;
            if m == margin::ACTIVE {
                active_len += cell_len_m(&grid, c, r);
            }
        }
        let near_coast = grid.nearest_source(&coast, PI as f32);
        let coast_margin: Vec<u8> =
            (0..n).map(|a| if near_coast.src[a] == u32::MAX { margin::NONE } else { node_margin[near_coast.src[a] as usize] }).collect();

        // --- 7. sutures and aprons, per landmass.
        let mut rng = Rng::new(p.seed, STREAM_COAST + 2);
        let mut landmasses = Vec::with_capacity(n_land);
        for (l, &area) in areas.iter().enumerate() {
            let k = if rng.f64() < 0.5 { 2 } else { 3 };
            let suture = if rng.f64() < 0.7 { Kind::Collision } else { Kind::Rift };
            let apron_m = rng.range(1_500_000.0, 2_500_000.0);
            let mut split_seeds: Vec<DVec3> = Vec::new();
            if area > SUTURE_AREA_KM2 {
                for _ in 0..k {
                    let mut best: Option<(f64, DVec3)> = None;
                    for _ in 0..40 {
                        let v = rng.unit_vec();
                        let c = grid.nearest_node(v);
                        if landmass_of[c] != l as u16 {
                            continue;
                        }
                        let d = grid.dir[c];
                        let score = split_seeds.iter().map(|s| s.dot(d).clamp(-1.0, 1.0).acos()).fold(PI, f64::min);
                        if best.is_none_or(|b| score > b.0) {
                            best = Some((score, d));
                        }
                    }
                    if let Some((_, d)) = best {
                        split_seeds.push(d);
                    }
                }
                if split_seeds.len() < 2 {
                    split_seeds.clear();
                }
            }
            landmasses.push(Landmass { area_km2: area, split_seeds, suture, apron_m });
        }

        let mut cls = Classification {
            grid,
            crust,
            landmass_of,
            landmasses,
            seaways,
            facing,
            coast_margin,
            active_share: if total_coast > 0.0 { active_len / total_coast } else { 0.0 },
        };
        if p.conditioning.strict {
            // Strict: crust = land + class-dependent margins (needs the classes, so last).
            cls.crust = crust::strict(&cls.grid, p, &src, &cls);
        }
        cls
    }

    /// Landmass labels for another grid's continental crust: fine nodes take the label of the
    /// nearest classification node, then labels spread through the rest of each crust
    /// component. Crust the classification never saw (small islands) stays `NONE`.
    pub fn landmass_labels(&self, grid: &Grid, cont: &[bool]) -> Vec<u16> {
        let mut label = vec![NONE; grid.n];
        let mut stack: Vec<u32> = Vec::new();
        for a in 0..grid.n {
            if cont[a] {
                let l = self.landmass_of[self.grid.nearest_node(grid.dir[a])];
                if l != NONE {
                    label[a] = l;
                    stack.push(a as u32);
                }
            }
        }
        while let Some(a) = stack.pop() {
            let l = label[a as usize];
            for &b in grid.neighbors(a as usize) {
                if cont[b as usize] && label[b as usize] == NONE {
                    label[b as usize] = l;
                    stack.push(b);
                }
            }
        }
        label
    }

    /// The seaway faced at direction `d` (`u32::MAX`: none).
    pub fn facing_at(&self, d: DVec3) -> u32 {
        self.facing[self.grid.nearest_node(d)]
    }

    /// The margin class of the coast nearest to direction `d`.
    pub fn margin_at(&self, d: DVec3) -> u8 {
        self.coast_margin[self.grid.nearest_node(d)]
    }
}

/// A seed-keyed set of points roughly `spacing` radians apart (dart throwing).
fn segment_points(seed: u64, spacing: f64) -> Vec<DVec3> {
    let mut rng = Rng::new(seed, STREAM_COAST + 3);
    let target = (4.0 * PI / (spacing * spacing)).ceil() as usize;
    let min_dot = (0.8 * spacing).cos();
    let mut pts: Vec<DVec3> = Vec::with_capacity(target);
    for _ in 0..target * 30 {
        if pts.len() >= target {
            break;
        }
        let v = rng.unit_vec();
        if pts.iter().all(|q| q.dot(v) < min_dot) {
            pts.push(v);
        }
    }
    pts
}

fn nearest_point(pts: &[DVec3], d: DVec3) -> usize {
    (0..pts.len()).max_by(|&i, &j| pts[i].dot(d).total_cmp(&pts[j].dot(d)).then(j.cmp(&i))).unwrap_or(0)
}
