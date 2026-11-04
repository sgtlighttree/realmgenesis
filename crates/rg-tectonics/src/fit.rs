//! R3 steps 3–4: plates **fitted to the coast classes**, and their motion solved to match.
//!
//! Design: `docs/design/coastline-conditioning.md` §3–4. In short:
//! - **Continental plates**: each landmass is a plate (big ones are split along sutures). It
//!   keeps its ocean out to the axis of every opening seaway (South America / Africa) and out
//!   to the trench offset off every active coast.
//! - **Oceanic plates** fill what's left (closing seaways, open ocean). They grow from seeds
//!   and never cross a closing seaway's axis, so that axis is always a plate boundary (its ridge).
//! - **Kinds are prescribed** by the classification (active coast → subduction with the
//!   ocean diving, opening axis → ridge, suture → collision or rift). Only boundaries between
//!   oceanic plates inside one region are left to the motion.
//! - **Motion**: one rotation vector per plate, solved by least squares so the convergence
//!   rate across each prescribed boundary is near its target (+55 mm/yr subduction, −40
//!   ridge, ...). The solve sets rates only; the kind always comes from the classification.
//!
//! The partition is built twice: on the fixed classification grid (which picks the oceanic
//! seeds and solves the motion, so plates are identical at every resolution) and on the
//! generation grid (which only re-draws the same plates with finer borders).

use crate::coast::{self, Classification, NONE, SeawayClass, margin};
use crate::crust::Crust;
use crate::grid::Grid;
use crate::plates::{BoundaryNode, Kind, PlateModel, velocity};
use crate::rng::Rng;
use crate::util::*;
use crate::{Params, Plate};
use glam::DVec3;
use rg_noise::salt;
use std::f64::consts::PI;

/// Distance from an active coast (the crust edge) to the plate boundary; the trench axis lies
/// 30 km beyond it (`relief.rs`).
pub const TRENCH_OFFSET_M: f64 = 120_000.0;
/// A passive coast's plate meets an oceanic plate here as a transform (lateral boundary).
pub const LATERAL_TRANSFORM_M: f64 = 400_000.0;
/// Unclaimed ocean pieces smaller than this join a neighbouring plate instead of getting one.
const SMALL_OCEAN_KM2: f64 = 500_000.0;
/// The margin field is drawn this far either side of the coast.
const MARGIN_BAND_M: f64 = 300_000.0;
const STREAM_FIT: u64 = 65;
const STREAM_PLATE_COST: u32 = 66;
/// Rate floors where the solved rate has the wrong sign or is too small, mm/yr.
const MIN_SUBDUCTION: f64 = 15.0;
const MIN_COLLISION: f64 = 10.0;
const MIN_RIDGE: f64 = 10.0;
const MIN_RIFT: f64 = 5.0;

/// How well the solved motion matches the prescribed kinds (on the classification grid).
#[derive(Clone, Debug, Default)]
pub struct FitStats {
    /// Length-weighted RMS of (rate − target) over prescribed boundaries, mm/yr.
    pub rms_residual_mm_yr: f64,
    /// Length share of prescribed convergent/divergent boundaries whose solved rate has the
    /// wrong sign (reported, not gated: wide closing seaways have some).
    pub wrong_sign_share: f64,
    /// Share of coast length that is active.
    pub active_share: f64,
    pub seaways: [u32; 3],
    pub plates: u32,
}

/// Per-grid geometry derived from the classification.
struct Geometry {
    landmass_of: Vec<u16>,
    region: Vec<u16>,
    coast_dist_m: Vec<f32>,
    coast_src: Vec<u32>,
    /// Margin class of the nearest coast node.
    margin: Vec<u8>,
    facing: Vec<u32>,
}

fn geometry(grid: &Grid, crust: &Crust, cls: &Classification, seed: u64, r: f64) -> Geometry {
    let landmass_of = cls.landmass_labels(grid, &crust.land);
    let cost = coast::region_cost(grid, seed);
    let region = coast::regions(grid, &landmass_of, cls.landmasses.len(), &cost);
    let coast = coast::coast_nodes(grid, &crust.land, &landmass_of);
    let near = grid.nearest_source(&coast, PI as f32);
    let margin = par_map(grid.n, |a| {
        let s = near.src[a];
        if s == u32::MAX { margin::NONE } else { cls.margin_at(grid.dir[s as usize]) }
    });
    let facing = par_map(grid.n, |a| cls.facing_at(grid.dir[a]));
    let coast_dist_m = near.dist_rad.iter().map(|&d| (d as f64 * r) as f32).collect();
    Geometry { landmass_of, region, coast_dist_m, coast_src: near.src, margin, facing }
}

/// Plate-id layout, fixed by the classification: landmass sub-plates first, then oceans.
struct Layout {
    base: Vec<u16>,
    n_cont: usize,
    /// Landmass of each continental plate.
    plate_land: Vec<u16>,
}

fn layout(cls: &Classification) -> Layout {
    let mut base = Vec::new();
    let mut plate_land = Vec::new();
    for (l, lm) in cls.landmasses.iter().enumerate() {
        base.push(plate_land.len() as u16);
        for _ in 0..lm.split_seeds.len().max(1) {
            plate_land.push(l as u16);
        }
    }
    Layout { base, n_cont: plate_land.len(), plate_land }
}

fn plate_cost(grid: &Grid, seed: u64) -> Vec<f32> {
    let s = salt(seed, STREAM_PLATE_COST);
    grid.coarse_field(96, |d| {
        let q = d + warp3(d, s.wrapping_add(5), 1.5, 2) * 0.25;
        let v = fbm01(q, s, 2.4, 4);
        (0.4 + 2.0 * v * v) as f32
    })
}

/// The node of `allowed` nearest to direction `d` (brute force; rare).
fn nearest_allowed(grid: &Grid, d: DVec3, allowed: impl Fn(usize) -> bool) -> Option<usize> {
    let a = grid.nearest_node(d);
    if allowed(a) {
        return Some(a);
    }
    (0..grid.n).filter(|&b| allowed(b)).max_by(|&x, &y| grid.dir[x].dot(d).total_cmp(&grid.dir[y].dot(d)).then(y.cmp(&x)))
}

/// Plate id of every landmass node (sub-plates grown from the suture seeds); `NONE` elsewhere.
fn land_plates(grid: &Grid, geo: &Geometry, cls: &Classification, lay: &Layout, cost: &[f32]) -> Vec<u16> {
    let mut seeds = Vec::new();
    for (l, lm) in cls.landmasses.iter().enumerate() {
        for (j, &d) in lm.split_seeds.iter().enumerate() {
            if let Some(a) = nearest_allowed(grid, d, |b| geo.landmass_of[b] == l as u16) {
                seeds.push((a as u32, lay.base[l] + j as u16));
            }
        }
    }
    let allowed: Vec<bool> = geo
        .landmass_of
        .iter()
        .map(|&l| l != NONE && !cls.landmasses[l as usize].split_seeds.is_empty())
        .collect();
    let grown = grid.grow_labelled(&seeds, &allowed, cost, Some(&geo.landmass_of));
    (0..grid.n)
        .map(|a| {
            let l = geo.landmass_of[a];
            if l == NONE {
                NONE
            } else if grown[a] != NONE {
                grown[a]
            } else {
                lay.base[l as usize]
            }
        })
        .collect()
}

/// Continental claims on the ocean (design §3.1). Returns the plate per node (`NONE`:
/// unclaimed, for the oceanic plates).
fn claims(grid: &Grid, geo: &Geometry, cls: &Classification, lay: &Layout, land: &[u16], cost: &[f32]) -> Vec<u16> {
    let n = grid.n;
    let mut claim_land = vec![NONE; n]; // which landmass claims the node
    let mut plate = land.to_vec(); // known plate (NONE: unknown or unclaimed)
    for a in 0..n {
        if geo.landmass_of[a] != NONE {
            claim_land[a] = geo.landmass_of[a];
            continue;
        }
        let r = geo.region[a];
        if r == NONE {
            continue;
        }
        let src = geo.coast_src[a];
        let near = geo.coast_dist_m[a] as f64;
        let active = geo.margin[a] == margin::ACTIVE;
        // Sub-plate of the nearest coast, if that coast is the region's own landmass.
        let sub = if src != u32::MAX && geo.landmass_of[src as usize] == r { land[src as usize] } else { NONE };
        // Keyed on the margin of the nearest coast first, so painted overrides (R3.3) need no
        // special case: an active coast keeps only the trench offset; beyond it the ocean
        // belongs to the plate across the seaway (opening/one-sided: it dives under this
        // coast) or to the oceanic plates (closing, open ocean). A passive coast keeps its
        // ocean out to the axis (or its apron, facing open ocean).
        let s = geo.facing[a];
        let sw = (s != u32::MAX).then(|| &cls.seaways[s as usize]).filter(|sw| sw.a == r || sw.b == r);
        let open = sw.is_none_or(|sw| sw.b == NONE);
        let owner = if active {
            if near < TRENCH_OFFSET_M {
                Some(r)
            } else {
                match sw {
                    Some(sw) if !open && sw.class != SeawayClass::Closing => Some(if sw.a == r { sw.b } else { sw.a }),
                    _ => None,
                }
            }
        } else if open {
            (near < cls.landmasses[r as usize].apron_m).then_some(r)
        } else {
            Some(r)
        };
        if let Some(o) = owner {
            claim_land[a] = o;
            if o == r {
                plate[a] = sub;
            }
        }
    }
    // Claimed nodes whose sub-plate is unknown take the nearest known one of their landmass.
    let seeds: Vec<(u32, u16)> =
        (0..n).filter(|&a| claim_land[a] != NONE && plate[a] != NONE).map(|a| (a as u32, plate[a])).collect();
    let allowed: Vec<bool> = claim_land.iter().map(|&l| l != NONE).collect();
    let grown = grid.grow_labelled(&seeds, &allowed, cost, Some(&claim_land));
    for a in 0..n {
        if claim_land[a] != NONE && plate[a] == NONE {
            plate[a] = if grown[a] != NONE { grown[a] } else { lay.base[claim_land[a] as usize] };
        }
    }
    plate
}

/// Oceanic seed directions, picked on the classification grid: one per sizeable unclaimed
/// piece, then extra ones spread out until the plate-count target is met.
fn ocean_seeds(grid: &Grid, p: &Params, geo: &Geometry, claimed: &[u16], n_cont: usize, rng: &mut Rng) -> Vec<DVec3> {
    let n = grid.n;
    let sphere_km2 = 4.0 * PI * (p.radius_m / 1000.0).powi(2);
    let mut comp = vec![u32::MAX; n];
    let mut big = vec![false; n];
    let mut seeds: Vec<DVec3> = Vec::new();
    let mut n_comp = 0;
    for start in 0..n {
        if claimed[start] != NONE || comp[start] != u32::MAX {
            continue;
        }
        comp[start] = n_comp;
        let mut stack = vec![start];
        let mut members = Vec::new();
        let (mut area, mut centroid) = (0.0, DVec3::ZERO);
        while let Some(a) = stack.pop() {
            members.push(a);
            area += grid.weight[a] as f64 * sphere_km2;
            centroid += grid.dir[a];
            for &b in grid.neighbors(a) {
                let b = b as usize;
                if claimed[b] == NONE && comp[b] == u32::MAX && geo.region[b] == geo.region[a] {
                    comp[b] = n_comp;
                    stack.push(b);
                }
            }
        }
        n_comp += 1;
        if area >= SMALL_OCEAN_KM2 {
            for &a in &members {
                big[a] = true;
            }
            let c = centroid.normalize_or_zero();
            let a = *members.iter().max_by(|&&x, &&y| grid.dir[x].dot(c).total_cmp(&grid.dir[y].dot(c)).then(y.cmp(&x))).unwrap();
            seeds.push(grid.dir[a]);
        }
    }
    let target = (p.plate_count.clamp(2, 60) as usize).saturating_sub(n_cont).max(seeds.len());
    let mut tries = 0;
    while seeds.len() < target && tries < 400 {
        let mut best: Option<(f64, DVec3)> = None;
        for _ in 0..20 {
            tries += 1;
            let a = grid.nearest_node(rng.unit_vec());
            if !big[a] {
                continue;
            }
            let d = grid.dir[a];
            let score = seeds.iter().map(|s| s.dot(d).clamp(-1.0, 1.0).acos()).fold(PI, f64::min) * rng.range(0.8, 1.2);
            if best.is_none_or(|b| score > b.0) {
                best = Some((score, d));
            }
        }
        if let Some((_, d)) = best {
            seeds.push(d);
        }
    }
    seeds
}

/// The full partition on one grid: claims, then oceanic plates grown from `seeds`, then
/// anything left over joins its nearest plate.
fn partition(grid: &Grid, geo: &Geometry, claimed: &[u16], seeds: &[DVec3], n_cont: usize, cost: &[f32]) -> Vec<u16> {
    let n = grid.n;
    let unclaimed: Vec<bool> = claimed.iter().map(|&c| c == NONE).collect();
    let mut seed_nodes = Vec::new();
    for (i, &d) in seeds.iter().enumerate() {
        // The nearest unclaimed node within a few hops (borders differ a little by resolution).
        let start = grid.nearest_node(d);
        let mut found = unclaimed[start].then_some(start);
        let mut frontier = vec![start];
        let mut seen = vec![start];
        for _ in 0..8 {
            if found.is_some() {
                break;
            }
            let mut next = Vec::new();
            for &x in &frontier {
                for &y in grid.neighbors(x) {
                    let y = y as usize;
                    if !seen.contains(&y) {
                        seen.push(y);
                        next.push(y);
                        if found.is_none() && unclaimed[y] {
                            found = Some(y);
                        }
                    }
                }
            }
            frontier = next;
        }
        if let Some(a) = found {
            seed_nodes.push((a as u32, (n_cont + i) as u16));
        }
    }
    let grown = grid.grow_labelled(&seed_nodes, &unclaimed, cost, Some(&geo.region));
    let mut plate: Vec<u16> = (0..n).map(|a| if claimed[a] != NONE { claimed[a] } else { grown[a] }).collect();
    if plate.contains(&NONE) {
        // Unseeded ocean pieces (too small for a plate of their own) join the nearest
        // *oceanic* plate, even across a closing axis. Joining the continental plate would
        // erase the trench off its active coast.
        let seeds: Vec<(u32, u16)> =
            (0..n).filter(|&a| plate[a] != NONE && claimed[a] == NONE).map(|a| (a as u32, plate[a])).collect();
        let fill = grid.grow_labelled(&seeds, &unclaimed, cost, None);
        for a in 0..n {
            if plate[a] == NONE {
                plate[a] = fill[a];
            }
        }
    }
    if plate.contains(&NONE) {
        let seeds: Vec<(u32, u16)> = (0..n).filter(|&a| plate[a] != NONE).map(|a| (a as u32, plate[a])).collect();
        let fill = grid.grow_labelled(&seeds, &vec![true; n], cost, None);
        for a in 0..n {
            if plate[a] == NONE {
                plate[a] = if fill[a] != NONE { fill[a] } else { 0 };
            }
        }
    }
    plate
}

/// One border node, measured from its 3-hop neighbourhood (as in s01's plates.rs).
struct Border {
    node: usize,
    own: u16,
    other: u16,
    /// Unit normal from `own` toward `other`, tangent to the sphere.
    normal: DVec3,
    a_cont: bool,
    b_cont: bool,
    /// The two plates meet across a landmass-region border (a seaway axis).
    across: bool,
}

fn borders(grid: &Grid, plate_of: &[u16], cont: &[bool], region: &[u16]) -> Vec<Border> {
    let n = grid.n;
    let mut out = Vec::new();
    let mut stamp = vec![0u32; n];
    let mut now = 0u32;
    for a in 0..n {
        if !grid.neighbors(a).iter().any(|&b| plate_of[b as usize] != plate_of[a]) {
            continue;
        }
        now += 1;
        stamp[a] = now;
        let mut ball = vec![a as u32];
        let mut frontier = vec![a as u32];
        for _ in 0..3 {
            let mut next = Vec::new();
            for &x in &frontier {
                for &y in grid.neighbors(x as usize) {
                    if stamp[y as usize] != now {
                        stamp[y as usize] = now;
                        ball.push(y);
                        next.push(y);
                    }
                }
            }
            frontier = next;
        }
        let own = plate_of[a];
        let mut counts: Vec<(u16, u32)> = Vec::new();
        for &x in &ball {
            let pl = plate_of[x as usize];
            if pl != own {
                match counts.iter_mut().find(|c| c.0 == pl) {
                    Some(c) => c.1 += 1,
                    None => counts.push((pl, 1)),
                }
            }
        }
        counts.sort_by_key(|&(pl, c)| (std::cmp::Reverse(c), pl));
        let Some(&(other, _)) = counts.first() else { continue };
        let (mut sa, mut ca, mut cont_a) = (DVec3::ZERO, 0u32, 0u32);
        let (mut sb, mut cb, mut cont_b) = (DVec3::ZERO, 0u32, 0u32);
        for &x in &ball {
            let x = x as usize;
            if plate_of[x] == own {
                sa += grid.dir[x];
                ca += 1;
                cont_a += cont[x] as u32;
            } else if plate_of[x] == other {
                sb += grid.dir[x];
                cb += 1;
                cont_b += cont[x] as u32;
            }
        }
        let pa = grid.dir[a];
        let toward = sb / cb as f64 - sa / ca as f64;
        let t = toward - pa * toward.dot(pa);
        if t.length() < 1e-9 {
            continue;
        }
        // Direct neighbours only: a one-sided trench can lie within a few cells of the axis.
        let across = grid.neighbors(a).iter().any(|&b| plate_of[b as usize] != own && region[b as usize] != region[a]);
        out.push(Border {
            node: a,
            own,
            other,
            normal: t.normalize(),
            a_cont: cont_a * 2 >= ca,
            b_cont: cont_b * 2 >= cb,
            across,
        });
    }
    out
}

/// The kind the classification prescribes at a border, and (for subduction) the diving plate.
/// `None`: free, decided by the motion.
fn prescribe(b: &Border, geo: &Geometry, cls: &Classification, lay: &Layout, cont: &[bool]) -> Option<(Kind, u16)> {
    let land = |pl: u16| if (pl as usize) < lay.n_cont { lay.plate_land[pl as usize] } else { NONE };
    let (lp, lq) = (land(b.own), land(b.other));
    let a = b.node;
    match (lp != NONE, lq != NONE) {
        (true, true) if lp == lq => {
            Some((if cont[a] { cls.landmasses[lp as usize].suture } else { Kind::Transform }, b.own))
        }
        (true, true) => {
            if b.across {
                Some((if b.a_cont && b.b_cont { Kind::Rift } else { Kind::Ridge }, b.own))
            } else {
                // One-sided seaway: both plates in the active side's region; the other
                // landmass's plate dives.
                let x = geo.region[a];
                if lp != x {
                    Some((Kind::Subduction, b.own))
                } else if lq != x {
                    Some((Kind::Subduction, b.other))
                } else {
                    None
                }
            }
        }
        (true, false) | (false, true) => {
            let ocean = if lp != NONE { b.other } else { b.own };
            if b.across {
                Some((Kind::Ridge, b.own))
            } else if geo.margin[a] == margin::ACTIVE {
                Some((Kind::Subduction, ocean))
            } else if (geo.coast_dist_m[a] as f64) < LATERAL_TRANSFORM_M {
                Some((Kind::Transform, b.own))
            } else {
                Some((Kind::Ridge, b.own))
            }
        }
        (false, false) => b.across.then_some((Kind::Ridge, b.own)),
    }
}

/// Target convergence rate (mm/yr, > 0 closing) and row weight per kind.
fn target(kind: Kind) -> (f64, f64) {
    match kind {
        Kind::Subduction => (55.0, 1.0),
        Kind::Collision => (30.0, 1.0),
        Kind::Ridge => (-40.0, 1.0),
        Kind::Rift => (-10.0, 1.0),
        Kind::Transform => (0.0, 0.3),
    }
}

/// Solve `m x = rhs` (dense, partial pivoting). `m` is `n × n`, row-major.
fn solve_dense(mut m: Vec<f64>, mut rhs: Vec<f64>, n: usize) -> Vec<f64> {
    for col in 0..n {
        let piv = (col..n).max_by(|&i, &j| m[i * n + col].abs().total_cmp(&m[j * n + col].abs())).unwrap();
        if m[piv * n + col].abs() < 1e-300 {
            continue;
        }
        if piv != col {
            for k in 0..n {
                m.swap(piv * n + k, col * n + k);
            }
            rhs.swap(piv, col);
        }
        let d = m[col * n + col];
        for i in col + 1..n {
            let f = m[i * n + col] / d;
            if f != 0.0 {
                for k in col..n {
                    m[i * n + k] -= f * m[col * n + k];
                }
                rhs[i] -= f * rhs[col];
            }
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let mut s = rhs[i];
        for k in i + 1..n {
            s -= m[i * n + k] * x[k];
        }
        let d = m[i * n + i];
        x[i] = if d.abs() > 1e-300 { s / d } else { 0.0 };
    }
    x
}

/// Rotation vectors (mm/yr at the surface per unit radius: `v(p) = u × p`) fitted to the
/// prescribed boundaries, regularised toward rolled priors `u0`, with no net rotation.
fn solve_motion(grid: &Grid, rows: &[(usize, u16, u16, DVec3, Kind)], u0: &[DVec3], area: &[f64]) -> Vec<DVec3> {
    let k = u0.len();
    let n = 3 * k;
    let mut m = vec![0.0; n * n];
    let mut rhs = vec![0.0; n];
    let mut total_w = 0.0;
    for &(a, own, other, normal, kind) in rows {
        let (t, wk) = target(kind);
        let w = wk * (grid.weight[a] as f64).sqrt();
        total_w += w;
        let g = grid.dir[a].cross(normal); // c = (u_own − u_other) · (p × n)
        let blocks = [(own as usize, 1.0), (other as usize, -1.0)];
        for &(bi, si) in &blocks {
            for r in 0..3 {
                rhs[3 * bi + r] += w * t * si * g[r];
                for &(bj, sj) in &blocks {
                    for c in 0..3 {
                        m[(3 * bi + r) * n + 3 * bj + c] += w * si * sj * g[r] * g[c];
                    }
                }
            }
        }
    }
    let total_w = total_w.max(1e-9);
    let lambda = 0.02 * total_w / k as f64;
    let mu = 10.0 * total_w;
    for i in 0..k {
        for r in 0..3 {
            m[(3 * i + r) * n + 3 * i + r] += lambda;
            rhs[3 * i + r] += lambda * u0[i][r];
            for j in 0..k {
                m[(3 * i + r) * n + 3 * j + r] += mu * area[i] * area[j];
            }
        }
    }
    let x = solve_dense(m, rhs, n);
    (0..k).map(|i| DVec3::new(x[3 * i], x[3 * i + 1], x[3 * i + 2])).collect()
}

/// Smooth 0..1 margin activity for the shelf widths (`relief::finish`), replacing the old
/// noise: 1 off active coasts, 0 off passive ones, blended over ~200 km.
pub fn activity(grid: &Grid, cls: &Classification) -> Vec<f32> {
    grid.coarse_field(48, |d| if cls.margin_at(d) == margin::ACTIVE { 1.0 } else { 0.0 })
}

pub fn build(grid: &Grid, p: &Params, crust: &Crust, cls: &Classification) -> (PlateModel, FitStats) {
    let r = p.radius_m;
    let tm = std::time::Instant::now();
    let lap = |what: &str| if std::env::var_os("RG_TIMING").is_some() { eprintln!("    fit {what}: {} ms", tm.elapsed().as_millis()) };
    let mm_per_yr = r * 1e-3; // rad/Myr -> mm/yr
    let lay = layout(cls);
    let mut rng = Rng::new(p.seed, STREAM_FIT);

    // --- on the classification grid: oceanic seeds and the motion solve.
    let cg = &cls.grid;
    let cgeo = geometry(cg, &cls.crust, cls, p.seed, r);
    let ccost = plate_cost(cg, p.seed);
    let cland = land_plates(cg, &cgeo, cls, &lay, &ccost);
    let cclaim = claims(cg, &cgeo, cls, &lay, &cland, &ccost);
    let seeds = ocean_seeds(cg, p, &cgeo, &cclaim, lay.n_cont, &mut rng);
    let cplate = partition(cg, &cgeo, &cclaim, &seeds, lay.n_cont, &ccost);
    let k = lay.n_cont + seeds.len();

    let mut area = vec![0.0f64; k];
    let mut cont_area = vec![0.0f64; k];
    for a in 0..cg.n {
        let pl = cplate[a] as usize;
        area[pl] += cg.weight[a] as f64;
        if cls.crust.cont[a] {
            cont_area[pl] += cg.weight[a] as f64;
        }
    }
    // Priors, rolled as before R3: oceanic plates a little faster (slab pull).
    let mut u0 = Vec::with_capacity(k);
    let mut density = Vec::with_capacity(k);
    for pl in 0..k {
        let continental = cont_area[pl] > 0.5 * area[pl];
        let deg = rng.range(0.15, 0.7) * if continental { 1.0 } else { 1.35 };
        u0.push(rng.unit_vec() * deg.to_radians() * mm_per_yr);
        density.push(rng.f64() as f32);
    }
    let cborders = borders(cg, &cplate, &cls.crust.cont, &cgeo.region);
    let rows: Vec<(usize, u16, u16, DVec3, Kind)> = cborders
        .iter()
        .filter_map(|b| prescribe(b, &cgeo, cls, &lay, &cls.crust.cont).map(|(kind, _)| (b.node, b.own, b.other, b.normal, kind)))
        .collect();
    let u = solve_motion(cg, &rows, &u0, &area);
    lap("classification grid + solve");

    let mut stats = FitStats { active_share: cls.active_share, plates: k as u32, ..Default::default() };
    for sw in &cls.seaways {
        stats.seaways[match sw.class {
            SeawayClass::Opening => 0,
            SeawayClass::Closing => 1,
            SeawayClass::OneSided => 2,
        }] += 1;
    }
    let (mut sw_res, mut sw_len, mut wrong, mut signed_len) = (0.0, 0.0, 0.0, 0.0);
    for &(a, own, other, normal, kind) in &rows {
        let (t, _) = target(kind);
        let c = (u[own as usize] - u[other as usize]).dot(cg.dir[a].cross(normal));
        let w = (cg.weight[a] as f64).sqrt();
        sw_res += w * (c - t) * (c - t);
        sw_len += w;
        if kind != Kind::Transform {
            signed_len += w;
            if c * t < 0.0 {
                wrong += w;
            }
        }
    }
    stats.rms_residual_mm_yr = (sw_res / sw_len.max(1e-12)).sqrt();
    stats.wrong_sign_share = wrong / signed_len.max(1e-12);

    // --- on the generation grid: the same plates with finer borders.
    let geo = geometry(grid, crust, cls, p.seed, r);
    lap("geometry");
    let cost = plate_cost(grid, p.seed);
    let land = land_plates(grid, &geo, cls, &lay, &cost);
    lap("land plates");
    let claim = claims(grid, &geo, cls, &lay, &land, &cost);
    lap("claims");
    let mut plate_of = partition(grid, &geo, &claim, &seeds, lay.n_cont, &cost);
    lap("partition");

    // Compact away plates the generation grid didn't reproduce (rare; tiny pieces).
    let mut used = vec![false; k];
    for &pl in &plate_of {
        used[pl as usize] = true;
    }
    let mut remap = vec![NONE; k];
    let mut kept = Vec::new();
    for pl in 0..k {
        if used[pl] {
            remap[pl] = kept.len() as u16;
            kept.push(pl);
        }
    }
    for pl in plate_of.iter_mut() {
        *pl = remap[*pl as usize];
    }
    let mut fine_area = vec![0.0f64; kept.len()];
    let mut fine_cont = vec![0.0f64; kept.len()];
    for a in 0..grid.n {
        fine_area[plate_of[a] as usize] += grid.weight[a] as f64;
        if crust.cont[a] {
            fine_cont[plate_of[a] as usize] += grid.weight[a] as f64;
        }
    }
    let min_omega = 0.03f64.to_radians();
    let plates: Vec<Plate> = kept
        .iter()
        .enumerate()
        .map(|(i, &pl)| {
            let w = u[pl] / mm_per_yr; // rad/Myr
            let (pole, omega) = if w.length() > 1e-12 { (w.normalize(), w.length()) } else { (u0[pl].normalize(), 0.0) };
            Plate { euler_pole: pole, omega: omega.clamp(min_omega, 0.045), continental: fine_cont[i] > 0.5 * fine_area[i] }
        })
        .collect();
    let density: Vec<f32> = kept.iter().map(|&pl| density[pl]).collect();
    // The fine layout, after compaction (continental plates keep their order).
    let fine_lay = Layout {
        base: lay.base.iter().map(|&b| remap[b as usize]).collect(),
        n_cont: kept.iter().filter(|&&pl| pl < lay.n_cont).count(),
        plate_land: kept.iter().filter(|&&pl| pl < lay.n_cont).map(|&pl| lay.plate_land[pl]).collect(),
    };

    // --- boundaries: prescribed kinds, rates from the solved motion.
    let n = grid.n;
    let s_strength = salt(p.seed, crate::plates::STREAM_PLATES + 1);
    let mut boundary = vec![0u8; n];
    let mut convergence = vec![0f32; n];
    let mut bnodes = Vec::new();
    for b in borders(grid, &plate_of, &crust.cont, &geo.region) {
        let a = b.node;
        let pa = grid.dir[a];
        let (po, pt) = (&plates[b.own as usize], &plates[b.other as usize]);
        let rel = velocity(po, pa) - velocity(pt, pa);
        let speed_mm = rel.length() * mm_per_yr;
        let c = rel.dot(b.normal) * mm_per_yr; // > 0: converging
        let (kind, mut sub_plate, prescribed) = match prescribe(&b, &geo, cls, &fine_lay, &crust.cont) {
            Some((kind, sub)) => (kind, sub, true),
            None => {
                let ratio = if speed_mm > 0.0 { c / speed_mm } else { 0.0 };
                let kind = if speed_mm < 3.0 || ratio.abs() < 0.4 {
                    Kind::Transform
                } else if ratio > 0.0 {
                    if b.a_cont && b.b_cont { Kind::Collision } else { Kind::Subduction }
                } else if b.a_cont && b.b_cont {
                    Kind::Rift
                } else {
                    Kind::Ridge
                };
                (kind, b.own, false)
            }
        };
        let rate = match kind {
            Kind::Subduction => c.max(MIN_SUBDUCTION),
            Kind::Collision => c.max(MIN_COLLISION),
            Kind::Ridge => (-c).max(MIN_RIDGE),
            Kind::Rift => (-c).max(MIN_RIFT),
            Kind::Transform => c.abs(),
        };
        let mut over_cont = false;
        if kind == Kind::Subduction {
            if prescribed {
                let over = if sub_plate == b.own { b.other } else { b.own };
                over_cont = (over as usize) < fine_lay.n_cont;
            } else {
                let own_dives =
                    if b.a_cont != b.b_cont { !b.a_cont } else { density[b.own as usize] > density[b.other as usize] };
                sub_plate = if own_dives { b.own } else { b.other };
                over_cont = if own_dives { b.b_cont } else { b.a_cont };
            }
        }
        boundary[a] = match kind {
            Kind::Subduction | Kind::Collision => crate::boundary_class::CONVERGENT,
            Kind::Ridge | Kind::Rift => crate::boundary_class::DIVERGENT,
            Kind::Transform => crate::boundary_class::TRANSFORM,
        };
        if matches!(kind, Kind::Subduction | Kind::Collision) {
            convergence[a] = rate as f32;
        }
        let rate_factor = match kind {
            Kind::Ridge | Kind::Rift => (rate / 2.0 / 40.0).clamp(0.3, 1.2),
            _ => (rate / 55.0).clamp(0.3, 1.4),
        };
        let along_strike = 0.55 + 0.9 * fbm01(pa, s_strength, 4.0, 3);
        bnodes.push(BoundaryNode {
            node: a as u32,
            kind,
            strength: (rate_factor * along_strike) as f32,
            speed_mm_yr: speed_mm.max(rate) as f32,
            sub_plate,
            over_cont,
        });
    }
    lap("boundaries");
    let margin: Vec<u8> =
        (0..n).map(|a| if (geo.coast_dist_m[a] as f64) < MARGIN_BAND_M { geo.margin[a] } else { margin::NONE }).collect();
    stats.plates = plates.len() as u32;
    (PlateModel { plate_of, plates, boundary, bnodes, convergence, margin }, stats)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Two plates meeting along part of the equator, normal pointing north from plate 0 to plate 1. A
    /// prescribed subduction boundary must come out converging, a ridge diverging.
    #[test]
    fn solve_matches_prescribed_signs() {
        let g = Grid::new(16);
        // A 60° arc of the equator: around a whole great circle no rigid rotation converges
        // everywhere (the rate is sinusoidal along strike; design §4).
        let eq: Vec<usize> = (0..g.n).filter(|&a| g.dir[a].z.abs() < 0.05 && g.dir[a].x > 0.87).collect();
        for (kind, sign) in [(Kind::Subduction, 1.0), (Kind::Ridge, -1.0)] {
            let rows: Vec<(usize, u16, u16, DVec3, Kind)> = eq
                .iter()
                .map(|&a| {
                    let p = g.dir[a];
                    let north = (DVec3::Z - p * p.z).normalize();
                    (a, 0u16, 1u16, north, kind)
                })
                .collect();
            let u0 = vec![DVec3::new(0.0, 0.0, 1.0), DVec3::new(0.0, 1.0, 0.0)];
            let u = solve_motion(&g, &rows, &u0, &[0.5, 0.5]);
            let mean: f64 = rows.iter().map(|&(a, o, t, n, _)| (u[o as usize] - u[t as usize]).dot(g.dir[a].cross(n))).sum::<f64>()
                / rows.len() as f64;
            assert!(mean * sign > 10.0, "{kind:?}: mean rate {mean}");
        }
    }
}
