//! The drained detail network (`docs/design/drained-detail-network.md`, v1.1).
//!
//! Land height is the minimum of a cap `U` (the drained triangle-split macro surface) and a
//! family of **valley wedges**, one per channel segment:
//!
//! `W_s(p) = z_s(t) + H_s(d)`: channel elevation at the nearest point of the segment plus a
//! strictly increasing hillslope profile of the distance to it.
//!
//! Channels form a hierarchy. Level 0 is the macro river network (node → receiver). Level k
//! has one jittered point per cell of a lattice 3^k times finer than the macro grid. Each
//! point drains to the nearest point of a coarser channel, and its channel rises from that
//! junction following the stream-power law. Because a local minimum of a min of functions is
//! a minimum of the active one, and every wedge's only minimum is its junction (where the
//! parent wedge continues downhill), the surface drains to the sea by construction.
//!
//! Everything a channel is made of (its point, junction, profile and hillslope) depends only
//! on its own lattice cell, so any two patches agree exactly: no seams between tiles.
//!
//! Evaluation is per **patch**: gather the channels that can reach a region (a tile, a probe
//! box), bin them, then evaluate each point against its bin.

use crate::network::{MacroNet, Node, Segment};
use glam::DVec3;
use rg_geo::lattice::{self, Cell};
use rg_geo::{CubeFace, CubeField};
use std::collections::HashMap;
use std::sync::Arc;

/// Stream-power concavity θ = m/n (n = 1), as in the macro solve.
const THETA: f64 = 0.45;
/// k_sn at which hillslopes are ~63% of the way to threshold, m^0.9 (Ouimet et al. 2009).
const KSN_SATURATE: f64 = 100.0;
/// Threshold hillslope gradient (tan 33°): soil-mantled slopes stop steepening here
/// (Montgomery & Brandon 2002; Roering et al. 2007).
const S_CRIT: f64 = 0.649;
/// Steepest channel a profile may have (colluvial/debris-flow channels), m/m.
const S_CHANNEL_MAX: f64 = 0.5;
/// Hilltop slope as a share of the valley-wall slope: slowly eroding (convex divides) and
/// fast-eroding (near-planar threshold slopes) end members.
const TOP_SLOPE_SHARE: f64 = 0.3;
const TOP_SHARE_FAST: f64 = 0.9;
/// Smallest rise per metre along any channel (keeps every profile strictly monotone, and
/// above f32 resolution over a cell).
const EPS_SLOPE: f64 = 1.0e-3;
/// Lattice jitter, cell widths (< 0.5).
const JITTER: f64 = 0.4;
/// Finest level spacing we ever build, m.
const FINEST_M: f64 = 80.0;
/// Junction angles are drawn per channel in this range (degrees). The shift downstream is
/// d / tan α, at most d / tan 44° ≈ 1.04 d, inside the margin derivation's 1.5× length bound.
const JUNCTION_ANGLE_MIN_DEG: f64 = 44.0;
const JUNCTION_ANGLE_MAX_DEG: f64 = 80.0;
/// Share of lattice cells that get a channel (the rest are thinned out).
const KEEP_FRACTION: f64 = 0.85;
/// No channel closer to its parent than this many level spacings (a stub).
const MIN_ATTACH: f64 = 0.35;
/// Floodplain width per √(km²) of drainage (see `floor_width`), m.
const FLOOR_M_PER_SQRT_KM2: f64 = 30.0;
/// Where channel heads aim relative to the macro surface, as a share of the ridge
/// allowance (heads at T + share·R; ridges rise up to T + R). Calibrated so the detail is
/// mean-neutral against the macro, which already carries the hypsometry (probe_hypsometry).
const HEAD_SHARE: f64 = -0.33;
/// Macro height above which the ridge allowance stops growing, m.
const RIDGE_FULL_M: f64 = 3000.0;
/// Ridge allowance per metre of macro height (see `DetailNet::ridge`).
const RIDGE_ALPHA: f64 = 0.3;
/// Junction search radius, in parent spacings. Tributaries meet the nearest coarser channel
/// well within one parent spacing; a head with nothing in range simply has no channel.
/// Sets the gather margins (see `DetailNet::margin`).
const ATTACH_SEARCH: f64 = 1.2;
/// Sibling attach: search radius and the longest channel that may serve as (or be) a
/// sibling-attached channel, in level spacings. Bounds the rank steps in the margins.
const SIBLING_SEARCH: f64 = 1.5;
const SIBLING_MAX_LEN: f64 = 2.5;
/// Steepest the surface can be within reach of its active wedge (m/m): bounds the minimum
/// over a block or cell from its centre value (exact pruning).
const MAX_GRADIENT: f64 = 1.5;
/// Profiling counters (native, RG_NET_PROF): candidate lists built, total length.
#[cfg(not(target_arch = "wasm32"))]
static CAND_STATS: [std::sync::atomic::AtomicUsize; 2] = [std::sync::atomic::AtomicUsize::new(0), std::sync::atomic::AtomicUsize::new(0)];
/// Reach wall length, spacings. The wall is sized per channel to climb from the bed to an
/// upper bound of the local cap by `reach_end` (gentle, not a 10 km cliff: walls that won
/// briefly drew dark crescent scarps, s01). Wedge reach in total: 2 + WALL_LEN spacings.
const WALL_LEN: f64 = 1.5;
const REACH_TOTAL: f64 = 2.0 + WALL_LEN;
/// Hack's exponent (L ∝ A^h, h ≈ 0.57; Hack 1957), for the along-channel area decline.
const HACK_H: f64 = 0.57;
/// Lake depth (m) at a node above which it counts as under a lake.
const LAKE_MIN_DEPTH_M: f64 = 0.5;

/// The water level of a lake whose nodes include a corner of the macro cell around `p`
/// (the highest such level), if any.
fn lake_level_near(h: &CubeField<f32>, depth: &CubeField<f32>, p: DVec3) -> Option<f64> {
    let (face, u, v) = CubeFace::from_dir(p);
    let r = h.res as f64;
    let (x, y) = (((u + 1.0) * 0.5 * r).clamp(0.0, r), ((v + 1.0) * 0.5 * r).clamp(0.0, r));
    let (i0, j0) = ((x.floor() as u32).min(h.res - 1), (y.floor() as u32).min(h.res - 1));
    let mut level: Option<f64> = None;
    for (di, dj) in [(0, 0), (1, 0), (0, 1), (1, 1)] {
        let k = h.index(face, i0 + di, j0 + dj);
        let d = depth.data[k] as f64;
        if d > LAKE_MIN_DEPTH_M {
            let l = h.data[k] as f64 + d;
            level = Some(level.map_or(l, |m: f64| m.max(l)));
        }
    }
    level
}
/// Pieces per child channel (meander).
const CHILD_PIECES: usize = 3;
/// Pieces per curved level-0 edge.
const CURVE_PIECES: usize = 4;
/// Fade-out terms (see `wedge`): large enough that a faded channel never wins.
/// A faded channel's bed climbs this steeply (m/m) from its junction, and its walls this
/// steeply, so a channel at presence 0 never wins anywhere but its junction.
const FADE_RAMP_SLOPE: f64 = 2.0;
const FADE_SLOPE: f64 = 50.0;
/// Most ancestors a channel can have (levels ≥ 1 above it).
const MAX_ANCESTORS: usize = 7;

/// One shard of the channel memo: (level, cell) → channel (None = no channel there).
type MemoShard = std::sync::Mutex<HashMap<(u8, Cell), (u8, Vec<Channel>)>>;
type CandShard = std::sync::Mutex<HashMap<(u8, Cell), Arc<Vec<Channel>>>>;

/// A parent candidate: its wedge at the head, the distance and parameter there.
#[derive(Clone, Copy)]
struct Cand {
    wq: f64,
    d: f64,
    t: f64,
    ch: Channel,
}

/// Cell-local facts for building a channel (see `DetailNet::prelim`).
struct Prelim {
    q: DVec3,
    u: [f64; 2],
    t_q: f64,
    coarse: f64,
    parent: Cand,
    rank: u8,
}

/// One channel segment with everything its wedge needs. `down` is the junction (or the
/// macro receiver), `up` the head.
#[derive(Clone, Copy, Debug)]
pub struct Channel {
    pub level: u8,
    pub up: DVec3,
    pub down: DVec3,
    pub z_up: f64,
    pub z_down: f64,
    /// Normalised channel steepness k_sn = S·A^θ (m^(2θ)), inherited down the hierarchy.
    pub ksn: f64,
    /// Valley-wall gradient, hilltop gradient, convexity width, floor width (m).
    pub s_wall: f64,
    pub s_top: f64,
    pub w: f64,
    pub floor: f64,
    /// Beyond `reach` (m) the wedge rises steeply (`wall_k` per m²), so it never wins there.
    pub reach: f64,
    pub wall_k: f64,
    /// Where the wall has climbed past any land height the wedge could compete with
    /// (reach + WALL_LEN spacings): beyond it the wedge never wins.
    pub reach_end: f64,
    /// Effective level spacing for LOD presence (m): the level spacing × a hashed factor in
    /// [0.75, 1.33], so siblings fade in at different zooms. A parent's (≥ 3/1.33 × the
    /// child's level spacing) always exceeds its child's, so a parent is never less present.
    pub lod_spacing: f64,
    /// Length of this piece, m.
    pub len_m: f64,
    /// Arc length from the channel's junction to this piece's down end, m (the fade ramp
    /// runs along the whole channel, so it doesn't restart at each piece).
    pub s0: f64,
    /// Arc length of the whole channel (all its pieces), m.
    pub chan_len: f64,
    /// Fading ancestors (levels ≥ 1), each (LOD spacing, distance from its junction to where
    /// this lineage attaches, m). A faded ancestor's bed is raised by (1 − φ)·ramp·distance
    /// there, and descendants carry that offset, so every junction sits on its parent's
    /// effective bed at every zoom (otherwise fading made junction pits ~100 m deep, s01).
    pub anc: [[f32; 2]; MAX_ANCESTORS],
    pub n_anc: u8,
}

impl Channel {
    /// (distance m, parameter t from down = 0 to up = 1) of `p` to the segment.
    #[inline]
    fn nearest(&self, p: DVec3, radius_m: f64) -> (f64, f64) {
        let ab = self.up - self.down;
        let l2 = ab.length_squared();
        let t = if l2 > 0.0 { ((p - self.down).dot(ab) / l2).clamp(0.0, 1.0) } else { 0.0 };
        ((p - (self.down + ab * t)).length() * radius_m, t)
    }

    #[inline]
    fn z_at(&self, t: f64) -> f64 {
        self.z_down + (self.z_up - self.z_down) * t
    }

    /// Lowest channel elevation (the junction).
    pub fn z_min(&self) -> f64 {
        self.z_down
    }

    /// Presence at band-limit `cutoff`. Level 0 never fades: consecutive macro segments
    /// aren't parent and child, so a ramp would break their continuity (band-limited floors
    /// handle coarse views instead).
    #[inline]
    fn presence(&self, cutoff: f64) -> f64 {
        if self.level == 0 { 1.0 } else { presence_of(self.lod_spacing, cutoff) }
    }

    /// Fade offset inherited from ancestors at band-limit `cutoff` (m, ≥ 0).
    #[inline]
    fn inherited_offset(&self, cutoff: f64) -> f64 {
        let mut off = 0.0;
        for a in &self.anc[..self.n_anc as usize] {
            let phi = presence_of(a[0] as f64, cutoff);
            off += (1.0 - phi) * FADE_RAMP_SLOPE * a[1] as f64;
        }
        off
    }

    /// The wedge at distance `d`, parameter `t`, faded in by `phi` (1 = fully present), with
    /// the valley floor at least `min_floor` wide (band-limiting: a level of detail can't
    /// show a floor narrower than it resolves, and sampling a sharper V makes grid pits).
    /// Every term is strictly increasing in `d` and in `t`, at any `phi`; `min_floor` is
    /// fixed per evaluation, so each wedge's shape is fixed as the proof requires.
    #[inline]
    fn wedge(&self, d: f64, t: f64, phi: f64, min_floor: f64, offset: f64) -> f64 {
        let f = self.floor.max(min_floor);
        let g = floor_hinge(d, f);
        let mut h = self.z_at(t) + offset + (self.s_wall - self.s_top) * self.w * (g / self.w).tanh() + self.s_top * g;
        let over = d - self.reach;
        if over > 0.0 {
            h += self.wall_k * over * over;
        }
        if phi < 1.0 {
            h += (1.0 - phi) * (FADE_RAMP_SLOPE * (self.s0 + self.len_m * t) + FADE_SLOPE * g);
        }
        h
    }
}

/// How close hillslopes are to threshold (0 = diffusive lowland, 1 = threshold): mean
/// hillslope gradient rises with k_sn and saturates near 100–150 m^0.9 (Ouimet, Whipple &
/// Granger 2009, "Beyond threshold hillslopes", Geology; θ_ref = 0.45).
fn threshold_fraction(ksn: f64) -> f64 {
    1.0 - (-ksn / tune("RG_KSN_SAT", KSN_SATURATE)).exp()
}

/// Distance up the valley side: 0 on a flat floor (|d| ≤ f/2), a smooth quadratic hinge up
/// to 1.5 f, then linear (C1 everywhere). Flat, not rounded: a rounded floor varies by metres
/// across the thalweg, which made grid pits on gentle channels; on a flat floor the wedge is
/// exactly the channel elevation, strictly decreasing downstream. Drainage holds because a
/// parent's floor is never narrower than its child's (ℓ/8 grows with level; the band-limit
/// minimum is shared), so a child's floor ends inside its parent's.
#[inline]
fn floor_hinge(d: f64, f: f64) -> f64 {
    let (a, b) = (0.5 * f, 1.5 * f);
    if d <= a {
        0.0
    } else if d < b {
        (d - a) * (d - a) / (2.0 * (b - a))
    } else {
        d - f
    }
}

/// Valley-floor width for a drainage area (km²), m: floodplain width grows like √A, because
/// channel width grows like √Q (hydraulic geometry, Leopold & Maddock 1953) and Q ∝ A.
/// 30 m·√(A/km²) is a calibration: ~300 m at 100 km², ~3 km at 10⁴ km², a few metres for
/// first-order valleys (which are V-shaped). It grows with level, as drainage requires.
fn floor_width(area_km2: f64) -> f64 {
    FLOOR_M_PER_SQRT_KM2 * area_km2.max(0.0).sqrt()
}

/// Valley-wall gradient from k_sn.
fn wall_slope(ksn: f64) -> f64 {
    S_CRIT * threshold_fraction(ksn)
}

/// Hilltop gradient as a share of the wall gradient: convex, soil-mantled hills in slowly
/// eroding terrain (~0.3); nearly planar threshold slopes with a narrow convex cap where
/// erosion is fast (Roering et al. 2007, "Functional relationships between denudation and
/// hillslope form and relief").
fn top_share(ksn: f64) -> f64 {
    TOP_SLOPE_SHARE + (tune("RG_TOP", TOP_SHARE_FAST) - TOP_SLOPE_SHARE) * threshold_fraction(ksn)
}

/// Presence of a channel with LOD spacing `lod_spacing` at band-limit `cutoff`.
#[inline]
fn presence_of(lod_spacing: f64, cutoff: f64) -> f64 {
    smoothstep(2.0 * cutoff, 4.0 * cutoff, lod_spacing)
}

/// Linear smoothstep.
fn smoothstep(e0: f64, e1: f64, x: f64) -> f64 {
    let t = ((x - e0) / (e1 - e0)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// The network for one world: macro hydrology plus the lattice hierarchy.
#[derive(Clone)]
pub struct DetailNet {
    pub macro_net: MacroNet,
    pub radius_m: f64,
    salt: u64,
    /// Spacing of each level, m (level 0 = macro cell).
    pub spacing: Vec<f64>,
    /// Lattice cells per face edge of each level (level 0 = macro res).
    cells: Vec<u32>,
    /// Memo of built channels per (level, cell). Only a cache: a channel depends on its cell
    /// alone (given complete coarser sets, which the gather margins guarantee), so patches
    /// built in any order agree. Sharded to keep worker threads from queueing.
    memo: Arc<[MemoShard; 16]>,
    /// Memo of coarser-candidate lists per (level, parent-lattice cell) (see `candidates`).
    cand_memo: Arc<[CandShard; 16]>,
    ridge_alpha: f64,
    /// Erodibility K of the macro solve, m^(1−2m)/yr (k_sn = U/K).
    k_erodibility: f64,
    /// Upper bound of the cap's gradient (m/m), from the steepest macro edge.
    cap_gradient: f64,
}

impl DetailNet {
    pub fn new(macro_net: MacroNet, radius_m: f64, seed: u64) -> Self {
        let res = macro_net.res();
        let l0 = std::f64::consts::FRAC_PI_2 * radius_m / res as f64;
        let (mut spacing, mut cells) = (vec![l0], vec![res]);
        while spacing.last().unwrap() / 3.0 >= FINEST_M && spacing.len() < 8 {
            spacing.push(spacing.last().unwrap() / 3.0);
            cells.push(cells.last().unwrap() * 3);
        }
        let memo = Arc::new(std::array::from_fn(|_| std::sync::Mutex::new(HashMap::new())));
        let cand_memo = Arc::new(std::array::from_fn(|_| std::sync::Mutex::new(HashMap::new())));
        // Steepest macro edge (axis and diagonal neighbours): a triangle facet's gradient is at
        // most √2 × its steepest edge; the cap multiplies heights by ≤ 1 + α.
        let h = &macro_net.height;
        let mut steepest = 0.0f64;
        for (k, face, i, j, d) in h.iter_samples() {
            for (di, dj) in [(1i32, 0i32), (0, 1), (1, 1), (1, -1)] {
                let (ii, jj) = (i as i32 + di, j as i32 + dj);
                if ii < 0 || jj < 0 || ii > h.res as i32 || jj > h.res as i32 {
                    continue;
                }
                let m = h.index(face, ii as u32, jj as u32);
                let dist = (h.sample_dir(face, ii as u32, jj as u32) - d).length() * radius_m;
                steepest = steepest.max(((h.data[m] - h.data[k]) as f64).abs() / dist.max(1.0));
            }
        }
        let ridge_alpha = tune("RG_RIDGE", RIDGE_ALPHA);
        Self { macro_net, radius_m, salt: seed ^ 0xD371_0A11_C4A7_0001, spacing, cells, memo, cand_memo, ridge_alpha, k_erodibility: 1.0e-5, cap_gradient: steepest * std::f64::consts::SQRT_2 * (1.0 + RIDGE_ALPHA.max(0.0) + 0.5) + 0.01 }
    }

    /// Set the macro solve's erodibility K (default 1e-5), used for k_sn = U/K.
    pub fn with_erodibility(mut self, k: f64) -> Self {
        self.k_erodibility = k;
        self
    }

    /// Coarser channels (levels < k) that can be the lowest coarser wedge somewhere in the
    /// level-(k−1) lattice cell containing `q`. Exact pruning, as for patch blocks: anywhere
    /// in the cell a wedge is ≥ z_down + H(max(0, d_c − r_cell)), and the minimum is ≤ the best
    /// value at the centre + MAX_GRADIENT·r_cell. The list depends only on the cell (given
    /// complete coarser sets), so it is memoised across patches.
    fn candidates(&self, k: usize, q: DVec3, coarser: &[Binned]) -> Arc<Vec<Channel>> {
        self.cell_candidates(k, lattice::cell_of(q, self.cells[k - 1]), coarser)
    }

    /// The list for level-(k−1) cell `cell`, built hierarchically: a channel of levels < k−1
    /// that can be the minimum somewhere in this cell is also a minimum (over levels < k−1)
    /// somewhere in the enclosing level-(k−2) cell, so it is in that cell's list. Add the
    /// level-(k−1) channels near the cell and prune. Exact, and no rescan of every level.
    fn cell_candidates(&self, k: usize, cell: Cell, coarser: &[Binned]) -> Arc<Vec<Channel>> {
        let key = (k as u8, cell);
        let shard = &self.cand_memo[(lattice::hash(k as u64 ^ 0xCA, cell) % 16) as usize];
        if let Some(v) = shard.lock().ok().and_then(|m| m.get(&key).cloned()) {
            return v;
        }
        let n = self.cells[k - 1];
        let pc = lattice::centre(cell, n);
        // Cell half-diagonal (cells are up to ~1.4 spacings across near face corners).
        let r_cell = self.spacing[k - 1];
        let mut all: Vec<Channel> = if k >= 2 {
            let parent_cell = lattice::cell_of(pc, self.cells[k - 2]);
            self.cell_candidates(k - 1, parent_cell, coarser).as_ref().clone()
        } else {
            Vec::new()
        };
        let base = all.len();
        let b = &coarser[k - 1];
        b.for_each_near(pc, REACH_TOTAL * b.spacing_m + b.max_len + r_cell, |ch| all.push(*ch));
        // Only the new level's entries can repeat (seam windows); dedupe those.
        let mut fresh = all.split_off(base);
        fresh.sort_by(|a, b| a.up.x.total_cmp(&b.up.x).then(a.up.y.total_cmp(&b.up.y)).then(a.up.z.total_cmp(&b.up.z)).then(a.down.x.total_cmp(&b.down.x)));
        fresh.dedup_by(|a, b| a.up == b.up && a.down == b.down);
        all.extend(fresh);
        let t = self.macro_tri(pc);
        let mut ub = t + self.ridge(t);
        for ch in &all {
            let (d, tt) = ch.nearest(pc, self.radius_m);
            ub = ub.min(ch.wedge(d, tt, 1.0, 0.0, 0.0));
        }
        let ub = ub + MAX_GRADIENT * r_cell;
        all.retain(|ch| {
            let (d, _) = ch.nearest(pc, self.radius_m);
            ch.wedge((d - r_cell).max(0.0), 0.0, 1.0, 0.0, 0.0) <= ub
        });
        #[cfg(not(target_arch = "wasm32"))]
        {
            CAND_STATS[0].fetch_add(1, std::sync::atomic::Ordering::Relaxed);
            CAND_STATS[1].fetch_add(all.len(), std::sync::atomic::Ordering::Relaxed);
        }
        let v = Arc::new(all);
        if let Ok(mut m) = shard.lock() {
            // Budget by channels held (~200 B each): ≤ ~60k per shard, ~200 MB in all.
            let held: usize = m.values().map(|l| l.len()).sum();
            if held + v.len() > 60_000 {
                m.clear();
            }
            m.insert(key, v.clone());
        }
        v
    }

    fn memo_put(&self, k: usize, c: Cell, v: (u8, Vec<Channel>)) {
        let shard = &self.memo[(lattice::hash(k as u64, c) % 16) as usize];
        if let Ok(mut m) = shard.lock() {
            // Bound memory: drop the shard when it grows large (it's only a cache).
            if m.len() > 50_000 {
                m.clear();
            }
            m.insert((k as u8, c), v);
        }
    }

    /// The macro surface T: drained triangle split.
    pub fn macro_tri(&self, p: DVec3) -> f64 {
        self.macro_net.height.sample_tri(p)
    }

    /// Ridge allowance above the macro surface, a function of macro height only (so the cap
    /// stays a strictly increasing function of a drained surface, hence drained). Local
    /// relief grows with mean elevation in orogens (Ahnert 1970).
    pub fn ridge(&self, t: f64) -> f64 {
        self.ridge_alpha * t.clamp(0.0, RIDGE_FULL_M)
    }

    /// Wall stiffness so the wedge exceeds any cap value it could meet by `reach_end`:
    /// (upper bound of the cap within reach − bed) / wall length².
    fn wall_for(&self, head: DVec3, len: f64, reach: f64, spacing: f64, z_down: f64) -> (f64, f64) {
        let wall = WALL_LEN * spacing;
        let reach_end = reach + wall;
        let cap_max = self.cap(head) + self.cap_gradient * (len + reach_end);
        (((cap_max - z_down).max(1.0) + 1.0) / (wall * wall), reach_end)
    }

    /// The cap U = T + R(T).
    pub fn cap(&self, p: DVec3) -> f64 {
        let t = self.macro_tri(p);
        t + self.ridge(t)
    }

    /// How present level `k` is at band-limit `cutoff_m` (0 = absent, 1 = full): a level of
    /// spacing ℓ needs about ℓ/4 resolution (its convexity width).
    /// (Upper bound over channels: their LOD spacing is up to 1.33 × the level's.)
    pub fn presence(&self, k: usize, cutoff_m: f64) -> f64 {
        smoothstep(2.0 * cutoff_m, 4.0 * cutoff_m, 1.33 * self.spacing[k])
    }

    fn level0_channel(&self, s: &Segment) -> Channel {
        let l0 = self.spacing[0];
        let len = (s.up - s.down).length() * self.radius_m;
        // k_sn: at stream-power steady state k_sn = U/K exactly (n = 1), from the stored
        // uplift. Fallback (older files): S·A^θ from the macro edge, inflated on lowlands
        // because macro heights are hypsometry-remapped, not physical.
        let ksn = match s.uplift_mm_yr {
            Some(u) => (u * 1.0e-3 / self.k_erodibility).max(0.1),
            None => {
                let slope = ((s.z_up - s.z_down) / len.max(1.0)).max(1.0e-5);
                let area = (s.area_km2 * 1.0e6).max(l0 * l0);
                slope * area.powf(THETA)
            }
        };
        let s_wall = wall_slope(ksn).max(0.002);
        // Floodplains in rising mountains narrow to gorges.
        let floor = floor_width(s.area_km2.max(l0 * l0 * 1.0e-6)).min(0.25 * l0) * (1.0 - 0.8 * threshold_fraction(ksn));
        Channel {
            level: 0,
            up: s.up,
            down: s.down,
            z_up: s.z_up,
            z_down: s.z_down,
            ksn,
            s_wall,
            s_top: top_share(ksn) * s_wall,
            w: 0.25 * l0,
            floor,
            reach: 2.0 * l0,
            wall_k: self.wall_for(s.up, len, 2.0 * l0, l0, s.z_down).0,
            reach_end: self.wall_for(s.up, len, 2.0 * l0, l0, s.z_down).1,
            lod_spacing: l0,
            len_m: len,
            s0: 0.0,
            chan_len: len,
            anc: [[0.0; 2]; MAX_ANCESTORS],
            n_anc: 0,
        }
    }

    /// How far (m) around a region's edge level-`k` channels must be gathered so that every
    /// channel within reach is present, and so that finer levels find all their junctions.
    ///
    /// With search σ = ATTACH_SEARCH parent spacings (3σ ℓ_k), a channel is at most 1.5× its
    /// search long (junction shift included): 4.5σ ℓ_k. Reach is 3 ℓ_k, so heads within
    /// (3 + 4.5σ) ℓ_k matter. Junction search from a head gathered at margin c ℓ_k needs the
    /// coarser level out to (c + 3σ + 13.5σ) ℓ_k ≤ 3c ℓ_k, i.e. c ≥ 8.25σ (and ≥ 3 + 4.5σ).
    /// Level 0 (macro, edges ≤ √2 ℓ_0) needs ≥ c/3 + σ + 1.5 spacings.
    fn margin(&self, k: usize) -> f64 {
        // With sibling chaining, rank-0 cells are gathered 2·(SIBLING_SEARCH + SIBLING_MAX_LEN)
        // = 8 ℓ_k beyond the margin, and their prelim reads the coarser surface: coarser heads
        // within (3 + 4.5σ)·3 ℓ_k of them. So 3c ≥ c + 8 + (3 + 4.5σ)·3, i.e. c ≥ 16.6 (σ 1.2).
        let rank_extra = 2.0 * (SIBLING_SEARCH + SIBLING_MAX_LEN);
        // + one parent spacing (3 ℓ_k): candidate lists are read from the parent cell's centre.
        // + one spacing for the meander bulge (≤ 0.15 × a channel ≤ 5.4 ℓ_k).
        let c = ((rank_extra + (REACH_TOTAL + 4.5 * ATTACH_SEARCH) * 3.0 + 3.0) / 2.0).max(8.25 * ATTACH_SEARCH) + 1.5;
        // Level 0: level-1 rank-0 heads sit at (c + rank_extra)/3 ℓ_0 and read level-0 wedges
        // whose heads are within 3 ℓ_0 + an edge (≤ 1.5 ℓ_0) of them.
        if k == 0 { ((c + rank_extra) / 3.0 + REACH_TOTAL + 1.5 + 1.0 + 0.5) * self.spacing[0] } else { c * self.spacing[k] }
    }

    /// The cell-local facts that decide a level-k cell's channel: its point, the coarser
    /// parent candidate, the coarser surface there and its rank. `None`: no channel.
    fn prelim(&self, k: usize, c: Cell, coarser: &[Binned]) -> Option<Prelim> {
        let lk = self.spacing[k];
        let q = lattice::jittered(c, self.cells[k], self.salt.wrapping_add(k as u64), JITTER);
        // Per-cell randomness (three uniforms) breaks the lattice's regularity.
        let hsh = lattice::hash(self.salt ^ 0x5EED_0000_0000_0000 ^ k as u64, c);
        let u = |i: u32| ((hsh >> (i * 21)) & 0x1F_FFFF) as f64 / 0x1F_FFFF as f64;
        if u(0) > tune("RG_KEEP", KEEP_FRACTION) {
            return None; // thinning: drainage density varies; dropping a channel is always safe
        }
        let t_q = self.macro_tri(q);
        if t_q <= 0.0 {
            return None;
        }
        // The coarser parent: the lowest coarser wedge at q within the search radius (q
        // drains to the hillslope it sits on; the *nearest* channel can be buried). Both it
        // and the coarser surface come from the exact candidate list of q's parent cell.
        let search = ATTACH_SEARCH * self.spacing[k - 1];
        let mut best: Option<Cand> = None;
        let u_q = t_q + self.ridge(t_q);
        let mut coarse = u_q;
        for ch in self.candidates(k, q, coarser).iter() {
            let (d, t) = ch.nearest(q, self.radius_m);
            let wq = ch.wedge(d, t, 1.0, 0.0, 0.0);
            coarse = coarse.min(wq);
            if d > search || d < 1.0 || best.as_ref().is_some_and(|b| b.wq <= wq) {
                continue;
            }
            best = Some(Cand { wq, d, t, ch: *ch });
        }
        let parent = best?;
        // Rank from hillslope position: cells low on the coarser hillside (near its bed)
        // are collectors (rank 0, drain to the coarser channel); higher cells may also drain
        // to lower-rank siblings. That gives side-branching trees instead of combs.
        let s_wall = wall_slope(parent.ch.ksn).max(0.002);
        let rel = (coarse - parent.ch.z_at(parent.t)).max(0.0) / (s_wall * lk);
        let jitter = (lattice::hash(self.salt ^ 0x0A4C_0000_0000_0000 ^ k as u64, c) >> 11) as f64 / (1u64 << 53) as f64;
        let rank = ((2.0 * rel + jitter - 0.5).floor().clamp(0.0, 2.0)) as u8;
        Some(Prelim { q, u: [u(1), u(2)], t_q, coarse, parent, rank })
    }

    /// Build the channel of a cell from its prelim and the lower-rank siblings already
    /// built at this level (empty for rank 0).
    fn build_channel(&self, k: usize, pre: &Prelim, siblings: Option<&Binned>) -> Option<Channel> {
        let lk = self.spacing[k];
        let q = pre.q;
        let mut cand = pre.parent;
        let mut coarse = pre.coarse;
        let mut sibling = false;
        // Siblings: lower-rank, short (≤ SIBLING_MAX_LEN ℓ) channels of this level within
        // SIBLING_SEARCH ℓ; the lowest wedge at q wins, coarser or sibling.
        if let Some(sb) = siblings {
            sb.for_each_near(q, SIBLING_SEARCH * lk + SIBLING_MAX_LEN * lk, |ch| {
                if ch.chan_len > SIBLING_MAX_LEN * lk {
                    return;
                }
                let (d, t) = ch.nearest(q, self.radius_m);
                let wq = ch.wedge(d, t, 1.0, 0.0, 0.0);
                coarse = coarse.min(wq);
                if d > SIBLING_SEARCH * lk || d < 1.0 || cand.wq <= wq {
                    return;
                }
                cand = Cand { wq, d, t, ch: *ch };
                sibling = true;
            });
        }
        let Cand { d, t, ch: parent, .. } = cand;
        if d < MIN_ATTACH * lk {
            return None; // a stub: too close to its parent to make a valley of its own
        }
        // Junction: shifted downstream along the parent so the tributary meets it at α,
        // drawn per cell (Seybold et al. 2017: humid basins ~60–70°, wide scatter).
        let alpha = (JUNCTION_ANGLE_MIN_DEG + (JUNCTION_ANGLE_MAX_DEG - JUNCTION_ANGLE_MIN_DEG) * pre.u[0]).to_radians();
        let parent_len = (parent.up - parent.down).length() * self.radius_m;
        let mut tj = (t - d / alpha.tan() / parent_len.max(1.0)).max(0.0);
        let mut junction = (parent.down + (parent.up - parent.down) * tj).normalize();
        // Length bounds the gather margins assume (see `margin`).
        let max_len = if sibling { SIBLING_MAX_LEN * lk } else { 1.5 * ATTACH_SEARCH * self.spacing[k - 1] };
        if (q - junction).length() * self.radius_m > max_len {
            tj = t;
            junction = (parent.down + (parent.up - parent.down) * tj).normalize();
            if (q - junction).length() * self.radius_m > max_len {
                return None;
            }
        }
        let z_a = parent.z_at(tj);
        let len = (q - junction).length() * self.radius_m;
        // LOD: a sibling's child shares the sibling's spacing (never shows while its parent
        // is faded); a coarser parent's child draws its own.
        let lod_spacing = if sibling { parent.lod_spacing } else { lk * (0.75 + 0.58 * pre.u[1]) };
        // Fading ancestors: the parent's, plus the parent (unless level 0, which never fades)
        // at the junction's distance along it. Same-spacing entries merge by summing (exact:
        // equal presence), which keeps sibling chains within MAX_ANCESTORS.
        let (mut anc, mut n_anc) = (parent.anc, parent.n_anc);
        if parent.level > 0 {
            let entry = [parent.lod_spacing as f32, (parent.s0 + tj * parent.len_m) as f32];
            if let Some(e) = anc[..n_anc as usize].iter_mut().find(|e| e[0] == entry[0]) {
                e[1] += entry[1];
            } else if (n_anc as usize) < MAX_ANCESTORS {
                anc[n_anc as usize] = entry;
                n_anc += 1;
            } else {
                return None; // can't represent the fade exactly: skip (always safe)
            }
        }
        // Stream-power profile from the junction: S = k_sn·A^−θ with A ≈ one level-k cell.
        let ksn = parent.ksn;
        let s_wall = wall_slope(ksn).max(0.002);
        let w = 0.25 * lk;
        let s_channel = (ksn * (lk * lk).powf(-THETA)).min(S_CHANNEL_MAX);
        let hollow = 0.5 * s_wall * w;
        let head_target = pre.t_q + tune("RG_HEAD", HEAD_SHARE) * self.ridge(pre.t_q);
        // Incised into the hillslope it drains: the head sits a hollow below the coarser
        // surface (and below the sibling hillsides seen in the search).
        let z_head = (z_a + s_channel * len).min(head_target - hollow).min(coarse - hollow).max(z_a + EPS_SLOPE * len);
        Some(Channel {
            level: k as u8,
            up: q,
            down: junction,
            z_up: z_head,
            z_down: z_a,
            ksn,
            s_wall,
            s_top: top_share(ksn) * s_wall,
            w,
            floor: floor_width(lk * lk * 1.0e-6) * (1.0 - 0.8 * threshold_fraction(ksn)),
            reach: 2.0 * lk,
            wall_k: self.wall_for(q, len, 2.0 * lk, lk, z_a).0,
            reach_end: self.wall_for(q, len, 2.0 * lk, lk, z_a).1,
            lod_spacing,
            len_m: len,
            s0: 0.0,
            chan_len: len,
            anc,
            n_anc,
        })
    }

    /// Split a straight channel into `CHILD_PIECES` pieces along a gently meandering path:
    /// hashed lateral offsets at the inner vertices (≤ 0.08 L next to the junction, so the
    /// junction angle moves < ~15°; ≤ 0.15 L further up). z is linear by arc length (strictly
    /// rising); every piece shares the channel's LOD and ancestors, and piece i+1's junction
    /// is piece i's head, so the drainage argument holds piece by piece.
    fn curve_pieces(&self, k: usize, c: Cell, ch: Channel) -> Vec<Channel> {
        let (a, q) = (ch.down, ch.up);
        let chord = (q - a).length();
        if chord <= 0.0 {
            return vec![ch];
        }
        let lat = (q - a).cross(a + q).normalize_or_zero();
        let h = lattice::hash(self.salt ^ 0xC0_7E00_0000_0000 ^ k as u64, c);
        let u = |i: u32| ((h >> (i * 21)) & 0x1F_FFFF) as f64 / 0x1F_FFFF as f64 * 2.0 - 1.0;
        let offs = [0.0, 0.08 * u(0), 0.15 * u(1), 0.0];
        let pts: Vec<DVec3> = (0..=CHILD_PIECES)
            .map(|i| {
                let f = i as f64 / CHILD_PIECES as f64;
                (a + (q - a) * f + lat * (offs[i] * chord)).normalize()
            })
            .collect();
        let lens: Vec<f64> = pts.windows(2).map(|w| (w[1] - w[0]).length() * self.radius_m).collect();
        let total: f64 = lens.iter().sum::<f64>().max(1e-9);
        // Concave profile: along the channel the drainage area shrinks towards the head, and
        // the stream-power steady state gives S ∝ A^−θ, so the bed steepens upstream. With
        // Hack's law the area at distance u from the head is A ∝ u^(1/h); with a channel-head
        // area (hollow) of ~5% of the junction's: A(s) = A_j·(0.05 + 0.95·((L − s)/L)^(1/h)).
        // z(s) = z_a + Δz·∫₀ˢ A^−θ / ∫₀ᴸ A^−θ, piece by piece (midpoint rule, 8 sub-steps):
        // strictly rising, the same total rise.
        let a_rel = |s: f64| 0.05 + 0.95 * ((total - s).max(0.0) / total).powf(1.0 / HACK_H);
        let steps = 8;
        let piece_int: Vec<f64> = {
            let mut acc = 0.0;
            lens.iter()
                .map(|&l| {
                    let s_start = acc;
                    acc += l;
                    (0..steps).map(|k| a_rel(s_start + (k as f64 + 0.5) / steps as f64 * l).powf(-THETA) * l / steps as f64).sum::<f64>()
                })
                .collect()
        };
        let int_total: f64 = piece_int.iter().sum::<f64>().max(1e-12);
        let mut out = Vec::with_capacity(CHILD_PIECES);
        let mut s_acc = 0.0;
        let mut i_acc = 0.0;
        for (i, w) in pts.windows(2).enumerate() {
            let z0 = ch.z_down + (ch.z_up - ch.z_down) * (i_acc / int_total);
            i_acc += piece_int[i];
            let z1 = ch.z_down + (ch.z_up - ch.z_down) * (i_acc / int_total);
            out.push(Channel { up: w[1], down: w[0], z_up: z1, z_down: z0, len_m: lens[i], s0: s_acc, chan_len: total, ..ch });
            s_acc += lens[i];
        }
        out
    }

    /// Gather every channel that can matter inside the disc of `radius_m` around `center`,
    /// for levels present at band-limit `cutoff_m`.
    pub fn patch(&self, center: DVec3, radius_m: f64, cutoff_m: f64) -> Patch {
        let levels = (0..self.spacing.len()).take_while(|&k| self.presence(k, cutoff_m) > 0.0).count().max(1);
        let mut binned: Vec<Binned> = Vec::with_capacity(levels);
        // Level 0: macro nodes in a window around the centre.
        let g0 = radius_m + self.margin(0);
        let mut segs = Vec::new();
        let mut l0 = Vec::new();
        self.macro_net.segments_in_disc(center, g0, self.radius_m, &mut segs);
        for s in &segs {
            l0.push(self.level0_channel(s));
        }
        binned.push(Binned::new(l0, self.cells[0], self.radius_m));
        #[cfg(not(target_arch = "wasm32"))]
        let prof = std::env::var("RG_NET_PROF").is_ok();
        for k in 1..levels {
            #[cfg(not(target_arch = "wasm32"))]
            let t_level = std::time::Instant::now();
            let (mut n_cells, mut n_hits, mut n_pre) = (0usize, 0usize, 0usize);
            let lk = self.spacing[k];
            let g = radius_m + self.margin(k);
            let n = self.cells[k];
            // Rank r cells are needed out to g + (2 − r)·RANK_STEP (their children's searches).
            let step = (SIBLING_SEARCH + SIBLING_MAX_LEN) * lk;
            let g_max = g + 2.0 * step;
            let rcells = (g_max / lk).ceil() as i32 + 1;
            let mut by_rank: [Vec<Channel>; 3] = Default::default();
            let mut pending: Vec<(Cell, Prelim, f64)> = Vec::new();
            for c in lattice::window(center, n, rcells) {
                let dist = (lattice::centre(c, n) - center).length() * self.radius_m;
                if dist > g_max + lk {
                    continue;
                }
                n_cells += 1;
                let shard = &self.memo[(lattice::hash(k as u64, c) % 16) as usize];
                if let Some((r, pieces)) = shard.lock().ok().and_then(|m| m.get(&(k as u8, c)).cloned()) {
                    by_rank[r as usize].extend(pieces);
                    n_hits += 1;
                    continue;
                }
                n_pre += 1;
                match self.prelim(k, c, &binned) {
                    None => self.memo_put(k, c, (0, Vec::new())),
                    Some(pre) if pre.rank == 0 => {
                        let pieces = self.build_channel(k, &pre, None).map_or_else(Vec::new, |ch| self.curve_pieces(k, c, ch));
                        self.memo_put(k, c, (0, pieces.clone()));
                        by_rank[0].extend(pieces);
                    }
                    Some(pre) => pending.push((c, pre, dist)),
                }
            }
            for r in 1..3u8 {
                let lower: Vec<Channel> = by_rank[..r as usize].iter().flatten().copied().collect();
                let sib = Binned::new(lower, n, self.radius_m);
                let lim = g + (2 - r) as f64 * step + lk;
                for (c, pre, dist) in pending.iter().filter(|p| p.1.rank == r) {
                    if *dist > lim {
                        continue;
                    }
                    let pieces = self.build_channel(k, pre, Some(&sib)).map_or_else(Vec::new, |ch| self.curve_pieces(k, *c, ch));
                    self.memo_put(k, *c, (r, pieces.clone()));
                    by_rank[r as usize].extend(pieces);
                }
            }
            let all: Vec<Channel> = by_rank.into_iter().flatten().collect();
            #[cfg(not(target_arch = "wasm32"))]
            if prof {
                let (lists, total) = (CAND_STATS[0].swap(0, std::sync::atomic::Ordering::Relaxed), CAND_STATS[1].swap(0, std::sync::atomic::Ordering::Relaxed));
                eprintln!("  L{k}: {n_cells} cells, {n_hits} memo hits, {n_pre} prelims, {} channels, {:.1} ms; {lists} lists, mean len {:.0}", all.len(), t_level.elapsed().as_secs_f64() * 1e3, total as f64 / lists.max(1) as f64);
            }
            let _ = (n_cells, n_hits, n_pre);
            binned.push(Binned::new(all, n, self.radius_m));
        }
        #[cfg(not(target_arch = "wasm32"))]
        let t_blocks = std::time::Instant::now();
        let p = Patch::new(self, center, radius_m, cutoff_m, binned);
        #[cfg(not(target_arch = "wasm32"))]
        if prof {
            eprintln!("  blocks+prune {:.1} ms", t_blocks.elapsed().as_secs_f64() * 1e3);
        }
        p
    }
}

/// Channels of one level, binned by the lattice cell of their head.
struct Binned {
    chans: Vec<Channel>,
    bins: HashMap<Cell, Vec<u32>>,
    n: u32,
    spacing_m: f64,
    /// Longest channel, m (heads can be this far from any point the channel passes).
    max_len: f64,
}

impl Binned {
    fn new(chans: Vec<Channel>, n: u32, radius_m: f64) -> Self {
        let mut bins: HashMap<Cell, Vec<u32>> = HashMap::new();
        let mut max_len = 0.0f64;
        for (i, ch) in chans.iter().enumerate() {
            bins.entry(lattice::cell_of(ch.up, n)).or_default().push(i as u32);
            max_len = max_len.max((ch.up - ch.down).length() * radius_m);
        }
        let spacing_m = std::f64::consts::FRAC_PI_2 * radius_m / n as f64;
        Self { chans, bins, n, spacing_m, max_len }
    }

    /// Visit every channel whose head is within `dist_m` (+ one cell) of `p`. Near a seam a
    /// channel can be visited twice; callers are idempotent (nearest-point searches).
    fn for_each_near(&self, p: DVec3, dist_m: f64, mut f: impl FnMut(&Channel)) {
        let r = (dist_m / self.spacing_m).ceil() as i32 + 1;
        let c0 = lattice::cell_of(p, self.n);
        for dj in -r..=r {
            for di in -r..=r {
                if let Some(ix) = self.bins.get(&lattice::offset(c0, di, dj, self.n)) {
                    for &i in ix {
                        f(&self.chans[i as usize]);
                    }
                }
            }
        }
    }
}

/// The channels that can reach one region, in spatial blocks for fast evaluation.
pub struct Patch {
    center: DVec3,
    east: DVec3,
    north: DVec3,
    radius_m: f64,
    planet_r: f64,
    block_m: f64,
    blocks_per_side: usize,
    /// Per block: channels that can reach it, sorted by junction elevation.
    blocks: Vec<Vec<Channel>>,
    /// The three band-limits themselves, m (this, parent, child): minimum floor widths.
    cutoffs: [f64; 3],
    net: Arc<CubeField<f32>>,
    ridge_alpha: f64,
    lake_depth: Option<Arc<CubeField<f32>>>,
}

impl Patch {
    fn new(dn: &DetailNet, center: DVec3, radius_m: f64, cutoff_m: f64, binned: Vec<Binned>) -> Self {
        let east = DVec3::Z.cross(center).normalize_or(DVec3::X);
        let north = center.cross(east);
        // Blocks about two cut-offs across (≥ 150 m), at most 16 per side.
        let blocks_per_side = ((2.0 * radius_m / (2.0 * cutoff_m).max(150.0)).ceil() as usize).clamp(1, 16);
        let block_m = 2.0 * radius_m / blocks_per_side as f64;
        let mut blocks = vec![Vec::new(); blocks_per_side * blocks_per_side];
        let half_diag = block_m * std::f64::consts::FRAC_1_SQRT_2;
        let pr = dn.radius_m;
        let to_xy = |d: DVec3| {
            let rel = d - center;
            (rel.dot(east) * pr, rel.dot(north) * pr)
        };
        for b in &binned {
            for ch in &b.chans {
                // A block can feel this channel only within reach + its wall (1.5 × reach)
                // of the segment: test the blocks overlapping the segment's padded box.
                let lim = ch.reach_end + half_diag;
                let ((x0, y0), (x1, y1)) = (to_xy(ch.up), to_xy(ch.down));
                let bi = |v: f64| ((v + radius_m) / block_m).floor();
                let (bx0, bx1) = (bi(x0.min(x1) - lim).max(0.0), bi(x0.max(x1) + lim).min(blocks_per_side as f64 - 1.0));
                let (by0, by1) = (bi(y0.min(y1) - lim).max(0.0), bi(y0.max(y1) + lim).min(blocks_per_side as f64 - 1.0));
                if bx0 > bx1 || by0 > by1 {
                    continue;
                }
                for by in by0 as usize..=by1 as usize {
                    for bx in bx0 as usize..=bx1 as usize {
                        let x = -radius_m + (bx as f64 + 0.5) * block_m;
                        let y = -radius_m + (by as f64 + 0.5) * block_m;
                        let p = (center + east * (x / pr) + north * (y / pr)).normalize();
                        if ch.nearest(p, pr).0 <= lim {
                            blocks[by * blocks_per_side + bx].push(*ch);
                        }
                    }
                }
            }
        }
        let cutoffs = [cutoff_m, 2.0 * cutoff_m, 0.5 * cutoff_m];
        // Exact pruning per block. Anywhere in the block a wedge is ≥ its lower bound
        // z_down + H(max(0, d_c − hd)) (z ≥ z_down, H monotone), and the block's minimum is
        // ≤ the best value at the centre + MAX_GRADIENT·hd (no surface is steeper within
        // reach). A channel whose bound exceeds that can never be the minimum there.
        for by in 0..blocks_per_side {
            for bx in 0..blocks_per_side {
                let list = &mut blocks[by * blocks_per_side + bx];
                if list.is_empty() {
                    continue;
                }
                let x = -radius_m + (bx as f64 + 0.5) * block_m;
                let y = -radius_m + (by as f64 + 0.5) * block_m;
                let pc = (center + east * (x / pr) + north * (y / pr)).normalize();
                let t = dn.macro_tri(pc);
                let mut ub = [t + dn.ridge(t); 3];
                for ch in list.iter() {
                    let (d, tt) = ch.nearest(pc, pr);
                    for i in 0..3 {
                        let ph = ch.presence(cutoffs[i]);
                        if ph > 0.0 {
                            ub[i] = ub[i].min(ch.wedge(d, tt, ph, cutoffs[i], ch.inherited_offset(cutoffs[i])));
                        }
                    }
                }
                let ub_max = ub.iter().fold(f64::MIN, |a, &b| a.max(b)) + MAX_GRADIENT * half_diag;
                list.retain(|ch| {
                    let (d, _) = ch.nearest(pc, pr);
                    let dl = (d - half_diag).max(0.0);
                    // The lower bound: full presence (fading only raises a wedge), the junction
                    // elevation, and the widest floor of the three LODs (wider floor → lower H).
                    ch.wedge(dl, 0.0, 1.0, cutoffs[1], 0.0) <= ub_max
                });
                list.sort_by(|a, c| a.z_min().total_cmp(&c.z_min()));
            }
        }
        Self { center, east, north, radius_m, planet_r: pr, block_m, blocks_per_side, blocks, cutoffs, net: dn.macro_net.height.clone(), ridge_alpha: dn.ridge_alpha, lake_depth: dn.macro_net.lake_depth_m.clone() }
    }

    /// Heights at this patch's band-limit, its parent LOD (2× coarser) and child LOD
    /// (2× finer, as far as gathered levels allow). Ocean returns the macro (bathymetry).
    pub fn heights(&self, p: DVec3) -> [f64; 3] {
        let t = self.net.sample_tri(p);
        let u = t + self.ridge_alpha * t.clamp(0.0, RIDGE_FULL_M);
        if t <= 0.0 {
            let m = self.net.sample_drained(p);
            return [m; 3];
        }
        let rel = p - self.center;
        let (x, y) = (rel.dot(self.east) * self.planet_r, rel.dot(self.north) * self.planet_r);
        let bx = (((x + self.radius_m) / self.block_m).floor().max(0.0) as usize).min(self.blocks_per_side - 1);
        let by = (((y + self.radius_m) / self.block_m).floor().max(0.0) as usize).min(self.blocks_per_side - 1);
        let mut best = [u; 3];
        for ch in &self.blocks[by * self.blocks_per_side + bx] {
            if ch.z_min() >= best[0].max(best[1]).max(best[2]) {
                break;
            }
            let (d, t) = ch.nearest(p, self.planet_r);
            for (b, &cut) in best.iter_mut().zip(&self.cutoffs) {
                let pres = ch.presence(cut);
                if pres > 0.0 {
                    *b = b.min(ch.wedge(d, t, pres, cut, ch.inherited_offset(cut)));
                }
            }
        }
        // Land stays land: every wedge is ≥ 0 (z ≥ 0 at sea) and the cap → 0 at the coast.
        let mut out = best.map(|h| h.max(1.0e-3_f64.min(u)));
        // Lakes: the water surface is flat at the lake's level (bed + depth at its nodes).
        // Near lake nodes (the cell's corners), water shows wherever the terrain is below that
        // level, so shorelines follow the terrain's own contours.
        if let Some(ld) = &self.lake_depth {
            // Only inside the lake's macro footprint (bed below the level), so the shoreline
            // follows the bed's contour rather than cell edges.
            if let Some(level) = lake_level_near(&self.net, ld, p).filter(|&l| self.net.sample_tri(p) < l) {
                out = out.map(|h| h.max(level));
            }
        }
        out
    }

    /// Diagnostics: true where the active wedge is inside its reach wall (d > reach), the
    /// source of small crescent scarps.
    pub fn in_wall(&self, p: DVec3) -> bool {
        let t0 = self.net.sample_tri(p);
        let mut best = (t0 + self.ridge_alpha * t0.clamp(0.0, RIDGE_FULL_M), false);
        let rel = p - self.center;
        let (x, y) = (rel.dot(self.east) * self.planet_r, rel.dot(self.north) * self.planet_r);
        let bx = (((x + self.radius_m) / self.block_m).floor().max(0.0) as usize).min(self.blocks_per_side - 1);
        let by = (((y + self.radius_m) / self.block_m).floor().max(0.0) as usize).min(self.blocks_per_side - 1);
        for ch in &self.blocks[by * self.blocks_per_side + bx] {
            let (d, t) = ch.nearest(p, self.planet_r);
            let pres = ch.presence(self.cutoffs[0]);
            if pres > 0.0 {
                let w = ch.wedge(d, t, pres, self.cutoffs[0], ch.inherited_offset(self.cutoffs[0]));
                if w < best.0 {
                    best = (w, d > ch.reach);
                }
            }
        }
        best.1
    }

    /// Diagnostics: the active term at `p` for the patch's own band-limit: `None` if the cap
    /// is active, else (level, t along the channel from its junction, distance m).
    pub fn active(&self, p: DVec3) -> Option<(u8, f64, f64)> {
        let t0 = self.net.sample_tri(p);
        let mut best = (t0 + self.ridge_alpha * t0.clamp(0.0, RIDGE_FULL_M), None);
        let rel = p - self.center;
        let (x, y) = (rel.dot(self.east) * self.planet_r, rel.dot(self.north) * self.planet_r);
        let bx = (((x + self.radius_m) / self.block_m).floor().max(0.0) as usize).min(self.blocks_per_side - 1);
        let by = (((y + self.radius_m) / self.block_m).floor().max(0.0) as usize).min(self.blocks_per_side - 1);
        for ch in &self.blocks[by * self.blocks_per_side + bx] {
            let (d, t) = ch.nearest(p, self.planet_r);
            let pres = ch.presence(self.cutoffs[0]);
            if pres > 0.0 {
                let w = ch.wedge(d, t, pres, self.cutoffs[0], ch.inherited_offset(self.cutoffs[0]));
                if w < best.0 {
                    best = (w, Some((ch.level, t, d)));
                }
            }
        }
        best.1
    }

    /// Every channel held by the patch's blocks (a channel can appear in several blocks).
    pub fn channels(&self) -> impl Iterator<Item = &Channel> {
        self.blocks.iter().flatten()
    }

    /// Approximate heap size, bytes (for cache budgets).
    pub fn bytes(&self) -> usize {
        self.channel_count() * std::mem::size_of::<Channel>() + self.blocks.len() * std::mem::size_of::<Vec<Channel>>()
    }

    pub fn channel_count(&self) -> usize {
        self.blocks.iter().map(|b| b.len()).sum()
    }
}

impl MacroNet {
    /// Land segments whose upstream node is within `dist_m` of `center`.
    pub fn segments_in_disc(&self, center: DVec3, dist_m: f64, radius_m: f64, out: &mut Vec<Segment>) {
        out.clear();
        let (face, u, v) = CubeFace::from_dir(center);
        let r = self.res() as f64;
        let cell = std::f64::consts::FRAC_PI_2 * radius_m / r;
        let half = (dist_m / cell).ceil() as i64 + 2;
        let (ci, cj) = (((u + 1.0) * 0.5 * r).round() as i64, ((v + 1.0) * 0.5 * r).round() as i64);
        let mut seen = std::collections::HashSet::new();
        for j in cj - half..=cj + half {
            for i in ci - half..=ci + half {
                let n = self.node(face, i, j);
                if !seen.insert(n) {
                    continue;
                }
                let d = self.dir(n);
                if (d - center).length() * radius_m > dist_m {
                    continue;
                }
                if let Some(m) = self.receiver(n) {
                    self.curved_edge(n, m, out);
                }
            }
        }
    }

    /// The node draining into `n` with the largest area (the main stem upstream), if any.
    fn main_donor(&self, n: Node) -> Option<Node> {
        let mut best: Option<(f32, Node)> = None;
        for dj in -1..=1i64 {
            for di in -1..=1i64 {
                let m = self.node(n.face, n.i as i64 + di, n.j as i64 + dj);
                if m == n || self.receiver(m) != Some(n) {
                    continue;
                }
                let a = self.area_km2.data[self.height.index(m.face, m.i, m.j)];
                if best.is_none_or(|(ba, bm)| a > ba || (a == ba && m < bm)) {
                    best = Some((a, m));
                }
            }
        }
        best.map(|b| b.1)
    }

    /// The edge n → m as `CURVE_PIECES` segments along a centripetal Catmull-Rom curve
    /// through (main donor of n, n, m, receiver of m). The curve passes through the nodes, so
    /// junctions and node heights stay exact; z is linear by arc length (strictly monotone).
    /// Bulge ≤ 0.15 × the edge; an edge stays straight where the curve would dip under the
    /// macro cap (the river would vanish there).
    fn curved_edge(&self, n: Node, m: Node, out: &mut Vec<Segment>) {
        let (p1, p2) = (self.dir(n), self.dir(m));
        let p0 = self.main_donor(n).map_or(p1 * 2.0 - p2, |d| self.dir(d));
        let p3 = self.receiver(m).map_or(p2 * 2.0 - p1, |r| self.dir(r));
        let z_up = self.surface_at(n);
        let z_down = self.surface_at(m).max(0.0);
        let area_km2 = self.area_km2.data[self.height.index(n.face, n.i, n.j)] as f64;
        let uplift = self.uplift_at(n);
        // Centripetal knots (α = 0.5): no cusps or self-intersections.
        let knot = |a: DVec3, b: DVec3| (a - b).length().sqrt().max(1e-12);
        let (t0, t1) = (0.0, knot(p0, p1));
        let t2 = t1 + knot(p1, p2);
        let t3 = t2 + knot(p2, p3);
        let lerp = |a: DVec3, b: DVec3, ta: f64, tb: f64, t: f64| a * ((tb - t) / (tb - ta)) + b * ((t - ta) / (tb - ta));
        let at = |t: f64| {
            let a1 = lerp(p0, p1, t0, t1, t);
            let a2 = lerp(p1, p2, t1, t2, t);
            let a3 = lerp(p2, p3, t2, t3, t);
            let b1 = lerp(a1, a2, t0, t2, t);
            let b2 = lerp(a2, a3, t1, t3, t);
            lerp(b1, b2, t1, t2, t).normalize()
        };
        let chord = (p2 - p1).length();
        let mut pts: Vec<DVec3> = (0..=CURVE_PIECES).map(|i| at(t1 + (t2 - t1) * i as f64 / CURVE_PIECES as f64)).collect();
        pts[0] = p1;
        pts[CURVE_PIECES] = p2;
        // Bound the bulge, and keep the bed under the cap.
        let mut straight = false;
        for (i, q) in pts.iter_mut().enumerate().take(CURVE_PIECES).skip(1) {
            let f = i as f64 / CURVE_PIECES as f64;
            let on_chord = (p1 + (p2 - p1) * f).normalize();
            let dev = (*q - on_chord).length();
            if dev > 0.15 * chord {
                *q = (on_chord + (*q - on_chord) * (0.15 * chord / dev)).normalize();
            }
            let z = z_up + (z_down - z_up) * f;
            let tq = self.height.sample_tri(*q);
            // The bed must stay below the macro surface there (else the river vanishes under
            // the cap); allow it to graze it within a few metres.
            if tq + 5.0 < z {
                straight = true;
            }
        }
        if straight {
            for (i, q) in pts.iter_mut().enumerate() {
                *q = (p1 + (p2 - p1) * (i as f64 / CURVE_PIECES as f64)).normalize();
            }
        }
        // z by arc length.
        let lens: Vec<f64> = pts.windows(2).map(|w| (w[1] - w[0]).length()).collect();
        let total: f64 = lens.iter().sum::<f64>().max(1e-15);
        let mut s_acc = 0.0;
        for (i, w) in pts.windows(2).enumerate() {
            let za = z_up + (z_down - z_up) * (s_acc / total);
            s_acc += lens[i];
            let zb = z_up + (z_down - z_up) * (s_acc / total);
            out.push(Segment { up: w[0], down: w[1], z_up: za, z_down: zb, area_km2, uplift_mm_yr: uplift });
        }
    }
}

/// Native-only tuning override from an environment variable (sweeps in the container).
fn tune(var: &str, default: f64) -> f64 {
    #[cfg(not(target_arch = "wasm32"))]
    if let Some(v) = std::env::var(var).ok().and_then(|v| v.parse().ok()) {
        return v;
    }
    let _ = var;
    default
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn floor_hinge_is_flat_then_c1() {
        let f = 100.0;
        assert_eq!(floor_hinge(0.0, f), 0.0);
        assert_eq!(floor_hinge(50.0, f), 0.0);
        let e = 1e-6;
        for &x in &[50.0, 150.0] {
            let left = (floor_hinge(x, f) - floor_hinge(x - e, f)) / e;
            let right = (floor_hinge(x + e, f) - floor_hinge(x, f)) / e;
            assert!((left - right).abs() < 1e-3, "kink at {x}: {left} vs {right}");
        }
        assert!((floor_hinge(400.0, f) - 300.0).abs() < 1e-9);
        // Non-decreasing everywhere, strictly increasing past the floor.
        let mut prev = -1.0;
        for k in 0..1000 {
            let v = floor_hinge(k as f64, f);
            assert!(v >= prev);
            if k > 50 {
                assert!(v > prev);
            }
            prev = v;
        }
    }
}
