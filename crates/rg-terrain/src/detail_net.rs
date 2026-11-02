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

use crate::network::{MacroNet, Segment};
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
/// A junction moves downstream along its parent by this share of the attach distance, so
/// tributaries meet at ~60° instead of 90°.
const JUNCTION_SHIFT: f64 = 0.5;
/// Floodplain width per √(km²) of drainage (see `floor_width`), m.
const FLOOR_M_PER_SQRT_KM2: f64 = 30.0;
/// Macro height above which the ridge allowance stops growing, m.
const RIDGE_FULL_M: f64 = 3000.0;
/// Ridge allowance per metre of macro height (see `DetailNet::ridge`).
const RIDGE_ALPHA: f64 = 0.3;
/// Junction search radius, in parent spacings. Tributaries meet the nearest coarser channel
/// well within one parent spacing; a head with nothing in range simply has no channel.
/// Sets the gather margins (see `DetailNet::margin`).
const ATTACH_SEARCH: f64 = 1.2;
/// Fade-out terms (see `wedge`): large enough that a faded channel never wins.
const FADE_RAMP_M: f64 = 1.0e4;
const FADE_SLOPE: f64 = 50.0;

/// One shard of the channel memo: (level, cell) → channel (None = no channel there).
type MemoShard = std::sync::Mutex<HashMap<(u8, Cell), Option<Channel>>>;

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

    /// The wedge at distance `d`, parameter `t`, faded in by `phi` (1 = fully present), with
    /// the valley floor at least `min_floor` wide (band-limiting: a level of detail can't
    /// show a floor narrower than it resolves, and sampling a sharper V makes grid pits).
    /// Every term is strictly increasing in `d` and in `t`, at any `phi`; `min_floor` is
    /// fixed per evaluation, so each wedge's shape is fixed as the proof requires.
    #[inline]
    fn wedge(&self, d: f64, t: f64, phi: f64, min_floor: f64) -> f64 {
        let f = self.floor.max(min_floor);
        let g = floor_hinge(d, f);
        let mut h = self.z_at(t) + (self.s_wall - self.s_top) * self.w * (g / self.w).tanh() + self.s_top * g;
        let over = d - self.reach;
        if over > 0.0 {
            h += self.wall_k * over * over;
        }
        if phi < 1.0 {
            h += (1.0 - phi) * (FADE_RAMP_M * t + FADE_SLOPE * g);
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
    ridge_alpha: f64,
    /// Erodibility K of the macro solve, m^(1−2m)/yr (k_sn = U/K).
    k_erodibility: f64,
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
        let ridge_alpha = tune("RG_RIDGE", RIDGE_ALPHA);
        Self { macro_net, radius_m, salt: seed ^ 0xD371_0A11_C4A7_0001, spacing, cells, memo, ridge_alpha, k_erodibility: 1.0e-5 }
    }

    /// Set the macro solve's erodibility K (default 1e-5), used for k_sn = U/K.
    pub fn with_erodibility(mut self, k: f64) -> Self {
        self.k_erodibility = k;
        self
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

    /// The cap U = T + R(T).
    pub fn cap(&self, p: DVec3) -> f64 {
        let t = self.macro_tri(p);
        t + self.ridge(t)
    }

    /// How present level `k` is at band-limit `cutoff_m` (0 = absent, 1 = full): a level of
    /// spacing ℓ needs about ℓ/4 resolution (its convexity width).
    pub fn presence(&self, k: usize, cutoff_m: f64) -> f64 {
        smoothstep(2.0 * cutoff_m, 4.0 * cutoff_m, self.spacing[k])
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
            wall_k: 1.0e4 / (l0 * l0),
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
        // Also: a head reads the coarser surface (every coarser wedge within reach of it):
        // heads of level k−1 within (3 + 4.5σ) ℓ_{k−1} of it, so c ≥ (3 + 4.5σ)·3/2 ≈ 12.6.
        let c = (8.25 * ATTACH_SEARCH).max((3.0 + 4.5 * ATTACH_SEARCH) * 1.5) + 0.5;
        if k == 0 { (c / 3.0 + ATTACH_SEARCH + 1.5 + 0.5) * self.spacing[0] } else { c * self.spacing[k] }
    }

    /// Build the channel of lattice cell `c` at level `k ≥ 1`, given all coarser channels
    /// near it (binned per level). `None` at sea or if no coarser channel is in range.
    fn channel_of(&self, k: usize, c: Cell, coarser: &[Binned]) -> Option<Channel> {
        let lk = self.spacing[k];
        let q = lattice::jittered(c, self.cells[k], self.salt.wrapping_add(k as u64), JITTER);
        let t_q = self.macro_tri(q);
        if t_q <= 0.0 {
            return None;
        }
        let u_q = t_q + self.ridge(t_q);
        // Channel heads aim at the macro level plus half the ridge allowance, so valleys cut
        // below the macro surface about as much as ridges rise above it.
        let head_target = t_q + 0.5 * self.ridge(t_q);
        // The parent is the coarser channel whose hillslope q sits on: the lowest coarser
        // wedge at q among channels within the search radius. (The *nearest* channel can be
        // buried under a lower neighbour's wedge, and a child hung from it stays buried.)
        let search = ATTACH_SEARCH * self.spacing[k - 1];
        let mut best: Option<(f64, f64, f64, Channel)> = None; // (wedge at q, d, t, channel)
        for b in coarser {
            b.for_each_near(q, search + b.max_len, |ch| {
                let (d, t) = ch.nearest(q, self.radius_m);
                if d > search || d < 1.0 {
                    return;
                }
                let wq = ch.wedge(d, t, 1.0, 0.0);
                if best.as_ref().is_some_and(|(bw, ..)| *bw <= wq) {
                    return;
                }
                best = Some((wq, d, t, *ch));
            });
        }
        let (_, d, t, parent) = best?;
        // Junction: shifted downstream along the parent (towards t = 0).
        let parent_len = (parent.up - parent.down).length() * self.radius_m;
        let tj = (t - JUNCTION_SHIFT * d / parent_len.max(1.0)).max(0.0);
        let junction = (parent.down + (parent.up - parent.down) * tj).normalize();
        let z_a = parent.z_at(tj);
        let len = (q - junction).length() * self.radius_m;
        // Stream-power profile from the junction: S = k_sn·A^−θ with A ≈ one level-k cell.
        let ksn = parent.ksn;
        let s_wall = wall_slope(ksn).max(0.002);
        let w = 0.25 * lk;
        let s_channel = (ksn * (lk * lk).powf(-THETA)).min(S_CHANNEL_MAX);
        let hollow = 0.5 * s_wall * w;
        // A tributary is incised into the hillslope it drains: its head sits a hollow depth
        // below the coarser surface there (all coarser wedges at full presence, intrinsic
        // floors, so the value depends on this cell alone). Without this, steep stream-power
        // profiles climbed above the parent's hillside and never cut it (smooth mesas).
        let mut coarse = u_q;
        for b in coarser {
            b.for_each_near(q, 3.0 * b.spacing_m + b.max_len, |ch| {
                let (d, t) = ch.nearest(q, self.radius_m);
                coarse = coarse.min(ch.wedge(d, t, 1.0, 0.0));
            });
        }
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
            wall_k: 1.0e4 / (lk * lk),
        })
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
        for k in 1..levels {
            let g = radius_m + self.margin(k);
            let n = self.cells[k];
            let rcells = (g / self.spacing[k]).ceil() as i32 + 1;
            let mut out = Vec::new();
            for c in lattice::window(center, n, rcells) {
                // Cheap pre-check on the cell centre (the jittered point is within a cell).
                if (lattice::centre(c, n) - center).length() * self.radius_m > g + self.spacing[k] {
                    continue;
                }
                let shard = &self.memo[(lattice::hash(k as u64, c) % 16) as usize];
                let cached = shard.lock().ok().and_then(|m| m.get(&(k as u8, c)).copied());
                let ch = match cached {
                    Some(ch) => ch,
                    None => {
                        let ch = self.channel_of(k, c, &binned);
                        if let Ok(mut m) = shard.lock() {
                            // Bound memory: drop the shard when it grows large (it's only a cache).
                            // ~50k × 16 shards × ~170 B ≈ 140 MB at most.
                            if m.len() > 50_000 {
                                m.clear();
                            }
                            m.insert((k as u8, c), ch);
                        }
                        ch
                    }
                };
                if let Some(ch) = ch {
                    out.push(ch);
                }
            }
            binned.push(Binned::new(out, n, self.radius_m));
        }
        Patch::new(self, center, radius_m, cutoff_m, binned)
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
    /// Presence of each level at the patch's band-limit and the two neighbouring LODs.
    presence: Vec<[f64; 3]>,
    /// The three band-limits themselves, m (this, parent, child): minimum floor widths.
    cutoffs: [f64; 3],
    net: Arc<CubeField<f32>>,
    ridge_alpha: f64,
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
                let lim = ch.reach * 1.5 + half_diag;
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
        let presence: Vec<[f64; 3]> = (0..binned.len())
            .map(|k| [dn.presence(k, cutoff_m), dn.presence(k, 2.0 * cutoff_m), dn.presence(k, 0.5 * cutoff_m)])
            .collect();
        let cutoffs = [cutoff_m, 2.0 * cutoff_m, 0.5 * cutoff_m];
        // Exact pruning per block. Anywhere in the block a wedge is ≥ its lower bound
        // z_down + H(max(0, d_c − hd)) (z ≥ z_down, H monotone), and the block's minimum is
        // ≤ the best value at the centre + MAX_GRADIENT·hd (no surface is steeper within
        // reach). A channel whose bound exceeds that can never be the minimum there.
        const MAX_GRADIENT: f64 = 1.5;
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
                        let ph = presence[ch.level as usize][i];
                        if ph > 0.0 {
                            ub[i] = ub[i].min(ch.wedge(d, tt, ph, cutoffs[i]));
                        }
                    }
                }
                let ub_max = ub.iter().fold(f64::MIN, |a, &b| a.max(b)) + MAX_GRADIENT * half_diag;
                list.retain(|ch| {
                    let (d, _) = ch.nearest(pc, pr);
                    let dl = (d - half_diag).max(0.0);
                    // The lower bound: full presence (fading only raises a wedge), the junction
                    // elevation, and the widest floor of the three LODs (wider floor → lower H).
                    ch.wedge(dl, 0.0, 1.0, cutoffs[1]) <= ub_max
                });
                list.sort_by(|a, c| a.z_min().total_cmp(&c.z_min()));
            }
        }
        Self { center, east, north, radius_m, planet_r: pr, block_m, blocks_per_side, blocks, presence, cutoffs, net: dn.macro_net.height.clone(), ridge_alpha: dn.ridge_alpha }
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
            let pres = self.presence[ch.level as usize];
            for i in 0..3 {
                if pres[i] > 0.0 {
                    best[i] = best[i].min(ch.wedge(d, t, pres[i], self.cutoffs[i]));
                }
            }
        }
        // Land stays land: every wedge is ≥ 0 (z ≥ 0 at sea) and the cap → 0 at the coast.
        best.map(|h| h.max(1.0e-3_f64.min(u)))
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
            let pres = self.presence[ch.level as usize][0];
            if pres > 0.0 {
                let w = ch.wedge(d, t, pres, self.cutoffs[0]);
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
                    out.push(Segment {
                        up: d,
                        down: self.dir(m),
                        z_up: self.height_at(n),
                        z_down: self.height_at(m).max(0.0),
                        area_km2: self.area_km2.data[self.height.index(n.face, n.i, n.j)] as f64,
                        uplift_mm_yr: self.uplift_at(n),
                    });
                }
            }
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
