//! Terrain at any zoom: the committed macro height field plus procedural detail.
//!
//! See `docs/design/tile-streaming.md`. In short:
//!
//! `h(dir) = macro(dir) + detail(dir)`, where detail is fBm whose octaves start just above
//! the macro field's resolution and continue down to about 10 m. A tile only evaluates
//! the octaves its own sample spacing can represent (band-limiting), so a coarse tile is
//! a smoothed version of the same surface a fine tile shows, not a different surface.
//!
//! **Coast rule:** detail can never move the coastline. Near sea level its amplitude
//! fades out, and the result keeps the sign of the macro height. The coastline you see
//! at 10 m/px is exactly the committed one, so maps and vector layers stay honest.
//!
//! Precision: this CPU version evaluates noise in f64 at `dir · R / λ`, which is exact
//! enough down to centimetres. The GPU port must switch to integer lattice coordinates
//! from the tile address (ARCHITECTURE §4).

pub mod detail_net;
pub mod erosion_filter;
pub mod gradient;
pub mod network;
pub mod tile;

pub use tile::{tile_indices, TileData, TileGenerator, TILE_CELLS};

use glam::DVec3;
use rg_geo::CubeField;
use std::collections::HashMap;
use std::sync::Arc;

/// Version of the detail function. Bump on any change that alters heights; stored sculpt
/// deltas (M6) record the version they were made against.
pub const DETAIL_GENERATOR_VERSION: u32 = 0; // 0 = not persistent yet (no sculpt deltas exist)

/// Network patches by (block cell, cells per face, band-limit bits).
type PatchCache = std::sync::Mutex<HashMap<(rg_geo::lattice::Cell, u32, u64), Arc<detail_net::Patch>>>;

/// Samples the full-resolution terrain surface. Cheap to clone (the field is shared), and
/// `Send`, so native builds can generate tiles on worker threads.
#[derive(Clone)]
pub struct TerrainSampler {
    pub macro_h: Arc<CubeField<f32>>,
    /// Macro drainage area (km²) from the erosion stage, if any: rivers get corridors.
    pub drainage_km2: Option<Arc<CubeField<f32>>>,
    pub radius_m: f64,
    salt: u32,
    ridge_salt: u32,
    gully_salt: u32,
    /// Wavelength of the first detail octave (m): about two macro cells.
    base_wavelength_m: f64,
    /// Wavelength of the first erosion-filter octave (m).
    pub gully_wavelength_m: f64,
    /// Metres per erosion-filter unit (maps the filter's first octave to the above).
    gully_unit_m: f64,
    /// Amplitude tuning (see the constants of the same name).
    hill_slope: f64,
    gully_strength: f64,
    slope_gain: f64,
    /// Finest wavelength ever generated (m).
    floor_wavelength_m: f64,
    /// The drained detail network, when the world has hydrology at macro resolution. It
    /// replaces hills, gullies and river corridors (`docs/design/drained-detail-network.md`).
    pub net: Option<detail_net::DetailNet>,
    /// Patches of the network, shared by clones (tiles on worker threads, probes). Keyed by
    /// (block cell, cells per face, band-limit bits). Only a cache.
    patches: Arc<PatchCache>,
    /// Macro only: no streamed detail of any kind (network, hills, gullies, river
    /// corridors). Every query returns the drained macro height. The app's "streamed
    /// detail: off" (s02); tests and probes keep the default (false).
    pub macro_only: bool,
}



impl TerrainSampler {
    pub fn new(macro_h: Arc<CubeField<f32>>, radius_m: f64, seed: u64) -> Self {
        let macro_cell_m = macro_h.spacing_rad() * radius_m;
        Self {
            macro_h,
            drainage_km2: None,
            radius_m,
            salt: rg_noise::salt(seed, 0x7E44),
            ridge_salt: rg_noise::salt(seed, 0x7E45),
            gully_salt: rg_noise::salt(seed, 0x7E46),
            // Hills start at one macro cell: longer wavelengths are the macro surface's job
            // (and there they are river-eroded; noise there would ignore the drainage).
            base_wavelength_m: macro_cell_m,
            gully_wavelength_m: tune("RG_GULLY_WL", GULLY_WAVELENGTH_M),
            gully_unit_m: tune("RG_GULLY_WL", GULLY_WAVELENGTH_M) / (0.15 * 0.7),
            slope_gain: tune("RG_SLOPE_GAIN", SLOPE_GAIN),
            hill_slope: tune("RG_HILL_SLOPE", HILL_SLOPE),
            gully_strength: tune("RG_GULLY_STRENGTH", GULLY_STRENGTH),
            floor_wavelength_m: 12.0,
            net: None,
            patches: Arc::new(std::sync::Mutex::new(HashMap::new())),
            macro_only: false,
        }
    }

    /// Switch off all streamed detail (see [`Self::macro_only`]).
    pub fn with_macro_only(mut self) -> Self {
        self.macro_only = true;
        self
    }

    /// All four levels at the macro height (macro-only mode, and the open ocean).
    fn flat_levels(m: f64) -> HeightLevels {
        HeightLevels { macro_h: m, h: m, h_parent: m, h_child: m }
    }

    /// Use the drained detail network for everything above sea level.
    pub fn with_network(mut self, net: detail_net::DetailNet) -> Self {
        self.net = Some(net);
        self
    }

    /// A network patch covering the disc `radius_m` around `centre` at `min_wavelength_m`
    /// (`None` without a network). Tiles build one for their own footprint.
    pub fn network_patch(&self, centre: DVec3, radius_m: f64, min_wavelength_m: f64) -> Option<detail_net::Patch> {
        if self.macro_only {
            return None;
        }
        let cutoff = min_wavelength_m.max(self.floor_wavelength_m);
        self.net.as_ref().map(|n| n.patch(centre, radius_m, cutoff))
    }

    /// Heights from a patch made by [`Self::network_patch`] (`dir` must lie inside it).
    pub fn levels_in_patch(&self, patch: &detail_net::Patch, dir: DVec3) -> HeightLevels {
        let [h, hp, hc] = patch.heights(dir);
        HeightLevels { macro_h: self.macro_height(dir), h, h_parent: hp, h_child: hc }
    }

    /// Network heights at (this, parent, child) band-limits via a cached patch. A patch
    /// covers a block of about 32 cut-offs (one tile's extent at that level of detail).
    fn net_levels(&self, net: &detail_net::DetailNet, dir: DVec3, min_wavelength_m: f64) -> HeightLevels {
        let cutoff = min_wavelength_m.max(self.floor_wavelength_m);
        let block_m = (32.0 * cutoff).max(2000.0);
        let n = ((std::f64::consts::FRAC_PI_2 * self.radius_m / block_m).round() as u32).max(1);
        let cell = rg_geo::lattice::cell_of(dir, n);
        let key = (cell, n, cutoff.to_bits());
        let cached = self.patches.lock().ok().and_then(|m| m.get(&key).cloned());
        let patch = match cached {
            Some(p) => p,
            None => {
                let centre = rg_geo::lattice::centre(cell, n);
                // Cells are up to ~1.4 blocks across near face corners; cover the whole cell.
                let p = Arc::new(net.patch(centre, block_m, cutoff));
                if let Ok(mut m) = self.patches.lock() {
                    // Budget ~96 MB: a patch is a few KB to a few MB.
                    let bytes: usize = m.values().map(|p| p.bytes()).sum();
                    if bytes > 96 << 20 {
                        m.clear();
                    }
                    m.insert(key, p.clone());
                }
                p
            }
        };
        let [h, hp, hc] = patch.heights(dir);
        HeightLevels { macro_h: self.macro_height(dir), h, h_parent: hp, h_child: hc }
    }

    /// Attach the macro drainage-area field so rivers carve corridors through the detail.
    pub fn with_drainage(mut self, drainage_km2: Arc<CubeField<f32>>) -> Self {
        self.drainage_km2 = Some(drainage_km2);
        self
    }

    /// River corridor at `dir`: (corridor strength 0–1, channel depth m). Along macro
    /// rivers the detail is suppressed into a valley floor and a channel is cut, deeper
    /// for bigger rivers (depth ∝ √area, cf. hydraulic geometry), so detail relief can never
    /// dam a macro river.
    pub fn river_corridor(&self, dir: DVec3) -> (f64, f64) {
        let Some(d) = self.drainage_km2.as_ref() else { return (0.0, 0.0) };
        let a = d.sample_f64(dir).max(0.0);
        if a < RIVER_CORRIDOR_MIN_KM2 {
            return (0.0, 0.0);
        }
        let c = smoothstep(RIVER_CORRIDOR_MIN_KM2.ln(), RIVER_CORRIDOR_FULL_KM2.ln(), a.ln());
        let depth = (RIVER_DEPTH_M_PER_SQRT_KKM2 * (a / 1000.0).sqrt()).min(RIVER_DEPTH_MAX_M);
        (c, depth * c)
    }

    /// Macro (committed) height at a unit direction, metres above sea level. Drained
    /// (bilinear with valleys carved) normally; in macro-only mode a smooth bicubic for
    /// display, with the coastline kept exactly where the drained surface puts it.
    pub fn macro_height(&self, dir: DVec3) -> f64 {
        let drained = self.macro_h.sample_drained(dir);
        if !self.macro_only {
            return drained;
        }
        let smooth = self.macro_h.sample_smooth(dir);
        if (smooth > 0.0) == (drained > 0.0) { smooth } else { drained }
    }

    /// Full height (metres above sea level) at `dir`, keeping only detail octaves with a
    /// wavelength of at least `min_wavelength_m`. Returns (height, macro height).
    pub fn height(&self, dir: DVec3, min_wavelength_m: f64) -> (f64, f64) {
        let s = self.levels(dir, min_wavelength_m);
        (s.h, s.macro_h)
    }

    /// Heights at three neighbouring LOD levels at once (one pass over the octaves):
    /// - `h`: this level (octaves with wavelength ≥ `min_wavelength_m`),
    /// - `h_parent`: the parent level (≥ 2 × that), used to morph tiles in without popping,
    /// - `h_child`: the child level (≥ ½ × that), used to measure this tile's error.
    pub fn levels(&self, dir: DVec3, min_wavelength_m: f64) -> HeightLevels {
        if self.macro_only {
            return Self::flat_levels(self.macro_height(dir));
        }
        if let Some(net) = &self.net {
            return self.net_levels(net, dir, min_wavelength_m);
        }
        let m = self.macro_height(dir);
        if m <= 0.0 {
            return HeightLevels { macro_h: m, h: m, h_parent: m, h_child: m };
        }
        // Single-point queries estimate the slope by forward differences; tiles pass the
        // slope from their own grid instead (see `tile.rs`), which is much cheaper.
        let east = DVec3::Z.cross(dir).normalize_or(DVec3::X);
        let north = dir.cross(east);
        let step = tune("RG_STEER_STEP", STEER_STEP_M);
        let ang = step / self.radius_m;
        let b = |v: DVec3| self.base_surface((dir + v * ang).normalize());
        let grad3 = east * ((b(east) - b(-east)) / (2.0 * step)) + north * ((b(north) - b(-north)) / (2.0 * step));
        self.levels_with_slope(dir, min_wavelength_m, grad3)
    }

    /// As [`Self::levels`], with the slope (3D tangent gradient, m/m) of the base surface
    /// supplied by the caller.
    pub fn levels_with_slope(&self, dir: DVec3, min_wavelength_m: f64, grad3: DVec3) -> HeightLevels {
        if self.macro_only {
            return Self::flat_levels(self.macro_height(dir));
        }
        if let Some(net) = &self.net {
            return self.net_levels(net, dir, min_wavelength_m);
        }
        let m = self.macro_height(dir);
        // Oceans render flat at sea level, so skip detail there entirely (~70% of the
        // planet). The coast fade below makes this continuous at m = 0.
        if m <= 0.0 {
            return HeightLevels { macro_h: m, h: m, h_parent: m, h_child: m };
        }
        let lim = min_wavelength_m.max(self.floor_wavelength_m);
        let lim_parent = (2.0 * min_wavelength_m).max(self.floor_wavelength_m);
        let lim_child = (0.5 * min_wavelength_m).max(self.floor_wavelength_m);

        // 1. Hills: low-frequency relief between the macro cell size and the largest gully
        //    wavelength, band-limited per LOD.
        let hl = self.hills(dir, m, [lim, lim_parent, lim_child]);

        // 2. Gullies: the erosion filter carves branching drainage down the slope of the
        //    *full* hill surface (not band-limited, so gully directions don't shift between
        //    LOD levels). The macro surface is already river-eroded, so its slopes lead into
        //    the macro valleys and the gullies drain into them.
        let g = if lim_child <= self.gully_wavelength_m {
            // Ridge/valley position in the hill relief (high = ridge), for the filter's rounding.
            let hill_scale = self.hill_slope * HILL_ROLLOFF_M / std::f64::consts::TAU * 2.0;
            let fade = (hl.full / hill_scale.max(1.0)).clamp(-1.0, 1.0);
            let unit = self.gully_unit_m;
            let pos = dir * (self.radius_m / unit);
            let cut = [lim, lim_parent, lim_child].map(|c| c / unit);
            erosion_filter::erosion_filter_sphere(
                pos,
                dir,
                grad3 * self.slope_gain,
                (m + hl.full) / unit,
                fade,
                &erosion_filter::FilterParams::default(),
                cut,
                DIFFUSION_WAVELENGTH_M / unit,
                self.gully_salt,
            )
        } else {
            erosion_filter::Banded::default()
        };
        // Valley depth follows local relief: full in mountains, gentle on lowland plains.
        let relief = self.relief(dir, m);
        let u = self.gully_unit_m * self.gully_strength * relief;
        // River corridors: flatten the detail along macro rivers and cut the channel.
        let (corridor, cut) = self.river_corridor(dir);
        let keep = 1.0 - RIVER_FLATTEN * corridor;
        HeightLevels {
            macro_h: m,
            h: coast_rule(m, (hl.level + g.at_level * u) * keep - cut),
            h_parent: coast_rule(m, (hl.parent + g.at_parent * u) * keep - cut),
            h_child: coast_rule(m, (hl.child + g.at_child * u) * keep - cut),
        }
    }

    /// Hill octaves (fBm + ridges in mountains) from two macro cells down to twice the
    /// first gully wavelength: the full sum plus the three LOD-banded sums.
    fn hills(&self, dir: DVec3, m: f64, cutoffs: [f64; 3]) -> Hills {
        let [lim, lim_parent, lim_child] = cutoffs;
        // Uneroded relief below the macro grid's resolution: a constant-slope spectrum
        // (each octave adds slopes of about HILL_SLOPE) scaled by how mountainous the macro
        // terrain is. The gully filter then carves valleys into this, as in Johansen's
        // design, which assumes a base already rough at the gully scale.
        let relief = self.relief(dir, m);
        // No ridged noise: (1 − |n|)² kinks along every zero-contour of the noise, drawing long
        // crisp arcs (seen at 25 km, s01). The gully filter sharpens ridges physically instead.
        let ridged_w = 0.0;
        let mut out = Hills::default();
        let mut wl = self.base_wavelength_m;
        let mut octave = 0u32;
        while wl >= self.gully_wavelength_m && octave < 12 {
            // Constant slope below HILL_ROLLOFF_M; amplitude saturates above it (relief is
            // limited by rock strength, and the macro surface already carries the long
            // wavelengths).
            // Band-pass spectrum: constant slope below HILL_ROLLOFF_M, amplitude falling as
            // 1/λ above it, so the hills hand the large scales over to the river-eroded
            // macro surface instead of adding drainage-blind relief there.
            let amp = self.hill_slope * wl / std::f64::consts::TAU / (1.0 + (wl / HILL_ROLLOFF_M).powi(2)) * relief;
            let p = dir * (self.radius_m / wl);
            let off = DVec3::new(octave as f64 * 31.7, octave as f64 * -17.3, octave as f64 * 9.1);
            let n = rg_noise::gradient3(p + off, self.salt.wrapping_add(octave));
            let v = if ridged_w > 0.0 {
                // Classic ridge (1 − |n|)²: sharp crests.
                let r = 1.0 - rg_noise::gradient3(p - off, self.ridge_salt.wrapping_add(octave)).abs();
                n * (1.0 - ridged_w) + (r * r * 2.0 - 0.8) * ridged_w
            } else {
                n
            };
            let c = v * amp;
            out.full += c;
            if wl >= lim_child {
                out.child += c;
            }
            if wl >= lim {
                out.level += c;
            }
            if wl >= lim_parent {
                out.parent += c;
            }
            wl *= 0.5;
            octave += 1;
        }
        out
    }

    /// How rugged the terrain is here, 0–1: local relief follows how steep the committed
    /// macro surface is (active, incising belts are steep), **not** elevation (a high
    /// plateau like Tibet is gentle; a low coastal range can be brutally steep). Measured
    /// on the macro field alone by 10 km differences, so it's smooth and identical at every
    /// LOD. Fades to nothing at sea level.
    pub fn relief(&self, dir: DVec3, m: f64) -> f64 {
        let east = DVec3::Z.cross(dir).normalize_or(DVec3::X);
        let north = dir.cross(east);
        let step = 10_000.0;
        let ang = step / self.radius_m;
        let me = self.macro_height((dir + east * ang).normalize());
        let mw = self.macro_height((dir - east * ang).normalize());
        let mn = self.macro_height((dir + north * ang).normalize());
        let ms = self.macro_height((dir - north * ang).normalize());
        let slope = (((me - mw) / (2.0 * step)).powi(2) + ((mn - ms) / (2.0 * step)).powi(2)).sqrt();
        let belt = smoothstep(RELIEF_SLOPE_LO, RELIEF_SLOPE_HI, slope);
        (HILL_LOWLAND + (1.0 - HILL_LOWLAND) * belt) * smoothstep(0.0, 300.0, m)
    }

    /// The full hill surface (macro + all hill octaves), whose slope steers the gullies.
    pub fn base_surface(&self, dir: DVec3) -> f64 {
        let m = self.macro_height(dir);
        if m <= 0.0 {
            return m;
        }
        if self.hill_slope == 0.0 {
            return m;
        }
        m + self.hills(dir, m, [0.0, 0.0, 0.0]).full
    }
}

/// Native-only tuning override from an environment variable (for parameter sweeps in the
/// container; unset in normal use, so generation stays deterministic).
fn tune(var: &str, default: f64) -> f64 {
    #[cfg(not(target_arch = "wasm32"))]
    if let Some(v) = std::env::var(var).ok().and_then(|v| v.parse().ok()) {
        return v;
    }
    let _ = var;
    default
}

#[derive(Default)]
struct Hills {
    full: f64,
    level: f64,
    parent: f64,
    child: f64,
}

/// Hill relief: slope contributed by each octave in mountains (×relief), and the share of
/// that kept on lowland plains.
/// 0 since s01: hills were unrouted fBm in the 4–13 km band the drainage network owns, and
/// put 13–16% of mountain land in closed depressions 150–210 m deep (probe_realism). Kept as
/// a knob (`RG_HILL_SLOPE`) for comparison; the gullies carry the sub-macro relief.
const HILL_SLOPE: f64 = 0.0;
const HILL_LOWLAND: f64 = 0.12;
/// Macro slope range mapping to relief 0 → 1 (tan of slope, measured over 20 km).
const RELIEF_SLOPE_LO: f64 = 0.004;
const RELIEF_SLOPE_HI: f64 = 0.04;
/// River corridors (see `river_corridor`): start, full strength, depth scaling and cap,
/// and how much of the detail is flattened inside a full corridor.
const RIVER_CORRIDOR_MIN_KM2: f64 = 3_000.0;
const RIVER_CORRIDOR_FULL_KM2: f64 = 30_000.0;
const RIVER_DEPTH_M_PER_SQRT_KKM2: f64 = 20.0;
const RIVER_DEPTH_MAX_M: f64 = 300.0;
const RIVER_FLATTEN: f64 = 0.85;
/// Wavelength above which hill amplitude stops growing, metres.
const HILL_ROLLOFF_M: f64 = 4000.0;
/// Gain on the slope the filter's masks see (gullies switch on at gentler macro slopes).
const SLOPE_GAIN: f64 = 1.0;
/// First (largest) gully wavelength, metres.
const GULLY_WAVELENGTH_M: f64 = 4000.0;
/// Central-difference step for the slope that steers the gullies, m. The macro surface is
/// only C0 across its cell edges (~13 km), so a short step flips the gully direction within
/// a few hundred metres there and draws a line along every cell edge; 3 km blends it.
const STEER_STEP_M: f64 = 3000.0;

/// Hillslope-diffusion length: gully octaves shorter than this are progressively smoothed
/// (soil creep), as on real hillslopes.
const DIFFUSION_WAVELENGTH_M: f64 = 150.0;

/// Valley/gully depth multiplier. The filter's octave amplitude is strength/cell_scale·λ
/// (≈ 0.31 λ at 1.0), i.e. constant slope per octave (gain 0.5 with lacunarity 2): at
/// 0.15 each octave contributes slopes up to ~0.3, and several stack in steep terrain.
const GULLY_STRENGTH: f64 = 2.0;

/// Heights of one point at three LOD levels; see [`TerrainSampler::levels`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeightLevels {
    pub macro_h: f64,
    pub h: f64,
    pub h_parent: f64,
    pub h_child: f64,
}

/// Coast rule: fade detail out within 150 m of sea level and never flip the sign, so the
/// coastline is exactly the committed macro contour at every zoom.
fn coast_rule(m: f64, d: f64) -> f64 {
    let d = d * smoothstep(0.0, 150.0, m.abs());
    if m >= 0.0 { (m + d).max(0.5 * m) } else { (m + d).min(0.5 * m) }
}

fn smoothstep(a: f64, b: f64, x: f64) -> f64 {
    let t = ((x - a) / (b - a)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn field() -> Arc<CubeField<f32>> {
        // Land in the +Z hemisphere rising to 3000 m, ocean in −Z.
        Arc::new(CubeField::from_fn(32, |d| (d.z * 4000.0) as f32))
    }

    #[test]
    fn macro_only_returns_the_macro_height_at_every_band() {
        let f = field();
        let s = TerrainSampler::new(f.clone(), 6_371_000.0, 7).with_macro_only();
        for i in 0..500 {
            let t = i as f64 * 0.731;
            let d = DVec3::new(t.cos(), t.sin(), (i as f64 * 0.004) - 1.0).normalize();
            let m = s.macro_height(d);
            for wl in [0.0, 40.0, 1000.0, 50_000.0] {
                for l in [s.levels(d, wl), s.levels_with_slope(d, wl, DVec3::X * 0.3)] {
                    assert_eq!([l.macro_h, l.h, l.h_parent, l.h_child], [m; 4], "detail at {d}, wl {wl}");
                }
            }
            assert!(s.network_patch(d, 10_000.0, 40.0).is_none());
        }
    }

    #[test]
    fn coast_never_moves() {
        let f = field();
        let s = TerrainSampler::new(f.clone(), 6_371_000.0, 7);
        for i in 0..4000 {
            let t = i as f64 * 0.731;
            let d = DVec3::new(t.cos(), t.sin(), (i as f64 * 0.0007) - 1.4).normalize();
            let (h, m) = s.height(d, 0.0);
            assert!(h.signum() == m.signum() || m == 0.0, "coast moved at {d}: macro {m}, h {h}");
        }
    }

    #[test]
    fn band_limited_tiles_converge() {
        // Coarser evaluations differ from the full surface by a small, shrinking amount.
        let f = field();
        let s = TerrainSampler::new(f.clone(), 6_371_000.0, 7);
        let d = DVec3::new(0.3, 0.2, 0.9).normalize();
        let full = s.height(d, 0.0).0;
        let e1 = (s.height(d, 1000.0).0 - full).abs();
        let e2 = (s.height(d, 100.0).0 - full).abs();
        assert!(e2 <= e1 + 1e-9 || e2 < 5.0, "e1 {e1} e2 {e2}");
        assert!(e2 < 30.0, "dropping sub-100 m octaves changed height by {e2} m");
    }

    #[test]
    fn deterministic() {
        let f = field();
        let a = TerrainSampler::new(f.clone(), 6_371_000.0, 7);
        let b = TerrainSampler::new(f.clone(), 6_371_000.0, 7);
        let d = DVec3::new(0.1, -0.4, 0.8).normalize();
        assert_eq!(a.height(d, 20.0), b.height(d, 20.0));
    }
}
