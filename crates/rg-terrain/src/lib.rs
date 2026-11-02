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

pub mod tile;

pub use tile::{tile_indices, TileData, TileGenerator, TILE_CELLS};

use glam::DVec3;
use rg_geo::CubeField;
use std::sync::Arc;

/// Version of the detail function. Bump on any change that alters heights; stored sculpt
/// deltas (M6) record the version they were made against.
pub const DETAIL_GENERATOR_VERSION: u32 = 0; // 0 = not persistent yet (no sculpt deltas exist)

/// Samples the full-resolution terrain surface. Cheap to clone (the field is shared), and
/// `Send`, so native builds can generate tiles on worker threads.
#[derive(Clone)]
pub struct TerrainSampler {
    pub macro_h: Arc<CubeField<f32>>,
    pub radius_m: f64,
    salt: u32,
    ridge_salt: u32,
    /// Wavelength of the first detail octave (m): about two macro cells.
    base_wavelength_m: f64,
    /// Finest wavelength ever generated (m).
    floor_wavelength_m: f64,
}

/// Amplitude ratio between successive octaves. 0.58 ≈ a Hurst exponent of 0.8, close to
/// measured spectra of real terrain (rougher than the textbook 0.5).
const GAIN: f64 = 0.58;

impl TerrainSampler {
    pub fn new(macro_h: Arc<CubeField<f32>>, radius_m: f64, seed: u64) -> Self {
        let macro_cell_m = macro_h.spacing_rad() * radius_m;
        Self {
            macro_h,
            radius_m,
            salt: rg_noise::salt(seed, 0x7E44),
            ridge_salt: rg_noise::salt(seed, 0x7E45),
            base_wavelength_m: 2.0 * macro_cell_m,
            floor_wavelength_m: 12.0,
        }
    }

    /// Macro (committed) height at a unit direction, metres above sea level.
    pub fn macro_height(&self, dir: DVec3) -> f64 {
        self.macro_h.sample(dir) as f64
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
        let m = self.macro_height(dir);
        // Oceans render flat at sea level, so skip detail there entirely (~70% of the
        // planet). The coast fade below makes this continuous at m = 0.
        if m <= 0.0 {
            return HeightLevels { macro_h: m, h: m, h_parent: m, h_child: m };
        }
        let lim = min_wavelength_m.max(self.floor_wavelength_m);
        let lim_parent = (2.0 * min_wavelength_m).max(self.floor_wavelength_m);
        let lim_child = (0.5 * min_wavelength_m).max(self.floor_wavelength_m);

        // Context: rugged in mountains, gentle on lowlands.
        let amp0 = 90.0 + 0.16 * m;
        let ridged_w = smoothstep(400.0, 2800.0, m);

        let (mut d, mut d_parent, mut d_child) = (0.0, 0.0, 0.0);
        let mut wl = self.base_wavelength_m;
        let mut amp = amp0;
        let mut octave = 0u32;
        while wl >= lim_child && octave < 24 {
            let p = dir * (self.radius_m / wl);
            let off = DVec3::new(octave as f64 * 31.7, octave as f64 * -17.3, octave as f64 * 9.1);
            let n = rg_noise::gradient3(p + off, self.salt.wrapping_add(octave));
            // Ridged octaves for the first few (large) scales in mountains: crests and
            // valleys instead of blobs.
            let v = if ridged_w > 0.0 && octave < 6 {
                // Classic ridge (1 − |n|)²: a sharp crest where n = 0, smooth everywhere
                // else. (An earlier (1 − 1.6|n|)² creased along |n| = 0.625: scalloped arcs.)
                let r = 1.0 - rg_noise::gradient3(p - off, self.ridge_salt.wrapping_add(octave)).abs();
                n * (1.0 - ridged_w) + (r * r * 2.0 - 0.8) * ridged_w
            } else {
                n
            };
            let c = v * amp;
            d_child += c;
            if wl >= lim {
                d += c;
            }
            if wl >= lim_parent {
                d_parent += c;
            }
            amp *= GAIN;
            wl *= 0.5;
            octave += 1;
        }
        HeightLevels { macro_h: m, h: coast_rule(m, d), h_parent: coast_rule(m, d_parent), h_child: coast_rule(m, d_child) }
    }
}

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
