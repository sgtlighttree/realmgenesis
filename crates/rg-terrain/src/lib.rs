//! The displayed terrain surface and its render tiles (R2, `docs/design/raster-renderer.md`).
//!
//! The surface is the committed macro height field, smoothly interpolated: a Catmull-Rom
//! bicubic through the samples, with the coastline kept where the drained interpolant puts it
//! (s02: bilinear facets showed at high relief). Lighting detail comes from the baked gradient
//! (`gradient.rs`). Sharper detail is the job of regional bakes (R4), not of this crate.
//!
//! History: until s03 this crate also streamed procedural detail (hills, gully filter, then
//! the drained detail network). It drew artifacts from orbit and stalled the UI, and was
//! removed after Matt signed off on R2; see git history before the s03 removal commit.

pub mod gradient;
pub mod tile;

pub use tile::{tile_indices, TileData, TileGenerator, TILE_CELLS};

use glam::DVec3;
use rg_geo::CubeField;
use std::sync::Arc;

/// Samples the displayed surface. Cheap to clone (the field is shared) and `Send`, so tiles
/// and bakes can run on worker threads.
#[derive(Clone)]
pub struct TerrainSampler {
    pub macro_h: Arc<CubeField<f32>>,
    pub radius_m: f64,
}

impl TerrainSampler {
    pub fn new(macro_h: Arc<CubeField<f32>>, radius_m: f64) -> Self {
        Self { macro_h, radius_m }
    }

    /// Height of the displayed surface at a unit direction, metres above sea level.
    /// Bicubic, except where the bicubic and the drained (routing) interpolant disagree about
    /// land vs sea: there the drained one wins, so coastlines match the rivers and lakes.
    pub fn macro_height(&self, dir: DVec3) -> f64 {
        let drained = self.macro_h.sample_drained(dir);
        let smooth = self.macro_h.sample_smooth(dir);
        if (smooth > 0.0) == (drained > 0.0) { smooth } else { drained }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn passes_through_samples_and_keeps_the_coast() {
        let f = Arc::new(CubeField::from_fn(32, |d| ((d.x * 7.0).sin() * 2000.0 + d.z * 500.0) as f32));
        let s = TerrainSampler::new(f.clone(), 6_371_000.0);
        for (k, _, _, _, d) in f.iter_samples().step_by(13) {
            assert!((s.macro_height(d) - f.data[k] as f64).abs() < 1e-2);
        }
        // Between samples, land/sea always agrees with the drained interpolant.
        for i in 0..2000 {
            let t = i as f64 * 0.7311;
            let d = DVec3::new(t.sin(), (t * 1.3).cos(), (t * 0.7).sin()).normalize();
            assert_eq!(s.macro_height(d) > 0.0, f.sample_drained(d) > 0.0);
        }
    }
}
