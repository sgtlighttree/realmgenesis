//! The baked surface gradient used for per-pixel lighting on the globe (R2,
//! `docs/design/raster-renderer.md`).
//!
//! Tiles carry heights per vertex, but their vertex spacing is coarse from orbit, so vertex
//! normals would smear the macro relief. Instead the fragment shader reads the slope of the
//! displayed surface from a texture: for every sample of a cube-face raster, the gradient of
//! `TerrainSampler::macro_height` (m per m) in the local (east, north) frame.
//!
//! The frame matches the shader exactly: `east = normalize(Z × up)`, `north = up × east`.
//! It's undefined at the two poles, where east falls back to +X (a pixel or so).

use crate::TerrainSampler;
use glam::DVec3;
use rg_geo::{CubeFace, CubeField};

/// Raster size for a macro field of `macro_res` cells per face edge: twice the macro (the
/// display surface is a bicubic, so its slope varies inside a macro cell), within 512..2048.
pub fn gradient_res(macro_res: u32) -> u32 {
    (macro_res * 2).clamp(512, 2048)
}

/// The local (east, north) frame at a unit direction, as the shader computes it.
pub fn east_north(up: DVec3) -> (DVec3, DVec3) {
    let east = DVec3::Z.cross(up).normalize_or(DVec3::X);
    (east, up.cross(east))
}

/// Bake the gradient of the displayed surface at `res` cells per face edge (corner-aligned
/// like every `CubeField`). Each sample is `[d h / d east, d h / d north]`, unitless.
///
/// Heights are evaluated once per sample, then differentiated against the four grid
/// neighbours in 3D (so the cube warp's skew is handled exactly). Samples on a face edge use
/// heights just across the edge, taken from the sampler directly, so the two faces agree.
pub fn bake_gradient(sampler: &TerrainSampler, res: u32) -> CubeField<[f32; 2]> {
    let n = (res + 1) as usize;
    let radius = sampler.radius_m;
    let mut out: CubeField<[f32; 2]> = CubeField::new(res);
    let face_len = n * n;
    // One face at a time keeps the height table at (res + 1)² floats (17 MB at 2048).
    for face in CubeFace::ALL {
        let uv = |k: i64| -1.0 + 2.0 * k as f64 / res as f64;
        let mut h = vec![0f32; face_len];
        par_rows(&mut h, n, |j, row| {
            for (i, v) in row.iter_mut().enumerate() {
                *v = sampler.macro_height(face.uv_to_dir(uv(i as i64), uv(j as i64))) as f32;
            }
        });
        let start = face as usize * face_len;
        let height = |i: i64, j: i64| -> (DVec3, f64) {
            let dir = face.uv_to_dir(uv(i), uv(j));
            if (0..n as i64).contains(&i) && (0..n as i64).contains(&j) {
                (dir, h[j as usize * n + i as usize] as f64)
            } else {
                // Past the face edge: uv_to_dir extrapolates smoothly onto the next face.
                (dir.normalize(), sampler.macro_height(dir.normalize()))
            }
        };
        par_rows(&mut out.data[start..start + face_len], n, |j, row| {
            let j = j as i64;
            for (i, g) in row.iter_mut().enumerate() {
                let i = i as i64;
                let up = face.uv_to_dir(uv(i), uv(j));
                let (pe, he) = height(i + 1, j);
                let (pw, hw) = height(i - 1, j);
                let (pn, hn) = height(i, j + 1);
                let (ps, hs) = height(i, j - 1);
                *g = tangent_gradient(up, (pe - pw) * radius, he - hw, (pn - ps) * radius, hn - hs);
            }
        });
    }
    out
}

/// The tangent-plane gradient at `up` from two height differences `da` and `db` taken
/// along the (not necessarily orthogonal) steps `a` and `b` (metres), in (east, north).
fn tangent_gradient(up: DVec3, a: DVec3, da: f64, b: DVec3, db: f64) -> [f32; 2] {
    let a = a - up * a.dot(up);
    let b = b - up * b.dot(up);
    // g = α a + β b with g·a = da and g·b = db: a 2×2 Gram system.
    let (aa, ab, bb) = (a.dot(a), a.dot(b), b.dot(b));
    let det = aa * bb - ab * ab;
    if det.abs() < 1e-12 {
        return [0.0, 0.0];
    }
    let alpha = (da * bb - db * ab) / det;
    let beta = (db * aa - da * ab) / det;
    let g = a * alpha + b * beta;
    let (east, north) = east_north(up);
    [g.dot(east) as f32, g.dot(north) as f32]
}

/// Run `f(row index, row)` over the rows of a row-major `data` with rows of `width`,
/// spread over the available cores natively (one thread on the web).
fn par_rows<T: Send>(data: &mut [T], width: usize, f: impl Fn(usize, &mut [T]) + Sync) {
    #[cfg(not(target_arch = "wasm32"))]
    {
        let threads = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4);
        let rows = data.len() / width;
        let per = rows.div_ceil(threads).max(1);
        std::thread::scope(|s| {
            for (c, chunk) in data.chunks_mut(per * width).enumerate() {
                let f = &f;
                s.spawn(move || {
                    for (r, row) in chunk.chunks_mut(width).enumerate() {
                        f(c * per + r, row);
                    }
                });
            }
        });
    }
    #[cfg(target_arch = "wasm32")]
    for (r, row) in data.chunks_mut(width).enumerate() {
        f(r, row);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    /// h = 3000 · z (metres) on a sphere of radius R: the tangent gradient is
    /// 3000 / R · (Z − (Z·up) up), whose north component is 3000/R · cos(lat), east 0.
    #[test]
    fn gradient_matches_an_analytic_field() {
        let r = 6_371_000.0;
        let f = Arc::new(CubeField::from_fn(64, |d| (d.z * 3000.0) as f32));
        let s = TerrainSampler::new(f, r, 1).with_macro_only();
        let g = bake_gradient(&s, 64);
        let mut worst: f64 = 0.0;
        for (k, _, _, _, dir) in g.iter_samples() {
            let lat = dir.z.clamp(-1.0, 1.0).asin();
            if lat.abs() > 1.4 {
                continue; // the frame degenerates at the poles
            }
            let expect_n = 3000.0 / r * lat.cos();
            let [ge, gn] = g.data[k];
            let err = ((ge as f64).powi(2) + (gn as f64 - expect_n).powi(2)).sqrt() / (3000.0 / r);
            worst = worst.max(err);
        }
        assert!(worst < 0.02, "worst relative error {worst}");
    }

    #[test]
    fn face_edge_twins_agree() {
        let r = 6_371_000.0;
        let f = Arc::new(CubeField::from_fn(32, |d| (d.x * 2000.0 + d.y * d.z * 4000.0) as f32));
        let s = TerrainSampler::new(f, r, 1).with_macro_only();
        let g = bake_gradient(&s, 32);
        // Group samples by quantised direction; twins must have nearly equal gradients.
        let mut by_pos: std::collections::HashMap<[i64; 3], [f32; 2]> = Default::default();
        let mut worst: f32 = 0.0;
        for (k, _, _, _, d) in g.iter_samples() {
            let key = [(d.x * 1e6).round() as i64, (d.y * 1e6).round() as i64, (d.z * 1e6).round() as i64];
            if let Some(o) = by_pos.insert(key, g.data[k]) {
                let scale = 4000.0 / r as f32;
                worst = worst.max(((o[0] - g.data[k][0]).abs() + (o[1] - g.data[k][1]).abs()) / scale);
            }
        }
        assert!(worst < 0.1, "twin mismatch {worst} (relative)");
    }
}
