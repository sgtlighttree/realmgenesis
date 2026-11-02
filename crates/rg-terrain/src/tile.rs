//! Render tiles: 64 × 64 cells (65² vertices) plus a skirt, generated on the CPU.
//!
//! Vertices are **unrelieved**: positions lie on the sea-level sphere, and the shader
//! raises them by `max(h, 0) × exaggeration` along the local up. That way changing the
//! exaggeration (it fades as you zoom in) never invalidates the cache.

use crate::TerrainSampler;
use glam::DVec3;
use rg_geo::TileId;

/// Cells per render-tile edge. Data tiles are 256 (ARCHITECTURE §4); 256 = 4 × 64, so a
/// render tile at level L lines up with data tile level L − 2.
pub const TILE_CELLS: u32 = 64;
const N: usize = TILE_CELLS as usize + 1; // vertices per edge
/// Vertices per tile including the skirt.
pub const TILE_VERTS: usize = N * N + 4 * N;

/// Everything the renderer needs for one tile. Positions are **relative to `center`**
/// (metres) so they fit f32 without jitter; `center` is f64 ECEF metres on the sphere.
pub struct TileData {
    pub id: TileId,
    pub center: DVec3,
    /// Sea-level sphere positions: 65² grid vertices (row-major, v then u), then the skirt.
    pub positions: Vec<[f32; 3]>,
    /// Surface normals at exaggeration 1 (the shader steepens them by the exaggeration).
    pub normals: Vec<[f32; 3]>,
    /// Height above sea level (m, negative under the sea) at this tile's level.
    pub heights: Vec<f32>,
    /// The parent tile's *mesh* surface under each vertex (parent-level heights, linearly
    /// interpolated along the parent's triangles): the geomorph target.
    pub parent_heights: Vec<f32>,
    /// 1.0 for skirt vertices, 0.0 for grid vertices.
    pub skirt: Vec<f32>,
    /// Radius around `center` containing the unrelieved patch.
    pub patch_radius_m: f64,
    /// Highest point (m above sea level, ≥ 0).
    pub max_height_m: f64,
    /// Largest height gap between this tile and its children's surface (m, unexaggerated).
    /// Drives the split decision: flat or ocean tiles stop splitting early.
    pub child_error_m: f64,
    /// Largest gap between this tile and its parent's surface (m), sizes the skirt.
    pub parent_error_m: f64,
}

#[derive(Clone)]
pub struct TileGenerator {
    pub sampler: TerrainSampler,
}

impl TileGenerator {
    /// Ground cell size of a render tile at `level`, metres.
    pub fn cell_size_m(&self, level: u8) -> f64 {
        cell_size_m(level, self.sampler.radius_m)
    }

    /// Pure function of (sampler, id): safe to run on any thread.
    pub fn generate(&self, id: TileId) -> TileData {
        let r = self.sampler.radius_m;
        let cell = self.cell_size_m(id.level);
        // Band-limit: only octaves this tile can represent (≥ 2 cells per wavelength).
        let min_wl = 2.0 * cell;
        let (u0, v0, u1, v1) = id.uv_bounds();
        let du = (u1 - u0) / TILE_CELLS as f64;
        let dv = (v1 - v0) / TILE_CELLS as f64;

        // Sample an (N+2)² grid. The extra ring (apron) gives normals central differences
        // that match the neighbouring tile's, so lighting has no seams. Apron samples past
        // a face edge use uv_to_dir's smooth extrapolation, fine for a normal estimate.
        let m = N + 2;
        // Heights on the (N+2)² grid (one apron ring for normals). The gully filter's slope
        // comes from fixed 250 m finite differences inside `levels`, not from this grid, so
        // every LOD level sees the same slope and a child's parent-level heights match the
        // parent tile exactly (needed for seamless geomorphing).
        let mut surf = vec![DVec3::ZERO; m * m]; // true (exaggeration 1) surface
        let mut dirs = vec![DVec3::ZERO; m * m];
        let mut lv = vec![crate::HeightLevels { macro_h: 0.0, h: 0.0, h_parent: 0.0, h_child: 0.0 }; m * m];
        // With the drained detail network, one patch covers the whole tile (+ apron).
        let patch = {
            let c = id.face.uv_to_dir((u0 + u1) * 0.5, (v0 + v1) * 0.5);
            let corner = id.face.uv_to_dir(u0 - du, v0 - dv);
            let corner2 = id.face.uv_to_dir(u1 + du, v1 + dv);
            let rad = (c - corner).length().max((c - corner2).length()) * r * 1.05;
            self.sampler.network_patch(c, rad, min_wl)
        };
        for j in 0..m {
            for i in 0..m {
                let u = u0 + (i as f64 - 1.0) * du;
                let v = v0 + (j as f64 - 1.0) * dv;
                let dir = id.face.uv_to_dir(u, v);
                let l = match &patch {
                    Some(p) => self.sampler.levels_in_patch(p, dir),
                    None => self.sampler.levels(dir, min_wl),
                };
                let k = j * m + i;
                surf[k] = dir * (r + l.h.max(0.0));
                dirs[k] = dir;
                lv[k] = l;
            }
        }

        let center_dir = id.face.uv_to_dir((u0 + u1) * 0.5, (v0 + v1) * 0.5);
        let center = center_dir * r;
        let mut out = TileData {
            id,
            center,
            positions: Vec::with_capacity(TILE_VERTS),
            normals: Vec::with_capacity(TILE_VERTS),
            heights: Vec::with_capacity(TILE_VERTS),
            parent_heights: Vec::with_capacity(TILE_VERTS),
            skirt: Vec::with_capacity(TILE_VERTS),
            patch_radius_m: 0.0,
            max_height_m: 0.0,
            child_error_m: 0.0,
            parent_error_m: 0.0,
        };
        for j in 0..N {
            for i in 0..N {
                let k = (j + 1) * m + (i + 1);
                let tu = surf[k + 1] - surf[k - 1];
                let tv = surf[k + m] - surf[k - m];
                let n = if lv[k].h <= 0.0 { dirs[k] } else { tu.cross(tv).normalize_or(dirs[k]) };
                let p = dirs[k] * r;
                out.positions.push((p - center).as_vec3().to_array());
                out.normals.push(n.as_vec3().to_array());
                out.heights.push(lv[k].h as f32);
                out.parent_heights.push(lv[k].h_parent as f32);
                out.skirt.push(0.0);
                out.patch_radius_m = out.patch_radius_m.max((p - center).length());
                out.max_height_m = out.max_height_m.max(lv[k].h);
                out.child_error_m = out.child_error_m.max((lv[k].h_child.max(0.0) - lv[k].h.max(0.0)).abs());
                out.parent_error_m = out.parent_error_m.max((lv[k].h.max(0.0) - lv[k].h_parent.max(0.0)).abs());
            }
        }
        // Macro undersampling: a coarse tile's cells can be wider than the macro field's,
        // so mountains between vertices get cut off. Compare the macro height at each cell
        // centre with the average of its corners and fold that into the split error.
        for j in 0..N - 1 {
            for i in 0..N - 1 {
                let k = (j + 1) * m + (i + 1);
                let corners = [lv[k].macro_h, lv[k + 1].macro_h, lv[k + m].macro_h, lv[k + m + 1].macro_h];
                let avg = corners.iter().map(|h| h.max(0.0)).sum::<f64>() / 4.0;
                let mid = id.face.uv_to_dir(u0 + (i as f64 + 0.5) * du, v0 + (j as f64 + 0.5) * dv);
                let e = (self.sampler.macro_height(mid).max(0.0) - avg).abs();
                out.child_error_m = out.child_error_m.max(e);
            }
        }
        // Morph target: the *parent mesh's* surface at each vertex. Even (i, j) vertices
        // coincide with parent vertices, so they keep the parent-level height. Odd ones lie
        // on parent edges (average of the two ends) or on the parent quad's a–d diagonal,
        // which is how quads are split into triangles (average of those corners). Morphing
        // fully to this makes a fine tile identical to its coarser neighbour at the border.
        let hp: Vec<f32> = out.parent_heights.clone();
        let at = |i: usize, j: usize| hp[j * N + i];
        for j in 0..N {
            for i in 0..N {
                let v = match (i % 2, j % 2) {
                    (0, 0) => at(i, j),
                    (1, 0) => 0.5 * (at(i - 1, j) + at(i + 1, j)),
                    (0, 1) => 0.5 * (at(i, j - 1) + at(i, j + 1)),
                    _ => 0.5 * (at(i - 1, j - 1) + at(i + 1, j + 1)),
                };
                out.parent_heights[j * N + i] = v;
            }
        }
        // Skirt: copies of the border. The shader drops them below the surface by a
        // per-tile depth, hiding cracks against coarser neighbours.
        for (i, j) in skirt_border() {
            let k = j * N + i;
            out.positions.push(out.positions[k]);
            out.normals.push(out.normals[k]);
            out.heights.push(out.heights[k]);
            out.parent_heights.push(out.parent_heights[k]);
            out.skirt.push(1.0);
        }
        out
    }
}

/// Ground cell size of a render tile at `level` on a sphere of `radius_m`.
pub fn cell_size_m(level: u8, radius_m: f64) -> f64 {
    (std::f64::consts::FRAC_PI_2 * radius_m) / (TILE_CELLS as f64 * (1u64 << level) as f64)
}

/// (i, j) of the border vertices in skirt order: bottom (j=0), right, top, left.
fn skirt_border() -> impl Iterator<Item = (usize, usize)> {
    let e = N - 1;
    (0..N)
        .map(move |i| (i, 0))
        .chain((0..N).map(move |j| (e, j)))
        .chain((0..N).rev().map(move |i| (i, e)))
        .chain((0..N).rev().map(move |j| (0, j)))
}

/// The shared index buffer for every tile: the 64² grid plus the skirt strip.
pub fn tile_indices() -> Vec<u32> {
    let n = N as u32;
    let mut idx = Vec::with_capacity((TILE_CELLS * TILE_CELLS * 6) as usize + 4 * N * 12);
    for j in 0..n - 1 {
        for i in 0..n - 1 {
            let a = j * n + i;
            let (b, c, d) = (a + 1, a + n, a + n + 1);
            idx.extend_from_slice(&[a, b, d, a, d, c]);
        }
    }
    let border: Vec<u32> = skirt_border().map(|(i, j)| (j * N + i) as u32).collect();
    let base = n * n;
    for s in 0..border.len() - 1 {
        let (a, b) = (border[s], border[s + 1]);
        let (a2, b2) = (base + s as u32, base + s as u32 + 1);
        // Both windings: a skirt is seen from either side depending on the neighbour.
        idx.extend_from_slice(&[a, a2, b, b, a2, b2, a, b, a2, b, b2, a2]);
    }
    idx
}

#[cfg(test)]
mod tests {
    use super::*;
    use rg_geo::{CubeFace, CubeField};
    use std::sync::Arc;

    fn generator() -> TileGenerator {
        let f = Arc::new(CubeField::from_fn(32, |d| (d.z * 3000.0) as f32));
        TileGenerator { sampler: TerrainSampler::new(f, 6_371_000.0, 1) }
    }

    #[test]
    fn tile_has_expected_shape_and_bounds() {
        let t = generator().generate(TileId { face: CubeFace::PosZ, level: 3, x: 2, y: 5 });
        assert_eq!(t.positions.len(), TILE_VERTS);
        assert_eq!(t.heights.len(), TILE_VERTS);
        assert!(t.patch_radius_m > 0.0 && t.patch_radius_m < 2.0e6);
        assert!(t.child_error_m >= 0.0 && t.max_height_m > 0.0);
        let max_idx = *tile_indices().iter().max().unwrap() as usize;
        assert!(max_idx < TILE_VERTS);
    }

    #[test]
    fn shared_edges_match_between_neighbours() {
        // Same-level neighbours must produce identical edge vertices (no cracks).
        let g = generator();
        let a = g.generate(TileId { face: CubeFace::PosX, level: 4, x: 3, y: 7 });
        let b = g.generate(TileId { face: CubeFace::PosX, level: 4, x: 4, y: 7 });
        for j in 0..N {
            let ka = j * N + (N - 1);
            let kb = j * N;
            let pa = DVec3::from(a.positions[ka].map(f64::from)) + a.center;
            let pb = DVec3::from(b.positions[kb].map(f64::from)) + b.center;
            assert!((pa - pb).length() < 0.05, "edge gap {} m", (pa - pb).length());
            assert!((a.heights[ka] - b.heights[kb]).abs() < 1e-3);
        }
    }

    #[test]
    fn ocean_tiles_are_flat_and_need_no_split() {
        let t = generator().generate(TileId { face: CubeFace::NegZ, level: 5, x: 10, y: 10 });
        assert_eq!(t.max_height_m, 0.0);
        assert_eq!(t.child_error_m, 0.0);
    }
}
