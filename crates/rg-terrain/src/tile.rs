//! Render tiles: 64 × 64 cells (65² vertices) plus a skirt, generated on the CPU.

use crate::TerrainSampler;
use glam::DVec3;
use rg_geo::TileId;

/// Cells per render-tile edge. Data tiles are 256 (ARCHITECTURE §4); 256 = 4 × 64, so a
/// render tile at level L lines up with data tile level L − 2.
pub const TILE_CELLS: u32 = 64;
const N: usize = TILE_CELLS as usize + 1; // vertices per edge

/// Everything the renderer needs for one tile. Positions are **relative to `center`** and
/// in metres, so they fit f32 without jitter. `center` is f64 ECEF metres.
pub struct TileData {
    pub id: TileId,
    pub center: DVec3,
    /// 65² grid vertices (row-major, v then u), then 4 × 65 skirt vertices.
    pub positions: Vec<[f32; 3]>,
    pub normals: Vec<[f32; 3]>,
    /// True height above sea level (m), negative under the sea. Used for colouring.
    pub heights: Vec<f32>,
    /// Radius of a sphere around `center` that contains the tile.
    pub bound_radius_m: f64,
    pub max_height_m: f64,
}

#[derive(Clone)]
pub struct TileGenerator {
    pub sampler: TerrainSampler,
    /// Vertical exaggeration for display.
    pub exaggeration: f64,
}

impl TileGenerator {
    /// Ground cell size of a render tile at `level`, metres.
    pub fn cell_size_m(&self, level: u8) -> f64 {
        (std::f64::consts::FRAC_PI_2 * self.sampler.radius_m) / (TILE_CELLS as f64 * (1u64 << level) as f64)
    }

    pub fn generate(&self, id: TileId) -> TileData {
        let r = self.sampler.radius_m;
        let cell = self.cell_size_m(id.level);
        // Band-limit: only octaves this tile can represent (≥ 2 cells per wavelength).
        let min_wl = 2.0 * cell;
        let (u0, v0, u1, v1) = id.uv_bounds();
        let du = (u1 - u0) / TILE_CELLS as f64;
        let dv = (v1 - v0) / TILE_CELLS as f64;

        // Sample a (N+2)² grid: the extra ring (apron) lets normals use central differences
        // that match the neighbouring tile's, so lighting has no seams. Apron samples may
        // lie past the face edge; uv_to_dir extrapolates smoothly there, which is fine for
        // a normal estimate.
        let m = N + 2;
        let mut world = vec![DVec3::ZERO; m * m];
        let mut true_h = vec![0f32; m * m];
        let mut max_h = f64::MIN;
        for j in 0..m {
            for i in 0..m {
                let u = u0 + (i as f64 - 1.0) * du;
                let v = v0 + (j as f64 - 1.0) * dv;
                let dir = id.face.uv_to_dir(u, v);
                let (h, _) = self.sampler.height(dir, min_wl);
                // Oceans render at sea level; depth only shows through colour.
                let shown = h.max(0.0) * self.exaggeration;
                world[j * m + i] = dir * (r + shown);
                true_h[j * m + i] = h as f32;
                max_h = max_h.max(shown);
            }
        }

        let center_dir = id.face.uv_to_dir((u0 + u1) * 0.5, (v0 + v1) * 0.5);
        let center = center_dir * r;
        let mut positions = Vec::with_capacity(N * N + 4 * N);
        let mut normals = Vec::with_capacity(N * N + 4 * N);
        let mut heights = Vec::with_capacity(N * N + 4 * N);
        let mut bound = 0f64;
        for j in 0..N {
            for i in 0..N {
                let k = (j + 1) * m + (i + 1);
                let p = world[k];
                let tu = world[k + 1] - world[k - 1];
                let tv = world[k + m] - world[k - m];
                let n = tu.cross(tv).normalize_or(p.normalize());
                positions.push((p - center).as_vec3().to_array());
                normals.push(n.as_vec3().to_array());
                heights.push(true_h[k]);
                bound = bound.max((p - center).length());
            }
        }

        // Skirt: a copy of the border pushed down, hiding cracks against coarser
        // neighbours. Order: bottom (j=0), right (i=N-1), top (j=N-1), left (i=0).
        let drop = cell * 3.0 + 50.0;
        for (i, j) in skirt_border() {
            let k = j * N + i;
            let p = DVec3::from(positions[k].map(f64::from)) + center;
            let p = p - p.normalize() * drop;
            positions.push((p - center).as_vec3().to_array());
            normals.push(normals[k]);
            heights.push(heights[k]);
        }

        TileData { id, center, positions, normals, heights, bound_radius_m: bound + drop, max_height_m: max_h }
    }
}

/// (i, j) of the border vertices in skirt order.
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
    let mut idx = Vec::with_capacity((TILE_CELLS * TILE_CELLS * 6) as usize + 4 * N * 6);
    for j in 0..n - 1 {
        for i in 0..n - 1 {
            let a = j * n + i;
            let (b, c, d) = (a + 1, a + n, a + n + 1);
            idx.extend_from_slice(&[a, b, d, a, d, c]);
        }
    }
    // Skirt quads between each border vertex and its dropped copy. The border is walked
    // counter-clockwise seen from outside, so (top edge, dropped) forms outward faces.
    let border: Vec<u32> = skirt_border().map(|(i, j)| (j * N + i) as u32).collect();
    let base = n * n;
    for s in 0..border.len() - 1 {
        let (a, b) = (border[s], border[s + 1]);
        let (a2, b2) = (base + s as u32, base + s as u32 + 1);
        // Both windings: skirts are seen from either side depending on the neighbour.
        idx.extend_from_slice(&[a, a2, b, b, a2, b2, a, b, a2, b, b2, a2]);
    }
    idx
}

#[cfg(test)]
mod tests {
    use super::*;
    use rg_geo::{CubeField, CubeFace};

    #[test]
    fn tile_has_expected_shape_and_bounds() {
        let f = std::sync::Arc::new(CubeField::from_fn(32, |d| (d.z * 3000.0) as f32));
        let g = TileGenerator { sampler: TerrainSampler::new(f.clone(), 6_371_000.0, 1), exaggeration: 1.0 };
        let t = g.generate(TileId { face: CubeFace::PosZ, level: 3, x: 2, y: 5 });
        assert_eq!(t.positions.len(), N * N + 4 * N);
        assert!(t.bound_radius_m > 0.0 && t.bound_radius_m < 2.0e6);
        let max_idx = *tile_indices().iter().max().unwrap() as usize;
        assert!(max_idx < t.positions.len());
    }

    #[test]
    fn shared_edges_match_between_neighbours() {
        // Same-level neighbours must produce identical edge positions (no cracks).
        let f = std::sync::Arc::new(CubeField::from_fn(32, |d| (d.z * 3000.0) as f32));
        let g = TileGenerator { sampler: TerrainSampler::new(f.clone(), 6_371_000.0, 1), exaggeration: 10.0 };
        let a = g.generate(TileId { face: CubeFace::PosX, level: 4, x: 3, y: 7 });
        let b = g.generate(TileId { face: CubeFace::PosX, level: 4, x: 4, y: 7 });
        for j in 0..N {
            let pa = DVec3::from(a.positions[j * N + (N - 1)].map(f64::from)) + a.center;
            let pb = DVec3::from(b.positions[j * N].map(f64::from)) + b.center;
            assert!((pa - pb).length() < 0.05, "edge gap {} m", (pa - pb).length());
        }
    }
}
