//! R1 spike (s03): erode on a **jittered** grid, so rivers aren't bound to the lattice's
//! 8 directions, then resample the results back onto the regular cube-face raster.
//!
//! Why: on a regular 8-connected grid, steepest descent can only step along the axes or the
//! diagonals, so trunks run straight with 45° turns (s02: "griddy", "mechanistic"). Rho8
//! randomises the diagonal length but still picks among the same 8 directions. Moving every
//! node a little (Cordonnier et al. 2016 use irregular grids for the same reason) makes the
//! true distances, and so the steepest neighbour, vary node by node. Stored on the jittered
//! positions and interpolated back, a river's path is an irregular polyline.
//!
//! Rules (from the s03 review, `docs/design/v2-direction-s03.md`):
//! - jitter only the erosion grid, never tectonics (that would move coasts);
//! - key the offsets on the node's integer sample index (not on `dir`, which depends on libm);
//! - keep face-edge nodes unjittered (twins on two faces must stay one place);
//! - recompute `min_edge`/`max_edge` (bucket queues depend on them);
//! - resample every field written back, so nothing is misregistered by the jitter.
//!
//! Area weights stay the lattice's: each node still stands for its lattice cell (a perturbed
//! lattice, not a Voronoi mesh). Drainage areas sum many cells, so the local error averages out.

use glam::DVec3;
use rg_geo::{CubeFace, CubeField};
use rg_tectonics::grid::Grid;

/// Jitter `grid`'s node positions by up to ±`amount` cells in face u and v (`amount` ≤ 0.4
/// keeps the lattice's triangles from folding over). Returns each node's offset in cells.
pub(crate) fn jitter_grid(grid: &mut Grid, amount: f64, seed: u64) -> Vec<[f64; 2]> {
    let res = grid.res;
    let n_edge = (res + 1) as usize;
    let mut offsets = vec![[0.0f64; 2]; grid.n];
    for (node, offset) in offsets.iter_mut().enumerate() {
        let k = grid.node_sample[node] as usize;
        let (face, i, j) = sample_ijk(k, n_edge);
        if i == 0 || j == 0 || i == res || j == res {
            continue;
        }
        let ou = (unit_hash(seed ^ 0x0011_7E12, k as u64) - 0.5) * 2.0 * amount;
        let ov = (unit_hash(seed ^ 0x0022_7E12, k as u64) - 0.5) * 2.0 * amount;
        *offset = [ou, ov];
        grid.dir[node] = face.uv_to_dir(cell_to_uv(i as f64 + ou, res), cell_to_uv(j as f64 + ov, res));
    }
    // Edge-length range, radians of arc, for the bucket queues.
    let (mut lo, mut hi) = (f32::INFINITY, 0f32);
    for a in 0..grid.n {
        for &b in grid.neighbors(a) {
            let ang = grid.dir[a].angle_between(grid.dir[b as usize]) as f32;
            lo = lo.min(ang);
            hi = hi.max(ang);
        }
    }
    grid.min_edge = lo;
    grid.max_edge = hi;
    offsets
}

/// Resample per-node values (on jittered positions) onto the regular `CubeField` lattice:
/// each lattice sample takes the linear interpolant of the jittered triangle that contains
/// it. Face-edge samples are unjittered nodes and copy their value.
pub(crate) fn resample(grid: &Grid, offsets: &[[f64; 2]], vals: &[f64]) -> CubeField<f32> {
    let res = grid.res;
    let mut out = CubeField::<f32>::new(res);
    let n_edge = (res + 1) as usize;
    for k in 0..out.data.len() {
        let node = grid.sample_node[k] as usize;
        let (face, i, j) = sample_ijk(k, n_edge);
        out.data[k] = if i == 0 || j == 0 || i == res || j == res {
            vals[node] as f32
        } else {
            interpolate_at(grid, offsets, vals, face, i, j) as f32
        };
    }
    out
}

/// The value at lattice point (i, j) of `face` from the jittered triangles around it.
fn interpolate_at(grid: &Grid, offsets: &[[f64; 2]], vals: &[f64], face: CubeFace, i: u32, j: u32) -> f64 {
    let n_edge = (grid.res + 1) as usize;
    let node_at = |ci: u32, cj: u32| grid.sample_node[face as usize * n_edge * n_edge + cj as usize * n_edge + ci as usize] as usize;
    // Position of the node at lattice (ci, cj), in cells, after jitter.
    let pos = |ci: u32, cj: u32| {
        let o = offsets[node_at(ci, cj)];
        (ci as f64 + o[0], cj as f64 + o[1])
    };
    let p = (i as f64, j as f64);
    let mut best: Option<(f64, f64)> = None; // (worst barycentric weight, value)
    // The four lattice cells touching (i, j), two triangles each.
    for (ci, cj) in [(i - 1, j - 1), (i, j - 1), (i - 1, j), (i, j)] {
        let corners = [(ci, cj), (ci + 1, cj), (ci, cj + 1), (ci + 1, cj + 1)];
        for tri in [[0usize, 1, 3], [0, 3, 2]] {
            let [a, b, c] = tri.map(|t| corners[t]);
            let (pa, pb, pc) = (pos(a.0, a.1), pos(b.0, b.1), pos(c.0, c.1));
            let Some(w) = barycentric(p, pa, pb, pc) else { continue };
            let v = w[0] * vals[node_at(a.0, a.1)] + w[1] * vals[node_at(b.0, b.1)] + w[2] * vals[node_at(c.0, c.1)];
            let worst = w[0].min(w[1]).min(w[2]);
            if worst >= -1e-9 {
                return v;
            }
            // Numerical fallback: the triangle the point is least outside of.
            if best.is_none_or(|(bw, _)| worst > bw) {
                best = Some((worst, v));
            }
        }
    }
    best.map_or(vals[node_at(i, j)], |(_, v)| v)
}

fn barycentric(p: (f64, f64), a: (f64, f64), b: (f64, f64), c: (f64, f64)) -> Option<[f64; 3]> {
    let det = (b.1 - c.1) * (a.0 - c.0) + (c.0 - b.0) * (a.1 - c.1);
    if det.abs() < 1e-12 {
        return None;
    }
    let wa = ((b.1 - c.1) * (p.0 - c.0) + (c.0 - b.0) * (p.1 - c.1)) / det;
    let wb = ((c.1 - a.1) * (p.0 - c.0) + (a.0 - c.0) * (p.1 - c.1)) / det;
    Some([wa, wb, 1.0 - wa - wb])
}

/// (face, i, j) of `CubeField` sample index `k` (face-major, then row j, then column i).
fn sample_ijk(k: usize, n_edge: usize) -> (CubeFace, u32, u32) {
    let face = CubeFace::ALL[k / (n_edge * n_edge)];
    let r = k % (n_edge * n_edge);
    (face, (r % n_edge) as u32, (r / n_edge) as u32)
}

fn cell_to_uv(x: f64, res: u32) -> f64 {
    -1.0 + 2.0 * x / res as f64
}

fn unit_hash(seed: u64, i: u64) -> f64 {
    let mut z = seed ^ i.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^= z >> 31;
    (z >> 11) as f64 / (1u64 << 53) as f64
}

/// A unit direction's position as a stable key for position-keyed noise (`DVec3` → lattice
/// coordinates at `wavelength_m`); used by the roughness below.
pub(crate) fn roughness_at(dir: DVec3, radius_m: f64, wavelength_m: f64, min_wavelength_m: f64, salt: u32) -> f64 {
    // fBm of 3D gradient noise from `wavelength_m` down to `min_wavelength_m` (gain 0.5).
    let (mut sum, mut norm, mut wl, mut amp) = (0.0, 0.0, wavelength_m, 1.0);
    let mut octave = 0u32;
    while wl >= min_wavelength_m && octave < 12 {
        sum += amp * rg_noise::gradient3(dir * (radius_m / wl), salt.wrapping_add(octave));
        norm += amp;
        wl *= 0.5;
        amp *= 0.5;
        octave += 1;
    }
    if norm == 0.0 { 0.0 } else { (2.0 * sum / norm).clamp(-1.0, 1.0) }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resampling_a_smooth_field_is_nearly_exact() {
        let mut g = Grid::new(48);
        let f = |d: DVec3| d.x * 1000.0 + d.y * d.z * 500.0;
        let offs = jitter_grid(&mut g, 0.35, 7);
        // Values at the jittered positions, resampled back to the lattice.
        let vals: Vec<f64> = g.dir.iter().map(|&d| f(d)).collect();
        let out = resample(&g, &offs, &vals);
        let mut worst: f64 = 0.0;
        for (k, _, _, _, d) in out.iter_samples() {
            worst = worst.max((out.data[k] as f64 - f(d)).abs());
        }
        // Linear interpolation over ~0.03 rad cells of a field with |∇| ≤ ~1500: tiny error.
        assert!(worst < 5.0, "worst resample error {worst} m");
    }

    #[test]
    fn jitter_moves_interior_nodes_only_and_updates_edges() {
        let reference = Grid::new(32);
        let mut g = Grid::new(32);
        let offs = jitter_grid(&mut g, 0.3, 1);
        let moved = (0..g.n).filter(|&i| g.dir[i] != reference.dir[i]).count();
        assert!(moved > g.n / 2);
        for (i, o) in offs.iter().enumerate() {
            assert!(o[0].abs() <= 0.3 && o[1].abs() <= 0.3);
            if *o == [0.0, 0.0] {
                assert_eq!(g.dir[i], reference.dir[i]);
            }
        }
        assert!(g.min_edge < reference.min_edge && g.max_edge > reference.max_edge);
    }
}
