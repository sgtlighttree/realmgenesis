//! Lattice cells on the cube-sphere with one identity across face seams.
//!
//! Procedural detail that places one feature per cell (a jittered point, a channel head)
//! needs to find the cells around a point, including cells on the neighbouring face. A cell
//! is (face, i, j) at `n` cells per face edge. Stepping past a face edge goes through the
//! cell centre's direction (`uv_to_dir` extrapolates smoothly for |u| < 2), so the cell found
//! is the one that owns that place, and no feature is ever duplicated or missing at a seam.

use crate::CubeFace;
use glam::DVec3;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Cell {
    pub face: CubeFace,
    pub i: u32,
    pub j: u32,
}

/// The cell containing `dir`, at `n` cells per face edge.
pub fn cell_of(dir: DVec3, n: u32) -> Cell {
    let (face, u, v) = CubeFace::from_dir(dir);
    let idx = |t: f64| (((t + 1.0) * 0.5 * n as f64).floor().max(0.0) as u32).min(n - 1);
    Cell { face, i: idx(u), j: idx(v) }
}

/// Face (u, v) of the centre of cell (i, j); `i`, `j` may lie outside the face.
fn centre_uv(i: i64, j: i64, n: u32) -> (f64, f64) {
    let s = 2.0 / n as f64;
    (-1.0 + (i as f64 + 0.5) * s, -1.0 + (j as f64 + 0.5) * s)
}

/// Unit direction of a cell centre.
pub fn centre(c: Cell, n: u32) -> DVec3 {
    let (u, v) = centre_uv(c.i as i64, c.j as i64, n);
    c.face.uv_to_dir(u, v)
}

/// The cell `(di, dj)` steps from `c`, resolved through the owning face past a seam.
/// Steps must stay within half a face (|di|, |dj| < n/2).
pub fn offset(c: Cell, di: i32, dj: i32, n: u32) -> Cell {
    let (i, j) = (c.i as i64 + di as i64, c.j as i64 + dj as i64);
    if (0..n as i64).contains(&i) && (0..n as i64).contains(&j) {
        return Cell { face: c.face, i: i as u32, j: j as u32 };
    }
    let (u, v) = centre_uv(i, j, n);
    cell_of(c.face.uv_to_dir(u, v), n)
}

/// The cells in a (2r + 1)² window around the cell containing `dir`, without duplicates
/// (near a cube corner two window positions can resolve to the same cell). Sorted, so the
/// result is deterministic.
pub fn window(dir: DVec3, n: u32, r: i32) -> Vec<Cell> {
    let c = cell_of(dir, n);
    let mut out = Vec::with_capacity(((2 * r + 1) * (2 * r + 1)) as usize);
    for dj in -r..=r {
        for di in -r..=r {
            out.push(offset(c, di, dj, n));
        }
    }
    out.sort_unstable();
    out.dedup();
    out
}

/// A jittered point inside cell `c`: the centre moved by up to `jitter` (0..0.5) cell widths
/// in u and v, from a hash of (salt, cell). Deterministic.
pub fn jittered(c: Cell, n: u32, salt: u64, jitter: f64) -> DVec3 {
    let h = hash(salt, c);
    let ru = ((h & 0xFFFF_FFFF) as f64 / 4_294_967_296.0) - 0.5;
    let rv = ((h >> 32) as f64 / 4_294_967_296.0) - 0.5;
    let (u, v) = centre_uv(c.i as i64, c.j as i64, n);
    let s = 2.0 / n as f64;
    c.face.uv_to_dir(u + 2.0 * jitter * ru * s, v + 2.0 * jitter * rv * s)
}

/// 64-bit hash of (salt, cell) (splitmix64 finaliser over a packed key).
pub fn hash(salt: u64, c: Cell) -> u64 {
    let key = (c.face as u64) << 58 ^ (c.i as u64) << 29 ^ c.j as u64;
    let mut z = salt ^ key.wrapping_mul(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn dirs(count: usize) -> Vec<DVec3> {
        // Fibonacci sphere, plus points hugging every face edge and corner.
        let golden = std::f64::consts::PI * (3.0 - 5f64.sqrt());
        let mut v: Vec<DVec3> = (0..count)
            .map(|k| {
                let z = 1.0 - 2.0 * (k as f64 + 0.5) / count as f64;
                let r = (1.0 - z * z).sqrt();
                DVec3::new(r * (golden * k as f64).cos(), r * (golden * k as f64).sin(), z)
            })
            .collect();
        for face in CubeFace::ALL {
            for t in [-0.999, -0.5, 0.0, 0.5, 0.999] {
                v.push(face.uv_to_dir(0.9999, t));
                v.push(face.uv_to_dir(t, -0.9999));
            }
        }
        v
    }

    /// The nearest jittered point found through a window equals the brute-force nearest
    /// over every cell on the sphere: seams neither hide nor duplicate points.
    #[test]
    fn window_search_matches_brute_force_across_seams() {
        let n = 12;
        let all: Vec<Cell> = CubeFace::ALL
            .into_iter()
            .flat_map(|face| (0..n).flat_map(move |j| (0..n).map(move |i| Cell { face, i, j })))
            .collect();
        for d in dirs(400) {
            let brute = all.iter().copied().min_by(|a, b| {
                let da = jittered(*a, n, 7, 0.4).angle_between(d);
                let db = jittered(*b, n, 7, 0.4).angle_between(d);
                da.total_cmp(&db).then(a.cmp(b))
            });
            let win = window(d, n, 2).into_iter().min_by(|a, b| {
                let da = jittered(*a, n, 7, 0.4).angle_between(d);
                let db = jittered(*b, n, 7, 0.4).angle_between(d);
                da.total_cmp(&db).then(a.cmp(b))
            });
            assert_eq!(brute, win, "at {d:?}");
        }
    }

    #[test]
    fn offset_round_trips_away_from_corners() {
        let n = 16;
        for face in CubeFace::ALL {
            for j in 2..n - 2 {
                // Step off the +u edge and back.
                let c = Cell { face, i: n - 1, j };
                let out = offset(c, 1, 0, n);
                assert_ne!(out.face, face);
                let back: Vec<Cell> = window(centre(out, n), n, 1);
                assert!(back.contains(&c), "{c:?} → {out:?} doesn't see it back");
            }
        }
    }

    #[test]
    fn cell_of_centre_is_itself() {
        let n = 9;
        for face in CubeFace::ALL {
            for j in 0..n {
                for i in 0..n {
                    let c = Cell { face, i, j };
                    assert_eq!(cell_of(centre(c, n), n), c);
                }
            }
        }
    }
}
