//! Whole-planet scalar fields on the cube-sphere: the shared currency between sim crates.
//!
//! A `CubeField` stores `(res + 1)²` corner-aligned samples per face (same convention as
//! tiles), so face edges are sampled on both sides with identical directions.

use crate::CubeFace;
use glam::DVec3;

#[derive(Clone, Debug, PartialEq)]
pub struct CubeField<T> {
    /// Cells per face edge. Samples per edge = `res + 1`.
    pub res: u32,
    /// Face-major, then row (v) major, then column (u): `data[face][j][i]`.
    pub data: Vec<T>,
}

impl<T: Copy + Default> CubeField<T> {
    pub fn new(res: u32) -> Self {
        let n = (res + 1) as usize;
        Self { res, data: vec![T::default(); 6 * n * n] }
    }

    /// Build by evaluating `f(dir)` at every sample.
    pub fn from_fn(res: u32, mut f: impl FnMut(DVec3) -> T) -> Self {
        let mut out = Self::new(res);
        for face in CubeFace::ALL {
            for j in 0..=res {
                for i in 0..=res {
                    let d = out.sample_dir(face, i, j);
                    let k = out.index(face, i, j);
                    out.data[k] = f(d);
                }
            }
        }
        out
    }
}

impl<T> CubeField<T> {
    pub fn samples_per_edge(&self) -> u32 {
        self.res + 1
    }

    #[inline]
    pub fn index(&self, face: CubeFace, i: u32, j: u32) -> usize {
        let n = (self.res + 1) as usize;
        face as usize * n * n + j as usize * n + i as usize
    }

    /// Unit direction of sample (face, i, j).
    pub fn sample_dir(&self, face: CubeFace, i: u32, j: u32) -> DVec3 {
        let u = -1.0 + 2.0 * i as f64 / self.res as f64;
        let v = -1.0 + 2.0 * j as f64 / self.res as f64;
        face.uv_to_dir(u, v)
    }

    /// Iterate (index, face, i, j, dir) over every sample.
    pub fn iter_samples(&self) -> impl Iterator<Item = (usize, CubeFace, u32, u32, DVec3)> + '_ {
        let r = self.res;
        CubeFace::ALL.into_iter().flat_map(move |face| {
            (0..=r).flat_map(move |j| {
                (0..=r).map(move |i| (self.index(face, i, j), face, i, j, self.sample_dir(face, i, j)))
            })
        })
    }

    /// Approximate spacing between neighbouring samples, in radians of arc.
    pub fn spacing_rad(&self) -> f64 {
        std::f64::consts::FRAC_PI_2 / self.res as f64
    }
}

impl CubeField<f32> {
    /// Bilinear sample at a unit direction.
    pub fn sample(&self, dir: DVec3) -> f32 {
        let (face, u, v) = CubeFace::from_dir(dir);
        let r = self.res as f64;
        let x = ((u + 1.0) * 0.5 * r).clamp(0.0, r);
        let y = ((v + 1.0) * 0.5 * r).clamp(0.0, r);
        let (i0, j0) = ((x.floor() as u32).min(self.res - 1), (y.floor() as u32).min(self.res - 1));
        let (fx, fy) = ((x - i0 as f64) as f32, (y - j0 as f64) as f32);
        let a = self.data[self.index(face, i0, j0)];
        let b = self.data[self.index(face, i0 + 1, j0)];
        let c = self.data[self.index(face, i0, j0 + 1)];
        let d = self.data[self.index(face, i0 + 1, j0 + 1)];
        let top = a + (b - a) * fx;
        let bot = c + (d - c) * fx;
        top + (bot - top) * fy
    }

    /// Bilinear sample at a unit direction, interpolated in f64. Use this when the result
    /// is differentiated: f32 interpolation of values ~1e4 adds ~1e-3 rounding noise,
    /// which finite differences turn into slope noise.
    pub fn sample_f64(&self, dir: DVec3) -> f64 {
        let (face, u, v) = CubeFace::from_dir(dir);
        let r = self.res as f64;
        let x = ((u + 1.0) * 0.5 * r).clamp(0.0, r);
        let y = ((v + 1.0) * 0.5 * r).clamp(0.0, r);
        let (i0, j0) = ((x.floor() as u32).min(self.res - 1), (y.floor() as u32).min(self.res - 1));
        let (fx, fy) = (x - i0 as f64, y - j0 as f64);
        let a = self.data[self.index(face, i0, j0)] as f64;
        let b = self.data[self.index(face, i0 + 1, j0)] as f64;
        let c = self.data[self.index(face, i0, j0 + 1)] as f64;
        let d = self.data[self.index(face, i0 + 1, j0 + 1)] as f64;
        let top = a + (b - a) * fx;
        let bot = c + (d - c) * fx;
        top + (bot - top) * fy
    }

    /// Drained sample, in f64: bilinear, with any diagonal valley carved straight.
    /// Use this for drained heights (macro terrain).
    ///
    /// Why: on a surface routed with D8 (8 neighbours), a node can drain diagonally across a
    /// cell. Bilinear interpolation averages in the two other corners there, so a diagonal
    /// valley rises mid-cell and dams itself (measured: up to ~10% of mountain land in
    /// closed depressions, 150 m deep). If node d drains diagonally to a, then a is lower
    /// than both other corners (otherwise the shorter, steeper axis neighbour would have
    /// been d's receiver), so the drainage diagonal is the one through the lowest corner.
    /// See [`drained_interp`] for the carve. (A plain triangle split is also drained, but
    /// its flat facets and creases showed in the hillshade, s01.)
    pub fn sample_drained(&self, dir: DVec3) -> f64 {
        let (face, u, v) = CubeFace::from_dir(dir);
        let r = self.res as f64;
        let x = ((u + 1.0) * 0.5 * r).clamp(0.0, r);
        let y = ((v + 1.0) * 0.5 * r).clamp(0.0, r);
        let (i0, j0) = ((x.floor() as u32).min(self.res - 1), (y.floor() as u32).min(self.res - 1));
        let (fx, fy) = (x - i0 as f64, y - j0 as f64);
        let a = self.data[self.index(face, i0, j0)] as f64;
        let b = self.data[self.index(face, i0 + 1, j0)] as f64;
        let c = self.data[self.index(face, i0, j0 + 1)] as f64;
        let d = self.data[self.index(face, i0 + 1, j0 + 1)] as f64;
        drained_interp(a, b, c, d, fx, fy)
    }

    /// Smooth (C¹) sample for display: Catmull-Rom bicubic over the 4×4 neighbourhood. It
    /// passes through every sample. Indices are clamped at face edges (one-sided there, so a
    /// faint crease can remain along the 12 cube edges).
    ///
    /// Not drained (valleys can rise slightly mid-cell, and it can overshoot a little), so
    /// never use it for routing or as the detail network's cap. It's for the macro-only
    /// display, where bilinear facets and creases showed at high relief (s02).
    pub fn sample_smooth(&self, dir: DVec3) -> f64 {
        let (face, u, v) = CubeFace::from_dir(dir);
        let r = self.res as f64;
        let x = ((u + 1.0) * 0.5 * r).clamp(0.0, r);
        let y = ((v + 1.0) * 0.5 * r).clamp(0.0, r);
        let last = self.res as i64;
        let (i0, j0) = ((x.floor() as i64).min(last - 1), (y.floor() as i64).min(last - 1));
        let (fx, fy) = (x - i0 as f64, y - j0 as f64);
        let at = |i: i64, j: i64| self.data[self.index(face, i.clamp(0, last) as u32, j.clamp(0, last) as u32)] as f64;
        // Catmull-Rom through p1 (t = 0) and p2 (t = 1), with p0 and p3 as the outer points.
        let cr = |p0: f64, p1: f64, p2: f64, p3: f64, t: f64| {
            p1 + 0.5 * t * (p2 - p0 + t * (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3 + t * (3.0 * (p1 - p2) + p3 - p0)))
        };
        let row = |j: i64| cr(at(i0 - 1, j), at(i0, j), at(i0 + 1, j), at(i0 + 2, j), fx);
        cr(row(j0 - 1), row(j0), row(j0 + 1), row(j0 + 2), fy)
    }

    /// Piecewise-linear sample on two triangles per cell, split along the diagonal through
    /// the cell's lowest corner: provably drained for a pit-free field (every axis edge is
    /// a triangle edge, and a lower diagonal neighbour is the cell's unique minimum, so its
    /// edge is kept). Faceted, so it's used as the network detail's cap, not for display.
    pub fn sample_tri(&self, dir: DVec3) -> f64 {
        let (face, u, v) = CubeFace::from_dir(dir);
        let r = self.res as f64;
        let x = ((u + 1.0) * 0.5 * r).clamp(0.0, r);
        let y = ((v + 1.0) * 0.5 * r).clamp(0.0, r);
        let (i0, j0) = ((x.floor() as u32).min(self.res - 1), (y.floor() as u32).min(self.res - 1));
        let (fx, fy) = (x - i0 as f64, y - j0 as f64);
        let a = self.data[self.index(face, i0, j0)] as f64;
        let b = self.data[self.index(face, i0 + 1, j0)] as f64;
        let c = self.data[self.index(face, i0, j0 + 1)] as f64;
        let d = self.data[self.index(face, i0 + 1, j0 + 1)] as f64;
        tri_interp(a, b, c, d, fx, fy)
    }

    pub fn min_max(&self) -> (f32, f32) {
        self.data.iter().fold((f32::MAX, f32::MIN), |(lo, hi), &v| (lo.min(v), hi.max(v)))
    }
}

/// Two triangles per cell, split along the diagonal through the lowest corner (ties: a–d).
#[inline]
pub fn tri_interp(a: f64, b: f64, c: f64, d: f64, fx: f64, fy: f64) -> f64 {
    if a.min(d) <= b.min(c) {
        if fx >= fy { a + (b - a) * fx + (d - b) * fy } else { a + (c - a) * fy + (d - c) * fx }
    } else if fx + fy <= 1.0 {
        a + (b - a) * fx + (c - a) * fy
    } else {
        d + (c - d) * (1.0 - fx) + (b - d) * (1.0 - fy)
    }
}

/// Identifies [`drained_interp`] in world files. Bump if the interpolant changes: the
/// committed macro heights only mean the same surface under the same interpolant.
pub const MACRO_INTERP_ID: &str = "drained-bilinear-v1";

/// Interpolate a cell with corners a = (0,0), b = (1,0), c = (0,1), d = (1,1).
///
/// Bilinear is `plane(fx, fy) − e·fx·fy` with `e = b + c − a − d`. Along the a–d diagonal
/// (fx = fy = t) it exceeds the straight line from a to d by `e·t(1 − t)`: when `e > 0` (b and
/// c are ridges) that bump is the dam. We subtract
///
/// `e · fx·fy·(1 − fx)(1 − fy) / (t(1 − t))`, with `t = (fx + fy)/2`,
///
/// which equals the bump on the diagonal (so the valley is exactly linear there), is zero on
/// all four cell edges (so neighbouring cells still agree), and is smooth in between (a
/// rounded valley, not a crease). The diagonal through the lowest corner is the one carved
/// (mirrored for b–c). With `e ≤ 0` the diagonal is a ridge and bilinear is left alone.
#[inline]
pub fn drained_interp(a: f64, b: f64, c: f64, d: f64, fx: f64, fy: f64) -> f64 {
    let bilinear = {
        let top = a + (b - a) * fx;
        let bot = c + (d - c) * fx;
        top + (bot - top) * fy
    };
    // Mirror so the carved diagonal always runs from (0,0) to (1,1) in (x, fy).
    let (e, x) = if a.min(d) <= b.min(c) { (b + c - a - d, fx) } else { (a + d - b - c, 1.0 - fx) };
    if e <= 0.0 {
        return bilinear;
    }
    let t = 0.5 * (x + fy);
    let tt = t * (1.0 - t);
    if tt <= 1e-12 {
        return bilinear; // at the diagonal's end corners the carve is 0
    }
    bilinear - e * x * fy * (1.0 - x) * (1.0 - fy) / tt
}

impl<T: Copy> CubeField<T> {
    /// Nearest-sample lookup at a unit direction (for categorical fields like plate ids).
    pub fn nearest(&self, dir: DVec3) -> T {
        let (face, u, v) = CubeFace::from_dir(dir);
        let r = self.res as f64;
        let i = ((u + 1.0) * 0.5 * r).round().clamp(0.0, r) as u32;
        let j = ((v + 1.0) * 0.5 * r).round().clamp(0.0, r) as u32;
        self.data[self.index(face, i, j)]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sample_reproduces_smooth_function() {
        let f = CubeField::<f32>::from_fn(32, |d| d.z as f32);
        for &d in &[DVec3::new(0.2, 0.5, 0.8), DVec3::new(-0.9, 0.1, -0.3), DVec3::Z] {
            let d = d.normalize();
            assert!((f.sample(d) - d.z as f32).abs() < 0.01);
        }
    }

    #[test]
    fn smooth_sample_passes_through_samples_and_tracks_smooth_fields() {
        let f = CubeField::<f32>::from_fn(32, |d| (d.x * 3000.0 + d.y * d.z * 2000.0) as f32);
        for (k, face, i, j, d) in f.iter_samples() {
            if face == CubeFace::ALL[2] && (2..30).contains(&i) && (2..30).contains(&j) {
                assert!((f.sample_smooth(d) - f.data[k] as f64).abs() < 1e-2, "misses sample ({i}, {j})");
            }
        }
        let g = CubeField::<f32>::from_fn(32, |d| (d.z * 1000.0) as f32);
        let d = DVec3::new(0.21, 0.13, 0.97).normalize();
        assert!((g.sample_smooth(d) - d.z * 1000.0).abs() < 1.0);
    }

    #[test]
    fn drained_interp_hits_corners_and_is_continuous() {
        for &(a, b, c, d) in &[(0.0, 5.0, 7.0, 1.0), (4.0, 0.0, 1.0, 9.0), (3.0, 3.0, 3.0, 3.0), (2.0, 9.0, 0.0, 8.0)] {
            assert_eq!(drained_interp(a, b, c, d, 0.0, 0.0), a);
            assert_eq!(drained_interp(a, b, c, d, 1.0, 0.0), b);
            assert_eq!(drained_interp(a, b, c, d, 0.0, 1.0), c);
            assert_eq!(drained_interp(a, b, c, d, 1.0, 1.0), d);
            // Linear along each edge (so neighbouring cells agree).
            for k in 0..=20 {
                let t = k as f64 / 20.0;
                let lin = |p: f64, q: f64| p + (q - p) * t;
                assert!((drained_interp(a, b, c, d, t, 0.0) - lin(a, b)).abs() < 1e-9);
                assert!((drained_interp(a, b, c, d, 1.0, t) - lin(b, d)).abs() < 1e-9);
                assert!((drained_interp(a, b, c, d, t, 1.0) - lin(c, d)).abs() < 1e-9);
                assert!((drained_interp(a, b, c, d, 0.0, t) - lin(a, c)).abs() < 1e-9);
            }
            // No jumps inside the cell.
            for i in 0..50 {
                for j in 0..50 {
                    let (x, y) = (i as f64 / 50.0 + 0.01, j as f64 / 50.0 + 0.01);
                    let h = drained_interp(a, b, c, d, x, y);
                    let hx = drained_interp(a, b, c, d, x + 1e-7, y);
                    let hy = drained_interp(a, b, c, d, x, y + 1e-7);
                    assert!((h - hx).abs() < 1e-4 && (h - hy).abs() < 1e-4, "jump at {x},{y}");
                }
            }
        }
    }

    /// A valley draining diagonally (d → a, with b and c as ridges) is exactly linear along
    /// the diagonal; bilinear interpolation dams it mid-cell. Same for the mirrored b–c case.
    #[test]
    fn diagonal_valley_is_not_dammed() {
        let (a, b, c, d) = (0.0, 100.0, 120.0, 10.0);
        for k in 0..=100 {
            let t = k as f64 / 100.0;
            assert!((drained_interp(a, b, c, d, t, t) - 10.0 * t).abs() < 1e-9);
        }
        assert!((a + b + c + d) / 4.0 > d, "bilinear would dam this valley");
        let (a, b, c, d) = (100.0, 0.0, 10.0, 120.0);
        for k in 0..=100 {
            let t = k as f64 / 100.0;
            assert!((drained_interp(a, b, c, d, 1.0 - t, t) - 10.0 * t).abs() < 1e-9);
        }
    }

    #[test]
    fn nearest_picks_exact_sample() {
        let f = CubeField::<u16>::from_fn(8, |d| if d.x > 0.0 { 1 } else { 2 });
        assert_eq!(f.nearest(DVec3::X), 1);
        assert_eq!(f.nearest(DVec3::NEG_X), 2);
    }
}
