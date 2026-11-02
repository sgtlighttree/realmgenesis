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

    /// Piecewise-linear sample on triangles, in f64: each cell is split along the diagonal
    /// through its **lowest corner**. Use this for drained heights (macro terrain).
    ///
    /// Why: on a surface routed with D8 (8 neighbours), a node can drain diagonally across a
    /// cell. Bilinear interpolation averages in the two other corners there, so a diagonal
    /// valley rises mid-cell and dams itself (measured: up to ~10% of mountain land in
    /// closed depressions, 150 m deep). If node d drains diagonally to a, then a is lower
    /// than both other corners (otherwise the shorter, steeper axis neighbour would have
    /// been d's receiver), so a is the cell's lowest corner and the split keeps a–d as an
    /// edge. Every drainage edge is then a triangle edge, and a piecewise-linear surface has
    /// minima only at vertices, so a drained (filled) field interpolates without depressions.
    /// Continuous everywhere (linear along every cell edge, like bilinear).
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
        tri_interp(a, b, c, d, fx, fy)
    }

    pub fn min_max(&self) -> (f32, f32) {
        self.data.iter().fold((f32::MAX, f32::MIN), |(lo, hi), &v| (lo.min(v), hi.max(v)))
    }
}

/// Interpolate a cell with corners a = (0,0), b = (1,0), c = (0,1), d = (1,1), split into
/// two triangles along the diagonal through the lowest corner (see
/// [`CubeField::sample_drained`]). Ties go to the a–d diagonal.
#[inline]
pub fn tri_interp(a: f64, b: f64, c: f64, d: f64, fx: f64, fy: f64) -> f64 {
    if a.min(d) <= b.min(c) {
        // Diagonal a–d: triangles (a, b, d) below it, (a, c, d) above it.
        if fx >= fy { a + (b - a) * fx + (d - b) * fy } else { a + (c - a) * fy + (d - c) * fx }
    } else {
        // Diagonal b–c: triangles (a, b, c) and (b, c, d).
        if fx + fy <= 1.0 { a + (b - a) * fx + (c - a) * fy } else { d + (c - d) * (1.0 - fx) + (b - d) * (1.0 - fy) }
    }
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
    fn tri_interp_hits_corners_and_is_continuous() {
        for &(a, b, c, d) in &[(0.0, 5.0, 7.0, 1.0), (4.0, 0.0, 1.0, 9.0), (3.0, 3.0, 3.0, 3.0)] {
            assert_eq!(tri_interp(a, b, c, d, 0.0, 0.0), a);
            assert_eq!(tri_interp(a, b, c, d, 1.0, 0.0), b);
            assert_eq!(tri_interp(a, b, c, d, 0.0, 1.0), c);
            assert_eq!(tri_interp(a, b, c, d, 1.0, 1.0), d);
            // Linear along each edge (so neighbouring cells agree).
            let mid = |p: f64, q: f64| (p + q) * 0.5;
            assert!((tri_interp(a, b, c, d, 0.5, 0.0) - mid(a, b)).abs() < 1e-12);
            assert!((tri_interp(a, b, c, d, 1.0, 0.5) - mid(b, d)).abs() < 1e-12);
            assert!((tri_interp(a, b, c, d, 0.5, 1.0) - mid(c, d)).abs() < 1e-12);
            assert!((tri_interp(a, b, c, d, 0.0, 0.5) - mid(a, c)).abs() < 1e-12);
            // No jumps across the diagonals inside the cell.
            for k in 0..=100 {
                let t = k as f64 / 100.0;
                let e = 1e-9;
                assert!((tri_interp(a, b, c, d, t + e, t) - tri_interp(a, b, c, d, t, t + e)).abs() < 1e-6);
                assert!((tri_interp(a, b, c, d, t + e, 1.0 - t) - tri_interp(a, b, c, d, t - e, 1.0 - t)).abs() < 1e-6);
            }
        }
    }

    /// A valley draining diagonally (d → a, with b and c as ridges) stays monotone along the
    /// diagonal; bilinear interpolation dams it mid-cell.
    #[test]
    fn diagonal_valley_is_not_dammed() {
        let (a, b, c, d) = (0.0, 100.0, 100.0, 10.0);
        let mut prev = f64::MAX;
        for k in (0..=100).rev() {
            let t = k as f64 / 100.0;
            let h = tri_interp(a, b, c, d, t, t);
            assert!(h <= prev, "rises at t = {t}");
            prev = h;
        }
        let bilinear_mid = (a + b + c + d) / 4.0;
        assert!(bilinear_mid > d, "bilinear would dam this valley");
    }

    #[test]
    fn nearest_picks_exact_sample() {
        let f = CubeField::<u16>::from_fn(8, |d| if d.x > 0.0 { 1 } else { 2 });
        assert_eq!(f.nearest(DVec3::X), 1);
        assert_eq!(f.nearest(DVec3::NEG_X), 2);
    }
}
