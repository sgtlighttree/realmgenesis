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

    pub fn min_max(&self) -> (f32, f32) {
        self.data.iter().fold((f32::MAX, f32::MIN), |(lo, hi), &v| (lo.min(v), hi.max(v)))
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
    fn nearest_picks_exact_sample() {
        let f = CubeField::<u16>::from_fn(8, |d| if d.x > 0.0 { 1 } else { 2 });
        assert_eq!(f.nearest(DVec3::X), 1);
        assert_eq!(f.nearest(DVec3::NEG_X), 2);
    }
}
