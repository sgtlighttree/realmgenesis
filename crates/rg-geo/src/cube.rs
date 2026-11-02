//! Cube-sphere: six faces, each a quadtree of tiles.
//!
//! The face coordinates (u, v) lie in [-1, 1]. They're warped with `tan(u·π/4)` before
//! projecting onto the sphere, which keeps the cell-area ratio near 1.4 instead of the
//! raw gnomonic ≈ 5.2.

use glam::DVec3;
use std::f64::consts::FRAC_PI_4;

/// Cells per tile edge. Tiles store `TILE_SAMPLES + 1` corner-aligned vertices per edge
/// (257²), plus a 1-sample skirt, so parent samples coincide with every second child
/// sample. See ARCHITECTURE §4.
pub const TILE_SAMPLES: u32 = 256;
/// Vertices per tile edge.
pub const TILE_VERTS: u32 = TILE_SAMPLES + 1;
/// Face-warp id stored in world files. 1 = tangent warp `tan(u·π/4)`.
pub const FACE_WARP_ID: u8 = 1;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
#[repr(u8)]
pub enum CubeFace {
    PosX = 0,
    NegX = 1,
    PosY = 2,
    NegY = 3,
    PosZ = 4,
    NegZ = 5,
}

impl CubeFace {
    pub const ALL: [CubeFace; 6] = [
        CubeFace::PosX,
        CubeFace::NegX,
        CubeFace::PosY,
        CubeFace::NegY,
        CubeFace::PosZ,
        CubeFace::NegZ,
    ];

    /// (normal, u axis, v axis). Each frame is right-handed, so `u × v = normal` and
    /// triangles wound (u, v) face outward.
    pub fn frame(self) -> (DVec3, DVec3, DVec3) {
        match self {
            CubeFace::PosX => (DVec3::X, DVec3::Y, DVec3::Z),
            CubeFace::NegX => (DVec3::NEG_X, DVec3::NEG_Y, DVec3::Z),
            CubeFace::PosY => (DVec3::Y, DVec3::NEG_X, DVec3::Z),
            CubeFace::NegY => (DVec3::NEG_Y, DVec3::X, DVec3::Z),
            CubeFace::PosZ => (DVec3::Z, DVec3::Y, DVec3::NEG_X),
            CubeFace::NegZ => (DVec3::NEG_Z, DVec3::Y, DVec3::X),
        }
    }

    /// Face coords (u, v) in [-1, 1] → unit direction.
    pub fn uv_to_dir(self, u: f64, v: f64) -> DVec3 {
        let (n, eu, ev) = self.frame();
        let wu = (u * FRAC_PI_4).tan();
        let wv = (v * FRAC_PI_4).tan();
        (n + eu * wu + ev * wv).normalize()
    }

    /// Unit direction → (face, u, v). The inverse of `uv_to_dir`.
    pub fn from_dir(d: DVec3) -> (CubeFace, f64, f64) {
        let a = d.abs();
        let face = if a.x >= a.y && a.x >= a.z {
            if d.x > 0.0 { CubeFace::PosX } else { CubeFace::NegX }
        } else if a.y >= a.z {
            if d.y > 0.0 { CubeFace::PosY } else { CubeFace::NegY }
        } else if d.z > 0.0 {
            CubeFace::PosZ
        } else {
            CubeFace::NegZ
        };
        let (n, eu, ev) = face.frame();
        let p = d / d.dot(n);
        let u = p.dot(eu).atan() / FRAC_PI_4;
        let v = p.dot(ev).atan() / FRAC_PI_4;
        (face, u, v)
    }
}

/// A quadtree tile: level L splits each face into 2^L × 2^L tiles.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct TileId {
    pub face: CubeFace,
    pub level: u8,
    pub x: u32,
    pub y: u32,
}

impl TileId {
    pub fn root(face: CubeFace) -> Self {
        Self { face, level: 0, x: 0, y: 0 }
    }

    pub fn children(self) -> [TileId; 4] {
        let (l, x, y) = (self.level + 1, self.x * 2, self.y * 2);
        [
            TileId { face: self.face, level: l, x, y },
            TileId { face: self.face, level: l, x: x + 1, y },
            TileId { face: self.face, level: l, x, y: y + 1 },
            TileId { face: self.face, level: l, x: x + 1, y: y + 1 },
        ]
    }

    /// The tile's (u0, v0, u1, v1) bounds in face coords.
    pub fn uv_bounds(self) -> (f64, f64, f64, f64) {
        let n = (1u64 << self.level) as f64;
        let s = 2.0 / n;
        let u0 = -1.0 + self.x as f64 * s;
        let v0 = -1.0 + self.y as f64 * s;
        (u0, v0, u0 + s, v0 + s)
    }

    /// Approximate cell size in metres at this level on a sphere of `radius_m`:
    /// a face edge spans a quarter great circle.
    pub fn cell_size_m(level: u8, radius_m: f64) -> f64 {
        (std::f64::consts::FRAC_PI_2 * radius_m) / (TILE_SAMPLES as f64 * (1u64 << level) as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frames_are_right_handed() {
        for f in CubeFace::ALL {
            let (n, u, v) = f.frame();
            assert!((u.cross(v) - n).length() < 1e-12, "{f:?}");
        }
    }

    #[test]
    fn dir_roundtrip() {
        for f in CubeFace::ALL {
            for &(u, v) in &[(0.0, 0.0), (0.5, -0.25), (-0.99, 0.99), (0.3, 0.7)] {
                let d = f.uv_to_dir(u, v);
                let (f2, u2, v2) = CubeFace::from_dir(d);
                assert_eq!(f, f2);
                assert!((u - u2).abs() < 1e-12 && (v - v2).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn shared_edges_agree() {
        // Every face's u = +1 edge must land on some other face's edge: no gaps.
        for f in CubeFace::ALL {
            for i in 0..=8 {
                let v = -1.0 + i as f64 * 0.25;
                let d = f.uv_to_dir(1.0, v);
                let (_, u2, v2) = CubeFace::from_dir(d);
                assert!(u2.abs().max(v2.abs()) > 1.0 - 1e-9);
            }
        }
    }

    #[test]
    fn level_sizes_match_architecture_doc() {
        let r = 6_371_000.0;
        assert!((TileId::cell_size_m(0, r) - 39_092.0).abs() < 100.0);
        assert!((TileId::cell_size_m(4, r) - 2_443.0).abs() < 10.0);
        assert!((TileId::cell_size_m(10, r) - 38.2).abs() < 0.5);
        assert!((TileId::cell_size_m(12, r) - 9.5).abs() < 0.1);
    }

    #[test]
    fn warp_keeps_area_ratio_low() {
        // Compare the solid angle of a centre cell with a corner cell at 64 × 64.
        let f = CubeFace::PosX;
        let n = 64.0;
        let cell = |i: f64, j: f64| {
            let (u0, v0) = (-1.0 + 2.0 * i / n, -1.0 + 2.0 * j / n);
            let s = 2.0 / n;
            let a = f.uv_to_dir(u0, v0);
            let b = f.uv_to_dir(u0 + s, v0);
            let c = f.uv_to_dir(u0, v0 + s);
            (b - a).cross(c - a).length()
        };
        let ratio = cell(32.0, 32.0) / cell(0.0, 0.0);
        assert!(ratio < 1.5, "area ratio {ratio}");
    }
}
