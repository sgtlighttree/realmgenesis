//! Planet shape and grids: the ellipsoid, the cube-sphere, and tile addresses.
//!
//! Conventions (see docs/ARCHITECTURE.md §4):
//! - Latitudes are **geodetic**, in radians unless a name ends in `_deg`.
//! - Longitude is east-positive, in -π..π.
//! - ECEF: +Z is the north pole, +X is (lat 0, lon 0), +Y is (lat 0, lon 90°E). Metres, f64.

pub mod cube;
pub mod field;
pub mod ellipsoid;

pub use cube::{CubeFace, TileId};
pub use field::CubeField;
pub use ellipsoid::Ellipsoid;
