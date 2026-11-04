//! Step 2: plate and boundary types. Since R3 the plates are built in `fit.rs`, fitted to the
//! coast classes (`coast.rs`): kinds are prescribed by the classification, and the plates'
//! rotations are solved so the rates match.
//!
//! * **Motion.** Each plate rotates about an Euler pole: the velocity at point `p` is
//!   `omega * (pole x p)`.
//! * **Boundary normal.** At each boundary node the normal comes from the local plate layout
//!   (the mean positions of the two plates' nodes within 3 hops); the relative velocity splits
//!   into a part across the boundary (convergent / divergent) and a part along it (transform).

use crate::Plate;
use glam::DVec3;

pub const STREAM_PLATES: u32 = 20;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Kind {
    /// Continent meets continent: broad high mountains and plateau.
    Collision,
    /// Oceanic plate dives under another plate: trench plus arc or cordillera.
    Subduction,
    /// Oceanic spreading ridge.
    Ridge,
    /// Continental rift valley.
    Rift,
    /// Plates slide past each other.
    Transform,
}

#[derive(Clone, Copy, Debug)]
pub struct BoundaryNode {
    pub node: u32,
    pub kind: Kind,
    /// Relative strength of the tectonic activity here (about 0.3 to 1.4), combining plate
    /// speed with a slowly varying along-the-boundary factor.
    pub strength: f32,
    /// Full relative speed of the two plates here, mm/yr (smooth along a boundary).
    pub speed_mm_yr: f32,
    /// For subduction: the plate that dives.
    pub sub_plate: u16,
    /// For subduction: is the overriding plate continental there?
    pub over_cont: bool,
}

pub struct PlateModel {
    pub plate_of: Vec<u16>,
    pub plates: Vec<Plate>,
    /// 0 none, 1 convergent, 2 divergent, 3 transform; set on both sides of every border.
    pub boundary: Vec<u8>,
    pub bnodes: Vec<BoundaryNode>,
    /// Clamped convergence rate at convergent boundary nodes, mm/yr (0 elsewhere). R3: the
    /// rate the relief was built from, committed so erosion's uplift reads the same one.
    pub convergence: Vec<f32>,
    /// Coast class per node near coasts (`coast::margin`).
    pub margin: Vec<u8>,
}

pub fn velocity(plate: &Plate, p: DVec3) -> DVec3 {
    plate.euler_pole.cross(p) * plate.omega
}
