//! `rg-climate`: monthly climate and Köppen-Geiger classes for an arbitrary planet.
//!
//! Tier A ([`tier_a`]) is the fast, rule-based preview used while rerolling a world. Tier B (a
//! physics solver, later) will reuse the grid type ([`LatLonGrid`]), the output struct
//! ([`ClimateFields`]) and the classifier ([`Koppen`]).
//!
//! Everything is deterministic, single-threaded and allocation-light so it runs in wasm32.

pub mod water;
pub mod circulation;
pub mod ebm;
pub mod geography;
pub mod grid;
pub mod insolation;
pub mod koppen;
pub mod synthetic;
pub mod tier_a;
pub mod types;

pub use grid::LatLonGrid;
pub use koppen::Koppen;
pub use tier_a::tier_a;
pub use types::{ClimateFields, ClimateInput, PlanetParams};
