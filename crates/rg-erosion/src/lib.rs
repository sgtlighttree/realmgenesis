//! Continental-scale fluvial erosion: the **stream-power law** solved with the implicit,
//! O(n) scheme of Braun & Willett (2013), on a drainage network routed to the sea with
//! Priority-Flood+ε (Barnes et al. 2014).
//!
//! # The physics
//!
//! Over geological time a landscape obeys
//!
//! ```text
//! ∂h/∂t = U − K · A^m · S^n
//! ```
//!
//! - `U`: rock uplift rate (m/yr), from tectonics: fast in active mountain belts, slow on cratons.
//! - `K · A^m · S^n`: the rate rivers cut down (detachment-limited incision). `A` is the
//!   upstream drainage area (m²), `S` the channel slope towards the next node downstream.
//!   Bigger rivers (more area) and steeper channels cut faster.
//! - `n = 1` (linear in slope) makes the implicit update closed-form and unconditionally
//!   stable. `m/n = θ ≈ 0.45` is the concavity seen in real river profiles (Lague 2014).
//!
//! The result is what makes terrain look real: dendritic valley networks, ridges that are
//! divides between basins, concave river profiles, mountain fronts dissected into spurs,
//! and every drop of rain finding the sea.
//!
//! # Each time step
//!
//! 1. **Receivers:** every land node drains to its steepest-descent neighbour (D8 on the
//!    8-neighbour sphere graph). If a node has no lower neighbour (a pit), the surface is
//!    filled first with Priority-Flood+ε: closed basins fill to their spill point and get a
//!    tiny gradient towards it (geologically: basins silt up into flat plains).
//! 2. **Stack:** a topological order from the coast upstream (every node after its receiver).
//! 3. **Drainage area:** one pass down the reversed stack, accumulating cell areas.
//! 4. **Implicit incision:** one pass up the stack:
//!    `h ← (h + U·Δt + F·h_r) / (1 + F)`, with `F = K·Δt·A^m / L` (L = distance to receiver).
//!    Because the receiver is solved first, the update is implicit and stable at any Δt.
//!
//! Ocean nodes are fixed base level (0 m at the coast; bathymetry is left as tectonics made
//! it). Land never drops below its receiver, so the coastline cannot move.
//!
//! # Uplift
//!
//! `U` comes from tectonic *process*, not from elevation: it decays exponentially with
//! distance from convergent plate boundaries, is stronger on continental crust and scales
//! with the plates' closing speed (see `uplift.rs`). A high plateau far from a boundary gets
//! only the background rate.
//!
//! # Calibration
//!
//! At 10–50 km cells the absolute constants (K, U) are underdetermined, so land heights are
//! finally quantile-mapped onto a land hypsometric curve measured from smoothed global
//! topography ([`hypsometry`]). Dimensionless shape (drainage pattern, concavity, Hack's law)
//! comes from the physics; overall amplitude from the curve. The [`ErosionStats`] check the
//! shape against Earth's laws before the remap.

pub mod hypsometry;
pub mod lakes;
mod solver;
pub mod stats;
mod uplift;

pub use solver::{
    diagonal_mask, erode, erode_climate, erode_climate_with_progress, erode_with_runoff, implicit_step, implicit_step_varying, priority_flood_fill, solve, solve_full, solve_full_with_progress, solve_with_k, Column, ErosionParams, ErosionResult, MassBudget, Network,
    RiverSegment, Solved, UpliftSource,
};
pub use stats::ErosionStats;
