//! Orchestration: runs the simulations in order and bakes display layers.
//!
//! The app (UI) calls into this crate and never runs simulations itself. M1 turns these
//! calls into stepped `Job`s (ARCHITECTURE §5b). For the MVP they're plain functions
//! run between frames.

pub mod export;
pub mod file;
pub mod layers;
mod palette;

use rg_climate::{ClimateFields, ClimateInput, LatLonGrid, PlanetParams};
use rg_tectonics::{TectonicSnapshot, TectonicsParams};

pub use layers::{Equirect, Layer};

/// Everything the user sets before generating.
#[derive(Clone, Debug, Default)]
pub struct WorldParams {
    pub tectonics: TectonicsParams,
    pub planet: PlanetParams,
}

/// Wall-clock timings of the last run of each stage, for the HUD.
#[derive(Clone, Debug, Default)]
pub struct Timings {
    pub tectonics_ms: f64,
    pub climate_ms: f64,
}

/// A generated world: the committed fields plus derived data.
pub struct World {
    pub params: WorldParams,
    pub tectonics: TectonicSnapshot,
    /// Heights resampled to the climate grid. Kept so display-time downscaling can
    /// compare fine heights with what the climate model saw.
    pub climate_heights: LatLonGrid<f32>,
    pub climate: ClimateFields,
    pub timings: Timings,
}

/// Climate grid resolution (1°). See ARCHITECTURE §4.
pub const CLIMATE_W: u32 = 360;
pub const CLIMATE_H: u32 = 180;

/// Generate a world: tectonics → macro heights → tier A climate.
pub fn generate(params: &WorldParams) -> World {
    let mut params = params.clone();
    // One radius for everything.
    params.planet.radius_m = params.tectonics.radius_m;

    let t0 = web_time::Instant::now();
    let tectonics = rg_tectonics::generate(&params.tectonics);
    let tectonics_ms = t0.elapsed().as_secs_f64() * 1000.0;

    let (climate_heights, climate, climate_ms) = run_climate(&params.planet, &tectonics);
    World { params, tectonics, climate_heights, climate, timings: Timings { tectonics_ms, climate_ms } }
}

/// Re-run only the climate (after changing planet parameters).
pub fn rerun_climate(world: &mut World, planet: &PlanetParams) {
    let mut planet = *planet;
    planet.radius_m = world.params.tectonics.radius_m;
    world.params.planet = planet;
    let (h, c, ms) = run_climate(&planet, &world.tectonics);
    world.climate_heights = h;
    world.climate = c;
    world.timings.climate_ms = ms;
}

fn run_climate(planet: &PlanetParams, t: &TectonicSnapshot) -> (LatLonGrid<f32>, ClimateFields, f64) {
    let t0 = web_time::Instant::now();
    let heights = LatLonGrid::from_cube_field(&t.height_m, CLIMATE_W, CLIMATE_H);
    let climate = rg_climate::tier_a(planet, &ClimateInput { height_m: heights.clone() });
    (heights, climate, t0.elapsed().as_secs_f64() * 1000.0)
}
