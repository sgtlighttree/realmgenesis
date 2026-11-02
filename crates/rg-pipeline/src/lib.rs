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
use rg_erosion::{ErosionParams, ErosionStats, RiverSegment};
use rg_geo::CubeField;
use rg_tectonics::{TectonicSnapshot, TectonicsParams};

pub use layers::{Equirect, Layer};

/// Everything the user sets before generating.
#[derive(Clone, Debug, Default)]
pub struct WorldParams {
    pub tectonics: TectonicsParams,
    pub erosion: ErosionParams,
    pub planet: PlanetParams,
}

/// Wall-clock timings of the last run of each stage, for the HUD.
#[derive(Clone, Debug, Default)]
pub struct Timings {
    pub tectonics_ms: f64,
    pub erosion_ms: f64,
    pub climate_ms: f64,
}

/// A generated world: the committed fields plus derived data.
/// Hydrology from the erosion stage (committed with the world).
pub struct Hydrology {
    /// Upstream drainage area, km² (erosion grid resolution).
    pub drainage_km2: CubeField<f32>,
    /// Downstream unit flow direction (tangent vectors; zero at sea).
    pub flow_dir: CubeField<[f32; 3]>,
    /// Standing water depth in closed basins, m.
    pub lake_depth_m: CubeField<f32>,
    /// Rock uplift rate, mm/yr (= erosion rate at steady state). None in files before
    /// schema 4.
    pub uplift_mm_yr: Option<CubeField<f32>>,
    pub rivers: Vec<RiverSegment>,
    pub stats: ErosionStats,
}

pub struct World {
    pub params: WorldParams,
    /// `tectonics.height_m` holds the **eroded** macro heights (the committed surface).
    pub tectonics: TectonicSnapshot,
    /// None for worlds generated with erosion off or loaded from schema-1 files.
    pub hydrology: Option<Hydrology>,
    /// Heights resampled to the climate grid. Kept so display-time downscaling can
    /// compare fine heights with what the climate model saw.
    pub climate_heights: LatLonGrid<f32>,
    pub climate: ClimateFields,
    pub timings: Timings,
}

/// Climate grid resolution (1°). See ARCHITECTURE §4.
pub const CLIMATE_W: u32 = 360;
pub const CLIMATE_H: u32 = 180;

/// Generate a world: tectonics → fluvial erosion → tier A climate.
pub fn generate(params: &WorldParams) -> World {
    let mut params = params.clone();
    // One radius for everything; erosion's random streams follow the world seed.
    params.planet.radius_m = params.tectonics.radius_m;
    params.erosion.seed = params.tectonics.seed;

    let t0 = web_time::Instant::now();
    let mut tectonics = rg_tectonics::generate(&params.tectonics);
    let tectonics_ms = t0.elapsed().as_secs_f64() * 1000.0;

    let t1 = web_time::Instant::now();
    let hydrology = params.erosion.enabled.then(|| {
        // The "mountains" setting scales the target hypsometry (Earth = 1.0).
        let mut ep = params.erosion.clone();
        ep.hypsometry_scale = params.tectonics.mountain_scale as f64 * params.erosion.hypsometry_scale;
        let r = rg_erosion::erode(&tectonics, params.tectonics.radius_m, &ep);
        tectonics.height_m = r.height_m;
        Hydrology {
            drainage_km2: r.drainage_km2,
            flow_dir: r.flow_dir,
            lake_depth_m: r.lake_depth_m,
            uplift_mm_yr: Some(r.uplift_mm_yr),
            rivers: r.rivers,
            stats: r.stats,
        }
    });
    let erosion_ms = t1.elapsed().as_secs_f64() * 1000.0;

    let (climate_heights, climate, climate_ms) = run_climate(&params.planet, &tectonics);
    World {
        params,
        tectonics,
        hydrology,
        climate_heights,
        climate,
        timings: Timings { tectonics_ms, erosion_ms, climate_ms },
    }
}

/// The streamed-detail sampler for a world (macro heights + river corridors). Use this
/// everywhere (app, tests, probes) so they all see the same surface.
pub fn terrain_sampler(world: &World) -> rg_terrain::TerrainSampler {
    let s = rg_terrain::TerrainSampler::new(
        std::sync::Arc::new(world.tectonics.height_m.clone()),
        world.params.tectonics.radius_m,
        world.params.tectonics.seed,
    );
    match &world.hydrology {
        Some(h) => s.with_drainage(std::sync::Arc::new(h.drainage_km2.clone())),
        None => s,
    }
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
