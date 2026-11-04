//! Orchestration: runs the simulations in order and bakes display layers.
//!
//! The app (UI) calls into this crate and never runs simulations itself. M1 turns these
//! calls into stepped `Job`s (ARCHITECTURE §5b). For the MVP they're plain functions
//! run between frames.

pub mod export;
pub mod file;
pub mod import;
pub mod layers;
pub mod legend;
mod palette;
pub mod svg;
/// The coast mask types (R3), re-exported for the app.
pub use rg_tectonics::mask;

use rg_climate::{ClimateFields, ClimateInput, LatLonGrid, PlanetParams};
use rg_erosion::{ErosionParams, ErosionStats, RiverSegment};
use rg_geo::CubeField;
use rg_tectonics::{TectonicSnapshot, TectonicsParams};

pub use layers::{Equirect, Layer};
pub use legend::Legend;

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
    /// Solved → committed height factor under linear calibration (k_sn = factor·U/K).
    pub height_scale: Option<f64>,
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

/// A generation stage, for progress reporting (the app's status line and console, s02).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Stage {
    Tectonics,
    /// Tier A climate on the tectonic surface (feeds erosion's runoff and lakes).
    ClimatePre,
    /// The erosion solve, reported per implicit step.
    Erosion,
    /// Tier A climate on the eroded surface.
    Climate,
    /// Detailing a finer world from the canonical one (`rg_erosion::refine`).
    Refine,
}

impl Stage {
    pub fn label(self) -> &'static str {
        match self {
            Stage::Tectonics => "tectonics",
            Stage::ClimatePre => "climate (pre-erosion)",
            Stage::Erosion => "erosion",
            Stage::Climate => "climate",
            Stage::Refine => "refine (detail)",
        }
    }
}

/// Where generation is: `step` of `steps` within `stage` (both 0 when the stage has no
/// sub-steps).
#[derive(Clone, Copy, Debug)]
pub struct Progress {
    pub stage: Stage,
    pub step: u32,
    pub steps: u32,
}

// The app moves a finished World from its generation thread to the UI thread (s02).
#[allow(dead_code)]
fn assert_world_is_send() {
    fn is_send<T: Send>() {}
    is_send::<World>();
}

/// True when a world is coarser than the canonical resolution: a **seed scanner**
/// (coastlines, plates, biomes, climate; no erosion). Erosion there would be its own,
/// different solve and wouldn't match Normal (s04, Matt: Draft as a seed scanner).
/// Coastlines and plates already match at every resolution (s03).
pub fn is_seed_scanner(params: &WorldParams) -> bool {
    let canon = params.erosion.canonical_res;
    canon > 0 && params.tectonics.res < canon
}

/// Generate a world: tectonics → fluvial erosion → tier A climate.
pub fn generate(params: &WorldParams) -> World {
    generate_with_progress(params, &|_| true).expect("generation without a cancelling callback always finishes")
}

/// As [`generate`], calling `on` at the start of each stage and each erosion step, on this
/// thread. If `on` returns false, generation stops at that checkpoint and returns `None`
/// (cancelled). Tectonics and the post-solve part of erosion have no checkpoints inside.
pub fn generate_with_progress(params: &WorldParams, on: &dyn Fn(Progress) -> bool) -> Option<World> {
    let stage = |stage: Stage| on(Progress { stage, step: 0, steps: 0 });
    let mut params = params.clone();
    // One radius for everything; erosion's random streams follow the world seed.
    params.planet.radius_m = params.tectonics.radius_m;
    params.erosion.seed = params.tectonics.seed;

    // Finer than the canonical world (Normal): erode at the canonical resolution, then
    // refine. The large-scale shape is then the same at every preset (s04).
    let canon = params.erosion.canonical_res;
    if params.erosion.enabled && canon > 0 && params.tectonics.res > canon {
        return generate_refined(params, canon, on);
    }
    let erode_now = params.erosion.enabled && !is_seed_scanner(&params);

    let t0 = web_time::Instant::now();
    if !stage(Stage::Tectonics) {
        return None;
    }
    let mut tectonics = rg_tectonics::generate(&params.tectonics);
    let tectonics_ms = t0.elapsed().as_secs_f64() * 1000.0;

    let t1 = web_time::Instant::now();
    let hydrology = if erode_now {
        // The "mountains" setting scales the target hypsometry (Earth = 1.0).
        let mut ep = params.erosion.clone();
        ep.hypsometry_scale = params.tectonics.mountain_scale as f64 * params.erosion.hypsometry_scale;
        // Climate feeds erosion: a first tier-A pass on the tectonic surface gives annual
        // precipitation; rivers carve with precipitation-weighted drainage (a discharge
        // proxy), so rain shadows and deserts erode less. Climate is rerun on the result.
        if !stage(Stage::ClimatePre) {
            return None;
        }
        let (_, climate0, _) = run_climate(&params.planet, &tectonics);
        let wb = water_balance(&params.planet, &climate0, &tectonics.height_m);
        let mean = wb.mean_land_runoff_mm.max(1.0) as f32;
        let runoff = CubeField { res: wb.runoff_mm.res, data: wb.runoff_mm.data.iter().map(|v| v / mean).collect() };
        // Lakes and closed basins from the water balance: open-water net evaporation ≈ PET − P.
        let net_evap = CubeField { res: wb.pet_mm.res, data: wb.pet_mm.data.iter().zip(&wb.precip_mm.data).map(|(e, p)| e - p).collect() };
        let lakes = ep.lakes.then_some((&wb.runoff_mm, &net_evap));
        if !stage(Stage::Erosion) {
            return None;
        }
        let on_step = |step, steps| on(Progress { stage: Stage::Erosion, step, steps });
        let r = rg_erosion::erode_climate_with_progress(&tectonics, params.tectonics.radius_m, &ep, Some(&runoff), lakes, &on_step)?;
        tectonics.height_m = r.height_m;
        Some(Hydrology {
            drainage_km2: r.drainage_km2,
            flow_dir: r.flow_dir,
            lake_depth_m: r.lake_depth_m,
            uplift_mm_yr: Some(r.uplift_mm_yr),
            height_scale: r.height_scale,
            rivers: r.rivers,
            stats: r.stats,
        })
    } else {
        None
    };
    let erosion_ms = t1.elapsed().as_secs_f64() * 1000.0;

    if !stage(Stage::Climate) {
        return None;
    }
    let (climate_heights, climate, climate_ms) = run_climate(&params.planet, &tectonics);
    Some(World {
        params,
        tectonics,
        hydrology,
        climate_heights,
        climate,
        timings: Timings { tectonics_ms, erosion_ms, climate_ms },
    })
}

/// Tectonics → climate pre-pass → erosion, the canonical way: the same steps (and bytes)
/// as `generate_with_progress` for a world at this resolution. Returns the eroded snapshot,
/// the erosion result and the relative runoff the solve used.
fn erode_canonical(params: &WorldParams, tectonics: &TectonicSnapshot, on: &dyn Fn(Progress) -> bool) -> Option<(rg_erosion::ErosionResult, CubeField<f32>)> {
    let stage = |stage: Stage| on(Progress { stage, step: 0, steps: 0 });
    let mut ep = params.erosion.clone();
    ep.hypsometry_scale = params.tectonics.mountain_scale as f64 * params.erosion.hypsometry_scale;
    if !stage(Stage::ClimatePre) {
        return None;
    }
    let (_, climate0, _) = run_climate(&params.planet, tectonics);
    let wb = water_balance(&params.planet, &climate0, &tectonics.height_m);
    let mean = wb.mean_land_runoff_mm.max(1.0) as f32;
    let runoff = CubeField { res: wb.runoff_mm.res, data: wb.runoff_mm.data.iter().map(|v| v / mean).collect() };
    let net_evap = CubeField { res: wb.pet_mm.res, data: wb.pet_mm.data.iter().zip(&wb.precip_mm.data).map(|(e, p)| e - p).collect() };
    let lakes = ep.lakes.then_some((&wb.runoff_mm, &net_evap));
    if !stage(Stage::Erosion) {
        return None;
    }
    let on_step = |step, steps| on(Progress { stage: Stage::Erosion, step, steps });
    let r = rg_erosion::erode_climate_with_progress(tectonics, params.tectonics.radius_m, &ep, Some(&runoff), lakes, &on_step)?;
    Some((r, runoff))
}

/// A world finer than the canonical resolution: the canonical world (Normal) is generated
/// and eroded exactly as at its own preset, then detailed onto the fine grid
/// (`rg_erosion::refine`). See `docs/design/canonical-refine.md`.
fn generate_refined(params: WorldParams, canon: u32, on: &dyn Fn(Progress) -> bool) -> Option<World> {
    let stage = |stage: Stage| on(Progress { stage, step: 0, steps: 0 });
    let t0 = web_time::Instant::now();
    if !stage(Stage::Tectonics) {
        return None;
    }
    let canon_tp = TectonicsParams { res: canon, ..params.tectonics.clone() };
    let canon_tect = rg_tectonics::generate(&canon_tp);
    let mut fine = rg_tectonics::generate(&params.tectonics);
    let tectonics_ms = t0.elapsed().as_secs_f64() * 1000.0;

    let t1 = web_time::Instant::now();
    let canon_params = WorldParams { tectonics: canon_tp, ..params.clone() };
    let (r, runoff) = erode_canonical(&canon_params, &canon_tect, on)?;
    if !stage(Stage::Refine) {
        return None;
    }
    let input = rg_erosion::refine::RefineInput {
        canon_h: &r.height_m,
        canon_t: &canon_tect.height_m,
        canon_lake_m: &r.lake_depth_m,
        canon_closed: &r.closed_lake,
        fine_t: &fine.height_m,
        height_scale: r.height_scale.unwrap_or(1.0),
        runoff: Some(&runoff),
    };
    let mut ep = params.erosion.clone();
    ep.hypsometry_scale = params.tectonics.mountain_scale as f64 * params.erosion.hypsometry_scale;
    let refined = rg_erosion::refine::refine(&input, params.tectonics.radius_m, &ep);
    fine.height_m = refined.height_m;
    // Uplift is the canonical rate, sampled onto the hydrology grid.
    let mut uplift = CubeField::<f32>::new(refined.drainage_km2.res);
    let samples: Vec<_> = uplift.iter_samples().map(|(k, _, _, _, d)| (k, d)).collect();
    for (k, d) in samples {
        uplift.data[k] = r.uplift_mm_yr.sample(d);
    }
    let hydrology = Some(Hydrology {
        drainage_km2: refined.drainage_km2,
        flow_dir: refined.flow_dir,
        lake_depth_m: refined.lake_depth_m,
        uplift_mm_yr: Some(uplift),
        height_scale: r.height_scale,
        rivers: refined.rivers,
        stats: r.stats,
    });
    let erosion_ms = t1.elapsed().as_secs_f64() * 1000.0;

    if !stage(Stage::Climate) {
        return None;
    }
    let (climate_heights, climate, climate_ms) = run_climate(&params.planet, &fine);
    Some(World { params, tectonics: fine, hydrology, climate_heights, climate, timings: Timings { tectonics_ms, erosion_ms, climate_ms } })
}

/// The displayed terrain surface of a world (the committed macro heights, smoothly
/// interpolated; `rg_terrain::TerrainSampler`). Use this everywhere (app, tests, probes) so
/// they all see the same surface.
pub fn terrain_sampler(world: &World) -> rg_terrain::TerrainSampler {
    let height = std::sync::Arc::new(world.tectonics.height_m.clone());
    rg_terrain::TerrainSampler::new(height, world.params.tectonics.radius_m)
}

/// Land water balance on a climate run: annual precipitation, potential evapotranspiration
/// (Thornthwaite) and runoff (Budyko–Fu), each on a cube field (res 128), mm/yr.
pub struct WaterBalance {
    pub precip_mm: CubeField<f32>,
    pub pet_mm: CubeField<f32>,
    pub runoff_mm: CubeField<f32>,
    /// Mean land runoff, mm/yr (Earth ≈ 300: ~40,000 km³/yr over ~135 M km²).
    pub mean_land_runoff_mm: f64,
}

pub fn water_balance(planet: &PlanetParams, climate: &ClimateFields, height: &CubeField<f32>) -> WaterBalance {
    let (w, h) = (climate.width, climate.height);
    let mut pr = Vec::with_capacity((w * h) as usize);
    let mut pe = Vec::with_capacity((w * h) as usize);
    let mut ro = Vec::with_capacity((w * h) as usize);
    for y in 0..h {
        let lat = layers::lat_of_row(y, h);
        for x in 0..w {
            let k = climate.index(x, y);
            let p: f64 = climate.precip_mm[k].iter().map(|&v| v as f64).sum();
            let pet = rg_climate::water::thornthwaite_pet_mm(&climate.temp_c[k], lat, planet.axial_tilt_deg);
            pr.push(p as f32);
            pe.push(pet as f32);
            ro.push(rg_climate::water::budyko_fu_runoff_mm(p, pet, rg_climate::water::FU_OMEGA) as f32);
        }
    }
    let to_cube = |v: Vec<f32>| {
        let g = LatLonGrid::from_vec(w, h, v);
        CubeField::<f32>::from_fn(128, |d| g.sample(d.z.clamp(-1.0, 1.0).asin(), d.y.atan2(d.x)))
    };
    let (precip_mm, pet_mm, runoff_mm) = (to_cube(pr), to_cube(pe), to_cube(ro));
    // Land mean (sample density on the cube is uneven; the tangent warp keeps it within
    // ~30%, fine for a normaliser and a sanity check).
    let (mut sum, mut cnt) = (0.0f64, 0usize);
    for (k, _, _, _, d) in runoff_mm.iter_samples() {
        if height.sample(d) > 0.0 {
            sum += runoff_mm.data[k] as f64;
            cnt += 1;
        }
    }
    WaterBalance { precip_mm, pet_mm, runoff_mm, mean_land_runoff_mm: sum / cnt.max(1) as f64 }
}

/// Runoff on the tectonic surface relative to its land mean (1 = average): the erosion
/// solve's discharge proxy. Runoff (Budyko), not precipitation: deserts shed far less water
/// than their rainfall suggests.
fn relative_runoff(planet: &PlanetParams, t: &TectonicSnapshot) -> CubeField<f32> {
    let (_, climate, _) = run_climate(planet, t);
    let wb = water_balance(planet, &climate, &t.height_m);
    let mean = wb.mean_land_runoff_mm.max(1.0) as f32;
    CubeField { res: wb.runoff_mm.res, data: wb.runoff_mm.data.iter().map(|v| v / mean).collect() }
}

/// [`relative_runoff`] for probes.
pub fn relative_runoff_for_probe(planet: &PlanetParams, t: &TectonicSnapshot) -> CubeField<f32> {
    relative_runoff(planet, t)
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

#[cfg(test)]
mod progress_tests {
    use super::*;

    fn tiny() -> WorldParams {
        let mut p = WorldParams::default();
        p.tectonics.res = 32;
        // A tiny world is coarser than Normal, so it would be a seed scanner (no erosion,
        // s04). These tests are about erosion progress: erode it directly.
        p.erosion.canonical_res = 0;
        p
    }

    #[test]
    fn cancelling_at_an_erosion_step_stops_generation() {
        let seen = std::cell::RefCell::new(Vec::new());
        let world = generate_with_progress(&tiny(), &|p| {
            seen.borrow_mut().push(p.stage);
            // Stop at the second erosion step.
            !(p.stage == Stage::Erosion && p.step == 2)
        });
        assert!(world.is_none(), "cancelled generation must not return a world");
        let seen = seen.into_inner();
        assert_eq!(seen.first(), Some(&Stage::Tectonics));
        assert!(!seen.contains(&Stage::Climate), "stopped before the final climate pass");
    }

    #[test]
    fn progress_reports_every_erosion_step() {
        let steps = std::cell::RefCell::new(Vec::new());
        let world = generate_with_progress(&tiny(), &|p| {
            if p.stage == Stage::Erosion && p.steps > 0 {
                steps.borrow_mut().push((p.step, p.steps));
            }
            true
        });
        assert!(world.is_some());
        let steps = steps.into_inner();
        let n = steps.last().map(|s| s.1).unwrap_or(0);
        assert!(n > 0 && steps.len() == n as usize, "{steps:?}");
    }
}

/// Resolution coherence guard (s04, Matt: "permanently fix or avoid popping up again").
/// A world finer than the canonical resolution must keep the canonical world's large-scale
/// shape: averaged over each canonical cell it matches the canonical world, and it shares
/// its height calibration. Before s04, erosion ran at full resolution and Ultra lost whole
/// massifs that Normal had (seed 1: height factor 0.75 vs 0.57).
#[cfg(test)]
mod coherence_tests {
    use super::*;
    use glam::DVec3;

    fn world(res: u32, canonical_res: u32) -> World {
        let mut p = WorldParams::default();
        p.tectonics.seed = 7;
        p.tectonics.res = res;
        p.erosion.canonical_res = canonical_res;
        generate(&p)
    }

    #[test]
    fn coarser_than_canonical_is_a_seed_scanner_with_the_same_coast() {
        // 128/96, not 64/48: at ~200 km cells the coastal band alone is ~6% of samples.
        let canon = world(128, 128);
        let scan = world(96, 128);
        assert!(scan.hydrology.is_none(), "a scanner world doesn't erode");
        let (hc, hs) = (&canon.tectonics.height_m, &scan.tectonics.height_m);
        let (mut same, mut total) = (0usize, 0usize);
        for (k, _, _, _, d) in hs.iter_samples() {
            total += 1;
            same += ((hs.data[k] > 0.0) == (hc.sample(d) > 0.0)) as usize;
        }
        let share = same as f64 / total as f64;
        assert!(share > 0.95, "land/sea agreement {share:.3}");
    }

    #[test]
    fn refined_world_keeps_the_canonical_shape() {
        let canon = world(64, 64);
        let fine = world(128, 64);
        let (hc, hf) = (&canon.tectonics.height_m, &fine.tectonics.height_m);
        assert_eq!(hf.res, 128);
        let scale = |w: &World| w.hydrology.as_ref().and_then(|h| h.height_scale);
        assert_eq!(scale(&canon), scale(&fine), "refinement must not recalibrate");

        // Fine heights averaged over each canonical cell (5×5 points), on canonical land.
        let cell = hc.spacing_rad();
        let (mut xs, mut ys) = (Vec::new(), Vec::new());
        for (k, _, _, _, d) in hc.iter_samples() {
            if hc.data[k] <= 0.0 || hf.sample(d) <= 0.0 {
                continue;
            }
            let east = DVec3::Z.cross(d).normalize_or(DVec3::X);
            let north = d.cross(east);
            let mean = (0..25)
                .map(|q| {
                    let (x, y) = ((q % 5) as f64 / 4.0 - 0.5, (q / 5) as f64 / 4.0 - 0.5);
                    hf.sample_f64((d + (east * x + north * y) * cell).normalize())
                })
                .sum::<f64>()
                / 25.0;
            xs.push(hc.data[k] as f64);
            ys.push(mean);
        }
        let n = xs.len() as f64;
        let rms = (xs.iter().zip(&ys).map(|(a, b)| (a - b).powi(2)).sum::<f64>() / n).sqrt();
        let q = |v: &mut Vec<f64>, f: f64| {
            v.sort_by(f64::total_cmp);
            v[((v.len() - 1) as f64 * f) as usize]
        };
        let (c50, f50) = (q(&mut xs.clone(), 0.5), q(&mut ys.clone(), 0.5));
        let (c99, f99) = (q(&mut xs.clone(), 0.99), q(&mut ys.clone(), 0.99));
        eprintln!("land cells {n} · rms {rms:.0} m · p50 {c50:.0}/{f50:.0} · p99 {c99:.0}/{f99:.0}");
        assert!((f50 - c50).abs() <= 0.1 * c50.max(100.0), "p50 {c50:.0} vs {f50:.0}");
        assert!((f99 - c99).abs() <= 0.1 * c99, "p99 {c99:.0} vs {f99:.0}");
        assert!(rms <= 0.1 * c99, "rms {rms:.0} m vs p99 {c99:.0}");
    }

    /// R3.1b gate 1 (erosion half): with strict coasts, erosion's ocean is the mask at its
    /// own nodes, so the eroded surface's sign matches the mask everywhere, and no lake bed
    /// or valley drops to sea level on land.
    #[test]
    fn strict_coasts_survive_erosion() {
        let mut params = WorldParams::default();
        params.tectonics.res = 48;
        params.erosion.canonical_res = 0;
        params.erosion.steps = 8;
        params.tectonics.conditioning.strict = true;
        params.tectonics.conditioning.islands = false;
        let snap = rg_tectonics::generate(&params.tectonics);
        let cover = snap.land_cover.as_ref().unwrap();
        let w = crate::generate(&params);
        let hm = &w.tectonics.height_m;
        let mut checked = 0;
        for (idx, _, _, _, d) in hm.iter_samples() {
            let land = cover.sample(d) > 0.5;
            let h = hm.data[idx];
            assert_eq!(land, h > 0.0, "sample {idx}: mask land {land}, height {h}");
            checked += 1;
        }
        assert!(checked > 0);
    }

    /// The same gate through the refine path (High/XH/Ultra = canonical + refine): the
    /// refined surface's sign must still be the fine mask's (advisor, s05).
    #[test]
    fn strict_coasts_survive_refine() {
        let mut params = WorldParams::default();
        params.tectonics.res = 96;
        params.erosion.canonical_res = 48;
        params.erosion.steps = 8;
        params.tectonics.conditioning.strict = true;
        params.tectonics.conditioning.islands = false;
        let snap = rg_tectonics::generate(&params.tectonics);
        let cover = snap.land_cover.as_ref().unwrap();
        let w = crate::generate(&params);
        let hm = &w.tectonics.height_m;
        assert!(hm.res > params.erosion.canonical_res, "refine path taken");
        let mut miss = 0;
        for (idx, _, _, _, d) in hm.iter_samples() {
            if (cover.sample(d) > 0.5) != (hm.data[idx] > 0.0) {
                miss += 1;
            }
        }
        assert_eq!(miss, 0, "{miss} samples off the mask");
    }
}

