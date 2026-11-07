//! A present-day plate-tectonics **snapshot**: macro relief for a planet, not a history.
//!
//! [`generate`] takes a seed and a few knobs and produces, on the cube-sphere:
//! heights (metres, 0 = sea level), plate ids, crust type and boundary types.
//!
//! How it works, in the order the code runs (each step has its own module):
//!
//! 1. [`crust`]: continental vs oceanic crust from its own noise and a handful of
//!    continental "cores". Independent of plates, so continents are not plate silhouettes.
//! 2. [`plates`]: plates grow from seeds over the sample graph with noisy costs (irregular,
//!    connected, varied in size). Each plate gets an Euler pole and rotation rate. At every
//!    border we compare the two plates' velocities to decide: convergent, divergent or transform.
//! 3. [`relief`]: heights built from distance-to-boundary profiles (mountains, trenches, ridges,
//!    rifts), a seafloor-age depth law, shelves, old eroded belts, and context-scaled noise.
//! 4. Sea level: a vertical offset found by bisection so the land share matches the target.
//!
//! Everything is deterministic (same params, same bytes) and uses no platform clock or hash
//! maps. It uses threads only on native targets (see `util::par_map`),
//! where they change the speed but never the result, so it runs the same on wasm.

pub mod coast;
mod crust;
pub mod fit;
pub mod grid;
pub mod mask;
pub mod plates;
mod relief;
mod rng;
pub mod util;

use glam::DVec3;
use rg_geo::CubeField;

pub(crate) type Params = TectonicsParams;

/// Knobs for [`generate`].
#[derive(Clone, Debug, PartialEq)]
pub struct TectonicsParams {
    pub seed: u64,
    /// Cube-sphere cells per face edge (the output fields have `res + 1` samples per edge).
    pub res: u32,
    /// Planet radius in metres (sets the physical size of mountain belts, shelves, etc.).
    pub radius_m: f64,
    /// Number of tectonic plates.
    pub plate_count: u32,
    /// Target share of the surface above sea level, 0..1.
    pub land_fraction: f32,
    /// Roughly how many large continents to make (smaller landmasses are added on top).
    pub continent_count_hint: u32,
    /// Multiplier on mountain-building amplitudes.
    pub mountain_scale: f32,
    /// Multiplier on fine terrain noise.
    pub roughness: f32,
    /// Target share of continental coast length that is active (trench + cordillera or arc).
    /// Earth: ~0.23 (Bradley 2008; design doc). 0 = every seaway opening ("all passive").
    pub active_share: f32,
    /// Where the coasts come from (R3): procedural and emergent (default), procedural strict,
    /// or an imported mask; plus brush strokes.
    pub conditioning: mask::Conditioning,
}

impl Default for TectonicsParams {
    fn default() -> Self {
        Self {
            seed: 1,
            res: 320,
            radius_m: 6_371_000.0,
            plate_count: 12,
            land_fraction: 0.29,
            continent_count_hint: 5,
            mountain_scale: 1.0,
            roughness: 1.0,
            active_share: 0.25,
            conditioning: mask::Conditioning::default(),
        }
    }
}

/// A rigid plate rotating about an Euler pole.
#[derive(Clone, Debug, PartialEq)]
pub struct Plate {
    /// Unit vector: the axis the plate rotates about.
    pub euler_pole: DVec3,
    /// Angular rate, radians per million years. The velocity at unit-sphere point `p` is
    /// `omega * (euler_pole x p)` (radians of arc per Myr; times the radius for metres).
    pub omega: f64,
    /// True if most of the plate's area is continental crust.
    pub continental: bool,
}

/// Boundary classes stored in [`TectonicSnapshot::boundary`].
pub mod boundary_class {
    pub const NONE: u8 = 0;
    pub const CONVERGENT: u8 = 1;
    pub const DIVERGENT: u8 = 2;
    pub const TRANSFORM: u8 = 3;
}

/// The result of [`generate`].
#[derive(Clone, Debug)]
pub struct TectonicSnapshot {
    /// Height in metres relative to sea level (0).
    pub height_m: CubeField<f32>,
    pub plate_id: CubeField<u16>,
    /// 0 oceanic, 1 continental.
    pub crust: CubeField<u8>,
    /// 0 none, 1 convergent, 2 divergent, 3 transform. Set on both sides of each border.
    pub boundary: CubeField<u8>,
    pub plates: Vec<Plate>,
    /// The vertical shift that was applied to hit the land-fraction target (metres).
    pub sea_level_offset_m: f32,
    /// Achieved share of the sphere's area above sea level.
    pub land_fraction: f32,
    /// Coast class per sample near coasts: 0 none, 1 passive, 2 active (`coast::margin`).
    /// R3; absent in worlds saved before schema 5.
    pub margin: Option<CubeField<u8>>,
    /// Clamped convergence rate at convergent boundary samples, mm/yr (0 elsewhere). R3: the
    /// rate the relief used, committed so erosion's uplift reads the same one (absent in old
    /// worlds; uplift then falls back to the plates' rotations).
    pub convergence_mm_yr: Option<CubeField<f32>>,
    /// Strict coasts (R3.1b): 1 land, 0 sea, per sample. Erosion takes its ocean from this
    /// (sampled at its own nodes) instead of the sign of the resampled heights, so coastal
    /// nodes can't flip. Only an erosion input; not stored in world files.
    pub land_cover: Option<CubeField<f32>>,
    /// Continental dynamic topography to build during erosion (m), if generated in this run.
    /// Not stored in world files: only erosion uses it, and its result is committed.
    pub dynamic_topo_m: Option<CubeField<f32>>,
    /// Heightmap pin weight 0..1 per sample (s07): erosion keeps pinned land's heights out
    /// of calibration and inverts uplift from them. Derived from the conditioning (stored);
    /// an erosion input only, not stored in world files. `None` without heightmap layers.
    pub pin: Option<CubeField<f32>>,
}

/// Area-weighted fraction of a height field above 0 (cube-sphere cells near face corners are
/// smaller than central ones, so a plain sample count would be biased).
pub fn land_fraction_of(height: &CubeField<f32>) -> f32 {
    let g = grid::Grid::new(height.res);
    let mut f = 0.0f64;
    for a in 0..g.n {
        if height.data[g.node_sample[a] as usize] > 0.0 {
            f += g.weight[a] as f64;
        }
    }
    f as f32
}

/// Build a tectonic snapshot.
pub fn generate(p: &TectonicsParams) -> TectonicSnapshot {
    generate_with_model(p).0
}

/// [`generate`], plus the plate model and the fit statistics (for probes and tests).
pub fn generate_with_model(p: &TectonicsParams) -> (TectonicSnapshot, plates::PlateModel, fit::FitStats) {
    let p = &TectonicsParams { res: p.res.max(8), ..p.clone() };
    let grid = grid::Grid::new(p.res);
    let tm = std::time::Instant::now();
    let lap = |what: &str| if std::env::var_os("RG_TIMING").is_some() { eprintln!("  {what}: {} ms", tm.elapsed().as_millis()) };
    // R3.1a: classify the coasts on the fixed grid, then fit plates and motion to them.
    let cls = coast::Classification::build(p);
    lap("classification");
    let strict = p.conditioning.strict;
    let mut crust = if strict {
        // R3.1b: the land mask first; crust = land + class-dependent margins.
        let value = crust::continentalness(&grid, p);
        let land: Vec<bool> = mask::land_cover(&grid, p, &value).iter().map(|&c| c > 0.5).collect();
        crust::strict(&grid, p, &land, &cls)
    } else {
        crust::build(&grid, p)
    };
    lap("crust");
    let (pm, stats) = fit::build(&grid, p, &crust, &cls);
    lap("fit");
    crust.activity = fit::activity(&grid, &cls);
    let composed = relief::compose(&grid, p, &crust, &pm);
    lap("relief");

    // Sea level: the offset at which `land_fraction` of the area is above water.
    // Strict coasts: no bisection, the heights conform to the mask instead.
    let offset = if strict {
        0.0
    } else {
        let target = p.land_fraction.clamp(0.0, 1.0) as f64;
        util::solve_threshold(&composed.raw(), &grid.weight, target)
    };
    // Heightmap layers (s07): imported heights and their pin weights.
    let imported = if strict { mask::imported_heights(&grid, p) } else { None };
    let height_node = relief::finish(&grid, p, &crust, &composed, offset, strict, imported.as_ref());
    let achieved = grid.fraction(|a| height_node[a] > 0.0) as f32;

    // Copy node values out to every cube sample (twins on face edges get identical values).
    let mut snap = TectonicSnapshot {
        height_m: CubeField::new(p.res),
        plate_id: CubeField::new(p.res),
        crust: CubeField::new(p.res),
        boundary: CubeField::new(p.res),
        margin: Some(CubeField::new(p.res)),
        convergence_mm_yr: Some(CubeField::new(p.res)),
        land_cover: strict.then(|| CubeField::new(p.res)),
        plates: pm.plates.clone(),
        sea_level_offset_m: offset,
        land_fraction: achieved,
        // Static offset in the heights (see relief.rs); the rate mechanism is off.
        dynamic_topo_m: None,
        pin: imported.as_ref().map(|_| CubeField::new(p.res)),
    };
    for (s, &node) in grid.sample_node.iter().enumerate() {
        let a = node as usize;
        snap.height_m.data[s] = height_node[a];
        snap.plate_id.data[s] = pm.plate_of[a];
        snap.crust.data[s] = crust.cont[a] as u8;
        snap.boundary.data[s] = pm.boundary[a];
        if let Some(f) = snap.margin.as_mut() {
            f.data[s] = pm.margin[a];
        }
        if let Some(f) = snap.convergence_mm_yr.as_mut() {
            f.data[s] = pm.convergence[a];
        }
        if let Some(f) = snap.land_cover.as_mut() {
            f.data[s] = if height_node[a] > 0.0 { 1.0 } else { 0.0 };
        }
        if let Some(f) = snap.dynamic_topo_m.as_mut() {
            f.data[s] = composed.dyn_topo[a];
        }
        if let (Some(f), Some(im)) = (snap.pin.as_mut(), imported.as_ref()) {
            f.data[s] = im.pin[a];
        }
    }
    (snap, pm, stats)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::OnceLock;

    fn default_snapshot() -> &'static TectonicSnapshot {
        static SNAP: OnceLock<TectonicSnapshot> = OnceLock::new();
        SNAP.get_or_init(|| generate(&TectonicsParams::default()))
    }

    #[test]
    fn deterministic_same_params_same_bytes() {
        let p = TectonicsParams { res: 64, seed: 99, ..Default::default() };
        let a = generate(&p);
        let b = generate(&p);
        let bits = |f: &CubeField<f32>| f.data.iter().map(|v| v.to_bits()).collect::<Vec<_>>();
        assert_eq!(bits(&a.height_m), bits(&b.height_m));
        assert_eq!(a.plate_id, b.plate_id);
        assert_eq!(a.crust, b.crust);
        assert_eq!(a.boundary, b.boundary);
        assert_eq!(a.plates, b.plates);
        assert_eq!(a.sea_level_offset_m.to_bits(), b.sea_level_offset_m.to_bits());
        // A different seed must actually differ.
        let c = generate(&TectonicsParams { seed: 100, ..p });
        assert_ne!(bits(&a.height_m), bits(&c.height_m));
    }

    #[test]
    fn land_fraction_hits_target() {
        let s = default_snapshot();
        let measured = land_fraction_of(&s.height_m);
        assert!((measured - 0.29).abs() < 0.005, "measured {measured}");
        assert!((s.land_fraction - 0.29).abs() < 0.005, "reported {}", s.land_fraction);
        for target in [0.15f32, 0.45] {
            let q = generate(&TectonicsParams { res: 96, land_fraction: target, ..Default::default() });
            let m = land_fraction_of(&q.height_m);
            assert!((m - target).abs() < 0.005, "target {target} measured {m}");
        }
    }

    #[test]
    fn heights_are_plausible() {
        let s = default_snapshot();
        let (lo, hi) = s.height_m.min_max();
        assert!(s.height_m.data.iter().all(|h| h.is_finite()));
        assert!((4000.0..=9000.0).contains(&hi), "max {hi}");
        assert!((-11_500.0..=-6000.0).contains(&lo), "min {lo}");
        // Typical abyssal depth: median of the sea floor between -3500 and -6000.
        let mut sea: Vec<f32> = s.height_m.data.iter().copied().filter(|&h| h < 0.0).collect();
        sea.sort_by(|a, b| a.total_cmp(b));
        let median = sea[sea.len() / 2];
        assert!((-6000.0..=-3500.0).contains(&median), "median sea floor {median}");
        // Land should not be one flat plain nor all mountains.
        let mut land: Vec<f32> = s.height_m.data.iter().copied().filter(|&h| h > 0.0).collect();
        land.sort_by(|a, b| a.total_cmp(b));
        let p50 = land[land.len() / 2];
        let p99 = land[land.len() * 99 / 100];
        assert!((200.0..=1500.0).contains(&p50), "median land {p50}");
        assert!(p99 > 3000.0, "p99 land {p99}");
    }

    /// Fraction of a plate's samples in its largest connected component.
    fn largest_component_share(res: u32, plate_id: &CubeField<u16>, plate: u16) -> f64 {
        let g = grid::Grid::new(res);
        let id_of = |a: usize| plate_id.data[g.node_sample[a] as usize];
        let mut seen = vec![false; g.n];
        let (mut total, mut largest) = (0usize, 0usize);
        for start in 0..g.n {
            if id_of(start) != plate {
                continue;
            }
            total += 1;
            if seen[start] {
                continue;
            }
            seen[start] = true;
            let mut stack = vec![start];
            let mut size = 0;
            while let Some(a) = stack.pop() {
                size += 1;
                for &b in g.neighbors(a) {
                    let b = b as usize;
                    if !seen[b] && id_of(b) == plate {
                        seen[b] = true;
                        stack.push(b);
                    }
                }
            }
            largest = largest.max(size);
        }
        // `total` counts every node once; components cover the same nodes.
        largest as f64 / total.max(1) as f64
    }

    #[test]
    fn every_plate_is_connected() {
        let s = default_snapshot();
        for id in 0..s.plates.len() as u16 {
            let share = largest_component_share(s.plate_id.res, &s.plate_id, id);
            assert!(share > 0.95, "plate {id}: largest component {share}");
        }
        // And every id used is a real plate.
        assert!(s.plate_id.data.iter().all(|&i| (i as usize) < s.plates.len()));
    }

    #[test]
    fn all_boundary_classes_present() {
        let s = default_snapshot();
        for class in [boundary_class::CONVERGENT, boundary_class::DIVERGENT, boundary_class::TRANSFORM] {
            let n = s.boundary.data.iter().filter(|&&b| b == class).count();
            assert!(n > 20, "class {class}: only {n} samples");
        }
        // Both crust types present.
        assert!(s.crust.data.contains(&0) && s.crust.data.contains(&1));
    }

    #[test]
    fn boundary_classes_present_across_seeds() {
        for seed in 1..=5u64 {
            let s = generate(&TectonicsParams { seed, res: 96, ..Default::default() });
            for class in 1..=3u8 {
                assert!(s.boundary.data.contains(&class), "seed {seed} lacks class {class}");
            }
        }
    }

    #[test]
    fn plates_are_valid_rotations() {
        for pl in &default_snapshot().plates {
            assert!((pl.euler_pole.length() - 1.0).abs() < 1e-9);
            assert!(pl.omega > 0.0 && pl.omega < 0.05);
        }
    }

    #[test]
    fn face_edge_twins_agree() {
        // Samples on a shared cube edge are one place; they must hold identical values.
        let s = default_snapshot();
        let g = grid::Grid::new(s.height_m.res);
        let mut by_node = vec![f32::NAN; g.n];
        for (sample, &node) in g.sample_node.iter().enumerate() {
            let h = s.height_m.data[sample];
            let slot = &mut by_node[node as usize];
            assert!(slot.is_nan() || *slot == h);
            *slot = h;
        }
    }

    /// R3: every discrete coast decision comes from the fixed classification grid, so it
    /// can't differ between Normal and Ultra (`refined_world_keeps_the_canonical_shape`
    /// depends on it).
    #[test]
    fn classification_is_resolution_independent() {
        let a = coast::Classification::build(&TectonicsParams { res: 64, ..Default::default() });
        let b = coast::Classification::build(&TectonicsParams { res: 640, ..Default::default() });
        assert_eq!(a.coast_margin, b.coast_margin);
        let classes = |c: &coast::Classification| c.seaways.iter().map(|s| (s.a, s.b, s.class, s.active_side)).collect::<Vec<_>>();
        assert_eq!(classes(&a), classes(&b));
        // The active share lands near the target (whole seaways are taken, so not exact).
        assert!((a.active_share - 0.25).abs() < 0.1, "active share {}", a.active_share);
    }

    #[test]
    fn all_passive_has_no_active_margin() {
        let s = generate(&TectonicsParams { res: 64, active_share: 0.0, ..Default::default() });
        let m = s.margin.as_ref().unwrap();
        assert!(!m.data.contains(&coast::margin::ACTIVE));
        assert!(m.data.contains(&coast::margin::PASSIVE));
    }

    #[test]
    fn active_coasts_get_trenches_and_convergence() {
        let (s, pm, st) = generate_with_model(&TectonicsParams { res: 96, ..Default::default() });
        assert!(s.margin.as_ref().unwrap().data.contains(&coast::margin::ACTIVE));
        let sub = pm.bnodes.iter().filter(|b| b.kind == plates::Kind::Subduction).count();
        assert!(sub > 20, "only {sub} subduction nodes");
        // Every convergent boundary node carries its committed (clamped) rate.
        for b in &pm.bnodes {
            if matches!(b.kind, plates::Kind::Subduction | plates::Kind::Collision) {
                assert!(pm.convergence[b.node as usize] >= 10.0);
            }
        }
        assert!(st.rms_residual_mm_yr.is_finite() && st.wrong_sign_share < 0.4, "{st:?}");
    }

    /// R3.1b gate 1 (tectonic half): with strict coasts the committed height's sign is the
    /// mask's at every node, land ≥ 1 m, sea ≤ −1 m.
    #[test]
    fn strict_coasts_conform_to_the_mask() {
        let mut p = TectonicsParams { res: 64, ..Default::default() };
        p.conditioning.strict = true;
        p.conditioning.islands = false;
        let s = generate(&p);
        let g = grid::Grid::new(64);
        let value = crust::continentalness(&g, &p);
        let cover = mask::land_cover(&g, &p, &value);
        for (a, &c) in cover.iter().enumerate() {
            let h = s.height_m.data[g.node_sample[a] as usize];
            assert_eq!(c > 0.5, h > 0.0, "node {a}: cover {c} h {h}");
            assert!(h >= 1.0 || h <= -1.0, "node {a}: h {h}");
        }
        assert!((s.land_fraction - 0.29).abs() < 0.02, "land fraction {}", s.land_fraction);
        assert!(s.land_cover.is_some());
        // Strokes paint land into the sea and sea into land.
        let sea_dir = (0..g.n).find(|&a| cover[a] < 0.5 && value[a] < value.iter().copied().fold(f32::MAX, f32::min) + 0.3).unwrap();
        p.conditioning.strokes.push(mask::Stroke::new(vec![g.dir[sea_dir]], 600.0, mask::StrokeOp::Land));
        let s2 = generate(&p);
        assert!(s2.height_m.data[g.node_sample[sea_dir] as usize] > 0.0, "painted land stays land");
    }

    /// s07 heightmap gate at the tectonics level: inside a heightmap layer (pin 1) land takes
    /// the imported height exactly, a 6 km import isn't squashed by the soft peak limit, the
    /// coast signs follow the mask, and a flat-sea map leaves the sea floor to bathymetry.
    #[test]
    fn heightmap_layer_sets_land_heights() {
        use mask::{Elevation, Fidelity, Placement, RasterMask};
        let (w, h) = (64u32, 64u32);
        // Land in the middle columns, rising north → south from grey 0.4 to 1.0; sea grey 0.
        let data: Vec<f32> = (0..w * h).map(|k| if (16..48).contains(&(k % w)) { 0.4 + 0.6 * (k / w) as f32 / 63.0 } else { 0.0 }).collect();
        let e = Elevation { max_m: 6000.0, fidelity: Fidelity::Exact, ..Default::default() };
        let mut m = RasterMask {
            data: std::sync::Arc::new(data),
            w,
            h,
            place: Placement::Patch { center_deg: [0.0, 0.0], km_per_px: 60.0, rotation_deg: 0.0 },
            adjust: mask::MaskAdjust { threshold: 0.2, invert: false },
            decode: Default::default(),
            channel_used: mask::Channel::Brightness,
            enabled: true,
            name: String::new(),
            source: std::sync::Arc::new(vec![]),
            source_kind: "png".into(),
            step: 1,
            elevation: Some(e),
            heights: None,
        };
        m.build_heights();
        let mut p = TectonicsParams { res: 64, ..Default::default() };
        p.conditioning.strict = true;
        p.conditioning.outside_procedural = false;
        p.conditioning.islands = false;
        p.conditioning.masks = vec![std::sync::Arc::new(m.clone())];
        let s = generate(&p);
        let g = grid::Grid::new(64);
        let cover = mask::land_cover(&g, &p, &vec![0.0; g.n]);
        let pin = s.pin.as_ref().expect("pin weights with a heightmap layer");
        let spacing = g.spacing_rad() * p.radius_m;
        let (mut checked, mut top) = (0, 0.0f32);
        for (a, &c) in cover.iter().enumerate() {
            let smp = g.node_sample[a] as usize;
            let hh = s.height_m.data[smp];
            assert_eq!(c > 0.5, hh > 0.0, "coast sign at node {a}");
            if c > 0.5 && pin.data[smp] >= 1.0 {
                let want = m.height_m(g.dir[a], p.radius_m, spacing).unwrap_or(1.0).max(1.0);
                assert!((hh - want).abs() < 1e-2, "node {a}: {hh} vs import {want}");
                checked += 1;
                top = top.max(hh);
            }
        }
        assert!(checked > 20, "only {checked} pinned land nodes");
        assert!(top > 5000.0, "soft limit bypassed: top {top}");
        // Flat-sea map: no sea floor from the import, so sea stays below −1 m (bathymetry).
        assert!(s.height_m.data.iter().zip(&pin.data).all(|(&hh, _)| hh >= 1.0 || hh <= -1.0));
    }

    /// s06 (Matt: "reefs / barrier island chains" around every coast). With strict coasts and
    /// tectonic islands on, land off the mask may only be an island on oceanic crust (arcs,
    /// hotspots), never the continental shelf poking up at the island guard distance.
    #[test]
    fn strict_islands_only_on_oceanic_crust() {
        let mut p = TectonicsParams { res: 96, ..Default::default() };
        p.conditioning.strict = true;
        p.conditioning.islands = true;
        let s = generate(&p);
        let g = grid::Grid::new(96);
        let value = crust::continentalness(&g, &p);
        let cover = mask::land_cover(&g, &p, &value);
        let mut on_shelf = 0;
        for (a, &c) in cover.iter().enumerate() {
            let smp = g.node_sample[a] as usize;
            if c < 0.5 && s.height_m.data[smp] > 0.0 && s.crust.data[smp] != 0 {
                on_shelf += 1;
            }
        }
        assert_eq!(on_shelf, 0, "{on_shelf} off-mask land nodes on continental crust");
        // An erase (sea) stroke over an island leaves open water.
        let island = (0..g.n).find(|&a| cover[a] < 0.5 && s.height_m.data[g.node_sample[a] as usize] > 0.0);
        let a = island.expect("the default seed has an oceanic island to erase");
        p.conditioning.strokes.push(mask::Stroke::new(vec![g.dir[a]], 300.0, mask::StrokeOp::Sea));
        let s2 = generate(&p);
        assert!(s2.height_m.data[g.node_sample[a] as usize] < 0.0, "erased island stays sea");
    }

    #[test]
    fn scales_with_parameters() {
        let base = TectonicsParams { res: 96, ..Default::default() };
        let lo = generate(&TectonicsParams { mountain_scale: 0.4, ..base.clone() });
        let hi = generate(&TectonicsParams { mountain_scale: 1.6, ..base });
        assert!(hi.height_m.min_max().1 > lo.height_m.min_max().1);
    }
}
