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

mod crust;
pub mod grid;
mod plates;
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
    /// Continental dynamic topography to build during erosion (m), if generated in this run.
    /// Not stored in world files: only erosion uses it, and its result is committed.
    pub dynamic_topo_m: Option<CubeField<f32>>,
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
    let p = &TectonicsParams { res: p.res.max(8), ..p.clone() };
    let grid = grid::Grid::new(p.res);
    let crust = crust::build(&grid, p);
    let pm = plates::build(&grid, p, &crust);
    let composed = relief::compose(&grid, p, &crust, &pm);

    // Sea level: the offset at which `land_fraction` of the area is above water.
    let target = p.land_fraction.clamp(0.0, 1.0) as f64;
    let offset = util::solve_threshold(&composed.raw(), &grid.weight, target);
    let height_node = relief::finish(&grid, p, &crust, &composed, offset);
    let achieved = grid.fraction(|a| height_node[a] > 0.0) as f32;

    // Copy node values out to every cube sample (twins on face edges get identical values).
    let mut snap = TectonicSnapshot {
        height_m: CubeField::new(p.res),
        plate_id: CubeField::new(p.res),
        crust: CubeField::new(p.res),
        boundary: CubeField::new(p.res),
        plates: pm.plates.clone(),
        sea_level_offset_m: offset,
        land_fraction: achieved,
        // Static offset in the heights (see relief.rs); the rate mechanism is off.
        dynamic_topo_m: None,
    };
    for (s, &node) in grid.sample_node.iter().enumerate() {
        let a = node as usize;
        snap.height_m.data[s] = height_node[a];
        snap.plate_id.data[s] = pm.plate_of[a];
        snap.crust.data[s] = crust.cont[a] as u8;
        snap.boundary.data[s] = pm.boundary[a];
        if let Some(f) = snap.dynamic_topo_m.as_mut() {
            f.data[s] = composed.dyn_topo[a];
        }
    }
    snap
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

    #[test]
    fn scales_with_parameters() {
        let base = TectonicsParams { res: 96, ..Default::default() };
        let lo = generate(&TectonicsParams { mountain_scale: 0.4, ..base.clone() });
        let hi = generate(&TectonicsParams { mountain_scale: 1.6, ..base });
        assert!(hi.height_m.min_max().1 > lo.height_m.min_max().1);
    }
}
