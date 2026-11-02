//! Step 1: where is continental crust and where is oceanic crust?
//!
//! This is deliberately **independent of the plates** (the key lesson of V1: if continents are
//! plate silhouettes they look like Voronoi polygons). Plates are grown later and are told
//! about the crust, never the other way round.
//!
//! The recipe:
//! 1. Scatter a few **continental cores**. Each is a short arc on the sphere (a "spine") that
//!    has a strength. Big cores make continents, weak short ones make microcontinents.
//!    Elongated spines give N-S or E-W continents instead of round blobs.
//! 2. Continentalness `C(x) = max over cores of (strength - distance_to_spine / R0)`, measured
//!    after warping the position with smooth noise (so shapes swirl), plus multi-scale fBm
//!    (coastline detail, bays, offshore islands and island arcs).
//! 3. Choose the threshold so continental crust covers `land_fraction * 1.4` of the sphere
//!    (on Earth continental crust, shelves included, is ~41% against ~29% land).

use crate::Params;
use crate::grid::Grid;
use crate::rng::Rng;
use crate::util::*;
use glam::DVec3;
use rg_noise::salt;
use std::f64::consts::PI;

pub const STREAM_CRUST: u32 = 10;

/// One continental core: a great-circle arc from `a` to `b` (angular length `len`).
#[derive(Clone, Copy, Debug)]
pub struct Core {
    a: DVec3,
    b: DVec3,
    axis: DVec3,
    len: f64,
    mid: DVec3,
    strength: f64,
}

impl Core {
    fn new(a: DVec3, toward: DVec3, len: f64, strength: f64) -> Core {
        let t = (toward - a * toward.dot(a)).normalize();
        let b = a * len.cos() + t * len.sin();
        let axis = a.cross(t).normalize(); // normal of the arc's plane
        let mid = (a * (len * 0.5).cos() + t * (len * 0.5).sin()).normalize();
        Core { a, b, axis, len, mid, strength }
    }

    /// Angular distance from unit vector `p` to this arc.
    fn distance(&self, p: DVec3) -> f64 {
        let s = p.dot(self.axis);
        let q = p - self.axis * s;
        let ql = q.length();
        if ql > 1e-9 {
            let q = q / ql;
            let t = self.a.cross(q).dot(self.axis).atan2(self.a.dot(q));
            if (0.0..=self.len).contains(&t) {
                return s.abs().asin();
            }
        }
        p.dot(self.a).clamp(-1.0, 1.0).acos().min(p.dot(self.b).clamp(-1.0, 1.0).acos())
    }
}

pub struct Crust {
    /// True for continental crust.
    pub cont: Vec<bool>,
    /// Signed distance to the continental/oceanic boundary in metres: positive inside
    /// continental crust, negative in the ocean.
    pub edge_m: Vec<f32>,
    /// How "active" the continental margin is, 0 (passive, Atlantic-like) to 1 (Pacific-rim).
    pub activity: Vec<f32>,
}

fn make_cores(p: &Params) -> Vec<Core> {
    let mut rng = Rng::new(p.seed, 1);
    let hint = p.continent_count_hint.max(1) as usize;
    // (strength range, spine length range, count)
    type Class = ((f64, f64), (f64, f64), usize);
    let classes: [Class; 3] = [
        ((0.95, 1.45), (0.15, 1.1), hint),             // continents
        ((0.55, 0.85), (0.0, 0.45), hint),             // medium landmasses
        ((0.25, 0.5), (0.0, 0.18), 2 * hint + 2),      // microcontinents
    ];
    let mut cores: Vec<Core> = Vec::new();
    for (strength, len, count) in classes {
        for _ in 0..count {
            // Best-candidate sampling: spread cores out, but not perfectly evenly.
            let mut best: Option<(f64, Core)> = None;
            for _ in 0..14 {
                let a = rng.unit_vec();
                let toward = rng.unit_vec();
                let c = Core::new(a, toward, rng.range(len.0, len.1), rng.range(strength.0, strength.1));
                let nearest = cores
                    .iter()
                    .map(|o| {
                        // Weighted: large existing cores repel more strongly.
                        let d = c.mid.dot(o.mid).clamp(-1.0, 1.0).acos();
                        d - 0.35 * o.strength - 0.5 * o.len
                    })
                    .fold(f64::MAX, f64::min);
                let score = if cores.is_empty() { 0.0 } else { nearest } + rng.range(0.0, 0.5);
                if best.as_ref().is_none_or(|b| score > b.0) {
                    best = Some((score, c));
                }
            }
            cores.push(best.unwrap().1);
        }
    }
    cores
}

pub fn build(grid: &Grid, p: &Params) -> Crust {
    let cores = make_cores(p);
    let s_warp = salt(p.seed, STREAM_CRUST);
    let s_big = salt(p.seed, STREAM_CRUST + 1);
    let s_fine = salt(p.seed, STREAM_CRUST + 2);
    let s_act = salt(p.seed, STREAM_CRUST + 3);
    const R0: f64 = 0.42; // radians of arc per unit of core strength

    // Coarse part: warped core distance field plus large-scale noise. Smooth, so a 64-cell
    // cube grid interpolated up is plenty.
    let value_coarse = grid.coarse_field(64, |d| {
        let w = warp3(d, s_warp, 1.1, 2) * 0.34 + warp3(d, s_warp.wrapping_add(7), 3.3, 2) * 0.10;
        let q = (d + w).normalize();
        let mut best = -2.0f64;
        for c in &cores {
            let reach = c.len * 0.5 + 1.6;
            if q.dot(c.mid) < reach.min(PI).cos() {
                continue;
            }
            let v = c.strength - c.distance(q) / R0;
            if v > best {
                best = v;
            }
        }
        (best + 0.42 * fbm(d, s_big, 1.4, 4, 0.55) + 0.33 * fbm(d, s_big.wrapping_add(50), 5.0, 3, 0.6)) as f32
    });
    // Fine part: coastline-scale detail, evaluated per node.
    let value: Vec<f32> = (0..grid.n)
        .map(|a| {
            let d = grid.dir[a];
            value_coarse[a] + (0.2 * fbm(d, s_fine, 14.0, 3, 0.6)) as f32
        })
        .collect();

    let cont_target = (p.land_fraction as f64 * 1.4).clamp(0.05, 0.9);
    let thr = solve_threshold(&value, &grid.weight, cont_target);
    let mut cont: Vec<bool> = value.iter().map(|&v| v > thr).collect();

    // Majority filter: removes single-node specks and pinholes.
    let before = cont.clone();
    for a in 0..grid.n {
        let nb = grid.neighbors(a);
        let same = nb.iter().filter(|&&b| before[b as usize] == before[a]).count();
        if same * 3 <= nb.len() {
            cont[a] = !before[a];
        }
    }

    // Signed distance to the crust boundary, metres.
    let mut edge_nodes: Vec<u32> = Vec::new();
    for a in 0..grid.n {
        if grid.neighbors(a).iter().any(|&b| cont[b as usize] != cont[a]) {
            edge_nodes.push(a as u32);
        }
    }
    let near = grid.nearest_source(&edge_nodes, PI as f32);
    let half_cell = (0.5 * grid.spacing_rad() * p.radius_m) as f32;
    let edge_m: Vec<f32> = (0..grid.n)
        .map(|a| {
            let d = near.dist_rad[a] * p.radius_m as f32 + half_cell;
            if cont[a] { d } else { -d }
        })
        .collect();

    // Margin activity: where plate boundaries are encouraged to hug the coast.
    let activity: Vec<f32> = grid.coarse_field(48, |d| {
        let x = 0.5 + 0.5 * 2.2 * fbm(d, s_act, 1.8, 3, 0.5);
        smoothstep(0.45, 0.85, x) as f32
    });

    Crust { cont, edge_m, activity }
}
