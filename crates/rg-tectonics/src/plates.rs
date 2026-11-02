//! Step 2: plates, their motion, and the boundaries between them.
//!
//! * **Shapes.** Each plate grows outward from one seed node over the sample graph
//!   (multi-source Dijkstra). Stepping costs more through noisy "rough" terrain, so the borders
//!   come out irregular, and each plate has its own growth speed, so sizes vary a lot (a few
//!   huge plates, several small ones). Because plates grow as connected blobs, a plate can
//!   never have an exclave.
//! * **Coupling to crust.** Along "active" stretches of coast the growth cost jumps at the
//!   continental margin, so fronts stall there and plate borders tend to hug the coast (Andes,
//!   Japan). Along "passive" stretches (Atlantic coasts) there is no such effect, and borders
//!   fall wherever the fronts meet (mid-ocean). Continents are never whole plates, because the
//!   crust field is made first and independently.
//! * **Motion.** Each plate rotates about an Euler pole: the velocity at point `p` is
//!   `omega * (pole x p)`.
//! * **Boundary type.** At each boundary node we estimate the boundary's normal from the local
//!   plate layout, then split the *relative* velocity into a part across the boundary
//!   (convergent / divergent) and a part along it (transform).

use crate::crust::Crust;
use crate::grid::Grid;
use crate::rng::Rng;
use crate::util::*;
use crate::{Params, Plate};
use glam::DVec3;
use rg_noise::salt;

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
}

pub fn velocity(plate: &Plate, p: DVec3) -> DVec3 {
    plate.euler_pole.cross(p) * plate.omega
}

pub fn build(grid: &Grid, p: &Params, crust: &Crust) -> PlateModel {
    let mut rng = Rng::new(p.seed, STREAM_PLATES as u64);
    let k = p.plate_count.clamp(2, 60) as usize;

    // --- seeds: spread out, biased a little toward continents so each continent tends to sit
    // on (or straddle) its own plates. Some plates get extra seeds strung out along their long
    // axis, which makes them elongated and non-convex.
    let mut seeds: Vec<u32> = Vec::new();
    let mut seed_region: Vec<u16> = Vec::new();
    let mut axis: Vec<DVec3> = Vec::new();
    let mut stretch: Vec<f32> = Vec::new();
    for region in 0..k {
        let mut best: Option<(f64, u32)> = None;
        for _ in 0..40 {
            let c = rng.index(grid.n);
            let nearest = seeds
                .iter()
                .map(|&s| grid.dir[c].dot(grid.dir[s as usize]).clamp(-1.0, 1.0).acos())
                .fold(3.2f64, f64::min);
            let bias = if crust.cont[c] { 1.25 } else { 1.0 };
            let score = nearest * bias * rng.range(0.8, 1.2);
            if best.is_none_or(|b| score > b.0) {
                best = Some((score, c as u32));
            }
        }
        let first = best.unwrap().1;
        seeds.push(first);
        seed_region.push(region as u16);
        let sd = grid.dir[first as usize];
        let v = rng.unit_vec();
        let g = (v - sd * v.dot(sd)).normalize(); // pole of the plate's long great circle
        axis.push(g);
        stretch.push(rng.range(0.3, 1.5) as f32);
        let t = g.cross(sd); // direction along the great circle
        let extra = match rng.f64() {
            x if x < 0.35 => 0,
            x if x < 0.8 => 1,
            _ => 2,
        };
        for _ in 0..extra {
            let ang = rng.range(0.25, 0.6) * if rng.f64() < 0.5 { -1.0 } else { 1.0 };
            let target = sd * ang.cos() + t * ang.sin();
            let node = (0..grid.n)
                .max_by(|&x, &y| grid.dir[x].dot(target).total_cmp(&grid.dir[y].dot(target)))
                .unwrap();
            seeds.push(node as u32);
            seed_region.push(region as u16);
        }
    }
    let speed: Vec<f32> = (0..k).map(|_| (0.3 * rng.gauss()).exp().clamp(0.7, 1.6) as f32).collect();

    // --- growth cost: rough noise plus the continental-margin wall.
    let s_cost = salt(p.seed, STREAM_PLATES);
    let rough = grid.coarse_field(96, |d| {
        let q = d + warp3(d, s_cost.wrapping_add(5), 1.5, 2) * 0.25;
        (0.5 + 0.5 * 1.7 * fbm(q, s_cost, 2.4, 4, 0.72)).clamp(0.0, 1.0) as f32
    });
    let rough_fine = grid.coarse_field(96, |d| {
        (0.5 + 0.5 * 1.7 * fbm(d, s_cost.wrapping_add(8), 9.0, 3, 0.75)).clamp(0.0, 1.0) as f32
    });
    let cost: Vec<f32> = (0..grid.n)
        .map(|a| {
            let wall = 0.55 * crust.activity[a] as f64 * bell(crust.edge_m[a] as f64, 140_000.0);
            0.35 + 2.6 * rough[a] * rough[a] + 1.6 * rough_fine[a] * rough_fine[a] + wall as f32
        })
        .collect();
    let plate_of = grid.grow_regions(&seeds, &seed_region, &speed, &axis, &stretch, &cost);

    // --- plate properties
    let mut cont_area = vec![0f64; k];
    let mut area = vec![0f64; k];
    for (a, &plate) in plate_of.iter().enumerate() {
        let pl = plate as usize;
        area[pl] += grid.weight[a] as f64;
        if crust.cont[a] {
            cont_area[pl] += grid.weight[a] as f64;
        }
    }
    let mut density = Vec::with_capacity(k);
    let mut plates = Vec::with_capacity(k);
    for pl in 0..k {
        let continental = cont_area[pl] > 0.5 * area[pl];
        // Oceanic plates move faster than continental ones (slab pull). 0.15 to 0.9 degrees/Myr.
        let deg_per_myr = rng.range(0.15, 0.7) * if continental { 1.0 } else { 1.35 };
        plates.push(Plate {
            euler_pole: rng.unit_vec(),
            omega: deg_per_myr.to_radians(),
            continental,
        });
        density.push(rng.f64() as f32);
    }

    // --- boundaries
    let n = grid.n;
    let mut is_border = vec![false; n];
    for a in 0..n {
        is_border[a] = grid.neighbors(a).iter().any(|&b| plate_of[b as usize] != plate_of[a]);
    }
    let s_strength = salt(p.seed, STREAM_PLATES + 1);
    let mut boundary = vec![0u8; n];
    let mut bnodes: Vec<BoundaryNode> = Vec::new();
    let mut stamp = vec![0u32; n];
    let mut stamp_now = 0u32;
    let mm_per_yr = p.radius_m * 1e-3; // rad/Myr -> mm/yr
    for a in 0..n {
        if !is_border[a] {
            continue;
        }
        // Gather everything within 3 hops.
        stamp_now += 1;
        stamp[a] = stamp_now;
        let mut ball = vec![a as u32];
        let mut frontier = vec![a as u32];
        for _ in 0..3 {
            let mut next = Vec::new();
            for &x in &frontier {
                for &y in grid.neighbors(x as usize) {
                    if stamp[y as usize] != stamp_now {
                        stamp[y as usize] = stamp_now;
                        ball.push(y);
                        next.push(y);
                    }
                }
            }
            frontier = next;
        }
        let own = plate_of[a];
        // The most common other plate in the ball.
        let mut counts: Vec<(u16, u32)> = Vec::new();
        for &x in &ball {
            let pl = plate_of[x as usize];
            if pl != own {
                match counts.iter_mut().find(|c| c.0 == pl) {
                    Some(c) => c.1 += 1,
                    None => counts.push((pl, 1)),
                }
            }
        }
        counts.sort_by_key(|&(pl, c)| (std::cmp::Reverse(c), pl));
        let Some(&(other, _)) = counts.first() else { continue };

        let (mut sa, mut ca, mut cont_a) = (DVec3::ZERO, 0u32, 0u32);
        let (mut sb, mut cb, mut cont_b) = (DVec3::ZERO, 0u32, 0u32);
        for &x in &ball {
            let pl = plate_of[x as usize];
            if pl == own {
                sa += grid.dir[x as usize];
                ca += 1;
                cont_a += crust.cont[x as usize] as u32;
            } else if pl == other {
                sb += grid.dir[x as usize];
                cb += 1;
                cont_b += crust.cont[x as usize] as u32;
            }
        }
        let pa = grid.dir[a];
        let toward = sb / cb as f64 - sa / ca as f64;
        let n_ab = toward - pa * toward.dot(pa);
        if n_ab.length() < 1e-9 {
            continue;
        }
        let n_ab = n_ab.normalize();
        let a_cont = cont_a * 2 >= ca;
        let b_cont = cont_b * 2 >= cb;

        let rel = velocity(&plates[own as usize], pa) - velocity(&plates[other as usize], pa);
        let speed_mm = rel.length() * mm_per_yr;
        let across = rel.dot(n_ab); // > 0: converging
        let ratio = if rel.length() > 0.0 { across / rel.length() } else { 0.0 };
        let rate = (across.abs() * mm_per_yr) as f32;
        let along_strike = 0.55 + 0.9 * fbm01(pa, s_strength, 4.0, 3);

        let (class, kind) = if speed_mm < 3.0 || ratio.abs() < 0.4 {
            (3u8, Kind::Transform)
        } else if ratio > 0.0 {
            (1u8, if a_cont && b_cont { Kind::Collision } else { Kind::Subduction })
        } else {
            (2u8, if a_cont && b_cont { Kind::Rift } else { Kind::Ridge })
        };
        boundary[a] = class;

        let rate_factor = match kind {
            Kind::Ridge | Kind::Rift => (rate as f64 / 2.0 / 40.0).clamp(0.3, 1.2),
            _ => (rate as f64 / 55.0).clamp(0.3, 1.4),
        };
        let strength = (rate_factor * along_strike) as f32;

        let mut sub_plate = own;
        let mut over_cont = false;
        if kind == Kind::Subduction {
            // The oceanic plate dives; if both are oceanic, the denser one dives.
            let own_dives = if a_cont != b_cont {
                !a_cont
            } else {
                density[own as usize] > density[other as usize]
            };
            sub_plate = if own_dives { own } else { other };
            over_cont = if own_dives { b_cont } else { a_cont };
        }
        bnodes.push(BoundaryNode { node: a as u32, kind, strength, speed_mm_yr: speed_mm as f32, sub_plate, over_cont });
    }

    PlateModel { plate_of, plates, boundary, bnodes }
}
