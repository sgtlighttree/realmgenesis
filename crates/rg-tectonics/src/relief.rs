//! Step 3: turn crust + plate boundaries into heights (metres, before the sea-level shift).
//!
//! Every node's height is built from layers, all measured as *distance to the nearest
//! boundary of some kind* and shaped by smooth bell/easing profiles:
//!
//! 1. **Base.** Continental crust rises from a ~-130 m shelf edge to a few hundred metres
//!    inland. Oceanic crust sits on a continental shelf/slope near coasts and otherwise follows
//!    the seafloor age-depth law: the ocean gets deeper (like sqrt of age) away from the
//!    nearest spreading ridge (GDH1, Stein & Stein 1992). Age = distance / spreading rate.
//! 2. **Collision** (continent meets continent): a broad plateau with a high range along the
//!    suture, and a shallow foreland basin beyond (Himalaya / Tibet / Ganges).
//! 3. **Subduction**: on the diving plate a deep trench and a slight outer bulge; on the
//!    overriding plate either a coastal cordillera + high plateau (Andes) or a chain of
//!    volcanic island arcs with a back-arc basin (Japan, Aleutians).
//! 4. **Continental rift**: a lowered valley with raised shoulders (East Africa).
//! 5. **Old belts**: eroded, lower and broader ranges inside continents (Urals/Appalachians).
//! 6. **Detail**: fBm everywhere, plus sharp ridged noise where the tectonic uplift is high.
//!    Amplitude follows context: rough in mountains, smooth on cratons and abyssal plains.

use crate::Params;
use crate::crust::Crust;
use crate::grid::Grid;
use crate::plates::{Kind, PlateModel};
use crate::rng::Rng;
use glam::DVec3;
use crate::util::*;
use rg_noise::salt;
use std::f64::consts::PI;

pub const STREAM_RELIEF: u32 = 40;

/// GDH1 seafloor depth (metres, positive down) at age `t` Myr.
pub fn gdh1_depth(t: f64) -> f64 {
    if t < 20.0 { 2500.0 + 350.0 * t.sqrt() } else { 5651.0 - 2473.0 * (-0.0278 * t).exp() }
}

/// One volcanic seamount of a hotspot chain (Hawaii-style).
struct Seamount {
    dir: DVec3,
    height: f64,
    width: f64,
}

/// A few hotspots under oceanic plates. A hotspot stays put while its plate slides over it, so
/// it leaves a chain of volcanoes: youngest and tallest at the hotspot, older ones carried
/// away along the plate's motion (a rotation about its Euler pole), shrinking and sinking.
fn hotspot_chains(grid: &Grid, p: &Params, crust: &Crust, pm: &PlateModel) -> Vec<Seamount> {
    let mut rng = Rng::new(p.seed, 50);
    let chains = 4 + p.plate_count / 3;
    let mut out = Vec::new();
    for _ in 0..chains {
        // Pick an oceanic spot (a few tries; give up quietly otherwise).
        let mut node = rng.index(grid.n);
        for _ in 0..30 {
            if !crust.cont[node] && crust.edge_m[node] < -400_000.0 {
                break;
            }
            node = rng.index(grid.n);
        }
        if crust.cont[node] {
            continue;
        }
        let plate = &pm.plates[pm.plate_of[node] as usize];
        let k = plate.euler_pole;
        let v0 = grid.dir[node];
        let count = rng.range(6.0, 15.0) as usize;
        let power = rng.range(0.7, 1.2);
        let width = rng.range(70_000.0, 110_000.0);
        let dt_myr = 4.0;
        for i in 0..count {
            let ang = plate.omega * dt_myr * i as f64;
            let (s, c) = ang.sin_cos();
            let v = v0 * c + k.cross(v0) * s + k * (k.dot(v0) * (1.0 - c));
            let fade = (1.0 - i as f64 / count as f64).powf(1.3);
            let height = 5200.0 * power * fade * rng.range(0.75, 1.15);
            out.push(Seamount { dir: v.normalize(), height, width: width * (0.6 + 0.4 * fade) });
        }
    }
    out
}

/// Everything `compose` computes, kept in layers so `finish` can rebuild the sea floor near
/// the coast once the sea level (and therefore the coast) is known.
pub struct Composed {
    /// Height before the sea-level shift, metres. Good enough to decide land vs sea.
    pub raw: Vec<f32>,
    /// Oceanic basin height from the age-depth law (negative metres), for every node.
    pub basin: Vec<f32>,
    /// Sum of the tectonic features (trench, arc, cordillera, ...), metres.
    pub tect: Vec<f32>,
    /// Detail noise in its "land" flavour and its "sea floor" flavour, metres.
    pub det_land: Vec<f32>,
    pub det_sea: Vec<f32>,
    /// Slowly varying 0..1 noise used to vary shelf widths.
    pub margin_w: Vec<f32>,
}

pub fn compose(grid: &Grid, p: &Params, crust: &Crust, pm: &PlateModel) -> Composed {
    let n = grid.n;
    let r = p.radius_m;
    let mscale = p.mountain_scale as f64;
    let rough = p.roughness as f64;

    // --- per-node attribute tables for boundary nodes
    let mut strength = vec![0f32; n];
    let mut speed = vec![0f32; n];
    let mut sub_plate = vec![0u16; n];
    let mut over_cont = vec![false; n];
    let (mut cc, mut sub, mut ridge, mut rift) = (Vec::new(), Vec::new(), Vec::new(), Vec::new());
    for b in &pm.bnodes {
        let a = b.node as usize;
        strength[a] = b.strength;
        speed[a] = b.speed_mm_yr;
        sub_plate[a] = b.sub_plate;
        over_cont[a] = b.over_cont;
        match b.kind {
            Kind::Collision => cc.push(b.node),
            Kind::Subduction => sub.push(b.node),
            Kind::Ridge => ridge.push(b.node),
            Kind::Rift => rift.push(b.node),
            Kind::Transform => {}
        }
    }
    let m_to_rad = |m: f64| (m / r) as f32;
    let near_cc = grid.nearest_source(&cc, m_to_rad(1_900_000.0));
    let near_sub = grid.nearest_source(&sub, m_to_rad(1_300_000.0));
    let near_ridge = grid.nearest_source(&ridge, PI as f32);
    let near_rift = grid.nearest_source(&rift, m_to_rad(600_000.0));

    // --- smooth fields (coarse grid)
    let s_reg = salt(p.seed, STREAM_RELIEF);
    let regional = grid.coarse_field(64, |d| fbm(d, s_reg, 2.0, 3, 0.5) as f32);
    let old_belts = grid.coarse_field(64, |d| {
        let q = d + warp3(d, s_reg.wrapping_add(9), 1.6, 2) * 0.25;
        ridged(q, s_reg.wrapping_add(1), 2.3, 3, 0.5) as f32
    });
    let band_warp = grid.coarse_field(64, |d| fbm01(d, s_reg.wrapping_add(2), 5.0, 2) as f32);
    let swell = grid.coarse_field(64, |d| fbm(d, s_reg.wrapping_add(3), 2.5, 3, 0.5) as f32);
    let plateau_n = grid.coarse_field(64, |d| fbm(d, s_reg.wrapping_add(5), 3.2, 3, 0.5) as f32);
    let margin_w = grid.coarse_field(64, |d| fbm01(d, s_reg.wrapping_add(4), 3.0, 2) as f32);

    // --- fine noise
    let s_det = salt(p.seed, STREAM_RELIEF + 1);
    let s_rid = salt(p.seed, STREAM_RELIEF + 2);
    let s_abyss = salt(p.seed, STREAM_RELIEF + 3);
    let s_arc = salt(p.seed, STREAM_RELIEF + 4);
    let oct = ((grid.res as f64).log2().floor() as i32 - 3).clamp(3, 6) as u32;

    let half_pi = 0.0;
    let _ = half_pi;
    let seamounts = hotspot_chains(grid, p, crust, pm);
    let mut raw = vec![0f32; n];
    let mut basin_out = vec![0f32; n];

    let mut tect_out = vec![0f32; n];
    let mut det_land_out = vec![0f32; n];
    let mut det_sea_out = vec![0f32; n];
    // Oceanic basin depth from seafloor age, everywhere. Where two ridges' influence zones
    // meet, the nearest-ridge distance has a crease; a few smoothing passes soften it.
    for a in 0..n {
        let age = if near_ridge.src[a] != u32::MAX {
            let src = near_ridge.src[a] as usize;
            let half_rate_km_myr = (speed[src] as f64 * 0.5).clamp(14.0, 55.0);
            near_ridge.dist_rad[a] as f64 * r / 1000.0 / half_rate_km_myr
        } else {
            180.0
        };
        basin_out[a] = -(gdh1_depth(age.min(180.0)) as f32 + 450.0 * swell[a]);
    }
    let mut scratch = basin_out.clone();
    for _ in 0..4 {
        for a in 0..n {
            let nb = grid.neighbors(a);
            let sum: f32 = nb.iter().map(|&b| basin_out[b as usize]).sum();
            scratch[a] = 0.5 * basin_out[a] + 0.5 * sum / nb.len() as f32;
        }
        std::mem::swap(&mut basin_out, &mut scratch);
    }

    for a in 0..n {
        let d = grid.dir[a];
        let e = crust.edge_m[a] as f64;
        let cont = crust.cont[a];
        let reg = regional[a] as f64;
        let bw = 0.8 + 0.4 * band_warp[a] as f64;

        let basin = basin_out[a] as f64;

        // ---------- base
        let mut h = if cont {
            let ease = smoothstep(0.0, 300_000.0, e);
            let interior = 1200.0 + 450.0 * reg + 500.0 * smoothstep(300_000.0, 1_600_000.0, e);
            -130.0 + interior * ease
        } else {
            basin
        };

        // ---------- tectonic layers
        let mut tect = 0.0f64;
        let mut uplift = 0.0f64; // positive contributions (for the ridged-detail mask)
        let mut range_src = u32::MAX; // boundary node that built the dominant mountain range

        if near_cc.src[a] != u32::MAX {
            let src = near_cc.src[a] as usize;
            let dist = near_cc.dist_rad[a] as f64 * r * bw;
            let s = strength[src] as f64 * mscale;
            let plateau = 2000.0 * s * (1.0 - smoothstep(120_000.0, 750_000.0, dist));
            let range = 2500.0 * s * bell(dist, 140_000.0);
            let foreland = -300.0 * s * bell(dist - 650_000.0, 160_000.0);
            let w = 0.35 + 0.65 * smoothstep(-100_000.0, 100_000.0, e);
            tect += (plateau + range + foreland) * w;
            uplift += (plateau + range).max(0.0) * w;
            range_src = near_cc.src[a];
        }

        if near_sub.src[a] != u32::MAX {
            let src = near_sub.src[a] as usize;
            let dist = near_sub.dist_rad[a] as f64 * r * bw;
            let s = strength[src] as f64 * mscale.max(0.3);
            let diving = pm.plate_of[a] == sub_plate[src];
            if diving {
                let ocean_w = 1.0 - smoothstep(-50_000.0, 50_000.0, e);
                tect -= (2800.0 + 1300.0 * s) * bell(dist - 40_000.0, 75_000.0) * ocean_w;
                tect += 250.0 * bell(dist - 170_000.0, 70_000.0) * ocean_w;
            } else if over_cont[src] {
                let w = 0.4 + 0.6 * smoothstep(-80_000.0, 80_000.0, e);
                let cordillera = 2900.0 * s * bell(dist - 170_000.0, 100_000.0);
                let altiplano = 1100.0 * s * (1.0 - smoothstep(120_000.0, 650_000.0, dist));
                tect += (cordillera + altiplano) * w;
                uplift += (cordillera + altiplano) * w;
                range_src = near_sub.src[a];
            } else {
                // Island arc: discontinuous volcanic chain, strongest where the noise says so.
                let sd = grid.dir[src];
                let chain = 0.15 + 1.25 * fbm01(sd, s_arc, 16.0, 2);
                let arc = 4900.0 * s * chain * bell(dist - 160_000.0, 48_000.0);
                let back_arc = -500.0 * bell(dist - 380_000.0, 130_000.0);
                tect += arc + back_arc;
                uplift += arc;
            }
        }

        if near_rift.src[a] != u32::MAX && cont {
            let src = near_rift.src[a] as usize;
            let dist = near_rift.dist_rad[a] as f64 * r * bw;
            let s = strength[src] as f64;
            tect += -800.0 * s * bell(dist, 60_000.0) + 450.0 * s * bell(dist - 150_000.0, 70_000.0);
        }

        // Hotspot volcanoes (mostly on the sea floor).
        let mut hot = 0.0f64;
        for m in &seamounts {
            let c2 = (d - m.dir).length_squared();
            let reach = 3.0 * m.width / r;
            if c2 < reach * reach {
                hot += m.height * (-(c2.sqrt() * r / m.width).powi(2)).exp();
            }
        }
        tect += hot * if cont { 0.15 } else { 1.0 };

        // Old, eroded belts: wide and lower than active ranges.
        let belt = smoothstep(0.58, 0.9, old_belts[a] as f64) * smoothstep(0.0, 250_000.0, e);
        tect += 1300.0 * mscale * belt;
        let belt_uplift = 0.5 * belt * mscale;

        // ---------- detail
        let relief_w = ((uplift - 500.0) / 2200.0).clamp(0.0, 1.0).max(belt_uplift.min(0.7));
        let land_like = cont || uplift > 300.0;
        let (mut det_land, mut det_sea) = (0.0f64, 0.0f64);
        if land_like {
            let craton = smoothstep(500_000.0, 1_400_000.0, e) * (1.0 - relief_w);
            let soft = fbm(d, s_det, 12.0, oct, 0.55);
            det_land = rough * 850.0 * soft * (0.3 + 0.7 * (1.0 - craton));
            // Tablelands with steep edges, and low rolling uplands.
            let table = smoothstep(0.12, 0.3, plateau_n[a] as f64) * smoothstep(0.0, 300_000.0, e);
            det_land += 650.0 * table * (1.0 - 0.5 * craton);
            det_land += rough * 450.0 * (ridged(d, s_det.wrapping_add(77), 5.0, oct - 1, 0.6) - 0.35)
                * (1.0 - 0.6 * craton);
            if relief_w > 0.02 {
                // Stretch the noise lookup across the belt so ridges run along it.
                let q = if range_src != u32::MAX {
                    let sd = grid.dir[range_src as usize];
                    sd + (d - sd) * 4.0
                } else {
                    d
                };
                let rd = ridged(q, s_rid, 8.0, oct, 0.6);
                det_land += rough * relief_w * 2100.0 * (rd - 0.32);
            }
        } else {
            // Sea floor: abyssal hills, rougher near ridges and trenches.
            let near_ridge_w = if near_ridge.src[a] != u32::MAX {
                bell(near_ridge.dist_rad[a] as f64 * r, 600_000.0)
            } else {
                0.0
            };
            let amp = 80.0 + 380.0 * near_ridge_w;
            det_sea = rough * amp * fbm(d, s_abyss, 14.0, oct, 0.6);
        }

        tect_out[a] = tect as f32;
        det_land_out[a] = det_land as f32;
        det_sea_out[a] = det_sea as f32;
        raw[a] = (h + tect + det_land + det_sea) as f32;
        let _ = &mut h;
    }
    Composed { raw, basin: basin_out, tect: tect_out, det_land: det_land_out, det_sea: det_sea_out, margin_w }
}

/// Shift by the sea level and rebuild the sea floor near coasts.
///
/// Before the sea-level solve nobody knows where the coast is, so `compose` only sketches
/// the ocean. Now we know: shelves and slopes are measured from the real coastline. Shelves
/// are wide on passive margins and narrow on active ones; islands on oceanic crust get only
/// a small apron. Land is untouched, and submerged nodes stay submerged, so the land fraction
/// found by the sea-level solve is preserved exactly.
pub fn finish(grid: &Grid, p: &Params, crust: &Crust, c: &Composed, offset: f32) -> Vec<f32> {
    let n = grid.n;
    let r = p.radius_m;
    let h0: Vec<f32> = c.raw.iter().map(|h| h - offset).collect();
    let mut coast: Vec<u32> = Vec::new();
    for a in 0..n {
        if h0[a] > 0.0 && grid.neighbors(a).iter().any(|&b| h0[b as usize] <= 0.0) {
            coast.push(a as u32);
        }
    }
    let near = grid.nearest_source(&coast, PI as f32);
    let half_cell = 0.5 * grid.spacing_rad() * r;
    let mut out = vec![0f32; n];
    for a in 0..n {
        let mut h = h0[a] as f64;
        if h0[a] <= 0.0 {
            let dc = (near.dist_rad[a] as f64 * r - half_cell).max(0.0);
            let src = near.src[a];
            let cont_src = src != u32::MAX && crust.cont[src as usize];
            let act = crust.activity[a] as f64;
            let var = 0.7 + 0.6 * c.margin_w[a] as f64;
            let shelf_w = lerp(120_000.0, 30_000.0, act) * var * if cont_src { 1.0 } else { 0.25 };
            let slope_w = lerp(260_000.0, 100_000.0, act) * var * if cont_src { 1.0 } else { 0.5 };
            let basin_depth = (-(c.basin[a] as f64)).max(330.0);
            let depth = if dc < shelf_w {
                8.0 + 122.0 * (dc / shelf_w).powf(1.3)
            } else {
                let t = smoothstep(shelf_w, shelf_w + slope_w, dc);
                lerp(130.0, basin_depth, t.powf(0.9))
            };
            let det = if crust.cont[a] {
                0.3 * c.det_land[a] as f64
            } else {
                c.det_sea[a] as f64 + 0.5 * c.det_land[a] as f64
            };
            h = (-depth + c.tect[a] as f64 + det).min(-1.0);
        }
        // Soft limits so extremes stay Earth-like (peaks <~ 8.8 km, trenches >~ -11 km).
        if h > 3800.0 {
            h = 3800.0 + 5000.0 * ((h - 3800.0) / 5000.0).tanh();
        }
        if h < -7000.0 {
            h = -7000.0 - 4000.0 * ((-7000.0 - h) / 4000.0).tanh();
        }
        out[a] = h as f32;
    }
    out
}
