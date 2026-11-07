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

/// Standard deviation of residual ocean depth (depth minus the age–depth law), m: mantle
/// dynamic topography and swells reach ±1 km over ~1000 km wavelengths (Hoggard et al.
/// 2016, Nature Geoscience). Calibrated with `probe_bathymetry` against Earth's depth bands.
pub const DYNAMIC_TOPO_STD_M: f64 = 350.0;
/// Share of the ocean floor that is oceanic plateau (Harris et al. 2014: 184 plateaus,
/// 18.5 M km², 5.11% of the oceans), and how far plateaus stand above the floor (Ontong
/// Java, Kerguelen, Shatsky: tops at 1.5–3 km in 5 km deep basins).
pub const PLATEAU_AREA: f64 = 0.051;
/// Continental residual (dynamic) topography std, m: a little below the oceanic value, as
/// thick continental lithosphere damps it (cf. Hoggard et al. 2016; not a fitted value).
pub const CONT_DYNAMIC_TOPO_STD_M: f64 = 300.0;
pub const PLATEAU_LIFT_M: f64 = 2200.0;

fn mean_sd(v: &[f64]) -> (f64, f64) {
    let n = v.len().max(1) as f64;
    let m = v.iter().sum::<f64>() / n;
    (m, (v.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / n).sqrt())
}

/// Oldest seafloor, Myr (older crust has been subducted).
pub const SEAFLOOR_MAX_AGE_MYR: f64 = 180.0;

/// Area–age calibration of the seafloor. Distance-from-ridge ages depend on where the
/// procedural ridges happen to sit, so their distribution is arbitrary; Earth's is close to
/// triangular, area falling linearly with age to ~180 Myr (Parsons 1982; mean ~60 Myr, near
/// the ~62 Myr of modern age grids). Oceanic nodes are ranked by raw age (area-weighted)
/// and given the age at that rank: `F(t) = 1 − (1 − t/T)²`, so `t = T·(1 − √(1 − F))`.
/// Monotonic, so ridges, flanks and the oldest basins stay where tectonics put them.
/// Continental nodes keep their raw age (capped), only used near margins.
pub fn earth_seafloor_age(grid: &Grid, cont: &[bool], raw: &[f64]) -> Vec<f64> {
    let t_max = SEAFLOOR_MAX_AGE_MYR;
    let mut out: Vec<f64> = raw.iter().map(|&t| t.min(t_max)).collect();
    let mut ocean: Vec<usize> = (0..grid.n).filter(|&a| !cont[a]).collect();
    ocean.sort_by(|&x, &y| raw[x].total_cmp(&raw[y]).then(x.cmp(&y)));
    let total: f64 = ocean.iter().map(|&a| grid.weight[a] as f64).sum();
    let mut below = 0.0;
    for &a in &ocean {
        let w = grid.weight[a] as f64;
        let f = ((below + 0.5 * w) / total.max(1e-300)).clamp(0.0, 1.0);
        below += w;
        out[a] = t_max * (1.0 - (1.0 - f).sqrt());
    }
    out
}

/// Isostatic lift of the seafloor per metre of sediment: a sediment load displaces water but
/// also depresses the crust, so the water gets shallower by `(ρm − ρs)/(ρm − ρw)` ≈ 0.49 of
/// the sediment thickness (ρm 3300, ρs 2200 compacted mean, ρw 1030 kg/m³).
pub const SEDIMENT_LIFT: f64 = (3300.0 - 2200.0) / (3300.0 - 1030.0);

/// Total sediment thickness on oceanic crust, m, from crustal age (Myr) and distance from
/// the continent–ocean boundary (m, positive offshore). Fitted to GlobSed (Straume et al.
/// 2019, Table 3): **404 m** mean on oceanic crust > 200 km offshore (pelagic rain; zero at
/// the ridge), and **~2,050 m** mean within 200 km (the terrigenous wedge of the continental
/// rise: 7.6×10⁷ km³ over ~3.7×10⁷ km²).
/// - Pelagic: `S·(1 − exp(−t/50 Myr))`, saturating because old crust drifts under
///   unproductive gyres and below the carbonate compensation depth (old Pacific crust carries
///   only a few hundred metres). S = 680 m gives the 404 m mean over the triangular age
///   distribution. (Linear in age over-lifted the oldest floor: a cliff at 5.5 km, s01.)
/// - Wedge: `A·exp(−x/L)` with L = 100 km; its 0–200 km mean is `0.432·A`, so A = 4,700 m
///   at the boundary.
pub fn sediment_m(age_myr: f64, offshore_m: f64) -> f64 {
    let pelagic = 680.0 * (1.0 - (-age_myr / 50.0).exp());
    let wedge = 4700.0 * (-offshore_m.max(0.0) / 100_000.0).exp();
    pelagic + wedge
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
        // Random places, not node indices, so hotspots don't move with resolution (s03).
        let mut node = grid.nearest_node(rng.unit_vec());
        for _ in 0..30 {
            if !crust.cont[node] && crust.edge_m[node] < -400_000.0 {
                break;
            }
            node = grid.nearest_node(rng.unit_vec());
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
    /// Per node: `[tect, det_land, det_sea, raw]`, all metres.
    /// * `tect`: sum of the tectonic features (trench, arc, cordillera, ...).
    /// * `det_land`, `det_sea`: detail noise in its "land" flavour and its "sea floor" flavour.
    /// * `raw`: height before the sea-level shift. Good enough to decide land vs sea.
    pub parts: Vec<[f32; 4]>,
    /// Oceanic basin height from the age-depth law (negative metres), for every node.
    pub basin: Vec<f32>,
    /// Slowly varying 0..1 noise used to vary shelf widths.
    pub margin_w: Vec<f32>,
    /// Continental dynamic (residual) topography to build during erosion, m per node:
    /// broad sags and swells ~1000–3000 km across (Africa's basin-and-swell; Hoggard et al.
    /// 2016), from the same mantle field as the sea-floor residual.
    pub dyn_topo: Vec<f32>,
}

impl Composed {
    /// The raw (pre sea level) heights as one array.
    pub fn raw(&self) -> Vec<f32> {
        self.parts.iter().map(|p| p[3]).collect()
    }
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
    // Four independent searches; they run side by side on native targets.
    let ((near_cc, near_sub), (near_ridge, near_rift)) = join(
        || {
            join(
                || grid.nearest_source(&cc, m_to_rad(1_900_000.0)),
                || grid.nearest_source(&sub, m_to_rad(1_300_000.0)),
            )
        },
        || join(|| grid.nearest_source(&ridge, PI as f32), || grid.nearest_source(&rift, m_to_rad(600_000.0))),
    );

    // --- smooth fields (coarse grid)
    let s_reg = salt(p.seed, STREAM_RELIEF);
    let regional = grid.coarse(64, |d| fbm(d, s_reg, 2.0, 3, 0.5) as f32);
    let old_belts = grid.coarse(64, |d| {
        let q = d + warp3(d, s_reg.wrapping_add(9), 1.6, 2) * 0.25;
        ridged(q, s_reg.wrapping_add(1), 2.3, 3, 0.5) as f32
    });
    let band_warp = grid.coarse(64, |d| fbm01(d, s_reg.wrapping_add(2), 5.0, 2) as f32);
    let swell = grid.coarse(64, |d| fbm(d, s_reg.wrapping_add(3), 2.5, 3, 0.5) as f32);
    let plateau_n = grid.coarse(64, |d| fbm(d, s_reg.wrapping_add(5), 3.2, 3, 0.5) as f32);
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

    // Oceanic basin depth from seafloor age, everywhere. Where two ridges' influence zones
    // meet, the nearest-ridge distance has a crease; a few smoothing passes soften it.
    // 1. Raw age: distance from the nearest ridge over its half spreading rate.
    let raw_age: Vec<f64> = par_map(n, |a| {
        if near_ridge.src[a] != u32::MAX {
            let src = near_ridge.src[a] as usize;
            let half_rate_km_myr = (speed[src] as f64 * 0.5).clamp(14.0, 55.0);
            near_ridge.dist_rad[a] as f64 * r / 1000.0 / half_rate_km_myr
        } else {
            f64::MAX
        }
    });
    // 2. Calibrate to Earth's area–age distribution (see `earth_seafloor_age`).
    let age = earth_seafloor_age(grid, &crust.cont, &raw_age);
    // 3. Depth: GDH1 subsidence, minus the isostatic lift of the sediment pile.
    // 4. Residual depth (mantle dynamic topography): the swell noise standardised to
    //    DYNAMIC_TOPO_STD_M over the ocean. 5. Oceanic plateaus (LIPs) on the top
    //    PLATEAU_AREA of an independent noise field.
    let ocean_nodes: Vec<usize> = (0..n).filter(|&a| !crust.cont[a]).collect();
    let swell_v: Vec<f64> = ocean_nodes.iter().map(|&a| grid.lookup(&swell, a) as f64).collect();
    let (sw_mean, sw_sd) = mean_sd(&swell_v);
    let lip = grid.coarse(64, |d| fbm(d, s_reg.wrapping_add(11), 4.0, 3, 0.5) as f32);
    let mut lip_v: Vec<f32> = ocean_nodes.iter().map(|&a| grid.lookup(&lip, a)).collect();
    lip_v.sort_by(f32::total_cmp);
    let lip_lo = lip_v.get(((lip_v.len() as f64) * (1.0 - PLATEAU_AREA)) as usize).copied().unwrap_or(f32::MAX) as f64;
    let lip_hi = lip_v.last().copied().unwrap_or(f32::MAX) as f64;
    let mut basin_out: Vec<f32> = par_map(n, |a| {
        let t = age[a];
        let residual = DYNAMIC_TOPO_STD_M * (grid.lookup(&swell, a) as f64 - sw_mean) / sw_sd.max(1e-9);
        // Flat top, abrupt sides: full lift within the top ~85% of the band above threshold.
        let lift = PLATEAU_LIFT_M * smoothstep(lip_lo, lip_lo + 0.15 * (lip_hi - lip_lo), grid.lookup(&lip, a) as f64);
        let water = gdh1_depth(t) - residual - lift;
        let sed = sediment_m(t, -(crust.edge_m[a] as f64));
        -((water - SEDIMENT_LIFT * sed) as f32)
    });
    drop(speed);
    for _ in 0..4 {
        basin_out = par_map(n, |a| {
            let nb = grid.neighbors(a);
            let sum: f32 = nb.iter().map(|&b| basin_out[b as usize]).sum();
            0.5 * basin_out[a] + 0.5 * sum / nb.len() as f32
        });
    }

    let dyn_topo: Vec<f32> = (0..n).map(|a| (CONT_DYNAMIC_TOPO_STD_M * (grid.lookup(&swell, a) as f64 - sw_mean) / sw_sd.max(1e-9)) as f32).collect();

    // The per-node work below depends only on read-only tables, so it runs in parallel.
    let parts: Vec<[f32; 4]> = par_map(n, |a| {
        let d = grid.dir[a];
        let e = crust.edge_m[a] as f64;
        let cont = crust.cont[a];
        let reg = grid.lookup(&regional, a) as f64;
        let bw = 0.8 + 0.4 * grid.lookup(&band_warp, a) as f64;

        let basin = basin_out[a] as f64;

        // ---------- base
        let h = if cont {
            let ease = smoothstep(0.0, 300_000.0, e);
            let interior = 1200.0 + 450.0 * reg + 500.0 * smoothstep(300_000.0, 1_600_000.0, e);
            // Continental dynamic topography: the same mantle-driven residual as on the sea
            // floor (one field across the margin), giving broad interior sags and swells
            // ~1000–3000 km across (Africa's basin-and-swell; Hoggard et al. 2016). Faded in
            // from the coast so shorelines stay put. Applied as a static offset: as an uplift
            // *rate* it raised endorheic land a little (2.3 → 3.6%) but pushed hypsometry off
            // Earth's (> 4 km 2.0 → 3.9%); hypsometry is the resolution-robust check (s01).
            match &crust.coast_m {
                // Strict (R3.1b): ramp up from the mask coast, so the coast sits at ~0 by
                // construction (the crust edge is 100–250 km offshore there).
                Some(cm) if cm[a] > 0.0 => 1.0 + (interior + dyn_topo[a] as f64) * smoothstep(0.0, 300_000.0, cm[a] as f64),
                Some(_) => -130.0,
                None => -130.0 + (interior + dyn_topo[a] as f64) * ease,
            }
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
                // Trench: a narrow flexural V (Earth: ~50–100 km wide, > 7 km only along the
                // deepest tens of km; > 7 km is 0.13% of the ocean floor). A 75 km Gaussian
                // gave 6–8× too much > 7 km floor (s01 probe_bathymetry).
                tect -= (2800.0 + 1300.0 * s) * bell(dist - 30_000.0, 40_000.0) * ocean_w;
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
        let belt = smoothstep(0.58, 0.9, grid.lookup(&old_belts, a) as f64) * smoothstep(0.0, 250_000.0, e);
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
            let table = smoothstep(0.12, 0.3, grid.lookup(&plateau_n, a) as f64) * smoothstep(0.0, 300_000.0, e);
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

        [tect as f32, det_land as f32, det_sea as f32, (h + tect + det_land + det_sea) as f32]
    });
    Composed { parts, basin: basin_out, margin_w, dyn_topo }
}

/// Shift by the sea level and rebuild the sea floor near coasts.
///
/// Before the sea-level solve nobody knows where the coast is, so `compose` only sketches
/// the ocean. Now we know: shelves and slopes are measured from the real coastline. Shelves
/// are wide on passive margins and narrow on active ones; islands on oceanic crust get only
/// a small apron. Land is untouched, and submerged nodes stay submerged, so the land fraction
/// found by the sea-level solve is preserved exactly.
///
/// **Strict** (R3.1b): the land/sea split is the mask's (`crust.land`). Land nodes are lifted
/// to ≥ 1 m; sea nodes are pushed to ≤ −1 m, and sea heights above −200 m are compressed
/// softly into [−200, −1] m, so an arc or hotspot crest becomes a shallow bank rather than
/// a flat sheet. With `islands`, tectonic islands on oceanic crust may stand above sea level
/// farther than `islands_guard_km` from the drawn coast.
pub fn finish(grid: &Grid, p: &Params, crust: &Crust, c: &Composed, offset: f32, strict: bool) -> Vec<f32> {
    let n = grid.n;
    let r = p.radius_m;
    let cond = &p.conditioning;
    let guard_m = (cond.islands_guard_km * 1000.0) as f32;
    // Islands only on oceanic crust (arcs, hotspots). Strict heights aren't shifted to a
    // solved sea level, so the continental shelf is still above zero past the guard and
    // used to stand up as a ring of "barrier islands" around every coast (s06, Matt).
    // Sea ("erase") strokes clear islands too: painted water stays open water.
    let erased = |a: usize| {
        cond.strokes.iter().any(|st| st.op == crate::mask::StrokeOp::Sea && st.weight(grid.dir[a], r) > 0.0)
    };
    let island = |a: usize, h: f32| -> bool {
        cond.islands
            && h > 0.0
            && !crust.cont[a]
            && crust.coast_m.as_ref().is_some_and(|cm| cm[a] < -guard_m)
            && !erased(a)
    };
    let h0: Vec<f32> = c
        .parts
        .iter()
        .enumerate()
        .map(|(a, p)| {
            let h = p[3] - offset;
            if !strict {
                h
            } else if crust.land[a] {
                h.max(1.0)
            } else if island(a, h) {
                h
            } else {
                h.min(-1.0)
            }
        })
        .collect();
    let mut coast: Vec<u32> = Vec::new();
    for a in 0..n {
        if h0[a] > 0.0 && grid.neighbors(a).iter().any(|&b| h0[b as usize] <= 0.0) {
            coast.push(a as u32);
        }
    }
    let near = grid.nearest_source(&coast, PI as f32);
    let half_cell = 0.5 * grid.spacing_rad() * r;
    par_map(n, |a| {
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
                0.3 * c.parts[a][1] as f64
            } else {
                c.parts[a][2] as f64 + 0.5 * c.parts[a][1] as f64
            };
            h = (-depth + c.parts[a][0] as f64 + det).min(-1.0);
        }
        if strict && h0[a] <= 0.0 && h > -200.0 {
            h = -1.0 - 199.0 * (-(h + 200.0) / 199.0).exp();
        }
        // Soft limits so extremes stay Earth-like (peaks <~ 8.8 km, trenches >~ -11 km).
        if h > 3800.0 {
            h = 3800.0 + 5000.0 * ((h - 3800.0) / 5000.0).tanh();
        }
        if h < -7000.0 {
            h = -7000.0 - 4000.0 * ((-7000.0 - h) / 4000.0).tanh();
        }
        if strict && h0[a] > 0.0 {
            h = h.max(1.0);
        }
        h as f32
    })
}
