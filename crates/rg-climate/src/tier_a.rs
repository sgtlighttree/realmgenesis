//! Tier A climate: a fast, rule-based but physically motivated monthly climate.
//!
//! The pipeline, in the order it runs (see each function for the details):
//!
//! 1. **Zonal energy balance** ([`crate::ebm`]): one seasonal 1-D EBM with an ocean-like heat
//!    capacity and one with a land-like heat capacity give sea-level temperature for every
//!    latitude and month.
//! 2. **Geography**: land mask, distance to the coast, and the orientation of coasts along each
//!    row (for ocean currents).
//! 3. **Per month**:
//!    - provisional temperature (continentality from plain coast distance + ocean currents);
//!    - the *thermal equator* of each longitude, which positions the rain belt (ITCZ) and is
//!      pulled toward hot land (monsoons);
//!    - prevailing winds: three-cell belts shifted with the season, plus a "thermal low" wind
//!      that blows toward hot land and away from cold land;
//!    - an **upwind march** from every cell: it finds how far inland the air has travelled
//!      (for continentality and moisture) and removes moisture on windward slopes, leaving
//!      rain shadows behind mountains;
//!    - final temperature (continentality including the upwind fetch, lapse rate) and
//!      precipitation (circulation rain belts x moisture, plus orographic, convective and
//!      monsoon rain).
//! 4. **Köppen** classes for land cells.
//!
//! Everything is deterministic. All tuning constants are collected below with notes.

use crate::circulation::{Cells, RowBelts, band_wind, piecewise_cos};
use crate::ebm::{HEAT_CAPACITY_LAND, HEAT_CAPACITY_OCEAN, solve_zonal_ebm};
use crate::geography::{Geometry, box_blur, distance_to_ocean_km, land_mask, row_runs_km};
use crate::grid::LatLonGrid;
use crate::insolation::{declination, month_mid_day};
use crate::koppen::Koppen;
use crate::types::{ClimateFields, ClimateInput, PlanetParams};

// ---------------------------------------------------------------------------------------------
// Tuning constants
// ---------------------------------------------------------------------------------------------

/// Standard lapse rate, C per km of land elevation.
const LAPSE_C_PER_KM: f32 = 6.5;

/// Land seasonal swing relative to the thin-land EBM. The EBM has no marine air advection or
/// evaporative cooling, so its raw continental swing is about 1.6x too large; this is the
/// empirical correction (checked against Moscow, Yakutsk, Beijing, Chicago normals).
const LAND_AMPLITUDE_GAIN: f32 = 0.62;
/// Share of the land-minus-ocean annual-mean difference of the two EBM runs that is kept. The
/// thin-land EBM is far too cold on average (it has no ocean heat transport reaching inland), so
/// only a quarter of it is used.
const LAND_MEAN_SHARE: f32 = 0.25;
/// Extra warmth of land over the zonal-mean EBM, C (see [`Model::sea_level_temps`]).
const LAND_WARMING_C: f32 = 3.0;
/// How fast a land cell becomes "continental" with distance from the sea: `1 - exp(-d / L)`.
const CONTINENTALITY_KM: f32 = 1500.0;
/// Weight of the upwind sea distance versus the plain nearest-coast distance. Air arrives from
/// upwind, so a coast to the west matters more than a coast to the east in the westerlies.
const UPWIND_WEIGHT: f32 = 0.65;

/// Length scale (km) over which ocean-boundary-current effects fade away from a coast.
const CURRENT_SCALE_OCEAN_KM: f32 = 1500.0;
const CURRENT_SCALE_LAND_KM: f32 = 650.0;

/// Seasonal migration of the belts, as a fraction of the (lagged) solar declination.
const BELT_SHIFT_PER_DECLINATION: f64 = 0.30;
/// Lag between the sun and the thermal/rain-belt response, days.
const SEASON_LAG_DAYS: f64 = 30.0;

/// Strength of the thermal-low wind: m/s per (K per 1000 km) of smoothed temperature anomaly.
const MONSOON_WIND_GAIN: f32 = 0.9;
const MONSOON_WIND_MAX: f32 = 6.0;
/// Poleward/equatorward flow (m/s) around the subtropical ocean highs, see [`Model::winds`].
const GYRE_WIND_MS: f32 = 3.0;
/// Coastal moisture top-up: share of the gap to sea-air saturation filled, and its reach (km).
const COASTAL_MOISTURE_SHARE: f32 = 0.75;
const COASTAL_MOISTURE_KM: f32 = 900.0;
/// Peak monsoon rain (mm/month), and how far inland (km of upwind land) it reaches.
const MONSOON_RAIN_MM: f32 = 300.0;
const MONSOON_FETCH_KM: f32 = 1500.0;
/// Moisture the air carries when it begins its march over land (fraction of saturation).
const MOISTURE_START_OVER_LAND: f32 = 0.3;
/// Moisture the land itself evaporates back into passing air (soil, lakes, forests):
/// the air over land relaxes toward this value instead of drying out completely.
const LAND_EVAPORATION_FLOOR: f32 = 0.30;
/// Fraction of the gap to saturation an air parcel regains per grid step over ocean.
const OCEAN_RECHARGE: f32 = 0.30;
/// Rise (m) at which half of the moisture is wrung out by a single step of orographic lift.
const OROGRAPHIC_HALF_RISE_M: f32 = 700.0;
/// Rain (mm/month) delivered when a fully moist airstream is lifted completely.
const OROGRAPHIC_RAIN_MM: f32 = 280.0;
/// Convective summer rain: mm per month per degree above 15 C (8 C near the poles), for a fully moist parcel.
const CONVECTIVE_MM_PER_K: f32 = 7.0;

// ---------------------------------------------------------------------------------------------
// Public entry point
// ---------------------------------------------------------------------------------------------

/// Run Tier A. The output grid has the same size as `input.height_m`.
pub fn tier_a(planet: &PlanetParams, input: &ClimateInput) -> ClimateFields {
    let model = Model::new(planet, input);
    model.run()
}

// ---------------------------------------------------------------------------------------------
// Model state
// ---------------------------------------------------------------------------------------------

/// Everything that does not change from month to month.
struct Model<'a> {
    planet: &'a PlanetParams,
    geo: Geometry,
    cells: Cells,
    land: Vec<bool>,
    /// Land elevation in metres (0 over ocean).
    elev: Vec<f32>,
    /// Km to the nearest ocean cell (0 over ocean).
    dist_coast: Vec<f32>,
    /// "How much this cell sits on the EAST side of its ocean basin" (cold-current side) and the
    /// west side (warm-current side), each 0..1; see [`Model::current_anomaly`].
    east_side: Vec<f32>,
    west_side: Vec<f32>,
    /// Boundary-current anomaly (C) of each cell: `cur_base + winter * cur_winter`, where
    /// `winter` is +1 in deep winter and -1 in deep summer. See [`Model::current_bands`].
    cur_base: Vec<f32>,
    cur_winter: Vec<f32>,
    /// Per-row monthly zonal sea-level temperatures from the two EBM runs, plus their annual means.
    ocean_row: Vec<[f32; 12]>,
    land_row: Vec<[f32; 12]>,
    ocean_mean: Vec<f32>,
    land_mean: Vec<f32>,
}

/// Scratch arrays for one month.
struct Month {
    /// Sea-level temperature, provisional, then reused for the final value without lapse.
    t_sea: Vec<f32>,
    u: Vec<f32>,
    v: Vec<f32>,
    /// Moisture (0..1) in the air arriving at each cell.
    q: Vec<f32>,
    /// Orographic rain at each cell, mm.
    oro: Vec<f32>,
    /// Km to the nearest upwind ocean.
    fetch_km: Vec<f32>,
    /// Ocean moisture capacity (0 over land).
    qcap: Vec<f32>,
    /// Where the rain belt sits at each longitude, degrees.
    itcz_deg: Vec<f32>,
}

impl<'a> Model<'a> {
    fn new(planet: &'a PlanetParams, input: &ClimateInput) -> Self {
        let (w, h) = (input.height_m.width, input.height_m.height);
        assert!(w >= 8 && h >= 4, "climate grid too small");
        let geo = Geometry::new(w, h, planet.radius_m);
        let land = land_mask(&input.height_m);
        let elev: Vec<f32> = input.height_m.data.iter().zip(&land).map(|(&e, &l)| if l { e } else { 0.0 }).collect();
        let dist_coast = distance_to_ocean_km(&land, &geo);

        // Coast orientation along rows -> east/west side of basin.
        let (west_km, east_km) = row_runs_km(&land, &geo);
        let n = land.len();
        let mut east_side = vec![0.0f32; n];
        let mut west_side = vec![0.0f32; n];
        for i in 0..n {
            // Ocean cell: land to the east = eastern edge of the basin. Land cell: ocean to the
            // west = west coast of the continent = eastern edge of the ocean to its west.
            let (to_east_edge, to_west_edge, scale) = if land[i] {
                (west_km[i], east_km[i], CURRENT_SCALE_LAND_KM)
            } else {
                (east_km[i], west_km[i], CURRENT_SCALE_OCEAN_KM)
            };
            east_side[i] = (-to_east_edge / scale).exp();
            west_side[i] = (-to_west_edge / scale).exp();
        }

        let model_cells = Cells::new(planet);
        let mut cur_base = vec![0.0f32; n];
        let mut cur_winter = vec![0.0f32; n];
        for y in 0..h as usize {
            let (sub, pol) = Self::bands(&model_cells, geo.lat_deg[y].abs());
            for x in 0..w as usize {
                let i = y * w as usize + x;
                let (e, wst) = (east_side[i], west_side[i]);
                cur_base[i] = sub * (4.0 * wst - 5.0 * e) + pol * (3.0 + 6.0 * e - 3.0 * wst);
                cur_winter[i] = pol * (3.6 * e - 0.9 * wst);
            }
        }

        // Zonal EBMs, sampled at every grid row.
        let ocean_z = solve_zonal_ebm(planet, HEAT_CAPACITY_OCEAN);
        let land_z = solve_zonal_ebm(planet, HEAT_CAPACITY_LAND);
        let row_table = |z: &crate::ebm::ZonalClimate| -> Vec<[f32; 12]> {
            geo.lat.iter().map(|&lat| std::array::from_fn(|m| z.at(lat, m) as f32)).collect()
        };
        let ocean_row = row_table(&ocean_z);
        let land_row = row_table(&land_z);
        let mean = |rows: &[[f32; 12]]| -> Vec<f32> { rows.iter().map(|r| r.iter().sum::<f32>() / 12.0).collect() };
        let (ocean_mean, land_mean) = (mean(&ocean_row), mean(&land_row));

        Self {
            planet,
            cells: Cells::new(planet),
            geo,
            land,
            elev,
            dist_coast,
            east_side,
            west_side,
            cur_base,
            cur_winter,
            ocean_row,
            land_row,
            ocean_mean,
            land_mean,
        }
    }

    fn run(&self) -> ClimateFields {
        let (w, h) = (self.geo.w, self.geo.h);
        let n = w * h;
        let tilt = self.planet.axial_tilt_deg.to_radians();

        let mut temp = vec![[0.0f32; 12]; n];
        let mut precip = vec![[0.0f32; 12]; n];
        let mut wind_u = vec![[0.0f32; 12]; n];
        let mut wind_v = vec![[0.0f32; 12]; n];

        let mut s = Month {
            t_sea: vec![0.0; n],
            u: vec![0.0; n],
            v: vec![0.0; n],
            q: vec![0.0; n],
            oro: vec![0.0; n],
            fetch_km: vec![0.0; n],
            qcap: vec![0.0; n],
            itcz_deg: vec![0.0; w],
        };
        let mut t_final = vec![0.0f32; n];

        for m in 0..12 {
            // Season: the lagged declination drives temperature response and the rain belts.
            let decl_lag = declination(month_mid_day(m) - SEASON_LAG_DAYS, tilt);
            // +1 = deep winter in the northern hemisphere, -1 = deep winter in the south.
            let winter_north = if tilt.abs() < 1e-6 { 0.0 } else { (-decl_lag / tilt.abs()).clamp(-1.0, 1.0) as f32 };
            let shift_deg = BELT_SHIFT_PER_DECLINATION * decl_lag.to_degrees();

            // (a) provisional temperatures with plain-distance continentality.
            self.sea_level_temps(m, winter_north, None, &mut s.t_sea);
            // (b) rain belt position per longitude.
            self.thermal_equator(&s.t_sea, &mut s.itcz_deg);
            // (c) winds.
            self.winds(winter_north, shift_deg, &s.t_sea, &mut s.u, &mut s.v);
            // (d) upwind march.
            self.ocean_moisture_capacity(&s.t_sea, &mut s.qcap);
            self.upwind_march(&mut s);
            // (e) final temperature and precipitation.
            self.sea_level_temps(m, winter_north, Some(&s.fetch_km), &mut s.t_sea);
            let offset = self.planet.global_temp_offset_c as f32;
            for ((t, &sea), &e) in t_final.iter_mut().zip(&s.t_sea).zip(&self.elev) {
                *t = sea - LAPSE_C_PER_KM * e / 1000.0 + offset;
            }
            self.precipitation(shift_deg, &s, &t_final, &mut precip, m);
            for i in 0..n {
                temp[i][m] = t_final[i];
                wind_u[i][m] = s.u[i];
                wind_v[i][m] = s.v[i];
            }
        }

        // Köppen for land.
        let koppen = LatLonGrid::from_vec(
            w as u32,
            h as u32,
            (0..n)
                .map(|i| {
                    if self.land[i] {
                        let north = self.geo.lat[i / w] >= 0.0;
                        Koppen::classify(&temp[i], &precip[i], north) as u8
                    } else {
                        0
                    }
                })
                .collect(),
        );
        ClimateFields { width: w as u32, height: h as u32, temp_c: temp, precip_mm: precip, koppen, wind_u, wind_v }
    }

    // -----------------------------------------------------------------------------------------
    // Temperature
    // -----------------------------------------------------------------------------------------

    /// Latitude weights (0..1) of the two boundary-current bands, for the given cells.
    ///
    /// Boundary-current anomaly in C, from the coast orientation:
    /// - Subtropical gyre band (about 0.8 x the Hadley edge, 24 degrees on Earth): the eastern
    ///   side of an ocean is cooled by equatorward flow and upwelling (Canary, California,
    ///   Humboldt, Benguela, West Australian currents), the western side is warmed (Gulf Stream,
    ///   Kuroshio, Brazil, Agulhas, East Australian).
    /// - Subpolar band (about 1.7 x the Hadley edge): the pattern reverses, so west coasts of
    ///   continents get the warm North Atlantic Drift / Alaska current effect, strongest in
    ///   winter, and the western sides are cold (Labrador, Oyashio).
    ///
    /// Amplitudes are a few kelvin, as observed. The per-cell values are precomputed in
    /// [`Model::new`] as `cur_base + winter * cur_winter`.
    fn bands(cells: &Cells, abs_lat_deg: f64) -> (f32, f32) {
        let h = cells.hadley_deg;
        let g = |x: f64| (-x * x).exp() as f32;
        (g((abs_lat_deg - 0.8 * h) / (0.45 * h)), g((abs_lat_deg - 1.7 * h) / (0.40 * h)))
    }

    fn current_bands(&self, abs_lat_deg: f64) -> (f32, f32) {
        Self::bands(&self.cells, abs_lat_deg)
    }

    /// Sea-level monthly temperature for every cell.
    ///
    /// `fetch_km = None` uses the plain coast distance for continentality; with a fetch (km of
    /// land air travels over, see [`Model::upwind_march`]) the two distances are blended.
    ///
    /// For a land cell with continentality `c` (0 at the coast, 1 deep inland):
    /// `T = mean + (1-c) * ocean_swing + c * gain * land_swing`, where the means and swings come
    /// from the two zonal EBM runs. Ocean cells (`c = 0`) just follow the ocean EBM.
    fn sea_level_temps(&self, m: usize, winter_north: f32, fetch_km: Option<&[f32]>, out: &mut [f32]) {
        let (w, h) = (self.geo.w, self.geo.h);
        for y in 0..h {
            let lat_deg = self.geo.lat_deg[y];
            let winter = if lat_deg >= 0.0 { winter_north } else { -winter_north };
            // Land swing is damped less near the equator (where it is tiny anyway) and more at
            // the poles, where snow and ice melting caps the summer warmth.
            let fade = ((lat_deg.abs() - 55.0) / 30.0).clamp(0.0, 1.0) as f32;
            let gain = LAND_AMPLITUDE_GAIN * (1.0 - 0.35 * fade * fade * (3.0 - 2.0 * fade));
            // Land runs warmer than the zonal-mean (ocean-dominated) EBM: dry air, low albedo and
            // heat trapped over the continents. Strongest under the subtropical ridge (deserts).
            let ridge = (-((lat_deg.abs() - self.cells.hadley_deg * 0.85) / (0.35 * self.cells.hadley_deg)).powi(2)).exp() as f32;
            let tropics = (1.0 - (lat_deg.abs() / 50.0)).clamp(0.0, 1.0) as f32;
            let land_warm = LAND_WARMING_C * (0.5 * tropics + ridge);
            let (o_mean, l_mean) = (self.ocean_mean[y], self.land_mean[y]);
            let o_swing = self.ocean_row[y][m] - o_mean;
            let l_swing = self.land_row[y][m] - l_mean;
            for x in 0..w {
                let i = y * w + x;
                let current = self.cur_base[i] + winter * self.cur_winter[i];
                if !self.land[i] {
                    out[i] = self.ocean_row[y][m] + current;
                    continue;
                }
                let plain = self.dist_coast[i];
                let d_eff = match fetch_km {
                    Some(f) => (1.0 - UPWIND_WEIGHT) * plain + UPWIND_WEIGHT * f[i].max(plain * 0.5),
                    None => plain,
                };
                let c = 1.0 - (-d_eff / CONTINENTALITY_KM).exp();
                let mean = o_mean + LAND_MEAN_SHARE * c * (l_mean - o_mean);
                out[i] = mean + (1.0 - c) * o_swing + c * gain * l_swing + current + land_warm * (0.4 + 0.6 * c);
            }
        }
    }

    // -----------------------------------------------------------------------------------------
    // Rain belt position and winds
    // -----------------------------------------------------------------------------------------

    /// Latitude of the tropical rain belt at each longitude.
    ///
    /// The thermal equator of a longitude is the temperature-weighted centroid latitude of its
    /// warmest cells. Over open ocean it sits near the (lagged) sub-solar latitude; over a big
    /// continent it lies much further poleward in summer, and the rain belt follows it, though
    /// not all the way (a tanh limits it to about the tilt). That is the monsoon.
    fn thermal_equator(&self, t_sea: &[f32], itcz_deg: &mut [f32]) {
        let (w, h) = (self.geo.w, self.geo.h);
        let tilt_deg = self.planet.axial_tilt_deg.abs();
        let search = (tilt_deg + 14.0).clamp(14.0, 60.0);
        let limit = tilt_deg.clamp(8.0, 40.0);
        let mut raw = vec![0.0f32; w];
        for (x, r) in raw.iter_mut().enumerate() {
            let mut tmax = f32::MIN;
            for y in 0..h {
                if self.geo.lat_deg[y].abs() <= search {
                    tmax = tmax.max(t_sea[y * w + x]);
                }
            }
            let (mut sum_w, mut sum_wl) = (0.0f64, 0.0f64);
            for y in 0..h {
                if self.geo.lat_deg[y].abs() <= search {
                    let wgt = ((t_sea[y * w + x] - tmax) / 1.2).exp() as f64;
                    sum_w += wgt;
                    sum_wl += wgt * self.geo.lat_deg[y];
                }
            }
            *r = (sum_wl / sum_w.max(1e-9)) as f32;
        }
        // Smooth around the globe (about 600 km each side).
        let radius = ((600.0 / self.geo.dx_km[h / 2].max(1.0)).round() as isize).clamp(1, (w / 4) as isize);
        for (x, out) in itcz_deg.iter_mut().enumerate() {
            let mut acc = 0.0;
            for d in -radius..=radius {
                acc += raw[(x as isize + d).rem_euclid(w as isize) as usize];
            }
            let te = (acc / (2 * radius + 1) as f32) as f64;
            *out = (limit * (0.65 * te / limit).tanh()) as f32;
        }
    }

    /// Prevailing surface wind for every cell: the belt winds plus a thermal-low wind.
    ///
    /// Heated land lowers surface pressure, so air flows toward warm anomalies (summer monsoon
    /// onshore flow) and away from cold ones (winter monsoon). The anomaly is the sea-level
    /// temperature minus its latitude-circle mean, blurred over about 450 km; the wind is its
    /// gradient, turned 30 degrees by the Coriolis force (clockwise in the north).
    fn winds(&self, winter_north: f32, shift_deg: f64, t_sea: &[f32], u: &mut [f32], v: &mut [f32]) {
        let (w, h) = (self.geo.w, self.geo.h);
        let mut anomaly = vec![0.0f32; w * h];
        for y in 0..h {
            let row = &t_sea[y * w..(y + 1) * w];
            let mean = row.iter().sum::<f32>() / w as f32;
            for x in 0..w {
                anomaly[y * w + x] = row[x] - mean;
            }
        }
        let radius = ((450.0 / self.geo.ds_km).round() as usize).clamp(1, 8);
        let smooth = box_blur(&box_blur(&anomaly, w, h, radius), w, h, radius);

        for y in 0..h {
            let lat_deg = self.geo.lat_deg[y];
            let (u0, v0) = band_wind(lat_deg, shift_deg, &self.cells);
            let (sub, _) = self.current_bands(lat_deg.abs());
            let gyre_season = 1.0 - 0.5 * winter_north_for(lat_deg, winter_north);
            let lat_sign = lat_deg.signum() as f32;
            let dx = self.geo.dx_km[y].max(self.geo.ds_km * 0.1) as f32;
            let ds = self.geo.ds_km as f32;
            let fade = (1.0 - ((lat_deg.abs() - 60.0) / 25.0).clamp(0.0, 1.0)) as f32;
            let angle = (-30.0f64.to_radians() * lat_deg.signum() * lat_deg.to_radians().sin().abs().sqrt()) as f32;
            let (sa, ca) = angle.sin_cos();
            for x in 0..w {
                let i = y * w + x;
                let (xe, xw) = ((x + 1) % w, (x + w - 1) % w);
                let (yn, ys) = (y.saturating_sub(1), (y + 1).min(h - 1));
                let gx = (smooth[y * w + xe] - smooth[y * w + xw]) / (2.0 * dx);
                let gy = (smooth[yn * w + x] - smooth[ys * w + x]) / ((yn.abs_diff(ys)).max(1) as f32 * ds);
                let (mx, my) = (gx * 1000.0 * MONSOON_WIND_GAIN, gy * 1000.0 * MONSOON_WIND_GAIN);
                let (mx, my) = (mx * ca - my * sa, mx * sa + my * ca);
                let mag = (mx * mx + my * my).sqrt();
                let k = if mag > MONSOON_WIND_MAX { MONSOON_WIND_MAX / mag } else { 1.0 } * fade;
                // Oceanic subtropical highs: poleward flow on the western flank (the east coasts
                // of continents, onshore, bringing warm moist air) and equatorward flow on the
                // eastern flank. Weaker in winter, when the continents are cold.
                let flank = self.west_side[i] - self.east_side[i];
                let gyre_v = GYRE_WIND_MS * sub * flank * gyre_season * lat_sign;
                u[i] = u0 as f32 + mx * k;
                v[i] = v0 as f32 + my * k + gyre_v;
            }
        }
    }

    // -----------------------------------------------------------------------------------------
    // Moisture
    // -----------------------------------------------------------------------------------------

    /// How much water vapour air can hold over each ocean cell, relative to 22 C sea water
    /// (Clausius-Clapeyron via the Tetens formula, about +7 % per kelvin). Land cells get 0.
    fn ocean_moisture_capacity(&self, t_sea: &[f32], qcap: &mut [f32]) {
        let sat = |t: f32| (17.27 * t / (t + 237.3)).exp();
        let reference = sat(22.0);
        for i in 0..qcap.len() {
            qcap[i] = if self.land[i] { 0.0 } else { (sat(t_sea[i].max(-30.0)) / reference).min(1.0) };
        }
    }

    /// Walk upwind from every cell and track the moisture in the air that arrives.
    ///
    /// For each cell we follow the wind backwards along a straight line of up to about 2600 km
    /// over land (700 km over ocean, where the answer saturates quickly) and then run forward
    /// again to the cell:
    /// - over **ocean** the air picks up moisture, up to the capacity set by sea temperature
    ///   (cold seas give dry air);
    /// - over **land** it slowly loses moisture to rain (`exp(-step / L)`, with `L` shorter in
    ///   the dry subtropics and longer in the wet tropics);
    /// - when the ground **rises** between two steps a fraction `rise / (rise + 700 m)` of the
    ///   moisture is wrung out. At the last step that is the orographic rain the cell receives;
    ///   at earlier steps it just dries the air, so the lee of a range is a rain shadow;
    /// - the distance to the first upwind ocean cell is recorded as the **fetch** (a measure of
    ///   continentality that respects wind direction).
    fn upwind_march(&self, s: &mut Month) {
        let (w, h) = (self.geo.w as i32, self.geo.h as i32);
        let ds = self.geo.ds_km as f32;
        let steps_land = ((2600.0 / ds).ceil() as usize).clamp(6, 56);
        let steps_ocean = ((700.0 / ds).ceil() as usize).clamp(3, 14);
        // Positions are tracked in 1/256 cell units so the inner loop needs no float rounding.
        const FP: i32 = 256;

        for y in 0..h as usize {
            let a = self.geo.lat_deg[y].abs();
            let hd = self.cells.hadley_deg;
            // Moisture e-folding distance over land by latitude (km).
            let decay_km = piecewise_cos(
                a,
                &[(0.0, 3600.0), (0.6 * hd, 3000.0), (hd, 2200.0), (hd + 12.0, 3200.0), (90.0, 2000.0)],
            ) as f32;
            let keep = (-ds / decay_km).exp();
            // Under the subtropical ridge the sinking air is dry and land cannot moisten it much.
            let ridge = (-((a - 0.9 * hd) / (0.3 * hd)).powi(2)).exp() as f32;
            let floor = LAND_EVAPORATION_FLOOR * (1.0 - 0.85 * ridge);
            let dx = self.geo.dx_km[y].max(self.geo.ds_km * 0.25) as f32;
            for x in 0..w as usize {
                let i = y * w as usize + x;
                let is_land = self.land[i];
                let (u, v) = (s.u[i], s.v[i]);
                let speed = (u * u + v * v).sqrt().max(1e-3);
                // Upwind step in cells: against the wind, in longitude scaled by the row width.
                let dcol = (-(u / speed) * ds / dx * FP as f32).round() as i32;
                let drow = ((v / speed) * FP as f32).round() as i32;
                let steps = if is_land { steps_land } else { steps_ocean };

                // The march is evaluated nearest-step-first as a chain of affine maps
                // `q_before_step = a * q_upwind + b`: step 1 (into this cell) is kept apart because
                // its rain is the cell's orographic rain; steps 2.. are folded into (acc_b, acc_a).
                let (mut a1, mut b1, mut f1) = (1.0f32, 0.0f32, 0.0f32);
                let (mut acc_b, mut acc_a) = (0.0f32, 1.0f32);
                let mut prev = i;
                let mut far = i;
                let mut fetch = None;
                let (mut px, mut py) = (x as i32 * FP, y as i32 * FP);
                for k in 1..=steps {
                    px += dcol;
                    py += drow;
                    let yy = (py + FP / 2) >> 8;
                    if yy < 0 || yy >= h {
                        break;
                    }
                    let mut xx = (px + FP / 2) >> 8;
                    while xx < 0 {
                        xx += w;
                    }
                    while xx >= w {
                        xx -= w;
                    }
                    let j = (yy * w + xx) as usize;
                    // Map for the step from cell `j` (upwind) into cell `prev`.
                    let (am, bm, f) = if !self.land[prev] {
                        (1.0 - OCEAN_RECHARGE, OCEAN_RECHARGE * s.qcap[prev], 0.0)
                    } else {
                        let rise = self.elev[prev] - self.elev[j];
                        let f = if rise > 0.0 { rise / (rise + OROGRAPHIC_HALF_RISE_M) } else { 0.0 };
                        ((1.0 - f) * keep, (1.0 - keep) * floor, f)
                    };
                    if k == 1 {
                        (a1, b1, f1) = (am, bm, f);
                    } else {
                        acc_b += acc_a * bm;
                        acc_a *= am;
                    }
                    if fetch.is_none() && !self.land[j] {
                        fetch = Some(k as f32 * ds);
                    }
                    prev = j;
                    far = j;
                    if acc_a < 0.02 {
                        break; // air from further upwind no longer matters
                    }
                }
                let kmax_dist = ds * steps as f32;
                s.fetch_km[i] = match (is_land, fetch) {
                    (false, _) => 0.0,
                    (true, Some(f)) => f,
                    (true, None) => kmax_dist.max(self.dist_coast[i]),
                };
                let q_start = if self.land[far] { MOISTURE_START_OVER_LAND } else { s.qcap[far] };
                let q_before = acc_b + acc_a * q_start;
                let speed_factor = (speed / 5.0).clamp(0.5, 1.5);
                s.oro[i] = if is_land { q_before * f1 * OROGRAPHIC_RAIN_MM * speed_factor } else { 0.0 };
                s.q[i] = if far == i { q_before } else { a1 * q_before + b1 }.clamp(0.0, 1.2);
            }
        }
    }

    // -----------------------------------------------------------------------------------------
    // Precipitation
    // -----------------------------------------------------------------------------------------

    /// Combine the circulation rain belts with moisture, orography and convection.
    ///
    /// `P = [belts(lat, itcz, season) x moisture^0.5 + orographic x dryness + convective]
    ///      x temperature factor x (1 - upwelling drying) x global multiplier`
    fn precipitation(&self, shift_deg: f64, s: &Month, t_final: &[f32], out: &mut [[f32; 12]], m: usize) {
        let (w, h) = (self.geo.w, self.geo.h);
        let sat = |t: f32| (17.27 * t / (t + 237.3)).exp();
        let sat_ref = sat(22.0);
        for y in 0..h {
            let lat_deg = self.geo.lat_deg[y];
            let a = lat_deg.abs();
            let (sub, _) = self.current_bands(a);
            // Convective rain fades out poleward of the Ferrel cell, where it is too cold anyway.
            let conv_lat = (1.0 - ((a - 62.0) / 18.0).clamp(0.0, 1.0)) as f32;
            // Convection needs less warmth to start at high latitudes (long summer days).
            let conv_threshold = (15.0 - 7.0 * ((a - 35.0) / 30.0).clamp(0.0, 1.0)) as f32;
            // Strength of the subtropical ridge here (sinking air), in the season's shifted frame.
            let row_belts = RowBelts::new(lat_deg, shift_deg, &self.cells);
            let ridge = row_belts.ridge as f32;
            // Moisture that storms draw from the sea around a coast, whatever the wind: the
            // capacity of the zonal-mean ocean at this latitude.
            let q_sea = (sat(self.ocean_row[y][m].max(-30.0)) / sat_ref).min(1.0);
            // The monsoon term lives in the tropics and warm mid-latitudes.
            let mon_lat = (1.0 - ((a - 35.0) / 20.0).clamp(0.0, 1.0)) as f32;
            for x in 0..w {
                let i = y * w + x;
                let belts = row_belts.at(s.itcz_deg[x] as f64) as f32;
                let mut q = s.q[i];
                if self.land[i] {
                    // Coastal air is also fed by the sea around it (cyclones, sea breezes).
                    q += (q_sea - q).max(0.0) * COASTAL_MOISTURE_SHARE * (-self.dist_coast[i] / COASTAL_MOISTURE_KM).exp();
                }
                let mut p = belts * q.sqrt();
                if self.land[i] {
                    let dryness = (belts / 70.0).clamp(0.25, 1.0);
                    p += s.oro[i] * dryness;
                    // Air arriving from the tropics (poleward surface flow) is warm and moist and
                    // is not subject to the subtropical ridge's subsidence: east coasts of
                    // continents, monsoon regions. Equatorward or zonal flow under the ridge is
                    // dry: west coasts and deserts (Mediterranean summer, California, Sahara).
                    let inflow = (s.v[i] * lat_deg.signum() as f32 / 2.0).clamp(0.0, 1.0);
                    let sinking = (1.0 - 0.9 * ridge * (1.0 - inflow)).powi(2);
                    let warm = (t_final[i] - conv_threshold).max(0.0);
                    p += CONVECTIVE_MM_PER_K * warm * q * q.sqrt() * conv_lat * sinking;
                    // Monsoon: warm land with moist tropical air coming from the sea nearby.
                    let onshore = (-s.fetch_km[i] / MONSOON_FETCH_KM).exp();
                    let warmth = ((t_final[i] - 12.0) / 10.0).clamp(0.0, 1.0);
                    p += MONSOON_RAIN_MM * q * q * onshore * warmth * inflow * mon_lat;
                }
                // Cold air holds little water: damp the rain toward the poles in winter.
                let cold = 0.6 + 0.4 * ((t_final[i] + 15.0) / 30.0).clamp(0.0, 1.0);
                // Upwelling (cold-current) coasts are foggy deserts: little evaporation and rain.
                let upwell = 1.0 - 0.55 * sub * self.east_side[i];
                out[i][m] = (p * cold * upwell * self.planet.global_precip_mult as f32).max(0.0);
            }
        }
    }
}

/// Winter factor (+1 deep winter, -1 deep summer) at this latitude's hemisphere.
fn winter_north_for(lat_deg: f64, winter_north: f32) -> f32 {
    if lat_deg >= 0.0 { winter_north } else { -winter_north }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::synthetic::earth_like_height;

    fn earth_run() -> (LatLonGrid<f32>, ClimateFields) {
        let height = earth_like_height(360, 180);
        let fields = tier_a(&PlanetParams::default(), &ClimateInput { height_m: height.clone() });
        (height, fields)
    }

    fn cell(f: &ClimateFields, lon: f64, lat: f64) -> usize {
        let x = (((lon + 180.0) / 360.0 * f.width as f64) as u32).min(f.width - 1);
        let y = (((90.0 - lat) / 180.0 * f.height as f64) as u32).min(f.height - 1);
        f.index(x, y)
    }

    fn class_at(f: &ClimateFields, lon: f64, lat: f64) -> String {
        Koppen::from_u8(f.koppen.data[cell(f, lon, lat)]).map_or("ocean".to_string(), |k| k.code().to_string())
    }

    fn annual_t(f: &ClimateFields, i: usize) -> f32 {
        f.temp_c[i].iter().sum::<f32>() / 12.0
    }

    fn annual_p(f: &ClimateFields, i: usize) -> f32 {
        f.precip_mm[i].iter().sum()
    }

    #[test]
    fn earth_like_planet_is_plausible() {
        let start = std::time::Instant::now();
        let (height, f) = earth_run();
        let elapsed = start.elapsed();
        println!("tier_a 360x180 (incl. building the test planet): {elapsed:?}");
        // The target is < 300 ms in release on one core (the `earthlike` example prints the real
        // number, about 160-250 ms). Tests run in parallel on a shared machine, so this bound
        // only guards against order-of-magnitude regressions.
        let budget_ms = if cfg!(debug_assertions) { 20_000 } else { 1_500 };
        assert!(elapsed.as_millis() < budget_ms, "too slow: {elapsed:?}");
        assert_eq!((f.width, f.height), (360, 180));

        // Equatorial rainforest: hot, very wet, class A.
        let congo = cell(&f, 22.0, 0.0);
        assert!(annual_t(&f, congo) > 24.0 && annual_p(&f, congo) > 900.0);
        assert!(class_at(&f, 22.0, 0.0).starts_with('A'));
        let amazon = class_at(&f, -60.0, -3.0);
        assert!(amazon == "Af" || amazon == "Am", "Amazon is {amazon}");

        // Hot deserts on the subtropical ridge, including a west coast.
        for &(lon, lat) in &[(10.0, 24.0), (-8.0, 22.0), (133.0, -25.0), (17.0, -25.0)] {
            let c = class_at(&f, lon, lat);
            assert!(c.starts_with('B'), "({lon},{lat}) should be arid, got {c}");
            assert!(annual_p(&f, cell(&f, lon, lat)) < 350.0);
        }

        // Mediterranean climates on west coasts at 35-40 degrees: dry summer, wet winter.
        for &(lon, lat) in &[(-5.0, 38.0), (-121.0, 36.0)] {
            let c = class_at(&f, lon, lat);
            assert!(c.starts_with("Cs"), "({lon},{lat}) should be Mediterranean, got {c}");
            let i = cell(&f, lon, lat);
            assert!(f.precip_mm[i][6] < 0.35 * f.precip_mm[i][0], "dry summer expected");
        }

        // Cold high latitudes.
        let siberia = class_at(&f, 100.0, 65.0);
        assert!(siberia.starts_with('D'), "Siberia is {siberia}");
        assert!(annual_t(&f, cell(&f, 100.0, 65.0)) < -2.0);
        assert_eq!(class_at(&f, 0.0, -80.0), "EF");
        assert!(["EF", "ET"].contains(&class_at(&f, -40.0, 72.0).as_str()));

        // Seasons: northern land is warmest in July, southern in January, and mid-latitude
        // continental interiors swing far more than the nearby ocean.
        let chicago = cell(&f, -90.0, 42.0);
        assert!(f.temp_c[chicago][6] > f.temp_c[chicago][0] + 15.0);
        let outback = cell(&f, 133.0, -30.0);
        assert!(f.temp_c[outback][0] > f.temp_c[outback][6] + 8.0);
        let atlantic = cell(&f, -40.0, 42.0);
        let swing = |i: usize| f.temp_c[i][6] - f.temp_c[i][0];
        assert!(swing(chicago) > 2.0 * swing(atlantic), "continentality: {} vs {}", swing(chicago), swing(atlantic));

        // Rain shadow: the windward (west) side of the Rockies is much wetter than the lee.
        assert!(annual_p(&f, cell(&f, -124.0, 45.0)) > 2.5 * annual_p(&f, cell(&f, -105.0, 46.0)));

        // Lapse rate: Tibet is far colder than lowland at the same latitude.
        assert!(annual_t(&f, cell(&f, 88.0, 33.0)) < annual_t(&f, cell(&f, 118.0, 33.0)) - 10.0);

        // Zonal structure over the ocean: warm equator, cooler poleward, both hemispheres.
        let ocean_mean = |lat: f64| -> f32 {
            let y = ((90.0 - lat) / 180.0 * f.height as f64) as u32;
            let mut sum = 0.0;
            let mut n = 0.0;
            for x in 0..f.width {
                if *height.get(x, y) < 0.0 {
                    sum += annual_t(&f, f.index(x, y));
                    n += 1.0;
                }
            }
            sum / n
        };
        let (eq, n30, n60, s30, s60) = (ocean_mean(0.0), ocean_mean(30.0), ocean_mean(60.0), ocean_mean(-30.0), ocean_mean(-60.0));
        assert!((24.0..32.0).contains(&eq), "equatorial ocean {eq}");
        assert!(eq > n30 && n30 > n60 && eq > s30 && s30 > s60, "{eq} {n30} {n60} {s30} {s60}");

        // Köppen raster: zero over ocean, a real class over land.
        for (i, &h) in height.data.iter().enumerate() {
            assert_eq!(f.koppen.data[i] == 0, h < 0.0);
        }
    }

    #[test]
    fn winds_follow_the_three_cell_pattern() {
        let (height, f) = earth_run();
        let mean_u = |lat: f64| -> f32 {
            let y = ((90.0 - lat) / 180.0 * f.height as f64) as u32;
            let (mut sum, mut n) = (0.0, 0.0);
            for x in 0..f.width {
                if *height.get(x, y) < 0.0 {
                    let i = f.index(x, y);
                    sum += f.wind_u[i].iter().sum::<f32>() / 12.0;
                    n += 1.0;
                }
            }
            sum / n
        };
        assert!(mean_u(15.0) < -2.0, "trades");
        assert!(mean_u(-15.0) < -2.0, "trades");
        assert!(mean_u(45.0) > 2.0 && mean_u(-45.0) > 2.0, "westerlies");
        assert!(mean_u(75.0) < 0.0, "polar easterlies");
    }

    #[test]
    fn output_is_finite_in_range_and_deterministic() {
        let (_, a) = earth_run();
        let (_, b) = earth_run();
        assert_eq!(a.temp_c, b.temp_c);
        assert_eq!(a.precip_mm, b.precip_mm);
        assert_eq!(a.koppen, b.koppen);
        for i in 0..a.temp_c.len() {
            for m in 0..12 {
                assert!((-90.0..70.0).contains(&a.temp_c[i][m]), "T {}", a.temp_c[i][m]);
                assert!((0.0..1500.0).contains(&a.precip_mm[i][m]), "P {}", a.precip_mm[i][m]);
                assert!(a.wind_u[i][m].is_finite() && a.wind_v[i][m].is_finite());
            }
        }
    }

    #[test]
    fn ocean_world_is_hemispherically_symmetric() {
        let height = LatLonGrid::filled(120, 60, -4000.0f32);
        let f = tier_a(&PlanetParams::default(), &ClimateInput { height_m: height });
        assert!(f.koppen.data.iter().all(|&k| k == 0));
        for &lat in &[10.0, 30.0, 50.0, 70.0] {
            // Mirror rows exactly (cell() would round the two latitudes differently).
            let north = cell(&f, 0.0, lat);
            let south = f.index(f.width / 2, f.height - 1 - (north / f.width as usize) as u32);
            for m in 0..12 {
                let (tn, ts) = (f.temp_c[north][m], f.temp_c[south][(m + 6) % 12]);
                assert!((tn - ts).abs() < 2.0, "lat {lat} month {m}: {tn} vs {ts}");
            }
        }
    }

    #[test]
    fn planet_parameters_do_what_they_say() {
        let height = earth_like_height(180, 90);
        let run = |p: PlanetParams| tier_a(&p, &ClimateInput { height_m: height.clone() });
        let base = run(PlanetParams::default());
        let chicago = |f: &ClimateFields| {
            let i = cell(f, -90.0, 42.0);
            (f.temp_c[i][6] - f.temp_c[i][0], annual_t(f, i), annual_p(f, i))
        };
        let (swing0, t0, p0) = chicago(&base);

        // No tilt, no seasons.
        let flat = run(PlanetParams { axial_tilt_deg: 0.0, ..Default::default() });
        assert!(chicago(&flat).0.abs() < 3.0, "tilt 0 swing {}", chicago(&flat).0);
        // More tilt, bigger seasons.
        let tilted = run(PlanetParams { axial_tilt_deg: 40.0, ..Default::default() });
        assert!(chicago(&tilted).0 > swing0);
        // Offset and multiplier are plain shifts.
        let warm = run(PlanetParams { global_temp_offset_c: 5.0, ..Default::default() });
        assert!((chicago(&warm).1 - t0 - 5.0).abs() < 1e-3);
        let wet = run(PlanetParams { global_precip_mult: 2.0, ..Default::default() });
        assert!((chicago(&wet).2 / p0 - 2.0).abs() < 1e-3);
        // A dimmer sun is colder.
        let dim = run(PlanetParams { solar_constant_rel: 0.9, ..Default::default() });
        assert!(chicago(&dim).1 < t0 - 3.0);
        // Slower rotation widens the Hadley cell, so the subtropical desert belt moves poleward.
        let slow = run(PlanetParams { rotation_hours: 40.0, ..Default::default() });
        let dry_lat = |f: &ClimateFields| -> f64 {
            // Latitude of the driest cell in a north Pacific column (ocean from 5 to 54 N).
            let mut best = (f32::MAX, 0.0);
            for y in 0..f.height / 2 {
                let lat = 90.0 - (y as f64 + 0.5) * 2.0;
                if !(5.0..54.0).contains(&lat) {
                    continue;
                }
                let x = 45; // longitude -135
                let p = annual_p(f, f.index(x, y));
                if p < best.0 {
                    best = (p, lat);
                }
            }
            best.1
        };
        assert!(dry_lat(&slow) > dry_lat(&base) + 5.0, "dry belt {} vs {}", dry_lat(&slow), dry_lat(&base));
    }

    #[test]
    fn works_at_other_resolutions_and_agrees_roughly() {
        let mut means = Vec::new();
        for &(w, h) in &[(72u32, 36u32), (180, 90), (360, 180)] {
            let height = earth_like_height(w, h);
            let f = tier_a(&PlanetParams::default(), &ClimateInput { height_m: height });
            assert_eq!(f.temp_c.len(), (w * h) as usize);
            means.push(annual_t(&f, cell(&f, 22.0, 0.0)));
            assert!(f.temp_c.iter().all(|t| t.iter().all(|v| v.is_finite())));
        }
        assert!((means[0] - means[2]).abs() < 3.0 && (means[1] - means[2]).abs() < 2.0, "{means:?}");
    }

    #[test]
    fn all_land_planet_does_not_break() {
        let height = LatLonGrid::filled(90, 45, 100.0f32);
        let f = tier_a(&PlanetParams::default(), &ClimateInput { height_m: height });
        assert!(f.koppen.data.iter().all(|&k| k != 0));
        assert!(f.temp_c.iter().all(|t| t.iter().all(|v| v.is_finite())));
    }
}
