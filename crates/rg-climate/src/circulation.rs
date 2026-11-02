//! Idealised atmospheric circulation: where the Hadley / Ferrel / polar cells are, the
//! prevailing surface winds they imply, and the latitude profile of rainfall they produce.
//!
//! These are lookup-style rules rather than a dynamical model (Tier B will compute pressure and
//! winds properly). Sources for the ideas: the three-cell model of the general circulation
//! (e.g. Hartmann, *Global Physical Climatology*, ch. 6-7) and Geoff Eddy's "Climate Cookbook"
//! for the worldbuilding-style rain belts.
//!
//! Latitudes in this module are degrees, north positive.

use crate::types::PlanetParams;

/// Where the circulation cells end, in degrees of latitude (measured from the shifted equator).
#[derive(Clone, Copy, Debug)]
pub struct Cells {
    /// Edge of the Hadley cell = the subtropical high-pressure ridge. Earth about 30 degrees.
    pub hadley_deg: f64,
    /// Edge of the Ferrel cell = the polar front / storm-track boundary. Earth about 60 degrees.
    pub ferrel_deg: f64,
}

impl Cells {
    /// Cell widths scale with rotation: the Hadley edge goes roughly as `1/(Omega * radius)^0.7`
    /// (Held & Hou 1980 give an exponent of 1 for the idealised case; GCM experiments land
    /// nearer 0.6-0.8). A fast spinner gets narrow cells, a slow one a wide Hadley cell.
    pub fn new(planet: &PlanetParams) -> Self {
        let omega_rel = 24.0 / planet.rotation_hours.max(0.5);
        let radius_rel = planet.radius_m / 6_371_000.0;
        let spin = (omega_rel * radius_rel).clamp(0.2, 6.0);
        let hadley = (30.0 * spin.powf(-0.7)).clamp(10.0, 55.0);
        let ferrel = (2.0 * hadley).clamp(hadley + 10.0, 80.0);
        Self { hadley_deg: hadley, ferrel_deg: ferrel }
    }
}

/// Smooth (cosine) interpolation through `(x, y)` knots sorted by `x`; flat outside the range.
pub fn piecewise_cos(x: f64, knots: &[(f64, f64)]) -> f64 {
    if x <= knots[0].0 {
        return knots[0].1;
    }
    for pair in knots.windows(2) {
        let ((x0, y0), (x1, y1)) = (pair[0], pair[1]);
        if x <= x1 {
            let t = (x - x0) / (x1 - x0).max(1e-9);
            let s = 0.5 - 0.5 * (std::f64::consts::PI * t).cos();
            return y0 + (y1 - y0) * s;
        }
    }
    knots[knots.len() - 1].1
}

/// Prevailing surface wind `(u east, v north)` in m/s at `lat_deg`, with the whole belt system
/// shifted north by `shift_deg` (the seasonal migration of the thermal equator).
///
/// Easterly trades in the tropics, westerlies in mid-latitudes, polar easterlies, with the
/// matching surface flow toward the equator in the trades, poleward in the Ferrel zone and
/// equatorward again at the poles. A weak easterly floor at the equator keeps the wind
/// direction defined in the doldrums.
pub fn band_wind(lat_deg: f64, shift_deg: f64, cells: &Cells) -> (f64, f64) {
    let (h, f) = (cells.hadley_deg, cells.ferrel_deg);
    let rel = lat_deg - shift_deg;
    let a = rel.abs().min(90.0);
    let mid_ferrel = 0.5 * (h + f);
    let mid_polar = 0.5 * (f + 90.0);
    let u = piecewise_cos(
        a,
        &[(0.0, -2.0), (0.5 * h, -5.5), (h, 0.0), (mid_ferrel, 7.0), (f, 0.0), (mid_polar, -3.5), (90.0, -1.0)],
    );
    let v_poleward = piecewise_cos(
        a,
        &[(0.0, 0.0), (0.5 * h, -2.4), (h, 0.0), (mid_ferrel, 2.2), (f, 0.0), (mid_polar, -1.5), (90.0, 0.0)],
    );
    (u, v_poleward * rel.signum())
}

/// Peak monthly rain of the tropical convergence zone over open ocean, mm.
const RAIN_ITCZ: f64 = 300.0;
/// Peak monthly rain of the mid-latitude storm track over open ocean, mm.
const RAIN_STORM: f64 = 145.0;

#[inline]
fn gauss(x: f64) -> f64 {
    (-x * x).exp()
}

/// Rainfall (mm/month) a location would get from the large-scale circulation alone, if it had
/// unlimited moisture: a wet convergence zone at `itcz_deg`, a wet storm track on the poleward
/// side of the subtropical ridge, dry subtropics in between, and a thin polar background.
///
/// `itcz_deg` is where the tropical rain belt sits at this longitude and month (it wanders
/// further from the equator over hot land: monsoons). `shift_deg` is the migration of the
/// subtropical ridge and storm tracks (smaller than the ITCZ's wander, and the same at all
/// longitudes). Because the storm track moves poleward in summer, west coasts near 35-45
/// degrees get a dry summer and a wet winter: the Mediterranean regime emerges.
pub fn base_precip_mm(lat_deg: f64, itcz_deg: f64, shift_deg: f64, cells: &Cells) -> f64 {
    RowBelts::new(lat_deg, shift_deg, cells).at(itcz_deg)
}

/// The part of [`base_precip_mm`] that depends only on latitude and season (not on where the
/// rain belt sits at a given longitude), so a whole grid row can share it.
#[derive(Clone, Copy, Debug)]
pub struct RowBelts {
    lat_deg: f64,
    /// Half-width of the rain belt in degrees.
    itcz_width: f64,
    /// Storm-track plus background rain, mm/month.
    fixed_mm: f64,
    /// Strength (0..1) of the subtropical ridge here.
    pub ridge: f64,
}

impl RowBelts {
    pub fn new(lat_deg: f64, shift_deg: f64, cells: &Cells) -> Self {
        let (h, f) = (cells.hadley_deg, cells.ferrel_deg);
        let a = (lat_deg - shift_deg).abs();
        let storm_centre = h + 0.55 * (f - h);
        let storm = RAIN_STORM * gauss((a - storm_centre) / (0.40 * (f - h)).max(6.0));
        // Sinking air under the subtropical ridge suppresses the storms and the background drizzle.
        let ridge = gauss((a - h) / (0.35 * h).max(4.0));
        let cos_lat = lat_deg.to_radians().cos();
        let background = (8.0 + 12.0 * cos_lat * cos_lat) * (1.0 - 0.5 * ridge);
        Self {
            lat_deg,
            itcz_width: (0.30 * h).clamp(5.0, 16.0),
            fixed_mm: storm * (1.0 - RIDGE_DRYING * ridge) + background,
            ridge,
        }
    }

    /// Rain (mm/month) with the tropical rain belt centred at `itcz_deg`. The belt is
    /// lopsided: it fades faster on its poleward flank, where the sinking air of the subtropical
    /// ridge begins.
    pub fn at(&self, itcz_deg: f64) -> f64 {
        let d = self.lat_deg - itcz_deg;
        let width = self.itcz_width * if d * self.lat_deg.signum() > 0.0 { 0.6 } else { 1.0 };
        RAIN_ITCZ * gauss(d / width) + self.fixed_mm
    }
}

/// Fraction of storm-track rain removed right under the subtropical ridge.
const RIDGE_DRYING: f64 = 0.65;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn earth_cells_are_30_and_60() {
        let c = Cells::new(&PlanetParams::default());
        assert!((c.hadley_deg - 30.0).abs() < 1e-9);
        assert!((c.ferrel_deg - 60.0).abs() < 1e-9);
    }

    #[test]
    fn faster_rotation_narrows_cells() {
        let fast = Cells::new(&PlanetParams { rotation_hours: 12.0, ..Default::default() });
        let slow = Cells::new(&PlanetParams { rotation_hours: 72.0, ..Default::default() });
        assert!(fast.hadley_deg < 25.0 && slow.hadley_deg > 40.0);
    }

    #[test]
    fn winds_have_the_three_cell_pattern() {
        let c = Cells::new(&PlanetParams::default());
        let (u15, v15) = band_wind(15.0, 0.0, &c);
        assert!(u15 < -3.0 && v15 < 0.0, "NE trades: u {u15} v {v15}");
        let (u_s15, v_s15) = band_wind(-15.0, 0.0, &c);
        assert!(u_s15 < -3.0 && v_s15 > 0.0, "SE trades");
        let (u45, v45) = band_wind(45.0, 0.0, &c);
        assert!(u45 > 4.0 && v45 > 0.0, "SW westerlies");
        let (u75, v75) = band_wind(75.0, 0.0, &c);
        assert!(u75 < 0.0 && v75 < 0.0, "polar easterlies");
    }

    #[test]
    fn rain_belts_are_in_the_right_places() {
        let c = Cells::new(&PlanetParams::default());
        let r = |lat: f64| base_precip_mm(lat, 0.0, 0.0, &c);
        assert!(r(0.0) > r(28.0) * 5.0, "equator much wetter than the subtropics");
        assert!(r(47.0) > r(28.0) * 2.0, "storm track wetter than the subtropical ridge");
        assert!(r(85.0) < r(47.0) / 3.0);
        // Mediterranean: at 40N winter (storms shifted south) beats summer (shifted north).
        let winter = base_precip_mm(40.0, -3.0, -3.0, &c);
        let summer = base_precip_mm(40.0, 9.0, 9.0, &c);
        assert!(winter > 2.0 * summer);
    }
}
