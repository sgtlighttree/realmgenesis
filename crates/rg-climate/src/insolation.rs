//! Top-of-atmosphere insolation from orbital geometry.
//!
//! The planet's orbit is taken as circular (eccentricity 0): good enough for a believable
//! seasonal cycle and it keeps the year symmetric. Formulas are the standard ones, e.g.
//! Berger (1978) / Hartmann, "Global Physical Climatology" (1994) section 2.3, and are the same
//! maths used by climlab's `DailyInsolation`.

use std::f64::consts::PI;

/// Earth's mean solar constant, W/m2.
pub const SOLAR_CONSTANT: f64 = 1361.0;

/// Day of the year (0..365) at which the northern spring equinox falls.
const EQUINOX_DAY: f64 = 79.5;
const DAYS_PER_YEAR: f64 = 365.0;

/// Mid-point of calendar month `m` (0 = January), as a day of the year. All months are the
/// same length (365/12 days); the crate treats a "month" as one twelfth of the orbit.
pub fn month_mid_day(m: usize) -> f64 {
    (m as f64 + 0.5) * DAYS_PER_YEAR / 12.0
}

/// Solar declination (radians) on `day` of the year for an axial tilt of `obliquity_rad`.
/// With a circular orbit, `sin(delta) = sin(obliquity) * sin(orbital longitude)`.
pub fn declination(day: f64, obliquity_rad: f64) -> f64 {
    let orbital_longitude = 2.0 * PI * (day - EQUINOX_DAY) / DAYS_PER_YEAR;
    (obliquity_rad.sin() * orbital_longitude.sin()).clamp(-1.0, 1.0).asin()
}

/// Daily-mean top-of-atmosphere insolation in W/m2 at latitude `lat_rad` when the sun is at
/// declination `decl_rad`. `s0` is the solar constant (W/m2).
///
/// `Q = (S0/pi) * (h0 sin(lat) sin(decl) + cos(lat) cos(decl) sin(h0))`, where `h0` is the
/// sunset hour angle `cos(h0) = -tan(lat) tan(decl)`. Polar day gives `h0 = pi`, polar night
/// `h0 = 0`.
pub fn daily_mean_insolation(lat_rad: f64, decl_rad: f64, s0: f64) -> f64 {
    let lat = lat_rad.clamp(-PI / 2.0 + 1e-6, PI / 2.0 - 1e-6);
    let cos_h0 = -lat.tan() * decl_rad.tan();
    let h0 = if cos_h0 >= 1.0 {
        0.0
    } else if cos_h0 <= -1.0 {
        PI
    } else {
        cos_h0.acos()
    };
    let q = s0 / PI * (h0 * lat.sin() * decl_rad.sin() + lat.cos() * decl_rad.cos() * h0.sin());
    q.max(0.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn equinox_equator_is_s0_over_pi() {
        let q = daily_mean_insolation(0.0, 0.0, SOLAR_CONSTANT);
        assert!((q - SOLAR_CONSTANT / PI).abs() < 1e-9);
    }

    #[test]
    fn june_solstice_pole_beats_equator() {
        let tilt = 23.44f64.to_radians();
        let decl = declination(EQUINOX_DAY + DAYS_PER_YEAR / 4.0, tilt);
        assert!((decl - tilt).abs() < 1e-9, "declination at solstice equals the tilt");
        let pole = daily_mean_insolation(PI / 2.0 - 0.001, decl, SOLAR_CONSTANT);
        let equator = daily_mean_insolation(0.0, decl, SOLAR_CONSTANT);
        assert!(pole > equator, "pole {pole} should exceed equator {equator}");
        // Polar value is S0 sin(decl) ~ 541 W/m2.
        assert!((pole - SOLAR_CONSTANT * tilt.sin()).abs() < 2.0);
        // And the opposite pole is in darkness.
        assert!(daily_mean_insolation(-PI / 2.0 + 0.001, decl, SOLAR_CONSTANT) < 1e-6);
    }

    #[test]
    fn annual_mean_has_earth_like_pole_to_equator_ratio() {
        let tilt = 23.44f64.to_radians();
        let annual = |lat: f64| {
            (0..360).map(|d| daily_mean_insolation(lat, declination(d as f64 + 0.5, tilt), SOLAR_CONSTANT)).sum::<f64>() / 360.0
        };
        let ratio = annual(PI / 2.0 - 0.001) / annual(0.0);
        // Textbook value for Earth is about 0.4.
        assert!((0.35..0.5).contains(&ratio), "ratio {ratio}");
    }

    #[test]
    fn zero_tilt_has_no_seasons() {
        for m in 0..12 {
            assert!(declination(month_mid_day(m), 0.0).abs() < 1e-12);
        }
    }
}
