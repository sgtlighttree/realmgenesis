//! Land water balance from monthly climate: potential evapotranspiration (Thornthwaite 1948)
//! and runoff (Budyko curve, Fu 1981 form). Used to drive erosion with runoff and to find
//! closed (endorheic) basins, where a lake's evaporation balances its inflow.

/// Annual potential evapotranspiration, mm, from 12 monthly mean temperatures (°C) at a
/// latitude (radians), for a planet with `axial_tilt_deg` (Thornthwaite 1948; the
/// high-temperature branch above 26.5 °C after Willmott et al. 1985).
pub fn thornthwaite_pet_mm(temps_c: &[f32; 12], lat_rad: f64, axial_tilt_deg: f64) -> f64 {
    let heat: f64 = temps_c.iter().map(|&t| if t > 0.0 { (t as f64 / 5.0).powf(1.514) } else { 0.0 }).sum();
    if heat <= 0.0 {
        return 0.0;
    }
    let a = 6.75e-7 * heat.powi(3) - 7.71e-5 * heat.powi(2) + 1.792e-2 * heat + 0.49239;
    let days = [31.0, 28.0, 31.0, 30.0, 31.0, 30.0, 31.0, 31.0, 30.0, 31.0, 30.0, 31.0];
    let mut total = 0.0;
    let mut doy = 0.0;
    for (m, &t) in temps_c.iter().enumerate() {
        let t = t as f64;
        let mid = doy + days[m] / 2.0;
        doy += days[m];
        if t <= 0.0 {
            continue;
        }
        let unadjusted = if t < 26.5 { 16.0 * (10.0 * t / heat).powf(a) } else { -415.85 + 32.24 * t - 0.43 * t * t };
        // Day length (hours) from the solar declination at mid-month.
        let decl = (axial_tilt_deg.to_radians()) * (std::f64::consts::TAU * (284.0 + mid) / 365.0).sin();
        let cos_h = (-lat_rad.tan() * decl.tan()).clamp(-1.0, 1.0);
        let daylength = 24.0 / std::f64::consts::PI * cos_h.acos();
        total += unadjusted.max(0.0) * (daylength / 12.0) * (days[m] / 30.0);
    }
    total
}

/// Annual runoff, mm, from precipitation and potential evapotranspiration (mm), by the
/// Budyko curve in Fu's form: ET/P = 1 + φ − (1 + φ^ω)^(1/ω), φ = PET/P. ω ≈ 2.6 is a
/// typical global value (Zhang et al. 2004). Limits: R → P when PET → 0, R → 0 when P ≪ PET.
pub fn budyko_fu_runoff_mm(p_mm: f64, pet_mm: f64, omega: f64) -> f64 {
    if p_mm <= 0.0 {
        return 0.0;
    }
    let phi = (pet_mm / p_mm).max(0.0);
    let et_over_p = 1.0 + phi - (1.0 + phi.powf(omega)).powf(1.0 / omega);
    (p_mm * (1.0 - et_over_p)).max(0.0)
}

/// Typical Fu ω for global land (Zhang et al. 2004).
pub const FU_OMEGA: f64 = 2.6;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn budyko_limits() {
        // Energy-limited (no PET): everything runs off.
        assert!((budyko_fu_runoff_mm(1000.0, 0.0, FU_OMEGA) - 1000.0).abs() < 1e-6);
        // Very arid: almost nothing runs off.
        assert!(budyko_fu_runoff_mm(100.0, 2000.0, FU_OMEGA) < 1.0);
        // Humid temperate (P 1000, PET 700): a few hundred mm, ET below both P and PET.
        let r = budyko_fu_runoff_mm(1000.0, 700.0, FU_OMEGA);
        assert!((250.0..600.0).contains(&r), "{r}");
    }

    #[test]
    fn thornthwaite_plausible() {
        // A temperate site (~45° N, seasonal 0–22 °C): PET of several hundred mm.
        let t = [0.0, 2.0, 6.0, 11.0, 16.0, 20.0, 22.0, 21.0, 17.0, 11.0, 5.0, 1.0];
        let pet = thornthwaite_pet_mm(&t, 45f64.to_radians(), 23.44);
        assert!((500.0..900.0).contains(&pet), "{pet}");
        // A hot tropical site: well over 1200 mm.
        let pet_trop = thornthwaite_pet_mm(&[27.0; 12], 5f64.to_radians(), 23.44);
        assert!(pet_trop > 1200.0, "{pet_trop}");
        // Polar: little.
        assert!(thornthwaite_pet_mm(&[-20.0; 12], 80f64.to_radians(), 23.44) < 1.0);
    }
}
