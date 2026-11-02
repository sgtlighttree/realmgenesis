//! One-dimensional seasonal energy-balance model (EBM) in latitude.
//!
//! This is the classic Budyko-Sellers / North (1975) model with a seasonal cycle, the same
//! family as climlab's `EBM_seasonal`:
//!
//! ```text
//! C dT/dt = Q(lat, t) (1 - albedo(T)) - (A + B T) + D d/dx[(1 - x^2) dT/dx]      x = sin(lat)
//! ```
//!
//! - `Q` is the daily-mean insolation ([`crate::insolation`]).
//! - `A + B T` is the outgoing long-wave radiation, linearised (North & Coakley 1979).
//! - The diffusion term stands for heat carried poleward by the atmosphere and oceans.
//! - The albedo rises when `T` drops below freezing (ice and snow feedback).
//! - `C` is the heat capacity of the surface layer, in J/m2/K. It is the only difference
//!   between "land-like" and "ocean-like" behaviour: a deep ocean mixed layer answers slowly
//!   (small, late seasonal swing) and thin land answers quickly (large swing).
//!
//! Tier A runs the model twice, once with ocean heat capacity and once with land heat
//! capacity, and then blends the two per cell by continentality. The solver uses backward Euler
//! (unconditionally stable, so thin "land" layers are fine) and repeats the year until the
//! monthly means stop changing.

use crate::insolation::{SOLAR_CONSTANT, daily_mean_insolation, declination};
use crate::types::PlanetParams;

/// Number of equal-area latitude bands (equal steps in sin(lat)).
const BANDS: usize = 90;
/// Time steps per year (about 3 days each).
const STEPS_PER_YEAR: usize = 120;
const MAX_YEARS: usize = 80;

/// Outgoing long-wave radiation `A + B T` (W/m2, W/m2/K). Roughly North & Coakley 1979.
const OLR_A: f64 = 203.3;
const OLR_B: f64 = 2.09;
/// Meridional heat diffusion for Earth (W/m2/K).
const DIFFUSION_EARTH: f64 = 0.56;
/// Surface albedo with no ice, at the equator and added toward the poles (clouds, bright
/// land).
const ALBEDO_EQUATOR: f64 = 0.275;
const ALBEDO_POLE_EXTRA: f64 = 0.04;
/// Additional polar brightening that only matters close to the poles (x^4 dependence).
const ALBEDO_POLE_EXTRA4: f64 = 0.09;
/// Extra albedo of a fully frozen surface (ramps in from -5 C to -17 C).
const ALBEDO_ICE_EXTRA: f64 = 0.28;

/// Heat capacity of a well-mixed ocean layer plus the air above it, J/m2/K.
pub const HEAT_CAPACITY_OCEAN: f64 = 8.0e7;
/// Heat capacity of land (thin soil layer plus the air above it), J/m2/K.
pub const HEAT_CAPACITY_LAND: f64 = 1.2e7;

/// Zonal-mean (averaged around a latitude circle) monthly sea-level temperatures.
pub struct ZonalClimate {
    /// `sin(lat)` of each band centre, south to north.
    x: Vec<f64>,
    /// Monthly mean temperature of each band, degrees C.
    temp: Vec<[f64; 12]>,
}

impl ZonalClimate {
    /// Monthly mean temperature (C) at `lat_rad`, linearly interpolated between bands in sin(lat).
    pub fn at(&self, lat_rad: f64, month: usize) -> f64 {
        let dx = 2.0 / BANDS as f64;
        let f = ((lat_rad.sin() - self.x[0]) / dx).clamp(0.0, (BANDS - 1) as f64);
        let i0 = (f.floor() as usize).min(BANDS - 2);
        let t = f - i0 as f64;
        self.temp[i0][month] * (1.0 - t) + self.temp[i0 + 1][month] * t
    }
}

/// Solve the seasonal EBM for a surface layer of the given heat capacity (J/m2/K).
pub fn solve_zonal_ebm(planet: &PlanetParams, heat_capacity: f64) -> ZonalClimate {
    let tilt = planet.axial_tilt_deg.to_radians();
    let s0 = SOLAR_CONSTANT * planet.solar_constant_rel;

    // Faster rotation = narrower, weaker-mixing circulation cells = less meridional transport.
    let rot = (planet.rotation_hours / 24.0).clamp(0.2, 8.0);
    let diffusion = DIFFUSION_EARTH * rot.sqrt().clamp(0.5, 2.0);

    let dx = 2.0 / BANDS as f64;
    let x_centre: Vec<f64> = (0..BANDS).map(|i| -1.0 + (i as f64 + 0.5) * dx).collect();
    // Diffusion conductance across the interface between band i-1 and i (W/m2/K).
    let mut k = vec![0.0; BANDS + 1];
    for (i, ki) in k.iter_mut().enumerate().take(BANDS).skip(1) {
        let xe = -1.0 + i as f64 * dx;
        *ki = diffusion * (1.0 - xe * xe) / (dx * dx);
    }
    let albedo_base: Vec<f64> =
        x_centre.iter().map(|x| ALBEDO_EQUATOR + ALBEDO_POLE_EXTRA * x * x + ALBEDO_POLE_EXTRA4 * x.powi(4)).collect();

    // Insolation table [step][band].
    let dt = DAYS_YEAR * 86400.0 / STEPS_PER_YEAR as f64;
    let q: Vec<Vec<f64>> = (0..STEPS_PER_YEAR)
        .map(|s| {
            let decl = declination((s as f64 + 0.5) * DAYS_YEAR / STEPS_PER_YEAR as f64, tilt);
            x_centre.iter().map(|x| daily_mean_insolation(x.asin(), decl, s0)).collect()
        })
        .collect();

    let steps_per_month = STEPS_PER_YEAR / 12;
    let c_dt = heat_capacity / dt;
    // Start from a smooth guess.
    let mut temp: Vec<f64> = x_centre.iter().map(|x| 28.0 - 55.0 * x * x).collect();
    let mut monthly = vec![[0.0f64; 12]; BANDS];
    let mut prev_monthly = vec![[f64::MAX; 12]; BANDS];

    // Scratch for the tridiagonal solve.
    let (mut cp, mut dp) = (vec![0.0; BANDS], vec![0.0; BANDS]);
    for _year in 0..MAX_YEARS {
        for row in monthly.iter_mut() {
            *row = [0.0; 12];
        }
        for (step, q_step) in q.iter().enumerate() {
            // Build and solve (c/dt + B + k_i + k_{i+1}) T_i - k_i T_{i-1} - k_{i+1} T_{i+1} = rhs.
            for i in 0..BANDS {
                let ice = ((-temp[i] - 5.0) / 12.0).clamp(0.0, 1.0);
                let albedo = albedo_base[i] + ALBEDO_ICE_EXTRA * ice;
                let rhs = c_dt * temp[i] + q_step[i] * (1.0 - albedo) - OLR_A;
                let lower = -k[i]; // coupling to i-1
                let upper = -k[i + 1]; // coupling to i+1
                let diag = c_dt + OLR_B + k[i] + k[i + 1];
                if i == 0 {
                    cp[0] = upper / diag;
                    dp[0] = rhs / diag;
                } else {
                    let m = diag - lower * cp[i - 1];
                    cp[i] = upper / m;
                    dp[i] = (rhs - lower * dp[i - 1]) / m;
                }
            }
            temp[BANDS - 1] = dp[BANDS - 1];
            for i in (0..BANDS - 1).rev() {
                temp[i] = dp[i] - cp[i] * temp[i + 1];
            }
            let month = (step / steps_per_month).min(11);
            for (i, t) in temp.iter().enumerate() {
                monthly[i][month] += t / steps_per_month as f64;
            }
        }
        let change = monthly
            .iter()
            .zip(&prev_monthly)
            .flat_map(|(a, b)| a.iter().zip(b.iter()).map(|(x, y)| (x - y).abs()))
            .fold(0.0, f64::max);
        prev_monthly.clone_from(&monthly);
        if change < 0.02 {
            break;
        }
    }
    ZonalClimate { x: x_centre, temp: monthly }
}

const DAYS_YEAR: f64 = 365.0;

#[cfg(test)]
mod tests {
    use super::*;

    fn annual(z: &ZonalClimate, lat_deg: f64) -> f64 {
        (0..12).map(|m| z.at(lat_deg.to_radians(), m)).sum::<f64>() / 12.0
    }

    #[test]
    fn earth_zonal_means_are_plausible() {
        let p = PlanetParams::default();
        let ocean = solve_zonal_ebm(&p, HEAT_CAPACITY_OCEAN);
        let land = solve_zonal_ebm(&p, HEAT_CAPACITY_LAND);
        for (name, z) in [("ocean", &ocean), ("land", &land)] {
            let line: Vec<String> =
                [0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 89.0].iter().map(|&l| format!("{l}:{:.1}", annual(z, l))).collect();
            println!("{name} annual {}", line.join(" "));
            let jj: Vec<String> = [0.0, 30.0, 45.0, 60.0, 80.0]
                .iter()
                .map(|&l: &f64| format!("{l}: {:.1}/{:.1}", z.at(l.to_radians(), 0), z.at(l.to_radians(), 6)))
                .collect();
            println!("{name} Jan/Jul {}", jj.join(" "));
        }
        let eq = annual(&ocean, 0.0);
        let pole = annual(&ocean, 89.0);
        assert!((24.0..30.0).contains(&eq), "equator {eq}");
        assert!((-30.0..-10.0).contains(&pole), "pole {pole}");
        // Seasonal swing of land is larger than ocean's at mid latitudes.
        let swing = |z: &ZonalClimate| z.at(45f64.to_radians(), 6) - z.at(45f64.to_radians(), 0);
        assert!(swing(&land) > 1.5 * swing(&ocean));
    }
}
