//! Input and output types shared by Tier A and (later) Tier B.

use crate::grid::LatLonGrid;

/// Physical description of the planet. `Default` is Earth.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanetParams {
    /// Mean radius in metres. Only used to turn grid cells into kilometres.
    pub radius_m: f64,
    /// Obliquity (axial tilt) in degrees. Earth 23.44.
    pub axial_tilt_deg: f64,
    /// Length of one rotation (a "day") in hours. Earth 24. Faster = narrower circulation cells.
    pub rotation_hours: f64,
    /// Stellar flux relative to Earth (1.0 = 1361 W/m2).
    pub solar_constant_rel: f64,
    /// Artistic knob: added to every temperature at the end, in degrees C.
    pub global_temp_offset_c: f64,
    /// Artistic knob: multiplies every precipitation value.
    pub global_precip_mult: f64,
}

impl Default for PlanetParams {
    fn default() -> Self {
        Self {
            radius_m: 6_371_000.0,
            axial_tilt_deg: 23.44,
            rotation_hours: 24.0,
            solar_constant_rel: 1.0,
            global_temp_offset_c: 0.0,
            global_precip_mult: 1.0,
        }
    }
}

/// What the climate model needs to know about the surface.
#[derive(Clone, Debug)]
pub struct ClimateInput {
    /// Elevation in metres relative to sea level. Negative = ocean. Any resolution; the
    /// climate fields come out on the same grid.
    pub height_m: LatLonGrid<f32>,
}

/// Monthly climate on the input grid. Month 0 is January; "January" is just the first twelfth
/// of the year, so on a tilted planet it is simply the month starting at the calendar origin
/// (the model puts the northern spring equinox at about day 80, like Earth).
#[derive(Clone, Debug)]
pub struct ClimateFields {
    pub width: u32,
    pub height: u32,
    /// Monthly mean 2 m air temperature per cell, degrees C (`y * width + x` order).
    pub temp_c: Vec<[f32; 12]>,
    /// Monthly precipitation totals per cell, mm per month (oceans included).
    pub precip_mm: Vec<[f32; 12]>,
    /// Köppen-Geiger class as `Koppen as u8`; 0 for ocean.
    pub koppen: LatLonGrid<u8>,
    /// Prevailing surface wind, eastward component, m/s per month.
    pub wind_u: Vec<[f32; 12]>,
    /// Prevailing surface wind, northward component, m/s per month.
    pub wind_v: Vec<[f32; 12]>,
}

impl ClimateFields {
    #[inline]
    pub fn index(&self, x: u32, y: u32) -> usize {
        y as usize * self.width as usize + x as usize
    }

    /// One month of temperature as a grid (handy for display layers).
    pub fn temp_grid(&self, month: usize) -> LatLonGrid<f32> {
        LatLonGrid::from_vec(self.width, self.height, self.temp_c.iter().map(|t| t[month]).collect())
    }

    /// One month of precipitation as a grid.
    pub fn precip_grid(&self, month: usize) -> LatLonGrid<f32> {
        LatLonGrid::from_vec(self.width, self.height, self.precip_mm.iter().map(|p| p[month]).collect())
    }

    /// Annual mean temperature per cell.
    pub fn annual_temp_grid(&self) -> LatLonGrid<f32> {
        LatLonGrid::from_vec(self.width, self.height, self.temp_c.iter().map(|t| t.iter().sum::<f32>() / 12.0).collect())
    }

    /// Annual precipitation total per cell.
    pub fn annual_precip_grid(&self) -> LatLonGrid<f32> {
        LatLonGrid::from_vec(self.width, self.height, self.precip_mm.iter().map(|p| p.iter().sum()).collect())
    }
}
