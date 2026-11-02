//! Display layers baked to equirectangular images.
//!
//! Climate runs at 1°, which is coarse. For display it's **downscaled per pixel**
//! (ARCHITECTURE §7): monthly temperature and precipitation are interpolated, temperature
//! is corrected by the lapse rate for the difference between the pixel's true height and
//! the height the climate model saw, and Köppen is re-classified per pixel. That's what
//! gives sharp climate edges along coasts and mountain ranges.

use crate::palette::{self, Rgb};
use crate::World;
use glam::DVec3;
use rg_climate::Koppen;
use std::f64::consts::PI;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Layer {
    Satellite,
    Elevation,
    Plates,
    Koppen,
    TempAnnual,
    TempJanuary,
    TempJuly,
    PrecipAnnual,
}

impl Layer {
    pub const ALL: [Layer; 8] = [
        Layer::Satellite,
        Layer::Elevation,
        Layer::Plates,
        Layer::Koppen,
        Layer::TempAnnual,
        Layer::TempJanuary,
        Layer::TempJuly,
        Layer::PrecipAnnual,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Layer::Satellite => "Satellite",
            Layer::Elevation => "Elevation",
            Layer::Plates => "Tectonic plates",
            Layer::Koppen => "Köppen climate",
            Layer::TempAnnual => "Temperature · year",
            Layer::TempJanuary => "Temperature · January",
            Layer::TempJuly => "Temperature · July",
            Layer::PrecipAnnual => "Precipitation · year",
        }
    }
}

/// Lapse rate used for downscaling, °C per metre.
const LAPSE: f32 = 0.0065;

/// An equirectangular raster of the world at display resolution, with fine heights
/// cached (they're needed by every layer).
pub struct Equirect {
    pub width: u32,
    pub height: u32,
    /// Metres above sea level, row 0 = north.
    pub heights: Vec<f32>,
}

impl Equirect {
    pub fn new(world: &World, width: u32, height: u32) -> Self {
        let mut heights = Vec::with_capacity((width * height) as usize);
        for y in 0..height {
            let lat = lat_of_row(y, height);
            for x in 0..width {
                let lon = lon_of_col(x, width);
                heights.push(world.tectonics.height_m.sample(dir(lat, lon)));
            }
        }
        Self { width, height, heights }
    }

    /// Bake one layer to RGBA8 sRGB bytes.
    pub fn bake(&self, world: &World, layer: Layer) -> Vec<u8> {
        let (w, h) = (self.width, self.height);
        let mut out = vec![0u8; (w * h * 4) as usize];
        for y in 0..h {
            let lat = lat_of_row(y, h);
            for x in 0..w {
                let lon = lon_of_col(x, w);
                let i = (y * w + x) as usize;
                let hf = self.heights[i];
                let rgb: Rgb = match layer {
                    Layer::Elevation => palette::elevation(hf),
                    Layer::Plates => plates_color(world, lat, lon, hf),
                    _ => {
                        let s = ClimateSample::at(world, lat, lon, hf);
                        match layer {
                            Layer::Satellite => satellite(&s, hf, lat),
                            Layer::Koppen => {
                                if hf < 0.0 {
                                    if s.sea_ice() { [214, 222, 230] } else { [196, 210, 222] }
                                } else {
                                    s.koppen(lat).color()
                                }
                            }
                            Layer::TempAnnual => palette::temperature(s.t.iter().sum::<f32>() / 12.0),
                            Layer::TempJanuary => palette::temperature(s.t[0]),
                            Layer::TempJuly => palette::temperature(s.t[6]),
                            Layer::PrecipAnnual => palette::precipitation(s.p.iter().sum()),
                            Layer::Elevation | Layer::Plates => unreachable!(),
                        }
                    }
                };
                out[i * 4..i * 4 + 4].copy_from_slice(&[rgb[0], rgb[1], rgb[2], 255]);
            }
        }
        out
    }
}

/// Monthly climate at one pixel after downscaling.
pub struct ClimateSample {
    pub t: [f32; 12],
    pub p: [f32; 12],
}

impl ClimateSample {
    pub fn at(world: &World, lat: f64, lon: f64, fine_h: f32) -> Self {
        let c = &world.climate;
        let (w, h) = (c.width as i64, c.height as i64);
        // Cell-centred grid: row 0 centre at lat 90 - Δ/2, col 0 centre at lon -180 + Δ/2.
        let fx = (lon + PI) / (2.0 * PI) * w as f64 - 0.5;
        let fy = (PI / 2.0 - lat) / PI * h as f64 - 0.5;
        let (x0, y0) = (fx.floor() as i64, fy.floor() as i64);
        let (tx, ty) = ((fx - x0 as f64) as f32, (fy - y0 as f64) as f32);
        let idx = |x: i64, y: i64| (y.clamp(0, h - 1) * w + x.rem_euclid(w)) as usize;
        let ids = [idx(x0, y0), idx(x0 + 1, y0), idx(x0, y0 + 1), idx(x0 + 1, y0 + 1)];
        let ws = [(1.0 - tx) * (1.0 - ty), tx * (1.0 - ty), (1.0 - tx) * ty, tx * ty];
        let mut t = [0f32; 12];
        let mut p = [0f32; 12];
        let mut coarse_h = 0f32;
        for (&k, &wt) in ids.iter().zip(&ws) {
            for m in 0..12 {
                t[m] += c.temp_c[k][m] * wt;
                p[m] += c.precip_mm[k][m] * wt;
            }
            coarse_h += world.climate_heights.data[k].max(0.0) * wt;
        }
        if fine_h > 0.0 {
            let dt = -LAPSE * (fine_h - coarse_h);
            for v in &mut t {
                *v += dt;
            }
        }
        Self { t, p }
    }

    pub fn koppen(&self, lat: f64) -> Koppen {
        Koppen::classify(&self.t, &self.p, lat >= 0.0)
    }

    pub fn sea_ice(&self) -> bool {
        self.t.iter().copied().fold(f32::MAX, f32::min) < -6.0
    }
}

/// Satellite-like albedo: biome colour from the downscaled Köppen class, nudged by
/// rainfall, with snow, ice and bare rock.
fn satellite(s: &ClimateSample, h: f32, lat: f64) -> Rgb {
    if h < 0.0 {
        let water = palette::sat_ocean(h);
        let coldest = s.t.iter().copied().fold(f32::MAX, f32::min);
        // Seasonal sea ice fades in below about -2 °C coldest month; permanent pack ice
        // when even the warmest month is freezing.
        let warmest = s.t.iter().copied().fold(f32::MIN, f32::max);
        let ice = if warmest < -1.0 { 0.92 } else { ((-coldest - 2.0) / 14.0).clamp(0.0, 0.55) };
        return palette::mix(water, [222, 230, 238], ice);
    }
    let annual_p: f32 = s.p.iter().sum();
    let warmest = s.t.iter().copied().fold(f32::MIN, f32::max);
    let mean_t = s.t.iter().sum::<f32>() / 12.0;
    use Koppen::*;
    let k = s.koppen(lat);
    let base: Rgb = match k {
        Af | Am => [34, 74, 32],
        Aw => [104, 112, 58],
        BWh => [206, 174, 124],
        BWk => [170, 150, 116],
        BSh => [158, 140, 88],
        BSk => [140, 132, 92],
        Csa | Csb | Csc => [112, 116, 70],
        Cfa | Cwa => [62, 96, 44],
        Cfb | Cfc | Cwb | Cwc => [58, 94, 48],
        Dfa | Dfb | Dsa | Dsb | Dwa | Dwb => [54, 84, 44],
        Dfc | Dfd | Dsc | Dsd | Dwc | Dwd => [40, 62, 40],
        ET => [124, 120, 100],
        EF => [236, 240, 246],
    };
    // Wetter = greener within a class; dry edges fade towards the steppe tone.
    let wet = ((annual_p - 400.0) / 1600.0).clamp(-1.0, 1.0);
    let mut c = if wet < 0.0 { palette::mix(base, [150, 136, 96], -wet * 0.35) } else { palette::mix(base, [36, 72, 34], wet * 0.25) };
    // Bare rock above the tree line on cold, high ground.
    if h > 2200.0 && warmest < 12.0 {
        c = palette::mix(c, [120, 112, 106], ((h - 2200.0) / 1200.0).clamp(0.0, 0.8));
    }
    // Snow cover where the year is cold enough.
    if !matches!(k, EF) {
        let snow = ((-mean_t - 4.0) / 10.0).clamp(0.0, 0.85).max(((-warmest + 1.5) / 4.0).clamp(0.0, 1.0));
        c = palette::mix(c, [234, 238, 244], snow);
    }
    c
}

fn plates_color(world: &World, lat: f64, lon: f64, h: f32) -> Rgb {
    let d = dir(lat, lon);
    let t = &world.tectonics;
    let id = t.plate_id.nearest(d);
    let mut c = palette::plate(id);
    if t.crust.nearest(d) == 0 {
        c = palette::mix(c, [20, 30, 50], 0.35);
    }
    if h >= 0.0 {
        c = palette::mix(c, [255, 255, 255], 0.18);
    }
    match t.boundary.nearest(d) {
        1 => [214, 52, 52],  // convergent
        2 => [56, 120, 230], // divergent
        3 => [60, 190, 90],  // transform
        _ => c,
    }
}

pub fn lat_of_row(y: u32, h: u32) -> f64 {
    PI / 2.0 - (y as f64 + 0.5) / h as f64 * PI
}

pub fn lon_of_col(x: u32, w: u32) -> f64 {
    -PI + (x as f64 + 0.5) / w as f64 * 2.0 * PI
}

pub fn dir(lat: f64, lon: f64) -> DVec3 {
    let (sl, cl) = lat.sin_cos();
    let (so, co) = lon.sin_cos();
    DVec3::new(cl * co, cl * so, sl)
}
