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
use rg_noise::{salt, Fbm};
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
    /// Unit ice-edge noise per pixel (see `polar_noise`). It depends only on the position, so
    /// it's computed once here instead of in every bake. Empty rows near the equator are zero.
    polar: Vec<[f32; 2]>,
}

impl Equirect {
    pub fn new(world: &World, width: u32, height: u32) -> Self {
        let n = (width * height) as usize;
        let seed = world.params.tectonics.seed;
        let mut heights = vec![0f32; n];
        par_rows(&mut heights, width as usize, |y, row| {
            let lat = lat_of_row(y as u32, height);
            for (x, v) in row.iter_mut().enumerate() {
                *v = world.tectonics.height_m.sample_drained(dir(lat, lon_of_col(x as u32, width))) as f32;
            }
        });
        let mut polar = vec![[0f32; 2]; n];
        par_rows(&mut polar, width as usize, |y, row| {
            let lat = lat_of_row(y as u32, height);
            if lat.abs() >= POLAR_MIN_LAT {
                for (x, v) in row.iter_mut().enumerate() {
                    *v = polar_unit(seed, dir(lat, lon_of_col(x as u32, width)));
                }
            }
        });
        Self { width, height, heights, polar }
    }

    /// Bake one layer to RGBA8 sRGB bytes.
    pub fn bake(&self, world: &World, layer: Layer) -> Vec<u8> {
        let (w, h) = (self.width, self.height);
        let lakes = world.hydrology.as_ref().map(|hy| {
            let d = &hy.lake_depth_m;
            visible_lakes(d, min_lake_nodes(d.res, world.params.tectonics.radius_m))
        });
        let mut out = vec![0u8; (w * h * 4) as usize];
        par_rows(&mut out, w as usize * 4, |y, row| {
            let lat = lat_of_row(y as u32, h);
            for (x, px) in row.chunks_exact_mut(4).enumerate() {
                let lon = lon_of_col(x as u32, w);
                let i = y * w as usize + x;
                let hf = self.heights[i];
                let rgb: Rgb = match layer {
                    Layer::Elevation => palette::elevation(hf),
                    Layer::Plates => plates_color(world, lat, lon, hf),
                    _ => {
                        let s = ClimateSample::at_with_noise(world, lat, lon, hf, self.polar[i]);
                        match layer {
                            Layer::Satellite => satellite(&s, hf, lat),
                            Layer::Koppen => {
                                if hf < 0.0 {
                                    if s.sea_ice() { KOPPEN_SEA_ICE } else { KOPPEN_OCEAN }
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
                // Lakes: standing water in closed basins (from the erosion stage).
                let rgb = match (&lakes, layer) {
                    (Some(lk), Layer::Satellite | Layer::Elevation | Layer::PrecipAnnual) if hf > 0.0 => {
                        let depth = lk.sample(dir(lat, lon));
                        if depth > 2.0 { lake_color(layer) } else { rgb }
                    }
                    _ => rgb,
                };
                px.copy_from_slice(&[rgb[0], rgb[1], rgb[2], 255]);
            }
        });
        if let (Some(hy), Layer::Satellite | Layer::Elevation | Layer::PrecipAnnual) = (&world.hydrology, layer) {
            draw_rivers(&mut out, w, h, &hy.rivers);
        }
        out
    }
}

/// Lakes smaller than this many macro nodes aren't drawn: single- and few-node pits left by
/// the erosion solve read as speckle ("dotted lakes", Matt, s02). Display only: the
/// committed `lake_depth_m` is unchanged.
/// The floor is an **area** (s04): as a node count (4 before s04) the same lake was hidden
/// at Normal and shown at Ultra, where it covers ~10× the nodes. 4 Normal cells ≈ 3,900 km².
/// Shared with refinement, which turns lakes below it into playas.
const MIN_LAKE_KM2: f64 = rg_erosion::refine::MIN_VISIBLE_LAKE_KM2;

/// `MIN_LAKE_KM2` in nodes of a cube grid with `res` cells per face edge (mean cell size).
fn min_lake_nodes(res: u32, radius_m: f64) -> usize {
    let cell_km = std::f64::consts::FRAC_PI_2 * radius_m * 1e-3 / res.max(1) as f64;
    (MIN_LAKE_KM2 / (cell_km * cell_km)).ceil().max(1.0) as usize
}

/// `depth` with every lake of fewer than `min_nodes` nodes set to 0. Lakes are flooded per
/// cube face (4-neighbours), so a lake cut by a cube edge is judged per side.
fn visible_lakes(depth: &rg_geo::CubeField<f32>, min_nodes: usize) -> rg_geo::CubeField<f32> {
    let mut out = depth.clone();
    let n = depth.samples_per_edge() as i64;
    let mut seen = vec![false; depth.data.len()];
    let mut stack: Vec<(i64, i64)> = Vec::new();
    let mut lake: Vec<usize> = Vec::new();
    for face in rg_geo::CubeFace::ALL {
        for j in 0..n {
            for i in 0..n {
                let k = depth.index(face, i as u32, j as u32);
                if seen[k] || depth.data[k] <= 0.0 {
                    continue;
                }
                seen[k] = true;
                lake.clear();
                stack.push((i, j));
                while let Some((a, b)) = stack.pop() {
                    lake.push(depth.index(face, a as u32, b as u32));
                    for (x, y) in [(a + 1, b), (a - 1, b), (a, b + 1), (a, b - 1)] {
                        if x < 0 || y < 0 || x >= n || y >= n {
                            continue;
                        }
                        let kk = depth.index(face, x as u32, y as u32);
                        if !seen[kk] && depth.data[kk] > 0.0 {
                            seen[kk] = true;
                            stack.push((x, y));
                        }
                    }
                }
                if lake.len() < min_nodes {
                    for &kk in &lake {
                        out.data[kk] = 0.0;
                    }
                }
            }
        }
    }
    out
}

/// Run `f(row_index, row)` over the rows of a row-major image (`row_len` items per row).
/// On native targets the rows are shared out to threads in small blocks; every row is
/// computed independently, so the result is the same either way. wasm32 has no threads and
/// runs the rows in order.
fn par_rows<T: Send>(data: &mut [T], row_len: usize, f: impl Fn(usize, &mut [T]) + Sync) {
    #[cfg(target_arch = "wasm32")]
    {
        for (y, row) in data.chunks_mut(row_len).enumerate() {
            f(y, row);
        }
    }
    #[cfg(not(target_arch = "wasm32"))]
    {
        const ROWS_PER_BLOCK: usize = 4;
        let threads = std::thread::available_parallelism().map_or(1, |t| t.get()).min(16);
        let blocks = std::sync::Mutex::new(data.chunks_mut(row_len * ROWS_PER_BLOCK).enumerate());
        let work = || {
            loop {
                let next = blocks.lock().unwrap().next();
                let Some((b, block)) = next else { break };
                for (k, row) in block.chunks_mut(row_len).enumerate() {
                    f(b * ROWS_PER_BLOCK + k, row);
                }
            }
        };
        std::thread::scope(|s| {
            for _ in 1..threads {
                s.spawn(work);
            }
            work();
        });
    }
}

/// Monthly climate at one pixel after downscaling.
pub struct ClimateSample {
    pub t: [f32; 12],
    pub p: [f32; 12],
    /// Extra °C of fine, ragged noise that only the ice and snow thresholds use (floe edges
    /// and snow patches). It is not added to `t`, so temperature readouts stay smooth.
    pub ice_fine: f32,
}

/// Noise stream for polar ice and snow edges (see `polar_noise`).
const STREAM_ICE: u32 = 0x1CE;

/// Closer to the equator than this the ice-edge noise is never used.
const POLAR_MIN_LAT: f64 = 30.0 * PI / 180.0;

/// Tier A climate is almost zonal over open ocean, so the sea-ice and snow lines would be
/// perfect latitude circles. This gives them a natural, deterministic wobble. It returns
/// `[shift, fine]` in °C, before fading by temperature:
/// * `shift` is a smooth, large-scale change to the temperature (several °C, like warm and
///   cold currents bending the isotherms).
/// * `fine` is smaller and much finer, for ragged edges.
fn polar_unit(seed: u64, d: DVec3) -> [f32; 2] {
    let s = salt(seed, STREAM_ICE);
    // Wavelengths from ~10000 km down to ~1000 km.
    const BROAD: Fbm = Fbm { octaves: 4, frequency: 4.0, lacunarity: 2.2, gain: 0.62 };
    // Wavelengths ~400 km down to ~100 km.
    const FINE: Fbm = Fbm { octaves: 3, frequency: 100.0, lacunarity: 2.0, gain: 0.6 };
    [11.0 * BROAD.sample(d, s) as f32, 3.5 * FINE.sample(d, s.wrapping_add(0x100)) as f32]
}

/// How much of the ice-edge noise applies at a place with this annual mean temperature:
/// 1 below -5 °C, 0 above +10 °C. Warm places are never affected.
fn polar_weight(mean_t: f32) -> f32 {
    let x = ((10.0 - mean_t) / 15.0).clamp(0.0, 1.0);
    x * x * (3.0 - 2.0 * x)
}

impl ClimateSample {
    pub fn at(world: &World, lat: f64, lon: f64, fine_h: f32) -> Self {
        // The same noise a baked pixel uses (hover readouts must agree with the map).
        let unit = if lat.abs() >= POLAR_MIN_LAT {
            polar_unit(world.params.tectonics.seed, dir(lat, lon))
        } else {
            [0.0; 2]
        };
        Self::at_with_noise(world, lat, lon, fine_h, unit)
    }

    /// Like [`ClimateSample::at`], with the unit ice-edge noise already computed.
    fn at_with_noise(world: &World, lat: f64, lon: f64, fine_h: f32, unit_noise: [f32; 2]) -> Self {
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
        let mean_t = t.iter().sum::<f32>() / 12.0;
        let w = polar_weight(mean_t);
        for v in &mut t {
            *v += w * unit_noise[0];
        }
        Self { t, p, ice_fine: w * unit_noise[1] }
    }

    pub fn koppen(&self, lat: f64) -> Koppen {
        Koppen::classify(&self.t, &self.p, lat >= 0.0)
    }

    pub fn sea_ice(&self) -> bool {
        self.t.iter().copied().fold(f32::MAX, f32::min) + self.ice_fine < -6.0
    }
}

/// Köppen map colours for open water (sea ice and ice-free).
pub const KOPPEN_SEA_ICE: Rgb = [214, 222, 230];
pub const KOPPEN_OCEAN: Rgb = [196, 210, 222];

/// Base land colour of the satellite layer for a Köppen class, before the rainfall tint,
/// rock and snow are applied. The legend reads this too.
pub fn satellite_base(k: Koppen) -> Rgb {
    use Koppen::*;
    match k {
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
    }
}

/// Satellite-like albedo: biome colour from the downscaled Köppen class, nudged by
/// rainfall, with snow, ice and bare rock.
fn satellite(s: &ClimateSample, h: f32, lat: f64) -> Rgb {
    if h < 0.0 {
        let water = palette::sat_ocean(h);
        // The fine noise shifts the thresholds, which makes ragged floe edges.
        let coldest = s.t.iter().copied().fold(f32::MAX, f32::min) + s.ice_fine;
        // Seasonal sea ice fades in below about -2 °C coldest month; permanent pack ice
        // when even the warmest month is freezing.
        let warmest = s.t.iter().copied().fold(f32::MIN, f32::max) + s.ice_fine;
        let ice = if warmest < -1.0 { 0.92 } else { ((-coldest - 2.0) / 14.0).clamp(0.0, 0.55) };
        return palette::mix(water, [222, 230, 238], ice);
    }
    let annual_p: f32 = s.p.iter().sum();
    let warmest = s.t.iter().copied().fold(f32::MIN, f32::max);
    let mean_t = s.t.iter().sum::<f32>() / 12.0;
    // Snow patches get half the fine noise the sea ice gets.
    let snow_fine = 0.5 * s.ice_fine;
    let k = s.koppen(lat);
    let base = satellite_base(k);
    // Wetter = greener within a class; dry edges fade towards the steppe tone.
    let wet = ((annual_p - 400.0) / 1600.0).clamp(-1.0, 1.0);
    let mut c = if wet < 0.0 { palette::mix(base, [150, 136, 96], -wet * 0.35) } else { palette::mix(base, [36, 72, 34], wet * 0.25) };
    // Bare rock above the tree line on cold, high ground.
    if h > 2200.0 && warmest < 12.0 {
        c = palette::mix(c, [120, 112, 106], ((h - 2200.0) / 1200.0).clamp(0.0, 0.8));
    }
    // Snow cover where the year is cold enough.
    if !matches!(k, Koppen::EF) {
        let snow = ((-(mean_t + snow_fine) - 4.0) / 10.0)
            .clamp(0.0, 0.85)
            .max(((-(warmest + snow_fine) + 1.5) / 4.0).clamp(0.0, 1.0));
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

fn lake_color(layer: Layer) -> Rgb {
    match layer {
        Layer::Satellite => [44, 92, 124],
        _ => [120, 170, 214],
    }
}

/// The tint the renderer applies over river coverage (stored in the layer's alpha).
pub fn river_color(layer: Layer) -> Rgb {
    match layer {
        Layer::Satellite => [46, 96, 132],
        _ => [60, 110, 200],
    }
}

/// Rasterise river reaches into the **alpha** of an equirect RGBA image (anti-aliased
/// coverage: alpha = 255 − coverage), leaving the colour untouched. The renderer tints by
/// coverage and fades rivers out when zoomed in, where a texel spans kilometres and the 3D
/// valleys show the river (s01: baked blue bands hid the terrain close up). Width grows
/// with the square root of drainage area, in pixels relative to a 2048-wide image.
pub fn draw_rivers(img: &mut [u8], w: u32, h: u32, rivers: &[rg_erosion::RiverSegment]) {
    let scale = w as f64 / 2048.0;
    let to_px = |d: DVec3| -> (f64, f64) {
        let lat = d.z.clamp(-1.0, 1.0).asin();
        let lon = d.y.atan2(d.x);
        ((lon + PI) / (2.0 * PI) * w as f64, (PI / 2.0 - lat) / PI * h as f64)
    };
    for r in rivers {
        let width = ((r.area_km2 as f64 / 2.0e4).sqrt() * 0.45 * scale).clamp(0.5 * scale, 2.6 * scale);
        let (x0, y0) = to_px(r.from);
        let (mut x1, y1) = to_px(r.to);
        if (x1 - x0).abs() > w as f64 / 2.0 {
            x1 += if x1 < x0 { w as f64 } else { -(w as f64) };
        }
        let half = width * 0.5;
        let (minx, maxx) = (x0.min(x1) - half - 1.0, x0.max(x1) + half + 1.0);
        let (miny, maxy) = (y0.min(y1) - half - 1.0, y0.max(y1) + half + 1.0);
        let (dx, dy) = (x1 - x0, y1 - y0);
        let len2 = (dx * dx + dy * dy).max(1e-9);
        for py in (miny.floor() as i64)..=(maxy.ceil() as i64) {
            if py < 0 || py >= h as i64 {
                continue;
            }
            for px in (minx.floor() as i64)..=(maxx.ceil() as i64) {
                let (cx, cy) = (px as f64 + 0.5, py as f64 + 0.5);
                let t = (((cx - x0) * dx + (cy - y0) * dy) / len2).clamp(0.0, 1.0);
                let (qx, qy) = (x0 + dx * t - cx, y0 + dy * t - cy);
                let d = (qx * qx + qy * qy).sqrt();
                let cover = (half + 0.5 - d).clamp(0.0, 1.0);
                if cover <= 0.0 {
                    continue;
                }
                let xi = px.rem_euclid(w as i64) as u32;
                let k = ((py as u32 * w + xi) * 4) as usize;
                let a = 255.0 - (cover * 255.0).round();
                img[k + 3] = img[k + 3].min(a as u8);
            }
        }
    }
}

#[cfg(test)]
mod lake_tests {
    use super::visible_lakes;
    use rg_geo::{CubeFace, CubeField};

    #[test]
    fn small_lakes_are_hidden_and_big_ones_kept() {
        let mut f = CubeField::<f32>::from_fn(16, |_| 0.0);
        let face = CubeFace::ALL[0];
        let one = f.index(face, 2, 2);
        f.data[one] = 30.0;
        for j in 8..11 {
            for i in 8..11 {
                let k = f.index(face, i, j);
                f.data[k] = 50.0;
            }
        }
        let v = visible_lakes(&f, 4);
        assert_eq!(v.data[one], 0.0, "a one-node pit is hidden");
        assert_eq!(v.data[f.index(face, 9, 9)], 50.0, "a 9-node lake is kept");
    }
}
