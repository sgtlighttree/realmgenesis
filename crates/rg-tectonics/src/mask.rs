//! R3 step 0: the land mask. Coasts become an **input** (design §0).
//!
//! In strict mode the land/sea split comes from here, not from the heights: a procedural
//! mask (today's continent noise thresholded to the land fraction), an imported raster placed
//! on the sphere, or both (the import inside its patch, procedural continents outside, kept
//! a guard distance away). Brush strokes (R3.3) apply on top.
//!
//! Coverage is sampled **at node positions** from a continuous source (bilinear raster, noise
//! field, distance to a stroke), never rasterised at one fixed resolution, so coasts agree
//! between resolutions (the s03 lesson).

use crate::Params;
use crate::grid::Grid;
use crate::util::{par_map, solve_threshold};
use glam::DVec3;
use std::f64::consts::PI;
use std::sync::Arc;

/// Where an imported raster sits on the sphere.
#[derive(Clone, Debug, PartialEq)]
pub enum Placement {
    /// A 2:1 image covering the whole sphere, plate carrée (row 0 = north, column 0 = 180° W).
    WorldEquirect,
    /// A regional map in an azimuthal-equidistant projection about its centre.
    /// `rotation_deg`: degrees clockwise from north of the image's "up".
    Patch { center_deg: [f64; 2], km_per_px: f64, rotation_deg: f64 },
}

/// Which channel of a PNG is the land signal.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum Channel {
    /// Alpha if it varies across the image, otherwise brightness.
    #[default]
    Auto,
    Brightness,
    Alpha,
}

/// How the source file was turned into raw values (PNG only; stored so the world file
/// rebuilds the same raster).
#[derive(Clone, Debug, PartialEq, Default)]
pub struct DecodeOpts {
    pub channel: Channel,
    /// Source-pixel rectangle `[x, y, w, h]` kept from the image; `None` = all of it.
    pub crop: Option<[u32; 4]>,
}

/// The land test applied to the raw values (s06, Matt: greyscale coastline maps).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaskAdjust {
    /// Raw value that becomes the coast (0..1). Above it is land.
    pub threshold: f32,
    /// Swap land and sea (dark land on a light sea).
    pub invert: bool,
}

impl Default for MaskAdjust {
    fn default() -> Self {
        Self { threshold: 0.5, invert: false }
    }
}

impl MaskAdjust {
    /// Coverage (0.5 = coast) from a raw, already-interpolated value. Shifting after the
    /// bilinear step puts the coast exactly on the threshold contour. The UI preview calls
    /// this too, so it can't drift from what Generate does.
    pub fn apply(&self, raw: f32) -> f32 {
        let v = if self.invert { 1.0 - raw } else { raw };
        if self.threshold == 0.5 {
            // Exactly the pre-s06 value, so older imported worlds stay byte-identical.
            v
        } else {
            (v - self.threshold) + 0.5
        }
    }
}

/// How much erosion may move a heightmap layer's heights (s07, `docs/design/heightmap-import.md`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Default)]
pub enum Fidelity {
    /// The import is the canonical surface (drained); refine adds sub-cell detail only.
    Exact,
    /// Erosion runs on uplift inverted from the import, so real valleys form close to it.
    #[default]
    Guided,
}

/// A heightmap layer's height mapping (s07). The coast grey is the layer's
/// `MaskAdjust::threshold` (g_sea); greys above it are land:
/// `h = max_m · ((g − g_sea)/(1 − g_sea))^gamma`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Elevation {
    /// Height of white, m. The panel sets the scale (Matt, s06: white has no fixed meaning).
    pub max_m: f32,
    /// Curve exponent on land (> 1 flattens the lowlands).
    pub gamma: f32,
    /// Use the map's sub-sea greys as depth (decided at import from their spread); otherwise
    /// the sea floor is R3's bathymetry.
    pub sea_floor: bool,
    /// Depth of black, m, when `sea_floor`: `depth = sea_depth_m · (g_sea − g)/g_sea`.
    pub sea_depth_m: f32,
    /// Gaussian σ in **source** pixels, applied to the heights only (Azgaar cell steps,
    /// 8-bit terraces). 0 = off. The coast comes from the unsmoothed raster.
    pub smooth_px: f32,
    pub fidelity: Fidelity,
}

impl Default for Elevation {
    fn default() -> Self {
        Self { max_m: 4000.0, gamma: 1.0, sea_floor: false, sea_depth_m: 4500.0, smooth_px: 0.0, fidelity: Fidelity::Guided }
    }
}

impl Elevation {
    /// Height (m) from a raw grey `g` (already inverted if the layer inverts) and the coast
    /// grey `g_sea`. `None` at sea without `sea_floor` (R3 bathymetry applies there).
    pub fn height_m(&self, g: f32, g_sea: f32) -> Option<f32> {
        if g > g_sea {
            let t = ((g - g_sea) / (1.0 - g_sea).max(1e-6)).clamp(0.0, 1.0);
            Some((self.max_m * t.powf(self.gamma.max(0.05))).max(1.0))
        } else if self.sea_floor {
            let t = ((g_sea - g) / g_sea.max(1e-6)).clamp(0.0, 1.0);
            Some((-self.sea_depth_m * t).min(-1.0))
        } else {
            None
        }
    }
}

/// A decoded raster: raw values 0..1 per pixel, turned into coverage (0 sea .. 1 land,
/// 0.5 = coast) by `adjust`.
#[derive(Clone, PartialEq)]
pub struct RasterMask {
    /// Raw values. Behind an `Arc` so editing the placement or the land test (a cheap clone
    /// of the mask) never copies the raster.
    pub data: Arc<Vec<f32>>,
    pub w: u32,
    pub h: u32,
    pub place: Placement,
    pub adjust: MaskAdjust,
    pub decode: DecodeOpts,
    /// The channel the decoder actually used (`Auto` resolved), for the UI.
    pub channel_used: Channel,
    /// Switched off in the edit list: kept (non-destructive) but not applied (s06).
    pub enabled: bool,
    /// The file's name, for the edit list.
    pub name: String,
    /// The original file (PNG/SVG bytes) and its kind, kept so the world file can store it
    /// and the world regenerates anywhere. Opaque to this crate.
    pub source: Arc<Vec<u8>>,
    pub source_kind: String,
    /// Source pixels per raster pixel (box downsample at import; s07). 1 for coastline
    /// layers. `place.km_per_px` is per **raster** pixel.
    pub step: u32,
    /// `Some` for a heightmap layer (s07).
    pub elevation: Option<Elevation>,
    /// The smoothed raw greys as a pyramid, for heights ([`RasterMask::build_heights`]).
    pub heights: Option<Arc<Pyramid>>,
}

impl std::fmt::Debug for RasterMask {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "RasterMask {{ {}x{}, {:?} }}", self.w, self.h, self.place)
    }
}

/// What a brush stroke does (R3.3).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StrokeOp {
    Land,
    Sea,
    Active,
    Passive,
}

/// One brush stroke: a polyline on the sphere with a radius. Resolution-independent (a
/// distance test), small enough to store as JSON.
#[derive(Clone, Debug, PartialEq)]
pub struct Stroke {
    pub path: Vec<DVec3>,
    pub radius_km: f64,
    pub op: StrokeOp,
    /// The seed it was painted on: it applies only there (s06, Matt). `None` = every seed
    /// (worlds saved before s06).
    pub seed: Option<u64>,
    /// Switched off in the edit list: kept (non-destructive) but not applied.
    pub enabled: bool,
}

impl Stroke {
    /// A stroke that applies on every seed (tests, old files).
    pub fn new(path: Vec<DVec3>, radius_km: f64, op: StrokeOp) -> Self {
        Self { path, radius_km, op, seed: None, enabled: true }
    }

    /// Does this stroke apply to a world with `seed`?
    pub fn applies(&self, seed: u64) -> bool {
        self.enabled && self.seed.is_none_or(|s| s == seed)
    }

    /// Great-circle distance (radians) from `d` to the stroke's polyline.
    pub fn distance_rad(&self, d: DVec3) -> f64 {
        let mut best = PI;
        for (i, &a) in self.path.iter().enumerate() {
            best = best.min(a.dot(d).clamp(-1.0, 1.0).acos());
            if let Some(&b) = self.path.get(i + 1) {
                let axis = a.cross(b);
                let l = axis.length();
                if l < 1e-12 {
                    continue;
                }
                let axis = axis / l;
                let s = d.dot(axis);
                let q = d - axis * s;
                if q.length() < 1e-12 {
                    continue;
                }
                let q = q.normalize();
                // q lies on the arc between a and b?
                if a.cross(q).dot(axis) >= 0.0 && q.cross(b).dot(axis) >= 0.0 {
                    best = best.min(s.abs().clamp(-1.0, 1.0).asin());
                }
            }
        }
        best
    }

    /// Soft inside weight 0..1 at `d` (1 inside, fading over the outer 20% of the radius).
    pub fn weight(&self, d: DVec3, radius_m: f64) -> f64 {
        let dist = self.distance_rad(d) * radius_m / 1000.0;
        let r = self.radius_km.max(1.0);
        (1.0 - (dist - 0.8 * r) / (0.2 * r)).clamp(0.0, 1.0)
    }
}

/// Where the coasts come from (R3). The default is today's behaviour: procedural crust,
/// coast emerging from the heights by the sea-level bisection.
#[derive(Clone, Debug, PartialEq)]
pub struct Conditioning {
    /// Coasts are drawn first and the heights conform to them (no bisection).
    pub strict: bool,
    /// Imported rasters (strict only), applied in order: a later layer wins inside its own
    /// patch, like image layers (s06: several coastlines on one world).
    pub masks: Vec<Arc<RasterMask>>,
    /// Where no import covers: procedural continents (kept `guard_km` from the imports' land),
    /// or open ocean: a blank world (s06).
    pub outside_procedural: bool,
    pub guard_km: f64,
    /// Tectonic islands (arcs, hotspots) may rise above sea level this far from drawn coasts.
    pub islands: bool,
    pub islands_guard_km: f64,
    pub strokes: Vec<Stroke>,
}

impl Default for Conditioning {
    fn default() -> Self {
        Self {
            strict: false,
            masks: Vec::new(),
            outside_procedural: true,
            guard_km: 800.0,
            islands: true,
            islands_guard_km: 100.0,
            strokes: Vec::new(),
        }
    }
}

pub fn dir_of(lat_deg: f64, lon_deg: f64) -> DVec3 {
    let (sl, cl) = lat_deg.to_radians().sin_cos();
    let (so, co) = lon_deg.to_radians().sin_cos();
    DVec3::new(cl * co, cl * so, sl)
}

impl RasterMask {
    /// Fractional pixel coordinates of direction `d`, or `None` outside a patch.
    pub fn pixel(&self, d: DVec3, radius_m: f64) -> Option<(f64, f64)> {
        let (w, h) = (self.w as f64, self.h as f64);
        match self.place {
            Placement::WorldEquirect => {
                let lat = d.z.clamp(-1.0, 1.0).asin();
                let lon = d.y.atan2(d.x);
                Some(((lon + PI) / (2.0 * PI) * w - 0.5, (PI / 2.0 - lat) / PI * h - 0.5))
            }
            Placement::Patch { center_deg, km_per_px, rotation_deg } => {
                let c = dir_of(center_deg[0], center_deg[1]);
                let (sp, cp) = center_deg[0].to_radians().sin_cos();
                let (sl, cl) = center_deg[1].to_radians().sin_cos();
                let north = DVec3::new(-sp * cl, -sp * sl, cp);
                let east = DVec3::new(-sl, cl, 0.0);
                let theta = c.dot(d).clamp(-1.0, 1.0).acos();
                let t = d - c * c.dot(d);
                let az = if t.length() < 1e-12 { 0.0 } else { t.dot(east).atan2(t.dot(north)) };
                let dist_km = theta * radius_m / 1000.0;
                let a = az - rotation_deg.to_radians();
                let (x_km, y_km) = (dist_km * a.sin(), dist_km * a.cos());
                let px = w * 0.5 + x_km / km_per_px.max(1e-9) - 0.5;
                let py = h * 0.5 - y_km / km_per_px.max(1e-9) - 0.5;
                let inside = (-0.5..=w - 0.5).contains(&px) && (-0.5..=h - 0.5).contains(&py) && theta < PI * 0.95;
                inside.then_some((px, py))
            }
        }
    }

    /// Bilinear coverage at `d`, or `None` outside a patch.
    pub fn coverage(&self, d: DVec3, radius_m: f64) -> Option<f32> {
        let (px, py) = self.pixel(d, radius_m)?;
        let (w, h) = (self.w as i64, self.h as i64);
        let wrap = matches!(self.place, Placement::WorldEquirect);
        let at = |x: i64, y: i64| -> f32 {
            let x = if wrap { x.rem_euclid(w) } else { x.clamp(0, w - 1) };
            let y = y.clamp(0, h - 1);
            self.data[(y * w + x) as usize]
        };
        let (x0, y0) = (px.floor(), py.floor());
        let (fx, fy) = ((px - x0) as f32, (py - y0) as f32);
        let (x0, y0) = (x0 as i64, y0 as i64);
        let top = at(x0, y0) * (1.0 - fx) + at(x0 + 1, y0) * fx;
        let bot = at(x0, y0 + 1) * (1.0 - fx) + at(x0 + 1, y0 + 1) * fx;
        Some(self.adjust.apply(top * (1.0 - fy) + bot * fy))
    }

    /// (Re)build the height pyramid from the raw greys, smoothed by `elevation.smooth_px`
    /// (source pixels → raster pixels via `step`). No-op on a coastline layer.
    pub fn build_heights(&mut self) {
        let Some(e) = self.elevation else {
            self.heights = None;
            return;
        };
        let sigma = e.smooth_px as f64 / self.step.max(1) as f64;
        let wrap = matches!(self.place, Placement::WorldEquirect);
        let base = if sigma > 0.05 { gaussian_blur(&self.data, self.w as usize, self.h as usize, sigma, wrap) } else { (*self.data).clone() };
        self.heights = Some(Arc::new(Pyramid::build(&base, self.w, self.h)));
    }

    /// Imported height (m) at `d` for a sample `spacing_m` wide, or `None` outside the patch,
    /// on a coastline layer, or at sea without the map's sea floor.
    pub fn height_m(&self, d: DVec3, radius_m: f64, spacing_m: f64) -> Option<f32> {
        let (e, pyr) = (self.elevation?, self.heights.as_ref()?);
        let (px, py) = self.pixel(d, radius_m)?;
        let px_m = match self.place {
            Placement::WorldEquirect => 2.0 * PI * radius_m / self.w as f64,
            Placement::Patch { km_per_px, .. } => km_per_px * 1000.0,
        };
        let raw = pyr.sample(px, py, spacing_m / px_m.max(1e-9), matches!(self.place, Placement::WorldEquirect));
        let g = if self.adjust.invert { 1.0 - raw } else { raw };
        e.height_m(g, self.adjust.threshold)
    }

    /// The direction at fractional pixel `(px, py)`: the inverse of [`Self::pixel`].
    pub fn dir_of_pixel(&self, px: f64, py: f64, radius_m: f64) -> DVec3 {
        let (w, h) = (self.w as f64, self.h as f64);
        match self.place {
            Placement::WorldEquirect => {
                let lon = (px + 0.5) / w * (2.0 * PI) - PI;
                let lat = PI / 2.0 - (py + 0.5) / h * PI;
                DVec3::new(lat.cos() * lon.cos(), lat.cos() * lon.sin(), lat.sin())
            }
            Placement::Patch { center_deg, km_per_px, rotation_deg } => {
                let c = dir_of(center_deg[0], center_deg[1]);
                let (sp, cp) = center_deg[0].to_radians().sin_cos();
                let (sl, cl) = center_deg[1].to_radians().sin_cos();
                let north = DVec3::new(-sp * cl, -sp * sl, cp);
                let east = DVec3::new(-sl, cl, 0.0);
                let x_km = (px + 0.5 - w * 0.5) * km_per_px;
                let y_km = (h * 0.5 - py - 0.5) * km_per_px;
                let theta = x_km.hypot(y_km) * 1000.0 / radius_m;
                let az = x_km.atan2(y_km) + rotation_deg.to_radians();
                c * theta.cos() + (north * az.cos() + east * az.sin()) * theta.sin()
            }
        }
    }
}

/// Separable Gaussian blur (σ in pixels, kernel ±3σ), edges clamped (x wraps if `wrap`).
/// Sequential, so deterministic.
pub fn gaussian_blur(data: &[f32], w: usize, h: usize, sigma: f64, wrap: bool) -> Vec<f32> {
    let r = (3.0 * sigma).ceil().max(1.0) as i64;
    let k: Vec<f64> = (-r..=r).map(|i| (-(i * i) as f64 / (2.0 * sigma * sigma)).exp()).collect();
    let ks: f64 = k.iter().sum();
    let (wi, hi) = (w as i64, h as i64);
    let mut tmp = vec![0.0f32; w * h];
    for y in 0..h {
        for x in 0..wi {
            let mut acc = 0.0f64;
            for (j, kv) in k.iter().enumerate() {
                let xx = x + j as i64 - r;
                let xx = if wrap { xx.rem_euclid(wi) } else { xx.clamp(0, wi - 1) };
                acc += kv * data[y * w + xx as usize] as f64;
            }
            tmp[y * w + x as usize] = (acc / ks) as f32;
        }
    }
    let mut out = vec![0.0f32; w * h];
    for y in 0..hi {
        for x in 0..w {
            let mut acc = 0.0f64;
            for (j, kv) in k.iter().enumerate() {
                let yy = (y + j as i64 - r).clamp(0, hi - 1);
                acc += kv * tmp[yy as usize * w + x] as f64;
            }
            out[y as usize * w + x] = (acc / ks) as f32;
        }
    }
    out
}

/// A mip pyramid of a raster: level 0 is the raster, each next level is 2× coarser by box
/// averaging (a partial block at the right or bottom edge averages what it holds). Heights
/// are sampled at a node's footprint, not at a point, so a 0.6 km/px map doesn't alias
/// differently at Normal and Ultra (s07 review: that would break the s04 coherence guard).
/// Built sequentially, so it's deterministic.
#[derive(Clone, Debug, PartialEq)]
pub struct Pyramid {
    /// (values, width, height) per level, finest first.
    pub levels: Vec<(Vec<f32>, u32, u32)>,
}

impl Pyramid {
    pub fn build(data: &[f32], w: u32, h: u32) -> Self {
        let mut levels = vec![(data.to_vec(), w, h)];
        while let Some((d, w, h)) = levels.last().filter(|(_, w, h)| (*w).max(*h) > 2) {
            let (w, h) = (*w as usize, *h as usize);
            let (nw, nh) = (w.div_ceil(2), h.div_ceil(2));
            let mut out = Vec::with_capacity(nw * nh);
            for y in 0..nh {
                for x in 0..nw {
                    let (mut sum, mut n) = (0.0f64, 0u32);
                    for yy in 2 * y..(2 * y + 2).min(h) {
                        for xx in 2 * x..(2 * x + 2).min(w) {
                            sum += d[yy * w + xx] as f64;
                            n += 1;
                        }
                    }
                    out.push((sum / n as f64) as f32);
                }
            }
            levels.push((out, nw as u32, nh as u32));
        }
        Self { levels }
    }

    /// Bilinear value at one level, at fractional **level-0** pixel `(px, py)` (pixel centres
    /// at integers). `wrap`: the x axis wraps (a whole-world map).
    fn bilinear(&self, level: usize, px: f64, py: f64, wrap: bool) -> f32 {
        let (d, w, h) = &self.levels[level];
        let (w, h) = (*w as i64, *h as i64);
        let f = (1u64 << level) as f64;
        let (x, y) = ((px + 0.5) / f - 0.5, (py + 0.5) / f - 0.5);
        let at = |x: i64, y: i64| -> f32 {
            let x = if wrap { x.rem_euclid(w) } else { x.clamp(0, w - 1) };
            d[(y.clamp(0, h - 1) * w + x) as usize]
        };
        let (x0, y0) = (x.floor(), y.floor());
        let (fx, fy) = ((x - x0) as f32, (y - y0) as f32);
        let (x0, y0) = (x0 as i64, y0 as i64);
        let top = at(x0, y0) * (1.0 - fx) + at(x0 + 1, y0) * fx;
        let bot = at(x0, y0 + 1) * (1.0 - fx) + at(x0 + 1, y0 + 1) * fx;
        top * (1.0 - fy) + bot * fy
    }

    /// Trilinear value at fractional level-0 pixel `(px, py)` for a sample `footprint_px`
    /// level-0 pixels wide: bilinear on the two levels around log2(footprint), blended.
    pub fn sample(&self, px: f64, py: f64, footprint_px: f64, wrap: bool) -> f32 {
        let top = self.levels.len() - 1;
        let lf = footprint_px.max(1.0).log2().min(top as f64);
        let l0 = lf.floor() as usize;
        let t = (lf - l0 as f64) as f32;
        let a = self.bilinear(l0, px, py, wrap);
        if t <= 0.0 || l0 >= top {
            return a;
        }
        a * (1.0 - t) + self.bilinear(l0 + 1, px, py, wrap) * t
    }
}

/// Land coverage per node, 0..1 (land = > 0.5), for strict mode. `value` is the procedural
/// continentalness (`crust::continentalness`), thresholded straight to the land fraction.
pub fn land_cover(grid: &Grid, p: &Params, value: &[f32]) -> Vec<f32> {
    let c = &p.conditioning;
    let r = p.radius_m;
    let thr = solve_threshold(value, &grid.weight, p.land_fraction.clamp(0.0, 1.0) as f64);
    let procedural = |a: usize| if value[a] > thr { 1.0f32 } else { 0.0 };
    // Imports, in order: a later layer replaces the earlier ones inside its own patch.
    let mut cover = vec![0.0f32; grid.n];
    let mut covered = vec![false; grid.n];
    for m in c.masks.iter().filter(|m| m.enabled) {
        let raster: Vec<Option<f32>> = par_map(grid.n, |a| m.coverage(grid.dir[a], r));
        for (a, v) in raster.into_iter().enumerate() {
            if let Some(v) = v {
                cover[a] = v.clamp(0.0, 1.0);
                covered[a] = true;
            }
        }
    }
    // Elsewhere: procedural continents, kept `guard_km` from the imports' land, or open sea.
    if c.outside_procedural {
        let import_land: Vec<u32> = (0..grid.n).filter(|&a| covered[a] && cover[a] > 0.5).map(|a| a as u32).collect();
        let guard = grid.nearest_source(&import_land, ((c.guard_km * 1000.0 / r) as f32).max(1e-6));
        for a in 0..grid.n {
            if !covered[a] && guard.src[a] == u32::MAX {
                cover[a] = procedural(a);
            }
        }
    }
    for st in c.strokes.iter().filter(|st| st.applies(p.seed)) {
        if !matches!(st.op, StrokeOp::Land | StrokeOp::Sea) || st.path.is_empty() {
            continue;
        }
        // Cheap reject: a cap around the stroke's points.
        let centre = st.path.iter().copied().sum::<DVec3>().normalize_or_zero();
        let extent = st.path.iter().map(|q| q.dot(centre).clamp(-1.0, 1.0).acos()).fold(0.0, f64::max);
        let cap = (extent + st.radius_km * 1000.0 / r).min(PI).cos();
        let w: Vec<f32> = par_map(grid.n, |a| {
            let d = grid.dir[a];
            if centre != DVec3::ZERO && d.dot(centre) < cap { 0.0 } else { st.weight(d, r) as f32 }
        });
        for a in 0..grid.n {
            cover[a] = match st.op {
                StrokeOp::Land => cover[a].max(w[a]),
                _ => cover[a].min(1.0 - w[a]),
            };
        }
    }
    cover
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_mask(w: u32, h: u32, place: Placement) -> RasterMask {
        RasterMask {
            data: Arc::new(vec![0.0; (w * h) as usize]),
            w,
            h,
            place,
            adjust: MaskAdjust::default(),
            decode: DecodeOpts::default(),
            channel_used: Channel::Brightness,
            enabled: true,
            name: String::new(),
            source: Arc::new(vec![]),
            source_kind: String::new(),
            step: 1,
            elevation: None,
            heights: None,
        }
    }

    /// `pixel` and `dir_of_pixel` invert each other: the gate for dragging a patch.
    #[test]
    fn pixel_round_trips() {
        let r = 6.371e6;
        let places = [
            Placement::WorldEquirect,
            Placement::Patch { center_deg: [10.0, 20.0], km_per_px: 10.0, rotation_deg: 0.0 },
            Placement::Patch { center_deg: [-35.0, 170.0], km_per_px: 25.0, rotation_deg: 37.0 },
            Placement::Patch { center_deg: [88.0, -60.0], km_per_px: 15.0, rotation_deg: -120.0 },
        ];
        for place in places {
            let m = test_mask(200, 100, place.clone());
            for (px, py) in [(0.0, 0.0), (199.0, 99.0), (57.3, 12.8), (100.0, 50.0), (3.0, 97.0)] {
                let d = m.dir_of_pixel(px, py, r);
                assert!((d.length() - 1.0).abs() < 1e-9);
                let (qx, qy) = m.pixel(d, r).unwrap();
                assert!((qx - px).abs() < 1e-6 && (qy - py).abs() < 1e-6, "{place:?} ({px},{py}) -> ({qx},{qy})");
            }
        }
    }

    #[test]
    fn adjust_threshold_and_invert() {
        let a = MaskAdjust::default();
        assert_eq!(a.apply(0.3), 0.3);
        let a = MaskAdjust { threshold: 0.2, invert: false };
        assert!((a.apply(0.2) - 0.5).abs() < 1e-6 && a.apply(0.3) > 0.5 && a.apply(0.1) < 0.5);
        let a = MaskAdjust { threshold: 0.5, invert: true };
        assert!(a.apply(0.9) < 0.5 && a.apply(0.1) > 0.5);
    }

    #[test]
    fn equirect_maps_corners() {
        let m = test_mask(4, 2, Placement::WorldEquirect);
        let (x, y) = m.pixel(dir_of(0.0, 0.0), 6.371e6).unwrap();
        assert!((x - 1.5).abs() < 1e-9 && (y - 0.5).abs() < 1e-9, "{x} {y}");
    }

    #[test]
    fn patch_centre_and_distance() {
        let m = test_mask(100, 100, Placement::Patch { center_deg: [10.0, 20.0], km_per_px: 10.0, rotation_deg: 0.0 });
        let r = 6.371e6;
        let (x, y) = m.pixel(dir_of(10.0, 20.0), r).unwrap();
        assert!((x - 49.5).abs() < 1e-6 && (y - 49.5).abs() < 1e-6);
        // 1° north ≈ 111.2 km ≈ 11.1 px up.
        let (x, y) = m.pixel(dir_of(11.0, 20.0), r).unwrap();
        assert!((x - 49.5).abs() < 1e-6 && (y - (49.5 - 11.12)).abs() < 0.02, "{x} {y}");
        // Far outside the patch.
        assert!(m.pixel(dir_of(-60.0, 20.0), r).is_none());
    }

    #[test]
    fn stroke_distance() {
        let s = Stroke::new(vec![dir_of(0.0, 0.0), dir_of(0.0, 10.0)], 100.0, StrokeOp::Land);
        let d = s.distance_rad(dir_of(1.0, 5.0)).to_degrees();
        assert!((d - 1.0).abs() < 1e-6, "{d}");
        assert!(s.weight(dir_of(0.0, 5.0), 6.371e6) > 0.99);
        assert!(s.weight(dir_of(5.0, 5.0), 6.371e6) == 0.0);
    }

    /// Heights: the curve on land, the map's sea floor only when asked, and pyramid sampling
    /// through the placement.
    #[test]
    fn heightmap_layer_heights() {
        let e = Elevation { max_m: 4000.0, ..Default::default() };
        assert_eq!(e.height_m(1.0, 0.2), Some(4000.0));
        assert!((e.height_m(0.6, 0.2).unwrap() - 2000.0).abs() < 1e-3);
        assert_eq!(e.height_m(0.1, 0.2), None);
        let es = Elevation { sea_floor: true, sea_depth_m: 4000.0, ..e };
        assert!((es.height_m(0.1, 0.2).unwrap() + 2000.0).abs() < 1e-3);
        let mut m = test_mask(8, 8, Placement::Patch { center_deg: [0.0, 0.0], km_per_px: 50.0, rotation_deg: 0.0 });
        m.data = Arc::new(vec![0.6; 64]);
        m.adjust.threshold = 0.2;
        m.elevation = Some(e);
        m.build_heights();
        let h = m.height_m(dir_of(0.0, 0.0), 6.371e6, 20_000.0).unwrap();
        assert!((h - 2000.0).abs() < 1e-2, "{h}");
        assert!(m.height_m(dir_of(40.0, 0.0), 6.371e6, 20_000.0).is_none());
        // Blur keeps a constant field constant.
        let b = gaussian_blur(&[3.0; 30], 6, 5, 1.5, true);
        assert!(b.iter().all(|v| (v - 3.0).abs() < 1e-5));
    }

    /// The pyramid: level sizes halve (rounding up), each level's mean matches level 0 on an
    /// even-sized raster, a 1-px footprint is plain bilinear, and a wide footprint sees the
    /// block mean rather than one aliased pixel.
    #[test]
    fn pyramid_levels_and_sampling() {
        let (w, h) = (16u32, 8u32);
        let data: Vec<f32> = (0..w * h).map(|k| if (k % w + k / w) % 2 == 0 { 1.0 } else { 0.0 }).collect();
        let p = Pyramid::build(&data, w, h);
        let sizes: Vec<(u32, u32)> = p.levels.iter().map(|l| (l.1, l.2)).collect();
        assert_eq!(sizes, vec![(16, 8), (8, 4), (4, 2), (2, 1)]);
        for (d, _, _) in &p.levels {
            let mean = d.iter().sum::<f32>() / d.len() as f32;
            assert!((mean - 0.5).abs() < 1e-6);
        }
        assert_eq!(p.sample(3.0, 2.0, 1.0, false), data[(2 * w + 3) as usize]);
        // A checkerboard seen through a 4-px footprint is grey, not 0 or 1.
        assert!((p.sample(5.3, 3.7, 4.0, false) - 0.5).abs() < 1e-5);
        // Odd sizes round up and wrap works.
        let p = Pyramid::build(&[1.0; 15], 5, 3);
        assert_eq!((p.levels[1].1, p.levels[1].2), (3, 2));
        assert!((p.sample(-0.4, 1.0, 1.0, true) - 1.0).abs() < 1e-6);
    }

    /// Two overlapping patches: the later one wins inside its own patch; outside both, a
    /// blank base is open sea; a disabled layer doesn't apply.
    #[test]
    fn layers_later_wins_and_blank_base() {
        let mut land = test_mask(10, 10, Placement::Patch { center_deg: [0.0, 0.0], km_per_px: 100.0, rotation_deg: 0.0 });
        land.data = Arc::new(vec![1.0; 100]);
        let mut sea = test_mask(4, 4, Placement::Patch { center_deg: [0.0, 0.0], km_per_px: 100.0, rotation_deg: 0.0 });
        sea.data = Arc::new(vec![0.0; 16]);
        let mut p = Params { res: 48, ..Default::default() };
        p.conditioning.strict = true;
        p.conditioning.outside_procedural = false;
        p.conditioning.masks = vec![Arc::new(land), Arc::new(sea.clone())];
        let g = Grid::new(48);
        let value = vec![0.0f32; g.n];
        let at = |cover: &[f32], lat: f64, lon: f64| cover[g.nearest_node(dir_of(lat, lon))];
        let cover = land_cover(&g, &p, &value);
        assert!(at(&cover, 0.0, 0.0) < 0.5, "the later (sea) layer wins at the centre");
        assert!(at(&cover, 0.0, 3.5) > 0.5, "the earlier land shows outside the small patch");
        assert!(at(&cover, 40.0, 120.0) < 0.5, "blank base: open sea elsewhere");
        let mut off = sea;
        off.enabled = false;
        p.conditioning.masks[1] = Arc::new(off);
        let cover = land_cover(&g, &p, &value);
        assert!(at(&cover, 0.0, 0.0) > 0.5, "a disabled layer doesn't apply");
    }
}
