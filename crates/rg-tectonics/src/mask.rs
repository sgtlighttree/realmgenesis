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
    /// The original file (PNG/SVG bytes) and its kind, kept so the world file can store it
    /// and the world regenerates anywhere. Opaque to this crate.
    pub source: Arc<Vec<u8>>,
    pub source_kind: String,
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
}

impl Stroke {
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
    /// An imported raster (strict only).
    pub mask: Option<Arc<RasterMask>>,
    /// Outside an imported patch: procedural continents (kept `guard_km` away) or ocean.
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
            mask: None,
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

/// Land coverage per node, 0..1 (land = > 0.5), for strict mode. `value` is the procedural
/// continentalness (`crust::continentalness`), thresholded straight to the land fraction.
pub fn land_cover(grid: &Grid, p: &Params, value: &[f32]) -> Vec<f32> {
    let c = &p.conditioning;
    let r = p.radius_m;
    let thr = solve_threshold(value, &grid.weight, p.land_fraction.clamp(0.0, 1.0) as f64);
    let procedural = |a: usize| if value[a] > thr { 1.0f32 } else { 0.0 };
    let mut cover: Vec<f32> = match &c.mask {
        None => (0..grid.n).map(procedural).collect(),
        Some(m) => {
            let raster: Vec<Option<f32>> = par_map(grid.n, |a| m.coverage(grid.dir[a], r));
            // Procedural continents outside the patch keep `guard_km` from its land.
            let patch_land: Vec<u32> =
                (0..grid.n).filter(|&a| raster[a].is_some_and(|v| v > 0.5)).map(|a| a as u32).collect();
            let guard = grid.nearest_source(&patch_land, ((c.guard_km * 1000.0 / r) as f32).max(1e-6));
            (0..grid.n)
                .map(|a| match raster[a] {
                    Some(v) => v.clamp(0.0, 1.0),
                    None if c.outside_procedural && guard.src[a] == u32::MAX => procedural(a),
                    None => 0.0,
                })
                .collect()
        }
    };
    for st in &c.strokes {
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
            source: Arc::new(vec![]),
            source_kind: String::new(),
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
        let s = Stroke { path: vec![dir_of(0.0, 0.0), dir_of(0.0, 10.0)], radius_km: 100.0, op: StrokeOp::Land };
        let d = s.distance_rad(dir_of(1.0, 5.0)).to_degrees();
        assert!((d - 1.0).abs() < 1e-6, "{d}");
        assert!(s.weight(dir_of(0.0, 5.0), 6.371e6) > 0.99);
        assert!(s.weight(dir_of(5.0, 5.0), 6.371e6) == 0.0);
    }
}
