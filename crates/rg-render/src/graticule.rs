//! The latitude/longitude grid: how far apart to draw the lines at the current zoom, and
//! where the degree labels go on screen. The lines themselves are drawn in `planet.wgsl`
//! (they get the spacing through `Globals::lod.w`); the labels are drawn by the UI.
//!
//! Screen positions here are viewport fractions: x and y in 0..1, y pointing down.

use crate::earth_camera::dir;
use crate::{EarthCamera, MapView, ViewMode};

/// Allowed line spacings in degrees, coarse to fine. All divide 360, so the shader's
/// longitude wrap-around lands exactly on a line.
pub const STEPS_DEG: [f64; 6] = [30.0, 15.0, 10.0, 5.0, 2.0, 1.0];

/// Closest we let two lines get, in screen pixels.
const MIN_LINE_GAP_PX: f64 = 80.0;

/// Pick the finest spacing whose lines stay at least `MIN_LINE_GAP_PX` apart.
/// If even the coarsest is tighter than that (tiny viewport), use the coarsest.
pub fn pick_step_deg(px_per_deg: f64) -> f64 {
    STEPS_DEG
        .iter()
        .rev()
        .copied()
        .find(|s| s * px_per_deg >= MIN_LINE_GAP_PX)
        .unwrap_or(STEPS_DEG[0])
}

/// Spacing for the current view. `size_px` is the viewport in pixels.
pub fn step_deg(mode: ViewMode, cam: &EarthCamera, map: &MapView, size_px: [u32; 2]) -> f64 {
    match mode {
        ViewMode::Globe => {
            // Metres per degree of arc, over ground metres per pixel at the view centre.
            let mpp = cam.metres_per_px(size_px[1] as f64);
            pick_step_deg(cam.radius_m.to_radians() / mpp)
        }
        ViewMode::Map => pick_step_deg(size_px[0] as f64 / map.width_rad.to_degrees()),
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LabelKind {
    /// A line of latitude. On the map it sits on the left edge (`fx` is 0).
    Parallel,
    /// A line of longitude. On the map it sits on the bottom edge (`fy` is 1).
    Meridian,
}

#[derive(Clone, Debug)]
pub struct GraticuleLabel {
    pub kind: LabelKind,
    pub text: String,
    pub fx: f64,
    pub fy: f64,
}

/// "30°N", "45°S", "0°".
pub fn format_lat(deg: i32) -> String {
    match deg {
        0 => "0°".to_string(),
        d if d > 0 => format!("{d}°N"),
        d => format!("{}°S", -d),
    }
}

/// "45°E", "45°W", "0°", "180°".
pub fn format_lon(deg: i32) -> String {
    // Fold into (-180, 180].
    let d = (deg + 180).rem_euclid(360) - 180;
    let d = if d == -180 { 180 } else { d };
    match d {
        0 | 180 => format!("{d}°"),
        d if d > 0 => format!("{d}°E"),
        d => format!("{}°W", -d),
    }
}

/// Project a point on the sea-level sphere to viewport fractions. `None` when it is on the
/// far side of the horizon, behind the camera or outside the viewport.
pub fn project_globe(cam: &EarthCamera, aspect: f64, lat: f64, lon: f64) -> Option<[f64; 2]> {
    let up = dir(lat, lon);
    let eye = cam.eye();
    // Visible side of a sphere: the point lies on the eye's side of its horizon circle.
    if up.dot(eye) <= cam.radius_m {
        return None;
    }
    let clip = cam.view_proj_f64(aspect, 1000.0) * (up * cam.radius_m - eye).extend(1.0);
    if clip.w <= 0.0 {
        return None;
    }
    let (x, y) = (clip.x / clip.w, clip.y / clip.w);
    if x.abs() > 1.0 || y.abs() > 1.0 {
        return None;
    }
    Some([(x + 1.0) * 0.5, (1.0 - y) * 0.5])
}

/// Labels for the globe: each parallel where it crosses the meridian through the view
/// centre, and each meridian where it crosses the parallel through the view centre.
pub fn globe_labels(cam: &EarthCamera, aspect: f64, step: f64) -> Vec<GraticuleLabel> {
    let mut out = Vec::new();
    let n = (90.0 / step).floor() as i32;
    for k in -n..=n {
        let deg = (k as f64 * step).round() as i32;
        if deg.abs() >= 90 {
            continue;
        }
        if let Some([fx, fy]) = project_globe(cam, aspect, (deg as f64).to_radians(), cam.lon) {
            out.push(GraticuleLabel { kind: LabelKind::Parallel, text: format_lat(deg), fx, fy });
        }
    }
    // Meridians bunch up at the poles, so their labels would pile on top of each other.
    if cam.lat.abs() < 80f64.to_radians() {
        let count = (360.0 / step).round() as i32;
        for k in 0..count {
            let deg = (k as f64 * step).round() as i32 - 180;
            if let Some([fx, fy]) = project_globe(cam, aspect, cam.lat, (deg as f64).to_radians()) {
                out.push(GraticuleLabel { kind: LabelKind::Meridian, text: format_lon(deg), fx, fy });
            }
        }
    }
    out
}

/// Labels for the map: parallels along the left edge, meridians along the bottom edge.
pub fn map_labels(map: &MapView, aspect: f64, step: f64) -> Vec<GraticuleLabel> {
    let mut out = Vec::new();
    let height = map.width_rad / aspect;
    let step_rad = step.to_radians();
    // Parallels: lat = center_lat - (fy - 0.5) * height.
    let lo = ((map.center_lat - height * 0.5) / step_rad).floor() as i32;
    let hi = ((map.center_lat + height * 0.5) / step_rad).ceil() as i32;
    for k in lo..=hi {
        let deg = (k as f64 * step).round() as i32;
        if deg.abs() > 90 {
            continue;
        }
        let fy = 0.5 - ((deg as f64).to_radians() - map.center_lat) / height;
        if (0.0..=1.0).contains(&fy) {
            out.push(GraticuleLabel { kind: LabelKind::Parallel, text: format_lat(deg), fx: 0.0, fy });
        }
    }
    // Meridians: lon = center_lon + (fx - 0.5) * width. The map can show more than 360°,
    // so the same meridian may appear twice.
    let lo = ((map.center_lon - map.width_rad * 0.5) / step_rad).floor() as i32;
    let hi = ((map.center_lon + map.width_rad * 0.5) / step_rad).ceil() as i32;
    for k in lo..=hi {
        let deg = (k as f64 * step).round() as i32;
        let fx = 0.5 + ((deg as f64).to_radians() - map.center_lon) / map.width_rad;
        if (0.0..=1.0).contains(&fx) {
            out.push(GraticuleLabel { kind: LabelKind::Meridian, text: format_lon(deg), fx, fy: 1.0 });
        }
    }
    out
}

/// Labels for the current view at the spacing the shader uses.
pub fn labels(mode: ViewMode, cam: &EarthCamera, map: &MapView, size_px: [u32; 2]) -> Vec<GraticuleLabel> {
    let step = step_deg(mode, cam, map, size_px);
    let aspect = size_px[0] as f64 / (size_px[1] as f64).max(1.0);
    match mode {
        ViewMode::Globe => globe_labels(cam, aspect, step),
        ViewMode::Map => map_labels(map, aspect, step),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn step_gets_finer_as_you_zoom_in() {
        assert_eq!(pick_step_deg(0.5), 30.0); // too small even for the coarsest
        assert_eq!(pick_step_deg(3.0), 30.0); // 15 deg would be 45 px
        assert_eq!(pick_step_deg(6.0), 15.0);
        assert_eq!(pick_step_deg(100.0), 1.0);
        assert_eq!(pick_step_deg(1e6), 1.0);
    }

    #[test]
    fn formats() {
        assert_eq!(format_lat(30), "30°N");
        assert_eq!(format_lat(-45), "45°S");
        assert_eq!(format_lat(0), "0°");
        assert_eq!(format_lon(-45), "45°W");
        assert_eq!(format_lon(190), "170°W");
        assert_eq!(format_lon(180), "180°");
        assert_eq!(format_lon(-180), "180°");
        assert_eq!(format_lon(0), "0°");
    }

    #[test]
    fn globe_centre_projects_to_middle_and_far_side_is_hidden() {
        let cam = EarthCamera::new(6_371_000.0);
        let [fx, fy] = project_globe(&cam, 1.5, cam.lat, cam.lon).unwrap();
        assert!((fx - 0.5).abs() < 1e-6 && (fy - 0.5).abs() < 1e-6);
        // The antipode is behind the planet.
        assert!(project_globe(&cam, 1.5, -cam.lat, cam.lon + std::f64::consts::PI).is_none());
        // A point on the near side but far off to the side is outside a narrow view.
        let close = EarthCamera { distance_m: 20_000.0, ..cam };
        assert!(project_globe(&close, 1.5, cam.lat + 0.5, cam.lon).is_none());
    }

    #[test]
    fn globe_labels_include_equator_and_stay_on_screen() {
        let cam = EarthCamera::new(6_371_000.0);
        let labels = globe_labels(&cam, 1.5, 30.0);
        assert!(labels.iter().any(|l| l.kind == LabelKind::Parallel && l.text == "0°"));
        assert!(labels.iter().any(|l| l.kind == LabelKind::Meridian));
        assert!(labels.iter().all(|l| (0.0..=1.0).contains(&l.fx) && (0.0..=1.0).contains(&l.fy)));
    }

    #[test]
    fn map_labels_match_lonlat_at() {
        let map = MapView { center_lon: 0.3, center_lat: 0.2, width_rad: 2.0 };
        let aspect = 1.6;
        let labels = map_labels(&map, aspect, 15.0);
        for l in &labels {
            let (lon, lat) = map.lonlat_at(l.fx, l.fy, aspect);
            let deg = match l.kind {
                LabelKind::Parallel => lat.to_degrees(),
                LabelKind::Meridian => lon.to_degrees(),
            };
            assert!((deg / 15.0 - (deg / 15.0).round()).abs() < 1e-6, "{deg}");
        }
        assert!(labels.iter().any(|l| l.kind == LabelKind::Parallel));
        assert!(labels.iter().any(|l| l.kind == LabelKind::Meridian));
    }
}
