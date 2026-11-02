//! A Google Earth–style camera: it looks at a target point on the planet from a distance,
//! with a heading (compass direction the view faces) and a tilt (0 = straight down).
//!
//! Everything is f64 metres. Rendering is camera-relative: [`EarthCamera::view_proj`]
//! puts the eye at the origin, and tile offsets are computed as (centre − eye) in f64.

use crate::camera::{look_at_rh, perspective_infinite_reverse_rh};
use glam::{DMat4, DVec3, Mat4};
use std::f64::consts::{PI, TAU};

#[derive(Clone, Copy, Debug)]
pub struct EarthCamera {
    pub radius_m: f64,
    /// Target latitude/longitude (radians) on the sea-level sphere.
    pub lat: f64,
    pub lon: f64,
    /// Distance from the target to the eye, metres.
    pub distance_m: f64,
    /// 0 = looking north, π/2 = looking east.
    pub heading: f64,
    /// 0 = straight down; grows towards the horizon. Capped by [`Self::max_tilt`].
    pub tilt: f64,
    pub fov_y: f64,
    /// Displayed terrain height under the target (m above sea level, exaggeration
    /// applied). The app updates it every frame so the camera orbits a point on the
    /// ground, not on the sea-level sphere under a mountain.
    pub target_height_m: f64,
}

impl EarthCamera {
    pub fn new(radius_m: f64) -> Self {
        Self { radius_m, lat: 0.35, lon: 0.6, distance_m: 2.2 * radius_m, heading: 0.0, tilt: 0.0, fov_y: 40f64.to_radians(), target_height_m: 0.0 }
    }

    pub const MIN_DISTANCE_M: f64 = 25.0;

    pub fn max_distance(&self) -> f64 {
        5.0 * self.radius_m
    }

    /// Tilt allowed at the current distance: none from orbit, up to 78° near the ground.
    pub fn max_tilt(&self) -> f64 {
        let t = smoothstep(4.0e6, 3.0e4, self.distance_m);
        78f64.to_radians() * t
    }

    /// Unit up vector at the target, plus east and north.
    pub fn frame(&self) -> (DVec3, DVec3, DVec3) {
        let up = dir(self.lat, self.lon);
        let east = DVec3::new(-self.lon.sin(), self.lon.cos(), 0.0);
        let north = up.cross(east);
        (up, east, north)
    }

    /// Horizontal direction the camera faces at the target.
    fn forward(&self) -> DVec3 {
        let (_, east, north) = self.frame();
        north * self.heading.cos() + east * self.heading.sin()
    }

    pub fn target(&self) -> DVec3 {
        dir(self.lat, self.lon) * (self.radius_m + self.target_height_m)
    }

    pub fn eye(&self) -> DVec3 {
        let (up, _, _) = self.frame();
        let f = self.forward();
        self.target() + (up * self.tilt.cos() - f * self.tilt.sin()) * self.distance_m
    }

    /// Altitude of the eye above the sea-level sphere, metres.
    pub fn altitude_m(&self) -> f64 {
        self.eye().length() - self.radius_m
    }

    /// View-projection with the eye at the origin (rotation only) and an infinite,
    /// reversed-Z projection. `clearance_m` is the eye's height above the terrain; the near
    /// plane is a tenth of it.
    pub fn view_proj(&self, aspect: f64, clearance_m: f64) -> Mat4 {
        self.view_proj_f64(aspect, clearance_m).as_mat4()
    }

    pub fn view_proj_f64(&self, aspect: f64, clearance_m: f64) -> DMat4 {
        let (up, _, _) = self.frame();
        let f = self.forward();
        let eye = self.eye();
        let screen_up = f * self.tilt.cos() + up * self.tilt.sin();
        let view = look_at_rh(DVec3::ZERO, self.target() - eye, screen_up);
        let near = (clearance_m * 0.1).clamp(0.05, 1.0e5);
        perspective_infinite_reverse_rh(self.fov_y, aspect, near) * view
    }

    /// Ground metres per screen pixel at the target.
    pub fn metres_per_px(&self, viewport_h: f64) -> f64 {
        2.0 * self.distance_m * (0.5 * self.fov_y).tan() / viewport_h.max(1.0)
    }

    /// Drag the ground: the point under the cursor follows it.
    pub fn pan(&mut self, dx: f64, dy: f64, viewport_h: f64) {
        let mpp = self.metres_per_px(viewport_h).min(self.radius_m * 2.5 / viewport_h.max(1.0));
        let (up, east, north) = self.frame();
        let f = self.forward();
        let right = f.cross(up);
        let delta = f * (dy * mpp) - right * (dx * mpp);
        self.lat = (self.lat + delta.dot(north) / self.radius_m).clamp(-1.55, 1.55);
        let coslat = self.lat.cos().max(0.02);
        self.lon = wrap(self.lon + delta.dot(east) / (self.radius_m * coslat));
    }

    /// Zoom by a factor on distance.
    pub fn zoom(&mut self, factor: f64) {
        self.distance_m = (self.distance_m * factor).clamp(Self::MIN_DISTANCE_M, self.max_distance());
        self.tilt = self.tilt.min(self.max_tilt());
    }

    /// Rotate heading and tilt (radians).
    pub fn orbit(&mut self, d_heading: f64, d_tilt: f64) {
        self.heading = (self.heading + d_heading).rem_euclid(TAU);
        self.tilt = (self.tilt + d_tilt).clamp(0.0, self.max_tilt());
    }

    /// Ray-pick the sea-level sphere under a viewport point (fractions 0..1, y down).
    /// Returns (lat, lon) radians.
    pub fn pick(&self, fx: f64, fy: f64, aspect: f64) -> Option<(f64, f64)> {
        let (up, _, _) = self.frame();
        let f = self.forward();
        let eye = self.eye();
        let fwd = (self.target() - eye).normalize();
        let screen_up = f * self.tilt.cos() + up * self.tilt.sin();
        let s = fwd.cross(screen_up).normalize();
        let u = s.cross(fwd);
        let th = (0.5 * self.fov_y).tan();
        let d = (fwd + s * ((fx * 2.0 - 1.0) * th * aspect) + u * ((1.0 - fy * 2.0) * th)).normalize();
        let r = self.radius_m;
        let b = eye.dot(d);
        let c = eye.length_squared() - r * r;
        let disc = b * b - c;
        if disc < 0.0 {
            return None;
        }
        let t = -b - disc.sqrt();
        if t < 0.0 {
            return None;
        }
        let p = (eye + d * t).normalize();
        Some((p.z.clamp(-1.0, 1.0).asin(), p.y.atan2(p.x)))
    }
}

pub fn dir(lat: f64, lon: f64) -> DVec3 {
    let (sl, cl) = lat.sin_cos();
    let (so, co) = lon.sin_cos();
    DVec3::new(cl * co, cl * so, sl)
}

fn wrap(lon: f64) -> f64 {
    (lon + PI).rem_euclid(TAU) - PI
}

fn smoothstep(a: f64, b: f64, x: f64) -> f64 {
    let t = ((x - a) / (b - a)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn straight_down_looks_at_target() {
        let c = EarthCamera { distance_m: 1000.0, ..EarthCamera::new(6_371_000.0) };
        let (lat, lon) = c.pick(0.5, 0.5, 1.5).unwrap();
        assert!((lat - c.lat).abs() < 1e-9 && (lon - c.lon).abs() < 1e-9);
        assert!((c.altitude_m() - 1000.0).abs() < 1e-3);
    }

    #[test]
    fn tilted_centre_still_hits_target() {
        let mut c = EarthCamera { distance_m: 20_000.0, ..EarthCamera::new(6_371_000.0) };
        c.orbit(0.7, 1.0);
        assert!(c.tilt > 0.5);
        let (lat, lon) = c.pick(0.5, 0.5, 1.5).unwrap();
        assert!((lat - c.lat).abs() < 1e-7 && (lon - c.lon).abs() < 1e-7);
    }

    #[test]
    fn target_projects_to_screen_centre() {
        let mut c = EarthCamera { distance_m: 50_000.0, ..EarthCamera::new(6_371_000.0) };
        c.orbit(1.2, 0.8);
        let m = c.view_proj_f64(1.5, 1000.0);
        let p = m * (c.target() - c.eye()).extend(1.0);
        assert!((p.x / p.w).abs() < 1e-9 && (p.y / p.w).abs() < 1e-9);
        assert!(p.z / p.w > 0.0 && p.z / p.w < 1.0);
    }

    #[test]
    fn pan_moves_target_north_when_dragging_down() {
        let mut c = EarthCamera { distance_m: 100_000.0, ..EarthCamera::new(6_371_000.0) };
        let lat0 = c.lat;
        c.pan(0.0, 100.0, 800.0);
        assert!(c.lat > lat0);
    }
}
