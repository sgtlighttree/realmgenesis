use glam::{DMat4, DVec3, Mat4};


/// An orbit camera around the planet centre. Units are planet radii, so the surface sits
/// at distance 1. f64 throughout. M1 adds camera-relative tile origins for close zoom.
#[derive(Clone, Copy, Debug)]
pub struct OrbitCamera {
    /// Longitude of the camera, radians.
    pub yaw: f64,
    /// Latitude of the camera, radians, clamped short of the poles.
    pub pitch: f64,
    /// Distance from the planet centre, in radii (> 1).
    pub distance: f64,
    /// Vertical field of view, radians.
    pub fov_y: f64,
}

impl Default for OrbitCamera {
    fn default() -> Self {
        Self { yaw: 0.6, pitch: 0.35, distance: 3.2, fov_y: 40f64.to_radians() }
    }
}

impl OrbitCamera {
    pub const MIN_DISTANCE: f64 = 1.02;
    pub const MAX_DISTANCE: f64 = 12.0;

    pub fn eye(&self) -> DVec3 {
        let (sp, cp) = self.pitch.sin_cos();
        let (sy, cy) = self.yaw.sin_cos();
        DVec3::new(cp * cy, cp * sy, sp) * self.distance
    }

    /// Drag by a screen delta in pixels. Rotation slows near the surface so the ground
    /// tracks the cursor at roughly the same speed at any zoom.
    pub fn drag(&mut self, dx: f64, dy: f64, viewport_h: f64) {
        let altitude = (self.distance - 1.0).max(1e-4);
        let k = (altitude / self.distance.max(1.0)) * self.fov_y / viewport_h.max(1.0) * 1.6;
        self.yaw -= dx * k;
        self.pitch = (self.pitch + dy * k).clamp(-1.55, 1.55);
    }

    /// Zoom by a factor on altitude (not distance), so close zoom stays controllable.
    pub fn zoom(&mut self, factor: f64) {
        let alt = (self.distance - 1.0) * factor;
        self.distance = (1.0 + alt).clamp(Self::MIN_DISTANCE, Self::MAX_DISTANCE);
    }

    /// View-projection with an infinite far plane and **reversed Z** (near → 1, far → 0),
    /// which keeps depth precision across planetary distances. Depth test: `Greater`.
    pub fn view_proj(&self, aspect: f64) -> Mat4 {
        let eye = self.eye();
        let view = look_at_rh(eye, DVec3::ZERO, DVec3::Z);
        let near = ((self.distance - 1.0) * 0.1).max(1e-5);
        let proj = perspective_infinite_reverse_rh(self.fov_y, aspect, near);
        (proj * view).as_mat4()
    }
}

/// Right-handed look-at view matrix (glam 0.34 no longer ships one).
pub fn look_at_rh(eye: DVec3, target: DVec3, up: DVec3) -> DMat4 {
    let f = (target - eye).normalize();
    let s = f.cross(up).normalize();
    let u = s.cross(f);
    DMat4::from_cols(
        glam::DVec4::new(s.x, u.x, -f.x, 0.0),
        glam::DVec4::new(s.y, u.y, -f.y, 0.0),
        glam::DVec4::new(s.z, u.z, -f.z, 0.0),
        glam::DVec4::new(-s.dot(eye), -u.dot(eye), f.dot(eye), 1.0),
    )
}

/// Right-handed infinite perspective with reversed Z, mapping depth to [0, 1] as WebGPU
/// expects: near plane → 1, infinity → 0.
pub fn perspective_infinite_reverse_rh(fov_y: f64, aspect: f64, near: f64) -> DMat4 {
    let f = 1.0 / (0.5 * fov_y).tan();
    DMat4::from_cols(
        glam::DVec4::new(f / aspect, 0.0, 0.0, 0.0),
        glam::DVec4::new(0.0, f, 0.0, 0.0),
        glam::DVec4::new(0.0, 0.0, 0.0, -1.0),
        glam::DVec4::new(0.0, 0.0, near, 0.0),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reversed_z_maps_near_to_one_and_far_to_zero() {
        let p = perspective_infinite_reverse_rh(1.0, 1.0, 0.1);
        let near = p * glam::DVec4::new(0.0, 0.0, -0.1, 1.0);
        let far = p * glam::DVec4::new(0.0, 0.0, -1e9, 1.0);
        assert!((near.z / near.w - 1.0).abs() < 1e-9);
        assert!((far.z / far.w).abs() < 1e-9);
    }

    #[test]
    fn planet_centre_is_in_front() {
        let c = OrbitCamera::default();
        let v = look_at_rh(c.eye(), DVec3::ZERO, DVec3::Z) * glam::DVec4::new(0.0, 0.0, 0.0, 1.0);
        assert!(v.z < 0.0 && (v.z + c.distance).abs() < 1e-9);
    }
}
