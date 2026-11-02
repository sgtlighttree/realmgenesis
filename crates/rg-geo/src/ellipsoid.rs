use glam::DVec3;

/// An oblate ellipsoid of revolution. `flattening = 0` is a sphere.
///
/// One code path serves both: there's no separate "sphere mode".
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Ellipsoid {
    /// Equatorial radius, metres.
    pub a: f64,
    /// Flattening f = (a - b) / a. Earth ≈ 1/298.257.
    pub flattening: f64,
}

impl Ellipsoid {
    pub const fn sphere(radius_m: f64) -> Self {
        Self { a: radius_m, flattening: 0.0 }
    }

    /// Polar radius b.
    pub fn b(&self) -> f64 {
        self.a * (1.0 - self.flattening)
    }

    /// First eccentricity squared, e² = f(2 - f).
    pub fn e2(&self) -> f64 {
        self.flattening * (2.0 - self.flattening)
    }

    /// Geodetic (lat, lon, height above the ellipsoid) → ECEF metres.
    pub fn geodetic_to_ecef(&self, lat: f64, lon: f64, h: f64) -> DVec3 {
        let (sl, cl) = lat.sin_cos();
        let (so, co) = lon.sin_cos();
        // Prime vertical radius of curvature.
        let n = self.a / (1.0 - self.e2() * sl * sl).sqrt();
        DVec3::new(
            (n + h) * cl * co,
            (n + h) * cl * so,
            (n * (1.0 - self.e2()) + h) * sl,
        )
    }

    /// ECEF metres → geodetic (lat, lon, h). Bowring's method with three iterations,
    /// which is sub-millimetre for any height a planet generator cares about.
    pub fn ecef_to_geodetic(&self, p: DVec3) -> (f64, f64, f64) {
        let lon = p.y.atan2(p.x);
        let r = (p.x * p.x + p.y * p.y).sqrt();
        let e2 = self.e2();
        if e2 == 0.0 {
            let lat = p.z.atan2(r);
            return (lat, lon, p.length() - self.a);
        }
        let b = self.b();
        let ep2 = (self.a * self.a - b * b) / (b * b);
        let mut beta = (self.a * p.z).atan2(b * r);
        let mut lat = 0.0;
        for _ in 0..3 {
            let (sb, cb) = beta.sin_cos();
            lat = (p.z + ep2 * b * sb * sb * sb).atan2(r - e2 * self.a * cb * cb * cb);
            beta = ((1.0 - self.flattening) * lat.sin()).atan2(lat.cos());
        }
        let sl = lat.sin();
        let n = self.a / (1.0 - e2 * sl * sl).sqrt();
        let h = if lat.cos().abs() > 1e-10 {
            r / lat.cos() - n
        } else {
            p.z.abs() - b
        };
        (lat, lon, h)
    }

    /// Surface point (h = 0) under a geocentric unit direction. Used to place cube-sphere
    /// samples, whose grid is defined in direction space.
    pub fn surface_along(&self, dir: DVec3) -> DVec3 {
        let b = self.b();
        // Solve |(x/a, y/a, z/b)| = 1 for the point t·dir.
        let q = (dir.x * dir.x + dir.y * dir.y) / (self.a * self.a) + dir.z * dir.z / (b * b);
        dir / q.sqrt()
    }

    /// Geodetic latitude of the surface point under a geocentric direction.
    pub fn geodetic_lat_of_dir(&self, dir: DVec3) -> f64 {
        let (lat, _, _) = self.ecef_to_geodetic(self.surface_along(dir));
        lat
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EARTH: Ellipsoid = Ellipsoid { a: 6_378_137.0, flattening: 1.0 / 298.257_223_563 };

    #[test]
    fn roundtrip_ellipsoid() {
        for &(lat, lon, h) in &[(0.0, 0.0, 0.0), (0.7, -2.0, 1200.0), (-1.5, 3.0, -8000.0), (1.57, 0.1, 50.0)] {
            let p = EARTH.geodetic_to_ecef(lat, lon, h);
            let (lat2, lon2, h2) = EARTH.ecef_to_geodetic(p);
            assert!((lat - lat2).abs() < 1e-9, "lat {lat} vs {lat2}");
            assert!((lon - lon2).abs() < 1e-9);
            assert!((h - h2).abs() < 1e-3, "h {h} vs {h2}");
        }
    }

    #[test]
    fn sphere_is_flattening_zero() {
        let s = Ellipsoid::sphere(4_500_000.0);
        let p = s.geodetic_to_ecef(0.3, 1.0, 10.0);
        assert!((p.length() - 4_500_010.0).abs() < 1e-6);
        let (lat, _, h) = s.ecef_to_geodetic(p);
        assert!((lat - 0.3).abs() < 1e-12 && (h - 10.0).abs() < 1e-6);
    }

    #[test]
    fn surface_along_lies_on_ellipsoid() {
        let d = DVec3::new(0.3, -0.5, 0.8).normalize();
        let p = EARTH.surface_along(d);
        let (_, _, h) = EARTH.ecef_to_geodetic(p);
        assert!(h.abs() < 1e-3);
    }
}
