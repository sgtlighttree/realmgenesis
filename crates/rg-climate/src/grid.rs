//! `LatLonGrid<T>`: a plain equirectangular (latitude/longitude) raster.
//!
//! Layout, shared by every climate field:
//! - cell **centred**: a grid of `width x height` cells covers the whole sphere, so one cell is
//!   `360/width` degrees wide and `180/height` degrees tall;
//! - **row 0 is the north edge** (latitude `+90 - dlat/2`), rows run south;
//! - **column 0 is the west edge** (longitude `-180 + dlon/2`), columns run east and wrap.
//!
//! Latitude and longitude are in radians throughout the crate (see `rg-geo` conventions).

use glam::DVec3;
use rg_geo::CubeField;
use std::f64::consts::PI;

#[derive(Clone, Debug, PartialEq)]
pub struct LatLonGrid<T> {
    pub width: u32,
    pub height: u32,
    /// Row-major: `data[y * width + x]`.
    pub data: Vec<T>,
}

impl<T> LatLonGrid<T> {
    /// Build a grid by calling `f(lat_rad, lon_rad)` at every cell centre.
    pub fn from_fn(width: u32, height: u32, mut f: impl FnMut(f64, f64) -> T) -> Self {
        let mut data = Vec::with_capacity((width * height) as usize);
        for y in 0..height {
            let lat = Self::lat_of_row_for(height, y);
            for x in 0..width {
                data.push(f(lat, Self::lon_of_col_for(width, x)));
            }
        }
        Self { width, height, data }
    }

    /// Wrap an existing buffer (must hold `width * height` values).
    pub fn from_vec(width: u32, height: u32, data: Vec<T>) -> Self {
        assert_eq!(data.len(), (width * height) as usize, "grid data has the wrong length");
        Self { width, height, data }
    }

    #[inline]
    pub fn index(&self, x: u32, y: u32) -> usize {
        y as usize * self.width as usize + x as usize
    }

    #[inline]
    pub fn get(&self, x: u32, y: u32) -> &T {
        &self.data[self.index(x, y)]
    }

    /// Latitude (radians) of the centre of row `y`.
    pub fn lat_of_row(&self, y: u32) -> f64 {
        Self::lat_of_row_for(self.height, y)
    }

    /// Longitude (radians, -pi..pi) of the centre of column `x`.
    pub fn lon_of_col(&self, x: u32) -> f64 {
        Self::lon_of_col_for(self.width, x)
    }

    fn lat_of_row_for(height: u32, y: u32) -> f64 {
        let dlat = PI / height as f64;
        PI / 2.0 - (y as f64 + 0.5) * dlat
    }

    fn lon_of_col_for(width: u32, x: u32) -> f64 {
        let dlon = 2.0 * PI / width as f64;
        -PI + (x as f64 + 0.5) * dlon
    }
}

impl<T: Clone> LatLonGrid<T> {
    /// A grid with every cell set to `value`.
    pub fn filled(width: u32, height: u32, value: T) -> Self {
        Self { width, height, data: vec![value; (width * height) as usize] }
    }
}

impl LatLonGrid<f32> {
    /// Bilinear sample at (lat, lon) in radians. Longitude wraps; latitude clamps at the poles.
    pub fn sample(&self, lat: f64, lon: f64) -> f32 {
        let (w, h) = (self.width as i64, self.height as i64);
        let dlat = PI / h as f64;
        let dlon = 2.0 * PI / w as f64;
        // Fractional cell coordinates where integer values are cell centres.
        let fy = (PI / 2.0 - lat) / dlat - 0.5;
        let fx = (lon + PI) / dlon - 0.5;
        let y0 = fy.floor();
        let x0 = fx.floor();
        let ty = (fy - y0) as f32;
        let tx = (fx - x0) as f32;
        let row = |y: i64| y.clamp(0, h - 1) as usize;
        let col = |x: i64| x.rem_euclid(w) as usize;
        let at = |x: i64, y: i64| self.data[row(y) * self.width as usize + col(x)];
        let (x0, y0) = (x0 as i64, y0 as i64);
        let top = at(x0, y0) * (1.0 - tx) + at(x0 + 1, y0) * tx;
        let bot = at(x0, y0 + 1) * (1.0 - tx) + at(x0 + 1, y0 + 1) * tx;
        top * (1.0 - ty) + bot * ty
    }

    /// Resample a cube-sphere field (e.g. the tectonics height field) onto a lat/lon grid.
    pub fn from_cube_field(field: &CubeField<f32>, width: u32, height: u32) -> Self {
        Self::from_fn(width, height, |lat, lon| {
            let (sl, cl) = lat.sin_cos();
            let (so, co) = lon.sin_cos();
            field.sample(DVec3::new(cl * co, cl * so, sl))
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cell_centres_follow_the_documented_layout() {
        let g = LatLonGrid::<f32>::filled(360, 180, 0.0);
        assert!((g.lat_of_row(0).to_degrees() - 89.5).abs() < 1e-9);
        assert!((g.lat_of_row(179).to_degrees() + 89.5).abs() < 1e-9);
        assert!((g.lon_of_col(0).to_degrees() + 179.5).abs() < 1e-9);
        assert!((g.lon_of_col(359).to_degrees() - 179.5).abs() < 1e-9);
    }

    #[test]
    fn sample_is_exact_at_centres_and_wraps() {
        let g = LatLonGrid::from_fn(36, 18, |lat, lon| (lat.sin() + 0.3 * lon.cos()) as f32);
        for &(x, y) in &[(0u32, 0u32), (5, 7), (35, 17), (20, 9)] {
            let v = g.sample(g.lat_of_row(y), g.lon_of_col(x));
            assert!((v - g.get(x, y)).abs() < 1e-5);
        }
        // Just across the antimeridian we interpolate between the last and first column.
        let a = g.sample(0.1, PI - 0.01);
        let b = g.sample(0.1, -PI + 0.01);
        assert!((a - b).abs() < 0.05);
    }

    #[test]
    fn cube_field_resample_matches_direction() {
        let f = CubeField::<f32>::from_fn(32, |d| d.z as f32);
        let g = LatLonGrid::from_cube_field(&f, 72, 36);
        let v = *g.get(10, 0);
        assert!(v > 0.95, "near-north-pole cell should have z ~ 1, got {v}");
    }
}
