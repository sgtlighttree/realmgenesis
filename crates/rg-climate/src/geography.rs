//! Grid geometry and land/ocean geography helpers: cell sizes in kilometres, distance to the
//! coast, coast orientation along rows, and blurring.

use crate::grid::LatLonGrid;
use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::f64::consts::PI;

/// Cell sizes of a lat/lon grid on a planet of a given radius.
#[derive(Clone, Debug)]
pub struct Geometry {
    pub w: usize,
    pub h: usize,
    pub radius_km: f64,
    /// Cell height in radians of latitude.
    pub dlat: f64,
    /// Cell width in radians of longitude.
    pub dlon: f64,
    /// North-south size of a cell, km.
    pub ds_km: f64,
    /// Latitude of each row, radians.
    pub lat: Vec<f64>,
    /// Latitude of each row, degrees.
    pub lat_deg: Vec<f64>,
    /// East-west size of a cell on each row, km (shrinks toward the poles).
    pub dx_km: Vec<f64>,
}

impl Geometry {
    pub fn new(width: u32, height: u32, radius_m: f64) -> Self {
        let radius_km = radius_m / 1000.0;
        let dlat = PI / height as f64;
        let dlon = 2.0 * PI / width as f64;
        let lat: Vec<f64> = (0..height).map(|y| PI / 2.0 - (y as f64 + 0.5) * dlat).collect();
        let dx_km = lat.iter().map(|l| radius_km * l.cos() * dlon).collect();
        Self {
            w: width as usize,
            h: height as usize,
            radius_km,
            dlat,
            dlon,
            ds_km: radius_km * dlat,
            lat_deg: lat.iter().map(|l| l.to_degrees()).collect(),
            lat,
            dx_km,
        }
    }
}

/// Great-circle-ish distance (km) from every land cell to the nearest ocean cell; 0 for ocean.
/// Dijkstra over the 8-neighbourhood with edge lengths measured in km on the sphere, so cells
/// near the poles are correctly close together. The result is measured to the coast line,
/// i.e. half a cell is subtracted. If the planet has no ocean at all every cell gets 6000 km.
pub fn distance_to_ocean_km(land: &[bool], geo: &Geometry) -> Vec<f32> {
    let (w, h) = (geo.w, geo.h);
    let mut dist = vec![f32::INFINITY; w * h];
    // Keys are f32 bit patterns: for non-negative floats, integer order equals numeric order.
    let mut heap: BinaryHeap<Reverse<(u32, u32)>> = BinaryHeap::new();
    for (i, &is_land) in land.iter().enumerate() {
        if !is_land {
            dist[i] = 0.0;
            heap.push(Reverse((0f32.to_bits(), i as u32)));
        }
    }
    while let Some(Reverse((key, idx))) = heap.pop() {
        let idx = idx as usize;
        let d = f32::from_bits(key);
        if d > dist[idx] {
            continue;
        }
        let (x, y) = (idx % w, idx / w);
        for dy in -1i32..=1 {
            let ny = y as i32 + dy;
            if ny < 0 || ny >= h as i32 {
                continue;
            }
            for dx in -1i32..=1 {
                if dx == 0 && dy == 0 {
                    continue;
                }
                let nx = (x as i32 + dx).rem_euclid(w as i32) as usize;
                let ny = ny as usize;
                let row_dx = geo.dx_km[y].min(geo.dx_km[ny]).max(geo.ds_km * 0.02);
                let step = ((dx as f64 * row_dx).powi(2) + (dy as f64 * geo.ds_km).powi(2)).sqrt() as f32;
                let nd = d + step;
                let ni = ny * w + nx;
                if nd < dist[ni] {
                    dist[ni] = nd;
                    heap.push(Reverse((nd.to_bits(), ni as u32)));
                }
            }
        }
    }
    for (i, d) in dist.iter_mut().enumerate() {
        *d = if land[i] {
            if d.is_finite() { (*d - 0.5 * geo.ds_km as f32).max(0.0) } else { 6000.0 }
        } else {
            0.0
        };
    }
    dist
}

/// For every cell, the distance in km travelling **west** and **east** along its row until a
/// cell of the *opposite* kind is met (land if the cell is ocean, ocean if the cell is land).
/// Rows wrap around the globe. `f32::INFINITY` when the whole row is one kind.
///
/// This is the cheap "coast orientation" test: an ocean cell with land close to the east sits
/// on the eastern side of its basin (where subtropical upwelling and cold currents live); a
/// land cell with ocean close to the west is on a west coast.
pub fn row_runs_km(land: &[bool], geo: &Geometry) -> (Vec<f32>, Vec<f32>) {
    let (w, h) = (geo.w, geo.h);
    let mut west = vec![f32::INFINITY; w * h];
    let mut east = vec![f32::INFINITY; w * h];
    for y in 0..h {
        let row = &land[y * w..(y + 1) * w];
        let cell_km = geo.dx_km[y] as f32;
        // Sweep eastward twice around the circle; the second lap sees every cell's western
        // neighbours. `last[k]` is the unwrapped column of the last cell of kind k seen so far.
        let mut last: [Option<usize>; 2] = [None, None];
        for k in 0..2 * w {
            let kind = row[k % w] as usize;
            if k >= w
                && let Some(p) = last[1 - kind]
            {
                west[y * w + k - w] = (k - p) as f32 * cell_km;
            }
            last[kind] = Some(k);
        }
        let mut last: [Option<usize>; 2] = [None, None];
        for k in (0..2 * w).rev() {
            let kind = row[k % w] as usize;
            if k < w
                && let Some(p) = last[1 - kind]
            {
                east[y * w + k] = (p - k) as f32 * cell_km;
            }
            last[kind] = Some(k);
        }
    }
    (west, east)
}

/// Separable box blur: longitude wraps, latitude clamps at the poles. `radius` is in cells.
pub fn box_blur(src: &[f32], w: usize, h: usize, radius: usize) -> Vec<f32> {
    if radius == 0 {
        return src.to_vec();
    }
    let r = radius;
    let norm = 1.0 / (2 * r + 1) as f32;
    // Horizontal pass (wrapping): copy each row with `r` wrapped cells on both sides, then slide.
    let mut tmp = vec![0.0f32; w * h];
    let mut ext = vec![0.0f32; w + 2 * r];
    for y in 0..h {
        let row = &src[y * w..(y + 1) * w];
        for (k, e) in ext.iter_mut().enumerate() {
            *e = row[(k + w * (r / w + 1) - r) % w];
        }
        let mut sum: f32 = ext[..2 * r + 1].iter().sum();
        for x in 0..w {
            tmp[y * w + x] = sum * norm;
            if x + 1 < w {
                sum += ext[x + 2 * r + 1] - ext[x];
            }
        }
    }
    // Vertical pass (clamped), done row by row so memory access stays sequential.
    let mut out = vec![0.0f32; w * h];
    let mut sum = vec![0.0f32; w];
    let row_of = |j: isize| j.clamp(0, h as isize - 1) as usize;
    for j in -(r as isize)..=(r as isize) {
        let src_row = &tmp[row_of(j) * w..][..w];
        for (s, v) in sum.iter_mut().zip(src_row) {
            *s += v;
        }
    }
    for y in 0..h {
        for (o, s) in out[y * w..(y + 1) * w].iter_mut().zip(&sum) {
            *o = s * norm;
        }
        let add = &tmp[row_of(y as isize + r as isize + 1) * w..][..w];
        let sub = &tmp[row_of(y as isize - r as isize) * w..][..w];
        for ((s, a), b) in sum.iter_mut().zip(add).zip(sub) {
            *s += a - b;
        }
    }
    out
}

/// Land mask from a height grid: `height >= 0` is land.
pub fn land_mask(height_m: &LatLonGrid<f32>) -> Vec<bool> {
    height_m.data.iter().map(|&h| h >= 0.0).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn distance_to_ocean_grows_inland() {
        // A 20-column-wide continent in the middle of a 72x36 world.
        let (w, h) = (72usize, 36usize);
        let geo = Geometry::new(w as u32, h as u32, 6_371_000.0);
        let land: Vec<bool> = (0..w * h).map(|i| (26..46).contains(&(i % w))).collect();
        let d = distance_to_ocean_km(&land, &geo);
        let row = 18 * w;
        assert_eq!(d[row + 10], 0.0);
        assert!(d[row + 26] < d[row + 30] && d[row + 30] < d[row + 36]);
        // 10 cells of 5 degrees at the equator is about 5500 km.
        assert!((d[row + 36] - 10.0 * 555.0).abs() < 800.0, "got {}", d[row + 36]);
    }

    #[test]
    fn row_runs_find_west_and_east_coasts() {
        let (w, h) = (72usize, 4usize);
        let geo = Geometry::new(w as u32, h as u32, 6_371_000.0);
        let land: Vec<bool> = (0..w * h).map(|i| (26..46).contains(&(i % w))).collect();
        let (west, east) = row_runs_km(&land, &geo);
        // Ocean cell just west of the continent: land is 1 cell east of it.
        assert!((east[25] - geo.dx_km[0] as f32).abs() < 1.0);
        // Land cell 26 is the western-most: ocean is 1 cell to its west.
        assert!((west[26] - geo.dx_km[0] as f32).abs() < 1.0);
        assert!(west[40] > west[26]);
    }

    #[test]
    fn box_blur_preserves_constants_and_mean() {
        let (w, h) = (16usize, 8usize);
        let flat = vec![3.0f32; w * h];
        assert!(box_blur(&flat, w, h, 2).iter().all(|v| (v - 3.0).abs() < 1e-5));
        let mut spike = vec![0.0f32; w * h];
        spike[4 * w + 8] = 25.0;
        let b = box_blur(&spike, w, h, 2);
        assert!((b.iter().sum::<f32>() - 25.0).abs() < 1e-3);
    }
}
