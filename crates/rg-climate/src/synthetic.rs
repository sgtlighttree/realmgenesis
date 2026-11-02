//! A crude, hand-made Earth-like planet for tests, examples and the app's "does the climate
//! look sane" check. Continents are coarse lon/lat polygons (a few dozen points each, good to a
//! few degrees), mountain ranges are Gaussian ridges. It is *not* a data set: it only needs to
//! put land in roughly the right places so Köppen zones can be eyeballed against the real map.

use crate::grid::LatLonGrid;

type Poly = &'static [(f64, f64)];

const NORTH_AMERICA: Poly = &[
    (-168.0, 66.0), (-162.0, 70.0), (-141.0, 70.0), (-125.0, 70.0), (-110.0, 68.0), (-95.0, 68.0),
    (-85.0, 70.0), (-82.0, 67.0), (-78.0, 62.0), (-70.0, 60.0), (-65.0, 60.0), (-62.0, 57.0),
    (-56.0, 52.0), (-60.0, 47.0), (-66.0, 44.0), (-70.0, 43.0), (-74.0, 40.0), (-76.0, 35.0),
    (-81.0, 31.0), (-80.0, 25.5), (-82.5, 28.0), (-84.0, 30.0), (-89.0, 30.0), (-94.0, 29.5),
    (-97.0, 27.0), (-97.5, 22.0), (-95.0, 18.5), (-91.0, 18.5), (-90.0, 21.0), (-87.0, 21.5),
    (-88.0, 16.0), (-84.0, 15.5), (-83.0, 10.0), (-79.0, 9.0), (-77.5, 8.0), (-80.0, 7.5),
    (-85.5, 9.5), (-87.0, 13.0), (-92.0, 14.5), (-96.0, 15.8), (-105.0, 20.0), (-106.0, 23.0),
    (-110.0, 23.0), (-112.0, 26.0), (-115.0, 30.0), (-117.0, 32.5), (-121.0, 35.0), (-124.0, 40.0),
    (-124.5, 46.0), (-123.0, 48.5), (-128.0, 51.0), (-134.0, 56.0), (-140.0, 60.0), (-148.0, 60.5),
    (-152.0, 59.0), (-158.0, 56.5), (-163.0, 55.0), (-158.0, 58.5), (-162.0, 60.0), (-165.0, 62.0),
];

const SOUTH_AMERICA: Poly = &[
    (-77.5, 8.5), (-72.0, 12.0), (-65.0, 10.5), (-60.0, 8.0), (-52.0, 5.0), (-50.0, 0.0),
    (-44.0, -2.5), (-35.0, -5.5), (-35.0, -9.0), (-39.0, -14.0), (-39.0, -18.0), (-41.0, -22.0),
    (-48.0, -26.0), (-49.0, -29.0), (-53.0, -34.0), (-57.0, -35.0), (-57.0, -38.0), (-62.0, -39.0),
    (-65.0, -41.0), (-64.0, -46.0), (-67.0, -47.0), (-69.0, -51.0), (-68.5, -53.0), (-71.0, -54.0),
    (-74.0, -52.0), (-75.0, -47.0), (-73.0, -42.0), (-73.5, -37.0), (-71.5, -30.0), (-70.5, -18.0),
    (-76.0, -14.0), (-81.0, -6.0), (-80.0, -2.0), (-80.0, 0.0), (-78.0, 2.0), (-77.5, 6.0),
];

const AFRICA: Poly = &[
    (-17.0, 21.0), (-17.0, 14.5), (-13.0, 9.0), (-8.0, 4.5), (-4.0, 5.0), (2.0, 6.0), (6.0, 4.5),
    (9.0, 4.0), (9.5, -1.0), (12.0, -5.0), (13.0, -10.0), (12.5, -17.0), (14.5, -23.0),
    (18.0, -32.0), (20.0, -35.0), (25.0, -34.0), (30.0, -31.0), (32.5, -26.0), (35.0, -20.0),
    (40.5, -15.0), (40.0, -10.0), (39.0, -5.0), (41.0, -1.0), (44.0, 3.0), (48.0, 8.0),
    (51.0, 11.5), (43.5, 12.0), (39.0, 16.0), (37.0, 21.0), (33.0, 28.0), (32.5, 31.0),
    (25.0, 32.0), (20.0, 31.0), (15.0, 32.0), (10.0, 34.0), (10.0, 37.0), (0.0, 36.0),
    (-5.0, 35.7), (-10.0, 31.0), (-13.0, 27.5), (-16.0, 23.5),
];

const EURASIA: Poly = &[
    (-9.0, 37.0), (-9.0, 43.0), (-2.0, 43.5), (-1.5, 46.0), (-4.5, 48.5), (-1.5, 49.5), (2.0, 51.0),
    (5.0, 53.0), (8.0, 54.0), (8.0, 57.0), (5.5, 58.5), (5.0, 62.0), (10.0, 64.0), (14.0, 67.5),
    (19.0, 70.0), (26.0, 71.0), (31.0, 70.0), (41.0, 67.0), (33.0, 66.0), (40.0, 64.5), (44.0, 66.0),
    (54.0, 68.5), (60.0, 69.0), (68.0, 69.0), (73.0, 72.0), (80.0, 73.0), (88.0, 75.5), (100.0, 77.0),
    (105.0, 77.5), (113.0, 74.0), (128.0, 73.0), (140.0, 72.5), (150.0, 71.5), (160.0, 69.5),
    (170.0, 70.0), (180.0, 69.0), (180.0, 65.0), (178.0, 64.5), (172.0, 64.5), (170.0, 60.0),
    (163.0, 59.5), (160.0, 53.0), (156.0, 51.0), (156.0, 57.0), (162.0, 61.0), (155.0, 59.5),
    (142.0, 59.0), (136.0, 54.5), (141.0, 53.0), (140.0, 48.0), (135.0, 43.0), (132.0, 43.0),
    (129.0, 41.0), (129.0, 36.0), (127.0, 34.5), (126.0, 37.5), (125.0, 39.5), (121.5, 40.5),
    (118.0, 39.0), (122.0, 37.0), (119.0, 34.5), (122.0, 31.0), (122.0, 29.0), (119.5, 25.5),
    (114.0, 22.5), (110.0, 21.0), (108.0, 21.5), (106.0, 19.0), (109.0, 15.0), (109.0, 11.5),
    (105.0, 8.5), (103.0, 10.5), (100.5, 13.5), (99.5, 10.0), (101.0, 6.5), (103.5, 1.5),
    (100.5, 3.5), (98.5, 8.0), (98.0, 12.0), (97.0, 16.5), (94.0, 18.0), (92.0, 21.5), (90.0, 22.0),
    (87.0, 21.0), (80.0, 15.5), (80.0, 10.0), (77.5, 8.0), (75.0, 12.0), (73.0, 17.0), (72.5, 21.0),
    (68.5, 23.5), (66.5, 25.3), (61.0, 25.2), (57.0, 25.8), (56.5, 27.0), (51.0, 28.0), (48.5, 30.0),
    (50.0, 26.5), (51.5, 24.5), (54.0, 24.0), (56.5, 25.5), (59.5, 22.5), (55.0, 17.0), (52.0, 16.0),
    (44.0, 12.8), (43.0, 16.0), (39.0, 21.0), (35.0, 28.0), (34.5, 29.5), (35.0, 32.0), (36.0, 34.0),
    (36.5, 36.8), (32.0, 36.2), (30.0, 36.5), (27.0, 37.0), (26.5, 39.0), (26.0, 40.5), (23.0, 38.0),
    (21.5, 37.0), (20.0, 40.0), (19.5, 42.0), (16.0, 44.0), (13.5, 45.5), (12.3, 44.5), (14.0, 42.0),
    (16.5, 40.5), (18.5, 40.0), (16.0, 38.0), (15.7, 40.0), (12.0, 42.0), (10.5, 43.5), (8.5, 44.3),
    (6.0, 43.0), (3.0, 43.2), (3.0, 42.0), (0.0, 39.0), (-0.5, 38.2), (-2.0, 36.7), (-5.5, 36.0),
];

const AUSTRALIA: Poly = &[
    (114.0, -22.0), (114.0, -26.0), (115.0, -34.0), (118.0, -35.0), (123.0, -34.0), (129.0, -31.5),
    (135.0, -34.8), (138.0, -35.0), (140.0, -38.0), (146.0, -39.0), (150.0, -37.0), (153.0, -31.0),
    (153.0, -25.0), (146.0, -19.0), (145.0, -15.0), (142.5, -10.7), (141.0, -17.0), (137.0, -16.0),
    (135.0, -12.0), (131.0, -11.3), (127.0, -14.0), (122.0, -17.0), (118.0, -20.5),
];

const GREENLAND: Poly = &[
    (-73.0, 78.0), (-60.0, 82.0), (-30.0, 83.5), (-20.0, 81.0), (-18.0, 76.0), (-22.0, 70.0),
    (-32.0, 68.0), (-42.0, 60.0), (-50.0, 62.0), (-54.0, 67.0), (-58.0, 75.0), (-68.0, 76.5),
];

/// Smaller islands and inland seas as ellipses: (lon, lat, half-width lon deg, half-height lat deg).
const ISLANDS: &[(f64, f64, f64, f64)] = &[
    (-2.5, 54.0, 3.0, 4.5),    // Great Britain
    (-8.0, 53.2, 1.8, 1.6),    // Ireland
    (-19.0, 65.0, 3.0, 1.5),   // Iceland
    (138.0, 36.0, 3.5, 2.5),   // Honshu
    (143.0, 43.5, 2.0, 1.8),   // Hokkaido
    (101.5, -0.3, 5.5, 1.5),   // Sumatra
    (114.0, 0.5, 3.5, 3.5),    // Borneo
    (121.0, -2.0, 1.8, 2.0),   // Sulawesi
    (110.5, -7.5, 3.0, 0.8),   // Java
    (141.0, -5.5, 6.5, 2.0),   // New Guinea
    (122.0, 12.5, 1.8, 3.5),   // Philippines
    (47.0, -19.5, 2.2, 6.5),   // Madagascar
    (80.7, 7.8, 0.8, 1.2),     // Sri Lanka
    (-79.0, 21.8, 4.5, 1.0),   // Cuba
    (-72.5, 19.0, 2.5, 1.0),   // Hispaniola
    (172.5, -41.5, 2.0, 2.5),  // New Zealand north island + south-west
    (169.5, -44.5, 2.5, 2.2),  // New Zealand south island
    (146.5, -42.0, 1.5, 1.0),  // Tasmania
    (-75.0, 70.0, 8.0, 3.0),   // Baffin
    (-100.0, 73.0, 8.0, 3.0),  // Arctic archipelago
    (55.0, 74.0, 6.0, 2.0),    // Novaya Zemlya
    (142.0, 52.0, 1.0, 3.0),   // Sakhalin
    (12.0, 78.5, 4.0, 1.5),    // Svalbard
];

/// Inland seas and bays carved out of the continents (same tuple layout).
const SEAS: &[(f64, f64, f64, f64)] = &[
    (-85.0, 60.0, 7.0, 5.0),  // Hudson Bay
    (34.0, 43.3, 6.0, 1.8),   // Black Sea
    (51.0, 41.5, 2.0, 5.0),   // Caspian
    (20.0, 57.8, 3.5, 2.5),   // Baltic
    (21.5, 62.5, 1.5, 2.5),   // Gulf of Bothnia
    (-90.0, 25.0, 6.0, 3.0),  // central Gulf of Mexico (keeps the gulf open)
    (4.0, 54.5, 2.5, 1.2),    // North Sea south
];

/// Gaussian mountain ranges: polyline of (lon, lat) points, half-width in degrees, height m.
struct Ridge {
    path: &'static [(f64, f64)],
    width_deg: f64,
    height_m: f64,
}

const RIDGES: &[Ridge] = &[
    Ridge { path: &[(-148.0, 63.5), (-136.0, 60.5)], width_deg: 2.0, height_m: 2500.0 },
    Ridge { path: &[(-124.0, 58.0), (-117.0, 52.0), (-111.0, 44.0), (-106.0, 36.0)], width_deg: 3.0, height_m: 2000.0 },
    Ridge { path: &[(-120.5, 48.5), (-119.5, 38.0)], width_deg: 1.4, height_m: 1500.0 },
    Ridge { path: &[(-83.0, 35.0), (-72.0, 44.0)], width_deg: 1.5, height_m: 800.0 },
    Ridge { path: &[(-104.0, 28.0), (-101.0, 19.0)], width_deg: 2.5, height_m: 1700.0 },
    Ridge { path: &[(-72.0, 10.0), (-78.0, 0.0), (-77.0, -10.0), (-69.0, -18.0), (-69.5, -30.0), (-71.5, -40.0), (-73.0, -50.0)], width_deg: 2.2, height_m: 3600.0 },
    Ridge { path: &[(-47.0, -14.0), (-45.0, -22.0)], width_deg: 4.0, height_m: 800.0 },
    Ridge { path: &[(-3.0, 43.0), (3.0, 42.6)], width_deg: 1.0, height_m: 1500.0 },
    Ridge { path: &[(6.0, 45.5), (15.0, 47.0)], width_deg: 1.8, height_m: 2100.0 },
    Ridge { path: &[(8.0, 60.0), (14.0, 65.0), (20.0, 69.0)], width_deg: 2.0, height_m: 1200.0 },
    Ridge { path: &[(60.0, 50.0), (60.0, 68.0)], width_deg: 1.3, height_m: 700.0 },
    Ridge { path: &[(40.0, 43.0), (48.0, 41.0)], width_deg: 1.5, height_m: 2500.0 },
    Ridge { path: &[(46.0, 35.0), (56.0, 27.0)], width_deg: 2.0, height_m: 1900.0 },
    Ridge { path: &[(66.0, 36.0), (76.0, 37.0), (95.0, 28.5)], width_deg: 2.3, height_m: 4000.0 },
    Ridge { path: &[(78.0, 43.0), (90.0, 46.0), (100.0, 50.0)], width_deg: 3.0, height_m: 2400.0 },
    Ridge { path: &[(74.0, 19.0), (77.0, 9.0)], width_deg: 1.0, height_m: 1100.0 },
    Ridge { path: &[(-6.0, 32.0), (8.0, 36.0)], width_deg: 2.0, height_m: 1600.0 },
    Ridge { path: &[(38.0, 11.0), (37.0, 5.0), (35.0, -5.0), (30.0, -12.0)], width_deg: 2.5, height_m: 1700.0 },
    Ridge { path: &[(29.0, -26.0), (28.0, -32.0)], width_deg: 2.0, height_m: 1500.0 },
    Ridge { path: &[(146.0, -17.0), (150.0, -37.0)], width_deg: 1.8, height_m: 900.0 },
    Ridge { path: &[(137.0, 60.0), (145.0, 68.0), (160.0, 66.0)], width_deg: 3.0, height_m: 1000.0 },
];

/// Plateaus and domes: (lon, lat, half-width lon, half-height lat, height m).
const PLATEAUS: &[(f64, f64, f64, f64, f64)] = &[
    (88.0, 33.0, 12.0, 4.5, 4200.0),   // Tibet
    (-110.0, 39.0, 6.0, 4.0, 1300.0),  // Colorado plateau / Great Basin
    (-67.0, -18.0, 3.0, 4.0, 3000.0),  // Altiplano
    (-40.0, 72.0, 12.0, 9.0, 2600.0),  // Greenland ice sheet
    (45.0, 8.0, 5.0, 3.0, 1000.0),     // Horn
    (24.0, -22.0, 9.0, 6.0, 1000.0),   // southern African plateau
];

fn in_polygon(poly: &[(f64, f64)], lon: f64, lat: f64) -> bool {
    let mut inside = false;
    let mut j = poly.len() - 1;
    for i in 0..poly.len() {
        let (xi, yi) = poly[i];
        let (xj, yj) = poly[j];
        if (yi > lat) != (yj > lat) && lon < (xj - xi) * (lat - yi) / (yj - yi) + xi {
            inside = !inside;
        }
        j = i;
    }
    inside
}

fn in_ellipse(e: &(f64, f64, f64, f64), lon: f64, lat: f64) -> bool {
    let dx = ((lon - e.0 + 540.0).rem_euclid(360.0) - 180.0) / e.2;
    let dy = (lat - e.1) / e.3;
    dx * dx + dy * dy <= 1.0
}

fn dist_to_segment(p: (f64, f64), a: (f64, f64), b: (f64, f64)) -> f64 {
    let (abx, aby) = (b.0 - a.0, b.1 - a.1);
    let len2 = (abx * abx + aby * aby).max(1e-9);
    let t = (((p.0 - a.0) * abx + (p.1 - a.1) * aby) / len2).clamp(0.0, 1.0);
    ((p.0 - a.0 - t * abx).powi(2) + (p.1 - a.1 - t * aby).powi(2)).sqrt()
}

/// Is (lon, lat) in degrees land on the crude Earth?
pub fn earth_is_land(lon: f64, lat: f64) -> bool {
    if lat < -71.0 + 3.0 * (lon.to_radians() * 2.0).sin() {
        return true; // Antarctica
    }
    let land = [NORTH_AMERICA, SOUTH_AMERICA, AFRICA, EURASIA, AUSTRALIA, GREENLAND]
        .iter()
        .any(|p| in_polygon(p, lon, lat))
        || ISLANDS.iter().any(|e| in_ellipse(e, lon, lat));
    land && !SEAS.iter().any(|e| in_ellipse(e, lon, lat))
}

/// Elevation in metres above sea level for the crude Earth (negative = ocean).
pub fn earth_height_m(lon: f64, lat: f64) -> f32 {
    if !earth_is_land(lon, lat) {
        return -3000.0;
    }
    let mut h: f64 = 250.0;
    if lat < -71.0 {
        h = 400.0 + 2600.0 * ((-71.0 - lat) / 15.0).clamp(0.0, 1.0);
    }
    let mut ridge_h: f64 = 0.0;
    for r in RIDGES {
        let d = r.path.windows(2).map(|s| dist_to_segment((lon, lat), s[0], s[1])).fold(f64::MAX, f64::min);
        ridge_h = ridge_h.max(r.height_m * (-(d / r.width_deg).powi(2)).exp());
    }
    h += ridge_h;
    for p in PLATEAUS {
        let dx = (lon - p.0) / p.2;
        let dy = (lat - p.1) / p.3;
        h = h.max(p.4 * (-(dx * dx + dy * dy).powi(2)).exp());
    }
    h as f32
}

/// The crude Earth on a `width x height` lat/lon grid.
pub fn earth_like_height(width: u32, height: u32) -> LatLonGrid<f32> {
    LatLonGrid::from_fn(width, height, |lat, lon| earth_height_m(lon.to_degrees(), lat.to_degrees()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn known_places_are_land_or_sea() {
        for &(lon, lat) in &[(-100.0, 40.0), (30.0, 5.0), (-60.0, -10.0), (100.0, 50.0), (135.0, -25.0), (0.0, -80.0)] {
            assert!(earth_is_land(lon, lat), "({lon},{lat}) should be land");
        }
        for &(lon, lat) in &[(-30.0, 30.0), (-150.0, 0.0), (80.0, -30.0), (0.0, 80.0), (-85.0, 60.0), (34.0, 43.0)] {
            assert!(!earth_is_land(lon, lat), "({lon},{lat}) should be ocean");
        }
    }

    #[test]
    fn land_fraction_is_earthlike() {
        let g = earth_like_height(180, 90);
        let (mut land, mut total) = (0.0, 0.0);
        for y in 0..g.height {
            let w = g.lat_of_row(y).cos();
            for x in 0..g.width {
                total += w;
                if *g.get(x, y) >= 0.0 {
                    land += w;
                }
            }
        }
        let frac = land / total;
        assert!((0.25..0.38).contains(&frac), "land fraction {frac}");
    }

    #[test]
    fn mountains_are_high() {
        assert!(earth_height_m(88.0, 33.0) > 3500.0, "Tibet");
        assert!(earth_height_m(-70.0, -25.0) > 2500.0, "Andes");
        assert!(earth_height_m(-50.0, 12.0) < 0.0);
    }
}
