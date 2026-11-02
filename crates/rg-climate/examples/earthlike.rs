//! Run Tier A on the crude synthetic Earth and write PNGs to look at.
//!
//! `cargo run -p rg-climate --release --example earthlike [out_dir]`
//!
//! Writes `koppen.png`, `temp_jan.png`, `temp_jul.png`, `precip_annual.png`, plus
//! `precip_jan.png`, `precip_jul.png` and `height.png`.

use rg_climate::synthetic::earth_like_height;
use rg_climate::{ClimateFields, ClimateInput, Koppen, LatLonGrid, PlanetParams, tier_a};
use std::path::{Path, PathBuf};
use std::time::Instant;

const DEFAULT_OUT: &str = "/tmp/claude-0/-home-user-realmgenesis/a86a78b6-6f71-50e0-ac50-4fba2f8d0ce0/scratchpad/climate";

fn main() {
    let out_dir = PathBuf::from(std::env::args().nth(1).unwrap_or_else(|| DEFAULT_OUT.to_string()));
    std::fs::create_dir_all(&out_dir).expect("create output directory");

    let planet = PlanetParams::default();
    let height = earth_like_height(360, 180);
    let input = ClimateInput { height_m: height.clone() };

    let start = Instant::now();
    let fields = tier_a(&planet, &input);
    println!("tier_a at 360x180: {:.1} ms", start.elapsed().as_secs_f64() * 1000.0);

    write_koppen(&out_dir.join("koppen.png"), &fields, &height);
    write_scalar(&out_dir.join("temp_jan.png"), &fields.temp_grid(0), &height, temp_color);
    write_scalar(&out_dir.join("temp_jul.png"), &fields.temp_grid(6), &height, temp_color);
    write_scalar(&out_dir.join("precip_annual.png"), &fields.annual_precip_grid(), &height, |p| rain_color(p / 12.0));
    write_scalar(&out_dir.join("precip_jan.png"), &fields.precip_grid(0), &height, rain_color);
    write_scalar(&out_dir.join("precip_jul.png"), &fields.precip_grid(6), &height, rain_color);
    write_scalar(&out_dir.join("height.png"), &height, &height, |h| {
        let v = (h / 4000.0).clamp(0.0, 1.0);
        let g = (255.0 * (1.0 - v * 0.8)) as u8;
        [g, g, g]
    });
    print_stats(&fields);
    print_zonal(&fields, &height);
    for (n, lon, lat) in [("Delhi", 77.2, 28.6), ("Atlanta", -84.0, 33.7), ("Beijing", 116.4, 39.9), ("Sahel", 10.0, 14.0), ("NYC", -74.0, 40.7), ("Houston", -95.0, 30.0)] {
        println!("wind {n:8} {}", wind_at(&fields, lon, lat));
    }
    println!("wrote PNGs to {}", out_dir.display());
}

/// Blue ... white ... red ramp, -40..40 C.
fn temp_color(t: f32) -> [u8; 3] {
    let stops: [(f32, [f32; 3]); 7] = [
        (-40.0, [90.0, 0.0, 140.0]),
        (-25.0, [40.0, 60.0, 200.0]),
        (-10.0, [80.0, 170.0, 240.0]),
        (5.0, [200.0, 235.0, 200.0]),
        (18.0, [255.0, 240.0, 120.0]),
        (28.0, [240.0, 130.0, 40.0]),
        (40.0, [160.0, 10.0, 10.0]),
    ];
    ramp(t, &stops)
}

/// White to deep blue for mm per month (0..400).
fn rain_color(mm_month: f32) -> [u8; 3] {
    let stops: [(f32, [f32; 3]); 6] = [
        (0.0, [245.0, 225.0, 170.0]),
        (20.0, [250.0, 250.0, 200.0]),
        (60.0, [150.0, 220.0, 150.0]),
        (120.0, [60.0, 160.0, 220.0]),
        (250.0, [20.0, 60.0, 180.0]),
        (450.0, [60.0, 0.0, 100.0]),
    ];
    ramp(mm_month, &stops)
}

fn ramp(v: f32, stops: &[(f32, [f32; 3])]) -> [u8; 3] {
    if v <= stops[0].0 {
        return stops[0].1.map(|c| c as u8);
    }
    for pair in stops.windows(2) {
        if v <= pair[1].0 {
            let t = (v - pair[0].0) / (pair[1].0 - pair[0].0);
            return std::array::from_fn(|k| (pair[0].1[k] + (pair[1].1[k] - pair[0].1[k]) * t) as u8);
        }
    }
    stops[stops.len() - 1].1.map(|c| c as u8)
}

fn write_scalar(path: &Path, grid: &LatLonGrid<f32>, height: &LatLonGrid<f32>, color: impl Fn(f32) -> [u8; 3]) {
    let mut rgb = Vec::with_capacity(grid.data.len() * 3);
    for (v, h) in grid.data.iter().zip(&height.data) {
        let c = color(*v);
        // Darken ocean slightly so coastlines show on temperature maps.
        let c = if *h < 0.0 { c.map(|x| (x as f32 * 0.85) as u8) } else { c };
        rgb.extend_from_slice(&c);
    }
    write_png(path, grid.width, grid.height, &rgb);
}

fn write_koppen(path: &Path, fields: &ClimateFields, height: &LatLonGrid<f32>) {
    let mut rgb = Vec::with_capacity(fields.koppen.data.len() * 3);
    for (k, h) in fields.koppen.data.iter().zip(&height.data) {
        match Koppen::from_u8(*k) {
            Some(c) => rgb.extend_from_slice(&c.color()),
            None => {
                let _ = h;
                rgb.extend_from_slice(&[225, 235, 245]);
            }
        }
    }
    write_png(path, fields.koppen.width, fields.koppen.height, &rgb);
}

fn write_png(path: &Path, w: u32, h: u32, rgb: &[u8]) {
    let file = std::fs::File::create(path).expect("create png");
    let mut enc = png::Encoder::new(std::io::BufWriter::new(file), w * 2, h * 2);
    enc.set_color(png::ColorType::Rgb);
    enc.set_depth(png::BitDepth::Eight);
    let mut writer = enc.write_header().expect("png header");
    // Upscale 2x (nearest) so the images are easy to read.
    let mut big = Vec::with_capacity(rgb.len() * 4);
    for y in 0..h as usize {
        let row: Vec<u8> = (0..w as usize).flat_map(|x| {
            let p = &rgb[(y * w as usize + x) * 3..][..3];
            [p[0], p[1], p[2], p[0], p[1], p[2]]
        }).collect();
        big.extend_from_slice(&row);
        big.extend_from_slice(&row);
    }
    writer.write_image_data(&big).expect("png data");
}

/// Reference stations: name, lon, lat, real annual mean temperature (C), real annual
/// precipitation (mm) and real Köppen class (rounded values from national normals / Wikipedia).
/// The synthetic Earth is crude, so only broad agreement is expected.
const STATIONS: &[(&str, f64, f64, f32, f32, &str)] = &[
    ("Manaus", -60.0, -3.1, 27.0, 2300.0, "Af"),
    ("Kinshasa", 15.3, -4.3, 25.0, 1400.0, "Aw"),
    ("Khartoum", 32.5, 15.6, 30.0, 160.0, "BWh"),
    ("Cairo", 31.2, 30.0, 22.0, 25.0, "BWh"),
    ("Timbuktu", -3.0, 16.8, 29.0, 200.0, "BWh"),
    ("CapeTown", 18.4, -33.9, 17.0, 515.0, "Csb"),
    ("Joburg", 28.0, -26.0, 16.0, 700.0, "Cwb"),
    ("Perth", 115.9, -32.0, 18.5, 730.0, "Csa"),
    ("Sydney", 151.2, -33.9, 18.0, 1200.0, "Cfa"),
    ("BuenosAires", -58.4, -34.6, 18.0, 1100.0, "Cfa"),
    ("Valparaiso", -71.6, -33.0, 14.5, 370.0, "Csb"),
    ("Lima", -77.0, -12.0, 19.0, 15.0, "BWh"),
    ("Madrid", -3.7, 40.4, 14.0, 420.0, "BSk"),
    ("Rome", 12.5, 41.9, 15.5, 800.0, "Csa"),
    ("London", -0.1, 51.5, 11.0, 600.0, "Cfb"),
    ("Paris", 2.3, 48.9, 11.0, 640.0, "Cfb"),
    ("Moscow", 37.6, 55.8, 5.5, 700.0, "Dfb"),
    ("Yakutsk", 129.7, 62.0, -9.0, 240.0, "Dwd"),
    ("Delhi", 77.2, 28.6, 25.0, 790.0, "Cwa"),
    ("Mumbai", 72.8, 19.0, 27.0, 2100.0, "Aw"),
    ("Bangkok", 100.5, 13.7, 28.0, 1500.0, "Aw"),
    ("Beijing", 116.4, 39.9, 12.9, 580.0, "Dwa"),
    ("Shanghai", 121.5, 31.2, 16.0, 1200.0, "Cfa"),
    ("Tokyo", 139.7, 35.7, 16.0, 1500.0, "Cfa"),
    ("NewYork", -74.0, 40.7, 13.0, 1200.0, "Cfa"),
    ("Chicago", -87.6, 41.9, 10.0, 950.0, "Dfa"),
    ("Atlanta", -84.4, 33.7, 17.0, 1300.0, "Cfa"),
    ("LosAngeles", -118.2, 34.0, 18.0, 370.0, "BSh"),
    ("Seattle", -122.3, 47.6, 11.0, 950.0, "Csb"),
    ("Denver", -105.0, 39.7, 10.0, 400.0, "BSk"),
    ("Edmonton", -113.5, 53.5, 4.0, 450.0, "Dfb"),
    ("Anchorage", -150.0, 61.2, 2.5, 410.0, "Dfc"),
    ("MexicoCity", -99.0, 19.4, 16.5, 700.0, "Cwb"),
    ("Phoenix", -112.0, 33.4, 24.0, 200.0, "BWh"),
    ("Tehran", 51.4, 35.7, 17.0, 230.0, "BSk"),
    ("Riyadh", 46.7, 24.7, 26.0, 100.0, "BWh"),
    ("Ulaanbaatar", 106.9, 47.9, -1.0, 260.0, "BSk"),
];

/// Print model vs real values for the reference stations, and the Köppen first-letter match.
fn print_stats(f: &ClimateFields) {
    let (mut letter_ok, mut exact_ok, mut counted) = (0, 0, 0);
    println!("station       model  real   MAT  (real)   MAP  (real)  Tjan Tjul  Pjan Pjul");
    for &(name, lon, lat, real_t, real_p, real_k) in STATIONS {
        let x = (((lon + 180.0) / 360.0 * f.width as f64) as u32).min(f.width - 1);
        let y = (((90.0 - lat) / 180.0 * f.height as f64) as u32).min(f.height - 1);
        let i = f.index(x, y);
        let (t, p) = (&f.temp_c[i], &f.precip_mm[i]);
        let k = Koppen::from_u8(*f.koppen.get(x, y)).map_or("ocean", |k| k.code());
        counted += 1;
        letter_ok += (k.chars().next() == real_k.chars().next()) as u32;
        exact_ok += (k == real_k) as u32;
        println!(
            "{name:12} {k:5} {real_k:5} {:5.1} ({real_t:5.1}) {:5.0} ({real_p:5.0}) {:5.1} {:5.1} {:5.0} {:5.0}",
            t.iter().sum::<f32>() / 12.0, p.iter().sum::<f32>(), t[0], t[6], p[0], p[6]
        );
    }
    println!("Koppen first letter matches {letter_ok}/{counted}, exact {exact_ok}/{counted}");
}

/// Zonal means over land cells, per 10 degree band, to compare with textbook values
/// (land annual precipitation roughly: 0-10N 1700, 10-20N 1000, 20-30N 400, 30-40N 550,
/// 40-50N 600, 50-60N 550, 60-70N 350 mm).
fn print_zonal(f: &ClimateFields, height: &LatLonGrid<f32>) {
    println!("band     landcells  MAT   Tjan  Tjul  Pann");
    for b in (-8..8).rev() {
        let (lo, hi) = (b * 10, b * 10 + 10);
        let (mut n, mut mat, mut tj, mut tl, mut p) = (0.0f32, 0.0, 0.0, 0.0, 0.0);
        for y in 0..f.height {
            let lat = 90.0 - (y as f32 + 0.5) * 180.0 / f.height as f32;
            if lat < lo as f32 || lat >= hi as f32 { continue; }
            for x in 0..f.width {
                if *height.get(x, y) < 0.0 { continue; }
                let i = f.index(x, y);
                n += 1.0;
                mat += f.temp_c[i].iter().sum::<f32>() / 12.0;
                tj += f.temp_c[i][0];
                tl += f.temp_c[i][6];
                p += f.precip_mm[i].iter().sum::<f32>();
            }
        }
        if n > 0.0 {
            println!("{lo:4}..{hi:3} {n:8.0} {:6.1} {:5.1} {:5.1} {:6.0}", mat / n, tj / n, tl / n, p / n);
        }
    }
}

#[allow(dead_code)]
fn wind_at(f: &ClimateFields, lon: f64, lat: f64) -> String {
    let x = (((lon + 180.0) / 360.0 * f.width as f64) as u32).min(f.width - 1);
    let y = (((90.0 - lat) / 180.0 * f.height as f64) as u32).min(f.height - 1);
    let i = f.index(x, y);
    format!("u/v Jan {:.1}/{:.1} Jul {:.1}/{:.1}", f.wind_u[i][0], f.wind_v[i][0], f.wind_u[i][6], f.wind_v[i][6])
}
