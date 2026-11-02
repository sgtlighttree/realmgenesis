//! Earth reference: the same terrain metrics as `probe_realism`, run on real elevation data,
//! so our targets come from measurement, not guesses:
//! `cargo run -p rg-pipeline --release --example earth_reference -- [half_km]`
//!
//! Data: AWS Terrain Tiles (Mapzen "terrarium" encoding, h = R·256 + G + B/256 − 32768),
//! zoom 11 (~76 m/px at the equator), SRTM-derived on land. Tiles are cached in
//! `target/earth-cache/`. Fetched with `curl` (no new dependency).
//!
//! Notes on comparability: real DEMs carry noise, vegetation and buildings (SRTM is a
//! surface model) and real lakes, so their depression share is a *ceiling* for "natural".
#[path = "common/metrics.rs"]
mod metrics;

use std::collections::HashMap;

const ZOOM: u32 = 11;
const R_EARTH: f64 = 6_371_000.0;

struct Tiles {
    cache: HashMap<(u32, u32), Vec<f64>>,
    dir: std::path::PathBuf,
}

impl Tiles {
    fn tile(&mut self, x: u32, y: u32) -> &Vec<f64> {
        if !self.cache.contains_key(&(x, y)) {
            let path = self.dir.join(format!("{ZOOM}-{x}-{y}.png"));
            if !path.exists() {
                let url = format!("https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{ZOOM}/{x}/{y}.png");
                let ok = std::process::Command::new("curl").args(["-sS", "-f", "-o"]).arg(&path).arg(&url).status().is_ok_and(|s| s.success());
                assert!(ok, "fetching {url} failed");
            }
            let dec = png::Decoder::new(std::io::BufReader::new(std::fs::File::open(&path).unwrap()));
            let mut rd = dec.read_info().unwrap();
            let mut buf = vec![0u8; rd.output_buffer_size().unwrap()];
            let info = rd.next_frame(&mut buf).unwrap();
            let ch = info.color_type.samples();
            let h: Vec<f64> = buf[..info.buffer_size()]
                .chunks(ch)
                .map(|p| p[0] as f64 * 256.0 + p[1] as f64 + p[2] as f64 / 256.0 - 32768.0)
                .collect();
            self.cache.insert((x, y), h);
        }
        &self.cache[&(x, y)]
    }

    /// Height at (lat, lon) degrees, bilinear in Web-Mercator pixels.
    fn height(&mut self, lat: f64, lon: f64) -> f64 {
        let n = (1u64 << ZOOM) as f64 * 256.0;
        let px = (lon + 180.0) / 360.0 * n;
        let lr = lat.to_radians();
        let py = (1.0 - (lr.tan() + 1.0 / lr.cos()).ln() / std::f64::consts::PI) / 2.0 * n;
        let (x0, y0) = ((px - 0.5).floor(), (py - 0.5).floor());
        let (fx, fy) = (px - 0.5 - x0, py - 0.5 - y0);
        let mut at = |x: f64, y: f64| {
            let (xi, yi) = (x as i64, y as i64);
            let (tx, ty) = ((xi.div_euclid(256)) as u32, (yi.div_euclid(256)) as u32);
            let t = self.tile(tx, ty);
            t[(yi.rem_euclid(256) * 256 + xi.rem_euclid(256)) as usize]
        };
        let (a, b, c, d) = (at(x0, y0), at(x0 + 1.0, y0), at(x0, y0 + 1.0), at(x0 + 1.0, y0 + 1.0));
        let top = a + (b - a) * fx;
        let bot = c + (d - c) * fx;
        top + (bot - top) * fy
    }
}

fn main() {
    let half_km: f64 = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(40.0);
    let half = half_km * 1000.0;
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../target/earth-cache");
    std::fs::create_dir_all(&dir).unwrap();
    let tiles = std::cell::RefCell::new(Tiles { cache: HashMap::new(), dir });
    // (name, lat, lon, regime)
    let regions = [
        ("Taiwan Central Range", 23.70, 121.05, "steep fluvial, fast uplift"),
        ("Appalachian Plateau", 38.40, -81.40, "dendritic fluvial, slow uplift"),
        ("Himalaya (Khumbu)", 27.90, 86.90, "extreme relief, glacial + fluvial"),
        ("Alps (Valais)", 46.10, 7.60, "glaciated"),
        ("Great Plains (Kansas)", 38.80, -99.50, "lowland"),
        ("Ozark Plateau", 36.80, -92.50, "dissected low plateau"),
    ];
    println!("region                  depress  depth   coh   asp    β    θ(n)        slope° mean/p50/p90/p99   regime");
    for (name, lat0, lon0, regime) in regions {
        let h = |x: f64, y: f64| {
            let lat = lat0 + (y / R_EARTH).to_degrees();
            let lon = lon0 + (x / (R_EARTH * lat0.to_radians().cos())).to_degrees();
            tiles.borrow_mut().height(lat, lon)
        };
        let n = (2.0 * half / 250.0) as usize;
        let g250 = metrics::grid(&h, n, 250.0);
        let (dp, dd) = metrics::depressions(&g250, n);
        let g60 = metrics::grid(&h, 167, 60.0);
        let (coh, asp) = metrics::orientation(&g60, 167, 60.0);
        let beta = metrics::spectrum_beta(&g60, 167, 60.0);
        let (theta, tn) = metrics::slope_area_theta(&g60, 167, 60.0);
        let sl = metrics::slopes(&h, half);
        println!(
            "{name:<22} {dp:6.2}%  {dd:4.0} m  {coh:.2}  {asp:.2}  {beta:4.2}  {theta:4.2}({tn:5})   {:4.1} {:4.1} {:4.1} {:4.1}   {regime}",
            sl[0], sl[1], sl[2], sl[3]
        );
    }
}
