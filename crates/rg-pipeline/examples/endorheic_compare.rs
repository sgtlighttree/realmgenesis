//! Topographic endorheic share, Earth vs ours, measured the same way at the same coarse
//! resolution (0.25° lat-lon, ~28 km, close to the Normal macro cell):
//! `cargo run -p rg-pipeline --release --example endorheic_compare -- <seed>`
//!
//! Rule (identical for both): sea = h ≤ 0 connected to the open ocean; Priority-Flood fill
//! from the sea; a *sink* is a connected depression (fill − h > 0) at least `MIN_DEPTH_M`
//! deep; endorheic land = land whose steepest descent on the raw surface ends in a sink.
//! Earth: AWS Terrain Tiles (terrarium) at zoom 3 (SRTM-derived on land), ±85° latitude.
//! Climate isn't used: this checks whether our *topography* holds Earth-like closures.
#[path = "common/mod.rs"]
mod common;
use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap};

const W: usize = 1440;
const H: usize = 680; // ±85°
const MIN_DEPTH_M: f64 = 50.0;

fn lat_of(y: usize) -> f64 {
    85.0 - (y as f64 + 0.5) * 170.0 / H as f64
}
fn lon_of(x: usize) -> f64 {
    -180.0 + (x as f64 + 0.5) * 360.0 / W as f64
}

fn earth_grid() -> Vec<f64> {
    let zoom = 3u32;
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../target/earth-cache");
    std::fs::create_dir_all(&dir).unwrap();
    let mut tiles: HashMap<(u32, u32), Vec<f64>> = HashMap::new();
    let mut tile = |x: u32, y: u32| -> Vec<f64> {
        tiles
            .entry((x, y))
            .or_insert_with(|| {
                let path = dir.join(format!("{zoom}-{x}-{y}.png"));
                if !path.exists() {
                    let url = format!("https://s3.amazonaws.com/elevation-tiles-prod/terrarium/{zoom}/{x}/{y}.png");
                    assert!(std::process::Command::new("curl").args(["-sS", "-f", "-o"]).arg(&path).arg(&url).status().is_ok_and(|s| s.success()), "fetch {url}");
                }
                let dec = png::Decoder::new(std::io::BufReader::new(std::fs::File::open(&path).unwrap()));
                let mut rd = dec.read_info().unwrap();
                let mut buf = vec![0u8; rd.output_buffer_size().unwrap()];
                let info = rd.next_frame(&mut buf).unwrap();
                let ch = info.color_type.samples();
                buf[..info.buffer_size()].chunks(ch).map(|p| p[0] as f64 * 256.0 + p[1] as f64 + p[2] as f64 / 256.0 - 32768.0).collect()
            })
            .clone()
    };
    let n = (1u64 << zoom) as f64 * 256.0;
    (0..W * H)
        .map(|k| {
            let (lat, lon) = (lat_of(k / W), lon_of(k % W));
            let px = ((lon + 180.0) / 360.0 * n).clamp(0.0, n - 1.0);
            let lr = lat.to_radians();
            let py = ((1.0 - (lr.tan() + 1.0 / lr.cos()).ln() / std::f64::consts::PI) / 2.0 * n).clamp(0.0, n - 1.0);
            let (xi, yi) = (px as u32, py as u32);
            let t = tile(xi / 256, yi / 256);
            t[((yi % 256) * 256 + xi % 256) as usize]
        })
        .collect()
}

#[derive(PartialEq)]
struct It(f64, usize);
impl Eq for It {}
impl Ord for It {
    fn cmp(&self, o: &Self) -> Ordering {
        o.0.total_cmp(&self.0).then(o.1.cmp(&self.1))
    }
}
impl PartialOrd for It {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}

fn neighbours(k: usize) -> impl Iterator<Item = usize> {
    let (x, y) = ((k % W) as i64, (k / W) as i64);
    (-1..=1i64).flat_map(move |dy| (-1..=1i64).map(move |dx| (dx, dy))).filter_map(move |(dx, dy)| {
        if dx == 0 && dy == 0 {
            return None;
        }
        let ny = y + dy;
        if ny < 0 || ny >= H as i64 {
            return None;
        }
        let nx = (x + dx).rem_euclid(W as i64);
        Some(ny as usize * W + nx as usize)
    })
}

/// (endorheic share of land %, sink count, land share in sinks themselves %).
fn endorheic(h: &[f64]) -> (f64, usize, f64) {
    let n = W * H;
    let area = |k: usize| lat_of(k / W).to_radians().cos();
    // Open ocean: flood from the deepest node through h ≤ 0.
    let seed = (0..n).min_by(|&a, &b| h[a].total_cmp(&h[b])).unwrap();
    let mut sea = vec![false; n];
    let mut stack = vec![seed];
    sea[seed] = true;
    while let Some(k) = stack.pop() {
        for j in neighbours(k) {
            if !sea[j] && h[j] <= 0.0 {
                sea[j] = true;
                stack.push(j);
            }
        }
    }
    // Priority-Flood fill from the sea.
    let mut f = h.to_vec();
    let mut done = sea.clone();
    let mut heap: BinaryHeap<It> = (0..n).filter(|&k| sea[k]).map(|k| It(h[k].max(0.0), k)).collect();
    while let Some(It(hc, c)) = heap.pop() {
        for j in neighbours(c) {
            if !done[j] {
                done[j] = true;
                f[j] = f[j].max(hc);
                heap.push(It(f[j], j));
            }
        }
    }
    // Depression components and their max depth.
    let mut comp = vec![u32::MAX; n];
    let mut deep: Vec<bool> = Vec::new();
    for k in 0..n {
        if sea[k] || comp[k] != u32::MAX || f[k] - h[k] <= 0.0 {
            continue;
        }
        let id = deep.len() as u32;
        let mut maxd = 0.0f64;
        let mut st = vec![k];
        comp[k] = id;
        while let Some(c) = st.pop() {
            maxd = maxd.max(f[c] - h[c]);
            for j in neighbours(c) {
                if !sea[j] && comp[j] == u32::MAX && f[j] - h[j] > 0.0 {
                    comp[j] = id;
                    st.push(j);
                }
            }
        }
        deep.push(maxd >= MIN_DEPTH_M);
    }
    // Steepest descent on raw h to its end: inside a deep sink → endorheic.
    let rec: Vec<usize> = (0..n)
        .map(|k| {
            let mut best = (h[k], k);
            for j in neighbours(k) {
                if h[j] < best.0 {
                    best = (h[j], j);
                }
            }
            best.1
        })
        .collect();
    let (mut land, mut endo, mut in_sink) = (0.0, 0.0, 0.0);
    for k in 0..n {
        if sea[k] || h[k] <= -12000.0 {
            continue;
        }
        land += area(k);
        let mut c = k;
        for _ in 0..4000 {
            if rec[c] == c || sea[c] {
                break;
            }
            c = rec[c];
        }
        if comp[c] != u32::MAX && deep[comp[c] as usize] {
            endo += area(k);
        }
        if comp[k] != u32::MAX && deep[comp[k] as usize] {
            in_sink += area(k);
        }
    }
    (endo / land * 100.0, deep.iter().filter(|&&d| d).count(), in_sink / land * 100.0)
}

fn main() {
    let seed: u64 = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(3);
    let earth = earth_grid();
    let (e, es, ei) = endorheic(&earth);
    println!("Earth (0.25°, ±85°): endorheic {e:.1}% of land, {es} sinks ≥ {MIN_DEPTH_M} m, {ei:.1}% of land inside sinks");
    let w = common::world(seed);
    let hm = &w.tectonics.height_m;
    let lake = w.hydrology.as_ref().map(|h| &h.lake_depth_m);
    // Our committed surface: the bed (lakes are depressions, as on Earth's lake-bed DEMs the
    // water surface is what SRTM sees; use the water surface for a like-for-like test).
    let ours: Vec<f64> = (0..W * H)
        .map(|k| {
            let (lat, lon) = (lat_of(k / W).to_radians(), lon_of(k % W).to_radians());
            let d = glam::DVec3::new(lat.cos() * lon.cos(), lat.cos() * lon.sin(), lat.sin());
            // Drained (triangle-split) interpolant: bilinear resampling dams diagonal valleys and
            // would invent sinks. Water surface over lakes (what a DEM sees).
            hm.sample_tri(d) + lake.map_or(0.0, |l| l.sample(d).max(0.0) as f64)
        })
        .collect();
    let (o, os, oi) = endorheic(&ours);
    println!("ours  (seed {seed}, same rule): endorheic {o:.1}% of land, {os} sinks ≥ {MIN_DEPTH_M} m, {oi:.1}% of land inside sinks");
}
