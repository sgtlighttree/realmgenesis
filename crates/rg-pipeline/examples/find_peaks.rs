//! Print the highest points of a world: `cargo run -p rg-pipeline --release --example find_peaks -- <seed>`.
use rg_pipeline::{generate, WorldParams};

fn main() {
    let seed = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(3);
    let mut p = WorldParams::default();
    p.tectonics.seed = seed;
    let w = generate(&p);
    let f = &w.tectonics.height_m;
    let mut pts: Vec<(f32, f64, f64)> = f
        .iter_samples()
        .map(|(k, _, _, _, d)| (f.data[k], d.z.asin().to_degrees(), d.y.atan2(d.x).to_degrees()))
        .collect();
    pts.sort_by(|a, b| b.0.total_cmp(&a.0));
    for (h, lat, lon) in pts.iter().step_by(400).take(8) {
        println!("{h:7.0} m  lat {lat:7.2}  lon {lon:7.2}");
    }
    // Rugged spots: largest height difference to a neighbouring sample (steep macro relief).
    let res = f.res;
    let mut rugged: Vec<(f32, f64, f64)> = Vec::new();
    for (k, face, i, j, d) in f.iter_samples() {
        if i == 0 || j == 0 || i == res || j == res || f.data[k] <= 500.0 {
            continue;
        }
        let nb = [f.index(face, i + 1, j), f.index(face, i - 1, j), f.index(face, i, j + 1), f.index(face, i, j - 1)];
        let relief = nb.iter().map(|&n| (f.data[n] - f.data[k]).abs()).fold(0.0f32, f32::max);
        rugged.push((relief, d.z.asin().to_degrees(), d.y.atan2(d.x).to_degrees()));
    }
    rugged.sort_by(|a, b| b.0.total_cmp(&a.0));
    println!("-- rugged (max neighbour step, m):");
    for (r, lat, lon) in rugged.iter().step_by(300).take(6) {
        println!("{r:7.0} m  lat {lat:7.2}  lon {lon:7.2}");
    }
}
