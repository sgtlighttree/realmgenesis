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
}
