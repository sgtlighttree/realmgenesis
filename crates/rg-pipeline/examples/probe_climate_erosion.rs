//! Does climate shape erosion? Mean erosion (tectonic height − eroded height, before the
//! hypsometric remap would hide it) on the wettest vs driest third of land, with and without
//! runoff weighting: `cargo run -p rg-pipeline --release --example probe_climate_erosion -- <seed>`
use rg_pipeline::WorldParams;

fn main() {
    let seed: u64 = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(3);
    let mut p = WorldParams::default();
    p.tectonics.seed = seed;
    let tect = rg_tectonics::generate(&p.tectonics);
    let mut ep = p.erosion.clone();
    ep.hypsometry_scale = 0.0; // compare raw physics (p99 rescale only)
    ep.linear_mean_m = 0.0;
    ep.seed = seed;
    let runoff = rg_pipeline::relative_runoff_for_probe(&p.planet, &tect);
    let w = rg_pipeline::generate(&p);
    let wb = rg_pipeline::water_balance(&w.params.planet, &w.climate, &w.tectonics.height_m);
    println!("mean land runoff {:.0} mm/yr (Earth ≈ 300)", wb.mean_land_runoff_mm);
    let dry_wet = |with: bool| {
        let r = rg_erosion::erode_with_runoff(&tect, p.tectonics.radius_m, &ep, with.then_some(&runoff));
        let mut v: Vec<(f32, f64)> = Vec::new();
        for (k, _, _, _, d) in r.height_m.iter_samples() {
            let t = tect.height_m.sample(d) as f64;
            if t > 200.0 {
                v.push((runoff.sample(d), t - r.height_m.data[k] as f64));
            }
        }
        v.sort_by(|a, b| a.0.total_cmp(&b.0));
        let third = v.len() / 3;
        let mean = |s: &[(f32, f64)]| s.iter().map(|x| x.1).sum::<f64>() / s.len().max(1) as f64;
        (mean(&v[..third]), mean(&v[2 * third..]))
    };
    for with in [false, true] {
        let (dry, wet) = dry_wet(with);
        println!("runoff {:5}: mean lowering, driest third {dry:6.0} m, wettest third {wet:6.0} m (wet/dry {:.2})", with, wet / dry.max(1.0));
    }
}
