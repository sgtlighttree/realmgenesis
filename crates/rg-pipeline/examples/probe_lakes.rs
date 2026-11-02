//! Endorheic land and lakes from the water balance (`rg_erosion::lakes`):
//! `cargo run -p rg-pipeline --release --example probe_lakes -- <seed> [surface]`
//! surface: `tectonic` (raw, before erosion) or `final` (the committed macro). Earth: ~20%
//! of ice-free land drains internally, mostly in dry (Köppen B) climates.
use rg_pipeline::{generate, WorldParams};

fn main() {
    let seed: u64 = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(3);
    let surface = std::env::args().nth(2).unwrap_or_else(|| "tectonic".into());
    let mut p = WorldParams::default();
    p.tectonics.seed = seed;
    let w = generate(&p);
    let tect = rg_tectonics::generate(&p.tectonics);
    let field = if surface == "final" { &w.tectonics.height_m } else { &tect.height_m };
    let wb = rg_pipeline::water_balance(&w.params.planet, &w.climate, &w.tectonics.height_m);
    let grid = rg_tectonics::grid::Grid::new(field.res);
    let n = grid.n;
    let r = p.tectonics.radius_m;
    let h: Vec<f64> = (0..n).map(|i| field.data[grid.node_sample[i] as usize] as f64).collect();
    let ocean: Vec<bool> = h.iter().map(|&v| v <= 0.0).collect();
    let area: Vec<f64> = grid.weight.iter().map(|&w| w as f64 * 4.0 * std::f64::consts::PI * r * r).collect();
    let runoff: Vec<f64> = (0..n).map(|i| wb.runoff_mm.sample(grid.dir[i]) as f64 * 1e-3).collect();
    let net_evap: Vec<f64> = (0..n).map(|i| (wb.pet_mm.sample(grid.dir[i]) - wb.precip_mm.sample(grid.dir[i])) as f64 * 1e-3).collect();
    let t0 = std::time::Instant::now();
    let l = rg_erosion::lakes::compute(&grid, &h, &ocean, &runoff, &net_evap, &area);
    let land: f64 = (0..n).filter(|&i| !ocean[i]).map(|i| area[i]).sum();
    let endo: f64 = (0..n).filter(|&i| !ocean[i] && l.endorheic[i]).map(|i| area[i]).sum();
    let lake: f64 = (0..n).filter(|&i| !ocean[i] && !l.level[i].is_nan()).map(|i| area[i]).sum();
    let closed: f64 = (0..n).filter(|&i| l.closed_lake[i]).map(|i| area[i]).sum();
    println!(
        "{surface}: endorheic land {:.1}% (Earth ~20%), lakes {:.2}% of land (closed {:.2}%), {:.0} ms",
        endo / land * 100.0,
        lake / land * 100.0,
        closed / land * 100.0,
        t0.elapsed().as_secs_f64() * 1e3
    );
}
