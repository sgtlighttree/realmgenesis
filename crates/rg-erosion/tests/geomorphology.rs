//! Generated landscapes must obey the empirical laws of real ones, and the solver must
//! match the analytical steady state. Uses a real tectonic world (seed 3) at a small res.

use rg_erosion::{erode, ErosionParams};

fn run(seed: u64) -> rg_erosion::ErosionResult {
    let tp = rg_tectonics::TectonicsParams { seed, res: 128, ..Default::default() };
    let snap = rg_tectonics::generate(&tp);
    erode(&snap, tp.radius_m, &ErosionParams::default())
}

#[test]
fn obeys_geomorphic_laws() {
    for seed in [3, 11] {
        let s = run(seed).stats;
        // Hard invariant: all rain reaches the sea.
        assert_eq!(s.drained_fraction, 1.0, "seed {seed}: undrained land");
        // Solver check: with uplift factored out, steady-state concavity equals m.
        assert!((s.concavity_normalised - 0.45).abs() < 0.03, "seed {seed}: θ_norm {}", s.concavity_normalised);
        // Hack's law exponent in the natural range 0.45–0.7.
        assert!((0.45..=0.7).contains(&s.hack_h), "seed {seed}: Hack h {}", s.hack_h);
        // Apparent concavity within what real rivers show (Lague 2014: ~0.3–0.7+).
        assert!((0.3..=0.8).contains(&s.concavity), "seed {seed}: θ {}", s.concavity);
        // Per-basin hypsometry roughly natural (coarse test grid, so a little wider).
        assert!((0.2..=0.6).contains(&s.basin_hypsometric_median), "seed {seed}: HI {}", s.basin_hypsometric_median);
    }
}

#[test]
fn deterministic_and_keeps_the_coast() {
    let tp = rg_tectonics::TectonicsParams { seed: 5, res: 96, ..Default::default() };
    let snap = rg_tectonics::generate(&tp);
    let a = erode(&snap, tp.radius_m, &ErosionParams::default());
    let b = erode(&snap, tp.radius_m, &ErosionParams::default());
    assert!(a.height_m.data == b.height_m.data, "not deterministic");
    for (k, (&t, &e)) in snap.height_m.data.iter().zip(&a.height_m.data).enumerate() {
        assert_eq!(t > 0.0, e > 0.0, "coast moved at sample {k}: {t} → {e}");
    }
}

#[test]
fn hypsometry_is_earth_like() {
    let tp = rg_tectonics::TectonicsParams { seed: 3, res: 128, ..Default::default() };
    let snap = rg_tectonics::generate(&tp);
    let r = erode(&snap, tp.radius_m, &ErosionParams::default());
    let mut land: Vec<f32> = r.height_m.data.iter().copied().filter(|&h| h > 0.0).collect();
    land.sort_by(f32::total_cmp);
    let mean = land.iter().map(|&h| h as f64).sum::<f64>() / land.len() as f64;
    let median = land[land.len() / 2] as f64;
    // Earth: mean ≈ 840 m, median ≈ 440 m (ETOPO-class statistics).
    assert!((700.0..1000.0).contains(&mean), "mean land {mean}");
    assert!((330.0..560.0).contains(&median), "median land {median}");
}
