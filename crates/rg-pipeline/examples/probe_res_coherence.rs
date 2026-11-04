//! How much does the same seed change between resolutions? (Matt, s02/s03: "generating at
//! different resolutions results in very different terrains, even with the same coastlines".)
//!
//! `cargo run -p rg-pipeline --release --example probe_res_coherence -- <seed> <res_a> <res_b>`
//! (other `RG_*` overrides as in `common/mod.rs`, e.g. `RG_JITTER`, `RG_ROUGH_WL`; both
//! worlds use them). Reports, over land cells of the coarser world A:
//! - **height r**: correlation of A's height with B's height averaged over A's cell;
//! - **river overlap**: of A's river cells (drainage ≥ 20,000 km²), the share where B has a
//!   river of at least half that area within one A cell. 1.0 = same rivers in the same places.
#[path = "common/mod.rs"]
mod common;
use glam::DVec3;

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let seed: u64 = a[0].parse().unwrap();
    let (ra, rb): (u32, u32) = (a[1].parse().unwrap(), a[2].parse().unwrap());
    let world_at = |res: u32| {
        // SAFETY: single-threaded at this point (before generation starts its own threads).
        unsafe { std::env::set_var("RG_RES", res.to_string()) };
        common::world(seed)
    };
    let wa = world_at(ra);
    let wb = world_at(rb);
    let (ha, hb) = (&wa.tectonics.height_m, &wb.tectonics.height_m);
    // No hydrology with erosion off (RG_EROSION=0): heights only.
    let drain = wa.hydrology.as_ref().zip(wb.hydrology.as_ref()).map(|(x, y)| (&x.drainage_km2, &y.drainage_km2));
    let cell = ha.spacing_rad();

    // Points spread over A's cell around `d` (a 5×5 pattern in the local tangent plane).
    let around = |d: DVec3| {
        let east = DVec3::Z.cross(d).normalize_or(DVec3::X);
        let north = d.cross(east);
        (0..25).map(move |k| {
            let (x, y) = ((k % 5) as f64 / 4.0 - 0.5, (k / 5) as f64 / 4.0 - 0.5);
            (d + (east * x + north * y) * cell).normalize()
        })
    };
    let (mut xs, mut ys) = (Vec::new(), Vec::new());
    let (mut rivers, mut matched) = (0usize, 0usize);
    // Plates: share of land cells whose plate id is the same at both resolutions, and share
    // of convergent-boundary cells of A that B also has as convergent within one A cell.
    let (pa, pb) = (&wa.tectonics.plate_id, &wb.tectonics.plate_id);
    let (ba, bb) = (&wa.tectonics.boundary, &wb.tectonics.boundary);
    let (mut same_plate, mut conv_a, mut conv_hit) = (0usize, 0usize, 0usize);
    for (k, _, _, _, d) in ha.iter_samples() {
        let h = ha.data[k] as f64;
        if h <= 0.0 || hb.sample(d) <= 0.0 {
            continue;
        }
        let hb_mean = around(d).map(|p| hb.sample_f64(p)).sum::<f64>() / 25.0;
        xs.push(h);
        ys.push(hb_mean);
        if pa.data[k] == pb.nearest(d) {
            same_plate += 1;
        }
        if ba.data[k] == rg_tectonics::boundary_class::CONVERGENT {
            conv_a += 1;
            let hit = around(d).chain(around(d).map(|p| (p + (p - d)).normalize())).any(|p| bb.nearest(p) == rg_tectonics::boundary_class::CONVERGENT);
            if hit {
                conv_hit += 1;
            }
        }
        let Some((da, db)) = drain else { continue };
        let area = da.data[k] as f64;
        if area >= 20_000.0 {
            rivers += 1;
            // Within one A cell: the 5×5 pattern at twice the spread.
            let east = DVec3::Z.cross(d).normalize_or(DVec3::X);
            let north = d.cross(east);
            let best = (0..25)
                .map(|k| {
                    let (x, y) = ((k % 5) as f64 / 2.0 - 1.0, (k / 5) as f64 / 2.0 - 1.0);
                    db.sample_f64((d + (east * x + north * y) * cell).normalize())
                })
                .fold(0.0, f64::max);
            if best >= 0.5 * area {
                matched += 1;
            }
        }
    }
    let r = pearson(&xs, &ys);
    let rms = (xs.iter().zip(&ys).map(|(x, y)| (x - y).powi(2)).sum::<f64>() / xs.len() as f64).sqrt();
    println!(
        "seed {seed} res {ra} vs {rb}: land cells {} · height r {r:.3} · rms {rms:.0} m · river overlap {:.1}% ({matched}/{rivers})",
        xs.len(),
        100.0 * matched as f64 / rivers.max(1) as f64
    );
    println!(
        "  same plate id {:.1}% · convergent boundary cells of A matched in B {:.1}% ({conv_hit}/{conv_a})",
        100.0 * same_plate as f64 / xs.len().max(1) as f64,
        100.0 * conv_hit as f64 / conv_a.max(1) as f64
    );
}

fn pearson(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let (mx, my) = (x.iter().sum::<f64>() / n, y.iter().sum::<f64>() / n);
    let (mut sxy, mut sxx, mut syy) = (0.0, 0.0, 0.0);
    for (a, b) in x.iter().zip(y) {
        sxy += (a - mx) * (b - my);
        sxx += (a - mx).powi(2);
        syy += (b - my).powi(2);
    }
    sxy / (sxx * syy).sqrt().max(1e-12)
}
