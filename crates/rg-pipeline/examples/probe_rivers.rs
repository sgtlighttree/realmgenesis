//! Numbers on the drawn river network (s08b: "too wide, too many" on imported maps):
//! `cargo run -p rg-pipeline --release --example probe_rivers -- <seed>`
//! Same `RG_*` env overrides as the other probes (`RG_MASK…` for an import).
//!
//! Prints the density (segments and length per 10⁶ km² of land), the width the pre-s08b
//! drawing gave them on the ground, how much of the length runs on flat ground (a drop under
//! 1 m per segment: filled depressions and drawn plains) or beside the coast, and how many
//! rivers are short stubs.
#[path = "common/mod.rs"]
mod common;
use glam::DVec3;
use std::collections::HashMap;

fn key(d: DVec3) -> (i64, i64, i64) {
    ((d.x * 1e6).round() as i64, (d.y * 1e6).round() as i64, (d.z * 1e6).round() as i64)
}

fn main() {
    let seed: u64 = std::env::args().nth(1).and_then(|a| a.parse().ok()).unwrap_or(1);
    let w = common::world(seed);
    let r_km = w.params.tectonics.radius_m / 1000.0;
    let hm = &w.tectonics.height_m;
    let hy = w.hydrology.as_ref().expect("erosion on");
    let land_frac = hm.data.iter().filter(|&&h| h > 0.0).count() as f64 / hm.data.len() as f64;
    let land_mkm2 = land_frac * 4.0 * std::f64::consts::PI * r_km * r_km / 1.0e6;
    let cell_km = hy.drainage_km2.spacing_rad() * r_km;
    let rivers = &hy.rivers;
    let len_km = |s: &rg_erosion::RiverSegment| s.from.angle_between(s.to) * r_km;
    let total: f64 = rivers.iter().map(len_km).sum();

    // Width on the ground as rivers were drawn before s08b (px of a 2048-wide world image);
    // since then the renderer draws them about 1.2 px wide on screen at any zoom.
    let km_per_px = 2.0 * std::f64::consts::PI * r_km / 2048.0;
    let width_km = |a: f32| ((a as f64 / 2.0e4).sqrt() * 0.45).clamp(0.5, 2.6) * km_per_px;
    let mut widths: Vec<f64> = rivers.iter().map(|s| width_km(s.area_km2)).collect();
    widths.sort_by(f64::total_cmp);
    let pct = |v: &[f64], p: f64| if v.is_empty() { 0.0 } else { v[((v.len() - 1) as f64 * p) as usize] };

    // Flat: the routed surface drops under 1 m along the segment. Coast: sea within 1.5 cells
    // of the segment's upstream end (but not its own mouth).
    let (mut flat, mut coastal) = (0.0, 0.0);
    for s in rivers {
        let l = len_km(s);
        if hm.sample_drained(s.from) - hm.sample_drained(s.to) < 1.0 {
            flat += l;
        }
        let mouth = hm.sample(s.to) <= 0.0;
        let east = DVec3::Z.cross(s.from).normalize_or_zero();
        let north = s.from.cross(east);
        let ring = (0..8).any(|k| {
            let a = k as f64 * std::f64::consts::FRAC_PI_4;
            let d = (s.from + (east * a.cos() + north * a.sin()) * (1.5 * cell_km / r_km)).normalize();
            hm.sample(d) <= 0.0
        });
        if ring && !mouth {
            coastal += l;
        }
    }

    // Rivers as chains: a source is a segment nothing flows into. Its course is the run of
    // segments down to the sea or to a bigger river (the first node where another segment
    // with a larger area upstream joins is not tracked here: the course is to the mouth).
    let by_from: HashMap<_, usize> = rivers.iter().enumerate().map(|(i, s)| (key(s.from), i)).collect();
    let mut inflow = vec![0u32; rivers.len()];
    for s in rivers {
        if let Some(&j) = by_from.get(&key(s.to)) {
            inflow[j] += 1;
        }
    }
    let sources: Vec<usize> = (0..rivers.len()).filter(|&i| inflow[i] == 0).collect();
    // Tributary length: from a source down to the first confluence (a node with 2+ inflows).
    let mut trib: Vec<f64> = sources
        .iter()
        .map(|&i| {
            let (mut k, mut n) = (i, 1.0);
            while let Some(&j) = by_from.get(&key(rivers[k].to)) {
                if inflow[j] > 1 {
                    break;
                }
                k = j;
                n += 1.0;
            }
            n
        })
        .collect();
    trib.sort_by(f64::total_cmp);
    let stubs = trib.iter().filter(|&&n| n <= 2.0).count();

    // s09: where rivers end, how long they run, and how much land is lake.
    let lake = &hy.lake_depth_m;
    let lake_frac = lake.iter_samples().filter(|&(k, _, _, _, d)| lake.data[k] > 2.0 && hm.sample(d) > 0.0).count() as f64 / lake.data.len() as f64;
    let (mut end_sea, mut end_lake, mut end_dry) = (0, 0, 0);
    let (mut len_in_lake, mut len_desc) = (0.0, vec![0.0f64; rivers.len()]);
    for (i, s) in rivers.iter().enumerate() {
        if lake.sample(s.from) > 2.0 {
            len_in_lake += len_km(s);
        }
        if !by_from.contains_key(&key(s.to)) {
            if hm.sample(s.to) <= 0.0 {
                end_sea += 1;
            } else if lake.sample(s.to) > 2.0 || lake.sample(s.from) > 2.0 {
                end_lake += 1;
            } else {
                end_dry += 1;
            }
        }
        let _ = i;
    }
    // Course length from each source to where its water leaves the drawn network.
    let mut courses: Vec<f64> = sources
        .iter()
        .map(|&i| {
            let (mut k, mut l) = (i, len_km(&rivers[i]));
            while let Some(&j) = by_from.get(&key(rivers[k].to)) {
                k = j;
                l += len_km(&rivers[k]);
            }
            l
        })
        .collect();
    courses.sort_by(f64::total_cmp);
    len_desc.clear();
    println!("lakes: {:.2}% of the globe's nodes are lake on land ({:.1}% of land); river length under lakes {:.0}%", lake_frac * 100.0, lake_frac / land_frac.max(1e-9) * 100.0, 100.0 * len_in_lake / total.max(1.0));
    println!("river ends: sea {end_sea}, lake {end_lake}, on dry land {end_dry}");
    println!("source-to-end course km: median {:.0}, p90 {:.0}, max {:.0}; courses ≥ 1000 km: {}", pct(&courses, 0.5), pct(&courses, 0.9), pct(&courses, 1.0), courses.iter().filter(|&&l| l >= 1000.0).count());
    println!("seed {seed}: land {:.1}% = {land_mkm2:.1} Mkm², cell {cell_km:.1} km, threshold {:.0} km²", land_frac * 100.0, w.params.erosion.river_min_area_km2);
    println!("segments {} ({:.0} per Mkm² land), length {total:.0} km ({:.0} km per Mkm² land)", rivers.len(), rivers.len() as f64 / land_mkm2, total / land_mkm2);
    println!("pre-s08b drawn width on the ground: min {:.1} km, median {:.1} km, p90 {:.1} km, max {:.1} km", pct(&widths, 0.0), pct(&widths, 0.5), pct(&widths, 0.9), pct(&widths, 1.0));
    println!("length on flat ground (< 1 m drop per segment): {:.0}%", 100.0 * flat / total.max(1.0));
    println!("length beside the coast (sea within 1.5 cells, not the mouth): {:.0}%", 100.0 * coastal / total.max(1.0));
    println!("rivers (sources) {}, of which stubs (≤ 2 segments to a confluence) {} = {:.0}%; median course {:.0} segments", sources.len(), stubs, 100.0 * stubs as f64 / sources.len().max(1) as f64, pct(&trib, 0.5));
}
