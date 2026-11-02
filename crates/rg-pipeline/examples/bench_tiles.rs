//! Tile generation cost on a real world, network detail vs the old gully path:
//! `cargo run -p rg-pipeline --release --example bench_tiles -- <seed>` (same cache as the probes).
#[path = "common/mod.rs"]
mod common;
use rg_geo::{CubeFace, TileId};
use rg_terrain::TileGenerator;

fn main() {
    let seed: u64 = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(3);
    let w = common::world(seed);
    let net = rg_pipeline::terrain_sampler(&w);
    let mut old = net.clone();
    old.net = None;
    // A mountain spot: the highest macro sample.
    let f = &w.tectonics.height_m;
    let (k, _) = f.data.iter().enumerate().fold((0, f32::MIN), |b, (i, &v)| if v > b.1 { (i, v) } else { b });
    let d = f.iter_samples().nth(k).unwrap().4;
    let (face, u, v) = CubeFace::from_dir(d);
    // Phase split for one deep tile: patch build vs per-vertex evaluation.
    if let Some(dn) = &net.net {
        for (cutoff, half) in [(1000.0, 16_000.0), (40.0, 640.0)] {
            let t0 = std::time::Instant::now();
            let p = dn.patch(d, half, cutoff);
            let t1 = t0.elapsed().as_secs_f64() * 1e3;
            let t0 = std::time::Instant::now();
            let east = glam::DVec3::Z.cross(d).normalize();
            let north = d.cross(east);
            let mut acc = 0.0;
            for j in 0..67 {
                for i in 0..67 {
                    let (x, y) = ((i as f64 / 66.0 - 0.5) * 2.0 * half, (j as f64 / 66.0 - 0.5) * 2.0 * half);
                    acc += p.heights((d + east * (x / 6_371_000.0) + north * (y / 6_371_000.0)).normalize())[0];
                }
            }
            let t2 = t0.elapsed().as_secs_f64() * 1e3;
            println!("cutoff {cutoff} m: patch build {t1:.1} ms ({} channel refs), eval 67² {t2:.1} ms  ({acc:.0})", p.channel_count());
        }
    }
    for (name, s) in [("network", net), ("gullies", old)] {
        let g = TileGenerator { sampler: s };
        let mut line = format!("{name:<8}");
        for level in [4u8, 7, 10, 12, 14] {
            let n = 1u32 << level;
            let (x, y) = (((u + 1.0) * 0.5 * n as f64) as u32, ((v + 1.0) * 0.5 * n as f64) as u32);
            let reps = 6;
            let t0 = std::time::Instant::now();
            for k in 0..reps {
                let _ = g.generate(TileId { face, level, x: (x + k) % n, y });
            }
            line += &format!("  L{level}: {:6.1} ms", t0.elapsed().as_secs_f64() * 1000.0 / reps as f64);
        }
        println!("{line}");
    }
}
