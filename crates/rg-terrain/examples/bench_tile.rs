//! `cargo run -p rg-terrain --release --example bench_tile` — per-tile generation cost by level.
use rg_geo::{CubeFace, CubeField, TileId};
use rg_terrain::{TerrainSampler, TileGenerator};

fn main() {
    let f = std::sync::Arc::new(CubeField::from_fn(320, |d| ((d.z * 3.0).sin() * 3000.0) as f32));
    let g = TileGenerator { sampler: TerrainSampler::new(f, 6_371_000.0) };
    for level in [0u8, 2, 4, 6] {
        let n = 1u32 << level;
        let t0 = std::time::Instant::now();
        let reps = 20;
        for k in 0..reps {
            let _ = g.generate(TileId { face: CubeFace::PosZ, level, x: (n / 2 + k) % n, y: n / 3 });
        }
        println!("L{level:>2}: {:.2} ms/tile", t0.elapsed().as_secs_f64() * 1000.0 / reps as f64);
    }
}
