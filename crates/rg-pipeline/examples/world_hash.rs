//! Fingerprint of a world's committed fields, to prove a change left a preset untouched
//! (s04 gate: "Normal is byte-identical after the canonical/refine change").
//!
//! `cargo run -p rg-pipeline --release --example world_hash -- <seed>` (generates, `RG_*`
//! overrides as in `common/mod.rs`, honours `RG_NO_CACHE`), or `-- --file <path.rgw>`.
#[path = "common/mod.rs"]
mod common;

fn fnv(bytes: impl Iterator<Item = u8>) -> u64 {
    bytes.fold(0xcbf2_9ce4_8422_2325u64, |h, b| (h ^ b as u64).wrapping_mul(0x0100_0000_01b3))
}

fn main() {
    let a: Vec<String> = std::env::args().skip(1).collect();
    let w = if a[0] == "--file" {
        let bytes = std::fs::read(&a[1]).unwrap();
        let f = rg_world::WorldFile::read_from(&bytes[..]).unwrap();
        rg_pipeline::file::from_file(&f).unwrap()
    } else {
        common::world(a[0].parse().unwrap())
    };
    let f32s = |v: &[f32]| fnv(v.iter().flat_map(|x| x.to_bits().to_le_bytes()));
    println!("height_m   {:016x} (res {})", f32s(&w.tectonics.height_m.data), w.tectonics.height_m.res);
    if let Some(h) = &w.hydrology {
        println!("drainage   {:016x}", f32s(&h.drainage_km2.data));
        println!("lake_depth {:016x}", f32s(&h.lake_depth_m.data));
        println!("rivers     {}", h.rivers.len());
    }
}
