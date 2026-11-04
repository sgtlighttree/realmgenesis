//! Shared by the probe examples: build (or load from the probe cache) the world for a seed,
//! with the `RG_*` generation overrides.
#![allow(dead_code)]
use rg_pipeline::{generate, World, WorldParams};
use rg_world::WorldFile;

pub fn env<T: std::str::FromStr>(k: &str) -> Option<T> {
    std::env::var(k).ok().and_then(|v| v.parse().ok())
}

pub fn world(seed: u64) -> World {
    let mut p = WorldParams::default();
    p.tectonics.seed = seed;
    if let Some(r) = env("RG_RES") {
        p.tectonics.res = r;
    }
    if let Some(r) = env("RG_EROSION_RES") {
        p.erosion.res = r;
    }
    if let Some(s) = env("RG_EROSION_STEPS") {
        p.erosion.steps = s;
    }
    if let Some(v) = env("RG_K_VAR") {
        p.erosion.k_variability = v;
    }
    if let Some(v) = env("RG_ROUGH") {
        p.erosion.initial_roughness_m = v;
    }
    if let Some(v) = env("RG_LINEAR") {
        p.erosion.linear_mean_m = v;
    }
    if let Some(v) = env::<u8>("RG_RHO8") {
        p.erosion.rho8 = v != 0;
    }
    if let Some(v) = env("RG_LITH_KM") {
        p.erosion.lithology_scale_km = v;
    }
    if let Some(v) = env::<u8>("RG_CONV_SCALE") {
        p.erosion.scale_by_convergence_rate = v != 0;
    }
    if let Some(v) = env::<u8>("RG_EROSION") {
        p.erosion.enabled = v != 0;
    }
    if let Some(v) = env("RG_JITTER") {
        p.erosion.node_jitter = v;
    }
    if let Some(v) = env("RG_ROUGH_WL") {
        p.erosion.roughness_wavelength_km = v;
    }
    if let Some(v) = env("RG_CONV_SMOOTH") {
        p.erosion.convergence_smoothing_km = v;
    }
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../target/probe-cache");
    let path = dir.join(format!(
        "s{seed}-t{}-e{}-n{}-k{}-r{}-l{}-rho{}-lin{}-cv{}-cs{}-j{}-rw{}-en{}-g{}-b{}.rgw",
        p.tectonics.res,
        p.erosion.res,
        p.erosion.steps,
        p.erosion.k_variability,
        p.erosion.initial_roughness_m,
        p.erosion.lithology_scale_km,
        p.erosion.rho8 as u8,
        p.erosion.linear_mean_m,
        p.erosion.scale_by_convergence_rate as u8,
        p.erosion.convergence_smoothing_km,
        p.erosion.node_jitter,
        p.erosion.roughness_wavelength_km,
        p.erosion.enabled as u8,
        rg_pipeline::file::WORLD_SCHEMA,
        build_stamp()
    ));
    if env::<u8>("RG_NO_CACHE").is_none()
        && let Ok(bytes) = std::fs::read(&path)
        && let Ok(f) = WorldFile::read_from(&bytes[..])
        && let Ok(w) = rg_pipeline::file::from_file(&f)
    {
        eprintln!("loaded {}", path.display());
        return w;
    }
    let t0 = std::time::Instant::now();
    let w = generate(&p);
    eprintln!(
        "generate {:.1} s (tectonics {:.0} ms, erosion {:.0} ms, macro res {})",
        t0.elapsed().as_secs_f64(),
        w.timings.tectonics_ms,
        w.timings.erosion_ms,
        w.tectonics.height_m.res
    );
    let _ = std::fs::create_dir_all(&dir);
    let mut buf = Vec::new();
    rg_pipeline::file::to_file(&w).write_to(&mut buf).unwrap();
    let _ = std::fs::write(&path, buf);
    w
}


/// Build time of this probe binary, so a rebuilt generator never reads a stale cached world.
fn build_stamp() -> u64 {
    std::env::current_exe()
        .and_then(std::fs::metadata)
        .and_then(|m| m.modified())
        .ok()
        .and_then(|t| t.duration_since(std::time::UNIX_EPOCH).ok())
        .map_or(0, |d| d.as_secs())
}
