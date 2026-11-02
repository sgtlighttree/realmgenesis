//! `World` ↔ world file. Everything committed is stored as data, never re-derived on load
//! (ARCHITECTURE §5): loading a file reproduces the world exactly, even if the generators
//! change later.

use crate::{Timings, World, WorldParams};
use rg_climate::{ClimateFields, LatLonGrid, PlanetParams};
use rg_geo::CubeField;
use rg_tectonics::{Plate, TectonicSnapshot, TectonicsParams};
use rg_world::{Chunk, WorldFile};
use serde_json::json;

/// Version of the chunk set below. Bump on any change to names or meanings.
/// 2 (s01): erosion params + hydrology chunks. Schema-1 files still load (no hydrology).
pub const WORLD_SCHEMA: u32 = 2;

pub fn to_file(w: &World) -> WorldFile {
    let t = &w.params.tectonics;
    let p = &w.params.planet;
    let plates: Vec<_> = w
        .tectonics
        .plates
        .iter()
        .map(|pl| json!({"pole": [pl.euler_pole.x, pl.euler_pole.y, pl.euler_pole.z], "omega": pl.omega, "continental": pl.continental}))
        .collect();
    let header = json!({
        "app": "RealmGenesis",
        "schema": WORLD_SCHEMA,
        "face_warp": rg_geo::cube::FACE_WARP_ID,
        "tectonics": {
            "seed": t.seed, "res": t.res, "radius_m": t.radius_m, "plate_count": t.plate_count,
            "land_fraction": t.land_fraction, "continent_count_hint": t.continent_count_hint,
            "mountain_scale": t.mountain_scale, "roughness": t.roughness,
        },
        "planet": {
            "axial_tilt_deg": p.axial_tilt_deg, "rotation_hours": p.rotation_hours,
            "solar_constant_rel": p.solar_constant_rel, "global_temp_offset_c": p.global_temp_offset_c,
            "global_precip_mult": p.global_precip_mult,
        },
        "snapshot": {
            "sea_level_offset_m": w.tectonics.sea_level_offset_m,
            "land_fraction": w.tectonics.land_fraction,
            "plates": plates,
        },
        "climate": {"width": w.climate.width, "height": w.climate.height, "model": "tier_a"},
        "erosion": {
            "enabled": w.params.erosion.enabled, "max_res": w.params.erosion.max_res,
            "duration_myr": w.params.erosion.duration_myr, "steps": w.params.erosion.steps,
            "k": w.params.erosion.k, "m": w.params.erosion.m,
            "uplift_max_mm_yr": w.params.erosion.uplift_max_mm_yr,
            "uplift_base_mm_yr": w.params.erosion.uplift_base_mm_yr,
            "uplift_exponent": w.params.erosion.uplift_exponent,
            "river_min_area_km2": w.params.erosion.river_min_area_km2,
            "hypsometry_scale": w.params.erosion.hypsometry_scale,
            "hydrology_res": w.hydrology.as_ref().map(|h| h.drainage_km2.res),
        },
    });
    let c = &w.climate;
    let mut hydro_chunks = Vec::new();
    if let Some(h) = &w.hydrology {
        hydro_chunks.push(Chunk::f32s("hydrology.drainage_km2", &h.drainage_km2.data));
        let fd: Vec<f32> = h.flow_dir.data.iter().flat_map(|v| v.iter().copied()).collect();
        hydro_chunks.push(Chunk::f32s("hydrology.flow_dir", &fd));
        hydro_chunks.push(Chunk::f32s("hydrology.lake_depth_m", &h.lake_depth_m.data));
        let rv: Vec<f32> = h
            .rivers
            .iter()
            .flat_map(|r| [r.from.x as f32, r.from.y as f32, r.from.z as f32, r.to.x as f32, r.to.y as f32, r.to.z as f32, r.area_km2])
            .collect();
        hydro_chunks.push(Chunk::f32s("hydrology.rivers", &rv));
    }
    let flat = |v: &Vec<[f32; 12]>| v.iter().flat_map(|m| m.iter().copied()).collect::<Vec<f32>>();
    WorldFile {
        header,
        chunks: vec![
            Chunk::f32s("tectonics.height_m", &w.tectonics.height_m.data),
            Chunk::u16s("tectonics.plate_id", &w.tectonics.plate_id.data),
            Chunk::u8s("tectonics.crust", &w.tectonics.crust.data),
            Chunk::u8s("tectonics.boundary", &w.tectonics.boundary.data),
            Chunk::f32s("climate.heights", &w.climate_heights.data),
            Chunk::f32s("climate.temp_c", &flat(&c.temp_c)),
            Chunk::f32s("climate.precip_mm", &flat(&c.precip_mm)),
            Chunk::f32s("climate.wind_u", &flat(&c.wind_u)),
            Chunk::f32s("climate.wind_v", &flat(&c.wind_v)),
            Chunk::u8s("climate.koppen", &c.koppen.data),
        ]
        .into_iter()
        .chain(hydro_chunks)
        .collect(),
    }
}

pub fn from_file(f: &WorldFile) -> Result<World, String> {
    let h = &f.header;
    let schema = h["schema"].as_u64().unwrap_or(0);
    if !(1..=WORLD_SCHEMA as u64).contains(&schema) {
        return Err(format!("world schema {schema} isn't supported (this build reads 1..={WORLD_SCHEMA})"));
    }
    let num = |v: &serde_json::Value, k: &str| v[k].as_f64().ok_or(format!("missing {k}"));
    let th = &h["tectonics"];
    let tectonics = TectonicsParams {
        seed: th["seed"].as_u64().ok_or("missing seed")?,
        res: num(th, "res")? as u32,
        radius_m: num(th, "radius_m")?,
        plate_count: num(th, "plate_count")? as u32,
        land_fraction: num(th, "land_fraction")? as f32,
        continent_count_hint: num(th, "continent_count_hint")? as u32,
        mountain_scale: num(th, "mountain_scale")? as f32,
        roughness: num(th, "roughness")? as f32,
    };
    let ph = &h["planet"];
    let planet = PlanetParams {
        radius_m: tectonics.radius_m,
        axial_tilt_deg: num(ph, "axial_tilt_deg")?,
        rotation_hours: num(ph, "rotation_hours")?,
        solar_constant_rel: num(ph, "solar_constant_rel")?,
        global_temp_offset_c: num(ph, "global_temp_offset_c")?,
        global_precip_mult: num(ph, "global_precip_mult")?,
    };
    let get = |name: &str| f.chunk(name).ok_or(format!("missing chunk {name}"));
    let res = tectonics.res;
    let cube_f32 = |name: &str| -> Result<CubeField<f32>, String> {
        let data = get(name)?.as_f32s().ok_or(format!("{name}: wrong type"))?;
        check_len(name, data.len(), 6 * (res as usize + 1).pow(2))?;
        Ok(CubeField { res, data })
    };
    let cube_u8 = |name: &str| -> Result<CubeField<u8>, String> {
        let data = get(name)?.bytes.clone();
        check_len(name, data.len(), 6 * (res as usize + 1).pow(2))?;
        Ok(CubeField { res, data })
    };
    let plate_data = get("tectonics.plate_id")?.as_u16s().ok_or("plate_id: wrong type")?;
    check_len("tectonics.plate_id", plate_data.len(), 6 * (res as usize + 1).pow(2))?;
    let sh = &h["snapshot"];
    let plates = sh["plates"]
        .as_array()
        .ok_or("missing plates")?
        .iter()
        .map(|p| {
            let pole = p["pole"].as_array().ok_or("bad plate")?;
            let v = |i: usize| pole.get(i).and_then(|x| x.as_f64()).ok_or("bad pole");
            Ok(Plate {
                euler_pole: glam::DVec3::new(v(0)?, v(1)?, v(2)?),
                omega: p["omega"].as_f64().ok_or("bad omega")?,
                continental: p["continental"].as_bool().unwrap_or(false),
            })
        })
        .collect::<Result<Vec<_>, &str>>()?;
    let snapshot = TectonicSnapshot {
        height_m: cube_f32("tectonics.height_m")?,
        plate_id: CubeField { res, data: plate_data },
        crust: cube_u8("tectonics.crust")?,
        boundary: cube_u8("tectonics.boundary")?,
        plates,
        sea_level_offset_m: num(sh, "sea_level_offset_m")? as f32,
        land_fraction: num(sh, "land_fraction")? as f32,
    };
    let ch = &h["climate"];
    let (cw, chh) = (num(ch, "width")? as u32, num(ch, "height")? as u32);
    let n = (cw * chh) as usize;
    let monthly = |name: &str| -> Result<Vec<[f32; 12]>, String> {
        let flat = get(name)?.as_f32s().ok_or(format!("{name}: wrong type"))?;
        check_len(name, flat.len(), n * 12)?;
        Ok(flat.as_chunks::<12>().0.to_vec())
    };
    let heights = get("climate.heights")?.as_f32s().ok_or("climate.heights: wrong type")?;
    check_len("climate.heights", heights.len(), n)?;
    let koppen = get("climate.koppen")?.bytes.clone();
    check_len("climate.koppen", koppen.len(), n)?;
    let climate = ClimateFields {
        width: cw,
        height: chh,
        temp_c: monthly("climate.temp_c")?,
        precip_mm: monthly("climate.precip_mm")?,
        koppen: LatLonGrid::from_vec(cw, chh, koppen),
        wind_u: monthly("climate.wind_u")?,
        wind_v: monthly("climate.wind_v")?,
    };
    let eh = &h["erosion"];
    let erosion = if schema >= 2 {
        rg_erosion::ErosionParams {
            enabled: eh["enabled"].as_bool().unwrap_or(true),
            max_res: num(eh, "max_res")? as u32,
            duration_myr: num(eh, "duration_myr")?,
            steps: num(eh, "steps")? as u32,
            k: num(eh, "k")?,
            m: num(eh, "m")?,
            uplift_max_mm_yr: num(eh, "uplift_max_mm_yr")?,
            uplift_base_mm_yr: num(eh, "uplift_base_mm_yr")?,
            uplift_exponent: num(eh, "uplift_exponent")?,
            river_min_area_km2: num(eh, "river_min_area_km2")?,
            hypsometry_scale: eh["hypsometry_scale"].as_f64().unwrap_or(0.0),
        }
    } else {
        rg_erosion::ErosionParams { enabled: false, ..Default::default() }
    };
    let hydrology = match eh["hydrology_res"].as_u64() {
        Some(hres) if schema >= 2 => {
            let hres = hres as u32;
            let samples = 6 * (hres as usize + 1).pow(2);
            let f32_field = |name: &str| -> Result<CubeField<f32>, String> {
                let data = get(name)?.as_f32s().ok_or(format!("{name}: wrong type"))?;
                check_len(name, data.len(), samples)?;
                Ok(CubeField { res: hres, data })
            };
            let fd = get("hydrology.flow_dir")?.as_f32s().ok_or("flow_dir: wrong type")?;
            check_len("hydrology.flow_dir", fd.len(), samples * 3)?;
            let rv = get("hydrology.rivers")?.as_f32s().ok_or("rivers: wrong type")?;
            let v3 = |a: &[f32]| glam::DVec3::new(a[0] as f64, a[1] as f64, a[2] as f64);
            Some(crate::Hydrology {
                drainage_km2: f32_field("hydrology.drainage_km2")?,
                flow_dir: CubeField { res: hres, data: fd.as_chunks::<3>().0.to_vec() },
                lake_depth_m: f32_field("hydrology.lake_depth_m")?,
                rivers: rv
                    .as_chunks::<7>()
                    .0
                    .iter()
                    .map(|r| rg_erosion::RiverSegment { from: v3(&r[0..3]), to: v3(&r[3..6]), area_km2: r[6] })
                    .collect(),
                stats: Default::default(),
            })
        }
        _ => None,
    };
    Ok(World {
        hydrology,
        params: WorldParams { tectonics, erosion, planet },
        tectonics: snapshot,
        climate_heights: LatLonGrid::from_vec(cw, chh, heights),
        climate,
        timings: Timings::default(),
    })
}

fn check_len(name: &str, got: usize, want: usize) -> Result<(), String> {
    if got == want { Ok(()) } else { Err(format!("{name}: {got} values, expected {want}")) }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn save_load_save_is_byte_identical() {
        let mut params = WorldParams::default();
        params.tectonics.res = 32;
        params.tectonics.seed = 7;
        let w = crate::generate(&params);
        let mut a = Vec::new();
        to_file(&w).write_to(&mut a).unwrap();
        let loaded = from_file(&WorldFile::read_from(&a[..]).unwrap()).unwrap();
        let mut b = Vec::new();
        to_file(&loaded).write_to(&mut b).unwrap();
        if let Some(i) = a.iter().zip(&b).position(|(x, y)| x != y) {
            let lo = i.saturating_sub(60);
            panic!(
                "first difference at byte {i}:\n a: {}\n b: {}",
                String::from_utf8_lossy(&a[lo..(i + 60).min(a.len())]),
                String::from_utf8_lossy(&b[lo..(i + 60).min(b.len())])
            );
        }
        assert!(a.len() == b.len(), "length {} vs {}", a.len(), b.len());
        assert!(loaded.tectonics.height_m.data == w.tectonics.height_m.data);
    }

    /// Golden determinism: generating twice gives byte-identical files.
    #[test]
    fn generation_is_deterministic() {
        let mut params = WorldParams::default();
        params.tectonics.res = 32;
        let (mut a, mut b) = (Vec::new(), Vec::new());
        to_file(&crate::generate(&params)).write_to(&mut a).unwrap();
        to_file(&crate::generate(&params)).write_to(&mut b).unwrap();
        assert!(a == b, "two generations differ");
    }
}
