//! R3.2: importing coastlines. PNG first (Matt, s04), then SVG; GeoJSON later.
//!
//! A mask becomes a [`RasterMask`]: coverage 0 (sea) .. 1 (land) per pixel, placed on the
//! sphere as a whole-world equirect image or a regional patch. The original file bytes travel
//! with it, so the world file can store them (`conditioning.source`) and regenerate the world
//! anywhere.

use rg_tectonics::mask::{Conditioning, Placement, RasterMask, Stroke, StrokeOp};
use rg_world::WorldFile;
use serde_json::{Value, json};
use std::sync::Arc;

/// Decode a PNG into coverage: the alpha channel if the image has a varying one (a
/// transparent sea around opaque land), otherwise brightness (white land on black sea).
pub fn coverage_from_png(bytes: &[u8]) -> Result<(Vec<f32>, u32, u32), String> {
    let mut dec = png::Decoder::new(std::io::Cursor::new(bytes));
    dec.set_transformations(png::Transformations::normalize_to_color8());
    let mut rd = dec.read_info().map_err(|e| format!("PNG: {e}"))?;
    let mut buf = vec![0u8; rd.output_buffer_size().ok_or("PNG: image too large")?];
    let info = rd.next_frame(&mut buf).map_err(|e| format!("PNG: {e}"))?;
    let (w, h) = (info.width, info.height);
    let ch = info.color_type.samples();
    let px = &buf[..info.buffer_size()];
    let has_alpha = matches!(info.color_type, png::ColorType::GrayscaleAlpha | png::ColorType::Rgba);
    let alpha_varies = has_alpha && {
        let first = px[ch - 1];
        px.chunks(ch).any(|p| p[ch - 1] != first)
    };
    let cov: Vec<f32> = px
        .chunks(ch)
        .map(|p| {
            if alpha_varies {
                p[ch - 1] as f32 / 255.0
            } else if ch >= 3 {
                (0.2126 * p[0] as f32 + 0.7152 * p[1] as f32 + 0.0722 * p[2] as f32) / 255.0
            } else {
                p[0] as f32 / 255.0
            }
        })
        .collect();
    if w == 0 || h == 0 {
        return Err("PNG: empty image".into());
    }
    Ok((cov, w, h))
}

/// Build a mask from a file's bytes (`kind`: "png" or "svg").
pub fn mask_from_bytes(bytes: Vec<u8>, kind: &str, place: Placement) -> Result<RasterMask, String> {
    let (data, w, h, place) = match kind {
        "png" => {
            let (d, w, h) = coverage_from_png(&bytes)?;
            (d, w, h, place)
        }
        "svg" => crate::svg::coverage_from_svg(&bytes, &place)?,
        other => return Err(format!("unsupported mask format {other:?}")),
    };
    Ok(RasterMask { data, w, h, place, source: bytes, source_kind: kind.to_string() })
}

/// The placement as the user gave it (an SVG patch's raster scale converted back to km per
/// viewBox unit), for the world file: loading re-runs [`mask_from_bytes`] with it.
pub fn user_placement(m: &RasterMask) -> Placement {
    match (&m.place, m.source_kind.as_str()) {
        (Placement::Patch { center_deg, km_per_px, rotation_deg }, "svg") => {
            let text = std::str::from_utf8(&m.source).unwrap_or("");
            let side = crate::svg::view_side(text).unwrap_or(1.0);
            let scale = m.w.max(m.h) as f64 / side; // raster px per unit
            Placement::Patch { center_deg: *center_deg, km_per_px: km_per_px * scale, rotation_deg: *rotation_deg }
        }
        (p, _) => p.clone(),
    }
}

/// Load a mask file from disk; the format comes from the extension.
pub fn load_mask(path: &std::path::Path, place: Placement) -> Result<RasterMask, String> {
    let kind = match path.extension().and_then(|e| e.to_str()).map(|e| e.to_ascii_lowercase()).as_deref() {
        Some("png") => "png",
        Some("svg") => "svg",
        _ => return Err("pick a .png or .svg file".into()),
    };
    let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    mask_from_bytes(bytes, kind, place)
}

fn op_name(op: StrokeOp) -> &'static str {
    match op {
        StrokeOp::Land => "land",
        StrokeOp::Sea => "sea",
        StrokeOp::Active => "active",
        StrokeOp::Passive => "passive",
    }
}

fn placement_json(p: &Placement) -> Value {
    match p {
        Placement::WorldEquirect => json!({"kind": "world_equirect"}),
        Placement::Patch { center_deg, km_per_px, rotation_deg } => json!({
            "kind": "patch", "center_deg": center_deg, "km_per_px": km_per_px, "rotation_deg": rotation_deg,
        }),
    }
}

fn placement_from(v: &Value) -> Placement {
    if v["kind"].as_str() == Some("patch") {
        let c = |i: usize| v["center_deg"][i].as_f64().unwrap_or(0.0);
        Placement::Patch {
            center_deg: [c(0), c(1)],
            km_per_px: v["km_per_px"].as_f64().unwrap_or(10.0),
            rotation_deg: v["rotation_deg"].as_f64().unwrap_or(0.0),
        }
    } else {
        Placement::WorldEquirect
    }
}

/// The world-file header block for the conditioning (schema 5). The mask's source bytes go
/// in the `conditioning.source` chunk (see [`source_chunk`]).
pub fn conditioning_header(c: &Conditioning) -> Value {
    let mode = match (&c.mask, c.strict) {
        (Some(_), _) => "import",
        (None, true) => "procedural_strict",
        (None, false) => "procedural",
    };
    let strokes: Vec<Value> = c
        .strokes
        .iter()
        .map(|s| {
            json!({
                "op": op_name(s.op),
                "radius_km": s.radius_km,
                "path": s.path.iter().map(|d| [d.x, d.y, d.z]).collect::<Vec<_>>(),
            })
        })
        .collect();
    let mut v = json!({
        "mode": mode,
        "strict": c.strict,
        "outside_procedural": c.outside_procedural,
        "guard_km": c.guard_km,
        "islands": c.islands,
        "islands_guard_km": c.islands_guard_km,
        "strokes": strokes,
    });
    if let Some(m) = &c.mask {
        v["placement"] = placement_json(&user_placement(m));
        v["source_kind"] = json!(m.source_kind);
    }
    v
}

/// The `conditioning.source` chunk (the original mask file), if any.
pub fn source_chunk(c: &Conditioning) -> Option<rg_world::Chunk> {
    c.mask.as_ref().map(|m| rg_world::Chunk::u8s("conditioning.source", &m.source))
}

/// Read the conditioning back. Absent (schema < 5) = procedural, as before R3.
pub fn conditioning_from_file(v: &Value, f: &WorldFile) -> Result<Conditioning, String> {
    let d = Conditioning::default();
    if !v.is_object() {
        return Ok(d);
    }
    let strokes = v["strokes"]
        .as_array()
        .map(|a| {
            a.iter()
                .filter_map(|s| {
                    let op = match s["op"].as_str()? {
                        "land" => StrokeOp::Land,
                        "sea" => StrokeOp::Sea,
                        "active" => StrokeOp::Active,
                        "passive" => StrokeOp::Passive,
                        _ => return None,
                    };
                    let path = s["path"]
                        .as_array()?
                        .iter()
                        .filter_map(|p| Some(glam::DVec3::new(p[0].as_f64()?, p[1].as_f64()?, p[2].as_f64()?)))
                        .collect();
                    Some(Stroke { path, radius_km: s["radius_km"].as_f64().unwrap_or(100.0), op })
                })
                .collect()
        })
        .unwrap_or_default();
    let mask = if v["mode"].as_str() == Some("import") {
        let src = f.chunk("conditioning.source").ok_or("import world without its conditioning.source chunk")?;
        let kind = v["source_kind"].as_str().unwrap_or("png");
        Some(Arc::new(mask_from_bytes(src.bytes.clone(), kind, placement_from(&v["placement"]))?))
    } else {
        None
    };
    Ok(Conditioning {
        strict: v["strict"].as_bool().unwrap_or(mask.is_some()),
        mask,
        outside_procedural: v["outside_procedural"].as_bool().unwrap_or(d.outside_procedural),
        guard_km: v["guard_km"].as_f64().unwrap_or(d.guard_km),
        islands: v["islands"].as_bool().unwrap_or(d.islands),
        islands_guard_km: v["islands_guard_km"].as_f64().unwrap_or(d.islands_guard_km),
        strokes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A tiny PNG in memory: white land square on black sea.
    pub fn test_png(w: u32, h: u32, land: impl Fn(u32, u32) -> bool) -> Vec<u8> {
        let mut out = Vec::new();
        {
            let mut enc = png::Encoder::new(&mut out, w, h);
            enc.set_color(png::ColorType::Grayscale);
            enc.set_depth(png::BitDepth::Eight);
            let mut wr = enc.write_header().unwrap();
            let data: Vec<u8> = (0..h).flat_map(|y| (0..w).map(move |x| (x, y))).map(|(x, y)| if land(x, y) { 255 } else { 0 }).collect();
            wr.write_image_data(&data).unwrap();
        }
        out
    }

    #[test]
    fn png_decodes_to_coverage() {
        let bytes = test_png(8, 4, |x, _| x < 4);
        let (cov, w, h) = coverage_from_png(&bytes).unwrap();
        assert_eq!((w, h), (8, 4));
        assert_eq!(cov[0], 1.0);
        assert_eq!(cov[7], 0.0);
    }
}
