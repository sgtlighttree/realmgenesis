//! R3.2: importing coastlines. PNG first (Matt, s04), then SVG; GeoJSON later.
//!
//! A mask becomes a [`RasterMask`]: coverage 0 (sea) .. 1 (land) per pixel, placed on the
//! sphere as a whole-world equirect image or a regional patch. The original file bytes travel
//! with it, so the world file can store them (`conditioning.source`) and regenerate the world
//! anywhere.

use rg_tectonics::mask::{Channel, Conditioning, DecodeOpts, Elevation, Fidelity, MaskAdjust, Placement, RasterMask, Stroke, StrokeOp};
use rg_world::WorldFile;
use serde_json::{Value, json};
use std::sync::Arc;

/// Decode a PNG into raw values 0..1: alpha or brightness, per `channel` (`Auto` = alpha if
/// the image has a varying one, a transparent sea around opaque land, else brightness).
/// Returns the values, the size and the channel actually used.
pub fn coverage_from_png(bytes: &[u8], channel: Channel) -> Result<(Vec<f32>, u32, u32, Channel), String> {
    raw_from_png(bytes, channel, None, 1)
}

/// A PNG's size in pixels, from its header (no decode).
pub fn png_size(bytes: &[u8]) -> Result<(u32, u32), String> {
    let rd = png::Decoder::new(std::io::Cursor::new(bytes)).read_info().map_err(|e| format!("PNG: {e}"))?;
    Ok((rd.info().width, rd.info().height))
}

/// Calls `f(y, row)` for each decoded row (8-bit samples, `samples` per pixel) until it
/// returns false. Rows stream one at a time, so a huge map never needs a whole-frame buffer
/// (s07: Klerga V2 is 11200×5600 RGBA, 250 MB). Interlaced PNGs (rare) fall back to one
/// frame. Returns (width, height, samples per pixel, has alpha).
fn for_each_png_row(bytes: &[u8], mut f: impl FnMut(u32, &[u8]) -> bool) -> Result<(u32, u32, usize, bool), String> {
    let mut dec = png::Decoder::new(std::io::Cursor::new(bytes));
    dec.set_transformations(png::Transformations::normalize_to_color8());
    let mut rd = dec.read_info().map_err(|e| format!("PNG: {e}"))?;
    let (w, h) = (rd.info().width, rd.info().height);
    if w == 0 || h == 0 {
        return Err("PNG: empty image".into());
    }
    let color = rd.output_color_type().0;
    let ch = color.samples();
    let has_alpha = matches!(color, png::ColorType::GrayscaleAlpha | png::ColorType::Rgba);
    let line = w as usize * ch;
    if rd.info().interlaced {
        let mut buf = vec![0u8; rd.output_buffer_size().ok_or("PNG: image too large")?];
        let info = rd.next_frame(&mut buf).map_err(|e| format!("PNG: {e}"))?;
        for y in 0..h {
            let at = y as usize * info.line_size;
            if !f(y, &buf[at..at + line]) {
                break;
            }
        }
    } else {
        let mut y = 0;
        while let Some(row) = rd.next_row().map_err(|e| format!("PNG: {e}"))? {
            if !f(y, &row.data()[..line]) {
                break;
            }
            y += 1;
        }
    }
    Ok((w, h, ch, has_alpha))
}

/// Decode a PNG into raw values 0..1 (see [`coverage_from_png`] for `channel`), keeping the
/// source-pixel rectangle `crop` (clamped like [`crop_raster`]) and box-averaging blocks of
/// `step`×`step` pixels (1 = full resolution; partial blocks at the right and bottom edges
/// average what they hold). At step 1 the values are exactly those of the whole-frame decode
/// plus `crop_raster`. Sequential, so the same bytes always give the same values.
pub fn raw_from_png(bytes: &[u8], channel: Channel, crop: Option<[u32; 4]>, step: u32) -> Result<(Vec<f32>, u32, u32, Channel), String> {
    let step = step.max(1) as usize;
    // The header pass: size and layout (the closure stops at once).
    let (w, h, ch, has_alpha) = for_each_png_row(bytes, |_, _| false)?;
    let use_alpha = match channel {
        Channel::Alpha => has_alpha,
        Channel::Brightness => false,
        // Alpha only if it varies somewhere: a pass that stores nothing and stops early.
        Channel::Auto => {
            let mut first: Option<u8> = None;
            let mut varies = false;
            if has_alpha {
                for_each_png_row(bytes, |_, row| {
                    for p in row.chunks(ch) {
                        let a = p[ch - 1];
                        if *first.get_or_insert(a) != a {
                            varies = true;
                            return false;
                        }
                    }
                    true
                })?;
            }
            varies
        }
    };
    let (x0, y0, cw, crows) = match crop {
        Some(r) => {
            let x0 = r[0].min(w - 1);
            let y0 = r[1].min(h - 1);
            (x0 as usize, y0 as usize, r[2].clamp(1, w - x0) as usize, r[3].clamp(1, h - y0) as usize)
        }
        None => (0, 0, w as usize, h as usize),
    };
    let (ow, oh) = (cw.div_ceil(step), crows.div_ceil(step));
    let mut out: Vec<f32> = Vec::with_capacity(ow * oh);
    let mut acc = vec![0.0f64; ow];
    for_each_png_row(bytes, |y, row| {
        let Some(ry) = (y as usize).checked_sub(y0) else { return true };
        if ry >= crows {
            return false;
        }
        for rx in 0..cw {
            let p = &row[(x0 + rx) * ch..(x0 + rx + 1) * ch];
            let v = if use_alpha {
                p[ch - 1] as f32 / 255.0
            } else if ch >= 3 {
                (0.2126 * p[0] as f32 + 0.7152 * p[1] as f32 + 0.0722 * p[2] as f32) / 255.0
            } else {
                p[0] as f32 / 255.0
            };
            acc[rx / step] += v as f64;
        }
        // The last row of a block (or of the crop): write the block row out.
        if (ry + 1) % step == 0 || ry + 1 == crows {
            let rows = ry % step + 1;
            for (ox, a) in acc.iter_mut().enumerate() {
                let cols = step.min(cw - ox * step);
                out.push((*a / (rows * cols) as f64) as f32);
                *a = 0.0;
            }
        }
        true
    })?;
    Ok((out, ow as u32, oh as u32, if use_alpha { Channel::Alpha } else { Channel::Brightness }))
}

/// Keep the source-pixel rectangle `[x, y, w, h]` (clamped to the image).
pub fn crop_raster(data: &[f32], w: u32, h: u32, rect: [u32; 4]) -> (Vec<f32>, u32, u32) {
    let x0 = rect[0].min(w - 1);
    let y0 = rect[1].min(h - 1);
    let cw = rect[2].clamp(1, w - x0);
    let ch = rect[3].clamp(1, h - y0);
    let out = (y0..y0 + ch).flat_map(|y| (x0..x0 + cw).map(move |x| data[(y * w + x) as usize])).collect();
    (out, cw, ch)
}

/// Build a mask from a file's bytes (`kind`: "png" or "svg") with the default land test.
pub fn mask_from_bytes(bytes: Vec<u8>, kind: &str, place: Placement) -> Result<RasterMask, String> {
    mask_from_bytes_with(bytes, kind, place, DecodeOpts::default(), MaskAdjust::default())
}

/// Build a mask from a file's bytes. `decode` (channel, crop) and `adjust` (threshold,
/// invert) apply to PNGs; an SVG's fill is already a clean land/sea split.
pub fn mask_from_bytes_with(
    bytes: Vec<u8>,
    kind: &str,
    place: Placement,
    decode: DecodeOpts,
    adjust: MaskAdjust,
) -> Result<RasterMask, String> {
    layer_from_bytes(bytes, kind, place, decode, adjust, None, None)
}

/// Longest raster side kept for a heightmap layer (s07): the source is box-downsampled by an
/// integer step to at most this. 6,144 keeps Klerga V2 (11,200 px around the world) at
/// 7.1 km/px, finer than Ultra's ~10 km cell.
pub const HEIGHTMAP_MAX_SIDE: u32 = 6144;

/// The integer downsample step for a heightmap of this source (after `crop`).
pub fn heightmap_step(bytes: &[u8], crop: Option<[u32; 4]>) -> Result<u32, String> {
    let (w, h) = png_size(bytes)?;
    let (cw, ch) = match crop {
        Some(r) => {
            let (x0, y0) = (r[0].min(w - 1), r[1].min(h - 1));
            (r[2].clamp(1, w - x0), r[3].clamp(1, h - y0))
        }
        None => (w, h),
    };
    Ok(cw.max(ch).div_ceil(HEIGHTMAP_MAX_SIDE).max(1))
}

/// Build a layer from a file's bytes: a coastline (`elevation` None) or a heightmap (PNG
/// only; downsampled by `step`, or by [`heightmap_step`] when `None`, with its height
/// pyramid built). `place` is in the user's units (km per **source** pixel), as stored in
/// the world file; the mask holds it per raster pixel.
pub fn layer_from_bytes(
    bytes: Vec<u8>,
    kind: &str,
    place: Placement,
    decode: DecodeOpts,
    adjust: MaskAdjust,
    elevation: Option<Elevation>,
    step: Option<u32>,
) -> Result<RasterMask, String> {
    let (data, w, h, place, used, decode, adjust, step) = match kind {
        "png" => {
            let step = match (elevation, step) {
                (None, _) => 1,
                (Some(_), Some(s)) => s.max(1),
                (Some(_), None) => heightmap_step(&bytes, decode.crop)?,
            };
            let channel = effective_channel(elevation.is_some(), decode.channel);
            let (d, w, h, used) = raw_from_png(&bytes, channel, decode.crop, step)?;
            (d, w, h, place, used, decode, adjust, step)
        }
        "svg" if elevation.is_some() => return Err("a heightmap must be a PNG".into()),
        "svg" => {
            let (d, w, h, place) = crate::svg::coverage_from_svg(&bytes, &place)?;
            (d, w, h, place, Channel::Brightness, DecodeOpts::default(), MaskAdjust::default(), 1)
        }
        other => return Err(format!("unsupported mask format {other:?}")),
    };
    let mut m = RasterMask {
        data: Arc::new(data),
        w,
        h,
        place,
        adjust,
        decode,
        channel_used: used,
        enabled: true,
        name: String::new(),
        source: Arc::new(bytes),
        source_kind: kind.to_string(),
        step,
        elevation,
        heights: None,
    };
    if step > 1 {
        m.place = raster_placement(&m, m.place.clone());
    }
    m.build_heights();
    Ok(m)
}

/// The channel a layer actually decodes. A heightmap's greys are its brightness: Auto would
/// pick a varying alpha (Klerga V1 has one, s07) and flatten every height. Alpha only if
/// asked for. The preview uses this too, so it shows what Generate reads.
pub fn effective_channel(heightmap: bool, channel: Channel) -> Channel {
    if heightmap && channel == Channel::Auto { Channel::Brightness } else { channel }
}

/// Does a heightmap's sea hold bathymetry? True when the greys below the coast spread over
/// more than ~5 levels (p5..p95); a flat sea (Azgaar's 0) would otherwise turn into a cliff
/// to the full depth at the coast (s07 advisor). A 256-bin histogram, so the panel can
/// re-evaluate it whenever the sea-level grey moves.
pub fn has_bathymetry(m: &RasterMask) -> bool {
    let mut hist = [0u64; 256];
    for &v in m.data.iter() {
        let g = if m.adjust.invert { 1.0 - v } else { v };
        hist[(g.clamp(0.0, 1.0) * 255.0).round() as usize] += 1;
    }
    let top = (m.adjust.threshold.clamp(0.0, 1.0) * 255.0).floor() as usize;
    let sea: u64 = hist[..=top].iter().sum();
    if sea < 16 {
        return false;
    }
    let quantile = |f: f64| {
        let want = (sea as f64 * f) as u64;
        let mut acc = 0;
        hist[..=top].iter().position(|&c| {
            acc += c;
            acc > want
        }).unwrap_or(top)
    };
    quantile(0.95) - quantile(0.05) > 5
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
        // A downsampled heightmap: km per raster pixel → per source pixel (s07).
        (Placement::Patch { center_deg, km_per_px, rotation_deg }, _) if m.step > 1 => {
            Placement::Patch { center_deg: *center_deg, km_per_px: km_per_px / m.step as f64, rotation_deg: *rotation_deg }
        }
        (p, _) => p.clone(),
    }
}

/// The inverse of [`user_placement`]: the placement to store on the mask for a placement in
/// the user's units (an SVG patch's km per viewBox unit becomes km per raster pixel).
pub fn raster_placement(m: &RasterMask, user: Placement) -> Placement {
    match (user, m.source_kind.as_str()) {
        (Placement::Patch { center_deg, km_per_px, rotation_deg }, "svg") => {
            let text = std::str::from_utf8(&m.source).unwrap_or("");
            let side = crate::svg::view_side(text).unwrap_or(1.0);
            let scale = m.w.max(m.h) as f64 / side; // raster px per unit
            Placement::Patch { center_deg, km_per_px: km_per_px / scale, rotation_deg }
        }
        (Placement::Patch { center_deg, km_per_px, rotation_deg }, _) if m.step > 1 => {
            Placement::Patch { center_deg, km_per_px: km_per_px * m.step as f64, rotation_deg }
        }
        (p, _) => p,
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
    let mode = match (c.masks.is_empty(), c.strict) {
        (false, _) => "import",
        (true, true) => "procedural_strict",
        (true, false) => "procedural",
    };
    let strokes: Vec<Value> = c
        .strokes
        .iter()
        .map(|s| {
            let mut j = json!({
                "op": op_name(s.op),
                "radius_km": s.radius_km,
                "path": s.path.iter().map(|d| [d.x, d.y, d.z]).collect::<Vec<_>>(),
            });
            // s06: optional, absent = every seed / enabled (older files).
            if let Some(seed) = s.seed {
                j["seed"] = json!(seed);
            }
            if !s.enabled {
                j["enabled"] = json!(false);
            }
            j
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
    // s06: several imports, each with its own source chunk `conditioning.source.{i}`.
    if !c.masks.is_empty() {
        v["masks"] = Value::Array(c.masks.iter().map(|m| mask_json(m)).collect());
    }
    v
}

/// One import's header entry: placement, format, land test, decode options, name, on/off.
fn mask_json(m: &RasterMask) -> Value {
    let mut v = json!({"placement": placement_json(&user_placement(m)), "source_kind": m.source_kind});
    // Optional fields; absent = defaults.
    if m.adjust != MaskAdjust::default() {
        v["adjust"] = json!({"threshold": m.adjust.threshold, "invert": m.adjust.invert});
    }
    if m.decode != DecodeOpts::default() {
        let channel = match m.decode.channel {
            Channel::Auto => "auto",
            Channel::Brightness => "brightness",
            Channel::Alpha => "alpha",
        };
        v["decode"] = json!({"channel": channel, "crop": m.decode.crop});
    }
    if !m.enabled {
        v["enabled"] = json!(false);
    }
    if !m.name.is_empty() {
        v["name"] = json!(m.name);
    }
    if let Some(e) = m.elevation {
        let fidelity = match e.fidelity {
            Fidelity::Exact => "exact",
            Fidelity::Guided => "guided",
        };
        v["elevation"] = json!({
            "max_m": e.max_m, "gamma": e.gamma, "sea_floor": e.sea_floor, "sea_depth_m": e.sea_depth_m,
            "smooth_px": e.smooth_px, "fidelity": fidelity,
            // The raster the heights came from: the step is stored, so a later cap change
            // can't alter an old world (s07).
            "step": m.step, "downsample_ver": 1,
        });
    }
    v
}

/// Rebuild one import from its header entry and its source bytes. Also reads the s05/early
/// s06 single-import layout, whose fields sat at the top of the conditioning block.
fn mask_from_json(v: &Value, bytes: Vec<u8>) -> Result<RasterMask, String> {
    let kind = v["source_kind"].as_str().unwrap_or("png");
    let d = MaskAdjust::default();
    let adjust = MaskAdjust {
        threshold: v["adjust"]["threshold"].as_f64().map_or(d.threshold, |t| t as f32),
        invert: v["adjust"]["invert"].as_bool().unwrap_or(d.invert),
    };
    let decode = DecodeOpts {
        channel: match v["decode"]["channel"].as_str() {
            Some("brightness") => Channel::Brightness,
            Some("alpha") => Channel::Alpha,
            _ => Channel::Auto,
        },
        crop: v["decode"]["crop"].as_array().and_then(|a| {
            let n: Vec<u32> = a.iter().filter_map(|x| x.as_u64().map(|x| x as u32)).collect();
            (n.len() == 4).then(|| [n[0], n[1], n[2], n[3]])
        }),
    };
    let ev = &v["elevation"];
    let elevation = ev.is_object().then(|| {
        let d = Elevation::default();
        let f = |k: &str, def: f32| ev[k].as_f64().map_or(def, |x| x as f32);
        Elevation {
            max_m: f("max_m", d.max_m),
            gamma: f("gamma", d.gamma),
            sea_floor: ev["sea_floor"].as_bool().unwrap_or(d.sea_floor),
            sea_depth_m: f("sea_depth_m", d.sea_depth_m),
            smooth_px: f("smooth_px", d.smooth_px),
            fidelity: if ev["fidelity"].as_str() == Some("exact") { Fidelity::Exact } else { Fidelity::Guided },
        }
    });
    let step = ev["step"].as_u64().map(|s| s as u32);
    let mut m = layer_from_bytes(bytes, kind, placement_from(&v["placement"]), decode, adjust, elevation, step)?;
    m.enabled = v["enabled"].as_bool().unwrap_or(true);
    m.name = v["name"].as_str().unwrap_or("").to_string();
    Ok(m)
}

/// The source chunks (the original mask files), one per import.
pub fn source_chunks(c: &Conditioning) -> Vec<rg_world::Chunk> {
    c.masks.iter().enumerate().map(|(i, m)| rg_world::Chunk::u8s(&format!("conditioning.source.{i}"), &m.source)).collect()
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
                    Some(Stroke {
                        path,
                        radius_km: s["radius_km"].as_f64().unwrap_or(100.0),
                        op,
                        seed: s["seed"].as_u64(),
                        enabled: s["enabled"].as_bool().unwrap_or(true),
                    })
                })
                .collect()
        })
        .unwrap_or_default();
    let masks: Vec<Arc<RasterMask>> = if let Some(list) = v["masks"].as_array() {
        list.iter()
            .enumerate()
            .map(|(i, mv)| {
                let name = format!("conditioning.source.{i}");
                let src = f.chunk(&name).ok_or(format!("import world without its {name} chunk"))?;
                mask_from_json(mv, src.bytes.clone()).map(Arc::new)
            })
            .collect::<Result<_, String>>()?
    } else if v["mode"].as_str() == Some("import") {
        // s05 / early s06: one import, fields at the top level, chunk `conditioning.source`.
        let src = f.chunk("conditioning.source").ok_or("import world without its conditioning.source chunk")?;
        vec![Arc::new(mask_from_json(v, src.bytes.clone())?)]
    } else {
        Vec::new()
    };
    Ok(Conditioning {
        strict: v["strict"].as_bool().unwrap_or(!masks.is_empty()),
        masks,
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
        let (cov, w, h, used) = coverage_from_png(&bytes, Channel::Auto).unwrap();
        assert_eq!(used, Channel::Brightness);
        assert_eq!((w, h), (8, 4));
        assert_eq!(cov[0], 1.0);
        assert_eq!(cov[7], 0.0);
    }

    /// An RGBA PNG with pseudo-random colours; alpha 255 except where `alpha` says.
    fn rgba_png(w: u32, h: u32, alpha: impl Fn(u32, u32) -> u8) -> Vec<u8> {
        let mut out = Vec::new();
        {
            let mut enc = png::Encoder::new(&mut out, w, h);
            enc.set_color(png::ColorType::Rgba);
            enc.set_depth(png::BitDepth::Eight);
            let mut wr = enc.write_header().unwrap();
            let mut data = Vec::new();
            for y in 0..h {
                for x in 0..w {
                    let k = (x * 7919 + y * 104_729) as u64;
                    let hsh = k.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> 40;
                    data.extend_from_slice(&[hsh as u8, (hsh >> 8) as u8, (hsh >> 16) as u8, alpha(x, y)]);
                }
            }
            wr.write_image_data(&data).unwrap();
        }
        out
    }

    /// The pre-s07 whole-frame decoder, kept here as the reference for byte identity.
    fn reference_decode(bytes: &[u8], channel: Channel) -> (Vec<f32>, u32, u32, Channel) {
        let mut dec = png::Decoder::new(std::io::Cursor::new(bytes));
        dec.set_transformations(png::Transformations::normalize_to_color8());
        let mut rd = dec.read_info().unwrap();
        let mut buf = vec![0u8; rd.output_buffer_size().unwrap()];
        let info = rd.next_frame(&mut buf).unwrap();
        let ch = info.color_type.samples();
        let px = &buf[..info.buffer_size()];
        let has_alpha = matches!(info.color_type, png::ColorType::GrayscaleAlpha | png::ColorType::Rgba);
        let use_alpha = match channel {
            Channel::Alpha => has_alpha,
            Channel::Brightness => false,
            Channel::Auto => has_alpha && px.chunks(ch).any(|p| p[ch - 1] != px[ch - 1]),
        };
        let cov = px
            .chunks(ch)
            .map(|p| {
                if use_alpha {
                    p[ch - 1] as f32 / 255.0
                } else if ch >= 3 {
                    (0.2126 * p[0] as f32 + 0.7152 * p[1] as f32 + 0.0722 * p[2] as f32) / 255.0
                } else {
                    p[0] as f32 / 255.0
                }
            })
            .collect();
        (cov, info.width, info.height, if use_alpha { Channel::Alpha } else { Channel::Brightness })
    }

    /// s07: the streaming decoder gives exactly the old values (crop included), so imported
    /// worlds stay byte-identical; Auto still scans the whole image for varying alpha.
    #[test]
    fn streaming_decode_matches_whole_frame() {
        let opaque = rgba_png(37, 23, |_, _| 255);
        let late_alpha = rgba_png(37, 23, |x, y| if y == 22 && x == 36 { 10 } else { 255 });
        for bytes in [&opaque, &late_alpha] {
            for ch in [Channel::Auto, Channel::Brightness, Channel::Alpha] {
                let (r, rw, rh, ru) = reference_decode(bytes, ch);
                let (d, w, h, u) = raw_from_png(bytes, ch, None, 1).unwrap();
                assert_eq!((w, h, u), (rw, rh, ru));
                assert!(d.iter().zip(&r).all(|(a, b)| a.to_bits() == b.to_bits()), "{ch:?}");
                let rect = [5, 3, 20, 100];
                let (cr, cw, chh) = crop_raster(&r, rw, rh, rect);
                let (cd, w, h, _) = raw_from_png(bytes, ch, Some(rect), 1).unwrap();
                assert_eq!((w, h), (cw, chh));
                assert!(cd.iter().zip(&cr).all(|(a, b)| a.to_bits() == b.to_bits()));
            }
        }
        assert_eq!(raw_from_png(&late_alpha, Channel::Auto, None, 1).unwrap().3, Channel::Alpha);
        assert_eq!(raw_from_png(&opaque, Channel::Auto, None, 1).unwrap().3, Channel::Brightness);
    }

    /// The box downsample: each output pixel is the mean of its block (partial at the edges).
    #[test]
    fn streaming_decode_box_downsamples() {
        let bytes = rgba_png(37, 23, |_, _| 255);
        let rect = [2, 1, 33, 21];
        let (full, fw, fh, _) = raw_from_png(&bytes, Channel::Brightness, Some(rect), 1).unwrap();
        let step = 4usize;
        let (d, w, h, _) = raw_from_png(&bytes, Channel::Brightness, Some(rect), step as u32).unwrap();
        assert_eq!((w, h), (fw.div_ceil(4), fh.div_ceil(4)));
        for oy in 0..h as usize {
            for ox in 0..w as usize {
                let (mut sum, mut n) = (0.0f64, 0);
                for y in oy * step..((oy + 1) * step).min(fh as usize) {
                    for x in ox * step..((ox + 1) * step).min(fw as usize) {
                        sum += full[y * fw as usize + x] as f64;
                        n += 1;
                    }
                }
                let want = (sum / n as f64) as f32;
                assert!((d[oy * w as usize + ox] - want).abs() < 1e-6);
            }
        }
        assert_eq!(png_size(&bytes).unwrap(), (37, 23));
    }

    /// An imported world saves its mask source and placement, loads back to the same
    /// conditioning, and re-saves byte-identically; the drawn coast is where the mask says.
    #[test]
    fn imported_world_round_trips() {
        let svg = br#"<svg viewBox="0 0 100 100"><rect x="20" y="20" width="60" height="60"/></svg>"#.to_vec();
        let place = Placement::Patch { center_deg: [20.0, 30.0], km_per_px: 40.0, rotation_deg: 15.0 };
        let mask = mask_from_bytes(svg, "svg", place).unwrap();
        let mut params = crate::WorldParams::default();
        params.tectonics.res = 32;
        params.erosion.enabled = false;
        params.tectonics.conditioning.masks = vec![Arc::new(mask)];
        params.tectonics.conditioning.strict = true;
        params.tectonics.conditioning.outside_procedural = false;
        params.tectonics.conditioning.strokes.push(Stroke::new(vec![rg_tectonics::mask::dir_of(-40.0, 100.0)], 500.0, StrokeOp::Land));
        // s06: a stroke from another seed and a disabled one round-trip too (and don't apply).
        let mut other = Stroke::new(vec![rg_tectonics::mask::dir_of(-20.0, -150.0)], 800.0, StrokeOp::Land);
        other.seed = Some(params.tectonics.seed + 1);
        params.tectonics.conditioning.strokes.push(other);
        let mut off = Stroke::new(vec![rg_tectonics::mask::dir_of(-20.0, -150.0)], 800.0, StrokeOp::Land);
        off.enabled = false;
        params.tectonics.conditioning.strokes.push(off);
        let w = crate::generate(&params);
        // Patch centre is land (inside the drawn square); the antipode is sea.
        let at = |lat: f64, lon: f64| w.tectonics.height_m.sample(rg_tectonics::mask::dir_of(lat, lon));
        assert!(at(20.0, 30.0) > 0.0);
        assert!(at(-20.0, -150.0) < 0.0);
        assert!(at(-40.0, 100.0) > 0.0, "painted land");
        let mut a = Vec::new();
        crate::file::to_file(&w).write_to(&mut a).unwrap();
        let loaded = crate::file::from_file(&WorldFile::read_from(&a[..]).unwrap()).unwrap();
        assert_eq!(loaded.params.tectonics.conditioning, w.params.tectonics.conditioning);
        let mut b = Vec::new();
        crate::file::to_file(&loaded).write_to(&mut b).unwrap();
        assert!(a == b, "re-saved file differs");
    }

    /// A greyscale PNG with a threshold, invert and crop: the options survive a save/load
    /// and re-save byte-identically, and the land test follows them.
    #[test]
    fn adjusted_png_round_trips() {
        // Dark-grey land (0.25) on light sea (0.75), left half land.
        let bytes = {
            let mut out = Vec::new();
            let mut enc = png::Encoder::new(&mut out, 40, 20);
            enc.set_color(png::ColorType::Grayscale);
            enc.set_depth(png::BitDepth::Eight);
            let mut wr = enc.write_header().unwrap();
            let data: Vec<u8> = (0..20).flat_map(|_| (0..40u32).map(|x| if x < 20 { 64u8 } else { 191 })).collect();
            wr.write_image_data(&data).unwrap();
            drop(wr);
            out
        };
        let decode = DecodeOpts { channel: Channel::Brightness, crop: Some([4, 2, 32, 16]) };
        let adjust = MaskAdjust { threshold: 0.4, invert: true };
        let place = Placement::Patch { center_deg: [0.0, 0.0], km_per_px: 100.0, rotation_deg: 0.0 };
        let mask = mask_from_bytes_with(bytes, "png", place, decode, adjust).unwrap();
        assert_eq!((mask.w, mask.h), (32, 16));
        // Inverted, dark grey is land; the left of the crop (x < 16) is land.
        let r = 6.371e6;
        assert!(mask.coverage(mask.dir_of_pixel(4.0, 8.0, r), r).unwrap() > 0.5);
        assert!(mask.coverage(mask.dir_of_pixel(28.0, 8.0, r), r).unwrap() < 0.5);
        let mut params = crate::WorldParams::default();
        params.tectonics.res = 32;
        params.erosion.enabled = false;
        params.tectonics.conditioning.masks = vec![Arc::new(mask)];
        params.tectonics.conditioning.strict = true;
        let w = crate::generate(&params);
        let mut a = Vec::new();
        crate::file::to_file(&w).write_to(&mut a).unwrap();
        let loaded = crate::file::from_file(&WorldFile::read_from(&a[..]).unwrap()).unwrap();
        assert_eq!(loaded.params.tectonics.conditioning, w.params.tectonics.conditioning);
        let mut b = Vec::new();
        crate::file::to_file(&loaded).write_to(&mut b).unwrap();
        assert!(a == b, "re-saved file differs");
    }

    /// s07: a downsampled heightmap layer saves its elevation and step, loads back equal, and
    /// keeps its placement in source pixels (Move patch can't jump).
    #[test]
    fn heightmap_layer_round_trips() {
        // Greys: a sea ramp 0..~0.2 on the left, land ramp to white on the right.
        let bytes = {
            let mut out = Vec::new();
            let mut enc = png::Encoder::new(&mut out, 64, 32);
            enc.set_color(png::ColorType::Grayscale);
            enc.set_depth(png::BitDepth::Eight);
            let mut wr = enc.write_header().unwrap();
            let data: Vec<u8> = (0..32).flat_map(|_| (0..64u32).map(|x| (x * 4) as u8)).collect();
            wr.write_image_data(&data).unwrap();
            drop(wr);
            out
        };
        let place = Placement::Patch { center_deg: [10.0, 5.0], km_per_px: 20.0, rotation_deg: 10.0 };
        let adjust = MaskAdjust { threshold: 0.235, invert: false };
        let e = Elevation { max_m: 3000.0, gamma: 1.3, smooth_px: 3.0, fidelity: Fidelity::Exact, ..Default::default() };
        let m = layer_from_bytes(bytes.clone(), "png", place.clone(), DecodeOpts::default(), adjust, Some(e), Some(4)).unwrap();
        assert_eq!((m.w, m.h, m.step), (16, 8, 4));
        assert_eq!(user_placement(&m), place);
        assert!(m.heights.is_some());
        assert!(has_bathymetry(&m), "a graded sea counts as bathymetry");
        let flat = layer_from_bytes(test_png(32, 16, |x, _| x > 20), "png", place.clone(), DecodeOpts::default(), adjust, Some(e), None).unwrap();
        assert!(!has_bathymetry(&flat), "a flat sea doesn't");
        let mut moved = flat.clone();
        moved.adjust.threshold = 0.6;
        assert!(!has_bathymetry(&moved), "nor after the sea-level grey moves (land is grey 255)");
        let mut params = crate::WorldParams::default();
        params.tectonics.res = 32;
        params.erosion.enabled = false;
        params.tectonics.conditioning.masks = vec![Arc::new(m)];
        params.tectonics.conditioning.strict = true;
        let w = crate::generate(&params);
        let mut a = Vec::new();
        crate::file::to_file(&w).write_to(&mut a).unwrap();
        let loaded = crate::file::from_file(&WorldFile::read_from(&a[..]).unwrap()).unwrap();
        assert_eq!(loaded.params.tectonics.conditioning, w.params.tectonics.conditioning);
        let mut b = Vec::new();
        crate::file::to_file(&loaded).write_to(&mut b).unwrap();
        assert!(a == b, "re-saved file differs");
    }

    /// A greyscale dome island with radial ridges, sea grey 0.1 (no bathymetry).
    fn dome_png(side: u32) -> Vec<u8> {
        let mut out = Vec::new();
        {
            let mut enc = png::Encoder::new(&mut out, side, side);
            enc.set_color(png::ColorType::Grayscale);
            enc.set_depth(png::BitDepth::Eight);
            let mut wr = enc.write_header().unwrap();
            let c = side as f32 / 2.0;
            let data: Vec<u8> = (0..side * side)
                .map(|k| {
                    let (x, y) = ((k % side) as f32 - c, (k / side) as f32 - c);
                    let r = x.hypot(y) / (0.42 * side as f32);
                    if r >= 1.0 {
                        return 25;
                    }
                    let ridges = 0.08 * (6.0 * y.atan2(x)).cos() * r;
                    ((0.3 + 0.65 * (1.0 - r) + ridges).clamp(0.0, 1.0) * 255.0) as u8
                })
                .collect();
            wr.write_image_data(&data).unwrap();
        }
        out
    }

    /// s07 gate (Exact, blank ocean, so no free land at all): after erosion, lakes and
    /// calibration, pinned land keeps the imported heights (within 1% of max height outside
    /// lakes), and rivers still drain it.
    #[test]
    fn exact_heightmap_survives_erosion() {
        let place = Placement::Patch { center_deg: [5.0, 20.0], km_per_px: 50.0, rotation_deg: 0.0 };
        let adjust = MaskAdjust { threshold: 0.2, invert: false };
        let e = Elevation { max_m: 5000.0, fidelity: Fidelity::Exact, ..Default::default() };
        let m = layer_from_bytes(dome_png(96), "png", place, DecodeOpts::default(), adjust, Some(e), None).unwrap();
        let mut params = crate::WorldParams::default();
        params.tectonics.res = 64;
        params.erosion.canonical_res = 0; // erode at this resolution
        let c = &mut params.tectonics.conditioning;
        c.strict = true;
        c.outside_procedural = false;
        c.islands = false;
        c.masks = vec![Arc::new(m)];
        let snap = rg_tectonics::generate(&params.tectonics);
        let world = crate::generate(&params);
        let hy = world.hydrology.as_ref().expect("eroded");
        let pin = snap.pin.as_ref().unwrap();
        let (mut n, mut bad, mut worst, mut top) = (0, 0, 0.0f32, 0.0f32);
        for k in 0..snap.height_m.data.len() {
            let t = snap.height_m.data[k];
            if t <= 0.0 || pin.data[k] < 1.0 || hy.lake_depth_m.data[k] > 0.0 {
                continue;
            }
            let d = (world.tectonics.height_m.data[k] - t).abs();
            n += 1;
            worst = worst.max(d);
            top = top.max(world.tectonics.height_m.data[k]);
            if d > 0.01 * 5000.0 {
                bad += 1;
            }
        }
        assert!(n > 50, "only {n} pinned land samples");
        assert!(bad * 50 <= n, "{bad}/{n} samples off by > 1% (worst {worst} m)");
        assert!(top > 4000.0, "heights not rescaled: top {top}");
        assert!(!hy.rivers.is_empty(), "rivers drain the import");
        // Guided on the same map must move some pinned land (> 1%): the solve really runs
        // when all land is pinned (s07: it was skipped, and the Exact checks passed vacuously).
        let mut gm = (*params.tectonics.conditioning.masks[0]).clone();
        gm.elevation.as_mut().unwrap().fidelity = Fidelity::Guided;
        params.tectonics.conditioning.masks = vec![Arc::new(gm)];
        let guided = crate::generate(&params);
        let moved = (0..snap.height_m.data.len())
            .filter(|&k| snap.height_m.data[k] > 0.0 && pin.data[k] >= 1.0)
            .filter(|&k| (guided.tectonics.height_m.data[k] - snap.height_m.data[k]).abs() > 0.01 * 5000.0)
            .count();
        assert!(moved * 10 > n, "Guided moved only {moved}/{n} pinned samples");
    }

    /// A world saved with the s05 single-import layout (fields at the top of the block,
    /// chunk `conditioning.source`) still loads, as one import layer.
    #[test]
    fn legacy_single_import_file_loads() {
        let svg = br#"<svg viewBox="0 0 100 100"><rect x="20" y="20" width="60" height="60"/></svg>"#.to_vec();
        let place = Placement::Patch { center_deg: [20.0, 30.0], km_per_px: 40.0, rotation_deg: 15.0 };
        let mask = mask_from_bytes(svg, "svg", place).unwrap();
        let mut params = crate::WorldParams::default();
        params.tectonics.res = 32;
        params.erosion.enabled = false;
        params.tectonics.conditioning.masks = vec![Arc::new(mask)];
        params.tectonics.conditioning.strict = true;
        let w = crate::generate(&params);
        let mut f = crate::file::to_file(&w);
        // Rewrite to the old layout.
        let cond = &mut f.header["conditioning"];
        let entry = cond["masks"][0].clone();
        cond.as_object_mut().unwrap().remove("masks");
        for k in ["placement", "source_kind"] {
            cond[k] = entry[k].clone();
        }
        for c in f.chunks.iter_mut().filter(|c| c.name == "conditioning.source.0") {
            c.name = "conditioning.source".into();
        }
        let loaded = crate::file::from_file(&f).unwrap();
        assert_eq!(loaded.params.tectonics.conditioning, w.params.tectonics.conditioning);
    }
}
