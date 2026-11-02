//! Exports for other tools (Blender, GIS).

use crate::layers::{dir, lat_of_row, lon_of_col};
use crate::World;

/// A 16-bit greyscale equirectangular heightmap plus the metadata needed to turn grey
/// levels back into metres.
pub struct HeightmapExport {
    pub png: Vec<u8>,
    /// JSON sidecar: min/max metres, radius, projection.
    pub sidecar: String,
}

/// Encode heights linearly: grey 0 = `min_m`, grey 65535 = `max_m`. Sea level sits at
/// `-min_m / (max_m - min_m)` of the range (written to the sidecar as `sea_level_grey`).
pub fn heightmap_png(world: &World, width: u32) -> Result<HeightmapExport, String> {
    let height = width / 2;
    let field = &world.tectonics.height_m;
    let (min_m, max_m) = field.min_max();
    let span = (max_m - min_m).max(1.0);
    let mut px = Vec::with_capacity((width * height * 2) as usize);
    for y in 0..height {
        let lat = lat_of_row(y, height);
        for x in 0..width {
            let h = field.sample(dir(lat, lon_of_col(x, width)));
            let g = (((h - min_m) / span).clamp(0.0, 1.0) * 65535.0).round() as u16;
            px.extend_from_slice(&g.to_be_bytes()); // PNG is big-endian
        }
    }
    let mut png = Vec::new();
    {
        let mut enc = png::Encoder::new(&mut png, width, height);
        enc.set_color(png::ColorType::Grayscale);
        enc.set_depth(png::BitDepth::Sixteen);
        let mut w = enc.write_header().map_err(|e| e.to_string())?;
        w.write_image_data(&px).map_err(|e| e.to_string())?;
    }
    let sidecar = serde_json::to_string_pretty(&serde_json::json!({
        "projection": "equirectangular (plate carrée), lon -180..180 left to right, lat 90..-90 top to bottom",
        "width": width, "height": height,
        "min_m": min_m, "max_m": max_m,
        "sea_level_grey": ((-min_m / span) * 65535.0).round(),
        "radius_m": world.params.tectonics.radius_m,
        "seed": world.params.tectonics.seed,
        "blender_hint": "Displace modifier: strength = (max_m - min_m) / radius_m in sphere units; midlevel = sea_level_grey / 65535",
    }))
    .map_err(|e| e.to_string())?;
    Ok(HeightmapExport { png, sidecar })
}

#[cfg(test)]
mod tests {
    #[test]
    fn heightmap_is_16_bit_and_spans_range() {
        let mut p = crate::WorldParams::default();
        p.tectonics.res = 32;
        let w = crate::generate(&p);
        let ex = super::heightmap_png(&w, 256).unwrap();
        let dec = png::Decoder::new(std::io::Cursor::new(&ex.png));
        let mut r = dec.read_info().unwrap();
        let info = r.info();
        assert_eq!((info.width, info.height), (256, 128));
        assert_eq!(info.bit_depth, png::BitDepth::Sixteen);
        let mut buf = vec![0; r.output_buffer_size().unwrap()];
        r.next_frame(&mut buf).unwrap();
        let vals: Vec<u16> = buf.as_chunks::<2>().0.iter().map(|b| u16::from_be_bytes(*b)).collect();
        let (lo, hi) = (*vals.iter().min().unwrap(), *vals.iter().max().unwrap());
        assert!(lo < 2000 && hi > 60000, "{lo}..{hi}");
        assert!(ex.sidecar.contains("sea_level_grey"));
    }
}
