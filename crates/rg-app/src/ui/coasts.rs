//! R3 coast controls: where the coasts come from (procedural, procedural strict, or an
//! imported PNG/SVG mask placed on the sphere), and the brush (R3.3): paint land, sea, or
//! force a coast active or passive, then regenerate with the same seed.
//!
//! s06 (Matt): an imported PNG gets a quick adjust panel (preview, threshold, invert,
//! channel, crop), and a regional patch can be dragged across the globe or map.

use crate::app::RealmApp;
use glam::{DQuat, DVec3};
use rg_pipeline::import;
use rg_pipeline::mask::{Channel, Conditioning, DecodeOpts, MaskAdjust, Placement, RasterMask, Stroke, StrokeOp};
use std::sync::Arc;

/// A small copy of the whole (uncropped) imported image, raw values, for the preview.
struct Thumb {
    raw: Vec<f32>,
    w: usize,
    h: usize,
    /// The source image's size in pixels (crop rectangles are in these).
    src_w: u32,
    src_h: u32,
}

/// What the preview texture was built from, to know when to rebuild it.
type TexKey = (u32, bool, Option<[u32; 4]>, Channel);

/// UI state for the coast section (not part of the world params).
pub struct CoastUi {
    pub file_name: String,
    pub patch: bool,
    pub lat: f64,
    pub lon: f64,
    pub km_per_px: f64,
    pub rotation: f64,
    pub msg: String,
    /// The active brush, if painting.
    pub brush: Option<StrokeOp>,
    pub brush_km: f64,
    /// The stroke being drawn (directions on the unit sphere).
    pub drawing: Vec<DVec3>,
    /// "Move patch": a primary drag moves the imported patch (never with a brush).
    pub placing: bool,
    drag_prev: Option<DVec3>,
    /// Draw the imported mask over the globe or map.
    pub show_mask: bool,
    thumb: Option<Thumb>,
    tex: Option<(egui::TextureHandle, TexKey)>,
    /// A crop being dragged on the preview, in source pixels (start, current).
    crop_drag: Option<([f32; 2], [f32; 2])>,
}

impl Default for CoastUi {
    fn default() -> Self {
        Self {
            file_name: String::new(),
            patch: false,
            lat: 0.0,
            lon: 0.0,
            km_per_px: 10.0,
            rotation: 0.0,
            msg: String::new(),
            brush: None,
            brush_km: 150.0,
            drawing: Vec::new(),
            placing: false,
            drag_prev: None,
            show_mask: true,
            thumb: None,
            tex: None,
            crop_drag: None,
        }
    }
}

impl CoastUi {
    fn placement(&self) -> Placement {
        if self.patch {
            Placement::Patch { center_deg: [self.lat, self.lon], km_per_px: self.km_per_px, rotation_deg: self.rotation }
        } else {
            Placement::WorldEquirect
        }
    }
}

const SOURCES: [&str; 3] = ["Procedural (coast from heights)", "Procedural, strict coasts", "Imported mask (PNG / SVG)"];

/// (lat, lon) in degrees at the centre of the current view, where a new patch goes.
pub fn view_centre_deg(app: &RealmApp) -> (f64, f64) {
    match app.view {
        rg_render::ViewMode::Globe => (app.camera.lat.to_degrees(), app.camera.lon.to_degrees()),
        rg_render::ViewMode::Map => (app.map.center_lat.to_degrees(), app.map.center_lon.to_degrees()),
    }
}

pub fn coasts_section(app: &mut RealmApp, ui: &mut egui::Ui) {
    let centre = view_centre_deg(app);
    let c = &mut app.params.tectonics.conditioning;
    let st = &mut app.coasts;
    let mut source = if c.mask.is_some() { 2 } else if c.strict { 1 } else { 0 };
    let before = source;
    egui::ComboBox::from_label("coasts")
        .selected_text(SOURCES[source])
        .show_ui(ui, |ui| {
            for (i, s) in SOURCES.iter().enumerate() {
                ui.selectable_value(&mut source, i, *s);
            }
        })
        .response
        .on_hover_text("Procedural: the coast emerges from the tectonic heights (land % sets the sea level). Strict: the coastline is drawn first and the terrain conforms to it. Imported: your own coastline map.");
    if source != before {
        match source {
            0 => {
                c.strict = false;
                c.mask = None;
            }
            1 => {
                c.strict = true;
                c.mask = None;
            }
            _ => c.strict = true,
        }
    }
    // Import is always one click away (Matt, s06: the combo box hid it).
    ui.horizontal(|ui| {
        #[cfg(not(target_arch = "wasm32"))]
        if ui
            .button(if c.mask.is_some() { "Import another…" } else { "Import coastline…" })
            .on_hover_text("Open a .png (land brighter than sea, or opaque land on a transparent sea; greyscale is fine) or an .svg (every filled shape is land). You can also drop the file on the window.")
            .clicked()
            && let Some(path) = rfd::FileDialog::new().add_filter("Coastline map", &["png", "svg"]).pick_file()
        {
            import_file(c, st, &path, centre);
            source = 2;
        }
        if source != 2 {
            ui.label(egui::RichText::new("or drop a .png / .svg here").small().color(super::theme::MUTED));
        }
    });
    if source == 2 {
        mask_panel(c, st, ui);
    }
    if c.strict {
        ui.checkbox(&mut c.islands, "tectonic islands offshore")
            .on_hover_text("Island arcs and hotspot volcanoes may rise above the sea farther than 100 km from the drawn coasts. Drawn coasts never move.");
    }

    // ---- brush (R3.3)
    ui.horizontal(|ui| {
        ui.label("brush");
        for (op, name, tip) in [
            (StrokeOp::Land, "land", "Paint land (switches to strict coasts)."),
            (StrokeOp::Sea, "sea", "Paint sea (switches to strict coasts)."),
            (StrokeOp::Active, "active", "Force a coast active: trench and mountains."),
            (StrokeOp::Passive, "passive", "Force a coast passive: wide shelf, no trench."),
        ] {
            let on = st.brush == Some(op);
            if ui.selectable_label(on, name).on_hover_text(tip).clicked() {
                st.brush = if on { None } else { Some(op) };
                st.placing = false;
            }
        }
    });
    if st.brush.is_some() {
        ui.add(egui::Slider::new(&mut st.brush_km, 20.0..=1500.0).logarithmic(true).text("brush km"));
        ui.label(egui::RichText::new("Drag on the globe to paint; right-drag still orbits. Generate to apply.").small().color(super::theme::MUTED));
    }
    if !c.strokes.is_empty() {
        ui.horizontal(|ui| {
            ui.label(egui::RichText::new(format!("{} strokes", c.strokes.len())).small());
            if ui.small_button("Undo").clicked() {
                c.strokes.pop();
            }
            if ui.small_button("Clear").clicked() {
                c.strokes.clear();
            }
        });
    }
}

/// The imported mask: file, placement, and (PNG) the adjust panel.
fn mask_panel(c: &mut Conditioning, st: &mut CoastUi, ui: &mut egui::Ui) {
    let Some(m) = c.mask.clone() else {
        ui.label(egui::RichText::new("no mask loaded yet").small().color(super::theme::MUTED));
        return;
    };
    let png = m.source_kind == "png";
    let channel = match m.channel_used {
        Channel::Alpha => "alpha",
        _ => "brightness",
    };
    let info = if png { format!("PNG {}×{}, land from {channel}", m.w, m.h) } else { format!("SVG, raster {}×{}", m.w, m.h) };
    ui.label(egui::RichText::new(format!("{} · {info}", st.file_name)).small());

    // Placement.
    let was_patch = st.patch;
    ui.horizontal(|ui| {
        ui.radio_value(&mut st.patch, false, "whole world").on_hover_text("A 2:1 equirectangular map of the whole planet.");
        ui.radio_value(&mut st.patch, true, "regional patch").on_hover_text("A map of one region, placed by its centre, scale and rotation.");
    });
    if st.patch != was_patch {
        // An SVG's raster depends on the placement kind, so re-decode.
        reload(c, st, m.decode.clone(), m.adjust);
        st.placing &= st.patch;
    }
    if st.patch {
        let mut changed = false;
        ui.horizontal(|ui| {
            changed |= ui.add(egui::DragValue::new(&mut st.lat).range(-90.0..=90.0).speed(0.2).prefix("lat ")).changed();
            changed |= ui.add(egui::DragValue::new(&mut st.lon).range(-180.0..=180.0).speed(0.2).prefix("lon ")).changed();
        });
        ui.horizontal(|ui| {
            changed |= ui
                .add(egui::DragValue::new(&mut st.km_per_px).range(0.01..=500.0).speed(0.05).prefix("km/px "))
                .on_hover_text("PNG: km per pixel. SVG: km per drawing unit.")
                .changed();
            changed |= ui.add(egui::DragValue::new(&mut st.rotation).range(-180.0..=180.0).speed(0.5).prefix("rot° ")).changed();
        });
        if changed {
            apply_placement(c, st);
        }
        ui.horizontal(|ui| {
            if ui
                .selectable_label(st.placing, "✥ Move patch")
                .on_hover_text("Drag the patch on the globe or map. ⌥ + scroll scales it. Right-drag still orbits.")
                .clicked()
            {
                st.placing = !st.placing;
                if st.placing {
                    st.brush = None;
                    st.show_mask = true;
                }
            }
            let side = st.km_per_px * m.w.max(m.h) as f64;
            ui.label(egui::RichText::new(format!("{side:.0} km across")).small().color(super::theme::MUTED));
        });
        if st.placing {
            ui.label(egui::RichText::new("Drag on the view to move · ⌥ scroll to scale").small().color(super::theme::MUTED));
        }
    }

    // Adjust (PNG only: an SVG's fill is already a clean land/sea split).
    if png {
        adjust_panel(c, st, ui, &m);
    }
    ui.checkbox(&mut st.show_mask, "show mask on the view").on_hover_text("Tint the imported land over the globe or map (before you Generate).");
    if st.patch {
        ui.checkbox(&mut c.outside_procedural, "procedural land outside the patch")
            .on_hover_text("Other continents around the imported region, kept 800 km away. Off: open ocean.");
    }
    if !st.msg.is_empty() {
        ui.label(egui::RichText::new(&st.msg).small().color(super::theme::MUTED));
    }
}

/// Preview with crop, channel, threshold and invert for a PNG mask.
fn adjust_panel(c: &mut Conditioning, st: &mut CoastUi, ui: &mut egui::Ui, m: &Arc<RasterMask>) {
    if st.thumb.is_none() {
        st.thumb = make_thumb(m);
    }
    let Some(th) = st.thumb.as_ref() else { return };
    let (src_w, src_h) = (th.src_w, th.src_h);
    // Preview texture: raw grey with the land tinted green, rebuilt when the land test changes.
    let key: TexKey = (m.adjust.threshold.to_bits(), m.adjust.invert, m.decode.crop, m.decode.channel);
    if st.tex.as_ref().is_none_or(|(_, k)| *k != key) {
        let img = thumb_image(th, &m.adjust);
        let handle = ui.ctx().load_texture("coast-mask-preview", img, egui::TextureOptions::LINEAR);
        st.tex = Some((handle, key));
    }
    let tex = st.tex.as_ref().map(|(t, _)| t.clone()).unwrap();
    let avail = ui.available_width().min(300.0);
    let size = egui::vec2(avail, avail * src_h as f32 / src_w.max(1) as f32).min(egui::vec2(avail, 220.0));
    let size = egui::vec2(size.y * src_w as f32 / src_h.max(1) as f32, size.y);
    let resp = ui
        .add(egui::Image::new(&tex).fit_to_exact_size(size).sense(egui::Sense::drag()))
        .on_hover_text("Green = land. Drag a box to crop.");
    let to_src = |p: egui::Pos2| -> [f32; 2] {
        let f = (p - resp.rect.min) / resp.rect.size();
        [(f.x.clamp(0.0, 1.0) * src_w as f32), (f.y.clamp(0.0, 1.0) * src_h as f32)]
    };
    let to_ui = |q: [f32; 2]| resp.rect.min + egui::vec2(q[0] / src_w as f32, q[1] / src_h as f32) * resp.rect.size();
    if resp.drag_started()
        && let Some(p) = resp.interact_pointer_pos()
    {
        st.crop_drag = Some((to_src(p), to_src(p)));
    }
    if let (Some((_, cur)), Some(p)) = (st.crop_drag.as_mut(), resp.interact_pointer_pos()) {
        *cur = to_src(p);
    }
    let painter = ui.painter_at(resp.rect);
    // Current crop (or the one being dragged) as a box; outside it dimmed.
    let crop_px = st.crop_drag.map(|(a, b)| [a[0].min(b[0]), a[1].min(b[1]), a[0].max(b[0]), a[1].max(b[1])]).or_else(|| {
        m.decode.crop.map(|[x, y, w, h]| [x as f32, y as f32, (x + w) as f32, (y + h) as f32])
    });
    if let Some([x0, y0, x1, y1]) = crop_px {
        let r = egui::Rect::from_min_max(to_ui([x0, y0]), to_ui([x1, y1]));
        let dim = egui::Color32::from_black_alpha(140);
        for outside in [
            egui::Rect::from_min_max(resp.rect.min, egui::pos2(resp.rect.max.x, r.min.y)),
            egui::Rect::from_min_max(egui::pos2(resp.rect.min.x, r.max.y), resp.rect.max),
            egui::Rect::from_min_max(egui::pos2(resp.rect.min.x, r.min.y), egui::pos2(r.min.x, r.max.y)),
            egui::Rect::from_min_max(egui::pos2(r.max.x, r.min.y), egui::pos2(resp.rect.max.x, r.max.y)),
        ] {
            painter.rect_filled(outside, 0.0, dim);
        }
        painter.rect_stroke(r, 0.0, (1.5, egui::Color32::WHITE), egui::StrokeKind::Middle);
    }
    if resp.drag_stopped()
        && let Some((a, b)) = st.crop_drag.take()
    {
        let (x0, y0) = (a[0].min(b[0]).floor() as u32, a[1].min(b[1]).floor() as u32);
        let (x1, y1) = (a[0].max(b[0]).ceil() as u32, a[1].max(b[1]).ceil() as u32);
        if x1 > x0 + 4 && y1 > y0 + 4 {
            let decode = DecodeOpts { crop: Some([x0, y0, x1 - x0, y1 - y0]), ..m.decode.clone() };
            reload(c, st, decode, m.adjust);
        }
    }

    // Threshold and invert apply live (the raster is shared, so this is cheap).
    let mut adj = m.adjust;
    let mut pct = adj.threshold * 100.0;
    ui.horizontal(|ui| {
        if ui
            .add(egui::Slider::new(&mut pct, 1.0..=99.0).integer().suffix("%").text("land above"))
            .on_hover_text("The grey level that becomes the coast. Brighter than this is land (darker, with Invert).")
            .changed()
        {
            adj.threshold = pct / 100.0;
        }
    });
    ui.horizontal(|ui| {
        ui.checkbox(&mut adj.invert, "invert").on_hover_text("Dark land on a light sea.");
        let mut ch = m.decode.channel;
        egui::ComboBox::from_id_salt("mask-channel")
            .selected_text(match ch {
                Channel::Auto => "auto",
                Channel::Brightness => "brightness",
                Channel::Alpha => "alpha",
            })
            .width(90.0)
            .show_ui(ui, |ui| {
                ui.selectable_value(&mut ch, Channel::Auto, "auto").on_hover_text("Alpha if the image has transparency, else brightness.");
                ui.selectable_value(&mut ch, Channel::Brightness, "brightness");
                ui.selectable_value(&mut ch, Channel::Alpha, "alpha");
            });
        if ch != m.decode.channel {
            reload(c, st, DecodeOpts { channel: ch, ..m.decode.clone() }, adj);
            st.thumb = None;
        }
        if m.decode.crop.is_some() && ui.small_button("Reset crop").clicked() {
            reload(c, st, DecodeOpts { crop: None, ..m.decode.clone() }, adj);
        }
    });
    if adj != m.adjust
        && let Some(mm) = c.mask.as_mut()
    {
        Arc::make_mut(mm).adjust = adj;
    }
}

/// Downsample the whole source image (no crop) to at most 384 px, box-averaged so thin
/// coasts survive.
fn make_thumb(m: &RasterMask) -> Option<Thumb> {
    let (raw, sw, sh) = if m.source_kind == "png" {
        let (d, w, h, _) = import::coverage_from_png(&m.source, m.decode.channel).ok()?;
        (d, w, h)
    } else {
        ((*m.data).clone(), m.w, m.h)
    };
    let step = (sw.max(sh) as usize).div_ceil(384).max(1);
    let (w, h) = ((sw as usize).div_ceil(step), (sh as usize).div_ceil(step));
    let mut out = vec![0.0f32; w * h];
    for ty in 0..h {
        for tx in 0..w {
            let (mut sum, mut n) = (0.0, 0);
            for y in ty * step..((ty + 1) * step).min(sh as usize) {
                for x in tx * step..((tx + 1) * step).min(sw as usize) {
                    sum += raw[y * sw as usize + x];
                    n += 1;
                }
            }
            out[ty * w + tx] = sum / n.max(1) as f32;
        }
    }
    Some(Thumb { raw: out, w, h, src_w: sw, src_h: sh })
}

/// The preview: the raw grey, with land (by the same test Generate uses) tinted green.
fn thumb_image(th: &Thumb, adj: &MaskAdjust) -> egui::ColorImage {
    let mut rgba = Vec::with_capacity(th.w * th.h * 4);
    for &v in &th.raw {
        let g = (v.clamp(0.0, 1.0) * 255.0) as u8;
        let px = if adj.apply(v) > 0.5 {
            [(g as u16 * 2 / 5 + 50) as u8, (g as u16 * 2 / 5 + 140) as u8, (g as u16 * 2 / 5 + 50) as u8, 255]
        } else {
            [g / 2, g / 2, (g / 2).saturating_add(40), 255]
        };
        rgba.extend_from_slice(&px);
    }
    egui::ColorImage::from_rgba_unmultiplied([th.w, th.h], &rgba)
}

/// Re-decode the current mask's file with new decode options (crop, channel) and the
/// panel's placement.
fn reload(c: &mut Conditioning, st: &mut CoastUi, decode: DecodeOpts, adjust: MaskAdjust) {
    let Some(m) = c.mask.as_ref() else { return };
    match import::mask_from_bytes_with((*m.source).clone(), &m.source_kind.clone(), st.placement(), decode, adjust) {
        Ok(new) => {
            c.mask = Some(Arc::new(new));
            st.msg = "Generate to apply".into();
        }
        Err(e) => st.msg = e,
    }
}

/// Write the panel's placement to the mask (cheap: the raster is shared).
fn apply_placement(c: &mut Conditioning, st: &CoastUi) {
    if let Some(m) = c.mask.as_mut() {
        let p = import::raster_placement(m, st.placement());
        Arc::make_mut(m).place = p;
    }
}

/// Viewport hook: in "Move patch" mode a primary drag moves the patch so the grabbed point
/// follows the pointer (a rotation of the sphere, so it behaves at the poles too).
pub fn place_drag(app: &mut RealmApp, dragging_primary: bool, released: bool, at: Option<(f64, f64)>) -> bool {
    if !app.coasts.placing || !app.coasts.patch || app.params.tectonics.conditioning.mask.is_none() {
        return false;
    }
    if dragging_primary && let Some((lat, lon)) = at {
        let d = rg_render::earth_camera::dir(lat, lon);
        if let Some(prev) = app.coasts.drag_prev {
            let q = DQuat::from_rotation_arc(prev, d);
            let st = &mut app.coasts;
            let c = q * rg_pipeline::mask::dir_of(st.lat, st.lon);
            st.lat = c.z.clamp(-1.0, 1.0).asin().to_degrees();
            st.lon = c.y.atan2(c.x).to_degrees();
            apply_placement(&mut app.params.tectonics.conditioning, &app.coasts);
        }
        app.coasts.drag_prev = Some(d);
    }
    if released {
        app.coasts.drag_prev = None;
        app.coasts.msg = "Generate to apply".into();
    }
    dragging_primary || released
}

/// ⌥ + scroll in "Move patch" mode scales the patch. Returns true when it used the scroll.
pub fn place_scroll(app: &mut RealmApp, scroll: f32, alt: bool) -> bool {
    if !(alt && app.coasts.placing && app.coasts.patch && scroll != 0.0) {
        return false;
    }
    app.coasts.km_per_px = (app.coasts.km_per_px * (-scroll as f64 * 0.0025).exp()).clamp(0.01, 500.0);
    apply_placement(&mut app.params.tectonics.conditioning, &app.coasts);
    true
}

/// Viewport hook: with a brush active, a primary drag paints instead of panning. Returns
/// true when it consumed the drag.
pub fn paint(app: &mut RealmApp, dragging_primary: bool, released: bool, at: Option<(f64, f64)>) -> bool {
    let Some(op) = app.coasts.brush else { return false };
    if dragging_primary && let Some((lat, lon)) = at {
        let d = DVec3::new(lat.cos() * lon.cos(), lat.cos() * lon.sin(), lat.sin());
        // Keep points ~a fifth of the brush apart.
        let min_rad = app.coasts.brush_km * 0.2 * 1000.0 / app.params.tectonics.radius_m;
        if app.coasts.drawing.last().is_none_or(|q| q.dot(d).clamp(-1.0, 1.0).acos() > min_rad) {
            app.coasts.drawing.push(d);
        }
    }
    if released && !app.coasts.drawing.is_empty() {
        let c = &mut app.params.tectonics.conditioning;
        if matches!(op, StrokeOp::Land | StrokeOp::Sea) {
            c.strict = true;
        }
        c.strokes.push(Stroke { path: std::mem::take(&mut app.coasts.drawing), radius_km: app.coasts.brush_km, op });
        app.status = format!("{} strokes; Generate to apply", c.strokes.len());
    }
    dragging_primary || released
}

/// A picked or dropped file: a 2:1 image is taken as a whole-world map, anything else as a
/// ~3,000 km regional patch at the centre of the view (drag it into place).
fn import_file(c: &mut Conditioning, st: &mut CoastUi, path: &std::path::Path, centre: (f64, f64)) {
    st.file_name = path.file_name().map_or_else(|| path.display().to_string(), |n| n.to_string_lossy().into_owned());
    st.thumb = None;
    st.tex = None;
    st.crop_drag = None;
    st.patch = false;
    st.show_mask = true;
    match import::load_mask(path, Placement::WorldEquirect) {
        Ok(m) => {
            c.mask = Some(Arc::new(m));
            c.strict = true;
            st.msg = "Generate to apply".into();
        }
        Err(e) => {
            st.msg = e;
            return;
        }
    }
    let m = c.mask.clone().unwrap();
    let ratio = m.w as f64 / m.h.max(1) as f64;
    if (ratio - 2.0).abs() > 0.05 {
        st.patch = true;
        st.placing = true;
        st.brush = None;
        (st.lat, st.lon) = centre;
        st.rotation = 0.0;
        // SVG: km per drawing unit; aim for ~3,000 km across either way.
        let side_units = if m.source_kind == "svg" {
            rg_pipeline::svg::view_side(std::str::from_utf8(&m.source).unwrap_or("")).unwrap_or(m.w.max(m.h) as f64)
        } else {
            m.w.max(m.h) as f64
        };
        st.km_per_px = ((3000.0 / side_units) * 1000.0).round() / 1000.0;
        reload(c, st, DecodeOpts::default(), MaskAdjust::default());
        st.msg = "Not 2:1, so placed as a regional patch at the view centre: drag it into place, then Generate".into();
    }
}

/// Files dropped on the window: the first one is imported as the coast mask (a wrong file
/// type shows the loader's error). Also draws a hint while a file hovers over the window.
pub fn handle_drops(app: &mut RealmApp, ui: &egui::Ui, rect: egui::Rect) {
    let (hovering, dropped) = ui.ctx().input(|i| (!i.raw.hovered_files.is_empty(), i.raw.dropped_files.clone()));
    if hovering {
        let r = egui::Rect::from_center_size(rect.center(), egui::vec2(320.0, 48.0));
        ui.painter().rect_filled(r, 6.0, super::theme::PANEL.gamma_multiply(0.92));
        ui.painter().text(
            r.center(),
            egui::Align2::CENTER_CENTER,
            "Drop a .png or .svg to import a coastline",
            egui::FontId::proportional(15.0),
            super::theme::TEXT,
        );
    }
    let Some(path) = dropped.first().map(|f| f.path().to_path_buf()) else { return };
    let centre = view_centre_deg(app);
    import_file(&mut app.params.tectonics.conditioning, &mut app.coasts, &path, centre);
    app.status = app.coasts.msg.clone();
}

fn op_colour(op: StrokeOp) -> [u8; 3] {
    match op {
        StrokeOp::Land => [120, 200, 90],
        StrokeOp::Sea => [70, 130, 230],
        StrokeOp::Active => [230, 90, 60],
        StrokeOp::Passive => [60, 185, 205],
    }
}

/// Point at angular distance `a` from `c` in the direction `t` (unit, tangent at `c`).
fn offset(c: DVec3, t: DVec3, a: f64) -> DVec3 {
    c * a.cos() + t * a.sin()
}

/// Some unit vector perpendicular to `c`.
fn any_tangent(c: DVec3) -> DVec3 {
    let helper = if c.z.abs() < 0.9 { DVec3::Z } else { DVec3::X };
    c.cross(helper).normalize()
}

/// A stroke's footprint on the sphere: a band `a` radians wide on each side of the path,
/// with round caps, as triangles to fill and a ring to outline.
struct Band {
    /// Triangles on the sphere (three points each).
    tris: Vec<[DVec3; 3]>,
    /// The closed outline.
    ring: Vec<DVec3>,
}

fn band(path: &[DVec3], a: f64) -> Band {
    const CAP: usize = 12;
    let mut tris = Vec::new();
    if path.len() == 1 {
        let c = path[0];
        let e1 = any_tangent(c);
        let e2 = c.cross(e1);
        let ring: Vec<DVec3> = (0..CAP * 2)
            .map(|k| {
                let th = k as f64 / (CAP * 2) as f64 * std::f64::consts::TAU;
                offset(c, e1 * th.cos() + e2 * th.sin(), a)
            })
            .collect();
        for k in 0..ring.len() {
            tris.push([c, ring[k], ring[(k + 1) % ring.len()]]);
        }
        return Band { tris, ring };
    }
    // Per-point forward tangent and left normal.
    let n = path.len();
    let frames: Vec<(DVec3, DVec3)> = (0..n)
        .map(|i| {
            let d = path[i];
            let raw = path[(i + 1).min(n - 1)] - path[i.saturating_sub(1)];
            let t = raw - d * d.dot(raw);
            let t = if t.length_squared() > 1e-24 { t.normalize() } else { any_tangent(d) };
            (t, d.cross(t))
        })
        .collect();
    let left: Vec<DVec3> = (0..n).map(|i| offset(path[i], frames[i].1, a)).collect();
    let right: Vec<DVec3> = (0..n).map(|i| offset(path[i], -frames[i].1, a)).collect();
    for i in 0..n - 1 {
        tris.push([left[i], right[i], right[i + 1]]);
        tris.push([left[i], right[i + 1], left[i + 1]]);
    }
    // Half-disc caps, so they don't overlap the band (overlaps would read darker).
    let cap = |c: DVec3, t: DVec3, nrm: DVec3| -> Vec<DVec3> {
        (0..=CAP)
            .map(|k| {
                let th = k as f64 / CAP as f64 * std::f64::consts::PI;
                offset(c, nrm * th.cos() + t * th.sin(), a)
            })
            .collect()
    };
    let end = cap(path[n - 1], frames[n - 1].0, frames[n - 1].1); // left → forward → right
    let start = cap(path[0], -frames[0].0, -frames[0].1); // right → backward → left
    for w in end.windows(2) {
        tris.push([path[n - 1], w[0], w[1]]);
    }
    for w in start.windows(2) {
        tris.push([path[0], w[0], w[1]]);
    }
    let mut ring = left;
    ring.extend(end);
    ring.extend(right.into_iter().rev());
    ring.extend(start);
    Band { tris, ring }
}

/// Viewport overlay while a brush is active: every stroke (painted but maybe not generated
/// yet) as a translucent band, the stroke in progress, and a circle cursor at the brush
/// radius. `hover` is the (lat, lon) under the pointer.
pub fn overlay(app: &RealmApp, ui: &egui::Ui, rect: egui::Rect, aspect: f64, hover: Option<(f64, f64)>) {
    let Some(op) = app.coasts.brush else { return };
    let r_m = app.params.tectonics.radius_m;
    let painter = ui.painter().with_clip_rect(rect);
    let project = projector(app, rect, aspect);
    // On the map, an edge longer than a quarter of the view wraps around the antimeridian.
    let max_jump = rect.width() * 0.25;
    let ok = |p: egui::Pos2, q: egui::Pos2| (p.x - q.x).abs() < max_jump && (p.y - q.y).abs() < max_jump;

    let draw_band = |path: &[DVec3], radius_km: f64, op: StrokeOp, strong: bool| {
        if path.is_empty() {
            return;
        }
        let [cr, cg, cb] = op_colour(op);
        let fill = egui::Color32::from_rgba_unmultiplied(cr, cg, cb, if strong { 110 } else { 80 });
        let line = egui::Color32::from_rgba_unmultiplied(cr, cg, cb, 220);
        let b = band(path, radius_km * 1000.0 / r_m);
        let mut mesh = egui::Mesh::default();
        for t in &b.tris {
            let (Some(p0), Some(p1), Some(p2)) = (project(t[0]), project(t[1]), project(t[2])) else { continue };
            if !(ok(p0, p1) && ok(p1, p2) && ok(p0, p2)) {
                continue;
            }
            let i = mesh.vertices.len() as u32;
            for p in [p0, p1, p2] {
                mesh.colored_vertex(p, fill);
            }
            mesh.add_triangle(i, i + 1, i + 2);
        }
        painter.add(egui::Shape::mesh(mesh));
        let pts: Vec<Option<egui::Pos2>> = b.ring.iter().map(|&d| project(d)).collect();
        for k in 0..pts.len() {
            if let (Some(p), Some(q)) = (pts[k], pts[(k + 1) % pts.len()])
                && ok(p, q)
            {
                painter.line_segment([p, q], (1.2, line));
            }
        }
    };
    for st in &app.params.tectonics.conditioning.strokes {
        draw_band(&st.path, st.radius_km, st.op, false);
    }
    draw_band(&app.coasts.drawing, app.coasts.brush_km, op, true);

    // Circle cursor at the brush radius.
    if let Some((lat, lon)) = hover {
        let c = rg_render::earth_camera::dir(lat, lon);
        let a = app.coasts.brush_km * 1000.0 / r_m;
        let e1 = any_tangent(c);
        let e2 = c.cross(e1);
        let ring: Vec<Option<egui::Pos2>> = (0..64)
            .map(|k| {
                let th = k as f64 / 64.0 * std::f64::consts::TAU;
                project(offset(c, e1 * th.cos() + e2 * th.sin(), a))
            })
            .collect();
        let [cr, cg, cb] = op_colour(op);
        for k in 0..ring.len() {
            if let (Some(p), Some(q)) = (ring[k], ring[(k + 1) % ring.len()])
                && ok(p, q)
            {
                painter.line_segment([p, q], (3.0, egui::Color32::from_black_alpha(140)));
                painter.line_segment([p, q], (1.5, egui::Color32::from_rgb(cr, cg, cb)));
            }
        }
        if let Some(p) = project(c) {
            painter.circle_filled(p, 2.0, egui::Color32::WHITE);
        }
    }
}

/// Screen position of a direction in the current view (globe or map), `None` when hidden.
fn projector(app: &RealmApp, rect: egui::Rect, aspect: f64) -> impl Fn(DVec3) -> Option<egui::Pos2> + '_ {
    use rg_render::ViewMode;
    let to_screen = move |fx: f64, fy: f64| egui::pos2(rect.left() + fx as f32 * rect.width(), rect.top() + fy as f32 * rect.height());
    move |d: DVec3| {
        let lat = d.z.clamp(-1.0, 1.0).asin();
        let lon = d.y.atan2(d.x);
        match app.view {
            ViewMode::Globe => rg_render::graticule::project_globe(&app.camera, aspect, lat, lon).map(|[fx, fy]| to_screen(fx, fy)),
            ViewMode::Map => {
                use std::f64::consts::{PI, TAU};
                let m = &app.map;
                let dlon = (lon - m.center_lon + PI).rem_euclid(TAU) - PI;
                Some(to_screen(0.5 + dlon / m.width_rad, 0.5 - (lat - m.center_lat) / (m.width_rad / aspect)))
            }
        }
    }
}

/// The imported mask over the view, before Generate: land tinted green on a coarse grid
/// (the same land test Generate uses), and a patch's outline.
pub fn mask_overlay(app: &RealmApp, ui: &egui::Ui, rect: egui::Rect, aspect: f64) {
    let Some(m) = app.params.tectonics.conditioning.mask.as_ref() else { return };
    if !app.coasts.show_mask {
        return;
    }
    let r_m = app.params.tectonics.radius_m;
    let project = projector(app, rect, aspect);
    let painter = ui.painter().with_clip_rect(rect);
    let max_jump = rect.width() * 0.25;
    let ok = |p: egui::Pos2, q: egui::Pos2| (p.x - q.x).abs() < max_jump && (p.y - q.y).abs() < max_jump;
    let patch = matches!(m.place, Placement::Patch { .. });
    let (w, h) = (m.w as f64, m.h as f64);
    let (nx, ny) = if patch {
        let n = 72.0;
        if w >= h { (n as usize, ((n * h / w).round() as usize).max(2)) } else { (((n * w / h).round() as usize).max(2), n as usize) }
    } else {
        (120, 60)
    };
    let corner = |i: usize, j: usize| m.dir_of_pixel(-0.5 + i as f64 * w / nx as f64, -0.5 + j as f64 * h / ny as f64, r_m);
    let pts: Vec<Option<egui::Pos2>> = (0..=ny).flat_map(|j| (0..=nx).map(move |i| (i, j))).map(|(i, j)| project(corner(i, j))).collect();
    let at = |i: usize, j: usize| pts[j * (nx + 1) + i];
    let land = egui::Color32::from_rgba_unmultiplied(110, 210, 90, 110);
    let sea = egui::Color32::from_rgba_unmultiplied(40, 70, 140, 50);
    let mut mesh = egui::Mesh::default();
    for j in 0..ny {
        for i in 0..nx {
            let (Some(a), Some(b), Some(cc), Some(d)) = (at(i, j), at(i + 1, j), at(i + 1, j + 1), at(i, j + 1)) else { continue };
            if !(ok(a, b) && ok(b, cc) && ok(cc, d) && ok(d, a)) {
                continue;
            }
            let centre = m.dir_of_pixel(-0.5 + (i as f64 + 0.5) * w / nx as f64, -0.5 + (j as f64 + 0.5) * h / ny as f64, r_m);
            let is_land = m.coverage(centre, r_m).is_some_and(|v| v > 0.5);
            let col = if is_land {
                land
            } else if patch {
                sea
            } else {
                continue;
            };
            let k = mesh.vertices.len() as u32;
            for p in [a, b, cc, d] {
                mesh.colored_vertex(p, col);
            }
            mesh.add_triangle(k, k + 1, k + 2);
            mesh.add_triangle(k, k + 2, k + 3);
        }
    }
    painter.add(egui::Shape::mesh(mesh));
    if patch {
        // Outline, densely sampled (edges are curves on the sphere).
        let per_side = 48;
        let mut ring = Vec::new();
        for k in 0..per_side {
            let t = k as f64 / per_side as f64;
            ring.push((t * w - 0.5, -0.5));
        }
        for k in 0..per_side {
            let t = k as f64 / per_side as f64;
            ring.push((w - 0.5, t * h - 0.5));
        }
        for k in 0..per_side {
            let t = k as f64 / per_side as f64;
            ring.push(((1.0 - t) * w - 0.5, h - 0.5));
        }
        for k in 0..per_side {
            let t = k as f64 / per_side as f64;
            ring.push((-0.5, (1.0 - t) * h - 0.5));
        }
        let ring: Vec<Option<egui::Pos2>> = ring.iter().map(|&(x, y)| project(m.dir_of_pixel(x, y, r_m))).collect();
        let col = if app.coasts.placing { egui::Color32::from_rgb(255, 220, 90) } else { egui::Color32::from_white_alpha(200) };
        for k in 0..ring.len() {
            if let (Some(p), Some(q)) = (ring[k], ring[(k + 1) % ring.len()])
                && ok(p, q)
            {
                painter.line_segment([p, q], (1.5, col));
            }
        }
        // The top edge's midpoint gets a tick, so the rotation reads.
        if let (Some(p), Some(q)) = (project(m.dir_of_pixel(w * 0.5 - 0.5, -0.5, r_m)), project(m.dir_of_pixel(w * 0.5 - 0.5, h * 0.06, r_m))) {
            painter.line_segment([p, q], (2.0, col));
        }
    }
}
