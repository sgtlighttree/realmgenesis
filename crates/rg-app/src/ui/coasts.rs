//! R3 coast controls: where the coasts come from (procedural, procedural strict, or an
//! imported PNG/SVG mask placed on the sphere), and the brush (R3.3): paint land, sea, or
//! force a coast active or passive, then regenerate with the same seed.

use crate::app::RealmApp;
use glam::DVec3;
use rg_pipeline::mask::{Placement, Stroke, StrokeOp};
use std::sync::Arc;

/// UI state for the coast section (not part of the world params).
pub struct CoastUi {
    pub path: String,
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
}

impl Default for CoastUi {
    fn default() -> Self {
        Self {
            path: String::new(),
            patch: false,
            lat: 0.0,
            lon: 0.0,
            km_per_px: 10.0,
            rotation: 0.0,
            msg: String::new(),
            brush: None,
            brush_km: 150.0,
            drawing: Vec::new(),
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

pub fn coasts_section(app: &mut RealmApp, ui: &mut egui::Ui) {
    let c = &mut app.params.tectonics.conditioning;
    let ui_state = &mut app.coasts;
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
            .button("Import coastline…")
            .on_hover_text("Open a .png (white land on black sea, or opaque land on transparent sea) or an .svg (every filled shape is land). You can also drop the file on the window.")
            .clicked()
            && let Some(path) = rfd::FileDialog::new().add_filter("Coastline map", &["png", "svg"]).pick_file()
        {
            import_file(c, ui_state, &path);
            source = 2;
        }
        if source != 2 {
            ui.label(egui::RichText::new("or drop a .png / .svg here").small().color(super::theme::MUTED));
        }
    });
    if source == 2 {
        ui.horizontal(|ui| {
            ui.label("file");
            ui.text_edit_singleline(&mut ui_state.path).on_hover_text(
                "Path to a .png (white land on black sea, or opaque land on transparent sea) or an .svg (every filled shape is land).",
            );
        });
        ui.horizontal(|ui| {
            ui.radio_value(&mut ui_state.patch, false, "whole world").on_hover_text("A 2:1 equirectangular map of the whole planet.");
            ui.radio_value(&mut ui_state.patch, true, "regional patch").on_hover_text("A map of one region, placed by its centre, scale and rotation.");
        });
        if ui_state.patch {
            ui.horizontal(|ui| {
                ui.add(egui::DragValue::new(&mut ui_state.lat).range(-90.0..=90.0).speed(0.2).prefix("lat "));
                ui.add(egui::DragValue::new(&mut ui_state.lon).range(-180.0..=180.0).speed(0.2).prefix("lon "));
            });
            ui.horizontal(|ui| {
                ui.add(egui::DragValue::new(&mut ui_state.km_per_px).range(0.01..=500.0).speed(0.1).prefix("km/px "))
                    .on_hover_text("PNG: km per pixel. SVG: km per drawing unit.");
                ui.add(egui::DragValue::new(&mut ui_state.rotation).range(-180.0..=180.0).speed(0.5).prefix("rot° "));
            });
        }
        if ui.button(if c.mask.is_some() { "Reload with this placement" } else { "Load mask" }).clicked() {
            load(c, ui_state);
        }
        if c.mask.is_some() {
            ui.checkbox(&mut c.outside_procedural, "procedural land outside the patch")
                .on_hover_text("Other continents around the imported region, kept 800 km away. Off: open ocean.");
        } else {
            ui.label(egui::RichText::new("no mask loaded yet").small().color(super::theme::MUTED));
        }
        if !ui_state.msg.is_empty() {
            ui.label(egui::RichText::new(&ui_state.msg).small().color(super::theme::MUTED));
        }
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
            let on = ui_state.brush == Some(op);
            if ui.selectable_label(on, name).on_hover_text(tip).clicked() {
                ui_state.brush = if on { None } else { Some(op) };
            }
        }
    });
    if ui_state.brush.is_some() {
        ui.add(egui::Slider::new(&mut ui_state.brush_km, 20.0..=1500.0).logarithmic(true).text("brush km"));
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

/// Load `st.path` with the placement in the panel; on success the mask drives the coasts.
fn load(c: &mut rg_pipeline::mask::Conditioning, st: &mut CoastUi) {
    match rg_pipeline::import::load_mask(std::path::Path::new(st.path.trim()), st.placement()) {
        Ok(m) => {
            st.msg = format!("loaded {}×{} {}", m.w, m.h, m.source_kind.to_uppercase());
            if let Placement::Patch { .. } = m.place {
                let side_km = st.km_per_px * m.w.max(m.h) as f64;
                if side_km > 6000.0 && m.source_kind == "png" {
                    st.msg.push_str(" · wide patch: distances distort past ~3,000 km from its centre; a world map fits better");
                }
            }
            st.msg.push_str(" · Generate to apply");
            c.mask = Some(Arc::new(m));
            c.strict = true;
        }
        Err(e) => st.msg = e,
    }
}

/// A picked or dropped file: a 2:1 image is taken as a whole-world map, anything else as a
/// ~3,000 km regional patch at the current centre (adjust it, then reload).
fn import_file(c: &mut rg_pipeline::mask::Conditioning, st: &mut CoastUi, path: &std::path::Path) {
    st.path = path.display().to_string();
    st.patch = false;
    load(c, st);
    let Some(m) = c.mask.as_ref() else { return };
    let ratio = m.w as f64 / m.h.max(1) as f64;
    if (ratio - 2.0).abs() > 0.05 {
        st.patch = true;
        st.km_per_px = (3000.0 / m.w.max(m.h) as f64 * 100.0).round() / 100.0;
        load(c, st);
        st.msg.push_str(" · not 2:1, so placed as a regional patch: set its centre and scale");
    }
}

/// Files dropped on the window: the first one is imported as the coast mask (a wrong file
/// type shows the loader's error). Also draws a
/// hint while a file hovers over the window.
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
    import_file(&mut app.params.tectonics.conditioning, &mut app.coasts, &path);
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
    use rg_render::ViewMode;
    let Some(op) = app.coasts.brush else { return };
    let r_m = app.params.tectonics.radius_m;
    let painter = ui.painter().with_clip_rect(rect);
    let to_screen = |fx: f64, fy: f64| egui::pos2(rect.left() + fx as f32 * rect.width(), rect.top() + fy as f32 * rect.height());
    let project = |d: DVec3| -> Option<egui::Pos2> {
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
    };
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
