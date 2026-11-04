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
        if ui.button("Load mask").clicked() {
            match rg_pipeline::import::load_mask(std::path::Path::new(ui_state.path.trim()), ui_state.placement()) {
                Ok(m) => {
                    ui_state.msg = format!("loaded {}×{} {}", m.w, m.h, m.source_kind.to_uppercase());
                    if let Placement::Patch { .. } = m.place {
                        let side_km = ui_state.km_per_px * m.w.max(m.h) as f64;
                        if side_km > 6000.0 && m.source_kind == "png" {
                            ui_state.msg.push_str(" · wide patch: distances distort past ~3,000 km from its centre; a world map fits better");
                        }
                    }
                    c.mask = Some(Arc::new(m));
                    c.strict = true;
                }
                Err(e) => ui_state.msg = e,
            }
        }
        if c.mask.is_some() {
            ui.checkbox(&mut c.outside_procedural, "procedural land outside the patch")
                .on_hover_text("Other continents around the imported region, kept 800 km away. Off: open ocean.");
        } else if !ui_state.msg.is_empty() || source == 2 {
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
