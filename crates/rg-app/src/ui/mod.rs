//! The egui front end. It reads `RealmApp` state and pushes `Command`s; it never runs
//! generation itself.

pub mod theme;

use crate::app::RealmApp;
use crate::commands::Command;
use theme::section;

pub fn draw(app: &mut RealmApp, ui: &mut egui::Ui, rs: &egui_wgpu::RenderState) {
    egui::Panel::left("controls")
        .resizable(true)
        .default_size(270.0)
        .show(ui, |ui| controls(app, ui));

    egui::CentralPanel::default()
        .frame(egui::Frame::NONE.fill(theme::BG))
        .show(ui, |ui| viewport(app, ui, rs));
}

fn controls(app: &mut RealmApp, ui: &mut egui::Ui) {
    ui.add_space(8.0);
    ui.heading("RealmGenesis");
    ui.label(egui::RichText::new("v2 · M0 preview").small().color(theme::MUTED));
    ui.separator();

    section(ui, "World");
    let mut changed = false;
    ui.horizontal(|ui| {
        ui.label("Seed");
        changed |= ui.add(egui::DragValue::new(&mut app.params.seed).range(0..=999_999)).changed();
        if ui.button("Random").on_hover_text("New random seed (R)").clicked() {
            app.push(Command::RandomizeSeed);
        }
    });
    changed |= ui
        .add(egui::Slider::new(&mut app.params.radius_km, 1000.0..=12000.0).text("radius km").integer())
        .changed();
    changed |= ui.add(egui::Slider::new(&mut app.params.sea_level, -0.4..=0.4).text("sea level")).changed();

    section(ui, "Display");
    changed |= ui.add(egui::Slider::new(&mut app.params.exaggeration, 1.0..=60.0).text("relief ×")).changed();
    changed |= ui.add(egui::Slider::new(&mut app.params.mesh_res, 32..=320).text("mesh res")).changed();
    if changed {
        app.push(Command::Regenerate);
    }
    if ui.button("Reset view").clicked() {
        app.push(Command::ResetCamera);
    }

    section(ui, "Performance");
    let s = &app.stats;
    ui.monospace(format!("frame  {:5.1} ms  (worst {:4.1})", s.avg_ms(), s.worst_ms()));
    ui.monospace(format!("gen    {:5.0} ms", s.last_gen_ms));
    ui.monospace(format!("tris   {:>8}", s.triangles));
    ui.label(egui::RichText::new(&app.adapter_info).small().color(theme::MUTED));

    if ui.input(|i| i.key_pressed(egui::Key::R)) && !ui.ctx().egui_wants_keyboard_input() {
        app.push(Command::RandomizeSeed);
    }
}

fn viewport(app: &mut RealmApp, ui: &mut egui::Ui, rs: &egui_wgpu::RenderState) {
    let rect = ui.available_rect_before_wrap();
    let response = ui.allocate_rect(rect, egui::Sense::click_and_drag());
    if response.dragged() {
        let d = response.drag_delta();
        app.camera.drag(d.x as f64, d.y as f64, rect.height() as f64);
    }
    if response.hovered() {
        let (scroll, zoom) = ui.input(|i| (i.smooth_scroll_delta.y, i.zoom_delta()));
        if scroll != 0.0 {
            app.camera.zoom((-scroll as f64 * 0.0025).exp());
        }
        if zoom != 1.0 {
            app.camera.zoom(1.0 / zoom as f64);
        }
    }
    let ppp = ui.ctx().pixels_per_point();
    let size_px = [(rect.width() * ppp).round() as u32, (rect.height() * ppp).round() as u32];
    if let Some(id) = app.globe_texture(rs, size_px) {
        ui.painter().image(
            id,
            rect,
            egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0)),
            egui::Color32::WHITE,
        );
    }
    // The globe is drawn every frame; keep the loop going while interacting.
    if response.dragged() || response.hovered() {
        ui.ctx().request_repaint();
    }
}
