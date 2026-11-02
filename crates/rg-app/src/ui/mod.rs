//! The egui front end. It reads `RealmApp` state and pushes `Command`s; it never runs
//! generation itself.

pub mod theme;

use crate::app::RealmApp;
use crate::commands::Command;
use rg_pipeline::Layer;
use rg_render::ViewMode;
use theme::section;

pub fn draw(app: &mut RealmApp, ui: &mut egui::Ui, rs: &egui_wgpu::RenderState) {
    shortcuts(app, ui);
    egui::Panel::left("controls")
        .resizable(true)
        .default_size(290.0)
        .show(ui, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| controls(app, ui));
        });
    egui::Panel::bottom("status").show(ui, |ui| status_bar(app, ui));
    egui::CentralPanel::default()
        .frame(egui::Frame::NONE.fill(theme::BG))
        .show(ui, |ui| viewport(app, ui, rs));
}

fn shortcuts(app: &mut RealmApp, ui: &mut egui::Ui) {
    if ui.ctx().egui_wants_keyboard_input() {
        return;
    }
    let view = app.view;
    let keys = ui.input(|i| {
        let mut v = Vec::new();
        if i.key_pressed(egui::Key::G) {
            v.push(Command::Regenerate);
        }
        if i.key_pressed(egui::Key::R) {
            v.push(Command::RandomizeSeed);
        }
        if i.key_pressed(egui::Key::M) {
            v.push(Command::SetView(if view == ViewMode::Globe { ViewMode::Map } else { ViewMode::Globe }));
        }
        let nums = [
            egui::Key::Num1,
            egui::Key::Num2,
            egui::Key::Num3,
            egui::Key::Num4,
            egui::Key::Num5,
            egui::Key::Num6,
            egui::Key::Num7,
            egui::Key::Num8,
        ];
        for (k, l) in nums.iter().zip(Layer::ALL) {
            if i.key_pressed(*k) {
                v.push(Command::SetLayer(l));
            }
        }
        v
    });
    for c in keys {
        app.push(c);
    }
}

fn controls(app: &mut RealmApp, ui: &mut egui::Ui) {
    ui.add_space(8.0);
    ui.heading("RealmGenesis");
    ui.label(egui::RichText::new("v2 · MVP preview").small().color(theme::MUTED));
    ui.separator();

    // ---------------------------------------------------------------- world
    section(ui, "World");
    let mut random = false;
    ui.horizontal(|ui| {
        ui.label("Seed");
        ui.add(egui::DragValue::new(&mut app.params.tectonics.seed).range(0..=999_999));
        random = ui.button("Random").on_hover_text("New random seed (R)").clicked();
    });
    if random {
        app.push(Command::RandomizeSeed);
    }
    let t = &mut app.params.tectonics;
    let mut radius_km = t.radius_m / 1000.0;
    if ui.add(egui::Slider::new(&mut radius_km, 1500.0..=12000.0).text("radius km").integer()).changed() {
        t.radius_m = radius_km * 1000.0;
    }
    let mut land = t.land_fraction * 100.0;
    if ui.add(egui::Slider::new(&mut land, 5.0..=70.0).text("land %").integer()).changed() {
        t.land_fraction = land / 100.0;
    }
    ui.add(egui::Slider::new(&mut t.plate_count, 4..=30).text("plates"));
    ui.add(egui::Slider::new(&mut t.continent_count_hint, 1..=12).text("continents"));
    ui.add(egui::Slider::new(&mut t.mountain_scale, 0.3..=2.0).text("mountains"));
    ui.add(egui::Slider::new(&mut t.roughness, 0.3..=2.0).text("roughness"));
    ui.horizontal(|ui| {
        ui.label("Detail");
        for (label, res) in [("Draft", 96), ("Normal", 192), ("High", 288)] {
            ui.selectable_value(&mut t.res, res, label);
        }
    });
    let generate = egui::Button::new(egui::RichText::new("Generate world").strong())
        .fill(theme::ACCENT.linear_multiply(0.45));
    if ui.add_sized([ui.available_width(), 30.0], generate).on_hover_text("G").clicked() {
        app.push(Command::Regenerate);
    }

    // ---------------------------------------------------------------- climate
    section(ui, "Planet & climate");
    let p = &mut app.params.planet;
    let mut climate_changed = false;
    climate_changed |= ui.add(egui::Slider::new(&mut p.axial_tilt_deg, 0.0..=60.0).text("axial tilt °")).drag_stopped();
    climate_changed |= ui.add(egui::Slider::new(&mut p.rotation_hours, 6.0..=96.0).text("day length h")).drag_stopped();
    climate_changed |= ui.add(egui::Slider::new(&mut p.solar_constant_rel, 0.7..=1.3).text("sunlight ×")).drag_stopped();
    climate_changed |=
        ui.add(egui::Slider::new(&mut p.global_temp_offset_c, -15.0..=15.0).text("temp offset °C")).drag_stopped();
    climate_changed |= ui.add(egui::Slider::new(&mut p.global_precip_mult, 0.3..=2.5).text("rainfall ×")).drag_stopped();
    let rerun = ui.add_enabled(app.world.is_some(), egui::Button::new("Re-run climate"));
    if rerun.clicked() || (climate_changed && app.world.is_some()) {
        app.push(Command::RerunClimate);
    }
    ui.label(egui::RichText::new("Tier A (rule-based). The physics solver comes later.").small().color(theme::MUTED));

    // ---------------------------------------------------------------- view
    section(ui, "View");
    ui.horizontal(|ui| {
        for (label, v) in [("Globe", ViewMode::Globe), ("Map", ViewMode::Map)] {
            if ui.selectable_label(app.view == v, label).on_hover_text("M toggles").clicked() {
                app.push(Command::SetView(v));
            }
        }
    });
    let mut layer = app.layer;
    egui::ComboBox::from_id_salt("layer")
        .width(ui.available_width() - 8.0)
        .selected_text(layer.label())
        .show_ui(ui, |ui| {
            for (n, l) in Layer::ALL.into_iter().enumerate() {
                ui.selectable_value(&mut layer, l, format!("{}  {}", n + 1, l.label()));
            }
        });
    if layer != app.layer {
        app.push(Command::SetLayer(layer));
    }
    let mut mesh_changed = false;
    mesh_changed |= ui.add(egui::Slider::new(&mut app.display.exaggeration, 1.0..=80.0).text("relief ×")).drag_stopped();
    mesh_changed |= ui.add(egui::Slider::new(&mut app.display.mesh_res, 64..=320).text("globe mesh")).drag_stopped();
    if mesh_changed {
        app.push(Command::RebuildMesh);
    }
    ui.add(egui::Slider::new(&mut app.style.map_exaggeration, 1.0..=40.0).text("map hillshade"));
    if ui.button("Reset view").clicked() {
        app.push(Command::ResetCamera);
    }

    // ---------------------------------------------------------------- files
    section(ui, "File");
    ui.horizontal_wrapped(|ui| {
        if ui.button("Save").on_hover_text("~/RealmGenesis/world-<seed>.rgworld").clicked() {
            app.push(Command::SaveWorld);
        }
        if ui.button("Load latest").clicked() {
            app.push(Command::LoadWorld);
        }
        if ui.button("Export heightmap").on_hover_text("16-bit PNG + JSON sidecar, 4096×2048").clicked() {
            app.push(Command::ExportHeightmap);
        }
    });

    // ---------------------------------------------------------------- perf
    section(ui, "Performance");
    let s = &app.stats;
    ui.monospace(format!("frame    {:6.1} ms (worst {:.0})", s.avg_ms(), s.worst_ms()));
    if let Some(w) = app.world.as_ref() {
        ui.monospace(format!("tecton.  {:6.0} ms", w.timings.tectonics_ms));
        ui.monospace(format!("climate  {:6.0} ms", w.timings.climate_ms));
    }
    ui.monospace(format!("bake     {:6.0} ms", s.bake_ms));
    ui.monospace(format!("mesh     {:6.0} ms  {} tris", s.mesh_ms, s.triangles));
    ui.label(egui::RichText::new(&app.adapter_info).small().color(theme::MUTED));
    ui.add_space(4.0);
    ui.label(
        egui::RichText::new("G generate · R random · M map/globe · 1–8 layers · drag to move · scroll to zoom")
            .small()
            .color(theme::MUTED),
    );
}

fn status_bar(app: &RealmApp, ui: &mut egui::Ui) {
    ui.horizontal(|ui| {
        ui.label(egui::RichText::new(&app.status).small());
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            if let Some(h) = &app.hover {
                ui.monospace(h);
            }
        });
    });
}

fn viewport(app: &mut RealmApp, ui: &mut egui::Ui, rs: &egui_wgpu::RenderState) {
    let rect = ui.available_rect_before_wrap();
    let response = ui.allocate_rect(rect, egui::Sense::click_and_drag());
    let aspect = (rect.width() / rect.height().max(1.0)) as f64;
    if response.dragged() {
        let d = response.drag_delta();
        match app.view {
            ViewMode::Globe => app.camera.drag(d.x as f64, d.y as f64, rect.height() as f64),
            ViewMode::Map => app.map.drag(d.x as f64, d.y as f64, rect.width() as f64),
        }
    }
    if response.hovered() {
        let (scroll, zoom) = ui.input(|i| (i.smooth_scroll_delta.y, i.zoom_delta()));
        let mut factor = 1.0;
        if scroll != 0.0 {
            factor *= (-scroll as f64 * 0.0025).exp();
        }
        if zoom != 1.0 {
            factor /= zoom as f64;
        }
        if factor != 1.0 {
            match app.view {
                ViewMode::Globe => app.camera.zoom(factor),
                ViewMode::Map => app.map.zoom(factor),
            }
        }
    }
    app.hover = response.hover_pos().and_then(|p| {
        let fx = ((p.x - rect.left()) / rect.width()) as f64;
        let fy = ((p.y - rect.top()) / rect.height()) as f64;
        let ll = match app.view {
            ViewMode::Globe => app.camera.pick(fx, fy, aspect),
            ViewMode::Map => {
                use std::f64::consts::{FRAC_PI_2, PI, TAU};
                let (lon, lat) = app.map.lonlat_at(fx, fy, aspect);
                (lat.abs() <= FRAC_PI_2).then(|| (lat, (lon + PI).rem_euclid(TAU) - PI))
            }
        };
        ll.and_then(|(lat, lon)| app.describe(lat, lon))
    });

    let ppp = ui.ctx().pixels_per_point();
    let size_px = [(rect.width() * ppp).round() as u32, (rect.height() * ppp).round() as u32];
    if let Some(id) = app.planet_texture(rs, size_px) {
        ui.painter().image(
            id,
            rect,
            egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0)),
            egui::Color32::WHITE,
        );
    }
    if let Some(msg) = app.busy {
        let r = egui::Rect::from_center_size(rect.center(), egui::vec2(260.0, 48.0));
        ui.painter().rect_filled(r, 6.0, theme::PANEL.gamma_multiply(0.92));
        ui.painter().text(r.center(), egui::Align2::CENTER_CENTER, msg, egui::FontId::proportional(15.0), theme::TEXT);
    }
    if response.dragged() || response.hovered() || app.busy.is_some() {
        ui.ctx().request_repaint();
    }
}
