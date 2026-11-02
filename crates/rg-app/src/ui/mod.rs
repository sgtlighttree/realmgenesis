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
    let show = *app.show_panel.get_or_insert_with(|| ui.ctx().content_rect().width() >= 700.0);
    if show {
        egui::Panel::left("controls")
            .resizable(true)
            .default_size(290.0)
            .show(ui, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| controls(app, ui));
            });
    }
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
        if i.key_pressed(egui::Key::H) {
            v.push(Command::TogglePanel);
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
    detail_picker(ui, &mut t.res, t.radius_m);
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
    ui.add(egui::Slider::new(&mut app.display.exaggeration, 1.0..=60.0).text("relief × (orbit)"))
        .on_hover_text("Vertical exaggeration seen from orbit. It fades to 1.5× as you zoom in to the ground.");
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
    let fps = match s.fps(crate::app::now_s()) {
        Some(f) => format!("{f:6.0}"),
        None => "  idle".to_string(),
    };
    ui.monospace(format!("fps      {fps}"))
        .on_hover_text("Frames actually presented in the last second. 'idle' = not repainting (nothing is moving).");
    ui.monospace(format!("cpu      {:6.1} ms (worst {:.1})", s.cpu_avg_ms(), s.cpu_worst_ms()))
        .on_hover_text("Time our code spends per frame (UI + render submission), averaged over 2 s. GPU time isn't included.");
    if let Some(w) = app.world.as_ref() {
        ui.monospace(format!("tecton.  {:6.0} ms", w.timings.tectonics_ms));
        ui.monospace(format!("climate  {:6.0} ms", w.timings.climate_ms));
    }
    ui.monospace(format!("bake     {:6.0} ms", s.bake_ms));
    if app.view == ViewMode::Globe {
        ui.monospace(format!("tiles    {:>4} drawn  {:>4} cached", s.tiles_drawn, s.tiles_resident))
            .on_hover_text("Terrain tiles on screen / held in memory. Each is 64×64 cells.");
        ui.monospace(format!("loading  {:>4}", s.tiles_in_flight));
        let alt = app.camera.altitude_m();
        let alt = if alt >= 10_000.0 { format!("{:.0} km", alt / 1000.0) } else { format!("{alt:.0} m") };
        ui.monospace(format!("altitude {alt:>8}  relief ×{:.1}", s.exaggeration));
    }
    ui.label(egui::RichText::new(&app.adapter_info).small().color(theme::MUTED));
    ui.add_space(4.0);
    ui.label(
        egui::RichText::new(
            "G generate · R random · H panel · M map/globe · 1–8 layers · drag to move · scroll to zoom · \
             right-drag (or Shift+drag, or two fingers) to tilt and turn",
        )
            .small()
            .color(theme::MUTED),
    );
}

/// Base-grid resolution presets (cube-face cells per edge). Cell size is about
/// `pi/2 * radius / res`, so on an Earth-sized planet Draft is ~52 km and Ultra ~10 km.
const DETAIL_PRESETS: [(&str, u32); 5] =
    [("Draft", 192), ("Normal", 320), ("High", 448), ("XHigh", 640), ("Ultra", 1024)];
const DETAIL_MIN_RES: u32 = 32;
const DETAIL_MAX_RES: u32 = 1536;

/// Size of one grid cell in km for a given resolution and planet radius.
fn cell_km(res: u32, radius_m: f64) -> f64 {
    std::f64::consts::FRAC_PI_2 * radius_m / 1000.0 / res as f64
}

/// The "Detail" row: preset buttons, a custom slider in km per cell, and a readout.
fn detail_picker(ui: &mut egui::Ui, res: &mut u32, radius_m: f64) {
    // Whether "Custom" is the chosen mode. It's kept in egui's temporary memory, because
    // a custom value can equal a preset's value and the buttons then look the same.
    let mode_id = egui::Id::new("detail_custom_mode");
    let mut custom = ui.data(|d| d.get_temp::<bool>(mode_id)).unwrap_or(false);
    ui.horizontal_wrapped(|ui| {
        ui.label("Detail");
        for (label, res) in DETAIL_PRESETS {
            if ui.selectable_label(!custom && *res == res, label).clicked() {
                *res = res;
                custom = false;
            }
        }
        let is_custom_value = !DETAIL_PRESETS.iter().any(|&(_, r)| r == *res);
        if ui.selectable_label(custom || is_custom_value, "Custom").clicked() {
            custom = true;
        }
    });
    if custom || !DETAIL_PRESETS.iter().any(|&(_, r)| r == *res) {
        let mut km = cell_km(*res, radius_m);
        let slider = egui::Slider::new(&mut km, 5.0..=100.0).logarithmic(true).text("km / cell").max_decimals(0);
        if ui.add(slider).changed() {
            let res = (std::f64::consts::FRAC_PI_2 * radius_m / 1000.0 / km).round();
            *res = (res as u32).clamp(DETAIL_MIN_RES, DETAIL_MAX_RES);
        }
    }
    ui.data_mut(|d| d.insert_temp(mode_id, custom));

    let points = 6.0 * (*res as f64 + 1.0).powi(2);
    ui.label(
        egui::RichText::new(format!(
            "≈ {:.0} km cells · {:.1} M points · est. {}",
            cell_km(*res, radius_m),
            points / 1.0e6,
            fmt_seconds(estimate_generate_s(*res)),
        ))
        .small()
        .color(theme::MUTED),
    );
}

/// Rough generation time (tectonics + climate), fitted to release-build timings on a
/// 4-core x86 machine: about 16 µs per cell of `res²`, near-linear in the point count.
/// The web build is single-threaded, so it's several times slower.
fn estimate_generate_s(res: u32) -> f64 {
    let native = 1.6e-5 * (res as f64).powi(2) + 0.2;
    if cfg!(target_arch = "wasm32") { native * 3.0 } else { native }
}

fn fmt_seconds(s: f64) -> String {
    if s < 10.0 { format!("{s:.1} s") } else { format!("{s:.0} s") }
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
    let multi = ui.input(|i| i.multi_touch());
    if let (ViewMode::Globe, Some(mt)) = (app.view, multi) {
        // Two fingers: twist to turn, slide up/down to tilt (pinch zoom is handled below).
        app.camera.orbit(-mt.rotation_delta as f64, -mt.translation_delta.y as f64 * 0.005);
    } else if response.dragged() {
        let d = response.drag_delta();
        let rotate = response.dragged_by(egui::PointerButton::Secondary)
            || response.dragged_by(egui::PointerButton::Middle)
            || ui.input(|i| i.modifiers.shift);
        match app.view {
            ViewMode::Globe if rotate => app.camera.orbit(d.x as f64 * 0.006, -d.y as f64 * 0.006),
            ViewMode::Globe => app.camera.pan(d.x as f64, d.y as f64, rect.height() as f64),
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
    // Sidebar toggle, top-left of the viewport. Added after the viewport's drag area so
    // it wins the click.
    let open = app.show_panel.unwrap_or(true);
    let toggle = egui::Button::new(egui::RichText::new(if open { "«" } else { "»" }).size(18.0))
        .fill(theme::PANEL.gamma_multiply(0.85))
        .min_size(egui::vec2(36.0, 36.0));
    let toggle_rect = egui::Rect::from_min_size(rect.min + egui::vec2(8.0, 8.0), egui::vec2(36.0, 36.0));
    if ui
        .put(toggle_rect, toggle)
        .on_hover_text(if open { "Hide controls (H)" } else { "Show controls (H)" })
        .clicked()
    {
        app.push(Command::TogglePanel);
    }
    if let Some(msg) = app.busy {
        let r = egui::Rect::from_center_size(rect.center(), egui::vec2(260.0, 48.0));
        ui.painter().rect_filled(r, 6.0, theme::PANEL.gamma_multiply(0.92));
        ui.painter().text(r.center(), egui::Align2::CENTER_CENTER, msg, egui::FontId::proportional(15.0), theme::TEXT);
    }
    if response.dragged() || response.hovered() || app.busy.is_some() || app.streaming {
        ui.ctx().request_repaint();
    }
}
