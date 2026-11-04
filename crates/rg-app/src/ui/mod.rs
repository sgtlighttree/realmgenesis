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
    if app.console_open {
        egui::Panel::bottom("console").resizable(true).default_size(150.0).show(ui, |ui| {
            egui::ScrollArea::vertical().stick_to_bottom(true).auto_shrink([false, false]).show(ui, |ui| {
                for line in app.console.lines() {
                    ui.label(egui::RichText::new(line).monospace().small());
                }
            });
        });
    }
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
        if i.key_pressed(egui::Key::Escape) {
            v.push(Command::CancelGeneration);
        }
        if i.key_pressed(egui::Key::H) {
            v.push(Command::TogglePanel);
        }
        if i.key_pressed(egui::Key::F12) {
            v.push(Command::SaveScreenshot);
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
    let busy = app.generating();
    ui.add_enabled_ui(!busy, |ui| {
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
        // ---------------------------------------------------------------- erosion
        let e = &mut app.params.erosion;
        ui.checkbox(&mut e.enabled, "River erosion")
            .on_hover_text("Stream-power erosion (Braun & Willett 2013): rivers carve valleys and drainage networks into the tectonic relief.");
        ui.add_enabled_ui(e.enabled, |ui| {
            ui.add(egui::Slider::new(&mut e.duration_myr, 1.0..=60.0).logarithmic(true).text("maturity Myr"))
                .on_hover_text("Simulated time. Short: young, sharp ranges with shallow valleys. Long: mature, deeply dissected relief.");
            // (s03: this used to be "uplift focus" = uplift_exponent, which only applies to the
            // legacy elevation-based uplift, so it did nothing with tectonic uplift.)
            ui.add(egui::Slider::new(&mut e.belt_width_km, 150.0..=800.0).text("belt width km"))
                .on_hover_text("How far mountain-building uplift reaches from a converging plate boundary. Narrow: tight ranges along the boundaries. Wide: broad uplands.");
            ui.add(egui::Slider::new(&mut e.convergence_smoothing_km, 0.0..=600.0).text("uplift blend km"))
                .on_hover_text("How smoothly uplift blends between boundary stretches that close at different speeds. 0: terraced \"floors\" with straight escarpments between them. Higher: natural ramps.");
            ui.add(egui::Slider::new(&mut e.node_jitter, 0.0..=0.4).text("drainage jitter"))
                .on_hover_text("R1 experiment: nudges the erosion grid's points off the lattice so rivers and ridges aren't limited to 8 directions (removes straight, parallel gullies). 0 = the regular grid. About 20% slower erosion.");
            ui.add(egui::Slider::new(&mut e.river_min_area_km2, 2_000.0..=200_000.0).logarithmic(true).text("river km²"))
                .on_hover_text("Smallest drainage area drawn as a river.");
        });
        let generate = egui::Button::new(egui::RichText::new("Generate world").strong())
            .fill(theme::ACCENT.linear_multiply(0.45));
        if ui.add_sized([ui.available_width(), 30.0], generate).on_hover_text("G").clicked() {
            app.push(Command::Regenerate);
        }
    });

    // ---------------------------------------------------------------- climate
    section(ui, "Planet & climate");
    let mut climate_changed = false;
    ui.add_enabled_ui(!busy, |ui| {
        let p = &mut app.params.planet;
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
    });
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
        .on_hover_text("Vertical exaggeration seen from orbit. It fades to 1.5× as you zoom in to the ground, unless kept below.");
    ui.checkbox(&mut app.display.keep_relief, "Keep relief when zooming in")
        .on_hover_text("Use the same exaggeration at every distance instead of fading to 1.5× near the ground.");
    ui.add(egui::Slider::new(&mut app.style.map_exaggeration, 1.0..=40.0).text("map hillshade"));
    ui.checkbox(&mut app.style.graticule, "Lat/lon grid")
        .on_hover_text("Latitude and longitude lines with degree labels. The spacing adapts to the zoom.");
    ui.add(
        egui::Slider::new(&mut app.display.max_px_per_texel, 4.0..=512.0)
            .logarithmic(true)
            .text(format!("zoom floor ({:.0} km)", app.zoom_floor_m / 1000.0)),
    )
    .on_hover_text("How far you can zoom in: the screen pixels one texel of the stored terrain may cover. Higher goes closer (and blurrier).");
    if ui.button("Reset view").clicked() {
        app.push(Command::ResetCamera);
    }

    // ---------------------------------------------------------------- files
    section(ui, "File");
    ui.horizontal_wrapped(|ui| {
        if ui.button("Save").on_hover_text("~/RealmGenesis/world-<seed>.rgworld").clicked() {
            app.push(Command::SaveWorld);
        }
        ui.add_enabled_ui(!busy, |ui| {
            if ui.button("Load latest").clicked() {
                app.push(Command::LoadWorld);
            }
        });
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
        ui.monospace(format!("erosion  {:6.0} ms", w.timings.erosion_ms));
        if let Some(hy) = w.hydrology.as_ref().filter(|h| h.stats.grid_res > 0) {
            ui.monospace(format!("Hack h   {:6.2}   θ {:.2}", hy.stats.hack_h, hy.stats.concavity))
                .on_hover_text("Geomorphic checks. Hack's law exponent (real rivers 0.45–0.7) and channel concavity θ (real 0.3–0.7).");
        }
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
            "G generate · Esc cancel · R random · H panel · M map/globe · 1–8 layers · drag to move · scroll to zoom · \
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
// Native: up to res 2048 (~5 km cells on an Earth-size planet, ~4 GB peak estimated).
// Matt confirmed Ultra (1024) runs well on the M1 (2026-11-02). Web: capped at 1024,
// since wasm32 can't address more than 4 GB and 1024 measured 1.04 GB peak.
#[cfg(not(target_arch = "wasm32"))]
const DETAIL_MAX_RES: u32 = 2048;
#[cfg(target_arch = "wasm32")]
const DETAIL_MAX_RES: u32 = 1024;

/// Estimated peak memory for generation, GB, scaling with the point count (res²). Fitted
/// to the app's peak RSS on Matt's M1 at res 2002 (4.6 GB, s02, n=1; swap also grew by
/// ~4.4 GB there, so ~2048 is this machine's practical ceiling). s01's container: 1.04 GB at
/// 1024.
fn estimate_peak_gb(res: u32) -> f64 {
    1.2 * (res as f64 / 1024.0).powi(2)
}

/// Size of one grid cell in km for a given resolution and planet radius.
fn cell_km(res: u32, radius_m: f64) -> f64 {
    std::f64::consts::FRAC_PI_2 * radius_m / 1000.0 / res as f64
}

/// The "Detail" row: preset buttons, a custom slider in km per cell, and a readout.
fn detail_picker(ui: &mut egui::Ui, res: &mut u32, radius_m: f64) {
    // Whether "Custom" is the chosen mode. It's kept in egui's temporary memory, because
    // a custom value can equal a preset's value and the buttons then look the same.
    let mode_id = ui.make_persistent_id("detail_custom_mode");
    let mut custom = ui.data(|d| d.get_temp::<bool>(mode_id)).unwrap_or(false);
    ui.horizontal_wrapped(|ui| {
        ui.label("Detail");
        for (label, preset) in DETAIL_PRESETS {
            if ui.selectable_label(!custom && *res == preset, label).clicked() {
                *res = preset;
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
        let min_km = cell_km(DETAIL_MAX_RES, radius_m).ceil();
        let slider = egui::Slider::new(&mut km, min_km..=100.0).logarithmic(true).text("km / cell").max_decimals(1);
        if ui.add(slider).changed() {
            let cells = (std::f64::consts::FRAC_PI_2 * radius_m / 1000.0 / km).round();
            *res = (cells as u32).clamp(DETAIL_MIN_RES, DETAIL_MAX_RES);
        }
    }
    ui.data_mut(|d| d.insert_temp(mode_id, custom));

    let points = 6.0 * (*res as f64 + 1.0).powi(2);
    let gb = estimate_peak_gb(*res);
    let heavy = gb > 2.0;
    ui.label(
        egui::RichText::new(format!(
            "≈ {:.1} km cells · {:.1} M points · est. {} · ~{:.1} GB",
            cell_km(*res, radius_m),
            points / 1.0e6,
            fmt_seconds(estimate_generate_s(*res)),
            gb,
        ))
        .small()
        .color(if heavy { egui::Color32::from_rgb(222, 170, 80) } else { theme::MUTED }),
    )
    .on_hover_text("Estimated peak memory while generating. Above ~2 GB, close other heavy apps first.");
}

/// Rough generation time, fitted to the app on Matt's M1 Air (background thread at utility
/// QoS, s02, n=1 each): 3.0 s at 320, 40.5 s at 1024, 276 s at 2002. Grows faster than the
/// point count (memory pressure at high res), so the exponent is ~2.3. The web build is
/// single-threaded, so it's several times slower.
fn estimate_generate_s(res: u32) -> f64 {
    let native = 3.0 * (res as f64 / 320.0).powf(2.3);
    if cfg!(target_arch = "wasm32") { native * 3.0 } else { native }
}

fn fmt_seconds(s: f64) -> String {
    if s < 10.0 { format!("{s:.1} s") } else { format!("{s:.0} s") }
}

fn status_bar(app: &mut RealmApp, ui: &mut egui::Ui) {
    let activity = app.activity(crate::app::now_s());
    let mut toggle_console = false;
    ui.horizontal(|ui| {
        match &activity {
            Some(a) => ui.label(egui::RichText::new(a).monospace().small().color(theme::ACCENT)),
            None => ui.label(egui::RichText::new(&app.status).small()),
        };
        ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
            if ui.selectable_label(app.console_open, "Console").on_hover_text("Show what the app is doing (log)").clicked() {
                toggle_console = true;
            }
            if let Some(h) = &app.hover {
                ui.monospace(h);
            }
        });
    });
    if toggle_console {
        app.push(Command::ToggleConsole);
    }
}

fn viewport(app: &mut RealmApp, ui: &mut egui::Ui, rs: &egui_wgpu::RenderState) {
    let rect = ui.available_rect_before_wrap();
    let response = ui.allocate_rect(rect, egui::Sense::click_and_drag());
    let aspect = (rect.width() / rect.height().max(1.0)) as f64;
    let multi = ui.input(|i| i.multi_touch());
    let touching = multi.is_some();
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
    if app.style.graticule {
        graticule_labels(app, ui, rect, size_px);
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
    if app.view == ViewMode::Globe {
        stop_tiles_button(app, ui, rect);
    }
    if let Some(msg) = app.busy {
        let r = egui::Rect::from_center_size(rect.center(), egui::vec2(260.0, 48.0));
        ui.painter().rect_filled(r, 6.0, theme::PANEL.gamma_multiply(0.92));
        ui.painter().text(r.center(), egui::Align2::CENTER_CENTER, msg, egui::FontId::proportional(15.0), theme::TEXT);
    }
    if app.world.is_none()
        && let Some(a) = app.activity(crate::app::now_s())
    {
        let r = egui::Rect::from_center_size(rect.center(), egui::vec2(420.0, 48.0));
        ui.painter().rect_filled(r, 6.0, theme::PANEL.gamma_multiply(0.92));
        ui.painter().text(r.center(), egui::Align2::CENTER_CENTER, a, egui::FontId::monospace(14.0), theme::ACCENT);
    }
    // Repaint for real input only: hovering alone used to keep a 60 fps loop running.
    let input_moving = response.hovered()
        && ui.input(|i| i.pointer.is_moving() || i.smooth_scroll_delta != egui::Vec2::ZERO || i.zoom_delta() != 1.0);
    if response.dragged() || input_moving || touching || app.settling || app.busy.is_some() || app.streaming {
        ui.ctx().request_repaint();
    }
}

/// Degree labels for the lat/lon grid, drawn over the planet image. `size_px` is the planet
/// texture's size, which the lines' spacing was chosen for.
fn graticule_labels(app: &RealmApp, ui: &egui::Ui, rect: egui::Rect, size_px: [u32; 2]) {
    use rg_render::graticule::{LabelKind, labels};
    let painter = ui.painter().with_clip_rect(rect);
    let font = egui::FontId::proportional(11.0);
    let text_color = egui::Color32::from_white_alpha(225);
    let shadow = egui::Color32::from_black_alpha(190);
    // Labels whose anchor is within this many points of an earlier one are skipped.
    let min_gap = 34.0;
    let mut placed: Vec<egui::Pos2> = Vec::new();
    for l in labels(app.view, &app.camera, &app.map, size_px) {
        let at = egui::pos2(rect.left() + l.fx as f32 * rect.width(), rect.top() + l.fy as f32 * rect.height());
        if placed.iter().any(|p| p.distance(at) < min_gap) {
            continue;
        }
        placed.push(at);
        let (pos, align) = match (app.view, l.kind) {
            (ViewMode::Map, LabelKind::Parallel) => (egui::pos2(rect.left() + 5.0, at.y - 1.0), egui::Align2::LEFT_BOTTOM),
            (ViewMode::Map, LabelKind::Meridian) => (egui::pos2(at.x + 3.0, rect.bottom() - 4.0), egui::Align2::LEFT_BOTTOM),
            (_, LabelKind::Parallel) => (at + egui::vec2(4.0, -2.0), egui::Align2::LEFT_BOTTOM),
            (_, LabelKind::Meridian) => (at + egui::vec2(3.0, 3.0), egui::Align2::LEFT_TOP),
        };
        painter.text(pos + egui::vec2(1.0, 1.0), align, &l.text, font.clone(), shadow);
        painter.text(pos, align, &l.text, font.clone(), text_color);
    }
}

/// The red "stop background rendering" button, top-right of the viewport (the sidebar can
/// be hidden, so it isn't there). Stopping cancels queued terrain tiles; resuming picks up
/// whatever the view needs.
fn stop_tiles_button(app: &mut RealmApp, ui: &mut egui::Ui, rect: egui::Rect) {
    let paused = app.tiles_paused();
    let (label, fill, tip) = if paused {
        (
            "RESUME BACKGROUND RENDERING",
            theme::SURFACE.gamma_multiply(0.92),
            "Terrain tiles are stopped, so zooming in won't add detail. Click to resume.",
        )
    } else {
        (
            "STOP BACKGROUND RENDERING",
            theme::DANGER,
            "Stop generating terrain tiles in the background and free the CPU. \
             A tile already being built finishes first (under a second).",
        )
    };
    let size = egui::vec2(250.0, 32.0);
    // Right-aligned, but never over the sidebar toggle at the top-left.
    let x = (rect.right() - size.x - 8.0).max(rect.left() + 52.0);
    let r = egui::Rect::from_min_size(egui::pos2(x, rect.top() + 8.0), size);
    let button = egui::Button::new(egui::RichText::new(label).strong().color(egui::Color32::WHITE)).fill(fill);
    if ui.put(r, button).on_hover_text(tip).clicked() {
        app.push(Command::ToggleTileStreaming);
    }
}
