//! Visual theme. Every colour and spacing choice lives here, so a redesign starts with
//! this one file.

use egui::{Color32, CornerRadius, FontFamily, FontId, Stroke, TextStyle};

pub const BG: Color32 = Color32::from_rgb(14, 17, 20);
pub const PANEL: Color32 = Color32::from_rgb(21, 25, 29);
pub const SURFACE: Color32 = Color32::from_rgb(30, 35, 40);
pub const LINE: Color32 = Color32::from_rgb(44, 51, 58);
pub const TEXT: Color32 = Color32::from_rgb(218, 224, 228);
pub const MUTED: Color32 = Color32::from_rgb(132, 144, 152);
/// Patina green: a nod to old map inks.
pub const ACCENT: Color32 = Color32::from_rgb(110, 186, 154);

pub fn apply(ctx: &egui::Context) {
    let mut style = (*ctx.global_style()).clone();
    let v = &mut style.visuals;
    *v = egui::Visuals::dark();
    v.panel_fill = PANEL;
    v.window_fill = PANEL;
    v.extreme_bg_color = BG;
    v.faint_bg_color = SURFACE;
    v.override_text_color = Some(TEXT);
    v.selection.bg_fill = ACCENT.linear_multiply(0.35);
    v.selection.stroke = Stroke::new(1.0, ACCENT);
    v.hyperlink_color = ACCENT;
    let r = CornerRadius::same(5);
    for w in [
        &mut v.widgets.noninteractive,
        &mut v.widgets.inactive,
        &mut v.widgets.hovered,
        &mut v.widgets.active,
        &mut v.widgets.open,
    ] {
        w.corner_radius = r;
    }
    v.widgets.noninteractive.bg_stroke = Stroke::new(1.0, LINE);
    v.widgets.inactive.bg_fill = SURFACE;
    v.widgets.inactive.weak_bg_fill = SURFACE;
    v.widgets.hovered.bg_stroke = Stroke::new(1.0, ACCENT.linear_multiply(0.6));
    v.widgets.active.bg_stroke = Stroke::new(1.0, ACCENT);

    style.spacing.item_spacing = egui::vec2(8.0, 7.0);
    style.spacing.button_padding = egui::vec2(10.0, 5.0);
    style.spacing.slider_width = 150.0;
    style.text_styles.insert(TextStyle::Heading, FontId::new(17.0, FontFamily::Proportional));
    style.text_styles.insert(TextStyle::Body, FontId::new(13.5, FontFamily::Proportional));
    style.text_styles.insert(TextStyle::Button, FontId::new(13.5, FontFamily::Proportional));
    style.text_styles.insert(TextStyle::Small, FontId::new(11.0, FontFamily::Proportional));
    style.text_styles.insert(TextStyle::Monospace, FontId::new(12.0, FontFamily::Monospace));
    ctx.set_global_style(style);
}

/// A small uppercase section label.
pub fn section(ui: &mut egui::Ui, text: &str) {
    ui.add_space(6.0);
    ui.label(egui::RichText::new(text.to_uppercase()).small().color(MUTED).extra_letter_spacing(0.8));
}
