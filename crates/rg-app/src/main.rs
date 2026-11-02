//! RealmGenesis desktop + web entry points. All app logic lives in `app.rs`.

mod app;
mod commands;
mod files;
mod tiles;
mod ui;

#[cfg(not(target_arch = "wasm32"))]
fn main() -> eframe::Result {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("warn")).init();
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("RealmGenesis")
            .with_inner_size([1440.0, 900.0])
            .with_min_inner_size([800.0, 500.0]),
        renderer: eframe::Renderer::Wgpu,
        ..Default::default()
    };
    eframe::run_native("RealmGenesis", options, Box::new(|cc| Ok(Box::new(app::RealmApp::new(cc)))))
}

#[cfg(target_arch = "wasm32")]
fn main() {
    use wasm_bindgen::JsCast as _;
    eframe::WebLogger::init(log::LevelFilter::Warn).ok();
    wasm_bindgen_futures::spawn_local(async {
        let document = web_sys::window().expect("no window").document().expect("no document");
        let canvas = document
            .get_element_by_id("realmgenesis")
            .expect("missing #realmgenesis canvas")
            .dyn_into::<web_sys::HtmlCanvasElement>()
            .expect("#realmgenesis is not a canvas");
        let result = eframe::WebRunner::new()
            .start(canvas, eframe::WebOptions::default(), Box::new(|cc| Ok(Box::new(app::RealmApp::new(cc)))))
            .await;
        if let Some(el) = document.get_element_by_id("loading") {
            match result {
                Ok(()) => el.remove(),
                Err(e) => el.set_inner_html(&format!("<p>RealmGenesis couldn't start: {e:?}</p><p>It needs a browser with WebGPU.</p>")),
            }
        }
    });
}
