use crate::commands::Command;
use crate::ui;
use rg_pipeline::{Equirect, Layer, World, WorldParams};
use crate::tiles::TileFeeder;
use rg_render::{EarthCamera, HeightImage, LayerImage, MapView, PlanetRenderer, RenderStyle, SelectView, ViewMode};
use rg_terrain::{TerrainSampler, TileGenerator};

/// Display resolution of baked layers: about 4 pixels per base-grid cell across the
/// equator, rounded up to a multiple of 256 and kept within a range. The web range is lower
/// because everything runs on one thread there.
fn layer_size(res: u32) -> (u32, u32) {
    #[cfg(not(target_arch = "wasm32"))]
    const RANGE: (u32, u32) = (2048, 4096);
    #[cfg(target_arch = "wasm32")]
    const RANGE: (u32, u32) = (1440, 2048);
    let width = (4 * res).next_multiple_of(256).clamp(RANGE.0, RANGE.1);
    (width, width / 2)
}

/// Display-only settings (no regeneration needed).
#[derive(Clone, Debug)]
pub struct DisplayParams {
    /// Vertical exaggeration seen from orbit. It fades to [`Self::NEAR_EXAGGERATION`] as
    /// you zoom in, so mountains look right from both orbit and the ground.
    pub exaggeration: f64,
}

impl DisplayParams {
    pub const NEAR_EXAGGERATION: f64 = 1.5;

    /// Exaggeration for a camera `distance_m` from its target, on a planet of `radius_m`.
    /// Driven by distance (not height above terrain) to avoid a feedback loop.
    pub fn effective_exaggeration(&self, distance_m: f64, radius_m: f64) -> f64 {
        let t = smoothstep(50_000.0, 0.6 * radius_m, distance_m);
        Self::NEAR_EXAGGERATION + (self.exaggeration.max(Self::NEAR_EXAGGERATION) - Self::NEAR_EXAGGERATION) * t
    }
}

fn smoothstep(a: f64, b: f64, x: f64) -> f64 {
    let t = ((x - a) / (b - a)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

impl Default for DisplayParams {
    fn default() -> Self {
        Self { exaggeration: 20.0 }
    }
}

/// Performance stats for the HUD. Two different things are measured:
///
/// - **FPS**: frames actually presented per second (wall clock). egui only repaints on
///   input or when asked, so an idle app shows "idle", not a misleading low number.
/// - **CPU frame time**: how long our own `ui()` call took (UI + render submission).
///   It's what we control; GPU time needs timestamp queries (later).
#[derive(Default)]
pub struct FrameStats {
    /// Wall-clock start times of recent frames, seconds.
    frame_starts: std::collections::VecDeque<f64>,
    /// (wall-clock time, CPU ms) of recent frames.
    cpu: std::collections::VecDeque<(f64, f32)>,
    pub bake_ms: f64,
    last_log: f64,
    pub tiles_drawn: usize,
    pub tiles_resident: usize,
    pub tiles_in_flight: usize,
    pub exaggeration: f64,
}

impl FrameStats {
    const WINDOW_S: f64 = 2.0;

    pub fn frame_started(&mut self, now: f64) {
        self.frame_starts.push_back(now);
        while self.frame_starts.front().is_some_and(|&t| now - t > Self::WINDOW_S) {
            self.frame_starts.pop_front();
        }
    }

    pub fn frame_finished(&mut self, now: f64, cpu_ms: f32) {
        self.cpu.push_back((now, cpu_ms));
        while self.cpu.front().is_some_and(|&(t, _)| now - t > Self::WINDOW_S) {
            self.cpu.pop_front();
        }
    }

    /// Presented frames per second over the last second, or None when idle (no frame in
    /// the last 0.5 s, i.e. egui isn't repainting).
    pub fn fps(&self, now: f64) -> Option<f32> {
        let last = *self.frame_starts.back()?;
        if now - last > 0.5 {
            return None;
        }
        let n = self.frame_starts.iter().filter(|&&t| now - t <= 1.0).count();
        Some(n as f32)
    }

    /// True at most once per second (debug logging).
    pub fn log_due(&mut self, now: f64) -> bool {
        if now - self.last_log >= 1.0 {
            self.last_log = now;
            return true;
        }
        false
    }

    pub fn cpu_avg_ms(&self) -> f32 {
        if self.cpu.is_empty() {
            return 0.0;
        }
        self.cpu.iter().map(|c| c.1).sum::<f32>() / self.cpu.len() as f32
    }

    pub fn cpu_worst_ms(&self) -> f32 {
        self.cpu.iter().map(|c| c.1).fold(0.0, f32::max)
    }
}

/// Seconds on a monotonic clock (works on wasm too).
pub fn now_s() -> f64 {
    thread_local! { static T0: web_time::Instant = web_time::Instant::now(); }
    T0.with(|t0| t0.elapsed().as_secs_f64())
}

pub struct RealmApp {
    pub params: WorldParams,
    pub display: DisplayParams,
    pub world: Option<World>,
    pub layer: Layer,
    pub view: ViewMode,
    pub camera: EarthCamera,
    pub map: MapView,
    pub style: RenderStyle,
    pub stats: FrameStats,
    pub adapter_info: String,
    /// A message shown in the viewport while a long step runs (it runs on the next frame,
    /// so the message gets painted first).
    pub busy: Option<&'static str>,
    pub status: String,
    /// Info under the cursor, for the readout.
    pub hover: Option<String>,
    /// True while tiles are loading or morphing in (keep repainting).
    pub streaming: bool,
    /// Whether the controls sidebar is shown. `None` until the first frame, when it
    /// defaults to hidden on narrow (phone-sized) screens.
    pub show_panel: Option<bool>,
    renderer: Option<PlanetRenderer>,
    feeder: TileFeeder,
    /// Terrain sampler for the current world (camera clearance, hover heights).
    terrain: Option<TerrainSampler>,
    equirect: Option<Equirect>,
    texture_id: Option<egui::TextureId>,
    pending: Vec<Command>,
    deferred: Vec<Command>,
}

impl RealmApp {
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        ui::theme::apply(&cc.egui_ctx);
        let mut app = Self {
            params: WorldParams::default(),
            display: DisplayParams::default(),
            world: None,
            layer: Layer::Satellite,
            view: ViewMode::Globe,
            camera: EarthCamera::new(6_371_000.0),
            map: MapView::default(),
            style: RenderStyle::default(),
            stats: FrameStats::default(),
            adapter_info: String::new(),
            busy: None,
            status: String::new(),
            hover: None,
            streaming: false,
            show_panel: None,
            renderer: None,
            feeder: TileFeeder::new(),
            terrain: None,
            equirect: None,
            texture_id: None,
            pending: vec![Command::Regenerate],
            deferred: Vec::new(),
        };
        app.apply_env_overrides();
        if let Some(rs) = cc.wgpu_render_state.as_ref() {
            let info = rs.adapter.get_info();
            app.adapter_info = format!("{} · {:?}", info.name, info.backend);
            #[cfg(not(target_arch = "wasm32"))]
            let capacity = 900;
            #[cfg(target_arch = "wasm32")]
            let capacity = 500;
            app.renderer = Some(PlanetRenderer::new(
                &rs.device,
                &rs.queue,
                &rg_terrain::tile_indices(),
                rg_terrain::TILE_CELLS,
                capacity,
            ));
        }
        app
    }

    /// Native-only startup overrides for scripted screenshots and testing:
    /// `RG_SEED=<n>`, `RG_VIEW=globe|map`, `RG_LAYER=<1-8>`, `RG_RES=<cells per face>`,
    /// camera: `RG_LAT`/`RG_LON`/`RG_HEADING`/`RG_TILT` (degrees), `RG_DIST` (metres).
    fn apply_env_overrides(&mut self) {
        #[cfg(not(target_arch = "wasm32"))]
        {
            let var = |k: &str| std::env::var(k).ok();
            if let Some(s) = var("RG_SEED").and_then(|v| v.parse().ok()) {
                self.params.tectonics.seed = s;
            }
            if let Some(r) = var("RG_RES").and_then(|v| v.parse().ok()) {
                self.params.tectonics.res = r;
            }
            let deg = |k: &str| var(k).and_then(|v| v.parse::<f64>().ok()).map(f64::to_radians);
            if let Some(v) = deg("RG_LAT") {
                self.camera.lat = v;
            }
            if let Some(v) = deg("RG_LON") {
                self.camera.lon = v;
            }
            if let Some(v) = var("RG_DIST").and_then(|v| v.parse::<f64>().ok()) {
                self.camera.distance_m = v;
            }
            if let Some(v) = deg("RG_HEADING") {
                self.camera.heading = v;
            }
            if let Some(v) = deg("RG_TILT") {
                self.camera.tilt = v.min(self.camera.max_tilt());
            }
            if var("RG_VIEW").as_deref() == Some("map") {
                self.view = ViewMode::Map;
            }
            if let Some(n) = var("RG_LAYER").and_then(|v| v.parse::<usize>().ok())
                && let Some(l) = Layer::ALL.get(n.saturating_sub(1)) {
                    self.layer = *l;
                }
        }
    }

    pub fn push(&mut self, c: Command) {
        self.pending.push(c);
    }

    fn apply(&mut self, c: Command, rs: &egui_wgpu::RenderState) {
        match c {
            // Slow commands: show the busy message this frame, run next frame.
            Command::Regenerate | Command::RandomizeSeed | Command::RerunClimate if self.busy.is_none() => {
                self.busy = Some(match c {
                    Command::RerunClimate => "Simulating climate…",
                    _ => "Generating world…",
                });
                self.deferred.push(c);
            }
            Command::Regenerate | Command::RandomizeSeed | Command::RerunClimate => self.deferred.push(c),
            Command::SetLayer(l) => {
                self.layer = l;
                self.bake_layer(rs);
            }
            Command::SetView(v) => self.view = v,
            Command::TogglePanel => self.show_panel = Some(!self.show_panel.unwrap_or(true)),
            Command::ResetCamera => {
                self.camera = EarthCamera::new(self.camera.radius_m);
                self.map = MapView::default();
            }
            Command::SaveWorld | Command::LoadWorld | Command::ExportHeightmap => {
                crate::files::handle(self, c, rs);
            }
        }
    }

    fn run_deferred(&mut self, c: Command, rs: &egui_wgpu::RenderState) {
        match c {
            Command::RandomizeSeed => {
                self.params.tectonics.seed = random_seed();
                self.regenerate(rs);
            }
            Command::Regenerate => self.regenerate(rs),
            Command::RerunClimate => {
                if let Some(w) = self.world.as_mut() {
                    rg_pipeline::rerun_climate(w, &self.params.planet);
                }
                self.after_world_changed(rs, false);
            }
            _ => {}
        }
    }

    fn regenerate(&mut self, rs: &egui_wgpu::RenderState) {
        self.world = Some(rg_pipeline::generate(&self.params));
        self.after_world_changed(rs, true);
    }

    /// Re-derive everything shown from the current world.
    pub fn after_world_changed(&mut self, rs: &egui_wgpu::RenderState, heights_changed: bool) {
        let Some(world) = self.world.as_ref() else { return };
        self.style.radius_m = world.params.tectonics.radius_m;
        if heights_changed || self.equirect.is_none() {
            let t0 = web_time::Instant::now();
            let (lw, lh) = layer_size(world.tectonics.height_m.res);
            let eq = Equirect::new(world, lw, lh);
            if let Some(r) = self.renderer.as_mut() {
                r.set_heights(
                    &rs.device,
                    &rs.queue,
                    &HeightImage { width: eq.width, height: eq.height, metres: eq.heights.clone() },
                );
            }
            self.equirect = Some(eq);
            self.stats.bake_ms = t0.elapsed().as_secs_f64() * 1000.0;
            self.start_tiles(rs);
        }
        self.bake_layer(rs);
        let w = self.world.as_ref().unwrap();
        self.status = format!(
            "Seed {} · land {:.1}% · {} plates",
            w.params.tectonics.seed,
            w.tectonics.land_fraction * 100.0,
            w.tectonics.plates.len()
        );
    }

    /// Point the tile streamer at the current world (clears the tile cache).
    fn start_tiles(&mut self, rs: &egui_wgpu::RenderState) {
        let Some(world) = self.world.as_ref() else { return };
        let radius = world.params.tectonics.radius_m;
        let sampler = rg_pipeline::terrain_sampler(world);
        if (self.camera.radius_m - radius).abs() > 1.0 {
            let scale = radius / self.camera.radius_m;
            self.camera.radius_m = radius;
            self.camera.distance_m *= scale;
        }
        if let Some(r) = self.renderer.as_mut() {
            self.feeder.set_world(TileGenerator { sampler: sampler.clone() }, &mut r.tiles, &rs.device, now_s());
        }
        self.terrain = Some(sampler);
    }

    fn bake_layer(&mut self, rs: &egui_wgpu::RenderState) {
        let (Some(world), Some(eq), Some(r)) = (self.world.as_ref(), self.equirect.as_ref(), self.renderer.as_mut()) else {
            return;
        };
        let t0 = web_time::Instant::now();
        let rgba = eq.bake(world, self.layer);
        r.set_layer(&rs.device, &rs.queue, &LayerImage { width: eq.width, height: eq.height, rgba });
        self.stats.bake_ms = t0.elapsed().as_secs_f64() * 1000.0;
        // Data layers read better without strong relief shading or glare.
        self.style.relief_shading = match self.layer {
            Layer::Satellite | Layer::Elevation => 1.0,
            _ => 0.35,
        };
        self.style.ocean_specular = self.layer == Layer::Satellite;
    }

    /// Render the planet into its offscreen texture at `size_px` and return the egui id.
    pub fn planet_texture(&mut self, rs: &egui_wgpu::RenderState, size_px: [u32; 2]) -> Option<egui::TextureId> {
        let now = now_s();
        let radius = self.camera.radius_m;
        let exag = self.display.effective_exaggeration(self.camera.distance_m, radius);
        self.stats.exaggeration = exag;
        // Orbit a point on the ground: track the displayed height under the target
        // (smoothed, so panning over ridges doesn't jolt the view).
        if let Some(t) = self.terrain.as_ref() {
            let target_dir = rg_render::earth_camera::dir(self.camera.lat, self.camera.lon);
            let h = t.height(target_dir, 50.0).0.max(0.0) * exag;
            let k = if self.camera.target_height_m == 0.0 { 1.0 } else { 0.25 };
            self.camera.target_height_m += (h - self.camera.target_height_m) * k;
        }
        // Keep the eye above the displayed terrain (checked under the eye; good enough
        // for M1a, which has no flight-through).
        let eye = self.camera.eye();
        let ground = self
            .terrain
            .as_ref()
            .map(|t| t.height(eye.normalize(), 30.0).0.max(0.0) * exag)
            .unwrap_or(0.0);
        let min_clear = 20.0;
        let clearance = eye.length() - radius - ground;
        if clearance < min_clear {
            self.camera.distance_m += min_clear - clearance;
        }
        let clearance = (self.camera.eye().length() - radius - ground).max(1.0);
        let r = self.renderer.as_mut()?;
        let recreated = r.resize(&rs.device, size_px);
        if self.view == ViewMode::Globe {
            let aspect = size_px[0] as f64 / size_px[1].max(1) as f64;
            let sel = r.tiles.select(&SelectView {
                eye: self.camera.eye(),
                view_proj: self.camera.view_proj_f64(aspect, clearance),
                radius_m: radius,
                exaggeration: exag,
                viewport_h_px: size_px[1] as f64,
                fov_y: self.camera.fov_y,
                tau_px: 2.0,
                max_cell_px: r.max_cell_px,
                max_level: 14,
                now_s: now,
            });
            self.feeder.update(&sel.wanted, &mut r.tiles, &rs.device, now);
            self.stats.tiles_drawn = sel.drawn;
            self.stats.tiles_resident = r.tiles.len();
            self.stats.tiles_in_flight = self.feeder.in_flight();
            self.streaming = sel.animating || !sel.wanted.is_empty() || self.feeder.in_flight() > 0;
            #[cfg(not(target_arch = "wasm32"))]
            if std::env::var_os("RG_LOG_TILES").is_some() && self.stats.log_due(now) {
                eprintln!(
                    "[tiles] drawn {} cached {} wanted {} in-flight {} alt {:.0} m exag {:.2} gen {}",
                    sel.drawn,
                    r.tiles.len(),
                    sel.wanted.len(),
                    self.feeder.in_flight(),
                    self.camera.altitude_m(),
                    exag,
                    self.feeder.last_generated
                );
            }
        } else {
            self.streaming = false;
        }
        r.render(&rs.device, &rs.queue, self.view, &self.camera, clearance, exag, &self.map, &self.style);
        let view = r.output_view()?;
        let mut egui_r = rs.renderer.write();
        match self.texture_id {
            Some(id) if recreated => {
                egui_r.update_egui_texture_from_wgpu_texture(&rs.device, view, wgpu::FilterMode::Linear, id)
            }
            Some(_) => {}
            None => self.texture_id = Some(egui_r.register_native_texture(&rs.device, view, wgpu::FilterMode::Linear)),
        }
        self.texture_id
    }

    /// Describe the point under the cursor (lat/lon radians).
    pub fn describe(&self, lat: f64, lon: f64) -> Option<String> {
        let w = self.world.as_ref()?;
        let d = rg_pipeline::layers::dir(lat, lon);
        let h = w.tectonics.height_m.sample_drained(d) as f32;
        let s = rg_pipeline::layers::ClimateSample::at(w, lat, lon, h);
        let ns = if lat >= 0.0 { 'N' } else { 'S' };
        let ew = if lon >= 0.0 { 'E' } else { 'W' };
        let mut out = format!(
            "{:.2}°{ns} {:.2}°{ew} · {} {:.0} m",
            lat.to_degrees().abs(),
            lon.to_degrees().abs(),
            if h >= 0.0 { "elev" } else { "depth" },
            h.abs()
        );
        let mean_t = s.t.iter().sum::<f32>() / 12.0;
        let p: f32 = s.p.iter().sum();
        out += &format!(" · {mean_t:.1} °C · {p:.0} mm/yr");
        if h >= 0.0 {
            out += &format!(" · {}", s.koppen(lat).code());
        }
        Some(out)
    }
}

impl eframe::App for RealmApp {
    fn ui(&mut self, ui: &mut egui::Ui, frame: &mut eframe::Frame) {
        let frame_t0 = now_s();
        self.stats.frame_started(frame_t0);
        let Some(rs) = frame.wgpu_render_state().cloned() else {
            ui.label("RealmGenesis needs a GPU (wgpu) renderer.");
            return;
        };
        // Deferred (slow) commands queued last frame run now, after their busy message
        // has been painted once.
        let deferred = std::mem::take(&mut self.deferred);
        let had_deferred = !deferred.is_empty();
        for c in deferred {
            self.run_deferred(c, &rs);
        }
        if had_deferred {
            self.busy = None;
        }
        for c in std::mem::take(&mut self.pending) {
            self.apply(c, &rs);
        }
        if !self.deferred.is_empty() {
            ui.ctx().request_repaint();
        }
        // Presentation light: the sun sits up and to the left of the camera, so the view
        // always shows relief and a terminator. Real seasonal sun comes with climate.
        let (up, east, north) = self.camera.frame();
        let (sh, ch) = self.camera.heading.sin_cos();
        // Light from the upper left of the screen, wherever the camera faces.
        let screen_left = -(east * ch - north * sh);
        let screen_up = north * ch + east * sh;
        self.style.sun_dir = up * 1.2 + screen_left * 0.55 + screen_up * 0.7;
        ui::draw(self, ui, &rs);
        let end = now_s();
        self.stats.frame_finished(end, ((end - frame_t0) * 1000.0) as f32);
    }
}

fn random_seed() -> u64 {
    // Time-based on purpose: this only picks a new seed. Generation itself is deterministic.
    let t = web_time::SystemTime::now().duration_since(web_time::UNIX_EPOCH).unwrap_or_default();
    (t.as_nanos() as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) % 1_000_000
}
