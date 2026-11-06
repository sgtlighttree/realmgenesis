use crate::commands::Command;
use crate::console::Console;
use crate::worldgen::{GenJob, Poll};
use crate::ui;
use rg_pipeline::{Equirect, Layer, World, WorldParams};
use crate::tiles::TileFeeder;
use rg_render::{EarthCamera, GradientImage, HeightImage, LayerImage, MapView, PlanetRenderer, RenderStyle, SelectView, ViewMode};
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
    /// Keep the orbit exaggeration all the way down instead of fading it (Matt, s02).
    pub keep_relief: bool,
    /// Zoom floor: how many screen pixels one texel of the stored surface may span before
    /// the camera stops descending (s03; the data is finite, so infinite zoom only shows
    /// blur). Higher = closer.
    pub max_px_per_texel: f64,
    /// Show the colour legend for the active layer.
    pub legend: bool,
}

impl DisplayParams {
    pub const NEAR_EXAGGERATION: f64 = 1.5;

    /// Exaggeration for a camera `distance_m` from its target, on a planet of `radius_m`.
    /// Driven by distance (not height above terrain) to avoid a feedback loop.
    pub fn effective_exaggeration(&self, distance_m: f64, radius_m: f64) -> f64 {
        if self.keep_relief {
            return self.exaggeration;
        }
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
        // 35: Matt's pick on the M1 (s02); 20 made network valleys read as flat grooves.
        // 32 px per texel: ~1300 km above an Ultra world on a 1600 px tall view (s03 start
        // value; Matt tunes it).
        Self { exaggeration: 35.0, keep_relief: false, max_px_per_texel: 32.0, legend: true }
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
    /// True while the orbit point is still easing onto the ground height (keep repainting).
    pub settling: bool,
    /// Whether the controls sidebar is shown. `None` until the first frame, when it
    /// defaults to hidden on narrow (phone-sized) screens.
    pub show_panel: Option<bool>,
    /// What the app has been doing (the status-bar console).
    pub console: Console,
    pub console_open: bool,
    /// Coast source / import / brush UI state (R3).
    pub coasts: crate::ui::coasts::CoastUi,
    /// The surface-gradient bake in flight (per-pixel lighting; see `normals.rs`).
    #[cfg(not(target_arch = "wasm32"))]
    normals_job: Option<crate::normals::NormalsJob>,
    /// Deepest tile level worth drawing: its cells are half a macro cell (finer would only
    /// re-sample the same smooth surface; lighting comes from the gradient texture).
    max_tile_level: u8,
    /// Size of one texel of the stored surface (m), for the zoom floor. 0 = no world yet.
    surface_texel_m: f64,
    /// The zoom floor (camera distance, m) used last frame, for the HUD and settings.
    pub zoom_floor_m: f64,
    /// A screenshot to save after the next render: (path, wait until the world, normals and
    /// tiles have all settled). F12, or `RG_SHOT=<path>` at startup for scripted looks.
    pending_shot: Option<(std::path::PathBuf, bool)>,
    /// Quit after the scripted screenshot (`RG_SHOT_EXIT=1`).
    exit_after_shot: bool,
    /// `RG_SHOT_UI=1`: the scripted screenshot is the whole window (panels, overlays and all)
    /// instead of just the planet view. macOS blocks `screencapture` here, so this is how a
    /// session checks the UI.
    shot_ui: bool,
    /// A whole-window screenshot to request this frame, and the flag its callback sets.
    ui_shot_due: Option<std::path::PathBuf>,
    ui_shot_done: std::sync::Arc<std::sync::atomic::AtomicBool>,
    close_requested: bool,
    /// World generation running on a background thread (native only; see `worldgen.rs`).
    gen_job: Option<GenJob>,
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
            settling: false,
            show_panel: None,
            console: Console::new(),
            console_open: false,
            coasts: Default::default(),
            #[cfg(not(target_arch = "wasm32"))]
            normals_job: None,
            max_tile_level: 4,
            surface_texel_m: 0.0,
            zoom_floor_m: 0.0,
            pending_shot: None,
            exit_after_shot: false,
            close_requested: false,
            shot_ui: false,
            ui_shot_due: None,
            ui_shot_done: Default::default(),
            gen_job: None,
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
    /// `RG_SEED=<n>`, `RG_VIEW=globe|map`, `RG_LAYER=<1-9>`, `RG_RES=<cells per face>`,
    /// camera: `RG_LAT`/`RG_LON`/`RG_HEADING`/`RG_TILT` (degrees), `RG_DIST` (metres);
    /// `RG_EXAG=<×>` and `RG_KEEP_RELIEF=1` (relief); `RG_CONV_SMOOTH=<km>` (uplift blend); `RG_SHOT=<png>` saves the planet view once the
    /// world, normals and tiles have settled (`RG_SHOT_UI=1`: the whole window), and `RG_SHOT_EXIT=1` then quits.
    fn apply_env_overrides(&mut self) {
        #[cfg(not(target_arch = "wasm32"))]
        {
            let var = |k: &str| std::env::var(k).ok();
            if let Some(s) = var("RG_SEED").and_then(|v| v.parse().ok()) {
                self.params.tectonics.seed = s;
            }
            if let Some(v) = var("RG_EXAG").and_then(|v| v.parse().ok()) {
                self.display.exaggeration = v;
            }
            if var("RG_KEEP_RELIEF").is_some() {
                self.display.keep_relief = true;
            }
            if let Some(v) = var("RG_PLATEAU").and_then(|v| v.parse().ok()) {
                self.params.erosion.plateau_runoff_mm = v;
            }
            if let Some(v) = var("RG_MATURITY").and_then(|v| v.parse().ok()) {
                self.params.erosion.duration_myr = v;
            }
            if let Some(v) = var("RG_JITTER").and_then(|v| v.parse().ok()) {
                self.params.erosion.node_jitter = v;
            }
            if let Some(v) = var("RG_ROUGH_WL").and_then(|v| v.parse().ok()) {
                self.params.erosion.roughness_wavelength_km = v;
            }
            if let Some(v) = var("RG_CONV_SMOOTH").and_then(|v| v.parse().ok()) {
                self.params.erosion.convergence_smoothing_km = v;
            }
            if let Some(p) = var("RG_SHOT") {
                self.pending_shot = Some((p.into(), true));
                self.exit_after_shot = var("RG_SHOT_EXIT").is_some();
                self.shot_ui = var("RG_SHOT_UI").is_some();
            }
            if let Some(r) = var("RG_RES").and_then(|v| v.parse().ok()) {
                self.params.tectonics.res = r;
            }
            // R3: `RG_ACTIVE_SHARE=<0..1>`, `RG_STRICT=1`, `RG_MASK=<png|svg>` (whole world, or a
            // patch with `RG_MASK_PATCH=lat,lon,km_per_px,rotation`).
            if let Some(v) = var("RG_ACTIVE_SHARE").and_then(|v| v.parse().ok()) {
                self.params.tectonics.active_share = v;
            }
            if var("RG_STRICT").is_some() {
                self.params.tectonics.conditioning.strict = true;
            }
            if let Some(path) = var("RG_MASK") {
                let place = match var("RG_MASK_PATCH").map(|v| v.split(',').filter_map(|x| x.trim().parse::<f64>().ok()).collect::<Vec<_>>()) {
                    Some(v) if v.len() >= 3 => rg_pipeline::mask::Placement::Patch {
                        center_deg: [v[0], v[1]],
                        km_per_px: v[2],
                        rotation_deg: v.get(3).copied().unwrap_or(0.0),
                    },
                    _ => rg_pipeline::mask::Placement::WorldEquirect,
                };
                match rg_pipeline::import::load_mask(std::path::Path::new(&path), place.clone()) {
                    Ok(m) => {
                        // Keep the coasts panel in step with the mask.
                        let st = &mut self.coasts;
                        st.file_name = std::path::Path::new(&path).file_name().map_or(path.clone(), |n| n.to_string_lossy().into_owned());
                        if let rg_pipeline::mask::Placement::Patch { center_deg, km_per_px, rotation_deg } = place {
                            (st.patch, st.lat, st.lon, st.km_per_px, st.rotation) = (true, center_deg[0], center_deg[1], km_per_px, rotation_deg);
                            // `RG_PLACE=1`: start in "Move patch" (handles shown), for screenshots.
                            st.placing = var("RG_PLACE").is_some();
                        }
                        self.params.tectonics.conditioning.mask = Some(std::sync::Arc::new(m));
                        self.params.tectonics.conditioning.strict = true;
                    }
                    Err(e) => eprintln!("RG_MASK: {e}"),
                }
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

    /// Whether background tile generation is stopped (the red button).
    pub fn tiles_paused(&self) -> bool {
        self.feeder.paused()
    }

    /// Add a line to the console (and stderr natively).
    pub fn log(&mut self, msg: impl Into<String>) {
        self.console.log(now_s(), msg);
    }

    /// True while a world is being generated in the background (native).
    pub fn generating(&self) -> bool {
        self.gen_job.is_some()
    }

    /// A live one-line description of ongoing work for the status bar, or `None` when idle.
    /// The spinner frame and elapsed time come from the clock, so the line visibly ticks even
    /// inside a long step.
    pub fn activity(&self, now_s: f64) -> Option<String> {
        const SPIN: [char; 4] = ['|', '/', '-', '\\'];
        let frame = SPIN[((now_s * 8.0) as usize) % SPIN.len()];
        if let Some(job) = &self.gen_job {
            let what = match job.latest {
                Some(p) if p.steps > 0 => format!("{} step {}/{}", p.stage.label(), p.step, p.steps),
                Some(p) => p.stage.label().to_string(),
                None => "starting".to_string(),
            };
            let elapsed = now_s - job.started_s;
            if job.cancelling() {
                return Some(format!("{frame} Cancelling · finishing {what} · {elapsed:.1} s"));
            }
            return Some(format!("{frame} Generating world · {what} · {elapsed:.1} s · Esc to cancel"));
        }
        if self.view == ViewMode::Globe && self.stats.tiles_in_flight > 0 {
            return Some(format!(
                "{frame} Streaming terrain · loading {} · {} drawn",
                self.stats.tiles_in_flight, self.stats.tiles_drawn
            ));
        }
        None
    }

    pub fn push(&mut self, c: Command) {
        self.pending.push(c);
    }

    fn apply(&mut self, c: Command, rs: &egui_wgpu::RenderState) {
        match c {
            // Enforced here, not only by disabled buttons: the G/R shortcuts bypass buttons.
            Command::Regenerate | Command::RandomizeSeed | Command::RerunClimate | Command::LoadWorld
                if self.generating() =>
            {
                self.log("Busy generating a world; that will have to wait until it's done.");
            }
            Command::Regenerate | Command::RandomizeSeed => self.start_generation(c),
            // Slow commands: show the busy message this frame, run next frame.
            Command::RerunClimate if self.busy.is_none() => {
                self.busy = Some("Simulating climate…");
                self.deferred.push(c);
            }
            Command::RerunClimate => self.deferred.push(c),
            Command::SetLayer(l) => {
                self.layer = l;
                self.bake_layer(rs);
            }
            Command::SetView(v) => self.view = v,
            Command::TogglePanel => self.show_panel = Some(!self.show_panel.unwrap_or(true)),
            Command::ToggleTileStreaming => {
                self.feeder.set_paused(!self.feeder.paused());
                self.log(if self.feeder.paused() { "Background rendering stopped." } else { "Background rendering resumed." });
            }
            Command::ToggleConsole => self.console_open = !self.console_open,
            Command::CancelGeneration => {
                if let Some(job) = self.gen_job.as_ref().filter(|j| !j.cancelling()) {
                    job.cancel();
                    self.log("Cancelling: generation stops at the next step.");
                }
            }
            Command::ResetCamera => {
                self.camera = EarthCamera::new(self.camera.radius_m);
                self.map = MapView::default();
            }
            Command::SaveScreenshot => {
                #[cfg(not(target_arch = "wasm32"))]
                match crate::files::shot_path() {
                    Ok(p) => self.pending_shot = Some((p, false)),
                    Err(e) => self.log(format!("⚠ {e}")),
                }
            }
            Command::SaveWorld | Command::LoadWorld | Command::ExportHeightmap => {
                crate::files::handle(self, c, rs);
            }
        }
    }

    /// Start generating a world from the current parameters. Native: on a background thread,
    /// with progress in the status bar (`poll_generation` picks the result up). Web: no
    /// threads, so the old path (busy message this frame, generate next frame).
    fn start_generation(&mut self, c: Command) {
        if matches!(c, Command::RandomizeSeed) {
            self.params.tectonics.seed = random_seed();
        }
        let (seed, res) = (self.params.tectonics.seed, self.params.tectonics.res);
        #[cfg(not(target_arch = "wasm32"))]
        match GenJob::spawn(self.params.clone(), now_s()) {
            Ok(job) => {
                self.log(format!("Generating world: seed {seed}, {res} cells per face"));
                self.gen_job = Some(job);
            }
            Err(e) => self.log(format!("⚠ Couldn't start generation: {e}")),
        }
        #[cfg(target_arch = "wasm32")]
        {
            self.log(format!("Generating world: seed {seed}, {res} cells per face"));
            if self.busy.is_none() {
                self.busy = Some("Generating world…");
            }
            self.deferred.push(Command::Regenerate);
        }
    }

    /// Pick up progress and the finished world from the background job (native).
    fn poll_generation(&mut self, rs: &egui_wgpu::RenderState, ctx: &egui::Context) {
        let Some(job) = self.gen_job.as_mut() else { return };
        let (stages, state) = job.poll();
        let started = job.started_s;
        for p in stages {
            self.log(format!("  {} (+{:.1} s)", p.stage.label(), now_s() - started));
        }
        match state {
            // Tick the status line at 10 Hz, not 60: nothing else changes meanwhile.
            Poll::Working => ctx.request_repaint_after(std::time::Duration::from_millis(100)),
            Poll::Done(world) => {
                self.gen_job = None;
                let t = world.timings.clone();
                self.world = Some(*world);
                self.after_world_changed(rs, true);
                self.log(format!(
                    "World ready in {:.1} s (tectonics {:.0} ms, erosion {:.0} ms, climate {:.0} ms)",
                    now_s() - started,
                    t.tectonics_ms,
                    t.erosion_ms,
                    t.climate_ms
                ));
            }
            Poll::Cancelled => {
                self.gen_job = None;
                self.log(format!("Generation cancelled after {:.1} s. The previous world is still shown.", now_s() - started));
            }
            Poll::Failed => {
                self.gen_job = None;
                self.log("⚠ World generation failed (its thread stopped). The previous world is still shown.");
            }
        }
    }

    fn run_deferred(&mut self, c: Command, rs: &egui_wgpu::RenderState) {
        match c {
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
        // R2 (s03): the display is the committed macro surface (smooth bicubic), lit per
        // pixel from the baked gradient. The streamed detail network is gone from the app.
        let sampler = rg_pipeline::terrain_sampler(world);
        let macro_res = world.tectonics.height_m.res;
        // Tile cells of half a macro cell: 64 · 2^L ≥ 2 · res.
        self.max_tile_level = ((2.0 * macro_res as f64 / rg_terrain::TILE_CELLS as f64).log2().ceil().max(2.0)) as u8;
        self.surface_texel_m =
            std::f64::consts::FRAC_PI_2 * radius / rg_terrain::gradient::gradient_res(macro_res) as f64;
        if (self.camera.radius_m - radius).abs() > 1.0 {
            let scale = radius / self.camera.radius_m;
            self.camera.radius_m = radius;
            self.camera.distance_m *= scale;
        }
        if let Some(r) = self.renderer.as_mut() {
            self.feeder.set_world(TileGenerator { sampler: sampler.clone() }, &mut r.tiles, &rs.device, now_s());
        }
        self.start_normals(rs, &sampler);
        self.terrain = Some(sampler);
    }

    /// Reset the lighting gradient to flat and bake the new world's (in the background
    /// natively).
    fn start_normals(&mut self, rs: &egui_wgpu::RenderState, sampler: &TerrainSampler) {
        if let Some(r) = self.renderer.as_mut() {
            r.set_gradient(&rs.device, &rs.queue, &GradientImage { res: 1, data: vec![[0.0; 2]; 24] });
        }
        #[cfg(not(target_arch = "wasm32"))]
        match crate::normals::NormalsJob::spawn(sampler.clone()) {
            Ok(job) => {
                self.normals_job = Some(job);
                self.log("  Baking surface normals…");
            }
            Err(e) => self.log(format!("⚠ Couldn't start the normals bake: {e}")),
        }
        #[cfg(target_arch = "wasm32")]
        {
            let img = crate::normals::bake(sampler);
            if let Some(r) = self.renderer.as_mut() {
                r.set_gradient(&rs.device, &rs.queue, &img);
            }
        }
    }

    /// Upload the gradient once its bake finishes (native).
    #[cfg(not(target_arch = "wasm32"))]
    fn poll_normals(&mut self, rs: &egui_wgpu::RenderState, ctx: &egui::Context) {
        use crate::normals::NormalsPoll;
        let Some(job) = self.normals_job.as_ref() else { return };
        match job.poll() {
            NormalsPoll::Working => ctx.request_repaint_after(std::time::Duration::from_millis(100)),
            NormalsPoll::Done(img, secs) => {
                self.normals_job = None;
                if let Some(r) = self.renderer.as_mut() {
                    r.set_gradient(&rs.device, &rs.queue, &img);
                }
                self.log(format!("  Surface normals baked in {secs:.1} s ({} per face)", img.res));
                ctx.request_repaint();
            }
            NormalsPoll::Failed => {
                self.normals_job = None;
                self.log("⚠ The normals bake stopped (its thread ended). Lighting stays flat.");
            }
        }
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
        self.style.river_rgb = rg_pipeline::layers::river_color(self.layer);
    }

    /// Render the planet into its offscreen texture at `size_px` and return the egui id.
    pub fn planet_texture(&mut self, rs: &egui_wgpu::RenderState, size_px: [u32; 2]) -> Option<egui::TextureId> {
        let now = now_s();
        let radius = self.camera.radius_m;
        let exag = self.display.effective_exaggeration(self.camera.distance_m, radius);
        self.stats.exaggeration = exag;
        let ground_at = |t: &TerrainSampler, dir: glam::DVec3| t.macro_height(dir).max(0.0) * exag;
        // Orbit a point on the ground: track the displayed height under the target
        // (smoothed, so panning over ridges doesn't jolt the view).
        if let Some(t) = self.terrain.as_ref() {
            let target_dir = rg_render::earth_camera::dir(self.camera.lat, self.camera.lon);
            let h = ground_at(t, target_dir);
            let k = if self.camera.target_height_m == 0.0 { 1.0 } else { 0.25 };
            self.camera.target_height_m += (h - self.camera.target_height_m) * k;
            self.settling = (h - self.camera.target_height_m).abs() > 0.5;
        }
        // Keep the eye above the displayed terrain (checked under the eye; good enough
        // for M1a, which has no flight-through).
        let eye = self.camera.eye();
        let ground = self.terrain.as_ref().map(|t| ground_at(t, eye.normalize())).unwrap_or(0.0);
        let min_clear = 20.0;
        let clearance = eye.length() - radius - ground;
        if clearance < min_clear {
            self.camera.distance_m += min_clear - clearance;
        }
        let clearance = (self.camera.eye().length() - radius - ground).max(1.0);
        // Zoom floor: stop where one texel of the stored surface spans `max_px_per_texel`
        // screen pixels (eye-to-target distance, so tilt doesn't change it).
        if self.surface_texel_m > 0.0 {
            let focal_px = size_px[1].max(1) as f64 / (2.0 * (0.5 * self.camera.fov_y).tan());
            self.zoom_floor_m = self.surface_texel_m * focal_px / self.display.max_px_per_texel.max(1.0);
            if self.camera.distance_m < self.zoom_floor_m {
                self.camera.distance_m = self.zoom_floor_m;
                self.camera.tilt = self.camera.tilt.min(self.camera.max_tilt());
            }
        }
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
                max_level: self.max_tile_level,
                now_s: now,
            });
            // While a new world is generating, the old one's tiles wait (CPU for generation).
            let generating = self.gen_job.is_some();
            if !generating {
                self.feeder.update(&sel.wanted, &mut r.tiles, &rs.device, now);
            }
            self.stats.tiles_drawn = sel.drawn;
            self.stats.tiles_resident = r.tiles.len();
            self.stats.tiles_in_flight = self.feeder.in_flight();
            // While stopped, wanted tiles never arrive: repainting for them would keep a
            // 60 fps loop running, so only morph-ins of resident tiles count then.
            let loading =
                !generating && !self.feeder.paused() && (!sel.wanted.is_empty() || self.feeder.in_flight() > 0);
            self.streaming = sel.animating || loading;
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
        #[cfg(not(target_arch = "wasm32"))]
        if let Some((path, wait)) = self.pending_shot.clone() {
            let settled = self.world.is_some() && self.gen_job.is_none() && self.normals_job.is_none() && !self.streaming;
            if (!wait || settled) && self.shot_ui {
                self.pending_shot = None;
                self.ui_shot_due = Some(path);
            } else if !wait || settled {
                self.pending_shot = None;
                let msg = match r.read_output(&rs.device, &rs.queue) {
                    Some((w, h, rgba)) => match crate::files::write_png_rgba(&path, w, h, &rgba) {
                        Ok(()) => format!("Screenshot saved: {}", path.display()),
                        Err(e) => format!("⚠ screenshot: {e}"),
                    },
                    None => "⚠ screenshot: nothing rendered yet".into(),
                };
                eprintln!("{msg}");
                // (`self.log` would borrow all of self while `r` holds the renderer.)
                self.console.log(now_s(), msg);
                self.close_requested = self.exit_after_shot;
            }
        }
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
        self.poll_generation(&rs, ui.ctx());
        #[cfg(not(target_arch = "wasm32"))]
        self.poll_normals(&rs, ui.ctx());
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
        // Presentation light: the cartographic convention, from the north-west of the point
        // the camera orbits, 40° above the horizon. Tied to the ground, not the camera, so
        // turning the view doesn't relight the relief (Matt, s03: the old camera-relative,
        // mostly overhead light shifted as he rotated and flattened the contrast). Real
        // seasonal sun comes with climate.
        let (up, east, north) = self.camera.frame();
        let north_west = (north - east).normalize();
        let elevation = 40f64.to_radians();
        self.style.sun_dir = up * elevation.sin() + north_west * elevation.cos();
        ui::draw(self, ui, &rs);
        #[cfg(not(target_arch = "wasm32"))]
        if let Some(path) = self.ui_shot_due.take() {
            let (done, ctx) = (self.ui_shot_done.clone(), ui.ctx().clone());
            ui.ctx().send_viewport_cmd(egui::ViewportCommand::Screenshot(egui::ScreenshotCallback::new(move |img| {
                let [w, h] = img.size;
                match crate::files::write_png_rgba(&path, w as u32, h as u32, img.as_raw()) {
                    Ok(()) => eprintln!("Window screenshot saved: {}", path.display()),
                    Err(e) => eprintln!("⚠ window screenshot: {e}"),
                }
                done.store(true, std::sync::atomic::Ordering::Relaxed);
                ctx.request_repaint();
            })));
            ui.ctx().request_repaint();
        }
        if self.ui_shot_done.swap(false, std::sync::atomic::Ordering::Relaxed) {
            self.close_requested = self.exit_after_shot;
        }
        if self.close_requested {
            ui.ctx().send_viewport_cmd(egui::ViewportCommand::Close);
        }
        let end = now_s();
        self.stats.frame_finished(end, ((end - frame_t0) * 1000.0) as f32);
    }
}

fn random_seed() -> u64 {
    // Time-based on purpose: this only picks a new seed. Generation itself is deterministic.
    let t = web_time::SystemTime::now().duration_since(web_time::UNIX_EPOCH).unwrap_or_default();
    (t.as_nanos() as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) % 1_000_000
}
