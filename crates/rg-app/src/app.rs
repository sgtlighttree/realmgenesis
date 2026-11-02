use crate::commands::Command;
use crate::ui;
use rg_pipeline::{Equirect, Layer, World, WorldParams};
use rg_render::{HeightImage, LayerImage, MapView, OrbitCamera, PlanetMesh, PlanetRenderer, RenderStyle, ViewMode};

/// Display resolution of baked layers. Lower on the web, where everything runs on one
/// thread.
#[cfg(not(target_arch = "wasm32"))]
const LAYER_SIZE: (u32, u32) = (2048, 1024);
#[cfg(target_arch = "wasm32")]
const LAYER_SIZE: (u32, u32) = (1440, 720);

/// Display-only settings (no regeneration needed).
#[derive(Clone, Debug)]
pub struct DisplayParams {
    pub exaggeration: f64,
    pub mesh_res: u32,
}

impl Default for DisplayParams {
    fn default() -> Self {
        Self { exaggeration: 25.0, mesh_res: 192 }
    }
}

/// Rolling frame-time stats for the HUD.
#[derive(Default)]
pub struct FrameStats {
    samples: std::collections::VecDeque<f32>,
    pub bake_ms: f64,
    pub mesh_ms: f64,
    pub triangles: usize,
}

impl FrameStats {
    pub fn push(&mut self, dt: f32) {
        self.samples.push_back(dt);
        if self.samples.len() > 120 {
            self.samples.pop_front();
        }
    }
    pub fn avg_ms(&self) -> f32 {
        if self.samples.is_empty() {
            return 0.0;
        }
        1000.0 * self.samples.iter().sum::<f32>() / self.samples.len() as f32
    }
    pub fn worst_ms(&self) -> f32 {
        1000.0 * self.samples.iter().copied().fold(0.0, f32::max)
    }
}

pub struct RealmApp {
    pub params: WorldParams,
    pub display: DisplayParams,
    pub world: Option<World>,
    pub layer: Layer,
    pub view: ViewMode,
    pub camera: OrbitCamera,
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
    renderer: Option<PlanetRenderer>,
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
            camera: OrbitCamera::default(),
            map: MapView::default(),
            style: RenderStyle::default(),
            stats: FrameStats::default(),
            adapter_info: String::new(),
            busy: None,
            status: String::new(),
            hover: None,
            renderer: None,
            equirect: None,
            texture_id: None,
            pending: vec![Command::Regenerate],
            deferred: Vec::new(),
        };
        app.apply_env_overrides();
        if let Some(rs) = cc.wgpu_render_state.as_ref() {
            let info = rs.adapter.get_info();
            app.adapter_info = format!("{} · {:?}", info.name, info.backend);
            app.renderer = Some(PlanetRenderer::new(&rs.device, &rs.queue));
        }
        app
    }

    /// Native-only startup overrides for scripted screenshots and testing:
    /// `RG_SEED=<n>`, `RG_VIEW=globe|map`, `RG_LAYER=<1-8>`, `RG_RES=<cells per face>`.
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
            Command::RebuildMesh => self.rebuild_mesh(rs),
            Command::SetView(v) => self.view = v,
            Command::ResetCamera => {
                self.camera = OrbitCamera::default();
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
            let eq = Equirect::new(world, LAYER_SIZE.0, LAYER_SIZE.1);
            if let Some(r) = self.renderer.as_mut() {
                r.set_heights(
                    &rs.device,
                    &rs.queue,
                    &HeightImage { width: eq.width, height: eq.height, metres: eq.heights.clone() },
                );
            }
            self.equirect = Some(eq);
            self.stats.bake_ms = t0.elapsed().as_secs_f64() * 1000.0;
            self.rebuild_mesh(rs);
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

    fn rebuild_mesh(&mut self, rs: &egui_wgpu::RenderState) {
        let Some(world) = self.world.as_ref() else { return };
        let t0 = web_time::Instant::now();
        let field = &world.tectonics.height_m;
        let mesh = PlanetMesh::build(self.display.mesh_res, world.params.tectonics.radius_m, self.display.exaggeration, |d| {
            field.sample(d) as f64
        });
        self.stats.triangles = mesh.indices.len() / 3;
        if let Some(r) = self.renderer.as_mut() {
            r.upload_mesh(&rs.device, &mesh);
        }
        self.stats.mesh_ms = t0.elapsed().as_secs_f64() * 1000.0;
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
        let r = self.renderer.as_mut()?;
        let recreated = r.resize(&rs.device, size_px);
        r.render(&rs.device, &rs.queue, self.view, &self.camera, &self.map, &self.style);
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
        let h = w.tectonics.height_m.sample(d);
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
        let dt = ui.ctx().input(|i| i.unstable_dt);
        self.stats.push(dt);
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
        let (y, p) = (self.camera.yaw + 0.55, (self.camera.pitch + 0.35).min(1.4));
        self.style.sun_dir = glam::DVec3::new(p.cos() * y.cos(), p.cos() * y.sin(), p.sin());
        ui::draw(self, ui, &rs);
    }
}

fn random_seed() -> u64 {
    // Time-based on purpose: this only picks a new seed. Generation itself is deterministic.
    let t = web_time::SystemTime::now().duration_since(web_time::UNIX_EPOCH).unwrap_or_default();
    (t.as_nanos() as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) % 1_000_000
}
