use crate::commands::Command;
use crate::ui;
use rg_render::{GlobeRenderer, GlobeStyle, OrbitCamera, PlanetMesh};

/// Generation parameters the UI edits. Placeholder until rg-tectonics lands.
#[derive(Clone, Debug)]
pub struct GenParams {
    pub seed: u64,
    pub radius_km: f64,
    pub sea_level: f64,
    pub exaggeration: f64,
    pub mesh_res: u32,
}

impl Default for GenParams {
    fn default() -> Self {
        Self { seed: 1, radius_km: 6371.0, sea_level: 0.05, exaggeration: 20.0, mesh_res: 160 }
    }
}

/// Rolling frame-time stats for the HUD.
#[derive(Default)]
pub struct FrameStats {
    samples: std::collections::VecDeque<f32>,
    pub last_gen_ms: f64,
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
    pub params: GenParams,
    pub camera: OrbitCamera,
    pub style: GlobeStyle,
    pub stats: FrameStats,
    pub adapter_info: String,
    renderer: Option<GlobeRenderer>,
    texture_id: Option<egui::TextureId>,
    pending: Vec<Command>,
}

impl RealmApp {
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        ui::theme::apply(&cc.egui_ctx);
        let mut app = Self {
            params: GenParams::default(),
            camera: OrbitCamera::default(),
            style: GlobeStyle::default(),
            stats: FrameStats::default(),
            adapter_info: String::new(),
            renderer: None,
            texture_id: None,
            pending: vec![Command::Regenerate],
        };
        if let Some(rs) = cc.wgpu_render_state.as_ref() {
            let info = rs.adapter.get_info();
            app.adapter_info = format!("{} · {:?}", info.name, info.backend);
            app.renderer = Some(GlobeRenderer::new(&rs.device));
        }
        app
    }

    pub fn push(&mut self, c: Command) {
        self.pending.push(c);
    }

    fn apply(&mut self, c: Command, rs: &egui_wgpu::RenderState) {
        match c {
            Command::Regenerate => self.regenerate(rs),
            Command::RandomizeSeed => {
                self.params.seed = random_seed();
                self.regenerate(rs);
            }
            Command::ResetCamera => self.camera = OrbitCamera::default(),
        }
    }

    fn regenerate(&mut self, rs: &egui_wgpu::RenderState) {
        let t0 = web_time::Instant::now();
        let p = self.params.clone();
        let radius_m = p.radius_km * 1000.0;
        let fbm = rg_noise::Fbm { octaves: 9, frequency: 1.4, lacunarity: 2.03, gain: 0.5 };
        let salt = rg_noise::salt(p.seed, 0);
        let ridge_salt = rg_noise::salt(p.seed, 1);
        // Placeholder relief: continents from fBm, ridged mountains on land. Replaced by
        // rg-tectonics in M2.
        let mesh = PlanetMesh::build(p.mesh_res, radius_m, p.exaggeration, |d| {
            let c = fbm.sample(d, salt) * 2.2 - p.sea_level;
            let mountains = fbm.ridged(d * 1.7, ridge_salt);
            if c > 0.0 {
                c.min(0.6) * 2500.0 + (c * 4.0).min(1.0) * mountains.powi(3) * 6500.0
            } else {
                (c * 9000.0).max(-6500.0)
            }
        });
        self.stats.triangles = mesh.indices.len() / 3;
        if let Some(r) = self.renderer.as_mut() {
            r.upload_mesh(&rs.device, &mesh);
        }
        self.stats.last_gen_ms = t0.elapsed().as_secs_f64() * 1000.0;
    }

    /// Render the globe into its offscreen texture at `size_px` and return the egui id.
    pub fn globe_texture(&mut self, rs: &egui_wgpu::RenderState, size_px: [u32; 2]) -> Option<egui::TextureId> {
        let r = self.renderer.as_mut()?;
        let recreated = r.resize(&rs.device, size_px);
        r.render(&rs.device, &rs.queue, &self.camera, &self.style);
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
}

impl eframe::App for RealmApp {
    fn ui(&mut self, ui: &mut egui::Ui, frame: &mut eframe::Frame) {
        let dt = ui.ctx().input(|i| i.unstable_dt);
        self.stats.push(dt);
        let Some(rs) = frame.wgpu_render_state().cloned() else {
            ui.label("RealmGenesis needs a GPU (wgpu) renderer.");
            return;
        };
        for c in std::mem::take(&mut self.pending) {
            self.apply(c, &rs);
        }
        ui::draw(self, ui, &rs);
    }
}

fn random_seed() -> u64 {
    // Time-based on purpose: this only picks a new seed. Generation itself is deterministic.
    let t = web_time::SystemTime::now().duration_since(web_time::UNIX_EPOCH).unwrap_or_default();
    (t.as_nanos() as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15) % 1_000_000
}
