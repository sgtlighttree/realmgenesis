//! The planet renderer: a 3D globe or a 2D equirectangular map, drawn offscreen.
//!
//! Both views sample the same **layer texture** (an equirectangular sRGB albedo image the
//! app bakes for the active layer: satellite, Köppen, plates…) and **height texture**
//! (metres, used for 2D hillshade). The globe additionally uses a displaced mesh.

use crate::terrain::{TileCache, TileInstance, TileVertex};
use crate::EarthCamera;
use bytemuck::{Pod, Zeroable};
use glam::DVec3;

pub const COLOR_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Rgba8UnormSrgb;
pub const DEPTH_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Depth32Float;
const SAMPLES: u32 = 4;

/// An equirectangular RGBA8 (sRGB) image: row 0 = north pole, column 0 = lon -180°.
#[derive(Clone, Debug)]
pub struct LayerImage {
    pub width: u32,
    pub height: u32,
    pub rgba: Vec<u8>,
}

/// An equirectangular height image in metres above sea level.
#[derive(Clone, Debug)]
/// The displayed surface's gradient on the cube faces (`rg_terrain::gradient`): six faces
/// of `(res + 1)²` samples, face-major then row-major, each `[d h/d east, d h/d north]`.
pub struct GradientImage {
    pub res: u32,
    pub data: Vec<[f32; 2]>,
}

pub struct HeightImage {
    pub width: u32,
    pub height: u32,
    pub metres: Vec<f32>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ViewMode {
    Globe,
    Map,
}

/// The 2D map's viewport.
#[derive(Clone, Copy, Debug)]
pub struct MapView {
    pub center_lon: f64,
    pub center_lat: f64,
    /// Visible longitude span in radians (2π = whole world).
    pub width_rad: f64,
}

impl Default for MapView {
    fn default() -> Self {
        Self { center_lon: 0.0, center_lat: 0.0, width_rad: std::f64::consts::TAU * 1.05 }
    }
}

impl MapView {
    /// Pan by a screen delta in pixels.
    pub fn drag(&mut self, dx: f64, dy: f64, viewport_w: f64) {
        let rad_per_px = self.width_rad / viewport_w.max(1.0);
        self.center_lon = (self.center_lon - dx * rad_per_px + std::f64::consts::PI)
            .rem_euclid(std::f64::consts::TAU)
            - std::f64::consts::PI;
        self.center_lat = (self.center_lat + dy * rad_per_px).clamp(-1.5, 1.5);
    }

    pub fn zoom(&mut self, factor: f64) {
        self.width_rad = (self.width_rad * factor).clamp(0.002, std::f64::consts::TAU * 1.2);
    }

    /// (lon, lat) radians under a point given in viewport fractions (0..1, y down).
    pub fn lonlat_at(&self, fx: f64, fy: f64, aspect: f64) -> (f64, f64) {
        let lon = self.center_lon + (fx - 0.5) * self.width_rad;
        let lat = self.center_lat - (fy - 0.5) * self.width_rad / aspect;
        (lon, lat)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct RenderStyle {
    pub sun_dir: DVec3,
    /// 0 = flat colour, 1 = full relief shading.
    pub relief_shading: f32,
    pub ocean_specular: bool,
    pub radius_m: f64,
    /// Vertical exaggeration for 2D hillshade.
    pub map_exaggeration: f32,
    /// River tint (sRGB 0–255) applied over the river coverage in the layer's alpha.
    pub river_rgb: [u8; 3],
    /// Latitude/longitude grid (the degree labels are drawn by the UI).
    pub graticule: bool,
}

impl Default for RenderStyle {
    fn default() -> Self {
        Self {
            sun_dir: DVec3::new(1.0, 0.6, 0.35),
            relief_shading: 1.0,
            ocean_specular: true,
            radius_m: 6_371_000.0,
            map_exaggeration: 12.0,
            river_rgb: [46, 96, 132],
            graticule: true,
        }
    }
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct Globals {
    view_proj: [[f32; 4]; 4],
    sun_dir: [f32; 4],
    eye: [f32; 4],
    params: [f32; 4],
    map: [f32; 4],
    /// Globe LOD: x = pixels per radian, y = cell size (px) at which a tile shows its
    /// parent's surface (morph 0), z = morph range (px). See `vs_tile`.
    /// w = lat/lon grid spacing in degrees (0 = off).
    lod: [f32; 4],
    /// xyz: river tint (linear), w: viewport width (px).
    river: [f32; 4],
}

struct Targets {
    size: [u32; 2],
    msaa: wgpu::TextureView,
    resolve: wgpu::Texture,
    resolve_view: wgpu::TextureView,
    depth: wgpu::TextureView,
}

pub struct PlanetRenderer {
    tile_pipeline: wgpu::RenderPipeline,
    map_pipeline: wgpu::RenderPipeline,
    bgl: wgpu::BindGroupLayout,
    globals: wgpu::Buffer,
    sampler: wgpu::Sampler,
    albedo: wgpu::Texture,
    heights: wgpu::Texture,
    /// Surface gradient per cube face (`Rg16Float`, 6 layers) for per-pixel lighting.
    gradient: wgpu::Texture,
    clamp_sampler: wgpu::Sampler,
    bind_group: wgpu::BindGroup,
    /// The streamed globe terrain.
    pub tiles: TileCache,
    /// Split threshold (projected cell size, px) used by selection; geomorphing must use
    /// the same scale so a tile equals its parent exactly where the parent would be drawn.
    pub max_cell_px: f64,
    targets: Option<Targets>,
}

impl PlanetRenderer {
    /// `tile_indices` / `tile_cells`: the shared tile index buffer and cells per tile edge
    /// (from rg-terrain). `tile_capacity`: max resident tiles.
    pub fn new(device: &wgpu::Device, queue: &wgpu::Queue, tile_indices: &[u32], tile_cells: u32, tile_capacity: usize) -> Self {
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("planet"),
            source: wgpu::ShaderSource::Wgsl(include_str!("planet.wgsl").into()),
        });
        let tex_entry = |binding, sample_type| wgpu::BindGroupLayoutEntry {
            binding,
            visibility: wgpu::ShaderStages::FRAGMENT,
            ty: wgpu::BindingType::Texture {
                sample_type,
                view_dimension: wgpu::TextureViewDimension::D2,
                multisampled: false,
            },
            count: None,
        };
        let bgl = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("planet"),
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::VERTEX_FRAGMENT,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                tex_entry(1, wgpu::TextureSampleType::Float { filterable: true }),
                wgpu::BindGroupLayoutEntry {
                    binding: 2,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Sampler(wgpu::SamplerBindingType::Filtering),
                    count: None,
                },
                tex_entry(3, wgpu::TextureSampleType::Float { filterable: false }),
                wgpu::BindGroupLayoutEntry {
                    binding: 4,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Texture {
                        sample_type: wgpu::TextureSampleType::Float { filterable: true },
                        view_dimension: wgpu::TextureViewDimension::D2Array,
                        multisampled: false,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 5,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Sampler(wgpu::SamplerBindingType::Filtering),
                    count: None,
                },
            ],
        });
        let globals = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("planet globals"),
            size: std::mem::size_of::<Globals>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let sampler = device.create_sampler(&wgpu::SamplerDescriptor {
            label: Some("layer"),
            address_mode_u: wgpu::AddressMode::Repeat,
            address_mode_v: wgpu::AddressMode::ClampToEdge,
            mag_filter: wgpu::FilterMode::Linear,
            min_filter: wgpu::FilterMode::Linear,
            ..Default::default()
        });
        // Cube-face textures are corner-aligned and end exactly at the face edge: clamp.
        let clamp_sampler = device.create_sampler(&wgpu::SamplerDescriptor {
            label: Some("cube faces"),
            mag_filter: wgpu::FilterMode::Linear,
            min_filter: wgpu::FilterMode::Linear,
            ..Default::default()
        });
        let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("planet"),
            bind_group_layouts: &[Some(&bgl)],
            immediate_size: 0,
        });
        let vattrs = wgpu::vertex_attr_array![0 => Float32x3, 1 => Float32, 2 => Float32x3, 3 => Float32, 4 => Float32];
        let iattrs = wgpu::vertex_attr_array![5 => Float32x3, 6 => Float32, 7 => Float32, 8 => Float32];
        let tile_buffers = [
            Some(wgpu::VertexBufferLayout {
                array_stride: std::mem::size_of::<TileVertex>() as u64,
                step_mode: wgpu::VertexStepMode::Vertex,
                attributes: &vattrs,
            }),
            Some(wgpu::VertexBufferLayout {
                array_stride: std::mem::size_of::<TileInstance>() as u64,
                step_mode: wgpu::VertexStepMode::Instance,
                attributes: &iattrs,
            }),
        ];
        let make_pipeline = |label, vs, fs, buffers: &[Option<wgpu::VertexBufferLayout>], depth: bool, cull| {
            device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
                label: Some(label),
                layout: Some(&layout),
                vertex: wgpu::VertexState {
                    module: &shader,
                    entry_point: Some(vs),
                    compilation_options: Default::default(),
                    buffers,
                },
                primitive: wgpu::PrimitiveState { cull_mode: cull, ..Default::default() },
                depth_stencil: Some(wgpu::DepthStencilState {
                    format: DEPTH_FORMAT,
                    depth_write_enabled: Some(depth),
                    // Reversed Z: nearer = larger depth.
                    depth_compare: Some(if depth { wgpu::CompareFunction::Greater } else { wgpu::CompareFunction::Always }),
                    stencil: Default::default(),
                    bias: Default::default(),
                }),
                multisample: wgpu::MultisampleState { count: SAMPLES, ..Default::default() },
                fragment: Some(wgpu::FragmentState {
                    module: &shader,
                    entry_point: Some(fs),
                    compilation_options: Default::default(),
                    targets: &[Some(wgpu::ColorTargetState {
                        format: COLOR_FORMAT,
                        blend: None,
                        write_mask: wgpu::ColorWrites::ALL,
                    })],
                }),
                multiview_mask: None,
                cache: None,
            })
        };
        let tile_pipeline = make_pipeline("tiles", "vs_tile", "fs_tile", &tile_buffers, true, Some(wgpu::Face::Back));
        let map_pipeline = make_pipeline("map", "vs_map", "fs_map", &[], false, None);

        let (albedo, heights) = Self::make_textures(device, 1, 1, 1, 1);
        let gradient = Self::make_gradient_texture(device, 1);
        let bind_group = Self::make_bind_group(device, &bgl, &globals, &albedo, &sampler, &heights, &gradient, &clamp_sampler);
        let tiles = TileCache::new(device, tile_indices, tile_cells, tile_capacity);
        let mut r = Self {
            tile_pipeline,
            map_pipeline,
            bgl,
            globals,
            sampler,
            albedo,
            heights,
            gradient,
            clamp_sampler,
            bind_group,
            tiles,
            max_cell_px: 6.0,
            targets: None,
        };
        r.set_layer(device, queue, &LayerImage { width: 1, height: 1, rgba: vec![40, 80, 120, 255] });
        r.set_heights(device, queue, &HeightImage { width: 1, height: 1, metres: vec![0.0] });
        r.set_gradient(device, queue, &GradientImage { res: 1, data: vec![[0.0; 2]; 6 * 4] });
        r
    }

    fn make_gradient_texture(device: &wgpu::Device, res: u32) -> wgpu::Texture {
        device.create_texture(&wgpu::TextureDescriptor {
            label: Some("surface gradient"),
            size: wgpu::Extent3d { width: res + 1, height: res + 1, depth_or_array_layers: 6 },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: wgpu::TextureFormat::Rg16Float,
            usage: wgpu::TextureUsages::TEXTURE_BINDING | wgpu::TextureUsages::COPY_DST,
            view_formats: &[],
        })
    }

    fn make_textures(device: &wgpu::Device, aw: u32, ah: u32, hw: u32, hh: u32) -> (wgpu::Texture, wgpu::Texture) {
        let make = |label, w, h, format| {
            device.create_texture(&wgpu::TextureDescriptor {
                label: Some(label),
                size: wgpu::Extent3d { width: w, height: h, depth_or_array_layers: 1 },
                mip_level_count: 1,
                sample_count: 1,
                dimension: wgpu::TextureDimension::D2,
                format,
                usage: wgpu::TextureUsages::TEXTURE_BINDING | wgpu::TextureUsages::COPY_DST,
                view_formats: &[],
            })
        };
        (
            make("layer albedo", aw, ah, wgpu::TextureFormat::Rgba8UnormSrgb),
            make("layer heights", hw, hh, wgpu::TextureFormat::R32Float),
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn make_bind_group(
        device: &wgpu::Device,
        bgl: &wgpu::BindGroupLayout,
        globals: &wgpu::Buffer,
        albedo: &wgpu::Texture,
        sampler: &wgpu::Sampler,
        heights: &wgpu::Texture,
        gradient: &wgpu::Texture,
        clamp_sampler: &wgpu::Sampler,
    ) -> wgpu::BindGroup {
        let gradient_view = gradient.create_view(&wgpu::TextureViewDescriptor {
            dimension: Some(wgpu::TextureViewDimension::D2Array),
            ..Default::default()
        });
        device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("planet"),
            layout: bgl,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: globals.as_entire_binding() },
                wgpu::BindGroupEntry {
                    binding: 1,
                    resource: wgpu::BindingResource::TextureView(&albedo.create_view(&Default::default())),
                },
                wgpu::BindGroupEntry { binding: 2, resource: wgpu::BindingResource::Sampler(sampler) },
                wgpu::BindGroupEntry {
                    binding: 3,
                    resource: wgpu::BindingResource::TextureView(&heights.create_view(&Default::default())),
                },
                wgpu::BindGroupEntry { binding: 4, resource: wgpu::BindingResource::TextureView(&gradient_view) },
                wgpu::BindGroupEntry { binding: 5, resource: wgpu::BindingResource::Sampler(clamp_sampler) },
            ],
        })
    }

    fn rebind(&mut self, device: &wgpu::Device) {
        self.bind_group = Self::make_bind_group(
            device,
            &self.bgl,
            &self.globals,
            &self.albedo,
            &self.sampler,
            &self.heights,
            &self.gradient,
            &self.clamp_sampler,
        );
    }

    /// Replace the layer (albedo) texture.
    pub fn set_layer(&mut self, device: &wgpu::Device, queue: &wgpu::Queue, img: &LayerImage) {
        assert_eq!(img.rgba.len(), (img.width * img.height * 4) as usize);
        if self.albedo.width() != img.width || self.albedo.height() != img.height {
            self.albedo = Self::make_textures(device, img.width, img.height, 1, 1).0;
            self.rebind(device);
        }
        write_texture(queue, &self.albedo, img.width, img.height, 4, &img.rgba);
    }

    /// Replace the height texture used for 2D hillshade.
    pub fn set_heights(&mut self, device: &wgpu::Device, queue: &wgpu::Queue, img: &HeightImage) {
        assert_eq!(img.metres.len(), (img.width * img.height) as usize);
        if self.heights.width() != img.width || self.heights.height() != img.height {
            self.heights = Self::make_textures(device, 1, 1, img.width, img.height).1;
            self.rebind(device);
        }
        write_texture(queue, &self.heights, img.width, img.height, 4, bytemuck::cast_slice(&img.metres));
    }

    /// Replace the surface-gradient texture (per-pixel lighting on the globe). Stored as
    /// half floats: ~3 significant digits, plenty for slopes.
    pub fn set_gradient(&mut self, device: &wgpu::Device, queue: &wgpu::Queue, img: &GradientImage) {
        let n = img.res + 1;
        assert_eq!(img.data.len(), (6 * n * n) as usize);
        if self.gradient.width() != n {
            self.gradient = Self::make_gradient_texture(device, img.res);
            self.rebind(device);
        }
        let halves: Vec<u16> =
            img.data.iter().flat_map(|g| g.map(|x| half::f16::from_f32(x).to_bits())).collect();
        queue.write_texture(
            wgpu::TexelCopyTextureInfo {
                texture: &self.gradient,
                mip_level: 0,
                origin: wgpu::Origin3d::ZERO,
                aspect: wgpu::TextureAspect::All,
            },
            bytemuck::cast_slice(&halves),
            wgpu::TexelCopyBufferLayout { offset: 0, bytes_per_row: Some(n * 4), rows_per_image: Some(n) },
            wgpu::Extent3d { width: n, height: n, depth_or_array_layers: 6 },
        );
    }

    /// Ensure the offscreen targets match `size`. Returns true if they were recreated, in
    /// which case the caller must re-register `output_view()` with its UI.
    pub fn resize(&mut self, device: &wgpu::Device, size: [u32; 2]) -> bool {
        let size = [size[0].max(1), size[1].max(1)];
        if self.targets.as_ref().is_some_and(|t| t.size == size) {
            return false;
        }
        let extent = wgpu::Extent3d { width: size[0], height: size[1], depth_or_array_layers: 1 };
        let make = |label, format, samples, usage| {
            device.create_texture(&wgpu::TextureDescriptor {
                label: Some(label),
                size: extent,
                mip_level_count: 1,
                sample_count: samples,
                dimension: wgpu::TextureDimension::D2,
                format,
                usage,
                view_formats: &[],
            })
        };
        let msaa = make("planet msaa", COLOR_FORMAT, SAMPLES, wgpu::TextureUsages::RENDER_ATTACHMENT);
        let resolve = make(
            "planet output",
            COLOR_FORMAT,
            1,
            wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::TEXTURE_BINDING | wgpu::TextureUsages::COPY_SRC,
        );
        let depth = make("planet depth", DEPTH_FORMAT, SAMPLES, wgpu::TextureUsages::RENDER_ATTACHMENT);
        self.targets = Some(Targets {
            size,
            msaa: msaa.create_view(&Default::default()),
            resolve_view: resolve.create_view(&Default::default()),
            resolve,
            depth: depth.create_view(&Default::default()),
        });
        true
    }

    /// Read the last rendered frame back from the GPU: (width, height, RGBA8 sRGB rows).
    /// Blocks until the GPU is done, so it's for screenshots, not per frame. The app's own
    /// screenshot works without macOS screen-recording permission (s03).
    pub fn read_output(&self, device: &wgpu::Device, queue: &wgpu::Queue) -> Option<(u32, u32, Vec<u8>)> {
        let t = self.targets.as_ref()?;
        let [w, h] = t.size;
        // Rows in a texture-to-buffer copy must be padded to 256 bytes.
        let row = w * 4;
        let padded = row.div_ceil(wgpu::COPY_BYTES_PER_ROW_ALIGNMENT) * wgpu::COPY_BYTES_PER_ROW_ALIGNMENT;
        let buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("screenshot"),
            size: (padded * h) as u64,
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });
        let mut enc = device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("screenshot") });
        enc.copy_texture_to_buffer(
            wgpu::TexelCopyTextureInfo {
                texture: &t.resolve,
                mip_level: 0,
                origin: wgpu::Origin3d::ZERO,
                aspect: wgpu::TextureAspect::All,
            },
            wgpu::TexelCopyBufferInfo {
                buffer: &buf,
                layout: wgpu::TexelCopyBufferLayout { offset: 0, bytes_per_row: Some(padded), rows_per_image: Some(h) },
            },
            wgpu::Extent3d { width: w, height: h, depth_or_array_layers: 1 },
        );
        queue.submit([enc.finish()]);
        buf.map_async(wgpu::MapMode::Read, .., |_| {});
        device.poll(wgpu::PollType::Wait { submission_index: None, timeout: None }).ok()?;
        let data = buf.get_mapped_range(..).ok()?;
        let mut out = Vec::with_capacity((row * h) as usize);
        for y in 0..h as usize {
            let start = y * padded as usize;
            out.extend_from_slice(&data[start..start + row as usize]);
        }
        drop(data);
        buf.unmap();
        Some((w, h, out))
    }

    pub fn output_view(&self) -> Option<&wgpu::TextureView> {
        self.targets.as_ref().map(|t| &t.resolve_view)
    }

    /// Draw one frame into the offscreen target. Call `resize` first, and for the globe,
    /// `tiles.select(..)` first (it decides what gets drawn).
    #[allow(clippy::too_many_arguments)]
    pub fn render(
        &mut self,
        device: &wgpu::Device,
        queue: &wgpu::Queue,
        mode: ViewMode,
        cam: &EarthCamera,
        clearance_m: f64,
        exaggeration: f64,
        map: &MapView,
        style: &RenderStyle,
    ) {
        let Some(t) = &self.targets else { return };
        let aspect = t.size[0] as f64 / t.size[1] as f64;
        let eye = cam.eye();
        let s = style.sun_dir.normalize();
        let w = match mode {
            ViewMode::Globe => exaggeration as f32,
            ViewMode::Map => style.map_exaggeration,
        };
        let g = Globals {
            view_proj: cam.view_proj(aspect, clearance_m).to_cols_array_2d(),
            sun_dir: [s.x as f32, s.y as f32, s.z as f32, 0.0],
            eye: [eye.x as f32, eye.y as f32, eye.z as f32, 1.0],
            params: [style.relief_shading, if style.ocean_specular { 1.0 } else { 0.0 }, style.radius_m as f32, w],
            map: [map.center_lon as f32, map.center_lat as f32, map.width_rad as f32, aspect as f32],
            lod: [
                (t.size[1] as f64 / (2.0 * (0.5 * cam.fov_y).tan())) as f32,
                (self.max_cell_px * 0.5) as f32,
                (self.max_cell_px * 0.35) as f32,
                if style.graticule { crate::graticule::step_deg(mode, cam, map, t.size) as f32 } else { 0.0 },
            ],
            river: {
                let lin = |c: u8| (c as f32 / 255.0).powf(2.2);
                [lin(style.river_rgb[0]), lin(style.river_rgb[1]), lin(style.river_rgb[2]), t.size[0] as f32]
            },
        };
        queue.write_buffer(&self.globals, 0, bytemuck::bytes_of(&g));

        // Background: space from orbit, blending to a daytime sky as the camera descends
        // into the atmosphere (globe view only).
        let space = [0.006, 0.008, 0.014];
        let sky = [0.30, 0.48, 0.80];
        let t = match mode {
            ViewMode::Globe => {
                let alt = eye.length() - style.radius_m;
                let x = ((alt - 0.004 * style.radius_m) / (0.03 * style.radius_m - 0.004 * style.radius_m)).clamp(0.0, 1.0);
                1.0 - x * x * (3.0 - 2.0 * x)
            }
            ViewMode::Map => 0.0,
        };
        let mix = |a: f64, b: f64| a + (b - a) * t;
        let clear = wgpu::Color { r: mix(space[0], sky[0]), g: mix(space[1], sky[1]), b: mix(space[2], sky[2]), a: 1.0 };
        let mut enc = device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: Some("planet") });
        {
            let t = self.targets.as_ref().unwrap();
            let mut pass = enc.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("planet"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &t.msaa,
                    depth_slice: None,
                    resolve_target: Some(&t.resolve_view),
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(clear),
                        store: wgpu::StoreOp::Discard,
                    },
                })],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &t.depth,
                    depth_ops: Some(wgpu::Operations { load: wgpu::LoadOp::Clear(0.0), store: wgpu::StoreOp::Discard }),
                    stencil_ops: None,
                }),
                timestamp_writes: None,
                occlusion_query_set: None,
                multiview_mask: None,
            });
            pass.set_bind_group(0, &self.bind_group, &[]);
            match mode {
                ViewMode::Globe => {
                    pass.set_pipeline(&self.tile_pipeline);
                    self.tiles.draw(device, queue, &mut pass);
                }
                ViewMode::Map => {
                    pass.set_pipeline(&self.map_pipeline);
                    pass.draw(0..3, 0..1);
                }
            }
        }
        queue.submit([enc.finish()]);
    }
}

fn write_texture(queue: &wgpu::Queue, tex: &wgpu::Texture, w: u32, h: u32, bpp: u32, data: &[u8]) {
    queue.write_texture(
        wgpu::TexelCopyTextureInfo { texture: tex, mip_level: 0, origin: wgpu::Origin3d::ZERO, aspect: wgpu::TextureAspect::All },
        data,
        wgpu::TexelCopyBufferLayout { offset: 0, bytes_per_row: Some(w * bpp), rows_per_image: Some(h) },
        wgpu::Extent3d { width: w, height: h, depth_or_array_layers: 1 },
    );
}
