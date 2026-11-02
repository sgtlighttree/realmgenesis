//! Streaming terrain: the GPU tile cache and per-frame LOD selection.
//!
//! The renderer never generates tiles. Each frame, [`TileCache::select`] decides what to
//! draw from what's resident and returns the tiles it wants; the app generates those
//! (worker threads natively, a time budget on the web) and hands them back via
//! [`TileCache::insert`]. See `docs/design/tile-streaming.md`.
//!
//! Rules (from the design review):
//! - Draw a tile or all of its visible children, never a mix (no z-fighting overlaps).
//! - While children are missing, keep drawing the parent.
//! - L0–L2 are pinned. Selection touches every ancestor it walks through, so the LRU never
//!   evicts the parent of something on screen.

use bytemuck::{Pod, Zeroable};
use glam::{DMat4, DVec3, DVec4};
use rg_geo::{CubeFace, TileId};
use std::collections::HashMap;
use wgpu::util::DeviceExt;

/// One terrain vertex (36 bytes). Positions are relative to the tile centre, on the
/// sea-level sphere; the shader lifts them by height × exaggeration.
#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
pub struct TileVertex {
    pub pos: [f32; 3],
    pub height_m: f32,
    pub normal: [f32; 3],
    pub parent_height_m: f32,
    pub skirt: f32,
}

/// Per-draw instance data (32 bytes): where the tile sits relative to the eye, its
/// morph-in progress and its skirt depth.
#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
pub struct TileInstance {
    pub offset: [f32; 3],
    pub morph: f32,
    pub skirt_m: f32,
    /// Ground size of one of this tile's cells, metres (for distance-based geomorphing).
    pub cell_m: f32,
    pub _pad: [f32; 2],
}

/// A tile ready for upload.
pub struct TileUpload {
    pub id: TileId,
    pub center: DVec3,
    pub vertices: Vec<TileVertex>,
    pub patch_radius_m: f64,
    pub max_height_m: f64,
    pub child_error_m: f64,
    pub parent_error_m: f64,
}

struct GpuTile {
    vbuf: wgpu::Buffer,
    center: DVec3,
    patch_radius_m: f64,
    max_height_m: f64,
    child_error_m: f64,
    parent_error_m: f64,
    born_s: f64,
    last_used: u64,
}

/// What selection needs to know about the view.
pub struct SelectView {
    pub eye: DVec3,
    /// Rotation-only view-projection (eye at the origin), f64.
    pub view_proj: DMat4,
    pub radius_m: f64,
    pub exaggeration: f64,
    pub viewport_h_px: f64,
    pub fov_y: f64,
    /// Allowed screen-space error, physical pixels.
    pub tau_px: f64,
    /// Largest allowed projected cell size on land, physical pixels (lighting detail).
    pub max_cell_px: f64,
    pub max_level: u8,
    pub now_s: f64,
}

/// Per-frame result of selection.
#[derive(Default)]
pub struct Selection {
    /// Tiles to generate, most important first.
    pub wanted: Vec<TileId>,
    /// True while any drawn tile is still morphing in (keep repainting).
    pub animating: bool,
    pub drawn: usize,
}

pub struct TileCache {
    tiles: HashMap<TileId, GpuTile>,
    pub index_buf: wgpu::Buffer,
    pub index_count: u32,
    instance_buf: wgpu::Buffer,
    instance_capacity: usize,
    /// (tile, instance) pairs selected for this frame.
    draw: Vec<(TileId, TileInstance)>,
    frame: u64,
    pub capacity: usize,
    /// Cell size of a level-0 tile at radius 1, for the curvature term.
    cells_per_tile: f64,
}

const MORPH_S: f64 = 0.25;

impl TileCache {
    pub fn new(device: &wgpu::Device, indices: &[u32], cells_per_tile: u32, capacity: usize) -> Self {
        let index_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("tile indices"),
            contents: bytemuck::cast_slice(indices),
            usage: wgpu::BufferUsages::INDEX,
        });
        let instance_capacity = 2048;
        let instance_buf = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("tile instances"),
            size: (instance_capacity * std::mem::size_of::<TileInstance>()) as u64,
            usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        Self {
            tiles: HashMap::new(),
            index_buf,
            index_count: indices.len() as u32,
            instance_buf,
            instance_capacity,
            draw: Vec::new(),
            frame: 0,
            capacity,
            cells_per_tile: cells_per_tile as f64,
        }
    }

    pub fn clear(&mut self) {
        self.tiles.clear();
        self.draw.clear();
    }

    pub fn len(&self) -> usize {
        self.tiles.len()
    }

    pub fn is_empty(&self) -> bool {
        self.tiles.is_empty()
    }

    pub fn contains(&self, id: &TileId) -> bool {
        self.tiles.contains_key(id)
    }

    pub fn insert(&mut self, device: &wgpu::Device, t: TileUpload, now_s: f64) {
        let vbuf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("tile vertices"),
            contents: bytemuck::cast_slice(&t.vertices),
            usage: wgpu::BufferUsages::VERTEX,
        });
        // Roots appear instantly; everything else morphs in from its parent's surface.
        let born_s = if t.id.level == 0 { now_s - MORPH_S } else { now_s };
        self.tiles.insert(
            t.id,
            GpuTile {
                vbuf,
                center: t.center,
                patch_radius_m: t.patch_radius_m,
                max_height_m: t.max_height_m,
                child_error_m: t.child_error_m,
                parent_error_m: t.parent_error_m,
                born_s,
                last_used: self.frame,
            },
        );
    }

    /// Choose what to draw this frame and what to generate next.
    pub fn select(&mut self, v: &SelectView) -> Selection {
        self.frame += 1;
        self.draw.clear();
        let mut sel = Selection::default();
        let planes = frustum_planes(&v.view_proj);
        let mut wanted: Vec<(f64, TileId)> = Vec::new();
        for face in CubeFace::ALL {
            let root = TileId::root(face);
            if !self.tiles.contains_key(&root) {
                wanted.push((-1.0, root));
                continue;
            }
            self.visit(root, v, &planes, &mut wanted, &mut sel);
        }
        wanted.sort_by(|a, b| a.0.total_cmp(&b.0));
        sel.wanted = wanted.into_iter().map(|(_, id)| id).collect();
        sel.drawn = self.draw.len();
        self.evict();
        sel
    }

    /// Bounding sphere of a resident tile under the current exaggeration: (centre, radius).
    fn bounds(&self, t: &GpuTile, v: &SelectView) -> (DVec3, f64) {
        let lift = t.max_height_m * v.exaggeration;
        let up = t.center.normalize();
        (t.center + up * (lift * 0.5), t.patch_radius_m + lift * 0.5 + t.parent_error_m * v.exaggeration)
    }

    fn visible(&self, c: DVec3, r: f64, v: &SelectView, planes: &[DVec4; 5]) -> bool {
        let rel = c - v.eye;
        for p in planes {
            if p.truncate().dot(rel) + p.w < -r * p.truncate().length() {
                return false;
            }
        }
        !occluded_by_horizon(c, r, v.eye, v.radius_m)
    }

    /// Pixels of geometric error if `id` is drawn instead of its children.
    fn error_px(&self, id: TileId, t: &GpuTile, c: DVec3, r: f64, v: &SelectView) -> f64 {
        let cell = std::f64::consts::FRAC_PI_2 * v.radius_m / (self.cells_per_tile * (1u64 << id.level) as f64);
        let sag = cell * cell / (8.0 * v.radius_m);
        let err = t.child_error_m * v.exaggeration + sag;
        let dist = ((c - v.eye).length() - r).max(1.0);
        let px_per_rad = v.viewport_h_px / (2.0 * (0.5 * v.fov_y).tan());
        err / dist * px_per_rad
    }

    /// Projected size of one of this tile's cells, pixels.
    fn cell_px(&self, id: TileId, c: DVec3, r: f64, v: &SelectView) -> f64 {
        let cell = std::f64::consts::FRAC_PI_2 * v.radius_m / (self.cells_per_tile * (1u64 << id.level) as f64);
        let dist = ((c - v.eye).length() - r).max(1.0);
        cell / dist * v.viewport_h_px / (2.0 * (0.5 * v.fov_y).tan())
    }

    fn visit(&mut self, id: TileId, v: &SelectView, planes: &[DVec4; 5], wanted: &mut Vec<(f64, TileId)>, sel: &mut Selection) {
        let frame = self.frame;
        let Some(t) = self.tiles.get_mut(&id) else { return };
        t.last_used = frame;
        let t = self.tiles.get(&id).unwrap();
        let (c, r) = self.bounds(t, v);
        if !self.visible(c, r, v, planes) {
            return;
        }
        // Split when the geometry is off by more than τ px, or (on land) when a cell
        // spans more than `max_cell_px`: geometry can be accurate while lighting still
        // needs denser normals. All-ocean tiles render flat, so only geometry matters.
        let land = t.max_height_m > 0.0;
        let split = id.level < v.max_level
            && (self.error_px(id, t, c, r, v) > v.tau_px || (land && self.cell_px(id, c, r, v) > v.max_cell_px));
        if split {
            let kids = id.children();
            let mut all_ready = true;
            for k in kids {
                if !self.tiles.contains_key(&k) {
                    // Children inherit the parent's bounds for a priority estimate.
                    all_ready = false;
                    let dist = (c - v.eye).length();
                    wanted.push((k.level as f64 * 1.0e9 + dist, k));
                }
            }
            if all_ready {
                for k in kids {
                    self.visit(k, v, planes, wanted, sel);
                }
                return;
            }
        }
        let t = self.tiles.get(&id).unwrap();
        let age = v.now_s - t.born_s;
        let morph = (age / MORPH_S).clamp(0.0, 1.0);
        if morph < 1.0 {
            sel.animating = true;
        }
        let cell = std::f64::consts::FRAC_PI_2 * v.radius_m / (self.cells_per_tile * (1u64 << id.level) as f64);
        let skirt = (t.parent_error_m * v.exaggeration + cell * 0.25 + 5.0) as f32;
        let off = (t.center - v.eye).as_vec3();
        self.draw.push((id, TileInstance { offset: off.to_array(), morph: morph as f32, skirt_m: skirt, cell_m: cell as f32, _pad: [0.0; 2] }));
    }

    fn evict(&mut self) {
        if self.tiles.len() <= self.capacity {
            return;
        }
        let mut old: Vec<(u64, TileId)> = self
            .tiles
            .iter()
            .filter(|(id, t)| id.level > 2 && t.last_used < self.frame)
            .map(|(id, t)| (t.last_used, *id))
            .collect();
        old.sort();
        let excess = self.tiles.len() - self.capacity;
        for (_, id) in old.into_iter().take(excess) {
            self.tiles.remove(&id);
        }
    }

    /// Record draws for the tiles chosen by the last `select`.
    pub fn draw<'a>(&'a mut self, device: &wgpu::Device, queue: &wgpu::Queue, pass: &mut wgpu::RenderPass<'a>) {
        if self.draw.len() > self.instance_capacity {
            self.instance_capacity = self.draw.len().next_power_of_two();
            self.instance_buf = device.create_buffer(&wgpu::BufferDescriptor {
                label: Some("tile instances"),
                size: (self.instance_capacity * std::mem::size_of::<TileInstance>()) as u64,
                usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            });
        }
        let inst: Vec<TileInstance> = self.draw.iter().map(|d| d.1).collect();
        if inst.is_empty() {
            return;
        }
        queue.write_buffer(&self.instance_buf, 0, bytemuck::cast_slice(&inst));
        pass.set_index_buffer(self.index_buf.slice(..), wgpu::IndexFormat::Uint32);
        pass.set_vertex_buffer(1, self.instance_buf.slice(..));
        for (slot, (id, _)) in self.draw.iter().enumerate() {
            if let Some(t) = self.tiles.get(id) {
                pass.set_vertex_buffer(0, t.vbuf.slice(..));
                let s = slot as u32;
                pass.draw_indexed(0..self.index_count, 0, s..s + 1);
            }
        }
    }
}

/// Left, right, bottom, top, near planes (Gribb–Hartmann) for a reversed-Z, infinite
/// projection with the eye at the origin. Inside ⇔ dot(n, p) + d ≥ 0.
fn frustum_planes(m: &DMat4) -> [DVec4; 5] {
    let r0 = m.row(0);
    let r1 = m.row(1);
    let r2 = m.row(2);
    let r3 = m.row(3);
    [r3 + r0, r3 - r0, r3 + r1, r3 - r1, r3 - r2]
}

/// Is the sphere (c, r) entirely hidden behind the planet (radius R) from `eye`?
/// Occluded ⇔ fully behind the horizon plane AND fully inside the eye's tangent cone.
fn occluded_by_horizon(c: DVec3, r: f64, eye: DVec3, radius: f64) -> bool {
    let d = eye.length();
    if d <= radius {
        return false;
    }
    let e_hat = eye / d;
    // Horizon plane: points X with X·ê = R²/d.
    if c.dot(e_hat) + r >= radius * radius / d {
        return false;
    }
    // Tangent cone around −ê with half-angle asin(R/d).
    let to = c - eye;
    let dist = to.length();
    if dist <= r {
        return false;
    }
    let angle = (to / dist).dot(-e_hat).clamp(-1.0, 1.0).acos();
    let alpha = (radius / d).asin();
    angle + (r / dist).asin() < alpha
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn horizon_hides_far_side_but_not_tall_limb() {
        let r = 6_371_000.0;
        let eye = DVec3::new(0.0, 0.0, r + 10_000.0);
        // Antipode: hidden.
        assert!(occluded_by_horizon(DVec3::new(0.0, 0.0, -r), 1000.0, eye, r));
        // Directly below: visible.
        assert!(!occluded_by_horizon(DVec3::new(0.0, 0.0, r), 1000.0, eye, r));
        // A point just past the horizon but lifted high: visible.
        let horizon_angle = (r / (r + 10_000.0)).acos();
        let a = horizon_angle * 1.05;
        let p = DVec3::new(a.sin(), 0.0, a.cos()) * (r + 20_000.0);
        assert!(!occluded_by_horizon(p, 1000.0, eye, r));
    }
}
