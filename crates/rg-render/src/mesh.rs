use bytemuck::{Pod, Zeroable};
use glam::{DVec3, Vec3};
use rg_geo::CubeFace;

#[repr(C)]
#[derive(Clone, Copy, Debug, Pod, Zeroable)]
pub struct Vertex {
    /// Position in planet radii (surface ≈ 1, displaced by exaggerated height).
    pub pos: [f32; 3],
    /// Height above sea level in metres. The shader colours from this, so colour changes
    /// need no mesh rebuild.
    pub height_m: f32,
    pub normal: [f32; 3],
    pub _pad: f32,
}

/// A whole-planet preview mesh: `res × res` quads per cube face.
pub struct PlanetMesh {
    pub vertices: Vec<Vertex>,
    pub indices: Vec<u32>,
}

impl PlanetMesh {
    /// `height_m(dir)` returns height above sea level in metres for a unit direction.
    /// `exaggeration` scales relief for display only.
    pub fn build(
        res: u32,
        radius_m: f64,
        exaggeration: f64,
        height_m: impl Fn(DVec3) -> f64,
    ) -> Self {
        let n = res + 1;
        let mut vertices = Vec::with_capacity((6 * n * n) as usize);
        let mut indices = Vec::with_capacity((6 * res * res * 6) as usize);
        for face in CubeFace::ALL {
            let base = vertices.len() as u32;
            for j in 0..n {
                for i in 0..n {
                    let u = -1.0 + 2.0 * i as f64 / res as f64;
                    let v = -1.0 + 2.0 * j as f64 / res as f64;
                    let dir = face.uv_to_dir(u, v);
                    let h = height_m(dir);
                    // Oceans render at sea level; the seabed shows through colour only.
                    let shown = h.max(0.0) * exaggeration;
                    let p = dir * (1.0 + shown / radius_m);
                    vertices.push(Vertex {
                        pos: p.as_vec3().into(),
                        height_m: h as f32,
                        normal: [0.0; 3],
                        _pad: 0.0,
                    });
                }
            }
            for j in 0..res {
                for i in 0..res {
                    let a = base + j * n + i;
                    let b = a + 1;
                    let c = a + n;
                    let d = c + 1;
                    // Counter-clockwise seen from outside (u × v = outward normal).
                    indices.extend_from_slice(&[a, b, d, a, d, c]);
                }
            }
        }
        let mut mesh = Self { vertices, indices };
        mesh.compute_normals();
        mesh
    }

    fn compute_normals(&mut self) {
        let mut acc = vec![Vec3::ZERO; self.vertices.len()];
        for t in self.indices.as_chunks::<3>().0 {
            let [a, b, c] = [t[0] as usize, t[1] as usize, t[2] as usize];
            let pa = Vec3::from(self.vertices[a].pos);
            let pb = Vec3::from(self.vertices[b].pos);
            let pc = Vec3::from(self.vertices[c].pos);
            let n = (pb - pa).cross(pc - pa); // area-weighted
            acc[a] += n;
            acc[b] += n;
            acc[c] += n;
        }
        for (v, n) in self.vertices.iter_mut().zip(acc) {
            let fallback = Vec3::from(v.pos).normalize_or_zero();
            v.normal = n.try_normalize().unwrap_or(fallback).into();
        }
    }
}
