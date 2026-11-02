//! Level 0 of the drained detail network (`docs/design/drained-detail-network.md`): the
//! macro river network, as segments from each land node to its receiver.
//!
//! The erosion stage stores, per macro sample, the downstream flow direction (a tangent
//! vector) and the drainage area. The receiver is recovered as the neighbouring sample whose
//! direction best matches the flow vector, which is robust on the tangent-warped grid (a
//! fixed step along the vector can land on the wrong sample where cells are 1:√2).

use glam::DVec3;
use rg_geo::{CubeFace, CubeField};
use std::sync::Arc;

/// A channel segment from `up` to `down` (unit directions), with its channel elevation at
/// both ends (m, `z_up > z_down` for land segments) and the drainage area at `up` (km²).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Segment {
    pub up: DVec3,
    pub down: DVec3,
    pub z_up: f64,
    pub z_down: f64,
    pub area_km2: f64,
}

/// A macro sample address. `i`, `j` in `0..=res`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Node {
    pub face: CubeFace,
    pub i: u32,
    pub j: u32,
}

/// The macro river network: committed heights plus the erosion stage's hydrology, which must
/// share one resolution (the world file guarantees it since schema 3).
#[derive(Clone)]
pub struct MacroNet {
    pub height: Arc<CubeField<f32>>,
    pub flow_dir: Arc<CubeField<[f32; 3]>>,
    pub area_km2: Arc<CubeField<f32>>,
}

impl MacroNet {
    /// `None` if the fields don't share a resolution.
    pub fn new(height: Arc<CubeField<f32>>, flow_dir: Arc<CubeField<[f32; 3]>>, area_km2: Arc<CubeField<f32>>) -> Option<Self> {
        (height.res == flow_dir.res && height.res == area_km2.res).then_some(Self { height, flow_dir, area_km2 })
    }

    pub fn res(&self) -> u32 {
        self.height.res
    }

    /// The sample `(i, j)` on `face`, with indices possibly off the face: resolved through
    /// the direction to the owning face's nearest sample.
    pub fn node(&self, face: CubeFace, i: i64, j: i64) -> Node {
        let r = self.res() as i64;
        if (0..=r).contains(&i) && (0..=r).contains(&j) {
            return Node { face, i: i as u32, j: j as u32 };
        }
        let s = 2.0 / r as f64;
        let d = face.uv_to_dir(-1.0 + i as f64 * s, -1.0 + j as f64 * s);
        self.nearest_node(d)
    }

    pub fn nearest_node(&self, d: DVec3) -> Node {
        let (face, u, v) = CubeFace::from_dir(d);
        let r = self.res() as f64;
        let k = |t: f64| ((t + 1.0) * 0.5 * r).round().clamp(0.0, r) as u32;
        Node { face, i: k(u), j: k(v) }
    }

    pub fn dir(&self, n: Node) -> DVec3 {
        self.height.sample_dir(n.face, n.i, n.j)
    }

    fn index(&self, n: Node) -> usize {
        self.height.index(n.face, n.i, n.j)
    }

    pub fn height_at(&self, n: Node) -> f64 {
        self.height.data[self.index(n)] as f64
    }

    /// The node `n` drains to, or `None` at sea and on flats.
    pub fn receiver(&self, n: Node) -> Option<Node> {
        let k = self.index(n);
        if self.height.data[k] <= 0.0 {
            return None;
        }
        let f = self.flow_dir.data[k];
        let f = DVec3::new(f[0] as f64, f[1] as f64, f[2] as f64);
        if f.length_squared() < 0.25 {
            return None;
        }
        let d0 = self.dir(n);
        let mut best = (0.7, None); // within ~45° of the stored direction
        for dj in -1..=1i64 {
            for di in -1..=1i64 {
                if di == 0 && dj == 0 {
                    continue;
                }
                let m = self.node(n.face, n.i as i64 + di, n.j as i64 + dj);
                if m == n {
                    continue;
                }
                let t = self.dir(m) - d0;
                let t = (t - d0 * t.dot(d0)).normalize_or_zero();
                let c = t.dot(f);
                if c > best.0 {
                    best = (c, Some(m));
                }
            }
        }
        best.1
    }

    /// Land segments (node → receiver) whose upstream node lies in the 4×4 window of samples
    /// around `dir`. Every segment passing within one cell of `dir` has its upstream node
    /// there (segments join neighbouring samples).
    pub fn segments_near(&self, dir: DVec3, out: &mut Vec<Segment>) {
        out.clear();
        let (face, u, v) = CubeFace::from_dir(dir);
        let r = self.res() as f64;
        let (x, y) = ((u + 1.0) * 0.5 * r, (v + 1.0) * 0.5 * r);
        let (i0, j0) = (x.floor() as i64, y.floor() as i64);
        let mut seen: Vec<Node> = Vec::with_capacity(16);
        for j in j0 - 1..=j0 + 2 {
            for i in i0 - 1..=i0 + 2 {
                let n = self.node(face, i, j);
                if seen.contains(&n) {
                    continue;
                }
                seen.push(n);
                let Some(m) = self.receiver(n) else { continue };
                let z_down = self.height_at(m).max(0.0);
                out.push(Segment {
                    up: self.dir(n),
                    down: self.dir(m),
                    z_up: self.height_at(n),
                    z_down,
                    area_km2: self.area_km2.data[self.index(n)] as f64,
                });
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A tilted plane draining towards -x on every face's interior: receivers must exist,
    /// lie one step away and be lower.
    #[test]
    fn receivers_follow_the_stored_flow() {
        let res = 16;
        let h = CubeField::<f32>::from_fn(res, |d| (1000.0 + 500.0 * d.x) as f32);
        let mut flow = CubeField::<[f32; 3]>::new(res);
        for (k, _, _, _, d) in h.iter_samples().collect::<Vec<_>>() {
            // Downhill = -grad(x) projected on the tangent plane.
            let g = DVec3::NEG_X - d * (-d.x);
            let g = g.normalize_or_zero();
            flow.data[k] = [g.x as f32, g.y as f32, g.z as f32];
        }
        let area = CubeField::<f32>::new(res);
        let net = MacroNet::new(Arc::new(h), Arc::new(flow), Arc::new(area)).unwrap();
        let mut found = 0;
        for face in CubeFace::ALL {
            for j in 1..res {
                for i in 1..res {
                    let n = Node { face, i, j };
                    if net.dir(n).x.abs() > 0.95 {
                        continue; // gradient vanishes at the ±x poles of the field
                    }
                    let m = net.receiver(n).expect("receiver");
                    assert!(net.height_at(m) < net.height_at(n), "{n:?} → {m:?} not lower");
                    assert!(net.dir(m).angle_between(net.dir(n)) < 2.0 * std::f64::consts::FRAC_PI_2 / res as f64);
                    found += 1;
                }
            }
        }
        assert!(found > 1000);
    }
}
