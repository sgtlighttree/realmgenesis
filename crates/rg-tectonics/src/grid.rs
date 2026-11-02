//! The sample graph all the tectonics code works on.
//!
//! A `CubeField` stores samples on the edges of cube faces twice (once per face). For
//! simulation we want every *place* on the sphere to exist once, so this module merges those
//! twins into one **node**. Everything in the crate (noise, plate growth, distance fields)
//! is computed per node, and only at the very end copied back out to the `CubeField`s.
//!
//! The graph connects each node to its 8 surrounding nodes (including diagonals), also across
//! face edges, with edge length = straight-line distance on the unit sphere.

use glam::DVec3;
use rg_geo::{CubeFace, CubeField};
use crate::util::par_map;


pub struct Grid {
    /// Cells per cube-face edge (same meaning as `CubeField::res`).
    pub res: u32,
    /// Number of unique nodes (`6 * res^2 + 2`).
    pub n: usize,
    /// Unit direction of each node.
    pub dir: Vec<DVec3>,
    /// Share of the sphere's area belonging to each node (sums to 1).
    pub weight: Vec<f32>,
    /// CSR adjacency: neighbours of node `a` are `adj[adj_start[a]..adj_start[a + 1]]`.
    pub adj_start: Vec<u32>,
    pub adj: Vec<u32>,
    /// Shortest and longest edge in `adj`, radians of arc.
    pub min_edge: f32,
    pub max_edge: f32,
    /// For every `CubeField` sample index, the node it belongs to.
    pub sample_node: Vec<u32>,
    /// For every node, one representative `CubeField` sample index.
    pub node_sample: Vec<u32>,
}

impl Grid {
    pub fn new(res: u32) -> Grid {
        let probe = CubeField::<u8>::new(res);
        let per_edge = (res + 1) as usize;
        let total = 6 * per_edge * per_edge;

        // Direction of every sample.
        let mut sample_dir = vec![DVec3::ZERO; total];
        for (idx, _, _, _, d) in probe.iter_samples() {
            sample_dir[idx] = d;
        }

        // Find twin samples: samples on a face edge that sit at the same place as a sample on
        // another face. Sort them by quantised position and group equal keys.
        let is_edge = |i: u32, j: u32| i == 0 || j == 0 || i == res || j == res;
        let mut edge_keys: Vec<([i64; 3], u32)> = Vec::new();
        for face in CubeFace::ALL {
            for j in 0..=res {
                for i in 0..=res {
                    if is_edge(i, j) {
                        let idx = probe.index(face, i, j);
                        let d = sample_dir[idx];
                        let q = |x: f64| (x * 1.0e7).round() as i64;
                        edge_keys.push(([q(d.x), q(d.y), q(d.z)], idx as u32));
                    }
                }
            }
        }
        edge_keys.sort_unstable();
        // rep[idx] = lowest sample index in the twin group (itself for non-edge samples).
        let mut rep: Vec<u32> = (0..total as u32).collect();
        let mut k = 0;
        while k < edge_keys.len() {
            let mut e = k + 1;
            while e < edge_keys.len() && edge_keys[e].0 == edge_keys[k].0 {
                e += 1;
            }
            let lowest = edge_keys[k..e].iter().map(|x| x.1).min().unwrap();
            for item in &edge_keys[k..e] {
                rep[item.1 as usize] = lowest;
            }
            k = e;
        }

        // Number the nodes in sample order.
        let mut sample_node = vec![u32::MAX; total];
        let mut node_sample: Vec<u32> = Vec::new();
        for idx in 0..total {
            if rep[idx] as usize == idx {
                sample_node[idx] = node_sample.len() as u32;
                node_sample.push(idx as u32);
            }
        }
        for idx in 0..total {
            sample_node[idx] = sample_node[rep[idx] as usize];
        }
        let n = node_sample.len();
        let dir: Vec<DVec3> = node_sample.iter().map(|&s| sample_dir[s as usize]).collect();
        drop(sample_dir);
        drop(rep);
        drop(edge_keys);

        // Area weights: solid angle element of the tangent-warped cube face, with trapezoid
        // end weights on the face edges so twin samples share their cell areas correctly.
        let mut weight64 = vec![0f64; n];
        let quarter_pi = std::f64::consts::FRAC_PI_4;
        for face in CubeFace::ALL {
            for j in 0..=res {
                for i in 0..=res {
                    let u = -1.0 + 2.0 * i as f64 / res as f64;
                    let v = -1.0 + 2.0 * j as f64 / res as f64;
                    let (x, y) = ((u * quarter_pi).tan(), (v * quarter_pi).tan());
                    let w = (1.0 + x * x) * (1.0 + y * y) / (1.0 + x * x + y * y).powf(1.5);
                    let tw = |k: u32| if k == 0 || k == res { 0.5 } else { 1.0 };
                    let idx = probe.index(face, i, j);
                    weight64[sample_node[idx] as usize] += w * tw(i) * tw(j);
                }
            }
        }
        let wsum: f64 = weight64.iter().sum();
        let weight: Vec<f32> = weight64.iter().map(|w| (w / wsum) as f32).collect();
        drop(weight64);

        // Adjacency. Each sample contributes its in-face neighbours to its node; nodes made of
        // several twin samples therefore collect neighbours from every face they touch.
        let mut deg = vec![0u32; n];
        let offsets: [(i32, i32); 8] =
            [(-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1)];
        let in_range = |v: i32| v >= 0 && v <= res as i32;
        for face in CubeFace::ALL {
            for j in 0..=res as i32 {
                for i in 0..=res as i32 {
                    let a = sample_node[probe.index(face, i as u32, j as u32)] as usize;
                    for (di, dj) in offsets {
                        if in_range(i + di) && in_range(j + dj) {
                            deg[a] += 1;
                        }
                    }
                }
            }
        }
        let mut start = vec![0u32; n + 1];
        for a in 0..n {
            start[a + 1] = start[a] + deg[a];
        }
        let mut fill = start.clone();
        let mut raw = vec![0u32; start[n] as usize];
        for face in CubeFace::ALL {
            for j in 0..=res as i32 {
                for i in 0..=res as i32 {
                    let a = sample_node[probe.index(face, i as u32, j as u32)] as usize;
                    for (di, dj) in offsets {
                        if in_range(i + di) && in_range(j + dj) {
                            let b = sample_node[probe.index(face, (i + di) as u32, (j + dj) as u32)];
                            raw[fill[a] as usize] = b;
                            fill[a] += 1;
                        }
                    }
                }
            }
        }
        // Sort and de-duplicate each node's list (only twin-edge nodes have duplicates), moving
        // the survivors to the front of the same buffer to avoid a second copy.
        drop(fill);
        drop(deg);
        let mut adj_start = vec![0u32; n + 1];
        let mut w = 0usize;
        for a in 0..n {
            let (lo, hi) = (start[a] as usize, start[a + 1] as usize);
            raw[lo..hi].sort_unstable();
            let mut last = u32::MAX;
            for i in lo..hi {
                let b = raw[i];
                if b != last && b as usize != a {
                    raw[w] = b;
                    w += 1;
                    last = b;
                }
            }
            adj_start[a + 1] = w as u32;
        }
        drop(start);
        raw.truncate(w);
        raw.shrink_to_fit();
        let adj = raw;

        // Shortest and longest edge (bounds for the search queues).
        let edge_range = par_map(n, |a| {
            let mut r = [f32::INFINITY, 0f32];
            for &b in &adj[adj_start[a] as usize..adj_start[a + 1] as usize] {
                let l = (dir[a] - dir[b as usize]).length() as f32;
                r = [r[0].min(l), r[1].max(l)];
            }
            r
        });
        let min_edge = edge_range.iter().fold(f32::INFINITY, |m, r| m.min(r[0]));
        let max_edge = edge_range.iter().fold(0f32, |m, r| m.max(r[1]));
        Grid { res, n, dir, weight, adj_start, adj, min_edge, max_edge, sample_node, node_sample }
    }

    #[inline]
    pub fn neighbors(&self, a: usize) -> &[u32] {
        &self.adj[self.adj_start[a] as usize..self.adj_start[a + 1] as usize]
    }

    /// Typical spacing between neighbouring nodes, radians.
    pub fn spacing_rad(&self) -> f64 {
        std::f64::consts::FRAC_PI_2 / self.res as f64
    }

    /// Decode a node's representative sample into (face, i, j).
    pub fn node_face_ij(&self, node: usize) -> (CubeFace, u32, u32) {
        let per_edge = (self.res + 1) as usize;
        let s = self.node_sample[node] as usize;
        let face = CubeFace::ALL[s / (per_edge * per_edge)];
        let r = s % (per_edge * per_edge);
        (face, (r % per_edge) as u32, (r / per_edge) as u32)
    }

    /// Evaluate `f(dir)` on a coarser cube grid and interpolate (bilinear in face coordinates)
    /// up to every node. Use this for smooth, low-frequency fields: it is much cheaper than
    /// evaluating noise at every node.
    pub fn coarse_field(&self, coarse_res: u32, f: impl Fn(DVec3) -> f32 + Sync) -> Vec<f32> {
        let coarse = self.coarse(coarse_res, f);
        par_map(self.n, |node| self.lookup(&coarse, node))
    }

    /// Only the coarse grid of [`Grid::coarse_field`]. Use [`Grid::lookup`] to read it at a
    /// node. This saves memory when a field is needed in just one pass: no per-node array.
    pub fn coarse(&self, coarse_res: u32, f: impl Fn(DVec3) -> f32 + Sync) -> CubeField<f32> {
        let mut coarse = CubeField::<f32>::new(coarse_res);
        let per = (coarse_res + 1) as usize;
        let data = par_map(6 * per * per, |idx| {
            let (face, r) = (CubeFace::ALL[idx / (per * per)], idx % (per * per));
            f(coarse.sample_dir(face, (r % per) as u32, (r / per) as u32))
        });
        coarse.data = data;
        coarse
    }

    /// Value of a coarse grid, interpolated bilinearly (in face coordinates) at a node.
    #[inline]
    pub fn lookup(&self, coarse: &CubeField<f32>, node: usize) -> f32 {
        let rc = coarse.res as f64;
        let (face, i, j) = self.node_face_ij(node);
        let x = i as f64 * rc / self.res as f64;
        let y = j as f64 * rc / self.res as f64;
        let i0 = (x.floor() as u32).min(coarse.res - 1);
        let j0 = (y.floor() as u32).min(coarse.res - 1);
        let (fx, fy) = ((x - i0 as f64) as f32, (y - j0 as f64) as f32);
        let at = |ii: u32, jj: u32| coarse.data[coarse.index(face, ii, jj)];
        let top = at(i0, j0) + (at(i0 + 1, j0) - at(i0, j0)) * fx;
        let bot = at(i0, j0 + 1) + (at(i0 + 1, j0 + 1) - at(i0, j0 + 1)) * fx;
        top + (bot - top) * fy
    }

    /// Area-weighted fraction of nodes for which `pred` holds.
    pub fn fraction(&self, mut pred: impl FnMut(usize) -> bool) -> f64 {
        (0..self.n).filter(|&a| pred(a)).map(|a| self.weight[a] as f64).sum()
    }

    /// Grow regions outward from `seeds` at once (multi-source shortest paths). Several seeds may
    /// belong to one region (`seed_region[i]`). Every node ends up with the region that
    /// reaches it cheapest. Stepping onto a node
    /// costs `edge length * cost[node] / speed[seed]`. Because every node inherits its parent's
    /// region, each region is **connected by construction**.
    ///
    /// Each region can also be *elongated*: `axis[k]` is the pole of a great circle, and steps
    /// along that circle's direction cost less than steps across it by a factor of up to
    /// `1 + stretch[k]`. This makes plates long and thin rather than round.
    pub fn grow_regions(
        &self,
        seeds: &[u32],
        seed_region: &[u16],
        speed: &[f32],
        axis: &[DVec3],
        stretch: &[f32],
        cost: &[f32],
    ) -> Vec<u16> {
        // Bounds on one step's cost, to size the bucket queue.
        let cost_min = cost.iter().copied().fold(f32::INFINITY, f32::min);
        let cost_max = cost.iter().copied().fold(0.0, f32::max);
        let speed_min = speed.iter().copied().fold(f32::INFINITY, f32::min);
        let speed_max = speed.iter().copied().fold(0.0, f32::max);
        let stretch_max = stretch.iter().copied().fold(0.0, f32::max);
        let step_min = self.min_edge * cost_min / speed_max;
        let step_max = self.max_edge * cost_max / speed_min * (1.0 + stretch_max);

        let mut best = vec![f32::INFINITY; self.n];
        let mut done = Bits::new(self.n); // set once a node has been expanded with its current label
        let mut owner = vec![u16::MAX; self.n];
        let mut queue = BucketQueue::new(2.0 * step_min, step_max);
        for (&s, &k) in seeds.iter().zip(seed_region) {
            best[s as usize] = 0.0;
            owner[s as usize] = k;
            queue.push(0.0, s);
        }
        while let Some(a) = queue.pop() {
            let a = a as usize;
            if done.test_and_set(a) {
                continue; // already expanded with its current label
            }
            let d = best[a];
            let k = owner[a] as usize;
            let (lo, hi) = (self.adj_start[a] as usize, self.adj_start[a + 1] as usize);
            for e in lo..hi {
                let b = self.adj[e] as usize;
                let edge_len = (self.dir[b] - self.dir[a]).length() as f32;
                let mut step = edge_len * 0.5 * (cost[a] + cost[b]) / speed[k];
                if stretch[k] > 0.0 {
                    let t = axis[k].cross(self.dir[a]);
                    let tl = t.length();
                    if tl > 1e-6 {
                        let along = ((self.dir[b] - self.dir[a]).dot(t) / (tl * edge_len as f64)).abs();
                        step *= 1.0 + stretch[k] * (1.0 - along as f32);
                    }
                }
                let nd = d + step;
                if nd < best[b] {
                    best[b] = nd;
                    owner[b] = k as u16;
                    done.clear(b);
                    queue.push(nd, b as u32);
                }
            }
        }
        owner
    }

    /// For every node, the nearest of `sources` and the true great-circle distance to it.
    ///
    /// A plain graph shortest-path would measure distance along grid lines and be up to ~8% off
    /// (diamond-shaped isolines). Instead we propagate "which source is nearest" through the
    /// graph, and measure the exact straight-line distance to that source. Nodes farther than
    /// `max_rad` from every source are left with `src = u32::MAX` and `dist = max_rad`.
    pub fn nearest_source(&self, sources: &[u32], max_rad: f32) -> Nearest {
        let max_chord2 = {
            let c = 2.0 * (max_rad as f64 * 0.5).sin();
            (c * c * 1.0001 + 1e-6) as f32
        };
        let mut best = vec![f32::INFINITY; self.n]; // squared chord distance
        let mut done = Bits::new(self.n);
        let mut src = vec![u32::MAX; self.n];
        // The queue is ordered by chord length (not squared): one edge changes it by at most
        // the edge length.
        let mut queue = BucketQueue::new(2.0 * self.min_edge, self.max_edge);
        for &s in sources {
            best[s as usize] = 0.0;
            src[s as usize] = s;
            queue.push(0.0, s);
        }
        while let Some(a) = queue.pop() {
            let a = a as usize;
            if done.test_and_set(a) {
                continue;
            }
            let s = src[a];
            let sd = self.dir[s as usize];
            for &b in self.neighbors(a) {
                let b = b as usize;
                let c = self.dir[b] - sd;
                let nd = c.length_squared() as f32;
                if nd < best[b] && nd <= max_chord2 {
                    best[b] = nd;
                    src[b] = s;
                    done.clear(b);
                    queue.push(nd.sqrt(), b as u32);
                }
            }
        }
        let dist_rad = best
            .iter()
            .map(|&c2| if c2.is_finite() { 2.0 * (0.5 * c2.sqrt()).min(1.0).asin() } else { max_rad })
            .collect();
        Nearest { dist_rad, src }
    }
}

/// A fixed-size set of bit flags.
struct Bits(Vec<u64>);

impl Bits {
    fn new(n: usize) -> Self {
        Bits(vec![0; n.div_ceil(64)])
    }

    /// Set bit `i`; returns whether it was already set.
    #[inline]
    fn test_and_set(&mut self, i: usize) -> bool {
        let (w, m) = (i / 64, 1u64 << (i % 64));
        let was = self.0[w] & m != 0;
        self.0[w] |= m;
        was
    }

    #[inline]
    fn clear(&mut self, i: usize) {
        self.0[i / 64] &= !(1u64 << (i % 64));
    }
}

/// A priority queue that pops keys *approximately* in order: keys are cut into buckets of a
/// fixed width, and inside one bucket the order is arbitrary. That is much faster than a
/// binary heap on millions of nodes. It is only used for shortest-path style searches that
/// re-expand a node whenever it gets a better label, so the final labels are the same as with
/// an exact queue (up to ties); a small ordering error only costs a few repeated expansions.
///
/// `max_step` is the largest amount by which a newly pushed key can exceed the key of the
/// entry being expanded. The buckets form a ring of just over `max_step / width` slots.
struct BucketQueue {
    inv_width: f32,
    ring: Vec<Vec<u32>>,
    /// Index (not wrapped) of the bucket being drained.
    cur: usize,
    len: usize,
    /// Value of `cur` for which the slot has been sorted.
    sorted_for: usize,
}

impl BucketQueue {
    fn new(width: f32, max_step: f32) -> Self {
        let slots = (max_step / width).ceil() as usize + 3;
        Self { inv_width: 1.0 / width, ring: vec![Vec::new(); slots], cur: 0, len: 0, sorted_for: usize::MAX }
    }

    fn push(&mut self, key: f32, node: u32) {
        // Keys below the current bucket (rounding, or a source that gets closer) are
        // processed right away.
        let k = ((key * self.inv_width) as usize).max(self.cur);
        debug_assert!(k - self.cur < self.ring.len(), "key too far ahead");
        let slot = k % self.ring.len();
        self.ring[slot].push(node);
        self.len += 1;
    }

    fn pop(&mut self) -> Option<u32> {
        while self.len > 0 {
            let slot = self.cur % self.ring.len();
            if self.sorted_for != self.cur {
                // Handle a bucket in node order (descending, since we pop from the back):
                // neighbouring nodes are then expanded one after another, which keeps
                // the memory accesses local.
                self.ring[slot].sort_unstable_by(|a, b| b.cmp(a));
                self.sorted_for = self.cur;
            }
            if let Some(node) = self.ring[slot].pop() {
                self.len -= 1;
                return Some(node);
            }
            self.cur += 1;
        }
        None
    }
}

/// Result of [`Grid::nearest_source`].
pub struct Nearest {
    /// Distance to the nearest source, radians of arc.
    pub dist_rad: Vec<f32>,
    /// Node index of that source, or `u32::MAX` when out of range.
    pub src: Vec<u32>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn node_count_and_weights() {
        let g = Grid::new(16);
        assert_eq!(g.n, 6 * 16 * 16 + 2);
        let s: f32 = g.weight.iter().sum();
        assert!((s - 1.0).abs() < 1e-4);
        // Symmetric graph.
        for a in 0..g.n {
            for &b in g.neighbors(a) {
                assert!(g.neighbors(b as usize).contains(&(a as u32)));
            }
        }
    }

    #[test]
    fn nearest_source_is_accurate() {
        let g = Grid::new(48);
        let s = 1234u32;
        let near = g.nearest_source(&[s], std::f32::consts::PI);
        for a in (0..g.n).step_by(97) {
            let truth = g.dir[a].dot(g.dir[s as usize]).clamp(-1.0, 1.0).acos();
            assert!((near.dist_rad[a] as f64 - truth).abs() < 2e-3, "{a}");
        }
    }
}
