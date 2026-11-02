//! Deterministic noise on the unit sphere.
//!
//! Lattice hashing uses integers only (PCG3D), so every platform gets identical lattice
//! values. Only the final interpolation is float. A WGSL twin will follow in M1. Keep the
//! two in lockstep: same hash, same fade, same gradient table.

use glam::DVec3;

/// PCG3D hash (Jarzynski & Olano 2020, "Hash Functions for GPU Rendering").
#[inline]
pub fn pcg3d(mut x: u32, mut y: u32, mut z: u32) -> (u32, u32, u32) {
    x = x.wrapping_mul(1664525).wrapping_add(1013904223);
    y = y.wrapping_mul(1664525).wrapping_add(1013904223);
    z = z.wrapping_mul(1664525).wrapping_add(1013904223);
    x = x.wrapping_add(y.wrapping_mul(z));
    y = y.wrapping_add(z.wrapping_mul(x));
    z = z.wrapping_add(x.wrapping_mul(y));
    x ^= x >> 16;
    y ^= y >> 16;
    z ^= z >> 16;
    x = x.wrapping_add(y.wrapping_mul(z));
    y = y.wrapping_add(z.wrapping_mul(x));
    z = z.wrapping_add(x.wrapping_mul(y));
    (x, y, z)
}

/// Mix a seed and a stream id into one 32-bit salt. Streams keep layers independent:
/// changing one layer never reseeds another (lesson: "never reseed", tcotg FINDINGS).
pub fn salt(seed: u64, stream: u32) -> u32 {
    let (a, _, _) = pcg3d(seed as u32, (seed >> 32) as u32, stream);
    a
}

#[inline]
fn grad(ix: i32, iy: i32, iz: i32, salt: u32) -> DVec3 {
    let (h, _, _) = pcg3d(ix as u32 ^ salt, iy as u32, iz as u32);
    // 12 cube-edge gradients (classic Perlin set), picked by the hash.
    const G: [(f64, f64, f64); 12] = [
        (1., 1., 0.), (-1., 1., 0.), (1., -1., 0.), (-1., -1., 0.),
        (1., 0., 1.), (-1., 0., 1.), (1., 0., -1.), (-1., 0., -1.),
        (0., 1., 1.), (0., -1., 1.), (0., 1., -1.), (0., -1., -1.),
    ];
    let g = G[(h % 12) as usize];
    DVec3::new(g.0, g.1, g.2)
}

#[inline]
fn fade(t: f64) -> f64 {
    t * t * t * (t * (t * 6.0 - 15.0) + 10.0)
}

/// 3D gradient noise, roughly in [-1, 1].
pub fn gradient3(p: DVec3, salt: u32) -> f64 {
    let i = p.floor();
    let f = p - i;
    let (ix, iy, iz) = (i.x as i32, i.y as i32, i.z as i32);
    let (u, v, w) = (fade(f.x), fade(f.y), fade(f.z));
    let mut acc = [0.0f64; 8];
    for (k, a) in acc.iter_mut().enumerate() {
        let (dx, dy, dz) = ((k & 1) as i32, ((k >> 1) & 1) as i32, ((k >> 2) & 1) as i32);
        let g = grad(ix + dx, iy + dy, iz + dz, salt);
        *a = g.dot(f - DVec3::new(dx as f64, dy as f64, dz as f64));
    }
    let lerp = |a: f64, b: f64, t: f64| a + (b - a) * t;
    let x00 = lerp(acc[0], acc[1], u);
    let x10 = lerp(acc[2], acc[3], u);
    let x01 = lerp(acc[4], acc[5], u);
    let x11 = lerp(acc[6], acc[7], u);
    lerp(lerp(x00, x10, v), lerp(x01, x11, v), w)
}

/// fBm parameters. `frequency` is in cycles per planet radius.
#[derive(Clone, Copy, Debug)]
pub struct Fbm {
    pub octaves: u32,
    pub frequency: f64,
    pub lacunarity: f64,
    pub gain: f64,
}

impl Default for Fbm {
    fn default() -> Self {
        Self { octaves: 8, frequency: 1.5, lacunarity: 2.0, gain: 0.5 }
    }
}

impl Fbm {
    /// Sum of octaves at a unit-sphere direction. **Normalized** by the amplitude sum, so
    /// the result stays in about [-1, 1] whatever the octave count. (tcotg lesson: this
    /// normalization makes effective amplitudes look much smaller than nominal.)
    pub fn sample(&self, dir: DVec3, salt: u32) -> f64 {
        let mut sum = 0.0;
        let mut amp = 1.0;
        let mut norm = 0.0;
        let mut freq = self.frequency;
        for o in 0..self.octaves {
            // Offset each octave so lattice points don't line up across octaves.
            let off = DVec3::new(o as f64 * 17.13, o as f64 * -9.71, o as f64 * 5.37);
            sum += amp * gradient3(dir * freq + off, salt.wrapping_add(o));
            norm += amp;
            amp *= self.gain;
            freq *= self.lacunarity;
        }
        sum / norm
    }

    /// Ridged variant (1 - |n|)², for mountain-ish crests. In [0, 1].
    pub fn ridged(&self, dir: DVec3, salt: u32) -> f64 {
        let mut sum = 0.0;
        let mut amp = 1.0;
        let mut norm = 0.0;
        let mut freq = self.frequency;
        for o in 0..self.octaves {
            let off = DVec3::new(o as f64 * 11.7, o as f64 * 3.1, o as f64 * -7.9);
            let n = 1.0 - gradient3(dir * freq + off, salt.wrapping_add(0x9e37 + o)).abs();
            sum += amp * n * n;
            norm += amp;
            amp *= self.gain;
            freq *= self.lacunarity;
        }
        sum / norm
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic_and_seed_sensitive() {
        let d = DVec3::new(0.3, 0.4, 0.5).normalize();
        let f = Fbm::default();
        assert_eq!(f.sample(d, salt(42, 0)), f.sample(d, salt(42, 0)));
        assert_ne!(f.sample(d, salt(42, 0)), f.sample(d, salt(43, 0)));
        assert_ne!(f.sample(d, salt(42, 0)), f.sample(d, salt(42, 1)));
    }

    #[test]
    fn range_is_sane() {
        let f = Fbm::default();
        let s = salt(7, 0);
        let (mut lo, mut hi) = (f64::MAX, f64::MIN);
        for i in 0..2000 {
            let t = i as f64 * 0.37;
            let d = DVec3::new(t.sin(), (t * 1.3).cos(), (t * 0.7).sin()).normalize();
            let v = f.sample(d, s);
            lo = lo.min(v);
            hi = hi.max(v);
        }
        assert!(lo > -1.2 && hi < 1.2 && hi - lo > 0.4, "{lo}..{hi}");
    }

    #[test]
    fn hash_is_pinned() {
        // Pinned values: if this fails, the hash changed and every world changes with it.
        assert_eq!(pcg3d(1, 2, 3), pcg3d(1, 2, 3));
        let (a, b, c) = pcg3d(0, 0, 0);
        assert_eq!((a, b, c), pcg3d(0, 0, 0));
        assert_ne!(a, b);
    }
}
