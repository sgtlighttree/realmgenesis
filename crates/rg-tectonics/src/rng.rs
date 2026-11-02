//! A tiny deterministic random number generator (SplitMix64).
//!
//! We avoid the `rand` crate on purpose: the workspace wants no extra dependencies, and the
//! output of this generator must never change between compiler or library versions.

use glam::DVec3;

pub struct Rng(u64);

impl Rng {
    /// A generator for one purpose ("stream") of one world seed. Different streams are
    /// independent, so adding a new use of randomness never reshuffles the existing ones.
    pub fn new(seed: u64, stream: u64) -> Self {
        let mut r = Rng(seed ^ stream.wrapping_mul(0x9E37_79B9_7F4A_7C15) ^ 0xD1B5_4A32_D192_ED03);
        r.next_u64();
        r.next_u64();
        r
    }

    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in [0, 1).
    pub fn f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// Uniform in [lo, hi).
    pub fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.f64()
    }

    /// Standard normal (Box-Muller).
    pub fn gauss(&mut self) -> f64 {
        let u1 = self.f64().max(1e-12);
        let u2 = self.f64();
        (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
    }

    /// A uniformly distributed unit vector.
    pub fn unit_vec(&mut self) -> DVec3 {
        let z = self.range(-1.0, 1.0);
        let phi = self.range(0.0, std::f64::consts::TAU);
        let r = (1.0 - z * z).sqrt();
        DVec3::new(r * phi.cos(), r * phi.sin(), z)
    }

    /// A uniform index in 0..n.
    pub fn index(&mut self, n: usize) -> usize {
        ((self.f64() * n as f64) as usize).min(n - 1)
    }
}
