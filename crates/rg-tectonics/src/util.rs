//! Small shared helpers: easing curves, noise wrappers and a bisection solver.

use glam::DVec3;
use rg_noise::Fbm;

#[inline]
pub fn smoothstep(a: f64, b: f64, x: f64) -> f64 {
    let t = ((x - a) / (b - a)).clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

/// exp(-(x/w)^2): 1 at x = 0, ~0.37 at x = w. A smooth bell used for most profiles.
#[inline]
pub fn bell(x: f64, w: f64) -> f64 {
    let t = x / w;
    (-t * t).exp()
}

#[inline]
pub fn lerp(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}

/// Normalised fBm in about [-1, 1]. `freq` is in lattice cells per planet radius.
/// (`rg_noise` normalises by the amplitude sum, so the *effective* amplitude is much smaller
/// than the nominal one; callers scale by measured output.)
#[inline]
pub fn fbm(dir: DVec3, salt: u32, freq: f64, octaves: u32, gain: f64) -> f64 {
    Fbm { octaves, frequency: freq, lacunarity: 2.0, gain }.sample(dir, salt)
}

/// Ridged noise in [0, 1]: sharp crests where the underlying noise crosses zero.
#[inline]
pub fn ridged(dir: DVec3, salt: u32, freq: f64, octaves: u32, gain: f64) -> f64 {
    Fbm { octaves, frequency: freq, lacunarity: 2.0, gain }.ridged(dir, salt)
}

/// fBm mapped to roughly [0, 1] (clamped), handy as a multiplier.
#[inline]
pub fn fbm01(dir: DVec3, salt: u32, freq: f64, octaves: u32) -> f64 {
    (0.5 + 0.5 * 1.7 * fbm(dir, salt, freq, octaves, 0.5)).clamp(0.0, 1.0)
}

/// A smooth 3-vector displacement field (three independent noises).
pub fn warp3(dir: DVec3, salt: u32, freq: f64, octaves: u32) -> DVec3 {
    DVec3::new(
        fbm(dir, salt, freq, octaves, 0.5),
        fbm(dir + DVec3::new(31.7, -12.3, 7.1), salt.wrapping_add(101), freq, octaves, 0.5),
        fbm(dir + DVec3::new(-18.2, 22.9, 41.3), salt.wrapping_add(202), freq, octaves, 0.5),
    )
}

/// Find `t` such that the weighted share of `values > t` equals `target_above`, by bisection.
/// (Multiplicative steps oscillate on flat plateaus; bisection always converges. We return the
/// closest iterate.)
pub fn solve_threshold(values: &[f32], weights: &[f32], target_above: f64) -> f32 {
    let (mut lo, mut hi) = values
        .iter()
        .fold((f32::MAX, f32::MIN), |(l, h), &v| (l.min(v), h.max(v)));
    let mut best = (f64::MAX, lo);
    for _ in 0..48 {
        let mid = 0.5 * (lo + hi);
        let frac: f64 = values
            .iter()
            .zip(weights)
            .filter(|(v, _)| **v > mid)
            .map(|(_, w)| *w as f64)
            .sum();
        let err = (frac - target_above).abs();
        if err < best.0 {
            best = (err, mid);
        }
        if frac > target_above {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo < 1e-4 {
            break;
        }
    }
    best.1
}
