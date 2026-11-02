// This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
// If a copy of the MPL was not distributed with this file, You can obtain one at
// https://mozilla.org/MPL/2.0/.
//
// Port of the "Advanced Terrain Erosion Filter" and Phacelle noise by Rune Skovbo
// Johansen (2025; https://www.shadertoy.com/view/wXcfWn,
// https://blog.runevision.com/2026/03/fast-and-gorgeous-erosion-filter.html), via the
// Rust port in `bevy_erosion_filter` by Korbin de Man (MPL-2.0). Lineage: Fewes (2023),
// Clay John (2020). This file stays under MPL-2.0; the rest of RealmGenesis is unaffected.
//
// Changes from upstream for planet scale (RealmGenesis, 2026-11-02):
// - f64 throughout, and an integer PCG3D lattice hash instead of IQ's
//   `fract(x·y·(x+y))` hash, which loses all precision at planet-scale coordinates.
// - Octaves are band-limited for LOD: partial sums are reported at three cutoffs
//   (this level, parent, child) in one pass.
// - On the sphere, a 3D generalisation (`phacelle_3d`, `erosion_filter_sphere`): jittered
//   cells on a 3D lattice and 3D tangent vectors for slope and gully direction. A first
//   triplanar attempt produced cliffs where the planes' independent patterns met.
// - Continuity: sign(s) in the gully-steering update is replaced by tanh(6·s), and the
//   triplanar plane cut-off is 1e-4 instead of 0.01; both removed visible height steps.
//   The sphere version also uses a Phacelle phase offset of 0 (upstream 0.25 cycles), so
//   it is symmetric when the gully direction flips sign.
// - Hillslope diffusion: each octave's height contribution is damped by the steady-state
//   response of linear diffusion, 1 / (1 + (λ_d / λ)²), so short wavelengths are smoothed
//   as soil creep smooths real hillslopes. Gully *directions* are unaffected.

use glam::{DVec2, DVec3};

const TAU: f64 = std::f64::consts::TAU;

/// 3D Phacelle cell jitter and squared kernel support (see `phacelle_3d`).
const JITTER: f64 = 0.22;
const SUPPORT2: f64 = 1.25 * 1.25;
/// Steepness of the smooth sign used for gully steering (see `erosion_filter`).
const SIGN_SHARPNESS: f64 = 6.0;

/// Lattice jitter in [-1, 1]² for an integer cell (replaces IQ's hash2).
#[inline]
fn hash2(cell: DVec2, salt: u32) -> DVec2 {
    let (a, b, _) = rg_noise::pcg3d(cell.x.floor() as i64 as u32, cell.y.floor() as i64 as u32, salt);
    DVec2::new(a as f64 / u32::MAX as f64 * 2.0 - 1.0, b as f64 / u32::MAX as f64 * 2.0 - 1.0)
}

/// Phacelle noise: a field of stripes aligned with `norm_dir`, blended smoothly between
/// jittered cells. Returns (cos wave, sin wave, side direction × freq·τ).
pub fn phacelle(p: DVec2, norm_dir: DVec2, freq: f64, offset_cycles: f64, normalization: f64, salt: u32) -> (f64, f64, DVec2) {
    let side_dir = DVec2::new(-norm_dir.y, norm_dir.x) * freq * TAU;
    let offset = offset_cycles * TAU;
    let p_int = p.floor();
    let p_frac = p - p_int;
    let mut phase = DVec2::ZERO;
    let mut weight_sum = 0.0;
    for i in -1..=2 {
        for j in -1..=2 {
            let g = DVec2::new(i as f64, j as f64);
            let v = p_frac - g - hash2(p_int + g, salt) * 0.5;
            let d2 = v.dot(v);
            // exp(−2d²) − 1/90 reaches zero at d² ≈ 2.25: skip those cells before the trig.
            if d2 >= 2.25 {
                continue;
            }
            let w = ((-d2 * 2.0).exp() - 0.011_11).max(0.0);
            weight_sum += w;
            let a = v.dot(side_dir) + offset;
            let (sn, cs) = a.sin_cos();
            phase += DVec2::new(cs, sn) * w;
        }
    }
    let interp = phase / weight_sum.max(1e-12);
    let mag = (1.0 - normalization).max(interp.length());
    let n = interp / mag.max(1e-12);
    (n.x, n.y, side_dir)
}

/// Filter parameters (upstream defaults, see Johansen's post for their meaning).
#[derive(Clone, Copy, Debug)]
pub struct FilterParams {
    pub scale: f64,
    pub strength: f64,
    pub gully_weight: f64,
    pub detail: f64,
    pub rounding: [f64; 4],
    pub onset: [f64; 4],
    pub assumed_slope: [f64; 2],
    pub cell_scale: f64,
    pub normalization: f64,
    pub lacunarity: f64,
    pub gain: f64,
}

impl Default for FilterParams {
    fn default() -> Self {
        Self {
            scale: 0.15,
            strength: 0.22,
            gully_weight: 0.5,
            detail: 1.5,
            rounding: [0.1, 0.3, 1.0, 2.0],
            onset: [1.25, 1.25, 2.8, 1.5],
            assumed_slope: [0.7, 1.0],
            cell_scale: 0.7,
            normalization: 0.5,
            lacunarity: 2.0,
            gain: 0.5,
        }
    }
}

/// Height deltas at three LOD cutoffs (in filter units).
#[derive(Clone, Copy, Debug, Default)]
pub struct Banded {
    pub at_level: f64,
    pub at_parent: f64,
    pub at_child: f64,
}

fn clamp01(t: f64) -> f64 {
    t.clamp(0.0, 1.0)
}
fn pow_inv(t: f64, p: f64) -> f64 {
    1.0 - (1.0 - clamp01(t)).powf(p)
}
fn ease_out(t: f64) -> f64 {
    let v = 1.0 - clamp01(t);
    1.0 - v * v
}
fn smooth_start(t: f64, s: f64) -> f64 {
    if t >= s { t - 0.5 * s } else { 0.5 * t * t / s.max(1e-12) }
}
fn mix(a: f64, b: f64, t: f64) -> f64 {
    a + (b - a) * t
}

/// The erosion filter at 2D point `p` (filter units) on a base surface with height `h`
/// and gradient `slope`. Octave k has wavelength `scale·cell_scale / lacunarity^k`;
/// octaves are evaluated while that wavelength ≥ `min_wl_child` and summed into the three
/// bands by their own cutoffs. Returns height deltas only.
#[allow(clippy::too_many_arguments)]
pub fn erosion_filter(p: DVec2, h: f64, slope: DVec2, fade_target_in: f64, prm: &FilterParams, cutoffs: [f64; 3], diffusion_wl: f64, salt: u32) -> Banded {
    let [min_wl, min_wl_parent, min_wl_child] = cutoffs;
    let mut strength = prm.strength * prm.scale;
    let mut fade_target = fade_target_in.clamp(-1.0, 1.0);
    let mut freq = 1.0 / (prm.scale * prm.cell_scale);
    let slope_len = slope.length().max(1e-10);
    let mut rounding_mult = 1.0;
    let rounding_in = mix(prm.rounding[1], prm.rounding[0], clamp01(fade_target + 0.5)) * prm.rounding[2];
    let mut combi = ease_out(smooth_start(slope_len * prm.onset[0], rounding_in * prm.onset[0]));
    let mut gully = slope.lerp(slope / slope_len * prm.assumed_slope[0], prm.assumed_slope[1]);
    let mut out = Banded::default();
    let mut dh = 0.0;
    let _ = h;
    let mut octave = 0u32;
    while 1.0 / freq >= min_wl_child && octave < 24 {
        let wl_here = 1.0 / freq;
        let damp = 1.0 / (1.0 + (diffusion_wl / wl_here).powi(2));
        if damp < 0.05 {
            // Diffusion has smoothed everything this short and shorter away.
            break;
        }
        if combi < 1.0e-8 {
            // The slope mask has closed (flat ground): each further octave just adds
            // fade_target·strength (the "faded" value with combi = 0). Exact, no noise.
            let contrib = fade_target * strength * damp;
            dh += contrib;
            out.at_child = dh;
            if wl_here >= min_wl {
                out.at_level = dh;
            }
            if wl_here >= min_wl_parent {
                out.at_parent = dh;
            }
            strength *= prm.gain;
            freq *= prm.lacunarity;
            octave += 1;
            continue;
        }
        let g = gully.normalize_or(DVec2::X);
        let (c, s, side) = phacelle(p * freq, g, prm.cell_scale, 0.25, prm.normalization, salt.wrapping_add(octave));
        let p_zw = side * -freq;
        // Upstream uses sign(s), which flips the next octave's gully direction abruptly
        // wherever s crosses zero: a visible height discontinuity. tanh(k·s) is the same
        // away from zero and continuous through it.
        gully += (s * SIGN_SHARPNESS).tanh() * p_zw * strength * prm.gully_weight;
        let octave_h = c;
        let faded = mix(fade_target, octave_h * prm.gully_weight, combi);
        let contrib = faded * strength * damp;
        dh += contrib;
        let wl = 1.0 / freq;
        out.at_child = dh;
        if wl >= min_wl {
            out.at_level = dh;
        }
        if wl >= min_wl_parent {
            out.at_parent = dh;
        }
        fade_target = faded;
        let rounding_oct = mix(prm.rounding[1], prm.rounding[0], clamp01(c + 0.5)) * rounding_mult;
        let new_mask = ease_out(smooth_start(s.abs() * prm.onset[1], rounding_oct * prm.onset[1]));
        combi = pow_inv(combi, prm.detail) * new_mask;
        strength *= prm.gain;
        freq *= prm.lacunarity;
        rounding_mult *= prm.rounding[3];
        octave += 1;
    }
    out
}

/// 3D Phacelle: stripes aligned with the tangent direction `dir3`, from jittered cells on a
/// 3D lattice. On a sphere this is seamless (no projection planes). `up` is the surface
/// normal; the stripes vary across `side = up × dir3`. Returns (cos, sin, side·freq·τ).
pub fn phacelle_3d(p: DVec3, up: DVec3, dir3: DVec3, freq: f64, offset_cycles: f64, normalization: f64, salt: u32) -> (f64, f64, DVec3) {
    let side = up.cross(dir3).normalize_or(up.any_orthonormal_vector()) * freq * TAU;
    let offset = offset_cycles * TAU;
    // Kernel support R = 1.25, cell jitter ≤ 0.22 per axis:
    // - coverage: the farthest a point can be from its nearest cell centre is
    //   √3·(0.5 + 0.22) ≈ 1.247 < R, so the weight sum never vanishes (a vanishing sum made
    //   the output snap to 0: a discontinuity);
    // - window: a contributing cell is within R + 0.22 = 1.47 < 1.5 per axis of p, i.e.
    //   within ±1 of round(p), so a 3×3×3 window is exact (27 cells; ~8 contribute).
    let c0 = p.round();
    let mut phase = DVec2::ZERO;
    let mut weight_sum = 0.0;
    for i in -1..=1 {
        for j in -1..=1 {
            for k in -1..=1 {
                let c = c0 + DVec3::new(i as f64, j as f64, k as f64);
                let (a, b, cc) = rg_noise::pcg3d(c.x as i64 as u32, c.y as i64 as u32, (c.z as i64 as u32) ^ salt);
                let jitter = DVec3::new(a as f64, b as f64, cc as f64) / u32::MAX as f64 * 2.0 - DVec3::ONE;
                let v = p - c - jitter * JITTER;
                let d2 = v.dot(v);
                if d2 >= SUPPORT2 {
                    continue;
                }
                // Smooth compact weight (1 − d²/R²)³: zero with zero slope at the edge.
                let q = 1.0 - d2 / SUPPORT2;
                let w = q * q * q;
                weight_sum += w;
                let (sn, cs) = (v.dot(side) + offset).sin_cos();
                phase += DVec2::new(cs, sn) * w;
            }
        }
    }
    if weight_sum <= 0.0 {
        return (0.0, 0.0, side);
    }
    let interp = phase / weight_sum;
    // Smooth normalisation: upstream's max(1 − normalization, |i|) has a kink along
    // |i| = 1 − normalization, which shows as crease lines on km-scale octaves.
    // √(|i|² + c²) behaves the same at both ends and is smooth everywhere.
    let c = 1.0 - normalization;
    let mag = (interp.length_squared() + c * c).sqrt();
    let n = interp / mag.max(1e-12);
    (n.x, n.y, side)
}

/// The erosion filter on a sphere, in 3D: `pos` in filter units, `up` the surface normal,
/// `grad3` the base surface's tangent gradient (rise per run). Same structure as
/// [`erosion_filter`] with every 2D vector replaced by a 3D tangent vector.
#[allow(clippy::too_many_arguments)]
pub fn erosion_filter_sphere(pos: DVec3, up: DVec3, grad3: DVec3, h: f64, fade_target_in: f64, prm: &FilterParams, cutoffs: [f64; 3], diffusion_wl: f64, salt: u32) -> Banded {
    let _ = h;
    let [min_wl, min_wl_parent, min_wl_child] = cutoffs;
    let tangent = |v: DVec3| v - up * v.dot(up);
    let slope = tangent(grad3);
    let mut strength = prm.strength * prm.scale;
    let mut fade_target = fade_target_in.clamp(-1.0, 1.0);
    let mut freq = 1.0 / (prm.scale * prm.cell_scale);
    let slope_len = slope.length().max(1e-10);
    let mut rounding_mult = 1.0;
    let rounding_in = mix(prm.rounding[1], prm.rounding[0], clamp01(fade_target + 0.5)) * prm.rounding[2];
    let mut combi = ease_out(smooth_start(slope_len * prm.onset[0], rounding_in * prm.onset[0]));
    let mut gully = slope.lerp(slope / slope_len * prm.assumed_slope[0], prm.assumed_slope[1]);
    let mut out = Banded::default();
    let mut dh = 0.0;
    let mut octave = 0u32;
    while 1.0 / freq >= min_wl_child && octave < 24 {
        let wl_here = 1.0 / freq;
        let damp = 1.0 / (1.0 + (diffusion_wl / wl_here).powi(2));
        if damp < 0.05 {
            break;
        }
        // Threshold must be tiny: with octave amplitudes of ~1 km, 1e-4 left 0.14 m steps.
        let contrib = if combi < 1.0e-8 {
            fade_target * strength * damp
        } else {
            let dir3 = tangent(gully).normalize_or(up.any_orthonormal_vector());
            // Offset 0 (upstream: 0.25 cycles). With an offset, cos(±x + φ) differs when the
            // gully direction flips sign (across ridges and valleys), stepping the height.
            // With 0 the stripes are symmetric under the flip, so the result is continuous.
            let (c, s, side) = phacelle_3d(pos * freq, up, dir3, prm.cell_scale, 0.0, prm.normalization, salt.wrapping_add(octave));
            let p_zw = side * -freq;
            gully = tangent(gully + (s * SIGN_SHARPNESS).tanh() * p_zw * strength * prm.gully_weight);
            let faded = mix(fade_target, c * prm.gully_weight, combi);
            fade_target = faded;
            let rounding_oct = mix(prm.rounding[1], prm.rounding[0], clamp01(c + 0.5)) * rounding_mult;
            let new_mask = ease_out(smooth_start(s.abs() * prm.onset[1], rounding_oct * prm.onset[1]));
            combi = pow_inv(combi, prm.detail) * new_mask;
            faded * strength * damp
        };
        dh += contrib;
        out.at_child = dh;
        if wl_here >= min_wl {
            out.at_level = dh;
        }
        if wl_here >= min_wl_parent {
            out.at_parent = dh;
        }
        strength *= prm.gain;
        freq *= prm.lacunarity;
        rounding_mult *= prm.rounding[3];
        octave += 1;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finite_and_band_ordered() {
        let prm = FilterParams::default();
        for &(x, y) in &[(0.42, 0.31), (1234.5, -987.25), (3.0e5, 1.7e5)] {
            let b = erosion_filter(DVec2::new(x, y), 0.5, DVec2::new(0.3, -0.2), 0.0, &prm, [0.01, 0.02, 0.005], 0.0, 7);
            assert!(b.at_level.is_finite() && b.at_parent.is_finite() && b.at_child.is_finite());
        }
    }

    #[test]
    fn flat_ground_gets_no_gullies() {
        let prm = FilterParams::default();
        let b = erosion_filter(DVec2::new(10.3, 4.2), 0.5, DVec2::ZERO, 0.0, &prm, [0.0, 0.0, 0.0], 0.0, 7);
        assert!(b.at_child.abs() < 1e-6, "{}", b.at_child);
    }
}

#[cfg(test)]
mod continuity {
    use super::*;

    fn max_jump(f: impl Fn(f64) -> f64, x0: f64, x1: f64, step: f64) -> (f64, f64) {
        let mut worst = (0.0, x0);
        let mut prev = f(x0);
        let mut x = x0 + step;
        while x < x1 {
            let v = f(x);
            if (v - prev).abs() > worst.0 {
                worst = ((v - prev).abs(), x);
            }
            prev = v;
            x += step;
        }
        worst
    }

    #[test]
    fn phacelle_3d_is_continuous() {
        let up = DVec3::new(0.3, 0.4, 0.866).normalize();
        let t = up.any_orthonormal_vector();
        let d = up.cross(t);
        let f = |x: f64| phacelle_3d(DVec3::new(10.0, 20.0, 30.0) + t * x, up, d, 0.7, 0.0, 0.5, 9).0;
        let (j, at) = max_jump(f, 0.0, 20.0, 1.0e-4);
        assert!(j < 2.0e-3, "phacelle jump {j} at {at}");
    }

    #[test]
    fn sphere_filter_is_continuous_with_fixed_slope() {
        let up = DVec3::new(0.3, 0.4, 0.866).normalize();
        let t = up.any_orthonormal_vector();
        let prm = FilterParams::default();
        let f = |x: f64| {
            let p = DVec3::new(10.0, 20.0, 30.0) + t * x;
            erosion_filter_sphere(p, up, t * 0.4, 0.0, 0.2, &prm, [0.0, 0.0, 0.0001], 0.0, 9).at_child
        };
        let (j, at) = max_jump(f, 0.0, 3.0, 1.0e-5);
        assert!(j < 2.0e-3, "filter jump {j} at {at}");
    }
}
