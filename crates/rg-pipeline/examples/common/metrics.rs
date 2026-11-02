//! Terrain metrics on any height function over a local box, shared by `probe_realism` (our
//! surface) and `earth_reference` (real DEMs), so both are scored by identical code.
//!
//! `h(x, y)`: height in m at (x, y) m east/north of the box centre; values ≤ 0 are sea.
#![allow(dead_code)]
use std::cmp::Ordering;
use std::collections::BinaryHeap;

const NB: [(i64, i64); 8] = [(-1, -1), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1)];

#[derive(PartialEq)]
struct It(f64, usize);
impl Eq for It {}
impl Ord for It {
    fn cmp(&self, o: &Self) -> Ordering {
        o.0.total_cmp(&self.0).then(o.1.cmp(&self.1))
    }
}
impl PartialOrd for It {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}

/// Sample an n×n grid at spacing `sp` centred on the box.
pub fn grid(h: &dyn Fn(f64, f64) -> f64, n: usize, sp: f64) -> Vec<f64> {
    let half = 0.5 * (n as f64 - 1.0) * sp;
    (0..n * n).map(|k| h(-half + (k % n) as f64 * sp, -half + (k / n) as f64 * sp)).collect()
}

/// Priority-Flood fill with the grid edge (and the sea) as outlets.
pub fn fill(g: &[f64], n: usize) -> Vec<f64> {
    let mut filled = g.to_vec();
    let mut done = vec![false; n * n];
    let mut heap = BinaryHeap::new();
    for k in 0..n * n {
        let (i, j) = (k % n, k / n);
        if i == 0 || j == 0 || i == n - 1 || j == n - 1 || g[k] <= 0.0 {
            done[k] = true;
            heap.push(It(g[k], k));
        }
    }
    while let Some(It(hc, c)) = heap.pop() {
        let (ci, cj) = ((c % n) as i64, (c / n) as i64);
        for (di, dj) in NB {
            let (ni, nj) = (ci + di, cj + dj);
            if ni < 0 || nj < 0 || ni >= n as i64 || nj >= n as i64 {
                continue;
            }
            let k = nj as usize * n + ni as usize;
            if done[k] {
                continue;
            }
            done[k] = true;
            filled[k] = filled[k].max(hc);
            heap.push(It(filled[k], k));
        }
    }
    filled
}

/// Share of land (%) in closed depressions deeper than 0.5 m, and their mean depth.
pub fn depressions(g: &[f64], n: usize) -> (f64, f64) {
    let f = fill(g, n);
    let land = g.iter().filter(|&&x| x > 0.0).count().max(1);
    let depths: Vec<f64> = (0..n * n).map(|k| f[k] - g[k]).filter(|&d| d > 0.5).collect();
    (depths.len() as f64 / land as f64 * 100.0, depths.iter().sum::<f64>() / depths.len().max(1) as f64)
}

/// 90 m central-difference slopes (degrees) at 2000 pseudo-random land points in ±half:
/// [mean, p50, p90, p99].
pub fn slopes(h: &dyn Fn(f64, f64) -> f64, half: f64) -> [f64; 4] {
    let mut rng = 12345u64;
    let mut next = || {
        rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (rng >> 11) as f64 / (1u64 << 53) as f64
    };
    let d = 45.0;
    let mut deg = Vec::new();
    for _ in 0..2000 {
        let (x, y) = ((next() * 2.0 - 1.0) * half, (next() * 2.0 - 1.0) * half);
        if h(x, y) <= 0.0 {
            continue;
        }
        let gx = (h(x + d, y) - h(x - d, y)) / (2.0 * d);
        let gy = (h(x, y + d) - h(x, y - d)) / (2.0 * d);
        deg.push((gx * gx + gy * gy).sqrt().atan().to_degrees());
    }
    deg.sort_by(f64::total_cmp);
    let q = |f: f64| deg.get(((deg.len().max(1) - 1) as f64 * f) as usize).copied().unwrap_or(0.0);
    [deg.iter().sum::<f64>() / deg.len().max(1) as f64, q(0.5), q(0.9), q(0.99)]
}

/// Structure-tensor coherence (mean over 25×25-sample windows) and gradient-weighted aspect
/// concentration, on an m×m grid at spacing sp.
pub fn orientation(g: &[f64], m: usize, sp: f64) -> (f64, f64) {
    let (mut coh_sum, mut coh_n) = (0.0, 0usize);
    let (mut ax, mut ay, mut aw) = (0.0, 0.0, 0.0);
    let win = 25;
    for wy in (1..m - win - 1).step_by(win) {
        for wx in (1..m - win - 1).step_by(win) {
            let (mut jxx, mut jyy, mut jxy) = (0.0, 0.0, 0.0);
            for j in wy..wy + win {
                for i in wx..wx + win {
                    let gx = (g[j * m + i + 1] - g[j * m + i - 1]) / (2.0 * sp);
                    let gy = (g[(j + 1) * m + i] - g[(j - 1) * m + i]) / (2.0 * sp);
                    jxx += gx * gx;
                    jyy += gy * gy;
                    jxy += gx * gy;
                    let mag = (gx * gx + gy * gy).sqrt();
                    if mag > 1e-6 {
                        let th = gy.atan2(gx);
                        ax += mag * (2.0 * th).cos();
                        ay += mag * (2.0 * th).sin();
                        aw += mag;
                    }
                }
            }
            let tr = jxx + jyy;
            if tr > 1e-12 {
                coh_sum += ((jxx - jyy).powi(2) + 4.0 * jxy * jxy).sqrt() / tr;
                coh_n += 1;
            }
        }
    }
    (coh_sum / coh_n.max(1) as f64, (ax * ax + ay * ay).sqrt() / aw.max(1e-12))
}

/// 1D power-spectrum exponent β (P ∝ k^−β) from the rows and columns of an m×m grid at
/// spacing sp: detrended, Hann-windowed, averaged, fitted over wavelengths 4 sp … m·sp/4.
pub fn spectrum_beta(g: &[f64], m: usize, sp: f64) -> f64 {
    spectrum(g, m, sp).0
}

/// β plus the log-residual (dB) of each octave-binned band around the power-law fit,
/// as (wavelength m, residual dB): lattice spikes show as positive residuals.
pub fn spectrum(g: &[f64], m: usize, sp: f64) -> (f64, Vec<(f64, f64)>) {
    let len = m.min(128);
    let mut power = vec![0.0f64; len / 2];
    let mut count = 0usize;
    let mut profile = |p: &[f64]| {
        // Detrend (least-squares line), Hann window.
        let nf = p.len() as f64;
        let (mx, my) = ((nf - 1.0) / 2.0, p.iter().sum::<f64>() / nf);
        let sxy: f64 = p.iter().enumerate().map(|(i, v)| (i as f64 - mx) * (v - my)).sum();
        let sxx: f64 = (0..p.len()).map(|i| (i as f64 - mx).powi(2)).sum();
        let b = sxy / sxx.max(1e-12);
        let w: Vec<f64> = p
            .iter()
            .enumerate()
            .map(|(i, v)| (v - my - b * (i as f64 - mx)) * (0.5 - 0.5 * (std::f64::consts::TAU * i as f64 / (nf - 1.0)).cos()))
            .collect();
        for (k, pw) in power.iter_mut().enumerate().skip(1) {
            let (mut re, mut im) = (0.0, 0.0);
            for (i, v) in w.iter().enumerate() {
                let a = std::f64::consts::TAU * (k * i) as f64 / nf;
                re += v * a.cos();
                im -= v * a.sin();
            }
            *pw += re * re + im * im;
        }
        count += 1;
    };
    for j in (0..m).step_by(4) {
        profile(&g[j * m..j * m + len]);
    }
    for i in (0..m).step_by(4) {
        let col: Vec<f64> = (0..len).map(|j| g[j * m + i]).collect();
        profile(&col);
    }
    let _ = count;
    // Fit log P vs log k for wavelengths between 4 sp and len·sp/4.
    let (mut xs, mut ys) = (Vec::new(), Vec::new());
    for (k, &pw) in power.iter().enumerate() {
        if k < 4 || k > len / 4 || pw <= 0.0 {
            continue;
        }
        let kk = k as f64 / (len as f64 * sp);
        xs.push(kk.ln());
        ys.push(pw.ln());
    }
    let slope = fit_slope(&xs, &ys);
    let n = xs.len().max(1) as f64;
    let (mx, my) = (xs.iter().sum::<f64>() / n, ys.iter().sum::<f64>() / n);
    // Residuals in half-octave bins.
    let mut bins: Vec<(f64, f64, usize)> = Vec::new();
    for (x, y) in xs.iter().zip(&ys) {
        let r = y - (my + slope * (x - mx));
        let b = (x / std::f64::consts::LN_2 * 2.0).floor();
        match bins.iter_mut().find(|e| e.0 == b) {
            Some(e) => {
                e.1 += r;
                e.2 += 1;
            }
            None => bins.push((b, r, 1)),
        }
    }
    let res = bins
        .iter()
        .map(|(b, r, c)| (1.0 / (2f64).powf((b + 0.5) / 2.0), 10.0 * (r / *c as f64) / std::f64::consts::LN_10))
        .collect();
    (-slope, res)
}

/// Slope–area concavity θ (S ∝ A^−θ) on an m×m grid at spacing sp: D8 on the filled
/// surface, drainage area by accumulation, fitted over channel cells with A in
/// [0.1, 10] km² (above the hillslope rollover). Also returns the cell count fitted.
pub fn slope_area_theta(g: &[f64], m: usize, sp: f64) -> (f64, usize) {
    let f = fill(g, m);
    // Tiny tie-break so flats in fills route (ε along the fill order is not available here;
    // use the raw surface for slopes and the filled one for routing).
    let mut recv = vec![usize::MAX; m * m];
    let mut drop = vec![0.0f64; m * m];
    for j in 1..m - 1 {
        for i in 1..m - 1 {
            let k = j * m + i;
            let mut best = (0.0, usize::MAX);
            for (di, dj) in NB {
                let kk = (j as i64 + dj) as usize * m + (i as i64 + di) as usize;
                let dist = sp * ((di * di + dj * dj) as f64).sqrt();
                let s = (f[k] - f[kk]) / dist;
                if s > best.0 {
                    best = (s, kk);
                }
            }
            recv[k] = best.1;
            drop[k] = best.0;
        }
    }
    let mut order: Vec<usize> = (0..m * m).collect();
    order.sort_by(|&a, &b| f[b].total_cmp(&f[a]).then(a.cmp(&b)));
    let mut area = vec![sp * sp; m * m];
    for &k in &order {
        if recv[k] != usize::MAX {
            area[recv[k]] += area[k];
        }
    }
    let (mut xs, mut ys) = (Vec::new(), Vec::new());
    for k in 0..m * m {
        let a_km2 = area[k] * 1e-6;
        if (0.1..=10.0).contains(&a_km2) && drop[k] > 1e-4 && g[k] > 0.0 {
            xs.push(a_km2.ln());
            ys.push(drop[k].ln());
        }
    }
    (-fit_slope(&xs, &ys), xs.len())
}

fn fit_slope(xs: &[f64], ys: &[f64]) -> f64 {
    let n = xs.len().max(1) as f64;
    let (mx, my) = (xs.iter().sum::<f64>() / n, ys.iter().sum::<f64>() / n);
    let sxy: f64 = xs.iter().zip(ys).map(|(x, y)| (x - mx) * (y - my)).sum();
    let sxx: f64 = xs.iter().map(|x| (x - mx).powi(2)).sum();
    sxy / sxx.max(1e-12)
}
