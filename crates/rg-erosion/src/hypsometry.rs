//! Hypsometric calibration: map the eroded land heights onto a target land hypsometric curve.
//!
//! # Where the target curve comes from (s02, 2026-11-02)
//!
//! The earlier curve was approximate ETOPO1-class numbers (mean 840 m, ~9% above 2 km). It
//! has two problems at the solver's 20-50 km cell size:
//!
//! 1. Those numbers come from ~1-2 km DEMs and include the **ice-sheet surfaces** of
//!    Antarctica and Greenland (mean ~2.2 km over 14 million km²), which a fantasy planet
//!    with the usual temperate land does not have;
//! 2. imposed on 20-50 km cells, a 1-2 km-resolution distribution has too heavy a tail,
//!    because averaging over a cell flattens peaks, and the streamed detail layer then adds
//!    its own variance on top.
//!
//! The new curve was **measured**, not recalled: GMT's `earth_relief_10m` grid, which is
//! SRTM15+ V2.7 (Tozer et al. 2019, doi:10.1029/2019EA000658) low-passed with a 52.4 km
//! Gaussian filter and sampled every 10 arc-minutes (~18 km). That is the same smoothing
//! scale as our cells. Land = z > 0, weighted by cos(latitude), Antarctica (south of 60° S)
//! and Greenland masked out. Results:
//!
//! | quantity | all land incl. ice | ex-Antarctica, ex-Greenland (used) |
//! |---|---|---|
//! | land area | 148.5 M km² | 132.0 M km² |
//! | mean | 794 m | 656 m |
//! | median | 423 m | 382 m |
//! | above 1 km | 25.7% | 19.9% |
//! | above 2 km | 10.7% | 5.6% |
//! | above 4 km | 1.7% | 1.9% |
//! | highest cell mean | 6118 m | 6118 m |
//!
//! The "all land" mean and median agree with the literature values (about 800 m and 420 m),
//! which is the cross-check that the grid and the weighting are right.
//!
//! **Caveats, stated plainly.** (a) This is the present-day ice-free *surface* of 87-90% of
//! the land. True ice-free Earth would add Antarctic and Greenland *bedrock* (mostly low,
//! with mountains in the Transantarctic and Gamburtsev ranges); we did not measure it
//! (Bedmap2 was not available offline), and it would lower the mean a little. (b) The curve
//! below is read off the measured quantiles and rounded; treat every number as ±10%.
//! (c) The smoothing scale (52 km filter) is somewhat coarser than the solver's finest grid
//! (res 512, ~20 km), so the tail is, if anything, slightly light at the finest resolution.

/// Target land hypsometric curve: `(fraction of land area above, height in m)`.
/// Read off quantiles of the smoothed ice-free-land measurement described in the module docs.
pub(crate) const LAND_HYPSOMETRY: [(f64, f64); 19] = [
    (1.000, 0.0),
    (0.900, 67.0),
    (0.800, 133.0),
    (0.700, 204.0),
    (0.600, 290.0),
    (0.500, 382.0),
    (0.400, 495.0),
    (0.300, 682.0),
    (0.200, 996.0),
    (0.150, 1180.0),
    (0.100, 1450.0),
    (0.070, 1760.0),
    (0.050, 2120.0),
    (0.030, 3015.0),
    (0.020, 3950.0),
    (0.010, 4690.0),
    (0.005, 5000.0),
    (0.001, 5350.0),
    (0.000, 6200.0),
];

/// Quantile-map land heights onto the target curve (× `scale`). `weight[i]` is node `i`'s
/// share of the sphere's area, so the *area* of land above a height (not the node count) is
/// what matches the curve. Monotonic, so the order of heights — and with it every drainage
/// direction, divide and valley — is kept; only how high things are is redistributed. Ties
/// break by node index (deterministic).
pub(crate) fn remap(h: &mut [f64], ocean: &[bool], weight: &[f64], scale: f64) {
    let mut land: Vec<usize> = (0..h.len()).filter(|&i| !ocean[i]).collect();
    if land.len() < 2 {
        return;
    }
    land.sort_by(|&a, &b| h[a].total_cmp(&h[b]).then(a.cmp(&b)));
    let total: f64 = land.iter().map(|&i| weight[i]).sum();
    let mut below = 0.0;
    for &i in &land {
        // Fraction of land *area* strictly above this node (midpoint rule within the node).
        let above = 1.0 - (below + 0.5 * weight[i]) / total;
        below += weight[i];
        h[i] = (height_at(above) * scale).max(1.0);
    }
}

/// Height at which `frac_above` of the land area lies above (linear between table rows).
pub(crate) fn height_at(frac_above: f64) -> f64 {
    let t = &LAND_HYPSOMETRY;
    for w in t.windows(2) {
        let (f0, h0) = w[0];
        let (f1, h1) = w[1];
        if frac_above <= f0 && frac_above >= f1 {
            let u = (f0 - frac_above) / (f0 - f1).max(1e-12);
            return h0 + (h1 - h0) * u;
        }
    }
    if frac_above > 1.0 { 0.0 } else { t[t.len() - 1].1 }
}

/// Mean height of the target curve (area-weighted), m. For tests and docs.
pub fn target_mean_m() -> f64 {
    LAND_HYPSOMETRY.windows(2).map(|w| (w[0].0 - w[1].0) * 0.5 * (w[0].1 + w[1].1)).sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn curve_is_monotonic_and_has_the_measured_statistics() {
        for w in LAND_HYPSOMETRY.windows(2) {
            assert!(w[0].0 > w[1].0 && w[0].1 < w[1].1);
        }
        let mean = target_mean_m();
        assert!((600.0..760.0).contains(&mean), "mean {mean}");
        assert!((height_at(0.5) - 382.0).abs() < 1.0);
        assert!(height_at(0.056) > 1900.0 && height_at(0.056) < 2200.0); // ~5.6% above 2 km
    }

    #[test]
    fn weighted_ranks_follow_area_not_count() {
        // Ten equal-height-spaced land nodes; the top node holds 40% of the area, so its
        // rank band is the top 40%: it must land near the 0.8-above quantile, not 0.05.
        let n = 10;
        let mut h: Vec<f64> = (0..n).map(|i| i as f64 + 1.0).collect();
        let mut w = vec![0.6 / 9.0; n];
        w[n - 1] = 0.4;
        let ocean = vec![false; n];
        remap(&mut h, &ocean, &w, 1.0);
        assert!((h[n - 1] - height_at(0.2)).abs() < 1.0, "top node {}", h[n - 1]);
        // Unweighted ranks would have given the top node height_at(0.05).
        assert!(h[n - 1] < height_at(0.05) * 0.6);
    }
}
