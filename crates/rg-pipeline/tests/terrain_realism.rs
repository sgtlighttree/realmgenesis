//! Final terrain (macro + streamed detail) must look like Earth's in distribution:
//! mountain belts steep, plains gentle, few closed depressions. See
//! docs/design/terrain-model.md for the reference numbers and sources.

use rg_pipeline::{generate, WorldParams};

struct Box2 {
    s: rg_terrain::TerrainSampler,
    r: f64,
}

impl Box2 {
    fn h(&self, lat0: f64, lon0: f64, x: f64, y: f64, wl: f64) -> f64 {
        let (la, lo) = (lat0 + y / self.r, lon0 + x / (self.r * lat0.cos()));
        self.s.height(glam::DVec3::new(la.cos() * lo.cos(), la.cos() * lo.sin(), la.sin()), wl).0
    }

    /// Mean 90 m slope (degrees) over a grid of points in a ±half box.
    fn mean_slope(&self, lat0: f64, lon0: f64, half: f64) -> f64 {
        let d = 45.0;
        let mut sum = 0.0;
        let k = 24;
        for j in 0..k {
            for i in 0..k {
                let x = -half + 2.0 * half * (i as f64 + 0.5) / k as f64;
                let y = -half + 2.0 * half * (j as f64 + 0.5) / k as f64;
                let gx = (self.h(lat0, lon0, x + d, y, 40.0) - self.h(lat0, lon0, x - d, y, 40.0)) / (2.0 * d);
                let gy = (self.h(lat0, lon0, x, y + d, 40.0) - self.h(lat0, lon0, x, y - d, 40.0)) / (2.0 * d);
                sum += (gx * gx + gy * gy).sqrt().atan().to_degrees();
            }
        }
        sum / (k * k) as f64
    }
}

#[test]
fn mountains_steep_plains_gentle_few_pits() {
    let mut p = WorldParams::default();
    p.tectonics.seed = 3;
    p.tectonics.res = 192;
    let w = generate(&p);
    let f = &w.tectonics.height_m;
    let b = Box2 { s: rg_pipeline::terrain_sampler(&w), r: p.tectonics.radius_m };

    // The most rugged land point (largest macro step to a neighbour) and a flat lowland one.
    let res = f.res;
    let (mut rugged, mut flat) = ((0.0f32, glam::DVec3::Z), (f32::MAX, glam::DVec3::Z));
    for (k, face, i, j, d) in f.iter_samples() {
        if i == 0 || j == 0 || i == res || j == res || f.data[k] <= 50.0 {
            continue;
        }
        let nb = [f.index(face, i + 1, j), f.index(face, i - 1, j), f.index(face, i, j + 1), f.index(face, i, j - 1)];
        let step = nb.iter().map(|&n| (f.data[n] - f.data[k]).abs()).fold(0.0f32, f32::max);
        if step > rugged.0 {
            rugged = (step, d);
        }
        if f.data[k] > 100.0 && f.data[k] < 400.0 && step < flat.0 {
            flat = (step, d);
        }
    }
    let ll = |d: glam::DVec3| (d.z.asin(), d.y.atan2(d.x));
    let (rl, ro) = ll(rugged.1);
    let (fl, fo) = ll(flat.1);
    let mountain = b.mean_slope(rl, ro, 20_000.0);
    let plain = b.mean_slope(fl, fo, 20_000.0);
    // Earth: mountain belts ~15–30° mean at 90 m; plains a few degrees.
    assert!((10.0..35.0).contains(&mountain), "mountain mean slope {mountain}°");
    assert!(plain < 4.0, "plain mean slope {plain}°");

    // Closed depressions at 500 m spacing in the rugged box: few (natural fluvial terrain).
    let (n, sp) = (60usize, 500.0);
    let hs: Vec<f64> = (0..n * n).map(|k| b.h(rl, ro, (k % n) as f64 * sp - 15_000.0, (k / n) as f64 * sp - 15_000.0, sp)).collect();
    let mut pits = 0;
    for j in 1..n - 1 {
        for i in 1..n - 1 {
            let c = hs[j * n + i];
            let lowest = [(-1i64, -1i64), (0, -1), (1, -1), (-1, 0), (1, 0), (-1, 1), (0, 1), (1, 1)]
                .iter()
                .all(|&(di, dj)| hs[(j as i64 + dj) as usize * n + (i as i64 + di) as usize] > c);
            if c > 0.0 && lowest {
                pits += 1;
            }
        }
    }
    let frac = pits as f64 / ((n - 2) * (n - 2)) as f64;
    assert!(frac < 0.01, "pit fraction {frac}");
}
