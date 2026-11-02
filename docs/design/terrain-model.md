# Terrain model (s01, 2026-11-02)

How RealmGenesis builds terrain, and how it's checked against Earth. Research behind it:
`docs/research/2026-11-02-terrain-erosion.md`.

## Pipeline

```
tectonics (plates, crust, boundary relief)                  rg-tectonics
  → fluvial erosion (stream-power law, implicit)             rg-erosion
  → hypsometric calibration (quantile map to Earth's curve)  rg-erosion
  → climate (tier A)                                         rg-climate
  → streamed detail per point (hills + erosion filter)       rg-terrain
```

## 1. Macro erosion (`rg-erosion`)

`∂h/∂t = U − K·A^m·S`, solved implicitly in O(n) (Braun & Willett 2013).
- **Routing:** D8 on the seam-free 8-neighbour sphere graph, with Priority-Flood+ε fills for closed basins (Barnes 2014).
- **Uplift U:** proportional to tectonic relief (`uplift_max` × (h/h95)^0.8 + base 0.15 mm/yr).
- **Grid:** runs at min(res, 512). Finer fields get the *change* upsampled (Schott 2023: blend changes, not surfaces).
- **Defaults:** 30 Myr, 50 steps, K = 1e-5, m = 0.45.
- **Verified (tests, `tests/geomorphology.rs`):**
  - all land drains to the sea;
  - uplift-normalised concavity = m to ±0.03 (0.4500 measured), the analytic steady state;
  - Hack's law h ∈ [0.45, 0.7] (0.54 measured; Schott 2023 reports 0.54–0.58);
  - apparent concavity ∈ [0.3, 0.8];
  - per-basin hypsometric integral ∈ [0.2, 0.6];
  - coast unchanged;
  - deterministic.

## 2. Hypsometric calibration

At 10–50 km cells, K and U are underdetermined, so absolute heights aren't meaningful. After erosion, land heights are **quantile-mapped onto Earth's land hypsometric curve** (× the "mountains" setting).
- The map is monotonic, so every drainage direction, divide and valley is unchanged. Solver stats are measured before the remap.
- **Measured (seed 3):** mean 841 m (Earth 840), median 429 m (~440), > 2 km 10.3% (~9%), > 4 km 2.8% (~2%).
- **Curve:** approximate ETOPO1-class numbers (table in `solver.rs`).

## 3. Detail (`rg-terrain`, per point, stateless, LOD-banded)

`h = macro + hills + gullies`, with the coast rule: no change within 150 m of sea level, sign preserved, so the coastline is exact at every zoom.
- **Hills:** fBm from 2 macro cells down to 4 km. Constant slope per octave below 4 km, amplitude saturating above it: `A = S·λ/2π / (1 + λ/4 km)`, S = 1.5.
  - Scaled by **relief**, a function of the **macro slope** measured over 20 km, not of elevation. Plateaus stay gentle; incising belts get rugged. Plains keep 12%.
- **Gullies:** Johansen's erosion filter (MPL-2.0 port, `erosion_filter.rs`), generalised to 3D on the sphere. It carves branching valleys down the slope of macro + hills, from 4 km down to about 30 m.
  - Hillslope diffusion damps short wavelengths: 1/(1 + (150 m/λ)²).
  - Strength 2.0 × relief.
- **Measured 90 m slopes** (`examples/probe_slope.rs`):
  - mountain belts: mean 15–20°, p90 30–37°, p99 48–54°;
  - plains: about 0.6°.
  - Earth's mountain belts typically sit at 15–30° mean (threshold hillslopes about 30–35°: Montgomery & Brandon 2002, Larsen & Montgomery 2012); plains a few degrees.
- **Continuity tests** (`erosion_filter.rs`, `continuity` module): Phacelle and the full filter, each to < 2e-3 per 1e-4 step.

### Artifacts found and fixed (keep these in mind)

| Symptom | Cause | Fix |
|---|---|---|
| Vertical cliffs along great circles | Triplanar planes with independent patterns | 3D Phacelle (3D lattice, tangent vectors) |
| Steps across ridges and valleys | `sign(s)` steering; Phacelle phase offset 0.25 isn't flip-symmetric | `tanh(6s)`; offset 0 |
| mm-scale steps of ~0.1 m | f32 bilinear macro sampling, ×filter sensitivity | `CubeField::sample_f64` |
| Output snapping to 0 in gaps | 3D kernel support didn't cover the jittered lattice | R = 1.25, jitter 0.22 (√3·0.72 < R), 3×3×3 window |
| Long sharp V-creases | Valley rounding 0, i.e. an `|slope|` kink in the mask | rounding [0.1, 0.3, 1.0, 2.0] |
| Long crisp arcs | Ridged noise `(1 − |n|)²` kinks at zero-contours | Ridged noise removed |
| Smooth sheets, sparse scratches | Base too smooth for the filter (macro slopes ~0.04) | Rough hill spectrum below the macro scale |
| Peaks too high, plateaus too rugged | Relief scaled by elevation | Relief from macro slope |
| ~2× Earth heights | p99 rescale to tectonic heights | Hypsometric calibration |
| Dark walls at LOD borders | Child tiles carry extra octaves | CDLOD geomorph to the parent mesh surface |

## Known limits

- **Hydrological consistency of detail:** hills and gullies follow local slopes, but between 5 and 30 km they aren't routed against the macro network, so closed depressions can exist in the detail. Rivers drawn on the layer come from the macro graph.
- **Cost:** about 40 ms per 65² tile at deep levels on the CPU, so the web build streams slowly. The GPU port (M1b) is the fix: WGSL with a 3D integer lattice.
- **Peaks:** detail can add 1–1.5 km on top of the highest macro peaks (max about 9–10 km against Earth's 8.8).
