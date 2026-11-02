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
- **Uplift U:** from tectonic process (distance to convergent boundaries, crust type, closing speed; `uplift.rs`), not from elevation.
- **Routing:** Rho8 (Fairfield & Leymarie 1991), diagonals found topologically; initial roughness ≤ 20 m (scaled by local relief); lithology log-normal K, σ_lnK 0.3, 250 km correlation.
- **Grid:** runs at the preset resolution (`ErosionParams::res`, 0 = tectonics res). If capped lower, the *change* is upsampled and PF-filled on the fine grid.
- **Interpolation:** `CubeField::sample_drained`: bilinear with any diagonal valley carved straight (bilinear dammed 9.5% of mountain land).
- **Defaults:** 30 Myr, 20 steps (same laws as 50), K = 1e-5, m = 0.45.
- **Verified (tests, `tests/geomorphology.rs`):**
  - all land drains to the sea;
  - uplift-normalised concavity = m to ±0.03 (0.4500 measured), the analytic steady state;
  - Hack's law h ∈ [0.45, 0.7] (0.54 measured; Schott 2023 reports 0.54–0.58);
  - apparent concavity ∈ [0.3, 0.8];
  - per-basin hypsometric integral ∈ [0.2, 0.6];
  - coast unchanged;
  - deterministic.

## 2. Hypsometric calibration

At 10–50 km cells, K and U are underdetermined, so absolute heights aren't meaningful. After erosion, land heights are **area-weighted quantile-mapped onto a measured, smoothed, ice-free land hypsometric curve** (× the "mountains" setting; mean 656 m, median 382 m), strictly increasing as stored f32.
- The map is monotonic, so every drainage direction, divide and valley is unchanged. Solver stats are measured before the remap.
- **Measured (seed 3):** mean 841 m (Earth 840), median 429 m (~440), > 2 km 10.3% (~9%), > 4 km 2.8% (~2%).
- **Curve:** approximate ETOPO1-class numbers (table in `solver.rs`).

## 3. Detail: the drained detail network (`rg-terrain/src/detail_net.rs`, s01)

Design and proofs: `docs/design/drained-detail-network.md` (v1.1 plus later changes in the s01 log). The old stateless detail (hills + Johansen gully filter) is kept only as the fallback for worlds without hydrology. It was "steep or drained, not both": hills put 13–16% of mountain land in closed depressions, and the gully filter's stripes scored coherence 0.48–0.91.

- **Surface.** Land height = min(cap, valley wedges).
  - Cap U = T + R(T): T is the triangle-split macro (provably drained), R = 0.3·min(T, 3 km).
  - Wedge = channel bed z(t) + strictly increasing hillslope H(d): a flat floor with a C1 hinge, tanh convexity, and an adaptive reach wall.
  - A local minimum of a min is a minimum of the active wedge, whose only minimum is its junction on a parent, so the surface drains to the sea **by construction**.
- **Hierarchy.**
  - Level 0 = macro river edges, curved by centripetal Catmull-Rom through the node chain.
  - Level k = one jittered point per cell of a lattice 3^k finer (15% thinned, no stubs).
  - Each point drains to the **lowest coarser wedge** within 1.2 parent spacings ("the hillslope it sits on"), or, by hillslope-position rank, to a lower-rank sibling (dendritic side-branching).
  - Junction angle 44–80° per channel (Seybold et al. 2017). Channels are 3-piece meanders.
- **Physics per channel.**
  - k_sn = U/K (exact at stream-power steady state, n = 1; U is the committed uplift).
  - Channel gradient S = k_sn·A^−θ with A ≈ one cell, capped at 0.5.
  - Head incised a hollow below the coarser surface, aiming at T − 0.33 R (mean-neutral against the macro).
  - Valley walls S_c·(1 − e^(−k_sn/100)) with S_c = tan 33°, after Ouimet et al. 2009; hilltop share 0.3 → 0.9 with erosion rate (Roering et al. 2007).
  - Floodplain width 30 m·√(A/km²) (hydraulic geometry), narrowed to gorges with uplift.
- **LOD.**
  - Per-channel presence (LOD spacing × [0.75, 1.33]); floors band-limited to the cut-off.
  - Ancestor-consistent fade: each channel carries its fading ancestors' bed offsets, so junctions sit on the parent's effective bed at every zoom.
  - Level 0 never fades.
- **Evaluation.**
  - Per patch (one per tile): gather with margins derived from reach, search and rank steps.
  - Exact block pruning; hierarchical exact candidate lists per parent cell.
  - Memos budgeted by contents. Tiles ~85–240 ms on CPU (web streams slowly; the GPU port is the fix).
- **Measured** (`probe_realism`, seed 3, res 320, final surface):
  - closed depressions **0.00%** in all boxes;
  - mountain 90 m slopes 13.8–20° mean, p99 33–54°; lowland 2.1° mean;
  - coherence 0.22–0.49;
  - mean detail − macro ≈ 0 (`probe_hypsometry`).
- **Refuted along the way:** smooth-min crease rounding (pits at every confluence); the triangle-split as display interpolant (facets); the nearest-parent attach (buried children, mesas); per-piece fade ramps (pits).

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

## 4. Bathymetry (`rg-tectonics/src/relief.rs`, s01)

Water depth on oceanic crust = GDH1 age–depth (Stein & Stein 1992) − dynamic topography − plateau lift − isostatic sediment lift, then trenches, arcs and seamounts from the plate boundaries.
- **Age:** distance from the nearest ridge over the half spreading rate, then area-weighted quantile-mapped onto Earth's triangular area–age distribution, F(t) = 1 − (1 − t/180 Myr)² (Parsons 1982; mean 60 Myr).
- **Sediment** (GlobSed, Straume et al. 2019, Table 3): pelagic 680·(1 − e^(−t/50 Myr)) m (404 m area mean, as GlobSed's deep ocean), plus a terrigenous wedge 4700·e^(−x/100 km) m offshore of the continent–ocean boundary (2,050 m mean within 200 km). Lift = (ρm − ρs)/(ρm − ρw) ≈ 0.49 per metre.
- **Dynamic topography:** residual depth with 350 m std (Hoggard et al. 2016: ±1 km at ~1000 km wavelengths).
- **Oceanic plateaus:** top 5.1% of an independent noise field (Harris et al. 2014: 5.11% of the oceans), raised 2.2 km, flat-topped.
- **Trenches:** Gaussian e-folding 40 km (a 75 km one gave 6–8× Earth's > 7 km area).
- **Measured** (`examples/probe_bathymetry.rs`, seeds 3 / 11, res 320) against Earth (Costello et al. 2010; ETOPO1 mean 3688 m):

  | | seed 3 | seed 11 | Earth |
  |---|---|---|---|
  | mean depth | 3652 m | 3687 m | 3688 m |
  | 0–200 m | 8.9% | 8.0% | 7.3% |
  | 200–1000 m | 4.7% | 4.5% | 4.4% |
  | 1–4 km | 31.2% | 34.1% | 35.2% |
  | 4–6 km | 54.3% | 52.3% | 52.1% |
  | 6–7 km | 0.5% | 0.6% | 0.9% |
  | > 7 km | 0.4% | 0.5% | 0.13% |

  > 7 km is resolution-limited (a trench is about one cell wide at 31 km cells).
- **Not yet:** abyssal-hill fabric aligned with the spreading direction, fracture zones, submarine canyons continuing rivers across the shelf, and streamed seafloor detail (oceans render flat at the surface).

## Known limits

- **Hydrological consistency of detail:** solved by the drained detail network (§3). The old fallback detail still dams.
- **Cost:** network tiles cost ~85–240 ms per 65² tile on the CPU, so the web build streams slowly. The GPU port (M1b) is the fix; the per-patch segment buffers port directly.
- **Peaks:** detail can add 1–1.5 km on top of the highest macro peaks (max about 9–10 km against Earth's 8.8).
