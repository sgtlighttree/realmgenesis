# Terrain erosion research (2026-11-02)

Deep web survey by a subagent (Haiku). Items marked verified were read at the source that session; the rest are marked. Full report condensed below. **Licence notes are load-bearing:** don't port CC BY-NC-SA shader code (Clay John, Gavoronoise, Shadertoy default) or GPL code (fastscapelib, richdem).

## Decision

- **Stage A (macro, at generation): Braun & Willett 2013 implicit stream-power law (SPL)** on a routed drainage tree, with D8 steepest-descent receivers + Priority-Flood+ε depression filling (Barnes 2014). O(n) per step. Expected: seconds at 6M cells, single-threaded; fine on wasm32. CPU first for determinism (committed fields). m/n ≈ 0.45–0.5 (Lague 2014 review).
  - GPU later: Schott et al. 2023 (TOG 42(5):162) iterative drainage area, MIT code at github.com/H-Schott/StreamPowerErosion; FastFlow (Jain et al., CGF 2024) for O(log n) routing.
  - **Pitfall (Schott 2023, Fig. 7):** blending two eroded fields destroys dendrites. Blend *uplift*, not heights.
  - **Pitfall:** flats after filling need ε gradients (fastscapelib and Landlab document this).
- **Stage B (detail, per point): Johansen's Advanced Terrain Erosion Filter / Phacelle noise** (2026 blog). Stateless, GPU-friendly, takes a per-point direction (the hook for macro flow). Simple Phacelle has discontinuities; the "sampled" version is continuous at 16× cost.
  - Rust port: `bevy_erosion_filter` 0.3 (MIT), with a pure-Rust `cpu` module when default features are off. **Reimplement from the blog maths; Johansen's shader licence is unverified.**
  - Alternatives: Grenier et al. 2024 (gradient-aligned ravines; amplitude from distance to rivers); Dendry 2019 (branching distance field, too slow for per-tile use); IQ derivative-damped "eroded fBm" (IQ code MIT; a cheap fallback).
  - **Conditioning:** commit per macro cell the receiver direction, log drainage area, height, slope. Per point: d(p) = the interpolated downslope vector (interpolate vectors, not angles); gate on log A; amplitude from distance to macro rivers; a **channel clamp** so gullies end at the channel bed. Noise alone can't form exact confluences with macro rivers (the bevy_erosion_filter README says so).

## Quality metrics (for tests)

- **Hack's law** l = c·a^h: natural h 0.45–0.7, c 1–6. Schott 2023's generated terrains landed at h 0.54–0.58, c 1.0–1.4.
- **Concavity θ** (log S against log A on channels) ≈ 0.45–0.5.
- **Horton ratios:** R_B 3–5, R_L 2–3, R_A 3–6.
- **Depression invariant:** every land cell drains to the sea after filling (a hard assert).
- **Seams:** continuous across face edges and tiles.
- **Macro–detail agreement:** recompute drainage on final heights and compare receivers.

## Sources (licence / code)

| Source | Year | What | Licence | Code |
|---|---|---|---|---|
| Braun & Willett, Geomorphology | 2013 | O(n) implicit SPL | paper | via fastscapelib |
| fastscapelib | 2025 | LEM toolkit | **GPL-3** | github.com/fastscape-lem/fastscapelib |
| Landlab | 2014– | Python LEM, Priority-Flood | MIT | github.com/landlab/landlab |
| Barnes et al., Priority-Flood | 2014 | Depression filling | paper (richdem **GPL-3**) | — |
| Schott et al., TOG | 2023 | Interactive SPL, GPU drainage area | MIT | github.com/H-Schott/StreamPowerErosion |
| Tzathas et al., CGF | 2024 | Analytical SPL + multigrid | unverified | gitlab.inria.fr/landscapes/analytical-terrains |
| FastFlow (Jain et al.), CGF | 2024 | GPU flow/depression routing | unverified | gitlab.inria.fr/landscapes/fastflow |
| Barnes, EMS | 2017 | Tiled accumulation | paper | — |
| Braun, CAGEO | 2016 | GPU topo-sort accumulation | paper | — |
| Cordonnier et al., CGF / TVCG | 2016 / 2018 | Uplift + SPL on graphs / plates | paper | — |
| Johansen (blog) | 2026 | Phacelle, erosion filter | unverified (Shadertoy default CC BY-NC-SA) | — |
| bevy_erosion_filter | 2026 | Rust/WGSL port of Johansen | MIT | github.com/korbindeman/bevy_erosion_filter |
| Grenier et al., CGF | 2024 | Gradient-aligned patterns | paper | unverified |
| Gaillard et al. (Dendry), I3D | 2019 | Branching distance fields | paper | unverified |
| Clay John / Gavoronoise | 2018–23 | Eroded noise | **CC BY-NC-SA** | do not port |
| IQ derivative noise | — | Analytic-derivative fBm | MIT | Shadertoy snippet |
| Lague 2014; Horton 1945 | — | SPL review; stream-order laws | paper | — |
