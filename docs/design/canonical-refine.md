# Canonical macro solve + fine refinement (s04)

Status: approach approved in direction by Matt (s04: "Ultra basically upscales Normal; Normal gives the grand scale"). Advisor-reviewed. Log: `docs/session-logs/2026-11-04-s04-r3-design.md`.

## Problem

Erosion is resolution-dependent:
- stream power treats every cell as a channel;
- the number of tiny hollows held as lakes grows with resolution;
- a global calibration (mean land height) spreads both effects to every region.

Seed 1, Normal vs Ultra: height factor 0.748 vs 0.565; in Matt's box, p99 4745 vs 3800 m. A snow massif at Normal is gone at Ultra.

Physical floors (a channel-head area, filling small basins) made the two coherent, but changed Normal's look and removed lakes. Rejected.

## Design

Normal (320) is the **canonical** world. Every preset above it is Normal plus detail.

1. **Canonical stage = Normal's pipeline, unchanged.** Tectonics at `canon = min(320, res)`, climate pre-pass, erosion with lakes and calibration. For a res-320 world nothing changes: **Normal stays byte-identical** (a gate).
2. **Fine tectonics** at the target res, for the fine fields (plates, crust, boundary, bathymetry) and the sub-cell detail `D = T_fine − up(T_canon)`. Tectonics is already coherent (r 0.98), so `D` is detail, not new large-scale shape.
3. **Refine** (`rg-erosion/src/refine.rs`, one function reused later by R4's regional bakes):
   - Land start: `h = up(H_canon) + s·D` (`up` = `sample_smooth`, Catmull-Rom; `s` = the canonical height factor, so the detail is in committed units). Floor 1 m. Ocean = `T_fine ≤ 0`, using `T_fine` heights (fine bathymetry).
   - Lakes: **the lake model does not run on the fine grid** (it multiplies tiny lakes there). Canonical lake levels are upsampled and held as base levels. Every other fine pit is plain-filled (priority flood to the ocean and the lakes).
   - Incision: a few implicit stream-power steps with **U = 0** (it only cuts sub-cell valleys and can't build large relief), with runoff weighting and lithology as in the canonical solve, over a short `refine_myr`. **No recalibration.**
   - Outputs on the fine grid: height, drainage, flow direction, rivers, lake depth (upsampled level − fine bed). Uplift is upsampled from canonical. The height factor and stats are the canonical ones.
4. Climate runs on the final heights as today.

## Gates

1. **Normal byte-identical** before and after (seeds 1 and 3, 320).
2. **Ultra low-passed to Normal's cell matches Normal**, regionally: box p50/p99 within a few %, plus the rms in metres (the probe's 5×5 cell average); global r is not enough. The strength (`refine_myr`) is the largest that keeps this passing.
3. Refined height factor = canonical, by construction.
4. **A `cargo test` guard:** canonical 64 → refine 128. It asserts gate 2 and that the canonical output is unchanged, so this can't come back silently (Matt asked).
5. Look: Matt in the app, Normal vs Ultra.

## Scope

- **Draft (192)** stays its own coarse preview: not coherent with Normal, documented. (Alternative: compute at 320 and downsample. Ask Matt.)
- **Custom:** `canon = min(320, res)`. High, XHigh and Ultra all go through refine.
- **World file:** the header gains `canonical_res` and `refine_myr`/`refine_steps`. Older files load as before (no refine).
- **Cost (Ultra):** about two tectonics runs (320 + 1024), one Normal erosion, and a few fine fill/route/step passes. Likely faster than today's full Ultra erosion (144 s capped).
- Pitfall: an upsampled surface is smooth, so filling makes ε-flats that route in parallel straight lines (the s01/s03 grain). `D` (position-keyed tectonic detail) supplies the roughness that breaks this. Watch for it in renders.
