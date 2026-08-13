# D7 polish — plate-shape realism without a substrate rebuild

**Date:** 2026-08-13 · **Milestone:** ROADMAP D7 (part 3 decision) · **Status:** design, approved in substance

## Decision that frames this: NO-GO on the Cortial rebuild

D7 part 3 was framed as the "Cortial boundary-curve rebuild" — plates as a topological
graph of boundary curves + terranes instead of seed-grown regions. Both the `advisor`
tool and the `fable-advisor` subagent independently ruled **NO-GO**, for the same reason:

- The rebuild replaces the plate-assignment substrate outright — killing
  `assignPlatesDijkstra`, `computeEdgeCosts`, `mergeSmallPlates`, `injectMicroplates`,
  and the **connectivity-by-construction 0-exclave invariant** earned in Session 9.
  Cortial's curve-graph has no equivalent guarantee; curve intersection + retriangulation
  on a sphere is weeks of geometric-robustness edge cases (degenerate junctions,
  self-intersection), re-proving the invariant from scratch.
- The commissioned research report's **own** engineering recommendation for a ~10k-cell
  browser budget (§4) is the heuristic tier — warped-distance Dijkstra region growth +
  distance-to-ridge GDH1 + shear microplates + abyssal noise — **which Sessions 9–10
  already shipped.** Cortial is listed as feasible but explicitly not recommended here.

**Empirical check (seed `realmgenesis`, current `main`, rendered 2026-08-13):** fine
boundaries are genuinely organic (Session 9's noisy edge costs work), but plate
*silhouettes* remain blobby-convex. That residual "Voronoi look" is a **seeding artifact**
(point-seed + Dijkstra → quasi-convex regions), not an edge-noise or substrate problem.

This spec closes the remaining gap with three additive, cheap items from the report's own
§4 tier — no substrate change.

## Goal

Make plates read as tectonically shaped, not Voronoi-shaped, while preserving:
- the 0-exclave / single-connected-region-per-plate invariant (verified over the display graph);
- run-to-run determinism;
- the existing param/erosion/climate pipeline downstream of `projectTectonicsToDisplay`.

All work is in `utils/tectonicsV3.ts` unless noted.

## Item 1 — Band / chain seeding *(dominant term)*

**Problem.** `assignPlatesDijkstra` (tectonicsV3.ts:296) seeds each plate at a **single**
macro cell (the one nearest its rotated Euler seed), then grows outward. Single-point
sources over a near-uniform cost field yield quasi-convex blobs regardless of boundary
noise — the Voronoi read.

**Change.** Seed each plate as a **connected chain** of macro cells, all pushed as dist-0
sources. Elongated seed → elongated plate. Because every chain cell is a graph-neighbor of
the previous one and all belong to the same plate, the plate stays a single connected
region — the 0-exclave invariant holds by construction.

**Mechanism** (in the seed-setup loop, tectonicsV3.ts:315–325):
1. Find nearest cell `best` to `rotatedSeeds[j]` as today.
2. Compute the plate's tangential velocity direction at `best` from its (rotated) Euler
   pole: `v = ω × r`, projected onto the local tangent plane, normalized. This is the
   chain axis — **derived from existing state, no RNG draw.**
3. Walk `K = 1 + round(plateElongation × 4)` cells: from the current tip, step to the
   unused neighbor whose direction from the tip best aligns (max dot) with the chain axis.
   Stop early if no forward neighbor exists.
4. Push every chain cell as `{cell, plate: j, dist: 0}`.

**Param.** New `plateElongation` (0–1, **default 0.4**). At 0, `K = 1` → identical to
today's single-source seeding, byte-identical. The default is non-zero on purpose: it
fixes the blobby look for all existing seeds (accepted re-baseline, below).

**Why velocity-aligned, not random:** physically, plates elongate along their spreading
direction; it needs zero new RNG (cleaner determinism); and the chain axis rotates
consistently with the seed each timestep. Random-axis (a `_elong_v3` stream) was the
alternative — rejected for the extra stream and no realism gain.

**Invariant gate:** the connected-components probe from Session 9 must still show every
macro plate = 1 component after this change.

## Item 3 — Seafloor age-band perturbation *(cosmetic, ~20 min)*

**Problem.** `computeSeafloorAge` sets `age = min(180, dist/spreadRate)` (tectonicsV3.ts:439).
Distance-to-ridge with a constant rate produces artificially clean, symmetric age bands —
the report's §4.2 named artifact.

**Change.** Perturb: `age = clamp(age × (1 + 0.1 × noise3D(macroPoints[i])), 0, 180)`,
where `noise3D` is a fresh simplex on side-stream `new RNG(seed + '_agenoise_v3')`.
Bathymetry (Height / seafloor) views only; does not touch land or plate layout.

**Risk:** none if the stream is fresh — never draw from `plateRng`/`microRng`.

## Item 2 — Transform-edge fracture *(conditional; attempt only after 1+3 are rendered)*

**Sequencing rule (advisor-endorsed):** build 1+3, render, judge. Band seeding is the
dominant term; if it alone kills the Voronoi read, item 2 is cheap polish. If item 2's
jaggedness turns out indistinguishable from `boundaryRoughness` at high values, **drop it**
— it is worth attempting only because approach (b) makes it ~an afternoon, not because the
deliverable needs it.

**Problem / why it's non-trivial.** A transform boundary is where relative plate velocity
is tangential to the boundary — classification requires plate **ownership**, but
`computeEdgeCosts` runs *before* `assignPlatesDijkstra`, so no edge can be known as a
transform at cost-computation time. (The existing `boundaryRoughness` noise is
ownership-agnostic, which is why it can live there.)

**Change — approach (b), ride the timestep loop.** `simulateTectonics` already re-runs
`assignPlatesDijkstra` every step and classifies boundaries right after. Feed step *k*'s
transform classification into step *k+1*'s cost field:

```
effectiveCost(edge) = staticCost(edge) × (isTransformLastStep(edge) ? fractureMul(edge) : 1)
```

- **Recompute as a SET each step, never accumulate.** If costs are multiplied again each
  step they compound and boundaries run away, and an edge that stops being a transform
  keeps a stale bump. Rebuilding from `staticCost` each step means edges cleanly revert.
- `fractureMul(edge)` from a seeded 1D noise keyed on edge index, side-stream
  `new RNG(seed + '_fracture_v3')`. Deterministic, bounded.
- Assignment stays region-growth → macro 0-exclave property holds; display projection
  inherits it with no cell swaps. (Display-space polyline offset was rejected — it is the
  peeling anti-pattern relocated, and risks display exclaves.)

**Accepted behavior:** the jagged pattern only fully settles if transforms are stable in
late steps. That's physically reasonable; do not force convergence.

## Determinism & re-baseline

- Three new side-streams: `_agenoise_v3`, `_fracture_v3` (and none for item 1 — it is
  RNG-free). **None touch `plateRng` / `microRng` draw order.**
- `plateElongation` default 0.4 changes plate shapes for **all** seeds → the
  `terrainSignature` baselines re-baseline. This is expected and has precedent (V3
  re-baseline, commit `7f572f5`). Re-scan lakes seeds if hydrology shifts (as Session 9/10 did).
- Reuse the MinHeap index/plate tie-break already in `assignPlatesDijkstra` — chain seeds
  must pop deterministically.

## New params

| Param | Range | Default | Notes |
|---|---|---|---|
| `plateElongation` | 0–1 | 0.4 | chain length K = 1 + round(e×4); 0 = today byte-identical |

`spreadRate`, `seafloorDetail`, `microplateIntensity` already exist (Session 10) and are
reused. Items 2–3 add **no** user-facing params: the age perturbation is fixed-amplitude
(0.1) and the transform `fractureMul` is a bounded noise multiplier — both always-on
textures, no sliders. Only `plateElongation` is new.

## Testing

- **paramLiveness:** add a `plateElongation` case (it must change `terrainSignature`).
  Items 2–3 ride existing terrain params; no new liveness cases required, but confirm the
  suite still passes after re-baseline.
- **0-exclave probe:** re-run the Session-9 connected-components check across
  `realmgenesis` / `route-test` / `abcxyz` after item 1 and after item 2.
- **Determinism:** same-seed-twice byte-identical height + plateId, after each item.
- **Gates:** typecheck 0, lint 0 errors / ≤30 warnings, full vitest, build. (Note the
  Session-10 parallel-load timeout flake in `paramLiveness` — real failures assert;
  timeouts are the flake.)

## Out of scope

- The Cortial rebuild (terranes, boundary-curve graph) — deliberate NO-GO, recorded above.
- Explicit triplet/cluster seeding for triple junctions — rejected: disjoint sources
  reintroduce exclaves; junctions emerge naturally where 3 plates meet. Bands (connected
  chains) are the safe elongation form.
