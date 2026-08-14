# Tectonics — the V3 terrain model

`utils/tectonicsV3.ts`. V3 is the **only** terrain model (V2 and the
`V3_ENABLED` flag were removed at Session 11). Its founding idea: **crust and
plates are independent fields.** Crust (continental vs oceanic) is seeded from
its own noise; plates deform it but do not decide where land is. This decouples
the coastline from the plate silhouette — the defect that motivated the rebuild
(V2 flipped whole plates to land/ocean, so continents *were* plates).

The simulation runs at coarse **macro resolution** (`simulationResolution`,
default 10,000 cells) and projects onto the display cells in one pass at the
end. It never runs at display resolution — that would be far out of a browser's
budget at the 200k cap.

## `simulateTectonics` — the steps

1. **Seed crust field** — `seedCrustField` (`utils/crust.ts`): per macro-cell
   continental/oceanic type + thickness, from the `_crust` stream. Independent
   of plates.
2. **Generate plates** — `generatePlates`: plate seeds with Euler poles
   (rotation axis + rate); `plateJitter` perturbs seed placement. Small plates
   are merged (`mergeSmallPlates`) into neighbors for a power-law size spread.
3. **Build macro neighbor graph** — `buildMacroNeighborGraph` (top-k, then
   **symmetrized** — Dijkstra needs an undirected graph).
4. **Assign plates by Dijkstra region-growth** — `assignPlatesDijkstra` grows
   each plate outward from the macro cell nearest its (rotated) seed, over
   **noisy static edge costs** (`computeEdgeCosts`: `chord × noiseMul ×
   marginMul`). `noiseMul` from `boundaryRoughness` (sampled at a warp-displaced
   edge midpoint, keeping `warpStrength` live) roughens fronts; `marginMul` from
   `marginCoupling` attracts boundaries to crust-type transitions. Then
   `mergeSmallPlates` again.
5. **Timestep loop** (`numTimesteps`, default 20): rotate seeds by their Euler
   poles, re-assign via Dijkstra each step, classify each boundary by relative
   velocity and crust pair, and accumulate uplift:
   - continental collision → broad symmetric uplift (×60)
   - oceanic subduction under continent → arc uplift (×30) + trench (−20)
   - transform → modest shear (×5)
   - continental–continental divergence → rift lowering
   - all scaled by `tectonicStrength` × per-boundary Simplex `noiseFactor`
     (segmented belts). **Oceanic divergence is NOT lowered here** — its
     elevation is GDH1's job (a corrected mistake; see
     [ENGINEERING-NOTES.md](ENGINEERING-NOTES.md)).
6. **Inject microplates** — `injectMicroplates` runs *after* the small-plate
   merge (so the cutoff doesn't eat them): finds top-shear boundary cells,
   appends a `PlateState` per pick with a low speed, plants **one** seed, and
   lets Dijkstra grow it. Seed injection, not post-hoc peeling — this preserves
   connectivity by construction. `microplateIntensity` controls the count.
7. **Seafloor age → bathymetry** — `computeSeafloorAge`: multi-source Dijkstra
   distance-to-ridge over **oceanic cells only** (continents block), `age =
   dist / spreadRate` capped at 180 Ma. `gdh1Depth` (Stein & Stein 1992 GDH1)
   maps age → metres; `depthToBandHeight` maps that into the oceanic **height
   band** (ridge ≈ −0.5, old floor ≈ −0.85) — *not* metres into the raw field.
   Empty-ridge worlds (Pangea) fall back to isostasy.
8. **Compose height** — oceanic-with-valid-age → GDH1 band + abyssal-hill noise;
   else `saturatedIsostasy(thickness)`; plus `upliftAccum·0.5`;
   `applyPassiveMargin` softens shelf transitions; normalize to [0,1].

## Projection to display

`projectTectonicsToDisplay`: each display cell copies its nearest macro cell's
`crustType`/`crustThickness`/`upliftAccum`/`plateId`/height, then adds
full-resolution structural noise (FBM blended with ridged via `ridgeBlend`,
scaled by `roughness`). Over deep-ocean cells the structural noise is **damped**
(baked factor) so the GDH1 gradient survives, plus low-amplitude abyssal-hill
noise. Final normalize.

## Load-bearing invariants

- **0 exclaves by construction.** Dijkstra region-growth makes each macro plate
  exactly one connected component; a shared `claimed` set keeps per-plate seed
  sets disjoint. Verified across seeds. This is *the* property the Cortial
  rebuild would have destroyed (deliberate NO-GO —
  [ENGINEERING-NOTES.md](ENGINEERING-NOTES.md)).
- **Macro connectivity is guaranteed; fine-mesh display connectivity is not.**
  The nearest-macro downsample can pinch a thin plate into disconnected display
  strays. Pre-existing; `tests/plateConnectivity.test.ts` guards it only for
  seed `realmgenesis`.
- **De-blob comes from `plateJitter` + `boundaryRoughness`** (both default 1.5),
  not `plateElongation` (near-inert at the macro silhouette). Full rationale in
  [ENGINEERING-NOTES.md](ENGINEERING-NOTES.md).
- **Determinism:** side-streams `_macro_v3`/`_crust`/`_plates_v3` (+ `_micro_v3`,
  `_abyssal_v3`, `_agenoise_v3` inside); the MinHeap carries a tiny index/plate
  tie-break so heap ordering is deterministic. Same seed → byte-identical.
