# s09 (2026-11-14, local): rivers on imported heightmaps

Goal (Matt): improve river generation on imports. Test maps: `~/Downloads/Klerga V2 Maps HeightmapMono 2.png` (whole world) and `Klerga (100k) Heightmap 2023-07-22-17-13.png` (patch 20,10 at 0.5 km/px). Matt: s08b's drawing is "still off"; the drop bug is fixed as long as the layer is checked.

**Stopped mid-work by the usage limit. The code is a WIP prototype with `std::env` knobs in `rg-erosion` (must go before this is done).**

## Findings
- **Density is fine, shape is not** (`probe_rivers`, extended this session; V2 freedom 0 vs procedural seed 1): same km of river per land area, but lakes cover 2.6% of land vs 0.1%, 15% of river length is under lakes, and trunks run ruler-straight for 500–1,000 km. Freedom 1 (Guided) is no better (lakes 3.0%, same flats): the inverted uplift makes the import a steady state, so erosion reorganises nothing.
- **Cause of the straight runs: steepest descent (D8) on a smooth drawn slope.** It only knows the lattice's 8 directions; the import has no cell-scale relief to break them. Verified by `river_map` renders (hillshade + full network). Refuted on the way: "the map is terraced" (the drop to the receiver is p10 6 m, median 80 m per cell: real slopes).
- **Tried, in order (all on the final routing only, heights untouched):**
  1. fBm added to the routing surface: bends courses but digs its own basins (lakes 6.5% of land). Rejected.
  2. Wander as a *bias* among truly lower neighbours (`Network::wander`): no new pits, but still lattice runs, plus zigzags where the bias wants uphill.
  3. Noise-weighted ε on filled flats (`priority_flood_fill_eps` `cost`): harmless, little visible effect on V2 (few true flats).
  4. Guided solve with initial roughness on pinned land (`RG_PIN_ROUGH`): slightly more organic, still herringbone. Not the answer alone.
  5. **Stochastic receiver: one of the truly lower neighbours drawn with odds ∝ slope^2** (`RG_COND_POW`, with wander strength 0.5): the best so far in `river_map`; trunks irregular and dendritic, a few lattice runs and zigzags left. n=1 region, unjudged by Matt.
- **Shallow drawn pits (≤ 40 m, `RG_COND_PIT`) filled before lakes; deeper basins kept:** lakes 2.6% → 1.7% of land. The remaining blocky lakes are the deep drawn basins.
- **The app showed none of the routing change:** suspected cause (unverified): the app's rivers come from `refine.rs` (finer preset than canonical 320), which has its own routing and is not conditioned yet. `TectonicsParams::default().res` is 320, so check what the app actually runs at before trusting this.
- **Second fault, separate from routing:** zoomed in on a small continent (V1 at 2,000 km) rivers are fat blurry strokes: the river field is 4.9 km per texel and the shader never draws under 1.2 texels. Needs real line geometry or a finer field near the camera.

## Next
1. Find why the app's rivers ignore the conditioning (refine path?) and condition that routing too.
2. Judge stochastic routing in the app on both maps; tune power / wander; remove zigzags.
3. Replace the env knobs with `ErosionParams` fields; delete the `RG_COND_DEBUG` blocks and `RG_PIN_ROUGH`.
4. Gates: procedural seed 1 `world_hash` identical; imports' height hash identical (drainage/lakes/rivers change); clippy; tests.
5. Zoomed-in river drawing (item above).
6. Fresh-context review before finalising (project rule), then Matt's verdict.
