# Regional detail ("Gaea mode"): design

Status: **proposed, M0 only** (2026-11-04, s02). Matt chose "A now, B next" (s02 log). This is B: ROADMAP-V2's regional detail DEMs (`ROADMAP-V2.md` on `main`: a world master raster plus regional DEMs at 100–200 m/px over ~500 km windows, cut into tiles).

## Goal

Pick a window on the globe (~300–600 km) and generate a high-resolution heightmap there (target 150–300 m/px) by a **coarse-to-fine erosion cascade** seeded and bounded by the global model. That's the workflow Gaea, World Machine, Houdini and Terragen all document (research digest below). It replaces the streamed detail network inside the window. Outside, the smooth macro surface stays.

## Order of work (advisor, s02): prove the look first

- **M0, headless spike (this doc's scope).** An example in `crates/rg-pipeline/examples/` (modelled on `river_map`) takes a cached world and a window, runs the cascade, and writes a hillshade PNG plus `examples/common/metrics.rs` stats, wall time and peak RSS (`/usr/bin/time -l`). **Matt judges the PNG before M1.** The s01 lesson is that the scorecard said Earth-like while the hillshade didn't.
- M1, storage: the region as optional schema-5 chunks (u16 heights with offset/scale, params in the header), plus a small golden test. **A fresh-context review is required here** (schema change).
- M2, rendering: `TerrainSampler` overrides `macro_height` inside regions. The region needs a mip pyramid so coarse tiles don't alias, and blends toward the macro across an edge band. Stored as the tile design's planned **data tiles** (256², sampled by render tiles at level L−2; `docs/design/tile-streaming.md`).
- M3, UI: pick a window, run it as a background job with progress and Esc cancel (the `worldgen.rs` pattern).

## M0 physics

**Reuse `rg-erosion` on a regional grid.** Verified in s02: `solve_full_with_progress`, `priority_flood_fill_to`, `lakes.rs`, `Network::route` and `diagonal_mask` use only `grid.dir`, `grid.n`, `grid.weight` and the neighbour lists. A regional raster can fill a `Grid` (8-neighbour CSR over a window of one cube face's (u, v); `res`/`sample_node` unused) without changing the solver.

**Window:** a rectangle of one cube face's (u, v), aligned to the tile quadtree (no windows across face edges in v1). The default is 2048² (about 4 M nodes).
- 4096² (16.7 M nodes) is **off-limits until M0's peak RSS is measured**. The global solve at ~25 M nodes peaked at 4.6 GB and pushed 4.4 GB into swap (s02).

**Cascade:** the macro is sampled onto the coarsest level (~1–2 km spacing). Each level then:
1. Upsample the previous level ×2 with the **drained** interpolant (`sample_drained` semantics), not the display bicubic. A non-drained start makes thousands of pits.
2. Add band-limited noise at the new octave, **keyed to position** (integer lattice coordinates of the global cube-face grid at that level, not node indices). Its amplitude is scaled by local relief.
3. Erode N implicit stream-power steps.

Levels double until the target. Steps per level go down as the resolution goes up (cost).

**Laws:**
- **Stream power with the macro's own constants:** `K_eff = K / height_scale` and the same m, n, so k_sn matches the committed surface.
- **Hillslopes: threshold slope (talus) at a repose angle of ~33°**, applied each step (Burbank et al. 1996 style; World Machine's thermal erosion is the same idea). Without it, stream power at 150 m gives knife ridges (advisor), which is the most likely difference between a Gaea-like and a network-like look. All parameters are in metres or degrees, so the result is resolution-independent (Houdini and Terragen document this).
- **No deposition in v1** (Yuan et al. 2019 is ROADMAP terrain #2).
- **No lakes in v1:** Priority-Flood fill, as in the no-lakes path. The macro's lakes stay as water bodies.

**Boundaries:**
- The edge ring is held at macro heights and acts as base level (outlets).
- **Inflow:** where the macro's `flow_dir` crosses into the window, the macro drainage area is injected at that edge node, as an extra runoff weight `1 + A_in / cell_area`. No solver change is needed. Without it, rivers entering the window starve and their valleys vanish.
- Ocean nodes are fixed from the macro, and land is clamped above 0 (the coast contract).

**Gates printed by M0:**
- the edge step (max |region − macro| on the border ring);
- the drift of the window's mean elevation against the macro;
- depressions, slopes, coherence and θ (shared metrics);
- wall time and peak RSS.

**Determinism:** grid methods only. Parallel droplet/particle erosion isn't bit-reproducible, and the region becomes committed data (ARCHITECTURE §5).

## Known limits to tell Matt

- **Colour inside the window stays blurry** (the global layer bake is ~10–20 km/px) until a regional albedo/biome pass. That's v1.1 and overlaps the M1c "detail albedo close up" item.
- Seams with the outside are blended, not physically continuous.

## Related finding: why resolution changes the terrain (Matt, s02)

"Generating at different resolutions also result[s] in very different terrains, even with the same coastlines and biomes."
- Coastlines and biomes come from position-keyed tectonics noise and the 1° climate.
- But the erosion's **initial roughness** (`solver.rs:561`, `unit_hash(seed, i)`) and **Rho8 diagonal randomisation** (`solver.rs:785`) are hashed by **node index**. A different resolution means a different random field, so different valleys capture different catchments.
- Stream-power erosion is also intrinsically resolution-dependent (drainage-area thresholds), as every tool's docs warn.
- Keying the roughness to position at a fixed physical wavelength would make major drainage much more consistent across presets. It changes every generated world, and may need a golden-test update, so it's **Matt's call** (a separate change from B). The cascade's new noise is position-keyed from the start.

## Research digest (Haiku + Exa, s02; licences checked by the agent from LICENSE files; re-verify before porting anything)

- **Tools:**
  - Gaea regions: build a window at higher resolution, where "almost all of the low-resolution detail is replaced" by erosion (docs.gaea.app, Regions).
  - World Machine and Gaea tiling: an extra border plus blending for simulation nodes.
  - Houdini: erode at low resolution, resample up, erode finer, with feature sizes in metres (sidefx.com heightfield docs).
  - Terragen: erode, resize, erode again (Planetside blog).
- **Algorithms:**
  - Braun & Willett 2013 (O(n) implicit stream power; we have it).
  - Schott et al. 2023 (GPU stream power, interactive up to 4096²; **MIT** code at H-Schott/StreamPowerErosion), the reference for a later wgpu port.
  - Yuan et al. 2019 and Davy & Lague 2009 (deposition, later).
  - Roering et al. nonlinear diffusion (hillslope alternative to talus).
  - Mei et al. 2007 pipe model and Beyer 2015 droplets (texture only; not deterministic in parallel).
- **Portable with attribution:** SebLague/Hydraulic-Erosion (MIT), H-Schott/StreamPowerErosion (MIT), henrikglass/erodr (MIT), symbios-ground (Rust, MIT), fabio-t/terrain-generator (Rust, Apache-2.0), mapgen4 (Apache-2.0).
- **Ideas only:** fastscapelib (GPL-3.0) and Badlands (GPL-3.0). Landlab, dandrino/terrain-erosion-3-ways, weigert's repos and rj00a/heightmap-erosion have unverified or missing licences.
- **Nested boundaries:** fixed-elevation open edges as base level (Terrainbento; Atkinson 2017). No source states an inflow-area rule, so ours above is a design choice.
