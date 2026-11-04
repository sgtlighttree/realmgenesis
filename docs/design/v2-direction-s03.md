# v2 direction after the s03 review (proposed, 2026-11-04)

Inputs: Matt's 12 answers (`docs/session-logs/2026-11-04-s03-project-review.md`), the s02 spike results, ROADMAP-V2 on `main` (raster-first). Status: **proposed**, awaiting a fresh-context review and Matt.

## What Matt wants (condensed)

- An all-in-one fantasy **worldbuilding** tool. Terrain first, then the "alive world" (Azgaar-style civs, borders, cities).
- **Naturalistic, photoreal-Earth terrain at modest resolution.** No infinite zoom. A finite raster world: ~5 km everywhere, sharp only where a region is baked.
- **Science first, but malleable:** import or edit coastlines, place them on the sphere, regenerate terrain that fits.
- **Native first.** Exports to Blender (DEM + albedo), Illustrator (vector), game engines (tiles).
- Reroll ≈ today's times (Ultra ~40 s); no overnight bakes.

## Diagnosis (from s01/s02, not new)

1. **The bones are fine:** tectonics → stream-power erosion + lakes → climate/Köppen; hypsometry and bathymetry calibrated to Earth. Matt: "the bones are there".
2. **Detail fails in two places:**
   - *Rendering*: CPU-built per-tile meshes with procedural detail (the drained detail network) give orbit artifacts, 56–240 ms/tile, UI stalls. Matt expected **raster displacement** all along.
   - *Generation*: the macro is D8-routed on an 8-connected cube-sphere grid (`rg-tectonics/src/grid.rs`), so its rivers run straight with 45° turns. The s02 lever round showed regional anchoring inherits that planform ("dumb upscale", "griddy", "mechanistic"). Trunks didn't reroute even at 4-cell anchoring.

## Direction

Raster-first data, raster-first rendering, un-gridded routing. Keep the Rust engine and the physics.

### R1: un-grid the macro (generation; small, judged by eye)
- Position-keyed node jitter on the global `Grid` (the solver only uses `dir/n/weight/adj`, verified s02), and/or MFD routing for the drainage-area pass. Fold in the parked **resolution-coherence** fix: position-keyed initial roughness/Rho8 at a fixed physical wavelength.
- Jittered nodes must be **resampled to the regular raster** before rendering or measuring (s02 artifact: unresampled jitter fakes isotropy).
- Gate: the 8-fold griddiness number G8 on the resampled raster, plus `river_map` PNGs judged by Matt, same seed and window as s02.
- Cost: **changes every world** (seed 1 Ultra won't reproduce); goldens updated deliberately.

### R2: raster renderer (display; biggest visible change)
- A committed **height raster** per cube face (the macro at the preset res; later baked regions as higher-res overlays) → GPU texture pyramid.
- Tiles become a fixed vertex grid displaced in the vertex shader; normals per pixel from the height texture (no facets, no CPU per-tile cost).
- Region rasters, when present, override the global one inside their footprint.
- The drained detail network and gully code leave the default path (code deleted once R2 lands and Matt signs off).
- **Zoom floor: 20 km camera altitude** to start; Matt adjusts.

### R3: coastline conditioning (the core interaction; design first, no code)
Today coastlines emerge from tectonics. Matt wants the generator to fit a land mask he imports (vector or raster, aligned/scaled on the sphere) or edits. Two broad options, not chosen yet:
- **(a) Mask constrains the tectonic seeding:** continental crust/plates grown to match the mask, so mountains, margins and bathymetry stay tectonically plausible.
- **(b) Mask imposed after tectonics:** uplift and sea level adjusted (or crust type overridden) to hit the mask, with erosion/climate run as usual. Simpler, less physical near imposed coasts.
Both then rerun erosion → climate. Editing = change the mask → regenerate (same seed).

### R4: regional bake ("Gaea mode")
- In-house (reading of Q9; Matt to confirm). Window solve on a jittered regional grid at ~150–300 m, **allowed to reroute** below a few macro cells; the region then **owns its rivers/lakes** inside its footprint (macro overlay suppressed there, edge crossings kept), plus hillslope diffusion/thermal for the finest texture.
- Output is stored rasters (determinism: only the stored data must be stable; GPU/particle methods allowed later).
- s02 budget: 25–60 s and ~2.7 GB for a 600 km window at 2049².

### R5: exports
Blender (16-bit PNG/EXR height + albedo, cube faces or equirect), Illustrator (SVG coastlines, rivers, later borders), engines (tiled heightmaps). The 16-bit heightmap export already exists.

Then the alive world (research doc `docs/research/2026-11-02-alive-world-features.md`).

## Order and why

R1 → R2 → R3 design → R4 → R5.
- R1 first: it's where the "mechanistic" look lives, it's small, it reuses the spike's jitter code, and it decides what the master raster contains. Matt can judge it from PNGs before any renderer work.
- R2 next: the display Matt expected; it also removes the tile-cost and UI-stall problems for good.
- R3 before R4: coastline editing is Matt's stated core interaction; regional detail is a refinement.
- Alternative considered: R2 first (most visible). Rejected: rendering the current griddy macro through a better renderer still looks griddy.

## Rules that change (proposed to Matt)

- Native-first: drop "CPU fallbacks single-threaded on wasm32" and "GPU-first *because the web has no threads*" (GPU stays where it pays). `rayon` allowed in core crates. The web build and artifact preview loop become optional.
- Settled decision "L12 (≈9.5 m) zoom floor" → replaced by the R2 floor.
- Rust stays; disk is managed by pruning `target/`, not by changing language.

## Open questions

- Q9 reading: in-house regional detail, external tools as export consumers only?
- R1: jitter alone, MFD alone, or both? (Decide from R1 PNGs.)
- R3: option (a) or (b), and the input formats (SVG/GeoJSON/PNG mask) to support first.

## Fresh-context review (s03): go with changes (adopted)

Reviewer read the code; I adopted all of these. They replace the sections above where they conflict.

- **Correction: the macro is already Rho8, not plain D8** (`solver.rs:121` `rho8: true`, applied at `:785`). Rho8 randomises the choice but still among the 8 lattice neighbours, so the planform stays lattice-bound. Jitter may help less than assumed, so R1 starts as a **jitter-only experiment**.
- **MFD is not cheap.** `Network` (`solver.rs:742–856`) is a single-receiver tree that stream power, lakes, `flow_dir`, rivers and `accumulate` all use. MFD is parked.
- **Jitter goes in the erosion grid only** (`solver.rs:545`), never the shared `Grid::new`, because tectonics uses that and coasts would move. The jittered grid also needs:
  - recomputed `weight` (else ±30% area noise);
  - `min_edge`/`max_edge` taken after jitter;
  - the hash keyed on the integer `node_sample` (not on `dir`, which depends on libm);
  - a check of the 8 valence-3 cube corners.
- **Resampling covers every field `to_field` writes** (`solver.rs:681–710`): height, drainage, lake depth, uplift, `flow_dir`. The fine Priority-Flood at `:703` needs rethinking too.
- **The raster must bake `sample_drained`** (`rg-geo/src/field.rs:121`). Plain GPU bilinear loses valley connectivity and pits the macro. Normals need cross-face aprons.
- **Store heights as R16 textures with manual bilinear** (vertex-texture float32 filtering needs a wgpu feature).
- **The zoom floor depends on bake coverage, not a flat 20 km.**
  - A 5 km texel at 20 km altitude is about 300 px.
  - Floor rule: one texel spans at most ~8 px. That gives **~800 km over unbaked land** and **~25 km over 150 m baked regions** (45° FOV, 1080p).
  - Alternative: synthetic detail on unbaked land, to be decided by eye.
- **Rivers need their own layer** once the detail network goes, since nothing else draws rivers under 5 km. A vector/overlay river layer also feeds the SVG export.
- **Schema:** region rasters need a new `rg-world` chunk (footprint + transform), and imported masks must be stored too (ARCHITECTURE §5).
- **Goldens:** there are no tectonic or erosion goldens, only the `terrain_realism` bounds and the lattice goldens. Re-baseline these when jitter lands.
- **Memory:** a native solve at 2048 peaks at ~4 GB and a region bake at ~2.7 GB. Never run both at once. The six R16 faces at 2049² are ~50 MB.

### Revised order

1. **R2: raster renderer on today's macro**, with drained interpolation baked in and a river overlay. Reason: it doesn't depend on jitter, and a 3D view settles whether the macro alone already reads as "decent terrain" (Matt, s02).
2. **R3: coastline-conditioning design** (doc only), before any R1 code. Option (b) feeds `tectonic_h` (`solver.rs:~552`), which decides where jitter can live.
3. **R1: jitter-only experiment** in the erosion grid, with full-field resampling. Judged by Matt from PNGs and the 3D view.
4. **R4: regional bake.** It reuses R2's override path.
5. **R5: exports.**
