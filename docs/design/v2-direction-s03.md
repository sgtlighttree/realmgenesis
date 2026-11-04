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
