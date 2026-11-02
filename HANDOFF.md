# HANDOFF

The current state, for the next session. History lives in `docs/session-logs/`. The newest log is the detail behind this file.

## State (2026-11-02, s01: MVP reached)

- **Branch:** `claude/peaceful-brahmagupta-rkl4y0`, an orphan history (v2 has no V1 ancestry). V1 lives on `main`.
- **What works** (verified by native Xvfb screenshots, plus all workspace tests and clean clippy):
  - Generate: tectonic snapshot (`rg-tectonics`, about 0.8 s at res 192) → tier A climate (`rg-climate`, about 0.2 s) → layers baked per pixel with climate downscaling (`rg-pipeline`).
  - Layers: satellite, elevation, plates, Köppen, temperature (year/Jan/Jul), precipitation. Globe + 2D equirect map, with a hover readout (lat/lon, height, °C, mm, Köppen).
  - Save/load `~/RealmGenesis/world-<seed>.rgworld` (`rg-world` container; byte-identical round trip; golden determinism test). 16-bit heightmap export + JSON sidecar.
  - Web preview artifact: https://claude.ai/artifact/1Tq1v7fz6zdbLTtUoLMPME. **Confirmed working by Matt on his Android phone** (claude.ai app, 2026-11-02): WebGPU in the artifact frame works on mobile too. GPU probe: https://claude.ai/artifact/64fD9SPXY3zpiZh8HjAURZ.
  - **M1a tile streaming on the globe** (`docs/design/tile-streaming.md`): zoom from orbit to the ground, Google Earth–style camera (drag pans, right-drag/Shift-drag/two fingers tilt and turn), worker-thread tile generation natively, frame-budgeted on the web. The web preview is not yet tested with tiles.
  - Quality presets Draft 192 / Normal 320 (default) / High 448 / XHigh 640 / Ultra 1024, plus Custom in km per cell (capped at 1024). Ultra ≈ 19 s tectonics, 1.04 GB peak (4-core container).
  - **Realistic terrain** (`docs/design/terrain-model.md`), checked against Earth:
    - tectonics → stream-power erosion (Braun–Willett + Priority-Flood; Hack h 0.54, concavity verified analytically) → hypsometric calibration (mean land 841 m vs Earth 840) → climate;
    - streamed detail = Johansen erosion filter (3D, MPL-2.0 port) on a drained macro interpolant (hills off since s01, see the current goal below);
    - rivers + lakes in the layers; world file schema 4;
    - CDLOD geomorph on the globe tiles.
  - Polar ice/snow edges are noise-perturbed (no longer latitude circles). The HUD shows fps vs cpu ms, tiles and altitude.
  - Controls sidebar toggles (`«`/`»` button at the top-left of the viewport, or **H**). It starts hidden when the window is under 700 px wide.
  - CI: `.github/workflows/build.yml` (test + clippy, macOS arm64 / Windows / Linux binaries, web bundle). **All green on 496a464** (run 41410046937). The toolchain is pinned to 1.97.0 in `rust-toolchain.toml`, because CI stable (1.99) added clippy lints that broke `-D warnings`.
- **Rebuild the preview:** `scripts/build-web-artifact.sh <out>`, then republish to the same URL (from a new session, pass `url`).
- **Visual checks in the container:** web screenshots don't work (WebGPU canvases capture blank). Use the native app under Xvfb + lavapipe:
  1. `apt-get install -y mesa-vulkan-drivers libxkbcommon-x11-0 xdotool`
  2. `Xvfb :99 &`
  3. `DISPLAY=:99 XDG_RUNTIME_DIR=/tmp RG_SEED=3 RG_LAYER=4 RG_VIEW=map ./target/debug/realmgenesis &`
  4. `PIL.ImageGrab.grab(xdisplay=':99')`

  xdotool keypresses don't reach the window (no window manager), so use the `RG_*` env vars.

## Current goal (Matt, s01): "very realistic, mathematically and scientifically accurate terrain", **plus bathymetry**

**Tools:**
- Scorecard: `cargo run -p rg-pipeline --release --example probe_realism -- 3`. `RG_PROBE_SURFACE=net` is the network via sub-patches; `=gully` the old detail; `=macro` the bare macro. It reports depressions, coherence/aspect (orientation), and 90 m slopes. Worlds are cached in `target/probe-cache/`, keyed by binary build time.
- Eyes: `--example river_map -- <seed> <lat> <lon> <half_km> <out.png>`, with `RG_MAP_FULL=1` for the final surface. Also `--example probe_bathymetry -- <seed>` and `--example bench_tiles`. `RG_NET_PROF=1` prints patch-build profiling.
- **Always check renders, not just numbers** (s01 lesson: the scorecard said Earth-like while the hillshade showed fishbones and mesas).

**State (all pushed):**
- **Macro:** erosion at the preset res; 20 steps; Rho8 with topological diagonals; relief-scaled roughness; lithology σ 0.3; k_sn = U/K from the committed uplift (schema 4).
- **Bathymetry:** Parsons age distribution, GlobSed sediment, Hoggard dynamic topography, Harris plateaus, narrow trenches. Mean ocean depth 3652–3687 m vs Earth 3688; bands within a few points.
- **Detail = drained detail network** (`rg-terrain/src/detail_net.rs`; design `docs/design/drained-detail-network.md` v1.1 plus the log):
  - min of valley wedges, drained by construction;
  - junction-up stream-power channels, incised below the coarser surface;
  - children attach to the lowest coarser wedge;
  - sibling chaining by hillslope rank;
  - per-channel junction angles, thinning, staggered LOD with an ancestor-consistent fade;
  - meandering 3-piece child channels; Catmull-Rom level-0 trunks;
  - flat hydraulic-geometry floors, band-limited;
  - exact block pruning and hierarchical candidate lists;
  - caches budgeted by bytes.
- **Scorecard** (seed 3, res 320): depressions **0.00%** in every box; mountains 15–23° mean, p99 31–35°; lowland 2.1°; coherence 0.28–0.49 (old gullies 0.48–0.91).
- **River overlay:** coverage in the layer alpha, faded when zoomed in. This fixes Matt's "terrain near rivers doesn't generate".
- **Web preview v7** at the URL above, with all of the above.

**Climate → erosion (done):** precipitation-weighted drainage; climate runs before and after erosion; wettest third of land erodes 7.1× the driest (it was 0.98×).

**Scientific review (s01 log, "Whole-pipeline scientific review"):** next steps in its order:
1. Earth reference harness (SRTM/GEBCO crops) + final-surface metrics (slope–area, Hack, Horton/Tokunaga, spectrum, shelf width, continental slope, ridge axis).
2. Endorheic basins/lakes from a water balance (climate is now in erosion).
3. Crustal thickness + isostasy instead of the hypsometric remap (hypsometry as a check only).
4. Sediment: transport-limited deposition (Yuan et al. 2019) and deltas/fans.
5. Glacial (ELA buzzsaw, U-valleys); ridge–transform segmentation; abyssal-hill fabric.

**Also next:**
1. **Tile cost:** 85–240 ms per tile natively (gullies: 9–33). Ideas: per-level window reuse across tiles, cheaper attach for the finest level, and eventually the GPU port. The web is single-threaded, so it streams slowly.
2. **Looks:** the uniform flat floodplain band of big rivers; faint arcs (reach walls); verify the globe app natively (Xvfb), including river areas.
3. **Final-surface hypsometry** (mean ~100 m below macro): a strictly increasing remap g with g(0) = 0. Stored receivers instead of flow_dir snapping.
4. Retire the gully/hill code once the app is happy; update `docs/design/terrain-model.md` §3 to the network.
5. **"Alive world" layer:** research in `docs/research/2026-11-02-alive-world-features.md`. Top 5: habitability/population → names → burgs → cultures → states. Licence flags are in the doc.
6. Reviewer leftovers: lithology tied to crust/terranes; χ-analysis.

## Settled decisions (don't re-litigate)

- Rust + wgpu + eframe/egui. The UI is isolated behind commands, and a redesign is planned later.
- GPU-first simulations, because the web preview has no threads.
- Sphere by default (`f = 0`), with the ellipsoid in the same code path. Geodetic latitudes.
- Cube-sphere tile quadtree with the tangent warp; 257² corner-aligned vertices per tile. Committed macro at L3 (≈4.9 km). L4 is sparse macro edits, L10 local edits, L12 (≈9.5 m) the zoom floor.
- Sims are stepped `Job`s (async web readback). High-octave noise uses integer lattice coordinates from the tile address. The world file + golden tests land at M2.
- Climate: tier A for the reroll preview, tier B for the one-shot commit. Committed fields are stored, never re-derived.
- Tectonics: a snapshot for the MVP; history mode parked.
- Borders: procedural by default, editable.
- Desktop/laptop only. Minimum spec: M1 Air.

## MVP definition (Matt's /goal: "until we have MVP", s01)

1. Native + web app; web preview artifact; CI binaries.
2. Generate: seed + params → tectonic snapshot → macro heights, fast reroll.
3. Tier A climate → monthly T/P → Köppen.
4. Layers: satellite-ish (biome colour + relief shading), height, plates, temperature, precipitation, Köppen.
5. 3D globe + 2D equirectangular map.
6. Save/load a versioned world file (native), and a 16-bit heightmap PNG export.

Post-MVP: LOD tile streaming, tier B climate, sculpting, the vector layer, erosion/rivers.

## Next (post-MVP, roughly by value)

1. Matt checks the web preview and the CI macOS binary on the M1; fix whatever he hits.
2. **Known quality gaps:**
   - Climate edges are blocky (1° grid; better downscaling or 0.5°).
   - Tier A weaknesses, listed in the rg-climate agent report (s01 log): dry Brazil-analogue interiors, weak monsoon, cold mid-latitude west coasts in winter.
   - Satellite palette tuning.
   - Tectonics: smooth Voronoi-ish plate borders, no fracture zones, land p99 about 6200 m.
3. **Terrain follow-ups:** detail isn't routed against the macro network between 5 and 30 km (closed depressions possible); peaks reach about 9–10 km with detail; tile CPU cost is about 40 ms deep, so the web build streams slowly (fix = GPU port).
3. **M1b:** `rg-gpu` compute tile generation, with a 3D integer-lattice noise port (see the review notes in the design doc); `Job` stepping for world generation (it still blocks with a busy overlay; Ultra ≈ 20 s).
   **M1c:** detail albedo close up (the global layer is ~10–20 km/px, so it's blurry); LOD normal seams; tile-height picking; a clearance check ahead of the camera; 2D map tiles.
4. Tier B climate on the GPU (see the research doc).
5. Hydrology/erosion/rivers → the vector layer (coastlines, rivers, borders, labels) → sculpting.
6. A file dialog (`rfd`) instead of the fixed folder; saving on the web.

## Open questions for Matt

- Paste one `mediasheet-for-obs` session log so the log style can match. It isn't in the public repo.
