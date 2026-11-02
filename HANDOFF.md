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
    - rivers + lakes in the layers; world file schema 3;
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

Scorecard: `cargo run -p rg-pipeline --release --example probe_realism -- 3` (scores the FINAL surface; `RG_RES`, `RG_EROSION_RES`, `RG_PROBE_SURFACE=macro`, detail `RG_*` knobs; worlds cached in `target/probe-cache/`). Eyes: `--example river_map -- <seed> <lat> <lon> <half_km> <out.png>` (`RG_MAP_FULL=1` for the detailed surface).

- **Done (s01):**
  - Erosion runs at the preset res (no 512 cap); 20 steps (same laws as 50).
  - Drained macro interpolant `sample_drained` (bilinear + diagonal-valley carve): bare-macro depressions 9.5% → 0.05–3%. A pure triangle split is fully drained but showed facets; rejected.
  - Rho8 (topological diagonals), initial roughness scaled by relief, lithology σ_lnK 0.3: straight channel steps 65% → 43–48%. This fixes Matt's "rivers too straight" combs.
  - Hills off (they caused 13–16% depressions). Gully steering by a 3 km central difference (no crease lines at macro cell edges).
  - Strict f32 remap; PF fill of upsampled macro; world schema 3 (height_res, macro_interp id, erosion seed).
- **Scorecard now (seed 3, mountain boxes):**
  - res 320: depressions 0.1–1.2%, slopes 6–9° mean (too gentle: coarse macro).
  - res 768: depressions 1.7–5%, slopes 14–17° mean, p99 55° (tail too heavy).
- **Next, in order:**
  1. **Dendry-style network detail** (Gaillard et al. 2019): hierarchical, point-evaluated sub-grid channel network; elevations from the SPL steady state with k_sn taken from the macro; hillslopes at threshold. Drained by construction (the min of drained surfaces is drained), replacing the gully filter's steep-or-drained trade-off and its "swirly" look. Write `docs/design/` first, then get a reviewer pass.
  2. **C1 macro:** square tonal blocks at the macro cell edges, since the shading of a C0 surface jumps there. A monotone bicubic plus the diagonal carve is one candidate; the network detail may cover it.
  3. **Bathymetry:** age–depth (Parsons & Sclater), shelf/slope/rise, trenches, abyssal hills, and Earth's ocean hypsometry as a check.
  4. Smooth river curves (map layer + corridor) instead of node-to-node straight segments.
  5. Reviewer leftovers: lithology tied to crust/terranes instead of isotropic noise; χ-analysis instead of the circular U/K concavity; validate the shipped (remapped) surface.

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
