# HANDOFF

The current state, for the next session. History lives in `docs/session-logs/`. The newest log is the detail behind this file.

## State (2026-11-02, s01: MVP reached, pending Matt's browser check)

- **Branch:** `claude/peaceful-brahmagupta-rkl4y0`, an orphan history (v2 has no V1 ancestry). V1 lives on `main`.
- **What works** (verified by native Xvfb screenshots, plus all workspace tests and clean clippy):
  - Generate: tectonic snapshot (`rg-tectonics`, about 0.8 s at res 192) → tier A climate (`rg-climate`, about 0.2 s) → layers baked per pixel with climate downscaling (`rg-pipeline`).
  - Layers: satellite, elevation, plates, Köppen, temperature (year/Jan/Jul), precipitation. Globe + 2D equirect map, with a hover readout (lat/lon, height, °C, mm, Köppen).
  - Save/load `~/RealmGenesis/world-<seed>.rgworld` (`rg-world` container; byte-identical round trip; golden determinism test). 16-bit heightmap export + JSON sidecar.
  - Web preview artifact (v2 = commit 7015550): https://claude.ai/artifact/1Tq1v7fz6zdbLTtUoLMPME. **Not yet confirmed loading in a real browser.** GPU probe: https://claude.ai/artifact/64fD9SPXY3zpiZh8HjAURZ.
  - CI: `.github/workflows/build.yml` (test + clippy, macOS arm64 / Windows / Linux binaries, web bundle). **All green on 496a464** (run 41410046937). The toolchain is pinned to 1.97.0 in `rust-toolchain.toml`, because CI stable (1.99) added clippy lints that broke `-D warnings`.
- **Rebuild the preview:** `scripts/build-web-artifact.sh <out>`, then republish to the same URL (from a new session, pass `url`).
- **Visual checks in the container:** web screenshots don't work (WebGPU canvases capture blank). Use the native app under Xvfb + lavapipe:
  1. `apt-get install -y mesa-vulkan-drivers libxkbcommon-x11-0 xdotool`
  2. `Xvfb :99 &`
  3. `DISPLAY=:99 XDG_RUNTIME_DIR=/tmp RG_SEED=3 RG_LAYER=4 RG_VIEW=map ./target/debug/realmgenesis &`
  4. `PIL.ImageGrab.grab(xdisplay=':99')`

  xdotool keypresses don't reach the window (no window manager), so use the `RG_*` env vars.

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
3. M1 proper: `rg-gpu`, `Job` stepping (generation is currently blocking between frames, with a busy overlay), tile quadtree LOD streaming with integer-lattice detail noise.
4. Tier B climate on the GPU (see the research doc).
5. Hydrology/erosion/rivers → the vector layer (coastlines, rivers, borders, labels) → sculpting.
6. A file dialog (`rfd`) instead of the fixed folder; saving on the web.

## Open questions for Matt

- Paste one `mediasheet-for-obs` session log so the log style can match. It isn't in the public repo.
