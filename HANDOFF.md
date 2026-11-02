# HANDOFF

The current state, for the next session. History lives in `docs/session-logs/`. The newest log is the detail behind this file.

## State (2026-11-02, s01)

- **Branch:** `claude/peaceful-brahmagupta-rkl4y0`, an **orphan** history (v2 has no V1 ancestry). V1 lives on `main`. Matt may rename the branch to `v2`.
- **Docs:** `docs/ARCHITECTURE.md` (draft), `docs/research/2026-11-02-climate-tier-b.md`, `docs/lineage/` (lessons carried from V1 and tcotg).
- **Code:** M0 runs natively: eframe shell, offscreen reversed-Z globe, placeholder fBm planet, orbit/zoom, perf HUD. `rg-geo` (ellipsoid, cube-sphere, `CubeField`), `rg-noise`, `rg-render`, `rg-app`. CI workflow `.github/workflows/build.yml`.
- **In flight (s01, parallel Sonnet agents, uncommitted):** `crates/rg-tectonics` (plate snapshot → macro heights) and `crates/rg-climate` (tier A + Köppen). If a session dies mid-way, check `git status` for these dirs: they may be partial.
- **Web preview artifact:** https://claude.ai/artifact/1Tq1v7fz6zdbLTtUoLMPME. Rebuild with `scripts/build-web-artifact.sh <out>` and republish to the same URL (from a new session, pass `url`). GPU probe: https://claude.ai/artifact/64fD9SPXY3zpiZh8HjAURZ.
- **Verifying visuals in the container:** web screenshots don't work. Headless and headed Chromium both give white or blank WebGPU canvases, whatever the flags. **Use the native app under Xvfb + lavapipe:** `apt-get install -y mesa-vulkan-drivers libxkbcommon-x11-0`, `Xvfb :99 &`, `DISPLAY=:99 ./target/debug/realmgenesis &`, then `PIL.ImageGrab.grab(xdisplay=':99')`.

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

## Next

1. Integrate rg-tectonics + rg-climate into the app (the `rg-pipeline` layer), with layer switching.
2. 2D equirect view. World file + export.

## Open questions for Matt

- Paste one `mediasheet-for-obs` session log so the log style can match. It isn't in the public repo.
