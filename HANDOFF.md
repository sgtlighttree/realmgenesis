# HANDOFF

The current state, for the next session. History lives in `docs/session-logs/`. The newest log is the detail behind this file.

## State (2026-11-02, s01)

- **Branch:** `claude/peaceful-brahmagupta-rkl4y0`, an **orphan** history (v2 has no V1 ancestry). V1 lives on `main`. Matt may rename the branch to `v2`.
- **Docs:** `docs/ARCHITECTURE.md` (draft), `docs/research/2026-11-02-climate-tier-b.md`, `docs/lineage/` (lessons carried from V1 and tcotg).
- **Code:** none yet. M0 scaffold in progress.
- **Preview artifact:** not published yet. GPU probe: https://claude.ai/artifact/64fD9SPXY3zpiZh8HjAURZ (WebGPU compute works in the frame; no threads).

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

## Next

1. M0: workspace + `rg-geo` + `rg-render` globe + `rg-app` (eframe, egui panel, orbit camera), building native and web.
2. A headless screenshot, then publish the trunk build as the preview artifact.
3. A CI workflow building macOS arm64 / Windows / Linux binaries and the web bundle.

## Open questions for Matt

- Paste one `mediasheet-for-obs` session log so the log style can match. It isn't in the public repo.
