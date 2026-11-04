# s02: background generation with a live progress console, and macro-only terrain

Status: **accepted with changes** (2026-11-04, s02; fresh-context review: "go, with changes", see §5). Matt chose "A now, B next" (see the s02 log): ship a macro-only look and a bigger monolithic basemap now, then design the regional "Gaea mode" (ROADMAP-V2's regional DEMs). He also asked for "an actively updating progress counter (not just a bar) so I know Something Is Being Done And The App Isn't Frozen. A bit like the System Console in V1."

## 1. Problem

- `RealmApp::run_deferred` calls `rg_pipeline::generate` **on the UI thread**. The window freezes for the whole run: ~3–6 s at Normal, ~17 s at 1024, > 1 min at 2048 (estimates). Nothing can update, so a live counter is impossible as things stand.
- The streamed detail network produces structural artifacts at every band (s02 log). The macro surface looks natural.
- The Custom detail cap is 2048 native (`ui/mod.rs: DETAIL_MAX_RES`).

## 2. Proposal

### 2.1 Generation off the UI thread (native), with progress

- New `rg_pipeline::generate_with_progress(params, on: &dyn Fn(Progress))`; `generate` becomes a wrapper with a no-op. `Progress { stage: Stage, step: u32, steps: u32 }`, where `Stage` is Tectonics, ClimatePre, Erosion (step k/N, from the solver loop), ClimateFinal.
- Erosion step progress: `rg_erosion::erode_climate_with_progress(.., on_step: &dyn Fn(u32, u32))` plus an internal `solve_full_with_progress`. The existing public functions keep their signatures and pass a no-op, so no test changes.
- Core crates stay UI-free: a plain callback, no egui types.
- The app spawns `std::thread` for Regenerate / RandomizeSeed. Progress and the finished `World` come back over an `mpsc` channel that the UI polls each frame. While generating:
  - the old world stays on screen and interactive (camera, layers);
  - Generate, Random, Load and Save are disabled, so only one generation runs at a time;
  - the UI repaints at ~10 Hz (`request_repaint_after(100 ms)`), not 60.
- When the `World` arrives: the same `after_world_changed` as today (equirect bake, tiles restart) on the UI thread. The bake costs ~35–40 ms, which is acceptable.
- **wasm32: unchanged** (no threads in the preview frame). It keeps the synchronous path with the busy overlay.
- RerunClimate (~80 ms) and LoadWorld stay synchronous.
- **No cancel in this step.** The thread can't be killed safely; a cancel flag checked between steps is a follow-up if Matt wants it.

### 2.2 Status-bar console (V1 "System Console" analogue)

- Status bar, left: a live line while work runs, e.g. `⠹ Generating · erosion step 7/20 · 4.3 s`. The spinner frame and elapsed time come from the clock each repaint, so they visibly tick at 10 Hz even inside a long single step.
- While tiles stream: `tiles · loading 12 · 220 drawn` (or `stopped` when paused).
- Idle: the existing world summary.
- A `Console ▴/▾` toggle on the status bar expands a panel of the last ~200 timestamped lines: stage starts and ends with ms, world summary, tile stop/resume, save/load/export results, errors. It's an in-memory ring buffer on `RealmApp`, also mirrored to stderr natively.

### 2.3 Macro-only terrain ("Streamed detail: off")

- A display setting, **not** stored in the world file (no schema change): `TerrainSampler` gets a detail mode {Network (today), Off}. With Off, `levels()` returns the drained macro height (`macro_height`) for `h`, `h_parent` and `h_child`. No patches, no gully path, no river corridors (rivers still show through the layer overlay).
- Toggling rebuilds the sampler and restarts tile streaming (`TileFeeder::set_world`).
- **Default: Off** for now. It's the natural look Matt approved, and it ends the per-frame heavy tile cost (a macro tile is a few ms). Network stays one click away.
- Tiles still refine for lighting (`max_cell_px`), but they're cheap.

### 2.4 Bigger monolith

- Raise the native Custom cap only after **measuring** the real peak RSS on the M1 at 1024 and 2048 (`estimate_peak_gb` is a 4-core-container fit: 1.04 GB at 1024, scaled by res²). Expected ceiling on 16 GB with other apps open: ~2048–2560. At 3072 the estimate is ~9.4 GB, too much.
- Keep the amber warning above 2 GB; add red above 6 GB with an "are you sure" confirmation.
- The real "deeper than the monolith" path is B (regional windows), not a bigger global grid.

## 3. Rejected / deferred

- **Cooperative `Job` stepping of generation (one stage per frame)** for native too: it would serve wasm as well, but the erosion step itself is seconds long at high res, so the UI would still freeze per step. Native threads now; wasm keeps the overlay.
- **Fixing the detail network now** (exaggerate only coarse relief, shallower fine levels): untested, and Matt's confidence is low. The code stays; it's off by default.
- **Raising the cap blind:** a wrong estimate means swap on a 19 GB-free disk, or an OOM kill.

## 4. Risks

- The `World` crosses threads: it must be `Send` (fields are `Vec`s and `CubeField`s; check at compile time).
- Two worlds in memory briefly (old displayed, new arriving): peak RAM roughly doubles at the swap. At 2048 that matters; the memory measurement should include it.
- `after_world_changed` stays on the UI thread: the equirect bake at 4096-wide is the main hitch (~40 ms measured at Normal; more at high res).

## 5. Review (fresh-context agent, 2026-11-04): go, with changes. Adopted:

1. **Macro-only seam:** one `macro_only` flag on `TerrainSampler`, checked at the top of `network_patch` (tiles call it directly, `tile.rs:87`), `levels` and `levels_with_slope`. Gating `levels()` alone would miss tiles, and "no network" alone falls into the old gully + river-corridor path. In macro-only mode the pipeline skips building the network (and its field clones).
2. **The core default stays Network:** `rg_pipeline::terrain_sampler` is unchanged, so tests and probes see today's surface; the app opts into Off through a new constructor.
3. **One generation at a time is enforced in `apply`**, not just by greyed-out buttons (the G/R shortcuts bypass buttons).
4. **Climate controls are disabled while generating** (otherwise RerunClimate hits the old world while the new one, with the old planet params, is on its way).
5. **CPU:** tile submission pauses while generating; the generation thread runs at utility QoS on macOS. Threads it spawns (tectonics/bake scoped pools) inherit the class (macOS QoS propagation; to verify).
6. **`World: Send`:** a compile-time assertion.
7. Progress callback `&dyn Fn` is called only from the orchestrating thread (the solver step loop); nothing inside `thread::scope` calls it.
8. wasm keeps the synchronous `busy` path; `std::thread::spawn` is behind `cfg(not(wasm32))`.
- **Faceting risk (to check, not fix yet):** macro-only tiles are piecewise-linear (`sample_drained`) with central-difference normals, so up close they may look low-poly with creases along the macro grid. Matt checks up close before Off stays the default; smoother normals in Off would be a `tile.rs` follow-up.
- **Console kept, but minimal:** status line + a collapsible panel of the last ~100 lines (Matt asked for the V1 console explicitly; the reviewer thought a 200-line panel was too much for a first step).
- **Unmetered phases:** tectonics and the post-solve part of erosion report no sub-steps; the clock-driven spinner and elapsed time cover them.
- Overstated: the bake is capped at 4096 wide (`layer_size`), so it doesn't grow past res 1024.
