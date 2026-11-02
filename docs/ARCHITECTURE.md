# RealmGenesis v2: Architecture

A fantasy planet generator: "Google Earth, but fantasy." A 3D globe and a flat 2D map; a raster view (satellite/data) and a vector view (political/cartographic); plate tectonics, a one-shot climate run, ocean currents; height detail generated as you zoom; sculpting.

Status: **draft, 2026-11-02, revised after the advisor-stand-in review (s01).** Decisions with their rationale and rejected alternatives also go in `HANDOFF.md`.

## 1. Stack

| Choice | Picked | Rejected, and why |
|---|---|---|
| Language | Rust (2024 edition) | Java: no web target, JVM to ship, GC pauses while streaming tiles. TS/Electron: heavy binary; the browser caps CPU/GPU use. |
| GPU | `wgpu` (Metal / Vulkan / DX12 natively, WebGPU on the web) | Raw Metal/Vulkan: one platform each. Bevy: heavy compile times, immature UI. |
| App shell + UI | `eframe` + `egui` (0.37) | Slint/Iced: kept as redesign candidates. Tauri: a webview defeats the point of native rendering. |
| Web build | `trunk` → wasm32 + WebGPU | — |

**The UI is replaceable.** Core crates never depend on egui or eframe. The UI reads state and sends `Command`s. A redesign rewrites `rg-app/src/ui/`, nothing else.

## 2. Targets and test loop

- **Minimum spec:** MacBook Air M1, 16 GB, 7-core GPU. Budgets are set against it. Secondary: RTX 3060.
- **Native:** macOS arm64, Windows x64, Linux. CI builds downloadable binaries.
- **Web preview:** the same app built with trunk and published as a private claude.ai artifact.
  - **Probed on the M1 in Vivaldi (2026-11-02):** WebGPU works inside the artifact frame (Metal 3, `shader-f16`, subgroups, 4 GiB storage buffers, max texture 16384). A compute shader did 1M cells × 64 iterations in 4 ms.
  - **There are no threads there:** no SharedArrayBuffer, because the frame isn't cross-origin isolated.
- **Consequence: simulations are GPU-first.** Compute shaders cover climate, erosion and detail generation. CPU code paths must run single-threaded on wasm. `rayon` is native-only, behind a cfg.

## 3. Crates (workspace; add each when its first code lands)

| Crate | Role | Deps |
|---|---|---|
| `rg-geo` | Ellipsoid (a, f; f = 0 is a sphere), geodetic ↔ ECEF, cube-sphere face/tile addressing, great-circle math | glam |
| `rg-gpu` | Device handle, compute and readback helpers, `Job` stepping, shared WGSL (plain `include_str!` concatenation; naga_oil only if that hurts) | wgpu |
| `rg-noise` | Deterministic integer-hashed noise sampled on the unit sphere, with a matching WGSL implementation | rg-geo |
| `rg-tectonics` | Plate snapshot → crust, boundaries, uplift (M2) | rg-geo, rg-noise |
| `rg-climate` | Tier A rules, tier B one-shot solver, Köppen (M3–M4) | rg-geo, wgpu |
| `rg-terrain` | Macro height composition, fine-detail amplification, erosion, hydrology (M2+) | rg-geo, rg-noise |
| `rg-world` | World document **types and serialization only**: params, committed fields, edit layers, the command/edit log (undo/redo) | rg-geo |
| `rg-pipeline` | Orchestration: drives the sims as `Job`s, writes results into `rg-world` | sims, rg-world, rg-gpu |
| `rg-render` | wgpu globe and flat renderers, tile quadtree, LOD, streaming. Owns offscreen colour + **reversed-Z Depth32Float** targets; takes `&Device`/`&Queue`, never sees eframe | rg-geo, rg-gpu, wgpu |
| `rg-app` | eframe app, egui UI, commands; native `main` and wasm entry point | everything |

Split for compile time: editing UI code rebuilds `rg-app` only. **The workspace `wgpu` version is pinned to the one egui-wgpu uses; upgrade them together.**

## 4. Geodesy and grids

- **One shape model:** an ellipsoid with semi-major axis `a` and flattening `f`. Default `f = 0`, a sphere, with a radius slider. "Make it oblate" sets `f` (later derived from rotation rate). **`f` is locked at generation;** changing it invalidates everything committed. Latitudes are **geodetic**. Heights are relative to the ellipsoid. Distances and bearings via geodesics (`geographiclib-rs` when needed).
- **Tile grid:** a cube-sphere quadtree with 6 faces. **Face warp: tangent (`tan(u·π/4)`), warp id 1, stored in the world file.** It inverts in closed form in shaders, and its area ratio is about 1.4. Rejected: equal-area warps that need iteration to invert.
- **Tile layout:** 256 cells = **257² corner-aligned vertices** (2^n + 1) plus a 1-sample skirt. Each parent sample coincides with every second child sample, so stitching is crack-free and delta downsampling is exact.
  - Cell size at level L ≈ `(π/2 · R) / (256 · 2^L)`. For R = 6371 km: L0 ≈ 39 km, **L3 ≈ 4.9 km (committed macro fields)**, **L4 ≈ 2.4 km (macro edit level, sparse deltas)**, **L10 ≈ 38 m (local edit level)**, **L12 ≈ 9.5 m (zoom floor)**.
  - **Sizes:** L3 is 6 × 2048² ≈ 25M samples per field, quantized to u16/f16 and compressed in the file. L4 dense would be about 100M samples (≈ 400 MB f32 per field), too big for the file, for wasm memory, and for serial hydrology.
- **Simulation grids:**
  - Tectonics: a graph of about 10k–100k cells (as in V1), upsampled to L3. Hydrology on L2–L3. Macro terrain committed at L3.
  - Climate runs on a regular lat-lon 1° grid (360×180, 0.5° later), with polar filtering. Every reference model uses lat-lon, and stencils stay trivial.
  - Results are resampled to the cube-sphere for display.
- **Precision:** world positions are f64 on the CPU. Rendering is camera-relative: tile origins relative to the camera in f64, uploaded as f32 offsets. Without this, jitter appears above about 1 m/px.
- **Noise precision (the `rg-noise` contract):** WGSL has no f64. An f32 unit-sphere position has about 0.4 m spacing at Earth radius, which stair-steps high octaves. **High-octave noise takes its lattice coordinates as integers from (face, level, tile x/y) plus an f32 fraction within the tile.** Only low octaves may use the raw direction.
- **2D mode:** the same tiles reprojected. Equirectangular first; others (equal-area, Dymaxion) later.

## 5. Determinism contract

Learned in tcotg: locked geography has to stay byte-identical. **The world file lands at M2**, not later, so this is testable from the start: a header, a schema version, field chunks. Golden tests: CPU paths byte-exact, GPU paths checked against a CPU reference within a tolerance.

- **Stored, never re-derived:** everything **committed** (macro heights, crust/plates, climate fields, sculpt deltas, vector features). The world file holds data, not just a seed. GPU floats are not bit-reproducible across vendors, so nothing gated may depend on GPU recomputation.
- **Re-derived:** fine procedural detail between the stored levels and the zoom floor. It's a pure function of (seed, committed coarse state, **committed vector features**, tile address, `detail_generator_version`).
  - **Detail is constrained:** amplitude is clamped near sea level using the committed coast signed-distance field, so the coastline never moves at zoom. It is carved along committed rivers, so rivers sit in their valleys.
  - **Versioned:** old `detail_generator_version`s stay callable, so stored sculpt deltas keep their meaning when the generator improves. It may differ in the last bits between GPUs. That's acceptable, because it's visual detail and no gate checks it.
- Noise uses integer hashing (no platform `sin` hashes), so CPU and WGSL agree up to float rounding.

## 5b. Jobs, not functions

GPU readback is async on the web: `map_async` resolves only after yielding to the browser. Long dispatches trip the GPU watchdog and lose the device. **So every simulation is a stepped `Job` state machine ticked once per frame.** Dispatches are chunked (under about 100 ms each), with a progress UI and device-lost recovery. Never `fn run_climate() -> Fields`.

## 6. Pipeline (user flow)

1. **Generate (reroll loop, seconds):** seed + params → tectonic snapshot → macro heights → tier A climate preview. Hit Generate until it looks close.
2. **Shape (hybrid):** paint land/sea intent and macro-sculpt at L4. Intent biases the crust field and re-runs macro terrain (lessons: `lineage/tcotg-v2-findings.md`).
3. **Commit:** tier B climate → hydrology/erosion → biomes → procedural borders and settlements. Committed fields are frozen into the world file. **A macro edit after commit marks climate stale**, and the UI offers a re-commit. Local (L10) sculpts never touch climate.
4. **Explore and edit:** stream tiles; local sculpt; edit borders, labels and settlements (vector layer).

## 7. Climate (see `research/2026-11-02-climate-tier-b.md`)

- **Tier A** (milliseconds): rules → monthly T/P → Köppen. Drives the reroll preview.
- **Tier B** (one-shot, target < 60 s on the M1): insolation → seasonal energy balance → pressure belts and winds → Sverdrup/Stommel gyres → moisture advection with orographic rain-out → sea ice → iterate 2–3× → Köppen. Steps 2, 4 and 5 run as compute shaders.
- **Display-time downscaling:** 1° misses ranges narrower than about 110 km. Temperature is downscaled via lapse rate from fine heights; precipitation via an orographic term (L3 slope · wind).

## 8. Layers and views

- **Raster:** satellite (biome albedo + hillshade + snow/ice + water), height, bathymetry, temperature, precipitation, wind, currents, Köppen, plates.
- **Vector:** coastlines (contoured from the committed macro heights), rivers, borders, settlements, roads, labels. Generated procedurally and editable. Rendered as screen-space-constant strokes in 2D and draped on the globe in 3D.

## 9. Milestones

| # | Deliverable |
|---|---|
| M0 | Workspace, eframe window (native + web), cube-sphere globe with placeholder fBm heights, orbit camera, egui panel (seed, radius, regenerate), **frame-time HUD**, CI binaries, artifact preview |
| M1 | `rg-gpu` + `Job` stepping; tile quadtree + LOD streaming with GPU-generated detail (integer-lattice noise); camera-relative precision; automated budget screenshots |
| M2 | Tectonic snapshot → macro heights (L3); fast reroll; **versioned world file + golden tests** |
| M3 | Tier A climate + Köppen; raster layers (height, climate, satellite v0); **flat 2D (equirect) mode** |
| M4 | Tier B climate on GPU |
| M5 | Hydrology, erosion, biomes; satellite shading |
| M6 | Sculpt: sparse L10 delta tiles with downsampled parents up to about L6, and sparse L4 macro deltas; intent painting; stale-climate flow |
| M7 | Vector layer: procedural borders/settlements/rivers/roads/labels, editable; 2D cartographic styles |
| M8 | Exports: 16-bit PNG/EXR height, equirect textures, GeoTIFF tags, glTF |

**Parked:** tectonic history (deep time), mobile, the custom-CRS WKT export, a UI redesign (Slint/Iced/themed egui).
