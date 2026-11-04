# HANDOFF

The current state, for the next session. History lives in `docs/session-logs/`. The newest log is the detail behind this file. **What to do next: `ROADMAP.md`.**

## s03 (2026-11-04): review in progress

- Matt answered the 12 questions: `docs/session-logs/2026-11-04-s03-project-review.md`.
- Proposed direction (**not yet approved**): `docs/design/v2-direction-s03.md`. In short: keep Rust; native-first; R1 un-grid macro routing (jitter/MFD + position-keyed roughness) → R2 raster-displacement renderer (zoom floor 20 km) → R3 coastline-conditioning design → R4 regional bake → R5 exports.
- Fresh-context review: **go with changes**, all adopted. Revised order: R2 raster renderer on today's macro (+ river overlay, coverage-based zoom floor) → R3 coastline design → R1 jitter-only experiment → R4 → R5.
- Matt approved (s03): go on R2; Q9 = in-house (FOSS Gaea alternative); native-first CLAUDE.md edits done; `node_modules/` → Trash, `target/` cleaned (40 GB free). `spike/` + `heightmap-1.*` → `~/.Trash/realmgenesis-s02-leftovers-20261104/` (Matt). Spike PNGs also in `~/RealmGenesis/spike/`.
- **R2 built and accepted by Matt** (`3b819c0`; clippy clean, workspace tests pass; zoom floor OK at 32): macro-only display (CPU bicubic heights per vertex) plus per-pixel lighting from a baked gradient texture, detail toggle gone, level cap, zoom-floor slider (default 32 px per texel). Design + review: `docs/design/raster-renderer.md` (the "revised design" section at the end).
- **Since R2 (s03):**
  - the sun is fixed to the ground (NW, 40°);
  - lat/lon grid with labels;
  - the detail-network code is deleted (−3k lines; tiles ~1 ms);
  - **uplift blend** slider: the "floors"/straight escarpments were Voronoi steps in the convergence-rate factor; default 200 km, 0 = old look;
  - "uplift focus" → **belt width** (the old slider did nothing);
  - in-app screenshots: F12, or `RG_SHOT=<png> RG_SHOT_EXIT=1` + `RG_EXAG`/`RG_KEEP_RELIEF`/camera envs. Use these instead of `screencapture`, which fails here.
- **Resolution coherence fixed at the root:**
  - `plates.rs`/`relief.rs` picked plate seeds and hotspots by random node *index*, so plates moved with resolution;
  - now `Grid::nearest_node(rng.unit_vec())`;
  - seed 1, 320 vs 640: height r 0.38 → 0.93, same plate 7.7% → 99.7% (`probe_res_coherence`);
  - every regenerated world changed.
- **R1 spike in the app:** "drainage jitter" slider (0–0.4, default 0).
  - At 0.35: the straight parallel gully grain goes away (Ultra, n=1); river coherence across resolutions 52% → 59%;
  - cost +22% erosion time, +0.3 GB at Ultra.
  - Position-keyed roughness (`roughness_wavelength_km`) had no measurable effect; keep it 0 (candidate for removal).
- **The app folder `~/RealmGenesis` is this repo** (case-insensitive macOS), so saves, exports and F12 shots land in the repo root (gitignored now). Proposed: move it to `~/Documents/RealmGenesis`.
- **Plateaus fixed and tunable** (s03): open lakes used to be fixed, never-eroding base levels, which carried dry interiors up as noisy plateaus. Now only closed lakes, plus open lakes drier than the **"plateau dryness mm"** slider (default 60), hold.
- **Legend** for each layer (bottom-left of the viewport, "Legend" checkbox; `rg-pipeline/src/legend.rs` holds the data, `ui/mod.rs` draws it). Unseen by me; Matt to check.
- **Decided (s03):**
  - **Continents first**, plates fitted to them by a plausibility pass (R3);
  - drainage jitter off and its slider removed (code kept for R4).
- **Next session: R3 design** (`docs/design/coastline-conditioning.md`, then a fresh review):
  - coastline import (vector/raster), align/scale on the sphere, edit, regenerate;
  - the plausibility pass: active vs passive coasts, boundaries + trench/arc/mountains on active coasts, shelf/rift geometry on passive ones;
  - **open question for Matt:** pick active coasts automatically (ocean-facing, a share per seed) or paint them on an imported map?
  - plate/coast alignment today: plates grown after the continents with only a soft margin wall (`plates.rs`, `crust.edge_m` cost).
- Parked: tall flat-topped ranges (`elevation_limit_m` taper?); coast staircase (coast SDF); erosion's remaining resolution dependence (rivers ~59% overlap).

## Previous: project-wide review brief (s02 wrap-up, 2026-11-04)

Matt wasn't satisfied with how terrain *detail* turned out ("the base maps plus the approach to geometry make decent terrains anyway. It's the extra details that we're running into problems"). He expected **raster displacement maps** for geometry. The next session is a **broad review of the whole project** before any further building.
- Start by reading the "Wrap-up" section at the end of `docs/session-logs/2026-11-03-s02-local-m1.md`: **12 grilling questions** (answers pending; ask Matt for them first), the B / M0 spike results, and the raster-displacement discussion.
- Context to weigh: `ROADMAP-V2.md` on `main` (Oct 7, raster-first paradigm from tcotg, Python + TS plan) vs what s01 built (Rust + wgpu, procedural per-tile detail); `docs/design/regional-detail.md`; `docs/design/s02-progress-console-and-macro-mode.md`.
- Get a fresh-context reviewer's verdict on whatever direction comes out (project rule).

## Local session quickstart (s01 wrap-up)

- Branch `claude/peaceful-brahmagupta-rkl4y0` is the whole of v2. It's an **orphan history**, so it shares no commits with `main` (V1, TypeScript). Don't PR it into `main`; check it out directly:
  - `git fetch origin claude/peaceful-brahmagupta-rkl4y0`
  - `git switch claude/peaceful-brahmagupta-rkl4y0`
  - To give it a friendlier name locally: `git switch -c v2 origin/claude/peaceful-brahmagupta-rkl4y0`
- Rust toolchain pinned to 1.97.0 (`rust-toolchain.toml`; rustup picks it up). Web builds need `trunk` and the `wasm32-unknown-unknown` target.
- Run: `cargo run -p rg-app --release` (native). Tests: `cargo test --workspace --release` (~2 min; `terrain_realism` takes ~45 s). CI (`.github/workflows/build.yml`) was green on the last s01 commit.
- Science probes (all in `crates/rg-pipeline/examples/`):
  - `probe_realism`: the scorecard (`RG_PROBE_SURFACE=net`);
  - `probe_hypsometry`;
  - `probe_bathymetry`;
  - `probe_lakes`;
  - `probe_climate_erosion`;
  - `earth_reference` (real-Earth reference; fetches AWS terrain tiles with curl);
  - `endorheic_compare`;
  - `river_map` (PNG hillshade; `RG_MAP_FULL=1`);
  - `bench_tiles`.
  - Probe worlds are cached in `target/probe-cache/`.
- **On Matt's M1 (s02):** see CLAUDE.md → "Local (M1) sessions": builds/tests capped with `taskpolicy -b -j 4`, timing runs uncapped, a lean offline toolchain (minimal profile + `rustup override set 1.97.0`, no wasm32 yet). The full test suite takes ~14 min capped.
- Web preview artifact (claude.ai): https://claude.ai/artifact/1Tq1v7fz6zdbLTtUoLMPME (v7, from before lakes/climate erosion; rebuild with `scripts/build-web-artifact.sh`).

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

**State of the physics (end of s01):**
- **Calibration:** one physical factor (`linear_mean_m` 656, i.e. fitting K); hypsometry emergent and checked. Uplift tapers to zero at 7 km (strength cap). Seed 3: median land 329 m (Earth ice-free ~330–400), > 2 km 10.7% (6–9), > 4 km 3.3% (~2), max 5.4 km (8.85).
- **Climate → erosion:** Thornthwaite PET + Budyko–Fu runoff (land mean 333 mm/yr vs Earth ~300); runoff-weighted incision (wet land erodes far more than dry).
- **Lakes / endorheic basins:** `rg-erosion/src/lakes.rs` (flux fill–spill–merge), run inside the solve as base level, no rock fill. Mass balance exact (~1e-14). Seed 3: endorheic **4.5%** of land (Earth ~20%), lakes 0.25%.
- **Committed fields:** uplift is the *effective* rate; `height_scale`; `lake_depth_m` = water depth (height = bed).
- **Earth reference** (`--example earth_reference`, AWS terrain tiles, same metric code as `probe_realism`):
  - coherence: Earth 0.23–0.59, ours 0.21–0.31;
  - β: Earth 2.3–4.7, ours 2.7–3.1;
  - θ: Earth 0.19–0.45, ours −0.04–0.14 (low);
  - mountain slopes: Earth 26° (steep ranges) / 16° (Appalachians), ours 8–13° in the rugged boxes (low).
- **Bathymetry:** mean depth 3626–3709 m on hold-out seeds (Earth 3688).

**Known gaps (ranked):**
0. **Performance.** Mobile is bad (Matt, web preview v7, end of s01); network tiles are 56–148 ms each on an M1 P-core, single-threaded on the web. **s02 (native, M1):** the app locked up the Mac and dragging crawled. Fixed in code, **not yet confirmed by Matt in the app**:
   - a red STOP BACKGROUND RENDERING button (top-right of the globe) that cancels queued tiles and stops the repaint loop;
   - camera height queries no longer build 240 ms detail patches on the UI thread mid-drag (macro height above 10 km clearance; one quantised cut-off below);
   - tile workers: half the cores (≤ 4), utility QoS on macOS, unwanted tiles skipped;
   - no repaint on mere hover.
   **s02 later:** world generation now runs on a background thread natively (live status line + console; the web build keeps the synchronous path), and **streamed detail is off by default** (macro-only tiles: idle at 0% CPU, 625 MB, n=1). The detail network is one toggle away ("Streamed detail (experimental)"); it draws structural artifacts from orbit at every tile band (s02 log). **Direction (Matt): A now (macro + bigger monolith), B next = ROADMAP-V2's regional "Gaea mode" DEMs** (design doc + fresh review first). Design: `docs/design/s02-progress-console-and-macro-mode.md`. On the web the button pauses the in-frame tile generation. See ROADMAP → Engineering.
1. **Leaky lake rims:** neighbour-basin wedges cross the macro divide below a closed lake's level; water stands as a wall on part of the footprint edge. Fix idea: partition wedges by terminal sink inside closed catchments, with a continuous rim (needs design; min/max compositions can make pits).
2. **Endorheic 4.5% vs ~20%:** Thornthwaite underestimates arid PET; no sediment infill; 1° climate.
3. **Mountains too gentle and θ low:** relief-limited channel heads, effective uplift; no glacial/threshold landscapes as steep as Taiwan/Himalaya.
4. **High ground ~1.5× Earth, peaks low:** no crustal thickening/isostasy; no deep incision above plateaus.
5. **Bathymetry partly circular:** the age remap and the dynamic-topography std were tuned on the validation bands.
6. Sediment transport/deposition, deltas; glacial; ridge–transform fabric (see the scientific review in the log).

**Also next:**
1. **Tile cost:** 85–240 ms per tile natively (gullies: 9–33). Ideas: per-level window reuse across tiles, cheaper attach for the finest level, and eventually the GPU port. The web is single-threaded, so it streams slowly.
2. **Looks:** the uniform flat floodplain band of big rivers; faint arcs (reach walls); verify the globe app natively (Xvfb), including river areas.
3. **Final-surface hypsometry** (mean ~100 m below macro): a strictly increasing remap g with g(0) = 0. Stored receivers instead of flow_dir snapping.
4. Retire the gully/hill code once the app is happy; update `docs/design/terrain-model.md` §3 to the network.
5. **"Alive world" layer:** research in `docs/research/2026-11-02-alive-world-features.md`. Top 5: habitability/population → names → burgs → cultures → states. Licence flags are in the doc.
6. Reviewer leftovers: lithology tied to crust/terranes; χ-analysis.

## Settled decisions (don't re-litigate)

- Rust + wgpu + eframe/egui. The UI is isolated behind commands, and a redesign is planned later.
- Native-first (s03): threads/rayon allowed in core; GPU where it pays. Web build optional.
- Sphere by default (`f = 0`), with the ellipsoid in the same code path. Geodetic latitudes.
- Cube-sphere tile quadtree with the tangent warp; 257² corner-aligned vertices per tile. Committed macro at L3 (≈4.9 km). L4 is sparse macro edits, L10 local edits, L12 (≈9.5 m) was the zoom floor; **s03 replaces it** with a coverage-based floor (one stored texel ≤ ~8 px: ~800 km unbaked, ~25 km over 150 m bakes).
- Sims are stepped `Job`s (async web readback). High-octave noise uses integer lattice coordinates from the tile address. The world file + golden tests land at M2.
- Climate: tier A for the reroll preview, tier B for the one-shot commit. Committed fields are stored, never re-derived.
- Tectonics: a snapshot for the MVP; history mode parked.
- Borders: procedural by default, editable.
- Desktop/laptop only. Minimum spec: M1 Air.
- Rust stays the engine (s03; rewrite rejected, disk managed by pruning `target/`).
- Product (s03): a fantasy-oriented, user-friendly FOSS alternative to Gaea; regional detail in-house; Gaea/Terragen parity is a long-term stretch. Blender/Illustrator/engines are export consumers.

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

- **The 12 review questions** at the end of `docs/session-logs/2026-11-03-s02-local-m1.md` (s02).
- Resolution coherence: key the erosion's initial roughness/Rho8 to position (agreed "after M0"; changes every generated world).

- Paste one `mediasheet-for-obs` session log so the log style can match. It isn't in the public repo.
