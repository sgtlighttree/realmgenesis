# s03: project-wide review (2026-11-04, local M1)

Context: s02 wrap-up (`2026-11-03-s02-local-m1.md`, "Wrap-up"). Matt unhappy with terrain *detail*; asked for a broad review before more building.

## Matt's answers to the 12 grilling questions (near-verbatim)

1. **Reference look:** "Basically a Photoreal Earth. I don't expect VFX/Hollywood grade levels of detail for terrain but I expect very naturalistic terrain but only at a lower resolution."
2. **Zoom range:** infinite zoom is "honestly just a habit. We can cap it to a reasonable limit. Your call, then I'll let you know where to push the floor."
3. **Product:** "Worldbuilding tool, for sure. Basically an all in one app for generating a fantasy world", hence the roadmap's civilizations and fantasy map.
4. **Science vs look:** "Science first, but also malleable." Wants to **import vector/raster coastlines** (e.g. from other fantasy map generators), align/transform/scale them on the sphere, and have the generator adapt.
5. **Interaction:** drawing isn't necessarily core. "Editing coastlines and regenerating terrain would be good."
6. **Finite resolution:** raster world (~5 km everywhere, sharp only in baked regions) is acceptable: "Also makes optimization easier."
7. **Stack:** offered the s01 cloud transcript (he asked Rust vs Java for cross-platform, forgot to ask about Python). Would like to save disk by dropping the Rust tooling, but sees it's a major rewrite.
8. **Web/mobile:** only mattered when the project was architected in the cloud session. **Native first for now.**
9. **External tools:** exports are for **Blender** (planetary rendering), **Illustrator** (vector maps) and **game engines** (interactive). (Didn't directly answer "Gaea/WM does regions"; my reading below.)
10. **Keep:** "the bones are there for the terrain generator, it just doesn't look right when we add more detail." Will flag other things as they come up.
11. **Waits:** overnight is overkill; an Ultra reroll on the M1 is fine as is (~40 s measured in s02).
12. **Usable:** "when I feel like it is." Terrain first: it's the first thing he sees and it feeds everything downstream (civs, borders, cultures, cities, à la Azgaar's).

## Findings this session

- **Disk (n=1, 2026-11-04):** 38 GB free on `/` (short of the 40 GB ideal). Rust footprint: `~/.rustup` 1.4 GB, `~/.cargo` 188 MB, `target/` 1.9 GB (≈3.5 GB total; `target/` is reclaimable any time). Untracked in the repo: `node_modules/` **372 MB** (contents look like V1's TS deps: `@dimforge`, `@google`, eslint; probably left over from a `main` checkout, unconfirmed), `spike/` 56 MB (s02 spike PNGs), `heightmap-1.png/json` 8.4 MB (an app heightmap export?). Asked Matt before touching them.
- **Macro grid confirmed:** `rg-tectonics/src/grid.rs`: cube-sphere nodes, 8-connected. So the macro planform's straight runs/45° turns are D8-on-cube-sphere, as s02 concluded.
- Committed the lever-round spike code (`02213a0`, was uncommitted at the s02 wrap).

## Decisions / readings (pending Matt's confirm where marked)

- **Q7 → keep Rust (recommendation).** Rewriting tectonics, erosion, lakes, climate and the app to save ~2 GB net is not worth it: a Python + numba + node/three.js stack is itself ~1–1.5 GB. The s01 transcript isn't needed: Rust was a deliberate choice over Java. Disk lever instead: `cargo clean` / prune `target/`, quarantine the stray `node_modules/`.
- **Q8 → native-first relaxes two settled rules:** "GPU-first because the web has no threads" and "CPU fallbacks single-threaded on wasm32; rayon native-only". The web-artifact preview loop becomes optional. CLAUDE.md edits **proposed to Matt, not made**.
- **Q9 reading (confirm):** regional detail stays in-house (all-in-one tool); Blender/Illustrator/engines are export *consumers*, not the regional generator.
- **Q2 → zoom floor proposal:** fixed **20 km camera altitude** to start (finest stored data will be ~150–300 m baked regions, ~5 km elsewhere; at 20 km a 150 m texel spans ~10 screen px, still readable with per-pixel normals). Matt pushes from there.
- **Q4/Q5 is the architecture-changing answer:** today coastlines *emerge* from tectonics; Matt wants terrain *conditioned on a land mask he supplies or edits*. Direction doc: `docs/design/v2-direction-s03.md`.

## Fresh-context review of the direction: go with changes

All points adopted; detail in `docs/design/v2-direction-s03.md` → "Fresh-context review".
- **Correction:** the macro is already Rho8 (`solver.rs:121`), not plain D8, so "D8 griddiness" was imprecise. The planform is still bound to the 8 lattice neighbours.
- **Refuted (my proposal): a flat 20 km zoom floor.** At 20 km a 5 km texel is about 300 px. New rule: one texel spans at most ~8 px, which gives ~800 km over unbaked land and ~25 km over 150 m baked regions.
- **Order changed** to R2 (raster renderer on today's macro + river overlay) → R3 design → R1 jitter-only experiment → R4 → R5. Rejected: R1 first. R2 doesn't depend on jitter, and a 3D view is the better judge.

## Matt's go (s03)

- **Go on R2** (raster renderer on today's macro).
- **Q9 confirmed:** "RG becomes a fantasy oriented, user friendly, but less feature rich FOSS alternative to Gaea." Feature parity with Gaea/Terragen is long-term.
- **Native-first:** CLAUDE.md "Code rules" and "Preview loop" updated; HANDOFF settled decisions updated (native-first, Rust stays, product line, zoom floor).
- **Disk:** `node_modules/` (V1 leftovers) moved to `~/.Trash/node_modules-realmgenesis-v1-20261104`; `cargo clean` freed 1.9 GB (`target/probe-cache` worlds included; they rebuild on the next probe run). Free space now **40 GB**.
- Side note (side agent, agrees with the s02/s03 finding): raster rendering alone won't fix the griddy planform. That's R1's job. Saved `.rgworld` files store committed fields, so an old seed-1 save should still load after R1 (unverified; schema permitting), but regenerating seed 1 won't reproduce it.

## R2 design (raster renderer)

- Draft: `docs/design/raster-renderer.md`. The key finding from reading the code: the quadtree, camera-relative instances, skirts, geomorph and exaggeration already exist (`rg-render/src/terrain.rs`, `rg-terrain/src/tile.rs`), so R2 only changes **where heights come from**: a baked 6×2049² `DisplayRaster` (R32F height + Rg16F gradient) on the GPU instead of per-vertex CPU sampling.
- **Expectation flag for Matt:** at the same data, R2 looks about like today's macro-only Ultra (which he already called "pretty decent" in s02). Its payoff is no stalls, per-pixel normals at every zoom, and the plumbing for R4/R5. Visible *detail* needs R1/R4.
- Fresh-context review: pending.
- **R2 review: go with changes.** Adopted a simpler design: heights stay per-vertex on the CPU but come from the `DisplayRaster` (cheap Catmull-Rom), plus a fragment-only baked gradient texture for per-pixel normals, plus one surface for tiles, camera and hover, plus a level cap at L6 and the zoom floor on `distance_m` using the raster texel.
  - **Rejected:** vertex-shader displacement for now. It needs the bind-group changes, a skirt-parity workaround, a corrected (odd, odd) diagonal geomorph, and a min/max pyramid for the culling/split errors, all for no visible gain at a max level of 6.
- Matt: "clear them" → `spike/` and `heightmap-1.*` moved to `~/.Trash/realmgenesis-s02-leftovers-20261104/`.
- ✓ R2 step 1: `rg-terrain/src/gradient.rs` (`bake_gradient`, `gradient_res` = 2× macro res within 512..2048). Tests pass: analytic field < 2% error, face-edge twins agree. Gotcha: the cargo registry cache was empty again (0 B; Mole?), re-fetched 21 MB.
- ✓ R2 app wiring `3b819c0`:
  - face uv per vertex and an Rg16Float gradient array in `fs_tile`;
  - normals bake job (`rg-app/src/normals.rs`; flat lighting until it lands);
  - detail toggle removed, river overlay fade removed;
  - level cap = ceil(log2(2·res/64));
  - zoom floor on `distance_m` (slider "zoom floor", default 32 px per texel).
- **Measured (n=1, seed 1 Normal, M1 on AC):** world 3.8 s, **normals bake 0.2 s** at 640 per face. Build at `-j 6` from a clean target: 1 m 38 s.
- `screencapture` fails from this terminal ("could not create image from display"; probably no Screen Recording permission). Matt checks the look.

## Matt's look at R2 (after the iTerm2 restart; screenshot 11:53, seed 1 Ultra, Elevation layer, ~36°N 150°W)

- "Ultra looks the same" (expected: same data). **Zoom floor OK** at the default (32 px per texel → 402 km on his view). No seams seen so far; he asked where to look (answered: the cube edges, see below).
- "Apart from some **oddly straight mountain ranges**, nothing else looks too off." His arrows mark an E–W ridge with parallel N–S gullies, and a ridge front running SW–NE.
  - My read (unconfirmed): two causes stack here.
    1. The plate-boundary uplift is straight and smooth (a known gap: "smooth Voronoi-ish plate borders").
    2. The gullies line up with the face grid axes. The spot is on the NegX face, where the grid's v axis runs ≈ north. That's the lattice-bound Rho8 planform, which is R1's target.
- **The resolution dependence is still there** (expected: R2 doesn't touch generation). Fix = position-keyed roughness/Rho8 (parked since s02), folded into R1.
- Cube edges (for the seam check): PosX is centred on lon 0 (lat = asin z, lon = atan2(y, x)). So the edges are the meridians **lon ±45° and ±135° between ~35°S and ~35°N**, plus a square around each pole whose corners sit at **lat ±35.3°, lon ±45°/±135°**.
- ✓ **Workspace tests pass after R2** (all suites ok; 142 s wall at `-j 6 --test-threads=2` on AC, with Matt's app open; peak RSS 4.0 GB per `/usr/bin/time -l`, n=1). Pending Matt: R1 before R3? Delete the detail-network code?

## Matt's round 2 on R2 (12:00–12:10)

- Lighting "changes slightly as I rotate… light comes from behind the camera… reduces contrast". **Cause:** `app.rs` tied the sun to the camera (screen upper-left, mostly overhead: up·1.2). **Fix (uncommitted at first, then in the next commit):** the sun now sits 40° up, from the north-west of the orbit point, independent of heading.
- Wants a **lat/lon grid with numbers** (status bar fine) → delegated (sonnet-medium, worktree).
- Approved: **R1 before R3**; **delete the detail-network code**.
- **Screenshot 2 (37°N 96°W, Ultra seed 1): a stair-stepped straight escarpment.**
  - New tool `examples/field_maps.rs` (6 panels: hillshade, uplift, runoff, lake depth, plates, crust; nearest-node sampling).
  - **Found:** the committed uplift field has a hard straight edge exactly where the texture changes. It doesn't follow the runoff cells (bilinear from res 128) or the plate boundary.
  - **Cause, verified (n=1):** `uplift.rs` scales the belt term by `factor(near.src)`, the convergence-rate factor of the *nearest* convergent seed (clamped 0.25–2.5). The nearest-seed regions are Voronoi cells, so the uplift jumps along straight bisectors. With `RG_CONV_SCALE=0` (new probe env), the edge disappears at the same spot.
- **Screenshot 3 (whole map): "floors" / massive plateaus, and some ranges form super-tall plateaus.** "Interesting, not natural… worth a slider… not something to fix, simply tune."
  - Likely the same Voronoi uplift steps (unconfirmed globally).
  - The tall flat-topped ranges are probably the uplift taper to 0 at the 7 km strength cap (`elevation_limit_m`; unconfirmed).
- **Plan:** smooth the convergence factor over a distance in km, exposed as a slider. At 0 it stays terraced (today's "floors"); higher values give natural ramps. Changes every world at non-zero values, so it goes in with R1.
- ✓ **Uplift blend:** `ErosionParams::convergence_smoothing_km`.
  - **Defaults:** 200 for new worlds; files saved before s03 load as 0, which reproduces their uplift.
  - **How it works:** `uplift::smooth_on_sphere` averages the field on a coarse grid with σ/2 cells, runs 6 neighbour passes, then reads it back bilinearly. Cost and result don't depend on the world's resolution.
  - **Unit test:** a step becomes a ramp and the mean is kept.
  - **Same spot, n=1:** the straight edge is gone and the uplift ramps; max uplift in the box went 1.54 → 1.18 mm/yr (before hypsometric calibration).
  - **UI:** "uplift blend km" (0–600).
- **Found:** the "uplift focus" slider (`uplift_exponent`) only affects the legacy elevation uplift source, so it did nothing with the default tectonic source. Repurposed as **"belt width km"** (`belt_width_km`, 150–800). Matt's earlier sessions' tuning of "uplift focus" had no effect.
- **Still to check:** whether the global "floors" go away at 200 (Matt's eye in the app), and whether the super-tall flat-topped ranges come from the 7 km `elevation_limit_m` taper (a candidate second slider).
- ✓ **Detail-network code deleted** (Matt: "go for it"). −2994 lines.
  - **`rg-terrain`:** now only the displayed surface + tiles + gradient. Gone: `detail_net.rs`, `erosion_filter.rs`, `network.rs`, hills/gullies/corridors, and the rg-noise dependency. `TerrainSampler::new(field, radius)` (no seed); `macro_height` = bicubic with the drained coast sign.
  - **`rg-pipeline`:** a single `terrain_sampler(world)`.
  - **Removed probes** (they measured the detail surface only): `bench_tiles`, `probe_drainage`, `probe_grid`, `probe_line`, `probe_pits`, `probe_slope`, `probe_hypsometry`. **Removed test:** `tests/terrain_realism.rs`. It asserted 90 m slopes of the detail surface, which the 10–50 km macro surface can't have. Its realism role moves to R4 (regional bakes), scored with `probe_realism`'s metrics (kept, now on the macro surface).
  - **`river_map`:** `RG_MAP_FULL` is gone.
  - **App:** camera height queries use `macro_height` directly; `camera_query_cutoff_m` and its tests removed.
  - **Tiles:** split error is now the bicubic-vs-flat-triangle gap at cell centres; skirt depth uses the parent-mesh gap (it was 0 in macro-only mode).
  - **Measured (`bench_tile`, synthetic field, n=1):** ~1 ms per tile at L0–L6 (detail network: 56–240 ms).
  - Clippy clean; workspace tests pass.
- ✓ **Lat/lon grid** (sonnet-medium, worktree → cherry-picked as `ef36072`'s copy):
  - `rg-render/src/graticule.rs` (step picker, globe projection, labels; 5 tests);
  - shader lines with adaptive 30/15/10/5/2/1° spacing;
  - a "Lat/lon grid" checkbox, on by default.
  - Gotcha: the agent's worktree started on a TS-era commit, and it reset to the branch itself.
- **Fresh-context review of uplift blend + deletion + graticule: no blocking bugs.**
  - **Low:** below res ~100 the blend's σ comes out wider than asked (the coarse grid is capped at the world res). The doc comment is fixed.
  - **Cosmetic:** globe labels sit at sea level, so with exaggeration they drift off lines over high ground; skirts get grid lines too.
- ✓ **In-app screenshots:**
  - `PlanetRenderer::read_output` copies the resolved frame back from the GPU;
  - **F12** saves to `~/RealmGenesis/shots/`;
  - `RG_SHOT=<png>` waits for world + normals + tiles to settle, `RG_SHOT_EXIT=1` then quits;
  - new envs `RG_EXAG`, `RG_KEEP_RELIEF`, `RG_CONV_SMOOTH`.
  - **Why:** `screencapture` fails from this shell even after Matt granted iTerm2 Screen Recording ("could not create image from window"), and CuaDriver's capture returns nothing despite both TCC grants (probably Tahoe's per-capture consent; not pushed, it would pop a dialog). The app's own readback needs no permission.
- **Verified by eye (n=1, seed 1 Ultra, 37.2°N 96.1°W, 600 km, relief 60× kept):**
  - blend 0 shows the straight E–W escarpment Matt saw; blend 200 removes it.
  - Shots in `~/RealmGenesis/shots/s03-blend{0,200}-37n96w.png`.
  - At 60×, a fine diamond crosshatch shows in places (the lattice/Rho8 grain, R1's target) and the coast steps along macro cells (drained-interpolant sign; a coast SDF is the known fix, parked).
- Perf note (n=1 each, unconfirmed cause): Ultra erosion took 58–63 s in this session's runs vs 28–32 s earlier. Matt's app instance may have been running concurrently. Normals bake at 2048: 12.6 s.

## R1 spike (s03)

- **Code** (`d…` spike commit, off by default):
  - `rg-erosion/src/jitter.rs`: `node_jitter` moves erosion nodes ±x cells (face-edge ring fixed, keyed on the integer sample index, min/max edge recomputed, Rho8 off);
  - every field written back is resampled by barycentric interpolation on the jittered triangles, then heights are Priority-Flood filled on the regular grid;
  - `roughness_wavelength_km` = position-keyed fBm roughness;
  - unit tests: resampling a smooth field < 5 m error; only interior nodes move.
- **Look test (n=1, seed 1, 36.3°N 150.2°W, 1500 km, relief 35× kept, blend 200; screenshots via `RG_SHOT`):**
  - **Ultra, jitter 0.35:** the straight, parallel NW–SE gully grain of the baseline (Matt's "oddly straight ranges" spot) is mostly gone; ridges and valleys curve and branch; rivers wander.
  - Adding roughness λ 150 km changes little at Ultra.
  - Normal is too coarse at this zoom to judge.
- **Cost (Ultra, n=1):** erosion 60.6 s → 74.0 s (+22%); peak RSS 2.71 → 3.03 GB.
- Gotcha: zsh doesn't word-split `$v` in `shot $v`, so the first batch ran six default worlds. Use `${=v}` or separate args.

## Resolution coherence: root cause found and fixed (s03)

- New probe `examples/probe_res_coherence.rs <seed> <res_a> <res_b>`. Over A's land it reports:
  - height r (B averaged over A's cell);
  - river overlap (A rivers ≥ 20,000 km² that have a B river ≥ half the area within one A cell);
  - same-plate share;
  - convergent-boundary match.
  - Probe env `RG_EROSION=0` turns erosion off.
- **Before (seed 1, n=1 each):**
  - full pipeline 320 vs 640: **height r 0.38**, rms 923 m, rivers 41%;
  - jitter or position-keyed roughness barely moved it (0.37–0.38; rivers 41–45%);
  - **tectonics alone (erosion off): r 0.25 (320/640), 0.19 (640/1024); same plate 7.7%, convergent boundaries 1.1%.**
- **Refuted (s02 hypothesis):** "the resolution dependence comes from index-hashed erosion roughness/Rho8". Not the main cause; roughness keyed to position changes nothing measurable.
- **Root cause:**
  - `plates.rs:84` picked plate seeds as `rng.index(grid.n)`, a random node *number*, which lands elsewhere at every resolution;
  - the same in `relief.rs` hotspot chains.
  - Coastlines matched because the continents come from position-keyed noise, but plates (so mountain belts, uplift and everything after) moved.
- **Fix:** `Grid::nearest_node(rng.unit_vec())` (random place → nearest node); `Rng::index` is unused now and annotated.
- **After (n=1):**
  - **tectonics r 0.98 / 0.99**, same plate 99.7%, convergent boundaries 97%;
  - **full pipeline 320 vs 640: height r 0.93 (rms 332 m), rivers 52%; with jitter 0.35: rivers 59%;**
  - roughness λ 150/50: 58–59% (no effect).
- **Consequence:** every regenerated world changes, seed 1 Ultra included. Saved `.rgworld` files load unchanged (committed fields).
- Remaining river incoherence (~41–48%) is erosion's own resolution dependence (stream power's intrinsic dependence and index-keyed Rho8/jitter); unconfirmed which dominates.
- **Found:** `~/RealmGenesis` (the app's save/export/screenshot folder, `files.rs::folder`) is the repo itself on macOS's case-insensitive filesystem. That's why `heightmap-1.*` and `spike/` showed up untracked. App output patterns are now in `.gitignore`. Proposed to Matt: move the app folder (e.g. `~/Documents/RealmGenesis`).

## Plateaus ("floors") investigation (s03, Matt: "toned down, need a way to tune how often they appear; flat areas are noisy")

Spot: seed 1 Ultra, 50°S 165°E, 2800 km box, `field_maps` (n=1 each).
- **Pre-erosion** (`RG_EROSION=0`): the tectonic surface has no plateaus or escarpments, just rolling interior with fine noise. So erosion makes them.
- Plateau = the original surface, standing above eroded basins behind sharp escarpments. Plateau-only speckle of tiny lakes (lake-depth panel).
- **Refuted:**
  1. "Small closed basins block drainage": a min-closed-basin knob (spill if catchment < X km²) changed nothing at 5k or 50k km². Reverted.
  2. "Fronts haven't had time": maturity 120 Myr leaves the fronts in place, and the plateau gets noisier.
- Drainage area panel: plateau and lowland are one integrated dendritic network. Not disconnected.
- **Confirmed: `RG_LAKES=0` removes plateaus, escarpments and speckle entirely** (erosion 8.5 s vs ~28 s).
- **Root cause:** in `solve_full_with_progress`, every lake node (`sub`: level > bed), *open* lakes included, was a fixed base level; its bed "rises with the rock but never erodes", and upstream land graded to it. The noisy, dip-riddled dry interior became lakes and was carried up as plateaus.
- **Fix:** only closed lakes, plus open lakes in dry ground (below), hold a base level; other open lakes route and erode normally. Seed 3 Normal: endorheic 2.4 → 2.0%; closed lakes 0.15 → 0.12%; open lakes 0.30% → ~0. Real closed basins survive.
- **Refuted as a knob:** scaling the open-lake erosion rate (plateau_persistence). Even ×0.001 drains them; plateaus are all-or-nothing per lake.
- **Knob adopted:** `plateau_runoff_mm`. Open lakes with local runoff below N mm/yr are held (dry interiors stay plateaus, cf. Altiplano/Iran). Default 30 (to calibrate); old files load as 1e9 (every open lake held = pre-s03).
- Also: Ultra erosion runs at ~25–28 s again when nothing else is running. The 58–88 s readings earlier were contention (Matt's app, agents' builds).
- **Calibration (seed 1 Ultra, whole map, n=1):** 0 → almost no floors; 30 → a few; 100 → a fair number; 1e9 → the pre-s03 abundance. **Default 60 mm/yr**; UI slider "plateau dryness mm" (0–300). Probe env `RG_PLATEAU`, app env too. `field_maps` panel 6 now shows drainage area (log) instead of crust; it works with erosion off.
