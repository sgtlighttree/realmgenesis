# s07 (2026-11-07, local M1): heightmap import, review and build

Continues HANDOFF → "a fresh-context review of `docs/design/heightmap-import.md`, then build".

## Fresh review (sonnet-medium, read-only): go with changes, all adopted
- Full list in the design doc, "Revised design (s07)". Key points:
  - **Blocker:** calibration (p99 rescale, mean fit or quantile remap), `finish`'s 3,800 m soft limit, the elevation taper and refine's global detail scale all rewrite imported heights. Fix: a per-node pin weight.
  - **Guided = inverted uplift** `U = K·A^m·S` on the import's filled network: the import becomes a fixed point of the implicit step. My own read before the review agreed: 30 Myr ≫ channel response time (~0.1–0.3 Myr at 30 km cells), so "start from the import" alone is forgotten. **Rejected:** "U ∝ imported height + pin ±x%" (no pin mechanism in the solver; steady state imposes its own χ-profile).
  - Exact and Guided are one mechanism (λ, steps). Loose is deferred.
  - Heights sample a mip pyramid (aliasing breaks Normal/Ultra coherence); coverage keeps level 0 (byte-identical old worlds).
  - Streaming row decode for V2 (~500 MB → tens of MB).
- **Checked:** both Klerga PNGs are **8-bit** (V1 7616² RGBA, V2 11200×5600 RGBA), so no 16-bit decode is needed. V1 also has a 32-bit TIF (deferred).

## Advisor (before build): plan holds, fold these in per step
1. Step 1: old imported worlds stay byte-identical: downsample only for the heightmap kind (with `downsample_ver` in the header); coastline masks keep the full-res path.
2. Auto channel needs a whole-image alpha check: two streaming passes (detect, then accumulate).
3. Integer step `ceil(max(w,h)/4096)`, plain box sum (V1 → 3808², V2 → 3734×1867).
4. The header placement stays in source px (the SVG `user_placement`/`raster_placement` pattern), or Move patch jumps.
5. Step 3: V1's flat-0 sea would become a cliff at D_max. "Has bathymetry" = sea greys spread above a small threshold, decided at import and stored per layer; otherwise R3 bathymetry.
6. Step 4: `height_scale` is one committed scalar: grep its consumers before choosing per-node vs derived (possible schema change; log it).
7. Blank ocean + one heightmap has ~no free land: fit p99/mean on free land only, skip when empty; test that setup.
8. Steps 5–6: the fixed point needs the solver's exact routing (same `route`, `rho8_seed`, L, `k_field`, runoff-weighted A, no roughness or taper on pinned nodes). **Write the one-step fixed-point test first** (|Δh| < 1e-6 m).
9. λ = 0.7 lets fitted plates add coastal ranges the drawing lacks: consider λ = 1 default, smoothing as the only deviation, λ a slider. Record the choice.
10. The ±10% clamp may cancel Guided's change: render Exact at XHigh before tuning; flip the default to Guided only once it exists.

## Step 1 built: streaming PNG decode + downsample + pyramid
- `import::raw_from_png(bytes, channel, crop, step)`: rows stream (`next_row`), crop then luma then an integer box step, f64 accumulators, sequential. Auto runs a storage-free alpha pass first (advisor 2). Coastline imports use step 1 → **byte-identical** (test `streaming_decode_matches_whole_frame` against the old whole-frame decoder kept in the test; crop too). `streaming_decode_box_downsamples` checks the block means.
- **Measured, V2 (11200×5600 RGBA, step 3):** 0.74 s, **36 MB peak RSS** (was ~500 MB: a 250 MB frame buffer + 250 MB f32). `examples/probe_import` prints size, timing and a grey histogram.
- `make_thumb` decodes straight into the 384 px thumbnail (still on the UI thread, now ~0.7 s for V2 rather than ~500 MB and several seconds).
- `mask::Pyramid` (2× box levels, trilinear `sample(px, py, footprint_px, wrap)`); test `pyramid_levels_and_sampling`. Not wired yet (step 3).
- **Side note (a side agent, relayed by Matt):** "upscale keeps the shape, not the detail; V2's sub-cell valleys get replaced". True of the original draft, mostly not of the revised one: the import goes into the tectonic snapshot at both resolutions, so refine's `D` term carries the file's own detail down to Ultra's ~10 km cell. Two real points adopted: (a) the 4,096 cap would have made V2 9.8 km/px, about Ultra's cell → **cap raised to 6,144** (V2 at 7.1 km/px); (b) below ~10 km only R4 bakes (re-decoding the stored source at full res) can carry it. Guided does move features within its clamp; Exact doesn't.

## Step 2 built: heightmap layer params
- `mask::Elevation {max_m, gamma, sea_floor, sea_depth_m, smooth_px, fidelity}`, `Fidelity {Exact, Guided}`; `RasterMask` gains `step`, `elevation`, `heights: Option<Arc<Pyramid>>` (smoothed greys; the coast still reads the unsmoothed level 0). `RasterMask::height_m(d, radius, spacing_m)` samples the pyramid at the node footprint.
- `import::layer_from_bytes(.., elevation, step)`; cap `HEIGHTMAP_MAX_SIDE` = 6,144; `has_bathymetry` (sea greys p5..p95 spread > 5 levels; V1's flat sea → false). The header stores `elevation` + `step` + `downsample_ver: 1`; placement stays in source px (`user_placement`/`raster_placement` divide/multiply by `step`). Test `heightmap_layer_round_trips` (byte-identical re-save, placement unchanged); coastline files unchanged.
- **Approximation (accepted):** a raster pixel = `step` source px exactly, except the partial last column/row, so a patch centre can shift < step/2 source px (< half a raster pixel).

## Step 3 built: imported heights in `relief::finish`
- `mask::imported_heights` → `Imported {h (NaN = none), pin}`: the topmost enabled layer covering a node decides (a later coastline layer removes heights under it); pin = smoothstep feather over min(150 km, 10% of the patch's shorter side) from the patch edge, 1 for a world map.
- `finish(.., imported)`: land = pin·import + (1 − pin)·procedural, **after** the soft peak limit (so pinned land skips it); the map's sea floor blends the same way with R3 bathymetry; strict signs re-applied. `TectonicSnapshot.pin` (erosion input, not stored; `file.rs` loads None).
- Gate `heightmap_layer_sets_land_heights` (res 64): pinned land equals the import to 1 cm, top > 5 km (soft limit bypassed), coast signs follow the mask.
- Not yet tested: the map-sea-floor blend (V2 path).

## Steps 4–5 built (merged): pinned calibration + inverted uplift (Exact)
- **Decision: Exact is not "no erosion steps".** In a mixed world (procedural continents + a heightmap patch) the free land still needs its solve. Instead the solve runs normally with `U = pin·U_inv + (1 − pin)·U_tect`; pinned land is a fixed point, free land erodes as before. Rejected: skipping the solve (breaks mixed worlds), or a separate code path per fidelity.
- `solver::invert_uplift` takes the routing from a one-step, zero-uplift run of the solver itself (same fill, lakes, held lakes, Rho8 seed, runoff-weighted area, L), so the fixed point is exact by construction. U may be negative (open-lake beds below their receiver); the **committed** uplift is clamped ≥ 0.
- **Gate (advisor 8):** `inverted_uplift_holds_the_surface`: 20 steps / 30 Myr move the surface < 1e-6 m, without and with lakes. ✓
- Solver: no roughness and no elevation taper on pinned nodes (`taper(h, lim, pin)`); p99, the mean fit and the quantile remap use free land only (pin < 0.5); each node's factor is `f + (1 − f)·pin`; no free land (blank ocean + heightmaps) → factors 1 (advisor 7).
- **`height_scale` decision (advisor 6):** it stays one committed scalar meaning the free-land factor (no schema change). Pinned land's factor is 1 by definition, and the pin is re-derivable from the stored conditioning. Refine takes `fine_pin` and scales detail by `s + (1 − s)·pin`. Consumers checked: refine, `file.rs`, `regional_spike`, `probe_res_coherence`.
- **Pipeline gate:** `exact_heightmap_survives_erosion` (res 64, blank ocean, dome with ridges): pinned land within 1% of max height outside lakes (≥ 98% of samples), top > 4 km (not rescaled), rivers present. ✓
- Non-heightmap worlds: every change multiplies by `1 − 0` or adds 0, so they should stay bit-identical; the full suite (golden hashes) is the check.
- Full workspace suite passed (capped, exit 0).

## Step 7 (UI) built, before Guided; V1 checked headless
- **Order change:** UI before Guided, so Matt can try Exact on his maps (advisor 10: look at Exact before tuning Guided).
- Import panel (PNG layers): **kind** coastline / heightmap; max height, curve, smooth px (pyramid rebuilt on release), fidelity (Guided labelled "not built yet: behaves as Exact"), "use the map's sea floor" (default from `has_bathymetry`) + depth of black. The threshold slider becomes **"sea-level grey" 0–255** for heightmaps. New heightmap layers default to **Exact** until Guided exists.
- App env for scripted runs: `RG_HEIGHT=max_m,sea_grey[,smooth_px]` (with `RG_MASK`), `RG_BLANK=1`. Probes (`examples/common`): `RG_MASK_HEIGHT=max_m[,smooth_px[,guided]]`, `RG_MASK_SEAFLOOR=0|1`; they print r and the share off by > 1% vs the import.
- **Bug found + fixed:** Klerga V1 has a varying alpha, so Auto picked alpha and every height came out flat. Heightmaps read brightness unless alpha is asked for explicitly.
- **V1 Exact (headless `river_map`, patch 20N 10E, 0.5 km/px, sea grey 50, smooth 3, blank ocean):** step 2 → 3808², r **1.0000**, 0% off by > 1%, 525 MB peak. Shape kept exactly. Rivers are a **straight D8 lattice on the smooth domes** (expected for Exact on a smooth drawing; Guided should fix it). PNG: scratchpad `v1_exact_normal.png` (not kept).
- App screenshot not possible: the screen was locked (`CGSSessionScreenIsLocked` true), and the `RG_SHOT` run idled at 0% CPU waiting for a frame.
