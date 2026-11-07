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
