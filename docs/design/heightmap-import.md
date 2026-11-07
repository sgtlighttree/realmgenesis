# Heightmap import and "upscale" (draft, s06, 2026-11-07)

Status: **reviewed (s07, go with changes, all adopted); see "Revised design (s07)" at the end, which supersedes the fidelity table and the Guided mechanism above.**

## Goal (Matt, s06)

- Import a low-res heightmap (from Azgaar, or a hand-made or eroded one) and get a RealmGenesis world that keeps its shape (coasts, ranges, basins) but has real drainage, erosion detail and climate.
- Make it work at every preset: "upscale" means Normal → Ultra and later R4 regional bakes.
- Example files (MATT 3, `32 Skybound Saga/5. Cartography Materials`):
  - **V1** `Klerga (100k) Heightmap`: 7616² RGBA greyscale, Azgaar style. Sea is 0 (71% of pixels); land starts at about grey 50 (Azgaar's 0–100 scale with sea level 20 ≈ grey 51); 171 grey levels, visibly cell-stepped.
  - **V2** `Klerga V2 Maps HeightmapMono`: 11200×5600 (2:1, a whole world). Already eroded-looking; sea level unknown (56% of pixels below grey 16).

## Model: a mask with an elevation channel

A heightmap import is a coastline import (R3.2) with heights:

- **Land/sea:** the existing threshold. The "land above" grey becomes the **sea-level grey**, so the import is a strict-coast layer: the plausibility pass, margins, trenches and shelves all work unchanged.
- **Elevation:** land heights come from grey above sea level:
  `h = max_m · ((g − g_sea)/(1 − g_sea))^γ`.
  - Inputs: `max_m` (default 4,000 m) and `γ` (default 1; > 1 flattens the lowlands).
  - **Presets:** Azgaar (g_sea = 20/100, the γ that matches Azgaar's own height display), Linear, Custom.
- **Sea floor:** by default our bathymetry (shelves, slopes, trenches from R3). Optionally the import's sub-sea greys as relative depth, if a map has them (V1 doesn't; sea is flat 0).
- **Layers:** a heightmap layer sits in the same layer list as coastline layers (later wins inside its patch). Its elevation applies only where it is the topmost enabled layer.

## Where it enters the pipeline

Today: tectonics composes relief, the sea-level solve (or the strict mask) sets the coast, uplift drives erosion at `canonical_res` (320), then `generate_refined` adds sub-cell detail for finer presets (`refine`: up(H_canon) + zero-mean tectonic detail + a short no-uplift incision).

Proposal: the imported heights replace the **composed land relief** inside the layer, before erosion. A **fidelity** control decides how much erosion may move them.

| fidelity | what runs | result |
|---|---|---|
| **Exact** | No macro erosion inside the layer: H_canon = the import (drained with priority-flood so rivers can route). Refine adds sub-cell detail and the short incision. | Shape preserved to the pixel at canonical res; valleys and detail only below a cell. |
| **Guided** (default) | Erosion runs with **uplift derived from the import** (U ∝ the imported height, so ranges hold), starting from the imported surface, then calibration is **pinned** so heights stay within ±x% of the import. | Real river networks and valley shapes, close to the drawn relief. |
| **Loose** | The import is only an initial surface plus the coast; the normal tectonic uplift applies. | Most natural; drifts furthest from the drawing. |

"Upscale" is then just the existing presets: the same import generated at High, XHigh or Ultra goes through `generate_refined`, and R4 regional bakes reuse `refine` on the imported macro.

What tectonics still does inside an imported layer:
- plates and margins from the mask (needed for trenches, bathymetry, and uplift outside the layer);
- **no** tectonic uplift or relief on land inside the layer in Exact/Guided;
- outside the layer, everything is unchanged.

## Data and memory

- **Downsample at import:** to ≤ 4,096 px on the long side (V2 becomes 4096×2048 f32 = 32 MB), box-filtered. Keep the original file bytes in the world file, as today. The 11200×5600 V2 decodes to ~250 MB RGBA + 250 MB f32 today, and `make_thumb` re-decoding it on the UI thread would freeze the app. Moving the decode to a worker thread is part of this work.
- **Cell steps:** Azgaar's stepped heights (171 levels on cells) need a light blur (≈ one source cell) before use, or the steps come out as terraces. Make it a "smooth" slider, default on for the Azgaar preset.
- **World file:** the layer's header entry gains `elevation: {max_m, gamma, g_sea, fidelity, smooth, sea_floor}`, all optional. Old files are unchanged.

## UI

In the import panel, a **kind** switch: Coastline / Heightmap. For Heightmap:
- a preset;
- the sea-level grey (the threshold slider, relabelled);
- max height and the curve;
- fidelity (Exact / Guided / Loose);
- smooth;
- the preview tinted as a hypsometric ramp instead of green/blue.

Placement, Move patch, crop, layers and undo all stay as they are.

## Gates (before calling it done)

- Exact at canonical res: |H − import| < 1% of `max_m` on land (excluding the drained pits); coast sign matches the mask exactly (existing strict gate).
- Guided: correlation with the import r > 0.95 at canonical res; rivers reach the sea (no undrained land outside lakes); hypsometry within the existing test bands.
- Upscale coherence: the s04 guard (Normal vs XHigh, r ≥ 0.99, same coast) holds for an imported world.
- V1 renders without terraces at XHigh (eyes, `river_map`).

## Matt's answers (s06, 2026-11-07)

- **Q1 fidelity: Guided** by default.
- **Q2 elevation:** the heightmap adjusts to the values set in the panel (max height, curve). There's no fixed meaning of white, so the panel's max height is the scale.
- **Q3 V2 sea level:** the grey Matt pointed at is the coast step. In his screenshot it measures ≈ grey 62–71 at the shelf edge against ≈ 40 beyond it. The file's histogram has a dip at grey 48–64 (0.4% of pixels), the sea-level step. So **g_sea ≈ 60/255 ≈ 0.235**; confirm on the file. **V2 includes bathymetry** below that grey.
- **Q4 sea floor: use the map's sub-sea greys when present**, else ours. (V1 has none; V2 has them.)
- **Q5: regional patch by default.** Azgaar exports have no projection, so a heightmap is placed like a patch (azimuthal equidistant about its centre), never as a 2:1 world unless the user picks that.

## Open questions (answered above, kept for the record)

1. Default fidelity: **Guided** (real valleys, close to the drawing), or **Exact**?
2. What max elevation did the Klerga maps mean? Azgaar's 100 is about 6–8 km in its own display, while ours defaults to 4 km.
3. V2's sea level: which grey is the coast? Or should it come from V1 or the Klerga coastline mask?
4. Sea floor: always ours, or use the map's sub-sea greys when they exist?
5. Is a heightmap always a whole continent or world, or do you also want small regional patches (which overlap R4)?

## Revised design (s07, 2026-11-07, after the fresh review)

Fresh-context review verdict: **go with changes**. Adopted:

1. **Calibration would erase the import (blocker).** `finish` soft-limits land above 3,800 m; the solver p99-rescales land, then fits one mean factor (`linear_mean_m`) or quantile-maps it; uplift tapers with `elevation_limit_m`; refine scales detail by a global factor. Fix: a per-node **pin weight** w ∈ [0, 1] (1 inside the heightmap layer, feathered at the patch edge). Pinned land skips the soft limit, the p99/mean/quantile calibration (fitted on free land only) and the taper; refine uses scale 1 there.
2. **Guided = inverted uplift**, not "U ∝ height". The 30 Myr solve reaches steady state (channel response ~0.1–0.3 Myr), so it forgets its initial surface. With n = 1 the implicit step's fixed point is `h − h_r = U·L / (K·A^m)`, so `U_inv = K_i·A_i^m·(h_i − h_r)/L_i` on the import's own filled network makes the import a steady state. One mechanism, two settings:
   - **Exact:** U = U_inv, no erosion steps (routing, lakes and stats only).
   - **Guided:** U_inv smoothed and clamped (≥ base rate, within [0.2, 5]× a local median), mixed with the tectonic U (λ ≈ 0.7), eroded from the import, then a post-solve clamp `h = import + clamp(h − import, ±max(10%·import, 30 m))`.
   - **Loose:** deferred (it's today's mask import plus an initial surface).
3. **Exact vs pits:** filled pits and playas are excluded from the 1% gate and reported separately; closed basins go through the lake model (bed kept).
4. **Sea floor:** `depth = D_max·((g_sea − g)/g_sea)^γb` (D_max default 4,500 m), blended with R3's shelf/slope formula by the pin weight; keep the strict `≤ −1 m` and the soft [−200, −1] compression.
5. **Determinism:** store params and the original bytes only (eroded heights are committed already). The downsample is a pure sequential function with a `downsample_ver` in the header. Heights sample a **mip pyramid** at the node spacing (point-bilinear on a 0.6 km/px raster aliases differently at Normal and Ultra and would break the s04 coherence guard). Coverage keeps sampling level 0, so existing coastline worlds stay byte-identical.
6. **Memory:** decode the PNG row by row into a box accumulator (crop, then luma, then downsample to ≤ 4,096); the peak drops from ~500 MB to tens of MB for V2. The thumbnail comes from the downsampled raster, built off the UI thread. Both Klerga PNGs are **8-bit** (checked s07), so the "smooth" slider matters for both (4 km over ~200 greys ≈ 20 m steps). V1's 32-bit TIF is deferred.
7. **Dropped gate:** "hypsometry within the existing bands" (the user sets `max_m`). Added Guided gates: rivers reach the sea, slope–area stats within the existing bands, r > 0.95.

**Build order** (small commits, each with a test): (1) streaming decode + downsample + pyramid; (2) `Elevation` params on the layer + header round trip; (3) imported heights in `finish` with the pin weight and the sea-floor blend; (4) pinned calibration in the solver + refine; (5) Exact; (6) Guided; (7) feather, UI, presets.
