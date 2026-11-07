# Heightmap import and "upscale" (draft, s06, 2026-11-07)

Status: **draft, not reviewed, no code.** It still needs a fresh-context review and Matt's answers to the open questions at the end.

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
