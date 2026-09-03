# F3 Phase 2 — WebGL fill surface for Map2D (Default style)

**Date:** 2026-09-04 · **Status:** design approved, pre-plan · **Milestone:** F3 (true 2D vector map), Phase 2

## Problem

The 2D map (`components/Map2D.tsx`) renders the whole scene — cell fills, all
overlays, and labels — into one offscreen canvas at a capped DPR, then blits it
to screen with the pan/zoom applied as an affine on `drawImage`. Zoom is
therefore **bitmap magnification**: past `MAX_SHARP_SCALE` the offscreen stops
gaining resolution and every pixel is stretched. At the 6× max zoom the surface
runs at ~0.5 device-px per CSS-px even on a retina M1 — soft fills, aliased and
"distorted" cell corners, fuzzy coastlines, blurred labels.

Phase 1 (S30) killed the 2.2s pan/zoom stall and built pre-projected geometry
caches (`utils/mapCache.ts`) sized for exactly this next step. Phase 2 makes the
Default-style 2D map a true vector surface: crisp at any zoom, no magnified
bitmap.

## Scope

**In:** Default map style, non-Dymaxion projections (Mercator, Equirectangular,
Winkel Tripel), all data view modes (biomes, political, province, culture,
religion, height, temperature, moisture, population).

**Out — keeps the current offscreen-blit path, unchanged:**
- **Parchment style.** Its look is Canvas2D hatching, paper grain, and procedural
  relief glyphs; a GPU port is a separate, larger effort (future "Phase 2b").
- **Dymaxion projection.** Stays on its raster source-buffer + pick-buffer route
  (F3 invariant 1).
- **3D globe** (`WorldViewer`). Separate renderer; the F2 ScreenOverlay work
  already covers it.

## Architecture

Two synchronized layers replace the offscreen-blit **only for the in-scope
path**. Everything out of scope, and the WebGL failure case, uses the existing
blit renderer.

```
┌──────────────────────────────────────────────┐
│ VectorOverlay  (Canvas2D, transparent, top)   │  per-frame, crisp
│  coast/border/lake Path2D, rivers, routes,    │  vector strokes + text
│  contours, grid, currents, markers, labels    │
├──────────────────────────────────────────────┤
│ GLFillSurface  (Three.js OrthographicCamera)  │  per-frame, crisp
│  earcut-triangulated cell fills, per-vertex   │  fills, 1 draw call
│  color; pan/zoom via one matrix uniform       │
└──────────────────────────────────────────────┘
        both driven by the SAME (scale, offset, dpr, mirror)
```

### Layer 1 — GLFillSurface (GPU fills)

- **Library:** Three.js (already a dependency via `WorldViewer`). Raw
  WebGL/regl is parked as a possible later swap. `WebGLRenderer` with
  **`antialias: true`** (Three's default is false — without MSAA the fills ship
  crisp-but-jaggy and lose the point).
- **Camera:** `OrthographicCamera` over CSS-px space. Pan/zoom and the
  horizontal mirror are expressed as the camera/model matrix; nothing is baked
  into vertices.
- **Geometry, built once per (world, projection, size):**
  - Earcut-triangulate each cell's projected **sub-rings** (see Unit
    `tessellate`). Voronoi cells are convex on the sphere but **not** after
    projection (Mercator polar stretch) or antimeridian clipping (a sub-ring
    bounded partly by the map edge). Centroid-fan triangulation of a non-convex
    ring produces overlapping/outside triangles precisely at cell corners — the
    exact artifact we are removing — so earcut per sub-ring, one path, no fan
    fast-path in v1.
  - `position`: Float32 `[x, y, 0]` per triangle vertex, in un-flipped CSS-px
    (same convention as `cellVerts`).
  - `color`: per-vertex RGB, filled from the per-cell color cache (Unit
    `mapFillColorBuffer`).
  - Store a **per-cell triangle-vertex range** (`cellTriRange: Uint32Array`,
    `[start_i, end_i)` into the vertex buffer) so paint strokes can rewrite one
    cell's colors without re-uploading the whole attribute.
- **Per frame:** update the camera/model matrix from `(scale, offset, dpr,
  mirror)`; one `renderer.render`. No geometry rebuild on pan/zoom.
- **Recolor (view mode / season / seaLevel / hillshade / paint change):** rewrite
  the `color` attribute only (topology unchanged) and flag it for GPU re-upload.
  Paint strokes rewrite only the affected cells' `cellTriRange` slices. The Three
  buffer partial-update API (`needsUpdate`, `addUpdateRange`/`updateRange` — it
  changed across recent Three versions) **must be pulled fresh from Context7**,
  not written from memory, and pinned to the installed version in `package.json`.
- **Disposal:** geometry, material, and renderer are disposed on unmount and on
  world/projection change, following the `WorldViewer` "every useMemo geometry
  has a matching disposal effect" invariant.

### Layer 2 — VectorOverlay (crisp linework + labels)

- A transparent Canvas2D layer stacked over the GL canvas, cleared and redrawn
  each frame.
- **Transform:** the SAME composition the current offscreen uses —
  `setTransform(dpr…)` → `translate(width, 0)` → `scale(-1, 1)` — **composed with
  the pan/zoom affine** `(scale, offset)`. Because that reproduces today's
  drawing context exactly (just at a live transform instead of baked-then-blitted),
  every existing draw call runs unchanged and lands crisp: `drawMapLabels`,
  `drawRiverPaths`, route/contour/grid/current/marker passes, and the cached
  `coast`/`borders`/`lakes` `Path2D`.
- **The mirror (load-bearing, see `docs/map-styles.md` §"The mirror trap"):**
  polygons and strokes come through the flipped context and land correctly; text
  and glyphs are compensated *inside* their own draw functions (double-flip →
  identity). Reusing the identical transform means labels stay forward-facing
  with zero new mirror logic. Do not hoist text into an unflipped group.
- **Boundary LOD:** the cached `coast`/`borders`/`lakes` paths are
  Douglas–Peucker-simplified at `tol = 0.5px` **at base fit** (`mapCache.ts`),
  so at 5× they carry ~2.5px of visible wobble — crisp fills over wobbly
  coastlines is still half the complaint. On **zoom-settle**, rebuild those three
  paths at `tol = 0.5 / scale`. This is normal vector-tile LOD (settle-refined
  *geometry*) and is categorically different from the raster-sharpness bug being
  removed (settle-refined *pixels*): the geometry is resolution-independent at
  every zoom; only its simplification budget tracks zoom. State this distinction
  in code comments so it does not read as reintroducing the bug.

### Layer 3 — Fallback

- Keep the current offscreen-blit renderer as the fallback path. Engage it on:
  1. WebGL context **init failure** (creation returns null / throws), and
  2. **`webglcontextlost`** mid-session (plausible on a fanless M1 under thermal
     load with After Effects running) — listen and fall back live.
- Out-of-scope paths (Parchment, Dymaxion) route to the blit renderer by the
  same selection logic; they are not "fallbacks," just the unchanged path.

## Units

Each is independently testable with a stated interface and dependency set.

### `utils/mapCache.ts` — add a subpath-boundary index (Task 1)

`cellVerts` concatenates every sub-ring's points with no boundary marker (its own
header comment says a renderer must add its own index first). Emit that index.
The data already exists in the `perCell: CapturedSubpath[][]` array built during
capture — expose it as flat offsets.

- **Add to `MapGeometryCache`:** `cellSubOffsets: Uint32Array` — flat, giving
  each cell's sub-ring start indices into `cellVerts` (structure to be finalized
  in the plan; e.g. per-cell run of sub-ring vertex-start offsets plus the cell's
  own `cellOffsets` bound as the terminator).
- **Contract:** lands under the existing **parity test** (cached vertices match a
  fresh `d3.geoPath` trace to sub-pixel precision) — extend that test to assert
  sub-ring boundaries line up with the capture; do **not** add a world-generating
  test.

### `utils/tessellate.ts` — pure triangulation

- **Input:** `cellVerts`, `cellOffsets`, `cellSubOffsets`.
- **Output:** `{ positions: Float32Array, cellTriRange: Uint32Array }`.
- **Method:** earcut per sub-ring. Uses the `earcut` library (small, standard;
  confirm/add dep). A cell that legitimately emits no sub-rings (pole-enclosing
  under a cylindrical projection — d3's clip stream drops it) contributes zero
  triangles; its `cellTriRange` is empty. **The parity test must not assert every
  cell produces triangles.**
- **Depends on:** `earcut` only. No world, no projection, no DOM.

### `utils/mapFillColorBuffer.ts` — per-vertex colors

- **Input:** the per-cell color source (`buildCellColorCache` output from
  `mapColorCache.ts`), optional `shadeMap`, and `cellTriRange`.
- **Output:** `Float32Array` RGB per vertex (or `Uint8Array` normalized —
  decide in plan), aligned to `positions`.
- **Recolor:** a `writeCellColors(buffer, cellIndices, …)` helper rewrites only
  the given cells' ranges for paint strokes.
- **Depends on:** the color cache and `cellTriRange`. Pure.

### `components/map2d/GLFillSurface.tsx`

Owns the Three renderer, camera, geometry, material, and the imperative draw +
disposal. Props: geometry buffers, color buffer, and the live `(scale, offset,
dpr, size)` transform inputs. Exposes an imperative recolor/redraw handle to
Map2D. No app state inside.

### `components/map2d/VectorOverlay.tsx`

Owns the transparent Canvas2D overlay and its per-frame redraw. Props: the same
transform inputs plus the world/overlay data the existing passes consume. Calls
the existing draw functions under the composed transform.

### `components/Map2D.tsx` — wiring

- Selects **GL path vs blit fallback** from: style === Default, projection ≠
  Dymaxion, WebGL available, context not lost.
- Drives both layers' transforms from the existing `scale`/`offset`/`size`/dpr.
- Retains the blit renderer intact for every other case.
- The `INTERACTION_DPR` / `qualityDpr` / `MAX_SHARP_*` raster-sharpness machinery
  is **bypassed on the GL path** (resolution is now zoom-independent) but stays
  for the blit path.

## Data flow

- **Build once** (world or projection or size changes): `mapCache` (+ new
  subOffsets) → `tessellate` → positions + cellTriRange; `mapFillColorBuffer` →
  color buffer. Upload to GPU. Measure build time at 200k — if it stalls the main
  thread, move `tessellate` into the worker (it is pure and serializable).
- **Per frame** (pan/zoom): matrix uniform update + one GL draw; overlay clear +
  redraw under the composed transform. No rebuild.
- **Recolor** (view mode / season / seaLevel / hillshade): rewrite color buffer,
  re-upload, redraw. Topology untouched.
- **Paint stroke:** rewrite affected cells' color slices only; re-upload with a
  bounded update range; redraw. Honors the "reuse geometry, mutate colors"
  invariant.
- **Zoom-settle:** rebuild boundary LOD paths at `tol = 0.5 / scale`; overlay
  picks them up on its next redraw.
- **Pick:** unchanged — `mapPick` quadtree over projected cell centres; the
  screen→world inverse of the pan/zoom affine is applied to the query point, as
  today.

## Correctness invariants

1. **Two canvases agree to the pixel** on CSS size, backing-store size (× dpr),
   and transform. Any disagreement makes fills and coastlines halo apart — which
   reads as the same "weird corners" artifact. Size them from one source.
2. **Mirror carried consistently:** GL matrix and overlay transform both apply
   the horizontal flip; text/glyphs stay compensated exactly as today (§mirror
   trap). Wrong here = backwards labels no unit test catches — verify in-browser.
3. **Lint ratchet is 30/30 with zero headroom** (HANDOFF S30b). New Three/GL code
   is the classic source of `any`. Type it properly from the first commit or CI
   fails.
4. **Disposal:** every GPU resource disposed on unmount and on world/projection
   change (WorldViewer pattern).
5. **Fallback parity:** the blit path remains byte-for-byte the current behavior;
   Phase 2 must not change what non-Default / Dymaxion / WebGL-absent users see.

## Testing strategy

- **`tessellate` (unit):** every triangle's vertices lie within their cell's
  projected bbox; total triangle area per cell ≈ the cell polygon area (within
  epsilon); a synthetic non-convex ring triangulates without outside triangles;
  a zero-sub-ring cell yields an empty range. No world generation.
- **`mapCache` subOffsets (unit):** extend the existing **parity** test — sub-ring
  boundaries match the capture; existing sub-pixel vertex parity still holds.
- **`mapFillColorBuffer` (unit):** per-vertex colors match the per-cell cache;
  `writeCellColors` touches only the named cells' ranges; `shadeMap` multiply
  matches the Canvas2D fill path.
- **Selection logic (unit):** Default+Mercator+WebGL → GL; Parchment / Dymaxion /
  no-WebGL / context-lost → blit.
- **Manual browser verification (required):** crisp fills, corners, coastlines,
  and labels at 1× and 6×; correct (non-mirrored) text; view-mode switch recolors
  correctly; paint stroke updates the right cells; context-loss falls back; pick
  hits the right cell at zoom. Rendering is verified in the browser per repo
  convention.
- **Test-running discipline:** touched-file runs by default; full suite only
  before merge, with Matt's speed choice; never alongside browser automation
  (per `CLAUDE.md`).

## Risks

- **200k triangulation build time / GPU memory** on the 16GB / 7-GPU-core M1.
  ~200k cells × ~4 triangles = ~800k tris ≈ 2.4M verts; position ~19MB + color.
  Mitigate: measure first; move `tessellate` to the worker if it stalls; consider
  a 16-bit position encoding only if memory bites.
- **Overlay per-frame cost at high coastline detail** (200k → long chains).
  Mitigate: boundary LOD already reduces vertices; if a live drag still drops
  frames, redraw the overlay at a lower cadence during drag and sharp on settle
  (GL fills stay live-crisp regardless). Note as a tactic, not v1 default.
- **Color-buffer rebuild on view-mode switch** (2.4M verts). Should be a fast
  typed loop; measure and, if needed, cache per-mode buffers.
- **Three version drift** in the buffer-update API — pin to installed version via
  Context7.

## Task ordering (for the plan)

1. `mapCache` subOffsets + parity test.
2. `utils/tessellate.ts` + unit tests (earcut dep).
3. `utils/mapFillColorBuffer.ts` + unit tests.
4. `GLFillSurface` component (build, draw, dispose) — static, no interaction yet.
5. `VectorOverlay` component + transform + boundary LOD.
6. Map2D wiring: selection, transform drive, fallback (init + context-loss).
7. Recolor + paint partial-update path.
8. Perf pass at 200k (build-time / worker decision / cadence).
9. Browser verification + docs (ROADMAP F3 Phase 2 → done, HANDOFF, map-styles
   note on the two-layer path).

## Out of scope (explicit)

Parchment GPU port (future 2b), Dymaxion changes, 3D globe changes, any change to
the generation pipeline, raw-WebGL/regl rewrite.
