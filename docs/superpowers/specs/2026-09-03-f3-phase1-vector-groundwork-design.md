# F3 Phase 1 — cached-Canvas2D map + Phase-2 vector groundwork

Status: DESIGN (approved in chat 2026-09-03, S30). Advisor-reviewed.
Roadmap: F3 (True 2D vector map). Follows the phased decision below.

## 1. Context and goal

The 2D map (`components/Map2D.tsx`) is Canvas2D + d3-geo. A pan or zoom at 30k
cells blocks the main thread ~2.2 s at **each end of the gesture** — the worst
number in the project (`docs/performance-findings.md` §2). The dragging itself is
free (a bitmap blit); the cost is the two full offscreen redraws triggered by the
DPR switch (interaction drops DPR to 1, settle raises it back), and each redraw:

- **re-projects every cell polygon vertex from scratch through d3's `geoPath`,
  several passes per cell** (land, hillshade, ocean fill, hatch, coastline) —
  there is no cached screen geometry; and
- **recomputes per-cell colour** (`getCellColor` + `seasonalTemperatureDelta` +
  a `THREE.Color` allocation) every redraw.

`computeShadeMap` is already memoised and is NOT the culprit.

**Decision (S30).** Industry map rendering is WebGL vector (Google 2019, Apple,
Bing/Azure, and OSM.org since July 2025); SVG-DOM is dead past ~10k elements, so
a one-node-per-cell DOM map at the 200k cap is a non-starter. The chosen path is
**phased**:

- **Phase 1 (this spec):** stay on Canvas2D, but project and colour **once** and
  cache, killing the stall, AND build the reusable geometry the WebGL phase will
  consume (edge-chained boundaries, a simplification primitive, a pick index).
- **Phase 2 (separate spec, later):** a WebGL vector surface on the Three.js
  already in the stack, eating Phase 1's cached vertex buffers.

Matt authorised the full-groundwork Phase 1 (the largest option): the two caches
PLUS edge-chaining, a simplification/LOD scaffold, and the picking quadtree.

## 2. Confirmed mechanics (grounding — verified in code, not assumed)

These decide the cache design; each was read, not inferred.

- **Projection is CSS-px.** `buildProjection(...).fitSize([size.width,
  size.height], Sphere)` (`utils/projections.ts:27,34`), memoised on
  `[size, projectionType]` (`Map2D.tsx:405-410`) — **independent of `scale`/zoom**.
- **DPR and the horizontal mirror both live in the CONTEXT transform, not the
  path coords.** The offscreen effect does
  `ctx.setTransform(renderDpr,0,0,renderDpr,0,0); ctx.translate(size.width,0);
  ctx.scale(-1,1)` (`Map2D.tsx:471-474`), where `renderDpr = qualityDpr`. So d3
  emits un-flipped CSS-px coordinates and the context applies both DPR scaling and
  the mirror.
- **Zoom is a blit-time transform + a DPR bump**, not a reprojection. The visible
  blit applies `setTransform(displayDpr*scale, …, offset)` (`Map2D.tsx:1056`);
  `qualityDpr` settles to `baseDpr * min(scale,2.5)` capped 3 (`Map2D.tsx:385-391`),
  so zooming re-enters the offscreen redraw at higher DPR to regain sharpness —
  but the base CSS-px projection is unchanged.

**Consequences that shape the design:**
1. The cache stores **un-flipped CSS-px geometry**. It is filled through the
   existing flipped + DPR-scaled context, so the mirror and DPR are handled for
   free and never baked into the cache.
2. The geometry cache is **independent of zoom and DPR** — its key is
   `(projectionType, width, height, world)` with no zoom dimension.
3. **Dymaxion is a different pipeline** (`Map2D.tsx:483-...`): it rasterises an
   equirectangular source per-pixel with its own flip and its own pick buffer.
   Phase 1 caching targets the d3-projection views only; Dymaxion is out of scope
   and untouched (see §6).

## 3. The five units

### Unit 1 — Screen-geometry cache (`utils/mapCache.ts`, new)

Project every cell ring to CSS-px **once** per `(projectionType, width, height,
world)`. Store:

- `verts: Float32Array` — projected ring vertices, ragged, addressed by
  `offsets: Uint32Array` (cell i's ring is `verts[offsets[i]*2 .. offsets[i+1]*2]`).
  This is the format Phase 2 WebGL consumes directly.
- `paths: Path2D[]` — one `Path2D` per cell, built from `verts`. Canvas2D fills
  these instead of re-running `geoPath(feature)`. Path2D construction is cheap;
  the expense being removed is the per-vertex projection.

Rebuilt only when the key changes. **NOT** rebuilt on pan/zoom (blit transform) or
DPR (context scale). Coordinates are un-flipped CSS-px (see §2). Boundary polylines
(coastlines/borders/rivers/routes, from Units 3/4) are projected into the same
cache in the same pass, as vertex arrays + `Path2D`.

Interface (pure, no React/DOM beyond `Path2D`):
```
buildMapGeometryCache(
  world, projection, width, height,
  chained: ChainedBoundaries,       // from Unit 3, in SPHERE space
): MapGeometryCache
// { cellVerts, cellOffsets, cellPaths, coastPaths, borderPaths, riverPaths, ... }
```

**Pipeline order (resolves the tolerance-units question):** Unit 3 produces
chains in **sphere** space → this build **projects** them to CSS-px → Unit 4
DP-simplifies the **projected** polylines (so the tolerance is genuinely in CSS
pixels) → store as vertex arrays + `Path2D`. Cell rings are projected but NOT
simplified in Phase 1 (fill LOD is Phase 2).

### Unit 2 — Colour cache

Per-cell fill colour as a `Uint32Array` (packed RGBA), computed once and read by
the fill passes. Key: `(viewMode, season, seaLevel, factionColors, cultureColors,
religionColors)` identities — recomputed only when one changes, NOT on pan/zoom.
`shadeMap` multiply is folded in at build time (it is already memoised upstream).

Interface:
```
buildCellColorCache(world, viewMode, colorCtx, shadeMap): Uint32Array
```
Passes read `colors[i]` instead of `getCellColor(cell,…) *
seasonalTemperatureDelta(...)` + `THREE.Color` alloc per redraw.

**Substrate primitive (the one sub-fork, decided).** Add
`fillCells(indices: Uint32Array | number[], colors: Uint32Array)` to the
`Substrate` interface, backed by the caches, because "fill all land / ocean /
hillshade cells" are the hot loops and a batched call keeps one code path.
- `Canvas2DSubstrate.fillCells`: loop indices, `ctx.fillStyle` from the packed
  colour, `ctx.fill(cache.cellPaths[i])`. Uses the geometry cache; no reprojection.
- `SvgSubstrate.fillCells`: **loops per-cell `<path>` exactly as its current
  `fillFeature` loop does — no aggregation, no change to export output.** The
  batching benefits only Canvas2D (and future WebGL). SVG-side aggregation
  (dissolving same-colour cells into multipolygons) is explicitly OUT of Phase 1
  and noted for later — the "parity" here is *behavioural sameness of output*, not
  a shared fast path. Precedent: `hatchFeatures` is already plural for cost while
  each substrate implements it differently.

Rejected alternative: making `fillFeature` cache-aware via a feature→index lookup.
Features carry no clean index; a batch primitive is cleaner and is the Phase-2 GPU
seam.

### Unit 3 — Edge-chaining (`utils/boundaries.ts` extension)

Today coastlines and borders are an **unordered soup of disjoint 2-point edges**
(`computeBoundarySegments`, `boundaries.ts:27-54`); rings are never assembled, and
matching is by **distance-threshold vertex comparison** (approximate).

Add `chainSegments(segments): Polyline[]` that assembles contiguous chains:
- **Quantise endpoints** to a fixed grid (snap to ~1e-6 on the unit sphere) so
  "shared endpoint" becomes an exact hash-map key rather than an approximate
  distance test — this is what makes the walk robust at T-junctions.
- Build an endpoint→segments adjacency map; walk each unused segment into a chain;
  at a branch point (endpoint shared by >2 segments — a T-junction where 3 cells
  meet) pick the continuation **deterministically** (e.g. smallest segment id) and
  leave the others to seed their own chains.
- A chain **closes** into a ring when its walk returns to its start endpoint.

Apply to: coastlines (`land≠land`), faction borders, province borders, and **lake
outlines** (boundary of each lake-cell cluster — new; lakes currently have no
outline geometry, only `LakeData.cellIds`).

**De-duplication:** `Map2D.tsx` has its own private `getFactionBorders`
(`Map2D.tsx:52-80`), a byte-equivalent second copy of the border scan. Phase 1
routes Map2D through the shared `utils/boundaries.ts` functions and deletes the
private copy (the `utils/ must not import components/` rule that motivated the copy
is satisfied because the shared code lives in `utils/`).

**Rivers and routes are already contiguous polylines** from generation
(`world.rivers: Point[][]`, `world.routes[].path`, per the geometry inventory) —
they bypass chaining and go straight to projection + simplification in Unit 1.
Only coastlines, borders, and lake outlines need Unit 3.

Output: sphere `Point3` polylines, consumed by Unit 1 (projected) then Unit 4
(simplified).

### Unit 4 — Simplification / LOD scaffold (`utils/simplify.ts`, new)

A pure Douglas–Peucker `simplifyPolyline(points, tolerance): Point[]`.

**Applied in Phase 1 to polyline features only** (coastlines, borders, rivers,
lake outlines), **once at cache-build time, on the PROJECTED CSS-px polylines**
(inside Unit 1's build, after projection — see the pipeline note in Unit 1), so
the fixed tolerance ≈ 0.5 CSS-px is a genuine sub-pixel quantity. Near-lossless at
base resolution. It is NOT zoom-keyed and does NOT re-key Unit 1's cache —
resolving the Unit1↔Unit4 conflict the advisor flagged: the base CSS-px geometry
is zoom-independent, so there is one cached copy and no `zoomBucket`.

**Scaffold, honestly scoped:** Phase 1 delivers the DP primitive + the
integration seam (the cache-build call site) + conservative defaults. It does
**not** implement cell-fill LOD / aggregation for the 200k case (dissolving
same-value cells) — that is Phase 2 work, deliberately deferred (§6). "Scaffold"
means the mechanism and the seam exist; it does not mean 200k is solved.

### Unit 5 — Picking quadtree (`utils/mapPick.ts`, new)

Replace the O(cells) nearest-cell scan (`getNearestCellId`, `Map2D.tsx:31-50`,
~20M dot products per click/hover/paint-move at 200k) with a `d3-quadtree`
(already in the d3 dependency tree) over **projected cell centres**, rebuilt on
`(projection, size, world)`. `getCellIdAtMapPoint` uses `quadtree.find(x,y)` for
d3-projection views.

**Dymaxion keeps its RGB pick-buffer path untouched** — it must mirror the visible
rasterisation (hard invariant, `docs/invariants.md:145`, `CLAUDE.md`). The
quadtree replaces only the d3-projection scan.

**Correctness caveat (advisor):** nearest-*projected*-centre ≠ nearest-*geodesic*
at the antimeridian and high Mercator latitudes. The parity test MUST sample
deliberately near the seam and at high latitude, not uniformly (see §5).

## 4. Integration (`components/Map2D.tsx`) — kept in-house

The offscreen redraw effect (`Map2D.tsx:462-937`) is rewired to:
1. On `(projection, size, world)` change: build Unit 3 chains → Unit 4 simplify →
   Unit 1 geometry cache; build Unit 5 quadtree.
2. On `(viewMode, season, colormaps, seaLevel, shadeMap)` change: build Unit 2
   colour cache.
3. Per redraw: passes call `sub.fillCells(indices, colorCache)` against the cached
   `Path2D`s; boundary passes stroke the cached polylines. No `geoPath(feature)`
   per redraw.

Cache lifecycle lives in refs + `useMemo`/`useEffect` keyed as above, mirroring the
existing memo pattern (`shadeMap`, `coastlines`, etc. already work this way).

This unit holds every invariant and is serial through one large file — not
delegated.

## 5. Testing

- **Geometry-cache parity** (`tests/mapCache.test.ts`): cached fill vs a fresh
  `geoPath`-projected fill for the same world/projection/size — sub-pixel max
  error, per cell. Include a non-identity projection and a non-square size.
- **Colour-cache parity** (`tests/mapColorCache.test.ts`): `colorCache[i]` equals
  `getCellColor(...) * shadeMap[i]` packed, across viewModes and a non-neutral
  season.
- **Edge-chaining** (`tests/boundaries.test.ts`, extend): on a generated world,
  **every coastline chain closes into a ring — count of unclosed coastline chains
  == 0** (hard assertion, not a spot check); total segment count is conserved
  (sum of chain edges == input segment count); borders may be open but must not
  drop segments. T-junction case covered with a synthetic fixture.
- **Simplification** (`tests/simplify.test.ts`): DP correctness (endpoints
  preserved, monotone vertex reduction with tolerance, collinear points removed);
  max deviation ≤ tolerance.
- **Pick parity** (`tests/mapPick.test.ts`): `quadtree.find` returns the same cell
  as the O(cells) scan for sample points **biased to the antimeridian seam and
  high latitudes** (not uniform) — document any divergence as a known edge rather
  than asserting false equality.
- **Perf**: `node scripts/renderMap2DPerf.mjs --points=30000 --dpr=2` before/after;
  record the settle-redraw ms and the settled pixel hash (unchanged output).
- Gates: typecheck 0, lint within ratchet, touched-file test runs (not the full
  suite unless Map2D/colors/types change — then ask per CLAUDE.md).

## 6. Out of scope (Phase 2 and beyond)

- **WebGL vector surface** — the whole point of Phase 2; consumes Unit 1's
  `verts`/`offsets` buffers.
- **Cell-fill LOD / aggregation for 200k** — dissolving same-value cells into
  multipolygons; Unit 4 scaffolds the primitive only.
- **SVG-side `fillCells` aggregation** — SVG export stays one-path-per-cell.
- **Dymaxion caching** — different per-pixel pipeline; untouched.
- **Native/GPU picking** — Phase 1 keeps Canvas picking (quadtree); GPU pick pass
  is Phase 2.

## 7. Files

New: `utils/mapCache.ts`, `utils/simplify.ts`, `utils/mapPick.ts`, and their tests.
Extend: `utils/boundaries.ts` (chaining, lake outlines), `utils/mapStyle/substrate.ts`
+ `substrateCanvas.ts` + `substrateSvg.ts` (`fillCells`), `utils/mapStyle/passes.ts`
(call `fillCells` + colour cache). Integrate: `components/Map2D.tsx`.

## 8. Delegation plan (for the implementation plan)

Pure, independently-testable units are delegatable to Sonnet with written briefs
(NOT "go investigate"): Unit 3 (edge-chaining), Unit 4 (simplification), Unit 5
(pick quadtree), Unit 2 (colour cache) — mutually independent, parallelisable.
Unit 1's pure module is delegatable; the Substrate `fillCells` change needs close
review; the Map2D integration is kept in-house (invariant-dense, serial). Order:
utils in parallel → substrate primitive → Map2D integration last.
