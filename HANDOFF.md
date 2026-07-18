# RealmGenesis 3D - Agent Handoff

This file is the quick state transfer for the next session. Read
`ARCHITECTURE.md` for deeper technical context and `AGENTS.md` for repo
workflow/style rules.

---

## Project Snapshot

RealmGenesis 3D is a browser-only procedural fantasy world generator built with
React 19, Three.js/R3F, Canvas2D, d3 geo tooling, and Gemini BYOK lore support.
There is no backend.

- Dev server: `npm run dev` on port 3000
- Quality gates (all enforced in CI): `npm run lint`, `npm run typecheck`,
  `npm test`, `npm run build`
- Vitest engine suite in `tests/`; no formatter is configured
- Deployment target: Netlify static SPA

---

## Work Completed In This Session

### Audit Hardening Batch (remaining AUDIT.md items)

- **CI added**: `.github/workflows/ci.yml` runs lint + typecheck + test + build
  on pushes to main and all PRs.
- **TypeScript strict mode enabled** (`"strict": true`); zero errors. Added
  `@types/d3`, `@types/react`, `@types/react-dom` and a local `vendor.d.ts`
  shim for `d3-geo-projection` / `d3-geo-voronoi` (no registry types).
  `WorldData.geoJson` is now a typed `GeoJsonCollection` (d3-compatible).
- **Vitest suite added** (`tests/`, 35 tests, ~2 s): RNG determinism, biome
  classification table, generation determinism/structure/abort/progress,
  param liveness (fails if any tunable param stops affecting output), paint
  utils, and import validation. `npm test`.
- **Engine fixes surfaced by the tests**:
  - `civSizeVariance` reworked from expansion budgets (which never bound on
    small maps) to per-faction competitive movement-cost scaling — effective
    at any map resolution.
  - Negative cell populations fixed (high-elevation suitability went negative
    and silently deflated faction totals).
  - `capitalSpacing` threshold was resolution-dependent and almost never
    fired; now a scale-independent squared-chord minimum
    (`spacing^2 * 4 / numFactions`).
- **Tailwind moved from CDN to the build pipeline** (tailwindcss v3 + PostCSS,
  `index.css`, purged ~23 kB output). CSP no longer allows any external script
  host or `unsafe-eval`.
- **`npm audit fix` applied**: 0 vulnerabilities (was 1 critical / 5 high).
- **Lore errors surfaced properly**: `generateWorldLore` throws instead of
  returning sentinel "Error World" lore, and validates the model's JSON field
  by field before mutating civData. `@google/genai` and `GLTFExporter` are now
  dynamically imported (smaller main bundle).
- **Import hardening**: `validateCivData` shape-checks imported civData
  (malformed metadata degrades to terrain-only load); points input capped at
  200k (UI + validation) to match the slider and avoid main-thread freezes.
- **`detailLevel` implemented** as the FBM octave count with a "Detail
  Octaves" Geo slider (default 3 = historical hardcoded value, so default
  worlds are unchanged).
- **Misc**: progress bar reaches exactly 100% (7 ticks / TOTAL_STAGES 7);
  Map2D caption reflects Mercator vs Dymaxion; `ExportResolution` type matches
  the UI (2K/4K/8K); shared Dymaxion geo helpers consolidated into
  `utils/geo.ts`; unused imports removed; lint script has a
  `--max-warnings 30` ratchet (remaining warnings are documented-intentional
  R3F `any`s and hook-dep patterns).

Validation this session: `npm run lint`, `npm run typecheck`, `npm test`
(35/35), `npm run build`, plus a Playwright pass with the bundled Tailwind:
fully styled UI, globe + Dymaxion rendering, new Geo sliders present,
Dymaxion caption correct, no console errors.

---

## Previous Session



### New Features (from AUDIT.md open questions)

- **`civSizeVariance` is now implemented**: each faction draws a per-faction
  expansion budget from `civRng` (base 200, spread `1 ± variance`, clamped to
  0.25x–2x). The Dijkstra frontier stops once a faction's cost exceeds its
  budget. The "Country Size Variance" slider in the Civ tab now does something.
- **`population` view mode implemented**: log-scaled heat gradient on land
  (dark blue → green → bright yellow), uninhabited land dark grey, ocean navy.
  New "Population" View Layer button in the Sys tab.
- **`province` view mode implemented**: faction base color with amplified
  per-province shade variation (variant strength 1.8). New "Provinces" button.
- **`Pangea` land style** now has a real engine branch in `worldGen.ts`
  (`landChance 0.6, landLevel 0.25, oceanLevel -0.45`), on top of the mask the
  preset already set.
- **Cell Jitter slider added to the Geo tab** (was documented but missing).
- **Build-time Gemini key fixed**: `services/gemini.ts` now reads
  `process.env.GEMINI_API_KEY`, matching the `vite.config.ts` define and the
  README's `.env.local` instructions (previously read `API_KEY`, which was
  always undefined). Verified via sentinel-key build.

### 3D Rendering Optimizations (`WorldViewer.tsx`)

- **World mesh geometry is now allocated once per world structure and refilled
  in place** on paint strokes / view changes (`useLayoutEffect`), instead of
  rebuilding + leaking a new ~45k-triangle `BufferGeometry` on every stroke
  event. Attributes use `DynamicDrawUsage`; bounding sphere is fixed (r = 1.1).
- **Dropped `computeVertexNormals` and the normal attribute** from the world
  mesh, selection overlay, and curved labels — the unlit basic material ignores
  normals and the standard material's `flatShading` derives them in-shader.
- **Every `useMemo` geometry now has a disposal effect** (world mesh, rivers,
  faction borders, brush ring, selection overlay, highlight outline, lat/long
  grid, Dymaxion overlay, curved label geometry). GPU memory no longer grows
  across regenerations and paint sessions.
- **Rivers keyed on `world.rivers`** (stable across strokes) so painting never
  re-runs CatmullRom smoothing; **curved faction labels keyed on position
  values** so unchanged centroids don't rebuild patches; the inline cell
  highlight was converted to a memoized `CellHighlightOutline` component.
- **`CityMarkers` no longer allocates Vector3s per instance per update.**
- **Map2D Dymaxion pick buffer keyed on `world.cells`** instead of world
  identity — skips a full-canvas per-pixel reprojection per stroke event.
- **MiniMap redraw debounced (150 ms)** and now passes live faction colors.

### Color Consistency

- New `buildFactionColorMap(civData)` helper in `colors.ts`; `MiniMap`,
  `export.ts` (both render loops), and `exportGLB.ts` now pass live faction
  colors, so user-edited faction colors appear in the minimap, PNG exports,
  and GLB vertex colors (closes the gap invariant 21 warned about).

### Correctness / Type Fixes

- `npx tsc --noEmit` is now clean (was 6 errors): `Controls.tsx` `setParams`
  prop typed as `React.Dispatch<React.SetStateAction<WorldParams>>`; dead
  `dymaxion` switch branch removed from `exportMap` (unreachable after the
  early return); dymaxion projection nodes typed with `children`; dead
  comparison removed in `Inspector.tsx`.
- The `G` keyboard shortcut no longer generates from stale params (and no
  longer reverts recent slider edits): the keydown effect re-registers on
  `params`/lock changes.
- Auto-update now also triggers on `mountainHeight`, `oceanDepth`, and
  `cellJitter` changes (previously only the first change fired, via the
  `landStyle → 'Custom'` side effect).

---

## Files Most Recently Touched

| File | Notes |
|------|-------|
| `components/WorldViewer.tsx` | Geometry reuse + in-place refill, disposal everywhere, no normals, rivers/labels keyed on stable refs, `CellHighlightOutline` component. |
| `utils/colors.ts` | `population` + `province` view modes, `buildFactionColorMap`, `getProvinceVariant` strength param. |
| `utils/worldGen.ts` | `civSizeVariance` expansion budgets; `Pangea` land style branch. |
| `components/Controls.tsx` | Provinces/Population view buttons, Cell Jitter slider, `setParams` Dispatch type, keyboard shortcut staleness fix, auto-update deps. |
| `components/Map2D.tsx` | Pick buffer keyed on `world.cells` (skips per-stroke rebuild). |
| `components/MiniMap.tsx` | Debounced redraw, live faction colors. |
| `utils/export.ts` | Live faction colors in both render loops; dead dymaxion branch removed. |
| `utils/exportGLB.ts` | Live faction colors in GLB vertex colors. |
| `services/gemini.ts` | Reads `process.env.GEMINI_API_KEY`. |
| `components/Inspector.tsx` | Dead comparison removed. |
| `utils/dymaxion.ts` | Projection tree nodes typed (`children`). |
| `README.md`, `ARCHITECTURE.md` | New view modes, civSizeVariance behavior, env var name, invariants 26–28 (geometry reuse, disposal, pick-buffer keying). |

---

## Validation

Last successful checks:

- `npm run build`
- `npm run lint` (zero errors; warnings only)
- `npx tsc --noEmit` (zero errors — new since this session)
- Headless-browser pass (Playwright + SwiftShader): generation completes with
  no console errors; Biomes/Provinces/Population layers render on the 3D globe,
  2D Mercator, and minimap; Cell Jitter slider present. (Tailwind CDN is
  blocked in the sandbox, so layout was stubbed — visual styling not covered.)
- `.env.local` sentinel key verified to land in the built bundle.

---

## Important Invariants

- All app-level state remains in `App.tsx` and is prop-drilled.
- 3D inspection should stay click-only unless `inspectMode === 'hover'`.
- 3D paint raycasting should only run during active strokes to avoid idle hover
  lag.
- Dymaxion picking must use the same raster pipeline as the visible Dymaxion map.
- Political brush undo must preserve `height`, `biome`, `regionId`, and
  `provinceId`.
- Political brush should not create provinces; it assigns cells to an existing
  selected-faction province or clears ownership when using the eraser.
- `getCellColor(..., factionColors)` must receive the live `factionColors` map
  on render paths that display political ownership.
- 3D faction labels are curved textured meshes, not HTML overlays or billboard
  sprites.

---

## Feature Roadmap (shortlist)

See `ROADMAP.md` for the full detail. Themes, by leverage:

- **A. Cartographic presentation** — map labels/typography, offline Markov
  namebases (Gemini becomes optional), fantasy map styles, hillshading +
  contours, great-circle ruler + scale bar, geodesic hex grid.
- **B. Physical geography** — lakes as first-class entities (the Priority-
  Flood fill in `generateRivers` already computes them and throws them away),
  river/lake editing, auto-detected + named ranges/seas/deserts/islands.
- **C. Worldbuilding depth** — cultures layer, religions, roads + sea trade
  routes (reuse civ Dijkstra costs), markers/POIs, editor completeness,
  diplomacy later.
- **D. Planet-scale simulation** — seasonal cycle from axial tilt, ocean
  currents feeding climate, ice caps, regional submap generation, planetary
  parameters. These are the sphere-native differentiators vs. flat-map tools.
- **E. Interoperability** — SVG export, GeoJSON export, Azgaar `.map` import
  (stretch).

**Recommended first milestone — "Map identity"**: labels + offline namebases
+ lakes + named geographic features (A1, A2, B1, B3). Turns the simulation
viewer into something that reads as a map; consumes data the pipeline already
computes. A phased implementation design is summarized in `ROADMAP.md`.

---

## Potential Next Tasks

- Manual browser regression pass across 3D, 2D Mercator, and 2D Dymaxion with
  real styling (the sandbox pass could not load Tailwind).
- Remaining items from `AUDIT.md`: CI workflow (lint + typecheck + build),
  Vitest suite over the pure engine, staged TypeScript strict mode,
  `npm audit fix`, Tailwind via build pipeline instead of CDN, code splitting.
- Optional next-level paint optimization: per-cell dirty ranges so strokes
  refill only affected vertices instead of all cells (current in-place refill
  is already allocation-free).
- Tune 3D curved faction label sizing/placement if very long faction names are
  introduced.
- Broaden editor features later: merge/split factions, province management,
  bulk rename, town/capital relocation.

---

## Workflow Notes

- Do not push without explicit user request.
- Git commit messages are strongly recommended to follow the 50/72 rule:
  subject line at or under 50 characters, body wrapped near 72 characters.
- Keep commits focused and imperative.
- The user normally runs the dev server; only start it when useful for manual
  verification or when asked.
