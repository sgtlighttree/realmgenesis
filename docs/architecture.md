# Architecture

RealmGenesis 3D is a browser-only procedural fantasy world engine. It simulates
planetary geography, climate, biomes, and civilizations on a sphere and renders
them as an interactive 3D globe, a 2D map, or a Dymaxion projection. It runs
entirely client-side — no backend — with seeded RNG for reproducibility.

Stack: React 19 + TypeScript (strict) + Vite 6; Three.js via
`@react-three/fiber`/`drei` for 3D; Canvas2D + d3 (`d3-geo-voronoi`,
`d3-geo-projection`) for 2D; `@google/genai` for optional lore. Tailwind CSS
(build-time, purged). Deployed as a static SPA on Netlify.

## Entry routing

`index.tsx` picks a root component from the `?shell` query param:

| `?shell=` | Root | Purpose |
|-----------|------|---------|
| *(none)* / `1` | `components/shell/ShellApp.tsx` | **Default.** The F1 redesign on real data. |
| `classic` | `App.tsx` | Legacy pre-redesign fork, kept until fully retired. |
| `stub` | `components/shell/DesignShell.tsx` | F1 layout prototype with stub panels. |

> The default route is the **shell**, not `App.tsx`. `App.tsx` is legacy. Any
> doc or instinct that says "start at `App.tsx`" is describing the classic
> route.

## State ownership

All application state lives in one hook: **`hooks/useWorldEngine.ts`**. It holds
`params`, `world`, view/display/inspect modes, edit state, overlay toggles
(grid, rivers, routes, hillshade, contours, labels), ruler and marker state, and
every handler (generate/cancel/load, paint/undo, civ/province recalc, faction
edit, lore). It returns them as one object; `ShellApp` consumes the hook and
prop-drills into the components. No Context, Redux, or Zustand — trace any value
to `useWorldEngine`.

`DEFAULT_PARAMS` is defined in `utils/defaultParams.ts` (a types-only module, so it
imports cleanly outside the worker chain). `useWorldEngine.ts` seeds its `params`
state from it, and `tests/helpers.ts` `makeParams` spreads it (overriding ~7 keys for
test speed) — one source of truth, so the two cannot silently diverge.

## The generation boundary — a Web Worker

Generation does **not** run on the main thread. `useWorldEngine` calls
`generateWorldInWorker()` (`utils/worldGenClient.ts`), which posts the params to
`workers/worldGen.worker.ts`; the worker runs `generateWorld` and returns a
Structure-of-Arrays serialization (`utils/worldTransfer.ts`) that the main
thread rehydrates into `Cell[]`.

```
ShellApp
  └── useWorldEngine.handleGenerate()
        └── generateWorldInWorker(params, onLog, signal, onProgress)   [worldGenClient.ts]
              └── postMessage → worldGen.worker.ts
                    └── generateWorld(params, …)                       [worldGen.ts]
              └── worker posts serialized world (SoA)                  [worldTransfer.ts]
        └── rehydrate → Cell[]  → setWorld(world)
              └── re-render WorldViewer / Map2D
```

This buys **responsiveness, not throughput** — the worker path is ~20% slower in
wall-clock (serialize + transfer + rehydrate) but keeps the UI at 60fps during a
multi-second generation. Abort is `worker.terminate()` (one worker per
generation); a running synchronous loop can't process an abort *message*. See
[invariants.md](invariants.md) §1 for the consequences.

Civilization and province recalculation run on the **main thread** (they're fast
and must not go through the rehydration path, which would replace `world.cells`
and force a geometry rebuild — see [invariants.md](invariants.md) §1).

## Module map

Purposes, not line numbers (line anchors rot fastest). See each doc for depth.

### Engine (`utils/`)
| Module | Purpose |
|--------|---------|
| `worldGen.ts` | The generation pipeline — `generateWorld`, `recalculateCivs`, `recalculateProvinces`, `determineBiome`, erosion, rivers/lakes. See [generation-pipeline.md](generation-pipeline.md). |
| `tectonicsV3.ts` | The V3 terrain model: crust fields, Euler-pole kinematics, Dijkstra plate growth, GDH1 bathymetry, microplates. See [tectonics-v3.md](tectonics-v3.md). |
| `crust.ts` | Independent crust-field seeding (`seedCrustField`) — continental vs oceanic, decoupled from plates. |
| `spherical.ts` | Sphere math: quaternion rotation, Euler poles, vector ops. |
| `rng.ts` | Mulberry32 seeded PRNG + 3D Simplex noise; isolated side-streams per subsystem. |
| `pathfinding.ts` | `MinHeap` + water-cost helpers, shared by rivers, civ Dijkstra, and routes. |
| `namegen.ts` | Offline seeded char-level Markov namebases (feature A2). |
| `features.ts` | Auto-detected + named geographic features: ranges, peaks, deserts, oceans, islands (B3). |
| `routes.ts` | Roads & sea trade routes over the terrain-cost model (`computeRoutes`, C3). |
| `civEdit.ts` | Faction/province editing: merge, rename, relocate capital (C5). |
| `measure.ts` | Great-circle geodesic distance + arc sampling for the ruler (A5). |
| `shading.ts` | Hillshade + contour derivation from per-cell height (A4). |
| `labels.ts` | The map-label model (`MapLabel`, styling, declutter) for 2D/3D labels (A1). |
| `geo.ts` | Shared geometry helpers (`insideTri`, `barycentric`, projection) — single source so export and Map2D project identically. |
| `dymaxion.ts` | Icosahedral geometry + Dymaxion net math + d3 polyhedral projection. |
| `colors.ts` | `getCellColor` — the single per-cell coloring function for every view mode + `buildFactionColorMap`. |
| `palette.ts` | Leaf palette constants; imports nothing but types, to keep `three` out of the worker bundle. |
| `vec.ts` | Minimal vector math; same worker-bundle-hygiene purpose as `palette.ts`. |
| `export.ts` | PNG raster export, JSON save/load + validation, `localStorage` saves. |
| `exportVector.ts` | SVG + GeoJSON vector export (E1/E2). |
| `exportGLB.ts` | GLB 3D export (mesh + rivers + city markers). |
| `worldGenClient.ts` | Worker client + `handleWorkerMessage` seam (testable). |
| `worldTransfer.ts` | Structure-of-Arrays (de)serialization across the worker boundary. |
| `paintUtils.ts` | Brush BFS + terrain/biome/political stroke functions + undo snapshots. |

### UI (`components/`)
| Module | Purpose |
|--------|---------|
| `shell/ShellApp.tsx` | The live app root: consumes `useWorldEngine`, wires everything. |
| `shell/{WideShell,NarrowShell}.tsx` | Responsive layout frames (desktop / mobile). |
| `shell/{DesignShell,shellKit}.tsx` | F1 layout prototype + shared shell primitives. |
| `Controls.tsx` | The parameter sidebar (System / Geo / Climate / Civ / Export tabs). |
| `WorldViewer.tsx` | The R3F 3D globe scene. See [rendering.md](rendering.md). |
| `Map2D.tsx` | Canvas2D Mercator + Dymaxion renderer with a color-ID pick buffer. |
| `Inspector.tsx` | Clicked-cell HUD + world-edit inputs. |
| `EditToolbar.tsx` | Floating paint/edit mode HUD. |
| `{MiniMap,Legend,ViewControls,DymaxionPreview2D,Select,ConfirmDialog}.tsx` | Overlays and shared UI. |

### Services / state
| Module | Purpose |
|--------|---------|
| `hooks/useWorldEngine.ts` | The single state owner (above). |
| `services/gemini.ts` | Google Gemini lore wrapper; mutates `world.civData` in place. |

## Data flow

1. User moves a control → `setParams(prev => ({ ...prev, … }))` in the hook.
2. `handleGenerate` → worker → `generateWorld` → serialized world → rehydrate →
   `setWorld`.
3. `world` prop-drills to `WorldViewer` (globe) or `Map2D` (2D).
4. Each cell is colored by `getCellColor(cell, viewMode, seaLevel, factionColors?)`.
5. Click a cell → `handleInspect` → `Inspector` reads `world.cells[id]`.

Civ/lore/edit paths mutate `world` in place and `setWorld({ ...world })` to
re-render (see [invariants.md](invariants.md) §7–§9).
