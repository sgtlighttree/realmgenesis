# RealmGenesis 3D -- Architecture Documentation

## Table of Contents

- [Overview](#overview)
- [LLM Quick-Navigation Guide](#llm-quick-navigation-guide)
- [System Architecture](#system-architecture)
- [Technology Stack](#technology-stack)
- [Project Structure](#project-structure)
- [Module Dependency Graph](#module-dependency-graph)
- [Module API Reference](#module-api-reference)
- [Data Model](#data-model)
  - [Core Types](#core-types-typests)
  - [WorldParams Reference](#worldparams-reference)
  - [BiomeType Classification Rules](#biometype-classification-rules)
  - [Cell Lifecycle](#cell-lifecycle)
- [Core Systems](#core-systems)
  - [World Generation Pipeline](#world-generation-pipeline)
  - [Climate Simulation](#climate-simulation)
  - [Civilization Engine](#civilization-engine)
  - [AI Lore Service](#ai-lore-service)
- [Rendering Architecture](#rendering-architecture)
  - [3D Globe Viewer](#3d-globe-viewer)
  - [2D Map Viewer](#2d-map-viewer)
  - [Dymaxion Projection](#dymaxion-projection)
- [State Management](#state-management)
- [Component Hierarchy](#component-hierarchy)
- [Data Flow](#data-flow)
- [Export & Persistence](#export--persistence)
- [Build & Deployment](#build--deployment)
- [Key Invariants & Gotchas](#key-invariants--gotchas)

---

## Overview

RealmGenesis 3D is a browser-based procedural fantasy world engine that simulates planetary geography, climate, biomes, and political systems on a sphere. It runs entirely client-side with no backend, using seeded random number generation for reproducibility.

The application generates worlds through a multi-stage pipeline: tectonic plate simulation → height map generation → hydraulic/thermal erosion → climate modeling → biome classification → river formation → civilization expansion. Results are visualized as an interactive 3D globe (Three.js), a 2D Mercator map, or an experimental Dymaxion (icosahedral) projection.

---

## LLM Quick-Navigation Guide

This section maps common tasks to the exact files and functions you should read first. The codebase is a client-side SPA with no backend — all logic is in-browser.

| Task | File | Key Symbol |
|------|------|-----------|
| Understand world generation | `utils/worldGen.ts` | `generateWorld()` (line 491) |
| Understand all data types | `types.ts` | `Cell`, `WorldData`, `WorldParams`, `BiomeType` |
| Understand all app state | `App.tsx` | `DEFAULT_PARAMS` + 15 `useState` calls (lines 13–65) |
| Understand cell color logic | `utils/colors.ts` | `getCellColor(cell, mode, seaLevel)` |
| Understand 3D rendering | `components/WorldViewer.tsx` | Full file — React Three Fiber scene |
| Understand 2D rendering | `components/Map2D.tsx` | Full file — Canvas2D with d3 projections |
| Understand Dymaxion math | `utils/dymaxion.ts` | `buildDymaxionNet()`, `createDymaxionProjection()` |
| Understand AI lore | `services/gemini.ts` | `generateWorldLore()` |
| Understand save/export | `utils/export.ts` | `exportMap()`, `saveMapToBrowser()`, `saveMapConfig()` |
| Understand 3D GLB export | `utils/exportGLB.ts` | `exportGLB(world, viewMode)` |
| Understand civilization generation | `utils/worldGen.ts` | `recalculateCivs()` (line 805), `recalculateProvinces()` (line 922) |
| Understand RNG | `utils/rng.ts` | `RNG` class (Mulberry32), `SimplexNoise` class |
| Understand biome rules | `utils/worldGen.ts` | `determineBiome()` (line 373) |

**Mental model of the data flow:**
1. User sets `WorldParams` in `Controls.tsx` → passed up to `App.tsx`
2. `App.tsx` calls `generateWorld(params)` → returns `WorldData` (cells + rivers + civData)
3. `WorldData` flows down via props to `WorldViewer` or `Map2D`
4. Each cell is colored by `getCellColor(cell, viewMode, seaLevel)` in `colors.ts`
5. User can click a cell → `Inspector` reads from `world.cells[cellId]`

**Key architectural constraint:** All state lives in `App.tsx` and is prop-drilled. There is no Context, Redux, or Zustand. If you need to trace a value, follow it up to `App.tsx`.

---

## System Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│                        Browser (Client)                         │
│                                                                 │
│  ┌──────────────┐  ┌──────────────────┐  ┌───────────────────┐  │
│  │   React UI   │  │  World Engine    │  │  Rendering Engine │  │
│  │  (Controls,  │  │  (worldGen.ts,   │  │  (Three.js +      │  │
│  │  Inspector,  │◄─┤  rng.ts, colors, │◄─┤   Canvas 2D)      │  │
│  │  Legend)     │  │  dymaxion.ts)    │  │                   │  │
│  └──────┬───────┘  └────────┬─────────┘  └─────────┬─────────┘  │
│         │                   │                       │            │
│  ┌──────▼───────────────────▼───────────────────────▼─────────┐  │
│  │                    App.tsx (Orchestrator)                   │  │
│  │  - State management (useState)                              │  │
│  │  - Generation lifecycle (AbortController)                   │  │
│  │  - View mode routing                                        │  │
│  │  - Lore API integration                                     │  │
│  └──────────────────────────┬─────────────────────────────────┘  │
│                             │                                    │
│  ┌──────────────────────────▼─────────────────────────────────┐  │
│  │              External Services                              │  │
│  │  - Google Gemini AI (lore generation, BYOK)                 │  │
│  │  - localStorage (map persistence)                           │  │
│  │  - File system (import/export JSON, image download)         │  │
│  └─────────────────────────────────────────────────────────────┘  │
└─────────────────────────────────────────────────────────────────┘
```

---

## Technology Stack

| Layer | Technology | Purpose |
|-------|-----------|---------|
| **Framework** | React 19 | UI component library |
| **Language** | TypeScript 5.8 | Type-safe development |
| **Build** | Vite 6 | Dev server, bundler, HMR |
| **3D Rendering** | Three.js 0.182 | WebGL scene graph |
| **3D React Bridge** | @react-three/fiber 9.5 | React renderer for Three.js |
| **3D Helpers** | @react-three/drei 10.7 | OrbitControls, Stars, Text |
| **Geo Math** | d3 7.9, d3-geo-voronoi 2.1, d3-geo-projection 4.0 | Voronoi tessellation, map projections |
| **AI** | @google/genai 1.38 | Gemini API client |
| **Icons** | lucide-react 0.563 | SVG icon library |
| **Styling** | Tailwind CSS (CDN) | Utility-first CSS |

---

## Project Structure

```
realmgenesis/
├── index.html                  # HTML shell (Tailwind CDN, CSP, root mount)
├── index.tsx                   # React DOM entry point (StrictMode mount)
├── App.tsx                     # Root component: state, orchestration, layout
├── types.ts                    # All TypeScript interfaces and enums
│
├── components/
│   ├── Controls.tsx            # Sidebar panel: 5 tabs (Sys/Geo/Clim/Civ/Exp)
│   ├── WorldViewer.tsx         # 3D globe: Three.js scene with mesh, markers, rivers
│   ├── Map2D.tsx               # 2D canvas: Mercator + Dymaxion projections
│   ├── DymaxionPreview2D.tsx   # Small preview for Dymaxion orientation
│   ├── Inspector.tsx           # HUD overlay: clicked cell data display
│   ├── MiniMap.tsx             # Bottom-right equirectangular overview
│   └── Legend.tsx              # Bottom-left biome color legend
│
├── utils/
│   ├── worldGen.ts             # Core engine: 12-stage generation pipeline
│   ├── rng.ts                  # Mulberry32 PRNG + 3D Simplex noise
│   ├── colors.ts               # View-mode color mapping (10 modes)
│   ├── dymaxion.ts             # Icosahedral geometry + Dymaxion net math
│   ├── export.ts               # Image export (PNG/raster), save/load, localStorage
│   └── exportGLB.ts            # 3D export: world mesh + rivers + city markers → GLB
│
├── services/
│   └── gemini.ts               # Google Gemini AI wrapper (lore generation)
│
├── vite.config.ts              # Vite config (port 3000, env injection)
├── tsconfig.json               # TypeScript config (ES2022, path alias @/*)
├── package.json                # Dependencies and scripts
├── public/
│   └── _redirects              # Netlify SPA fallback
└── .codacy/                    # Codacy static analysis configuration
```

---

## Module Dependency Graph

```
index.tsx
  └── App.tsx
        ├── types.ts
        ├── utils/worldGen.ts
        │     ├── utils/rng.ts
        │     ├── utils/colors.ts   (BIOME_COLORS only)
        │     └── types.ts
        ├── utils/export.ts
        │     ├── utils/dymaxion.ts
        │     ├── utils/colors.ts
        │     └── types.ts
        ├── utils/exportGLB.ts
        │     ├── utils/colors.ts
        │     └── types.ts
        ├── services/gemini.ts
        │     └── types.ts
        ├── components/Controls.tsx
        │     ├── utils/export.ts
        │     ├── utils/exportGLB.ts
        │     ├── services/gemini.ts
        │     └── types.ts
        ├── components/WorldViewer.tsx
        │     ├── utils/colors.ts
        │     ├── utils/dymaxion.ts
        │     └── types.ts
        ├── components/Map2D.tsx
        │     ├── utils/colors.ts
        │     ├── utils/dymaxion.ts
        │     └── types.ts
        ├── components/Inspector.tsx
        │     └── types.ts
        ├── components/Legend.tsx
        │     └── utils/colors.ts
        ├── components/MiniMap.tsx
        │     ├── utils/colors.ts
        │     └── types.ts
        └── components/DymaxionPreview2D.tsx
              ├── utils/dymaxion.ts
              └── types.ts
```

`utils/colors.ts` and `types.ts` are the most widely imported modules — nearly every component depends on them. `utils/worldGen.ts` is the only module that imports `utils/rng.ts`.

---

## Module API Reference

This section documents the public API surface of each module. Internal helpers (not exported) are noted separately where they are architecturally significant.

### `utils/worldGen.ts`

**Exported functions:**

| Function | Signature | Description |
|----------|-----------|-------------|
| `generateWorld` | `(params: WorldParams, onLog?: (msg: string) => void, signal?: AbortSignal) => Promise<WorldData>` | Runs the full 12-stage async pipeline; cancellable via AbortSignal |
| `recalculateCivs` | `(world: WorldData, params: WorldParams, onLog?: (msg: string) => void) => WorldData` | Replaces faction/territory data on an existing world; called independently without re-generating terrain |
| `recalculateProvinces` | `(world: WorldData, params: WorldParams) => WorldData` | Subdivides existing factions into provinces and places towns |

**Key internal functions (not exported, but architecturally significant):**

| Function | Description |
|----------|-------------|
| `generateFibonacciSphere(samples, rng, jitter)` | Distributes N points evenly on a unit sphere using the golden angle |
| `fbm(simplex, x, y, z, octaves, persistence, lacunarity)` | Fractal Brownian Motion — layered simplex noise for terrain height |
| `applyHydraulicErosion(cells, iterations, seaLevel, signal?)` | Simulates water flow, sediment transport, and deposition |
| `applyThermalErosion(cells, iterations, signal?)` | Smooths steep slopes by talus angle |
| `determineBiome(height, temp, moisture, seaLevel)` | Maps physical params to a `BiomeType` value (see classification table below) |

---

### `utils/rng.ts`

**Exported classes:**

| Class | Constructor | Key Methods | Description |
|-------|-------------|-------------|-------------|
| `RNG` | `new RNG(seedStr: string)` | `next(): number`, `range(min, max): number` | Mulberry32 seeded PRNG; FNV-1a hash converts the string seed to a 32-bit integer |
| `SimplexNoise` | `new SimplexNoise(rng: RNG)` | `noise3D(x, y, z): number` | 3D Simplex noise seeded from an RNG instance; returns values in approximately [-1, 1] |

Multiple independent `RNG` instances are created inside `generateWorld` with different seed suffixes (e.g., `seed + '_macro'`, `seed + '_plates_h'`, `seed + '_civs'`) to keep subsystems independent.

---

### `utils/colors.ts`

**Exported symbols:**

| Symbol | Type | Description |
|--------|------|-------------|
| `BIOME_COLORS` | `Record<BiomeType, string>` | Hex color string for each of the 15 biome types; used by `Legend` and as a fallback in export |
| `getCellColor` | `(cell: Cell, mode: ViewMode, seaLevel: number) => THREE.Color` | Primary color-mapping function; switches on `ViewMode` to produce a `THREE.Color` for a cell |

**ViewMode color logic summary:**

| ViewMode | Logic |
|----------|-------|
| `biome` | Direct lookup in `BIOME_COLORS` |
| `satellite` | Biome-aware realistic coloring with rock/snow blending at high elevation |
| `height` | HSL gradient: deep blue (ocean) → green (low land) → grey/white (peaks) |
| `height_bw` | Greyscale by normalized height |
| `temperature` | HSL hue sweep: blue (cold) → red (hot), range −30°C to 50°C |
| `moisture` | Saturation-based: dark blue (ocean) → light blue/white (dry land) |
| `plates` | 18-color palette indexed by `cell.plateId` |
| `political` | 18-color palette by `cell.regionId`; uses `getProvinceVariant()` if `provinceId` is set; territorial waters blend faction color toward deep blue |
| `population` | Falls through to `biome` default (not separately implemented) |
| `province` | Falls through to `biome` default (province coloring is handled within `political` mode) |

---

### `utils/dymaxion.ts`

**Exported symbols:**

| Symbol | Signature | Description |
|--------|-----------|-------------|
| `DymaxionNetFace` | `type` | Describes a single face of the unfolded icosahedron net (vertices in 2D, face index, rotation) |
| `buildDymaxionNet` | `(layout: DymaxionLayout) => DymaxionNetFace[]` | Constructs the 2D net for the given layout. `'classic'` builds the Fuller/Dymaxion diagonal net via a spanning-tree unfold; `'blender'` returns a hardcoded net whose UV coordinates match Blender's default icosphere UV unwrap exactly |
| `createDymaxionProjection` | `(layout: DymaxionLayout) => GeoProjection` | Returns a d3-compatible `geoPolyhedral` projection using `geoGnomonic` per-face; used by `WorldViewer` and `DymaxionPreview2D` |

---

### `utils/exportGLB.ts`

**Exported functions:**

| Function | Signature | Description |
|----------|-----------|-------------|
| `exportGLB` | `(world: WorldData, viewMode: ViewMode) => void` | Builds a Three.js scene from `WorldData` and triggers a GLB binary download. Produces three named objects: `World` (mesh with per-vertex colors), `Rivers` (line geometry, present if `world.rivers` is non-empty), and `Capitals`/`Towns` (merged low-poly cylinder geometry, present if `world.civData` exists). Uses `GLTFExporter` from `three/examples/jsm`. |

**GLB scene structure:**
- `World` — `MeshStandardMaterial` with `vertexColors: true`; per-vertex colors match `getCellColor(cell, viewMode, seaLevel)`; same `hMult = 1 + height * 0.05` elevation relief as the 3D globe viewer
- `Rivers` — `LineSegments` with `LineBasicMaterial` (colour `#38bdf8`); each river path stored as consecutive vertex pairs
- `Capitals` — merged `MeshBasicMaterial` (colour `#ef4444`); 6-sided cylinders, radius 0.008, height 0.08
- `Towns` — merged `MeshBasicMaterial` (colour `#ffffff`); 5-sided cylinders, radius 0.005, height 0.04

Vertex colors are exported as the `COLOR_0` attribute in GLTF. In Blender, set the material's Base Color to **Vertex Color** (or use an Attribute node with name `COLOR_0`) to display them.

---

### `utils/export.ts`

**Exported symbols:**

| Symbol | Signature | Description |
|--------|-----------|-------------|
| `ExportResolution` | `type = 4096 \| 8192 \| 16384 \| 32768` | Valid pixel widths for image export |
| `ProjectionType` | `type` | `'equirectangular' \| 'mercator' \| 'winkeltripel' \| 'orthographic' \| 'robinson' \| 'mollweide' \| 'dymaxion'` |
| `exportMap` | `async (world, viewMode, resolution, projection, dymaxionSettings?) => void` | Renders world to canvas at target resolution and triggers a file download |
| `saveMapConfig` | `(params: WorldParams, world?: WorldData) => void` | Downloads a JSON config file (params + optional civData) |
| `loadMapConfig` | `async (file: File) => Promise<LoadedMap \| null>` | Parses and validates a JSON config file; returns `null` on error |
| `getSavedMaps` | `() => SavedMapEntry[]` | Returns all maps persisted in `localStorage` |
| `saveMapToBrowser` | `(name: string, params: WorldParams, civData?: CivData) => void` | Serializes params + civData to `localStorage` |
| `deleteSavedMap` | `(name: string) => void` | Removes a saved map entry from `localStorage` |
| `SavedMapEntry` | `interface` | `{ name, timestamp, params, civData? }` — the localStorage record structure |
| `LoadedMap` | `interface` | `{ params, civData? }` — result of a successful `loadMapConfig` call |

---

### `services/gemini.ts`

**Exported functions:**

| Function | Signature | Description |
|----------|-----------|-------------|
| `setRuntimeApiKey` | `(key: string) => void` | Sets the Gemini API key for the current session; resets the `GoogleGenAI` client instance so the new key is picked up; key is module-level, not persisted |
| `generateWorldLore` | `(world: WorldData) => Promise<LoreData>` | Calls Gemini with a structured prompt based on `world.params.loreLevel`; **mutates `world.civData` in-place** with generated names; returns `{ name, description }` |

The model used is `gemini-3-flash-preview` with `responseMimeType: "application/json"`. Key is sourced from `runtimeKey` (set via UI) or `process.env.API_KEY` (build-time env var).

---

## Data Model

### Core Types (`types.ts`)

#### Point
```typescript
{ x: number; y: number; z: number }
```
A 3D Cartesian coordinate on the unit sphere. All cell centers and vertices are Points.

#### Cell
```typescript
{
  id: number,           // Index into WorldData.cells array (stable)
  center: Point,        // 3D position of cell centroid on unit sphere
  vertices: Point[],    // Voronoi polygon vertices (variable count)
  neighbors: number[],  // IDs of adjacent cells (graph edges)

  // Physical (populated by stages 3–11)
  height: number,       // 0–1 normalized; seaLevel ~0.55 by default
  plateId: number,      // Tectonic plate index
  temperature: number,  // Celsius; −50 to 50 typical range
  moisture: number,     // 0–1; 0 = arid, 1 = saturated
  biome: BiomeType,     // Determined from height + temp + moisture
  flux?: number,        // Water flux accumulation (erosion stage)

  // Political (populated after recalculateCivs / recalculateProvinces)
  regionId?: number,    // Faction ID; undefined for unclaimed cells
  provinceId?: number,  // Province ID local to faction
  isCapital?: boolean,
  isTown?: boolean,
  population?: number,
}
```
The world is composed of N cells (default 5,000). Cell ID equals its index in `WorldData.cells`.

#### BiomeType (enum with 15 values)
```typescript
enum BiomeType {
  OCEAN, DEEP_OCEAN,                              // water
  ICE_CAP, TUNDRA,                               // polar (E)
  HOT_DESERT, COLD_DESERT, STEPPE,               // dry (B)
  TROPICAL_RAINFOREST, TROPICAL_SAVANNA,         // tropical (A)
  MEDITERRANEAN, TEMPERATE_FOREST,               // temperate (C)
  TEMPERATE_RAINFOREST, BOREAL_FOREST,           // continental (D)
  BEACH, VOLCANIC,                               // special
}
```

#### WorldData
```typescript
{
  cells: Cell[],
  params: WorldParams,           // The params used to generate this world
  geoJson: Record<string, unknown>, // d3 GeoJSON FeatureCollection, cached for export
  civData?: CivData,
  rivers?: Point[][],            // Array of smoothed river paths (Point arrays)
}
```

#### Political Hierarchy
```
WorldData.civData
  └── factions: FactionData[]
        ├── id, name, color, capitalId, totalPopulation, description?
        └── provinces: ProvinceData[]
              ├── id, name, totalPopulation, color?
              └── towns: TownData[]
                    └── name, cellId, population, isCapital
```

#### View & Display Modes
```typescript
type DisplayMode    = 'globe' | 'mercator' | 'dymaxion'
type ViewMode       = 'biome' | 'height' | 'height_bw' | 'temperature' |
                      'moisture' | 'plates' | 'political' | 'population' |
                      'province' | 'satellite'
type InspectMode    = 'click' | 'hover' | 'off'
type DymaxionLayout = 'classic' | 'blender'
// 'classic' → Fuller diagonal net (default live view + raster export)
// 'blender' → Blender icosphere UV net (export-only; square image, net in lower ~47%)
```

---

### WorldParams Reference

All parameters live in `types.ts`. Defaults are set in `App.tsx` (`DEFAULT_PARAMS`).

#### System
| Parameter | Type | Default | Range/Options | Controls |
|-----------|------|---------|---------------|---------|
| `mapName` | `string` | `'map'` | Any string | Display name and export filename |
| `points` | `number` | `5000` | 500–20,000+ | Number of Voronoi cells; higher = more detail, slower |
| `seed` | `string` | `'realmgenesis'` | Any string | Terrain RNG seed (hashed to uint32) |
| `planetRadius` | `number` | `6371` | km | Display only; affects no simulation logic |
| `axialTilt` | `number` | `23.5` | 0–90° | Modulates temperature latitudinal gradient |

#### Geography
| Parameter | Type | Default | Range/Options | Controls |
|-----------|------|---------|---------------|---------|
| `landStyle` | `LandStyle` | `'Continents'` | `'Continents' \| 'Archipelago' \| 'Islands' \| 'Pangea' \| 'Custom'` | Sets terrain preset (adjusts noise/mask params) |
| `cellJitter` | `number` | `0.5` | 0–1 | Randomizes Fibonacci sphere points; 0 = regular grid |

#### Advanced Terrain
| Parameter | Type | Default | Range/Options | Controls |
|-----------|------|---------|---------------|---------|
| `noiseScale` | `number` | `0.4` | 0.1–2.0 | FBM feature frequency; lower = broader continents |
| `ridgeBlend` | `number` | `0.1` | 0–1 | 0 = smooth FBM, 1 = sharp ridged noise mountains |
| `maskType` | `MaskType` | `'None'` | `'None' \| 'Pangea'` | Optional supercontinent height mask |
| `warpStrength` | `number` | `0.5` | 0–2 | Domain warp intensity for organic shapes |
| `plateInfluence` | `number` | `0.5` | 0–2 | Weight of tectonic stress on height |
| `erosionIterations` | `number` | `2` | 0–50 | Hydraulic + thermal erosion pass count |
| `plates` | `number` | `12` | 2–30 | Number of tectonic plates |
| `seaLevel` | `number` | `0.55` | 0–1 | Height threshold separating ocean from land |
| `roughness` | `number` | `0.5` | 0–1 | FBM persistence (controls terrain roughness) |
| `detailLevel` | `number` | `2` | 1–8 | FBM octave count |

#### Climate
| Parameter | Type | Default | Range/Options | Controls |
|-----------|------|---------|---------------|---------|
| `baseTemperature` | `number` | `30` | −30 to 60°C | Equatorial temperature before elevation adjustment |
| `poleTemperature` | `number` | `−30` | −60 to 10°C | Polar temperature |
| `rainfallMultiplier` | `number` | `1.0` | 0.1–3.0 | Scales moisture values globally |
| `moistureTransport` | `number` | `0.5` | 0–1 | How far wind carries moisture inland |
| `temperatureVariance` | `number` | `5` | 0–20 | Noise added to temperature for local variation |

#### Political
| Parameter | Type | Default | Range/Options | Controls |
|-----------|------|---------|---------------|---------|
| `numFactions` | `number` | `6` | 1–20 | Number of political factions |
| `civSeed` | `string` | `'realmgenesis_civs'` | Any string | Separate RNG seed for faction placement |
| `borderRoughness` | `number` | `0.2` | 0–1 | Noise on Dijkstra costs for irregular borders |
| `civSizeVariance` | `number` | `0.5` | 0–1 | How unequal faction sizes can be |
| `waterCrossingCost` | `number` | `0.8` | 0–1 | Dijkstra cost multiplier for crossing water |
| `territorialWaters` | `number` | `0.15` | 0–1 | Max graph distance from land to claim water cells |
| `capitalSpacing` | `number` | `0.5` | 0–1 | Minimum angular distance between faction capitals |
| `provinceSize` | `number` | `0.5` | 0.1–1.0 | Province target size (0.1 = small/many, 1.0 = large/few) |

#### Meta
| Parameter | Type | Default | Options | Controls |
|-----------|------|---------|---------|---------|
| `loreLevel` | `1 \| 2 \| 3` | `1` | 1/2/3 | Gemini prompt depth: 1 = names only, 2 = + provinces, 3 = + backstories |

---

### BiomeType Classification Rules

`determineBiome(height, temp, moisture, seaLevel)` in `utils/worldGen.ts` (line 373):

```
height < seaLevel:
  height < seaLevel × 0.6  →  DEEP_OCEAN
  otherwise                →  OCEAN

height ≥ seaLevel (land):
  landH = (height − seaLevel) / (1 − seaLevel)   [0–1 normalized land elevation]

  landH < 0.02 AND temp > 15°C   →  BEACH     (coastal fringe)
  landH > 0.85 AND temp > −5°C   →  VOLCANIC  (extreme high elevation)

  temp < −10°C                   →  ICE_CAP
  temp < 0°C                     →  TUNDRA

  moisture < 0.15:
    temp > 25°C   →  HOT_DESERT
    temp > 10°C   →  STEPPE
    otherwise     →  COLD_DESERT

  moisture < 0.40:
    temp > 25°C   →  TROPICAL_SAVANNA
    temp > 10°C   →  MEDITERRANEAN
    otherwise     →  STEPPE

  moisture ≥ 0.40:
    temp > 25°C   →  TROPICAL_RAINFOREST
    temp > 15°C   →  TEMPERATE_RAINFOREST
    temp > 5°C    →  TEMPERATE_FOREST
    otherwise     →  BOREAL_FOREST
```

Note: BEACH and VOLCANIC are checked before temperature/moisture rules, acting as overrides for extreme elevation bands.

---

### Cell Lifecycle

Each field on `Cell` is populated progressively through the generation pipeline:

| Stage | Pipeline Step | Fields Set |
|-------|--------------|-----------|
| **Stage 1** | Point distribution | (raw Point array, not yet Cell objects) |
| **Stage 2** | Voronoi tessellation | `id`, `center`, `vertices`, `neighbors` |
| **Stage 3** | Tectonic plate assignment | `plateId` |
| **Stage 4** | Connectivity enforcement | (re-assigns orphaned cells' `plateId`) |
| **Stage 5** | Stress calculation | (internal stress array, not stored on Cell) |
| **Stage 6** | Height generation | `height` (raw, unnormalized) |
| **Stage 7** | Normalization | `height` (scaled to 0–1) |
| **Stage 8** | Hydraulic erosion | `height` (modified in-place) |
| **Stage 9** | Thermal erosion | `height` (modified in-place) |
| **Stage 10** | Climate simulation | `temperature`, `moisture` |
| **Stage 11** | Biome assignment | `biome` |
| **Stage 12** | River generation | `flux` (water flux accumulation) |
| **`recalculateCivs`** | Faction expansion | `regionId` |
| **`recalculateProvinces`** | Province subdivision | `provinceId`, `isCapital`, `isTown`, `population` |

After Stage 12, `WorldData.rivers` is populated (array of Point-path arrays for smooth river rendering). `WorldData.geoJson` is built during Stage 2 and cached for all export operations.

---

## Core Systems

### World Generation Pipeline

The `generateWorld()` function in `utils/worldGen.ts` (line 491) executes a 12-stage async pipeline:

```
Stage 1:  Point Distribution
          └── generateFibonacciSphere(points, rng, cellJitter)
              Uses golden angle (φ = π(3−√5)) with optional jitter

Stage 2:  Voronoi Tessellation
          └── d3-geo-voronoi spherical Voronoi → Cell graph with neighbors
              Also builds WorldData.geoJson (FeatureCollection for export)

Stage 3:  Tectonic Plate Assignment
          └── K-means-like clustering with warp noise for organic boundaries
              Each cell assigned to nearest of N plate centers

Stage 4:  Connectivity Enforcement
          └── Flood-fill BFS per plate; orphaned cells reassigned to
              majority-neighbor plate (enforces single connected component)

Stage 5:  Stress Calculation
          └── Computes convergence/divergence at plate boundaries
              (dot product of plate motion vectors); used in Stage 6

Stage 6:  Height Generation
          └── fbm(simplex, x, y, z, octaves, persistence, lacunarity)
              + ridged noise (blended via ridgeBlend)
              + plate influence (plateInfluence weight on stress map)
              + tectonic stress (convergent = mountain, divergent = rift)
              + detail noise (high-frequency variation)
              + optional Pangea mask (maskType = 'Pangea')

Stage 7:  Normalization
          └── Scales heights to 0–1 range

Stage 8:  Hydraulic Erosion (async, cancellable)
          └── applyHydraulicErosion(cells, erosionIterations, seaLevel, signal)
              Water flow simulation: rainfall → flow → sediment transport → deposition

Stage 9:  Thermal Erosion (async, cancellable)
          └── applyThermalErosion(cells, erosionIterations, signal)
              Talus slope smoothing: steep slopes shed material to neighbors

Stage 10: Climate Simulation
          └── Wind vectors (latitude-based prevailing winds)
              → moisture transport from ocean (8 iterative passes)
              → temperature (latitude gradient + axial tilt + elevation lapse rate)

Stage 11: Biome Assignment
          └── determineBiome(height, temp, moisture, seaLevel) per cell
              See classification table above

Stage 12: River Generation
          └── Priority-Flood depression filling (MinHeap) → ensure drainage
              → flux accumulation from precipitation + upstream flow
              → path tracing from high-flux sources to ocean sinks
              → CatmullRomCurve3 smoothing for render-ready Point paths
```

Between stages, the pipeline calls `await new Promise(r => setTimeout(r, 0))` to yield to the browser event loop, and checks `signal.aborted` to support cancellation.

#### Key Algorithms

- **Fibonacci Sphere**: Distributes points evenly on a sphere using the golden angle, with configurable jitter for natural variation.
- **Spherical Voronoi**: Uses `d3-geo-voronoi` to compute Voronoi cells on a sphere, producing a dual graph with neighbor relationships.
- **Fractal Brownian Motion (fbm)**: Layered simplex noise at multiple octaves for terrain height generation.
- **Ridged Noise**: Absolute-value noise for sharp mountain ridges, blended with FBM via `ridgeBlend` parameter.
- **Priority-Flood**: Depression-filling algorithm using a MinHeap priority queue for realistic river drainage.
- **Hydraulic Erosion**: Simulates rainfall, water flow, sediment transport, and deposition over iterations.

### Climate Simulation

The climate system models:

1. **Wind Vectors**: Prevailing wind patterns based on latitude (trade winds, westerlies, polar easterlies).
2. **Moisture Transport**: Iterative moisture propagation from ocean cells inland over 8 passes, with orographic rain shadow effects.
3. **Temperature**: Latitude-based gradient from equator to poles, modulated by:
   - Axial tilt (seasonal variation)
   - Elevation (lapse rate: approximately −6°C per normalized unit of height above sea level, scaled to produce ~−60°C at max elevation)
   - `temperatureVariance` parameter (simplex noise offset)

### Civilization Engine

Two-phase political simulation, each independently callable:

#### Phase 1: Faction Expansion (`recalculateCivs`, line 805)
1. Places faction capitals using `capitalSpacing` as minimum angular separation
2. Expands territories outward using Dijkstra's algorithm
3. Terrain-dependent costs: ocean cells cost `waterCrossingCost × base`, mountains/deserts add penalties, `borderRoughness` injects random noise
4. Water cells within `territorialWaters` graph-distance of a land cell are claimed as territorial waters
5. `civSizeVariance` modulates how different faction sizes can be by adjusting initial cell budgets

#### Phase 2: Province Subdivision (`recalculateProvinces`, line 922)
1. Subdivides each faction into provinces based on `provinceSize` parameter
2. Places towns within each province
3. Assigns population based on biome suitability (fertile biomes = higher, deserts/tundra = lower)
4. Province 0 of each faction contains the capital town

### AI Lore Service

`services/gemini.ts` integrates Google Gemini for procedural world lore:

- **Model**: `gemini-3-flash-preview` with JSON response mode (`responseMimeType: "application/json"`)
- **API Key**: Ephemeral; set at runtime via `setRuntimeApiKey()` or baked into build as `process.env.API_KEY`; never persisted
- **Lore Levels**:
  - Level 1: World name, description, faction names, capital names
  - Level 2: + Province and town names
  - Level 3: + Faction backstories (~50 words each)
- `generateWorldLore()` **mutates `world.civData` in-place** — names are applied directly to `FactionData` and `ProvinceData` objects; the caller must `setWorld({ ...world })` to trigger a re-render

---

## Rendering Architecture

### 3D Globe Viewer

`components/WorldViewer.tsx` renders an interactive Three.js scene via `@react-three/fiber`:

| Element | Implementation |
|---------|---------------|
| **World Mesh** | Triangle-based geometry with vertex colors from `getCellColor(cell, viewMode, seaLevel)`. Each Voronoi cell is triangulated from its center to vertices. |
| **City Markers** | `InstancedMesh` cylinders: red for capitals, white for towns |
| **River Lines** | `LineSegments` with `CatmullRomCurve3` smoothing |
| **Faction Borders** | Line segments between adjacent cells of different regions |
| **Country Labels** | 3D `<Text>` components (drei) for faction names in political mode |
| **Lat/Long Grid** | 10-degree latitude/longitude grid lines |
| **Dymaxion Overlay** | Rotatable icosahedron wireframe |
| **Background** | `<Stars>` component (drei) |
| **Camera** | `OrbitControls` with auto-rotation (paused in overlay mode) |

R3F element names (e.g., `"bufferGeometry"`, `"lineSegments"`) are passed as strings to bypass TypeScript's JSX element type checking — this is intentional and documented in AGENTS.md.

Pointer interaction supports click and hover inspection, propagating cell IDs to the `Inspector` HUD.

### 2D Map Viewer

`components/Map2D.tsx` uses an offscreen Canvas2D for raster rendering:

- **Mercator Mode**: `d3.geoMercator` projection with GeoJSON polygon features
- **Dymaxion Mode**: Pixel-by-pixel reprojection from an equirectangular source through the icosahedral net (via `buildDymaxionNet`)
- **Adaptive DPR**: Reduces device pixel ratio to 1× during interaction for 60fps, sharpens to 2–3× when settled
- **Pan/Zoom**: Drag to pan, scroll wheel to zoom, throttled via `requestAnimationFrame`
- **Hit Detection**: Color-coded pick buffer maps screen pixels back to cell IDs
- **River Rendering**: Antimeridian crossing detection for correct line wrapping

### Dymaxion Projection

`utils/dymaxion.ts` implements the Buckminster Fuller Dymaxion map:

1. **Icosahedron Geometry**: 12 vertices, 20 triangular faces using the golden ratio (PHI = (1+√5)/2)
2. **Face Construction**: Oriented faces with correct winding order (cross-product normals)
3. **2D Net** (`buildDymaxionNet`): Unfolds the icosahedron into a flat layout with barycentric coordinate transforms for pixel-level reprojection
4. **D3 Integration** (`createDymaxionProjection`): Creates a `d3.geoPolyhedral` projection using `geoGnomonic` per-face, suitable for GeoJSON rendering
5. **Orientation**: Configurable `lon`/`lat`/`roll` for rotating the projection center (stored in `DymaxionSettings`)

The Dymaxion projection is available in both the 3D viewer (as an icosahedron wireframe overlay) and the 2D viewer (as a full raster map via `buildDymaxionNet`).

---

## State Management

The application uses React `useState` at the `App.tsx` level with prop drilling. No external state management library is used.

### App State
| State | Type | Purpose |
|-------|------|---------|
| `params` | `WorldParams` | Generation configuration |
| `world` | `WorldData \| null` | Generated world data |
| `viewMode` | `ViewMode` | Current visualization mode |
| `displayMode` | `DisplayMode` | Globe / Mercator / Dymaxion |
| `inspectMode` | `InspectMode` | Click / Hover / Off |
| `inspectedCellId` | `number \| null` | Selected cell for inspector |
| `isGenerating` | `boolean` | Loading state |
| `logs` | `string[]` | Generation console output |
| `lore` | `LoreData \| null` | AI-generated world lore |
| `isLoreLoading` | `boolean` | Lore API loading state |
| `showGrid` | `boolean` | Toggle lat/long grid |
| `showRivers` | `boolean` | Toggle river display |
| `sidebarOpen` | `boolean` | Mobile sidebar toggle |
| `dymaxionSettings` | `DymaxionSettings` | Dymaxion projection config |
| `apiKey` | `string` | Gemini API key (ephemeral) |

### Generation Lifecycle

World generation uses `AbortController` for cancellation:

```
handleGenerate()
  ├── Abort previous generation (if running)
  ├── Create new AbortController
  ├── Set isGenerating = true
  ├── await generateWorld(params, onLog, signal)
  │     ├── Async pipeline stages (1-12)
  │     ├── Check signal.aborted between stages
  │     └── Throw "Generation Cancelled" if aborted
  ├── Set world = newWorld
  └── Set isGenerating = false
```

The controller reference is stored in a `useRef` to persist across renders. The `onLog` callback calls `setLogs(prev => [...prev, msg])` which appends to the console panel in `Controls.tsx`.

---

## Component Hierarchy

```
<App>
  ├── <Controls>                          # Left sidebar (5 tabs)
  │     ├── System tab: seed, points, name, presets, console
  │     ├── Geography tab: plates, sea level, terrain params
  │     ├── Climate tab: temperature, rainfall, moisture
  │     ├── Civilizations tab: factions, provinces, AI lore
  │     └── Export tab: image export, save/load, import
  │
  ├── <WorldViewer> (displayMode = 'globe')
  │     ├── <Canvas>
  │     │     ├── <WorldMesh>             # Triangle mesh with vertex colors
  │     │     ├── <CityMarkers>           # InstancedMesh for capitals/towns
  │     │     ├── <RiverLines>            # Smoothed river paths
  │     │     ├── <FactionBorders>        # Border line segments
  │     │     ├── <CountryLabels>         # 3D text labels
  │     │     ├── <DymaxionOverlay>       # Icosahedron wireframe
  │     │     ├── <LatLongGrid>           # Grid lines
  │     │     └── <Stars>                 # Background stars
  │     └── <OrbitControls>               # Camera interaction
  │
  ├── <Map2D> (displayMode = 'mercator' | 'dymaxion')
  │     └── <canvas>                      # Offscreen 2D rendering
  │
  ├── <BiomeLegend>                       # Bottom-left overlay
  ├── <MiniMap>                           # Bottom-right overview (globe mode)
  └── <Inspector>                         # Floating HUD overlay
```

### Controls Component (~1,150 lines)

The largest UI component, organized into 5 tabs:

- **Sys**: Map name, resolution (points), seed (with lock/randomize), terrain presets (Continents, Pangea, Archipelago, Islands, Custom), auto-update toggle, generation console, generate/cancel buttons
- **Geo**: Plates, sea level, roughness, detail level, cell jitter, noise scale, ridge blend, mask type, warp strength, plate influence, erosion iterations
- **Clim**: Base temperature, pole temperature, rainfall multiplier, moisture transport, temperature variance, axial tilt
- **Civ**: Number of factions, civ seed, border roughness, size variance, water crossing cost, territorial waters, capital spacing, province size, AI lore generation, update civs/provinces buttons
- **Exp**: View mode selector, display mode selector, inspect mode, grid/rivers toggles, Dymaxion settings (lon/lat/roll), image export (resolution + projection), Dymaxion preview, browser storage manager, JSON import/export

---

## Data Flow

### World Generation
```
Controls (user input)
  └── setParams()
        └── handleGenerate()
              └── generateWorld(params, onLog, signal)
                    ├── Returns WorldData
                    └── Calls onLog() at each stage
                          └── addLog() → setLogs()
              └── setWorld(newWorld)
                    └── Re-renders WorldViewer / Map2D
```

### Cell Inspection
```
WorldViewer pointer event
  └── onInspect(cellId)
        └── setInspectedCellId(cellId)
              └── <Inspector> reads cell from world.cells[cellId]
                    └── Displays: biome, temperature, moisture,
                                  elevation, population, faction
```

### AI Lore Generation
```
Controls "Generate Lore" button
  └── handleGenerateLore()
        └── generateWorldLore(world)   ← mutates world.civData in-place
              └── Gemini API call (async)
                    └── Returns LoreData
                          └── setLore(newLore)
                          └── setWorld({ ...world })   ← shallow copy triggers re-render
```

### Map Save/Load
```
Save:
  Controls "Save" → saveMapToBrowser(name, params, civData?)
    └── Serialize params + civData → localStorage

Load:
  Controls "Load" → loadMapConfig(file)
    └── Parse JSON → validate params structure
    └── handleLoadWorld(params, savedCivData)
          └── Regenerate world (same seed → same terrain)
          └── Restore saved names/descriptions from civData
```

---

## Export & Persistence

### Image Export (`utils/export.ts`)

`exportMap()` renders the world to a canvas at configurable resolutions:
- **Resolutions**: 4K (4096px), 8K (8192px), 16K (16384px), 32K (32768px) width
- **Projections**: Equirectangular, Mercator, Winkel Tripel, Robinson, Mollweide, Orthographic, Dymaxion
- **Classic Dymaxion raster**: Pixel-by-pixel reprojection via `buildDymaxionNet('classic')`; auto-fit with padding; output is `width × round(width × 0.6)`
- **Blender Dymaxion raster**: Uses `buildDymaxionNet('blender')`; direct UV→pixel mapping (`px = u*W, py = (1-v)*H`); output is square (`width × width`). The net occupies the bottom ~47% of the image, matching Blender's UV space exactly for any icosphere subdivision level.

### 3D Export (`utils/exportGLB.ts`)

`exportGLB()` builds a Three.js scene from `WorldData` and exports it as a binary GLB file:
- **World mesh**: per-vertex colored `MeshStandardMaterial`; same fan-triangulation and `hMult` elevation as the live 3D globe
- **Rivers**: `LineSegments` geometry; GLTF LINES primitive; colour `#38bdf8`
- **City markers**: merged non-indexed cylinder geometry; capitals (6-sided, height 0.08) in red, towns (5-sided, height 0.04) in white; only included when `world.civData` is present
- Vertex colors are exported as GLTF attribute `COLOR_0`

### Browser Storage

Uses `localStorage` with a `SavedMapEntry[]` structure:
```typescript
{
  name: string,
  timestamp: number,
  params: WorldParams,
  civData?: CivData,
}
```

### JSON Config

`saveMapConfig` / `loadMapConfig` serialize/deserialize full world configuration as JSON files, including all parameters and civilization metadata. `loadMapConfig` validates the structure before use; malformed files return `null`.

---

## Build & Deployment

### Development
```bash
npm install        # Install dependencies
npm run dev        # Start Vite dev server on port 3000
```

### Production Build
```bash
npm run build      # Vite production build → dist/
npm run preview    # Preview production build locally
```

### Deployment Target: Netlify

- `public/_redirects` provides SPA fallback routing (`/* /index.html 200`)
- Environment variable `GEMINI_API_KEY` can be set at build time (baked in as `process.env.API_KEY`) or provided at runtime via the UI
- No server-side rendering; fully static SPA

### Build Configuration

- **Vite 6** with `@vitejs/plugin-react` for HMR
- **TypeScript**: ES2022 target, ESNext modules, `react-jsx` transform
- **Path alias**: `@/*` maps to project root (configured in tsconfig but intentionally unused — use relative imports)
- **CSP**: HTML meta tag allows self, Tailwind CDN, and Google Generative Language API

---

## Key Invariants & Gotchas

These are non-obvious facts that are critical for making correct changes:

1. **Single-threaded generation**: No Web Workers are used. The 12-stage pipeline runs on the main thread. UI responsiveness is maintained solely by `await new Promise(r => setTimeout(r, 0))` yields between stages. Long-running stage changes may introduce noticeable UI stalls.

2. **`@/` path alias is intentionally unused**: Although configured in `tsconfig.json` and `vite.config.ts`, all imports in the codebase use relative paths (`../types`, `./colors`). Never add `@/` imports — they work at runtime but contradict established convention.

3. **In-place civData mutation**: `generateWorldLore()` mutates `world.civData` objects directly (faction names, descriptions, province names). The caller must do `setWorld({ ...world })` (a shallow copy of `WorldData`) to trigger React re-renders. Do not replace `world.civData` with a new object — update the existing one in-place.

4. **R3F element names as strings**: In `WorldViewer.tsx`, Three.js elements are referenced as string literals (e.g., `<primitive object={...}>`, `"bufferGeometry"`) to bypass TypeScript's JSX type system. This is an established pattern; `@typescript-eslint/no-explicit-any` is set to `warn`, not `error`, for this reason.

5. **Ephemeral API key**: The Gemini API key stored in `apiKey` state and via `setRuntimeApiKey()` resets to empty on every page reload. It is never written to `localStorage`, cookies, or any persistent store. Do not add persistence for it.

6. **Cell ID stability**: A cell's `id` equals its index in `WorldData.cells`. This is stable within a single generated world (cells are never reordered), but is not stable across different generation runs — even with the same seed, layout changes could theoretically alter Voronoi topology.

7. **GeoJSON cache**: `WorldData.geoJson` is a `d3` GeoJSON `FeatureCollection` built during Stage 2 and cached for the lifetime of the world. Feature index `i` corresponds to `world.cells[i]`. Export functions rely on this cache — do not clear or replace it after generation.

8. **`seaLevel` is passed to `getCellColor`**: The function signature is `getCellColor(cell, viewMode, seaLevel)`. The third argument must be `world.params.seaLevel`, not a hardcoded value. Missing this causes incorrect ocean/land color boundaries in all rendering modes.

9. **Export resolutions are canvas-based**: Very large exports (16K, 32K) create large offscreen canvases. On low-memory devices or browsers with canvas size limits, these may fail silently or throw. 32K (32768px) exceeds most browser canvas limits and should be considered experimental.

10. **No test framework**: There are no automated tests. Quality gates are: `npm run build` succeeds, `npm run lint` has zero errors, TypeScript has zero type errors. All behavioral testing is manual via the browser.

11. **GLB vertex colors require a Blender material step**: `exportGLB` bakes cell colors into the GLTF `COLOR_0` vertex attribute. Blender does not display these automatically — the imported material must have its Base Color connected to an **Attribute** node (name: `COLOR_0`) or the viewport shading must be set to **Vertex Color**.

12. **Blender UV net is export-only**: `DymaxionLayout = 'blender'` produces a square image for use as a UV texture; it has no effect on the live 3D globe or 2D map views. The toggle in Dymaxion Controls only influences the raster export path in `exportMap()`.

13. **`exportGLB` builds geometry independently of the live scene**: The function creates fresh Three.js objects from `WorldData` — it does not capture the R3F canvas or scene ref. This means it works from any tab and any display mode, but always reflects the current `viewMode` color scheme, not any transient render state.
