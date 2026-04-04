# RealmGenesis 3D -- Architecture Documentation

## Table of Contents

- [Overview](#overview)
- [System Architecture](#system-architecture)
- [Technology Stack](#technology-stack)
- [Project Structure](#project-structure)
- [Data Model](#data-model)
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

---

## Overview

RealmGenesis 3D is a browser-based procedural fantasy world engine that simulates planetary geography, climate, biomes, and political systems on a sphere. It runs entirely client-side with no backend, using seeded random number generation for reproducibility.

The application generates worlds through a multi-stage pipeline: tectonic plate simulation → height map generation → hydraulic/thermal erosion → climate modeling → biome classification → river formation → civilization expansion. Results are visualized as an interactive 3D globe (Three.js), a 2D Mercator map, or an experimental Dymaxion (icosahedral) projection.

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
│   └── export.ts               # Image export, save/load, localStorage
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

## Data Model

### Core Types (`types.ts`)

#### Point
```
{x: number, y: number, z: number}
```
A 3D Cartesian coordinate on the unit sphere.

#### Cell
```
{
  id: number,
  center: Point,
  vertices: Point[],
  neighbors: number[],
  height: number,          // 0-1 normalized elevation
  plateId: number,
  temperature: number,     // Celsius
  moisture: number,        // 0-1
  biome: BiomeType,
  flux?: number,           // Water flux (erosion)
  regionId?: number,       // Faction ID
  provinceId?: number,     // Local province ID
  isCapital?: boolean,
  isTown?: boolean,
  population?: number,
}
```
A single Voronoi cell representing a discrete geographic unit. The world is composed of N cells (default 5000).

#### BiomeType (15 values)
Ocean, Deep Ocean, Ice Cap, Tundra, Hot Desert, Cold Desert, Steppe, Tropical Rainforest, Tropical Savanna, Mediterranean, Temperate Forest, Temperate Rainforest, Boreal Forest, Beach, Volcanic.

Classification follows a simplified Koppen system based on height, temperature, and moisture.

#### WorldParams (30+ configurable parameters)
Organized into categories:
- **System**: mapName, points, seed, planetRadius, axialTilt
- **Geography**: landStyle, cellJitter
- **Advanced Terrain**: noiseScale, ridgeBlend, maskType, warpStrength, plateInfluence, erosionIterations
- **Base**: plates, seaLevel, roughness, detailLevel
- **Climate**: baseTemperature, poleTemperature, rainfallMultiplier, moistureTransport, temperatureVariance
- **Political**: numFactions, civSeed, borderRoughness, civSizeVariance, waterCrossingCost, territorialWaters, capitalSpacing, provinceSize
- **Meta**: loreLevel

#### Political Hierarchy
```
WorldData
  └── civData
        └── factions[] (FactionData)
              └── provinces[] (ProvinceData)
                    └── towns[] (TownData)
```

#### WorldData
```
{
  cells: Cell[],
  params: WorldParams,
  geoJson: Record<string, unknown>,
  civData?: CivData,
  rivers?: Point[][],
}
```

#### View & Display Modes
- **DisplayMode**: `'globe' | 'mercator' | 'dymaxion'`
- **ViewMode**: `'biome' | 'height' | 'height_bw' | 'temperature' | 'moisture' | 'plates' | 'political' | 'population' | 'province' | 'satellite'`
- **InspectMode**: `'click' | 'hover' | 'off'`

---

## Core Systems

### World Generation Pipeline

The `generateWorld()` function in `utils/worldGen.ts` executes a 12-stage async pipeline:

```
Stage 1:  Point Distribution
          └── Fibonacci sphere with jitter (generateFibonacciSphere)

Stage 2:  Voronoi Tessellation
          └── d3-geo-voronoi spherical Voronoi → Cell graph with neighbors

Stage 3:  Tectonic Plate Assignment
          └── K-means-like clustering with warp noise for organic boundaries

Stage 4:  Connectivity Enforcement
          └── Ensures each plate is a single connected component (flood fill)

Stage 5:  Stress Calculation
          └── Computes tectonic stress at plate boundaries

Stage 6:  Height Generation
          └── FBM noise + ridged noise + plate influence + tectonic stress
              + detail noise + optional Pangea mask

Stage 7:  Normalization
          └── Scales heights to 0-1 range

Stage 8:  Hydraulic Erosion
          └── Water flow simulation with deposition (applyHydraulicErosion)

Stage 9:  Thermal Erosion
          └── Talus slope smoothing (applyThermalErosion)

Stage 10: Climate Simulation
          └── Wind vectors → moisture transport (8 passes) → temperature
              with axial tilt and latitude effects

Stage 11: Biome Assignment
          └── determineBiome() using height + temperature + moisture

Stage 12: River Generation
          └── Priority-Flood depression filling → flux accumulation
              → path tracing from sources to sinks
```

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
   - Elevation (lapse rate)
   - Temperature variance parameter

### Civilization Engine

Two-phase political simulation:

#### Phase 1: Faction Expansion (`recalculateCivs`)
1. Places faction capitals with spacing constraints
2. Expands territories using Dijkstra's algorithm with terrain-dependent costs
3. Water cells are claimed as territorial waters within a configurable distance
4. Border roughness adds natural-looking boundary irregularities

#### Phase 2: Province Subdivision (`recalculateProvinces`)
1. Subdivides each faction into provinces based on `provinceSize` parameter
2. Places towns within provinces
3. Calculates population based on biome suitability (higher for fertile biomes, lower for deserts/tundra)

### AI Lore Service

`services/gemini.ts` integrates Google Gemini for procedural world lore:

- **Model**: `gemini-3-flash-preview` with JSON response mode
- **API Key**: Ephemeral, set at runtime via UI input or `.env.local` (never persisted in app state)
- **Lore Levels**:
  - Level 1: World name, description, faction names, capital names
  - Level 2: + Province and town names
  - Level 3: + Faction backstories (~50 words each)
- Generated names are applied in-place to `WorldData.civData`

---

## Rendering Architecture

### 3D Globe Viewer

`components/WorldViewer.tsx` renders an interactive Three.js scene via `@react-three/fiber`:

| Element | Implementation |
|---------|---------------|
| **World Mesh** | Triangle-based geometry with vertex colors from `getCellColor()`. Each Voronoi cell is triangulated from its center to vertices. |
| **City Markers** | `InstancedMesh` cylinders: red for capitals, white for towns |
| **River Lines** | `LineSegments` with `CatmullRomCurve3` smoothing |
| **Faction Borders** | Line segments between adjacent cells of different regions |
| **Country Labels** | 3D `<Text>` components (drei) for faction names in political mode |
| **Lat/Long Grid** | 10-degree latitude/longitude grid lines |
| **Dymaxion Overlay** | Rotatable icosahedron wireframe |
| **Background** | `<Stars>` component (drei) |
| **Camera** | `OrbitControls` with auto-rotation (paused in overlay mode) |

Pointer interaction supports click and hover inspection, propagating cell IDs to the `Inspector` HUD.

### 2D Map Viewer

`components/Map2D.tsx` uses an offscreen Canvas2D for raster rendering:

- **Mercator Mode**: `d3.geoMercator` projection with GeoJSON polygon features
- **Dymaxion Mode**: Pixel-by-pixel reprojection from an equirectangular source through the icosahedral net
- **Adaptive DPR**: Reduces device pixel ratio during interaction for performance, sharpens when settled
- **Pan/Zoom**: Drag to pan, scroll wheel to zoom, throttled via `requestAnimationFrame`
- **Hit Detection**: Color-coded pick buffer maps screen pixels back to cell IDs
- **River Rendering**: Antimeridian crossing detection for correct line wrapping

### Dymaxion Projection

`utils/dymaxion.ts` implements the Buckminster Fuller Dymaxion map:

1. **Icosahedron Geometry**: 12 vertices, 20 triangular faces using the golden ratio
2. **Face Construction**: Oriented faces with correct winding order
3. **2D Net**: Unfolds the icosahedron into a flat layout with barycentric coordinate transforms
4. **D3 Integration**: Creates a `d3.geoPolyhedral` projection using `geoGnomonic` per-face
5. **Orientation**: Configurable lon/lat/roll for rotating the projection center

The Dymaxion projection is available in both the 3D viewer (as a wireframe overlay) and the 2D viewer (as a full raster map).

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

The controller reference is stored in a `useRef` to persist across renders.

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

### Controls Component (1150 lines)

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
        └── generateWorldLore(world)
              └── Gemini API call (async)
                    └── Returns LoreData
                          └── setLore(newLore)
                          └── Mutates world.civData in-place (names)
                          └── setWorld({ ...world })
```

### Map Save/Load
```
Save:
  Controls "Save" → saveMapToBrowser(name, world)
    └── Serialize params + civData → localStorage

Load:
  Controls "Load" → loadMapConfig(file)
    └── Parse JSON → validateWorldParams()
    └── handleLoadWorld(params, savedCivData)
          └── Regenerate world (same seed)
          └── Restore saved names/descriptions
```

---

## Export & Persistence

### Image Export (`utils/export.ts`)

`exportMap()` renders the world to a canvas at configurable resolutions:
- **Resolutions**: 2K (2048px), 4K (4096px), 8K (8160px), 16K (16384px)
- **Projections**: Equirectangular, Mercator, Winkel Tripel, Robinson, Mollweide, Orthographic, Dymaxion
- **Dymaxion raster**: Pixel-by-pixel reprojection from equirectangular source

### Browser Storage

Uses `localStorage` with a `SavedMapEntry[]` structure:
```
{
  name: string,
  timestamp: number,
  params: WorldParams,
  civData?: CivData,
}
```

### JSON Config

Save/load full world configuration as JSON files, including all parameters and civilization metadata.

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
- Environment variable `GEMINI_API_KEY` can be set at build time or provided at runtime
- No server-side rendering; fully static SPA

### Build Configuration

- **Vite 6** with `@vitejs/plugin-react` for HMR
- **TypeScript**: ES2022 target, ESNext modules, `react-jsx` transform
- **Path alias**: `@/*` maps to project root
- **CSP**: HTML meta tag allows self, Tailwind CDN, and Google Generative Language API
