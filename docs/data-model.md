# Data Model

Every type lives in `types.ts` (except `NameStyle` in `utils/namegen.ts`). This
is the shape of a world. Parameter semantics are in
[params-reference.md](params-reference.md).

## Point & Cell

```ts
Point = { x, y, z }   // 3D Cartesian on the unit sphere
```

A world is `N` cells (default 5,000; cap 200,000). **`cell.id` == its index in
`WorldData.cells`** — stable within a world, not across runs.

```ts
Cell {
  id, center: Point, vertices: Point[], neighbors: number[]   // structure (Stage 2)

  height, plateId, temperature, moisture, biome, flux?        // physical

  crustType?, crustThickness?, upliftAccum?                   // V3 tectonics

  regionId?, provinceId?, isCapital?, isTown?, population?,   // political
  cultureId?, religionId?                                     // cultural (undefined on water)
}
```

`height` is 0–1 normalized (`seaLevel` default 0.55). `crustType` is `0`=oceanic
/ `1`=continental. Optional fields are absent until their pipeline stage runs —
see [generation-pipeline.md](generation-pipeline.md) for which stage sets what.

## BiomeType — 17 values

```
water:       OCEAN, DEEP_OCEAN
polar (E):   ICE_CAP, TUNDRA
dry (B):     HOT_DESERT, COLD_DESERT, STEPPE
tropical (A):TROPICAL_RAINFOREST, TROPICAL_SAVANNA
temperate(C):MEDITERRANEAN, TEMPERATE_FOREST, TEMPERATE_RAINFOREST
continental: BOREAL_FOREST
special:     BEACH, VOLCANIC
hydrology:   LAKE, SALT_LAKE
```

> `LAKE`/`SALT_LAKE` are **derived** — assigned after biome classification from
> the drainage solve, never by `determineBiome`. (Legacy docs list 15 values;
> these two were added with the lakes feature, B1.)

## WorldData

```ts
WorldData {
  cells: Cell[]
  params: WorldParams
  geoJson: GeoJsonCollection   // d3 polygons, cached from Stage 2, feature[i] ↔ cells[i]
  civData?:  CivData           // factions → provinces → towns
  rivers?:   Point[][]         // smoothed river paths
  lakes?:    LakeData[]        // B1
  features?: GeoFeature[]      // named geography, B3
  markers?:  MarkerData[]      // user POIs, C4
  cultures?: CultureData[]     // C1
  religions?:ReligionData[]    // C2
  routes?:   RouteData[]       // C3
}
```

**What is serialized on save vs regenerated:** terrain, rivers, lakes, features,
routes, cultures, and religions all **regenerate deterministically from the
seed** and are *not* stored in save files — only names/edits are (see
[export.md](export.md)). `markers` are the exception: user-placed, so they're
persisted and anchored to a **sphere `position`, not a `cellId`**, to survive
regeneration.

## Political hierarchy

```
CivData.factions: FactionData[]
  ├── id, name, color, capitalId, totalPopulation, description?
  └── provinces: ProvinceData[]
        ├── id (local), name, totalPopulation, color?
        └── towns: TownData[]
              └── name, cellId, population, isCapital
```

## Worldbuilding entities

| Type | Feature | Key fields | Notes |
|------|---------|-----------|-------|
| `CultureData` | C1 | `homeCellId`, `nameStyle`, `color` | Terrain-affinity clusters that predate factions; a faction inherits its capital's culture naming style. |
| `ReligionData` | C2 | `kind: 'folk'\|'organized'`, `cultureId`, `holyCellId` | Folk = one per culture; organized = spreads from a holy town within a budget. |
| `LakeData` | B1 | `cellIds`, `surfaceLevel`, `outflowCellId`, `isEndorheic`, `isSalt` | A filled depression from the Priority-Flood solve; salt = hot+arid basin. |
| `GeoFeature` | B3 | `kind` (range/desert/forest/sea/ocean/island/lake), `anchor`, `size` | Terrain-derived clusters, named from the **terrain** seed; anchor is the label position. |
| `MarkerData` | C4 | `kind`, `note`, `position: Point` | User POI; sphere-anchored, persisted. |
| `RouteData` | C3 | `path: Point[]`, `kind: 'road'\|'searoute'`, from/to | Town-to-town over the terrain cost model; not serialized. |

## Modes & enums

```ts
DisplayMode = 'globe' | 'mercator' | 'dymaxion'
ViewMode    = 'biome' | 'height' | 'height_bw' | 'temperature' | 'moisture'
            | 'plates' | 'political' | 'population' | 'province' | 'satellite'
            | 'culture' | 'religion'                         // 12 modes
InspectMode = 'click' | 'hover' | 'off'
EditMode    = 'off' | 'terrain-raise' | 'terrain-lower' | 'terrain-flatten'
            | 'terrain-smooth' | 'biome' | 'political' | 'world-edit'
PaintStyle  = 'adaptive' | 'freeform'
```

> `ViewMode` gained `culture` and `religion` (legacy docs list 10). Each needs
> its color map passed to `getCellColor` — see [invariants.md](invariants.md) §11.

## Supporting types

- `GeoJsonFeature` / `GeoJsonCollection` — minimal d3-compatible shapes;
  `features[i]` corresponds to `cells[i]`; cached for the world's lifetime.
- `TectonicResult` — the V3 sim output (typed arrays per **macro**-cell),
  projected onto display cells in one pass. See [tectonics-v3.md](tectonics-v3.md).
- `UndoSnapshot` — `Map<cellId, {height, biome, regionId?, provinceId?}>`, the
  grow-in-place stroke record ([invariants.md](invariants.md) §8).
- `LabelVisibility` — per-kind label toggles (factions, capitals, towns,
  provinces, borders, geography, markers); `DEFAULT_LABEL_VISIBILITY` sets the
  initial state.
- `DymaxionSettings` — `layout` (classic/blender), `lon`/`lat`/`roll`, overlay
  toggle, control `mode`.
- `LoreData` — `{ name, description }` returned by the Gemini lore service.
