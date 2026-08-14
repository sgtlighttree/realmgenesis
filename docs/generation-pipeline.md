# Generation Pipeline

`generateWorld(params, onLog?, signal?, onProgress?)` in `utils/worldGen.ts`,
run inside the Web Worker (see [architecture.md](architecture.md)). It reports
**7 progress ticks** (`TOTAL_STAGES = 7`) and returns a fully-populated
`WorldData` including the first civ pass.

The old "12-stage" V2 description (K-means plates, separate connectivity/stress
stages) is gone — this is the current flow.

## Terrain half

```
1. Init grid          Fibonacci sphere (display points, jitter = cellJitter·0.8)
                      → spherical Voronoi (d3 geoVoronoi) → Cell[] + neighbors
                      → WorldData.geoJson (cached polygon FeatureCollection)

2. Simulate tectonics Fibonacci MACRO points (simulationResolution, default 10k)
                      → simulateTectonics(...)  [utils/tectonicsV3.ts]
                      RNG streams: _macro_v3, _crust, _plates_v3
                      → TectonicResult (per-macro heights/crust/uplift/plateId)
                      See tectonics-v3.md

3. Project to display projectTectonicsToDisplay: each display cell takes its
                      nearest macro cell's fields, + structural FBM/ridged noise
                      at full resolution, coastline flattening, optional Pangea
                      mask → normalize to [0,1]

4. Erosion (if >0)    hydraulic then thermal; step counts scale with
                      √(points/5000); heights re-normalized after

5. Stage-9b remap     mountainHeight / oceanDepth / seafloorDepth
                      (only runs if any ≠ default; coastline held fixed)
                      See params-reference.md and invariants.md §2

6. Climate            wind vectors (latitude bands) → 8-pass moisture transport
                      → temperature (latitude² blend + axial tilt + −60·elevation
                      lapse + variance noise) → determineBiome per cell

7. Rivers & lakes     generateRivers: Priority-Flood depression fill (MinHeap)
                      → flux accumulation → path trace to ocean → Catmull-Rom
                      smoothing; filled depressions surface as LakeData (B1)
                      → detectFeatures (named geography, B3)
                      → recalculateCivs (the civ half, below)
```

Determinism comes from isolated RNG side-streams per subsystem — never reorder
draws within a stream (see [invariants.md](invariants.md) §4).

### Neighbor derivation is coordinate-keyed
Cell adjacency is matched from `voronoi.links()` via string keys of lon/lat
rounded to **4 decimals**. Two generator points within ~1e-4° would collide and
mis-wire an edge — unlikely at default density, probability grows toward the
200k cap. Known fragility, not currently a defect.

### Biome classification
`determineBiome(height, temp, moisture, seaLevel)`: below `seaLevel` →
`OCEAN`/`DEEP_OCEAN` (< `seaLevel·0.6`); above, a Köppen-ish decision tree on
normalized land elevation + temperature + moisture, with `BEACH`/`VOLCANIC` as
elevation-band overrides. `LAKE`/`SALT_LAKE` are assigned **later** by the
drainage solve, never here.

## Civilization half

The terrain pipeline ends with `recalculateCivs`, which is **also independently
callable** — the hook re-runs it (and `recalculateProvinces`) without
regenerating terrain. The chain:

```
recalculateCivs(world, params)
  ├── recalculateCultures      terrain-affinity culture regions (C1), RNG-seeded
  │                            from the terrain seed; additive, mutates cells
  ├── faction expansion        capitals placed (capitalSpacing) → Dijkstra
  │                            territory growth over terrain costs (waterCrossingCost,
  │                            borderRoughness, civSizeVariance) → territorial waters
  └── return recalculateProvinces(world, params)
        ├── subdivide factions into provinces (provinceSize) + place towns
        └── recalculateReligions  folk (one per culture) + organized (spread from
                                  holy towns) — runs here so holy cities exist
```

Roads & sea routes (`computeRoutes`, C3) are **not** in the core pipeline —
they're computed lazily in `useWorldEngine` when the routes toggle is on,
because they're O(towns · A*). `recalculateProvinces` clears `world.routes` so
they recompute against fresh towns.

See [civilization.md](civilization.md) for the worldbuilding layers in depth and
[tectonics-v3.md](tectonics-v3.md) for stage 2.
