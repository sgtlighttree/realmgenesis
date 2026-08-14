# Export & Persistence

Five output paths — raster image, SVG, GeoJSON, GLB, and save files. The guiding
principle: **geometry regenerates from the seed; only names and edits are
stored.** Save files are tiny and version-tagged.

## Raster image — `utils/export.ts`

`exportMap(world, viewMode, resolution, projection, dymaxionSettings?)` renders
to an offscreen canvas and downloads a PNG.

- **Resolutions** (`ExportResolution`): `2048 | 4096 | 8192`. (16K/32K were
  removed — they exceeded browser canvas limits.)
- **Projections** (`ProjectionType`): equirectangular, mercator, winkeltripel,
  orthographic, robinson, mollweide, dymaxion.
- **Heightmap**: export in `height_bw` at equirectangular for a greyscale DEM
  suitable as a displacement map.
- **Dymaxion** raster reprojects pixel-by-pixel through `buildDymaxionNet` —
  `classic` gives the Fuller net (`width × ~0.6·width`), `blender` gives a square
  UV-space image (net in the lower ~47%) matching Blender's icosphere unwrap.
- Political/culture/religion exports build a `factionColors` map so user color
  edits appear (see [invariants.md](invariants.md) §11).

## Vector — `utils/exportVector.ts`

- `exportSVG(...)` / `downloadSVG(...)` — coastlines, borders, rivers, labels as
  vector for Inkscape/Illustrator (E1). Projection is any
  `VectorProjectionType` (every `ProjectionType` except dymaxion).
- `exportGeoJSON(world)` / `downloadGeoJSON(world)` — cells/coastlines/borders/
  rivers as GeoJSON (E2). The data is genuinely geodesic (lon/lat on a sphere),
  so it opens cleanly in QGIS with no fudging.

Both build from the cached `WorldData.geoJson` + d3-geo path generation.

## 3D — `utils/exportGLB.ts`

`exportGLB(world, viewMode)` builds a fresh Three.js scene (independent of the
live canvas) and downloads a binary GLB:
- `World` — per-vertex-colored `MeshStandardMaterial` (colors as GLTF `COLOR_0`),
  same fan-triangulation and `hMult` elevation as the globe.
- `Rivers` — `LineSegments`.
- `Capitals` / `Towns` — merged cylinders (present only if `civData` exists).

GLB vertex colors need a Blender material step to display — see
[invariants.md](invariants.md) §21.

## Save / load — `utils/export.ts`

- **JSON config**: `saveMapConfig(params, world?)` writes a **version `"1.4"`**
  file containing `params`, `civData`, and `markers` — **not** cell geometry.
  `loadMapConfig(file)` parses, validates, and returns `{ params, civData?,
  markers? }` or `null`.
- **Browser storage**: `saveMapToBrowser` / `getSavedMaps` / `deleteSavedMap`
  persist `SavedMapEntry` records (params + civData + markers) to `localStorage`.
- **Load flow**: params → regenerate world (same seed → same terrain) → restore
  saved names/markers on top.

### Validation
Imports are defended: `validateWorldParams` bounds-checks numeric params;
`validateCivData` shape-checks factions (a corrupt file degrades to terrain-only
rather than crashing); `validateMarkers` sanitizes the marker array. `civData`
and markers are metadata, so a failed check is non-fatal.

> The params validator bounds `tectonicStrength: [0, 2.0]` (renamed from the old
> dead `plateInfluence` key post-Session-13) — see
> [params-reference.md](params-reference.md).
