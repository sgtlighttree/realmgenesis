# Rendering

Three surfaces share one coloring function and one world: the 3D globe, the 2D
map, and the export paths. The single source of per-cell color is
`getCellColor(cell, viewMode, seaLevel, factionColors?, cultureColors?,
religionColors?)` in `utils/colors.ts` — change it once and globe, Map2D,
minimap, PNG, and GLB all follow.

The GPU-lifecycle and pick-buffer mechanics that break things if violated live
in [invariants.md](invariants.md) §10–§17; this doc is the map of the surfaces.

## View modes

`getCellColor` switches on `ViewMode` (12): `biome`, `satellite`, `height`,
`height_bw`, `temperature`, `moisture`, `plates`, `political`, `population`,
`province`, `culture`, `religion`. Political / province / culture / religion
modes need their color map passed (`buildFactionColorMap`, and the
`cultureColors` / `religionColors` maps from the hook) or user color edits won't
appear (invariant §11). `PLATE_COLORS` (plates mode) and `FACTION_COLORS`
(everything political) are deliberately separate palettes — don't cross them.

## 3D globe — `components/WorldViewer.tsx`

An R3F (`@react-three/fiber`) scene:
- **World mesh** — one `BufferGeometry` allocated per `world.cells` identity,
  fan-triangulated center→vertices, per-vertex colors refilled in place on
  paint/view change (no per-stroke realloc; no normals; fixed bounding sphere
  r = 1.1 — invariants §13/§14).
- **City markers** — `InstancedMesh` cylinders (capitals red, towns white).
- **Rivers** — smoothed `LineSegments`.
- **Faction borders** — line segments between differing regions; independent of
  the active view layer.
- **Country labels** — curved textured mesh patches just above the surface, so
  they rotate with the planet and are occluded by it (invariants §16). Other
  label kinds (capitals, towns, provinces, geography, markers) toggle via
  `LabelVisibility`.
- **Lat/long grid**, **Dymaxion overlay** (rotatable icosahedron wireframe),
  **Stars** background, **OrbitControls** camera.

Pointer interaction: click inspection raycasts only on click/up (so dragging the
globe doesn't select); paint strokes use native pointer listeners to avoid idle
hover lag.

## 2D map — `components/Map2D.tsx`

Offscreen Canvas2D:
- **Mercator** — `d3.geoMercator` over the cached GeoJSON polygons.
- **Dymaxion** — pixel-by-pixel reprojection of an equirectangular source
  through the icosahedral net (`buildDymaxionNet`).
- **Adaptive DPR** — drops to 1× during interaction for 60fps, sharpens when
  settled.
- **Pick / hit-testing** — Mercator inverts the projection + nearest-cell
  lookup; Dymaxion uses a hidden **color-ID pick buffer** generated through the
  *same* pipeline as the visible raster (invariants §15).
- **Overlays** — borders, faction/geography labels, rivers (with antimeridian
  wrap), and routes draw over any base layer.
- **Seam prevention** — every cell is filled *and* stroked with its own hex
  (invariants §17).

## Derived overlays

- **Hillshade / contours** — `utils/shading.ts` derives a light-vector·normal
  hillshade and elevation isolines from per-cell height (A4); toggled via
  `showHillshade` / `showContours`.
- **Labels** — `utils/labels.ts` holds the label model (importance sizing,
  greedy priority declutter, elongated-feature angling) consumed by both 2D and
  3D (A1).
- **MiniMap** — bottom-right equirectangular overview; builds its own
  `factionColors` (it's outside the React color-map path).
- **Ruler** — `utils/measure.ts` great-circle distance between clicked points,
  scaled by `planetRadius` (A5).
