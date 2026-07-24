# RealmGenesis 3D — Feature Roadmap

This document lays out the feature direction for turning RealmGenesis from a
planetary-simulation viewer into a full fantasy map / planet creation tool.
Azgaar's Fantasy Map Generator is the reference point for scope; the goal is
not feature parity but complementarity — RealmGenesis is sphere-first and
simulation-first, and the roadmap leans into that.

**Where we stand vs. Azgaar.** RealmGenesis already leads on the planetary
side: true spherical Voronoi cells, plate tectonics with stress, hydraulic and
thermal erosion, latitude/tilt-aware climate, Dymaxion and multi-projection
export, GLB export. The gaps are in cartographic *presentation* (nothing on
the map is labeled or styled like a map), worldbuilding depth beyond
factions/provinces/towns, and a handful of simulation features that only make
sense on a sphere — which flat-map tools cannot do correctly at all.

Themes are ordered roughly by leverage. Within each theme, items are ordered
by suggested implementation order. File references point at the current code
that each feature would build on.

---

## A. Cartographic presentation

The biggest gap. A generated world currently *renders* but does not *read* as
a map: no names appear anywhere, and every view mode is a data visualization
rather than a cartographic style.

### A1. Map labels & typography

Render names — factions, towns, rivers, lakes, geographic features — on the
2D map and the 3D globe. Names already exist in `CivData`
(`types.ts`) but are only visible in the Inspector and the 3D curved faction
labels; nothing else is ever labeled.

- 2D (`components/Map2D.tsx`): canvas text drawn in the composite pass, sized
  by importance (faction > ocean > capital > range/sea > town > lake > river),
  with greedy priority decluttering and zoom-based level of detail. Elongated
  features (rivers, ranges) get angled labels along their principal axis;
  per-glyph curved text is a stretch goal.
- 3D (`components/WorldViewer.tsx`): billboarded sprite labels with
  back-of-globe and distance culling. Factions keep the existing curved
  textured-mesh labels.
- Toggleable per label kind; included in PNG export (`utils/export.ts`).
- The Dymaxion point projector currently local to `Map2D.tsx` gets promoted
  into shared geo code so export projects label anchors identically
  (single-source-of-truth invariant, `utils/geo.ts`).

This is the single most transformative feature for "proper fantasy map" feel.

### A2. Offline name generation (namebases)

Seeded Markov-chain / syllable namebases so every world gets names with no
API key. Azgaar's namebase model is the reference: several built-in styles
(generic fantasy, Norse, Latin-imperial, desert, elven, …), each a word list
that trains a small char-level Markov model.

- Deterministic via fresh RNG side-streams from `utils/rng.ts`
  (`seed + '_names'`), so existing seeds keep their terrain and civ geometry.
- Hooks where placeholder names are currently assigned in `recalculateCivs` /
  `recalculateProvinces` (`utils/worldGen.ts`); geographic names key off the
  terrain seed so re-rolling civs doesn't rename mountains.
- Gemini (`services/gemini.ts`) demotes to an optional lore *enhancer* that
  overwrites the offline names when invoked — no longer a requirement for a
  named world.
- New `nameStyle` param + Civ-tab selector.

### A3. Map style system (fantasy rendering)

A style layer over the existing view modes: parchment/hand-drawn preset with
coastline outline strokes, ocean hatching or stipple patterns, paper texture,
and either relief glyphs (mountain/hill/forest icons) or hillshading. Styles
apply to Map2D, the globe texture path, and PNG export. This is what makes
output pinboard-worthy rather than diagnostic.

### A4. Hillshading & contour lines

Cheap wins from existing per-cell height: a hillshade pass (light-vector dot
cell normal) usable in any view, and elevation isolines. Both are derived
data — no pipeline changes.

### A5. Scale bar & great-circle measurement

`planetRadius` is already a parameter, so true geodesic distance is one
formula away. A ruler tool measuring great-circle distance between clicked
points (with an optional travel-time estimate reusing the Dijkstra terrain
costs from civ expansion in `utils/worldGen.ts`), plus a projection-aware
scale bar on exports. Flat-map tools fundamentally cannot do this correctly —
it's a headline geodesic feature.

### A6. Grid overlays

The 10° graticule exists on the globe; add configurable density, 2D-view
graticules, and a geodesic hex grid (icosahedral subdivision) overlay for
hex-crawl RPG use.

---

## B. Physical geography completeness

### B1. Lakes as first-class features

The Priority-Flood depression fill in `generateRivers()`
(`utils/worldGen.ts`) already computes a per-cell `waterLevel`; land cells
where `waterLevel > height` are exactly the filled depressions — currently
discarded. Surfacing them:

- `Lake` entity: member cells, surface level, outflow river or endorheic
  flag, fresh vs. salt (arid + hot basins become endorheic salt lakes).
- New `LAKE` / `SALT_LAKE` biomes; lake cells render as water in all views
  (`utils/colors.ts` is the single coloring function, so one change
  propagates to globe, Map2D, minimap, PNG, and GLB).
- Rivers terminate at lake entry and restart at the outflow.
- Civ integration: no capitals/towns/population on lakes; lake adjacency
  counts as coast.

### B2. River & lake editing / naming

Extend the existing paint/edit system (`utils/paintUtils.ts`) to hydrology:
rename rivers and lakes, delete or redirect river segments, adjust lake
levels. Painting terrain over a lake clears it (snapshot semantics), with a
"resurvey" action to recompute hydrology-derived features.

### B3. Auto-detected & named geographic features

Cluster cells to detect mountain ranges (contiguous high elevation), peaks,
deserts, forests, oceans/seas (water bodies by size), and islands (small land
components); name them via A2 and store anchors (plus a principal axis for
elongated features) for labeling. All O(n) BFS over the existing cell graph —
fine at the 200k-cell cap. Surfaced in the Inspector and consumed by A1.

---

## C. Worldbuilding depth (civilization layer)

### C1. Cultures layer

Distinct from factions: culture regions spread by terrain affinity, each tied
to a namebase from A2. Factions then inherit (and can span) cultures —
enabling believable multi-ethnic states and consistent regional naming.

### C2. Religions layer

Folk religions derived per culture, plus organized religions spreading from
holy cities along the faction graph and (later) trade routes.

### C3. Roads & trade routes

A*/Dijkstra paths between towns over the existing terrain-cost model (reuse
the cost logic in `recalculateCivs`), plus sea routes between ports. Route
connectivity should feed back into town importance/population.

### C4. Markers & notes (POIs)

User-placed pins (dungeon, ruin, battlefield, portal, …) with free-text
notes, persisted in the JSON save; optional Gemini lore per marker.

### C5. Editor completeness

Already tracked in `HANDOFF.md` next tasks: merge/split factions, province
management, bulk rename, town/capital relocation.

### C6. Diplomacy & relations (later)

Inter-faction relations, vassals, alliances, war status — mostly flavor and
map overlays, cheap once C1–C3 exist.

---

## D. Planet-scale simulation (the geodesic differentiators)

Features that exploit the sphere-native architecture — the category where
RealmGenesis can do things no flat-map tool can.

### D1. Seasonal cycle

`axialTilt` currently produces a single static climate. Compute insolation
per orbital position for a season slider: temperature bands, ice extent, and
optionally biome edges shift through the year.

### D2. Ocean currents

Currents deflected by continents and Coriolis, feeding the existing 8-pass
moisture transport and temperature model (warm currents moderating coastal
climates — think Gulf Stream). A major realism upgrade to climate.

### D3. Ice caps & glaciers

Sea-ice extent from temperature (seasonal with D1), glacial land cover, and
optionally sea-level coupling.

### D4. Regional zoom / submap generation

Re-run generation at higher point density inside a selected spherical cap for
local/regional maps — Azgaar's "submap" concept, but projection-correct, with
parent-world values as boundary constraints.

### D5. Planetary parameters

Day length, gravity, star class, moons — mostly lore/export metadata at
first, but D1 (seasons) and tides give them mechanical hooks over time.

### D6.  (added by Matt/Maintainer)

Another overhaul of terrain generaton algorithm to make it more realistic
and get rid of seam lines on plate boundaries, and more detailed heightmap
rendering and calculation without increasing cell count

---

## E. Interoperability & export

### E1. SVG export

Vector export of coastlines, borders, rivers, and labels for post-processing
in Inkscape/Illustrator — a key Azgaar strength and a frequent request for
any mapping tool. The cached GeoJSON cell geometry (`WorldData.geoJson`) plus
d3-geo path generation makes this tractable.

### E2. GeoJSON export

Cells, coastlines, borders, and rivers as GeoJSON. The data is genuinely
geodesic (lon/lat on a sphere), so exports open cleanly in QGIS and web-GIS
tooling with zero fudging — worth advertising.

### E3. Azgaar `.map` import (stretch)

Project an Azgaar flat map onto the sphere (equirectangular assumption,
re-tessellate onto the cell graph). Lossy by nature; stretch goal.

---

## F. Frontend/Rendering/UI Overhaul

# F1. Redesigning the UI/Frontend
The UI has become a mess in an attempt to add more features, especially on mobile.
A full redesign and rearchitecture is warranted here.
Can come *before or alongside* D6 with the roadmap in mind.

# F2. 3D Mode Presentation
Part of redesign is figuring out how the planet is presented; overlays like borders, rivers and roads and routes and the lat/lon grid are also 3D objects, not 2D overlays simply composited over the 3D globe, which affects visibility and accuracy. Or maybe make the globe entirely smooth by default, instead of applying height per cell. Perhaps 3D mode should more like Google Earth Pro in this respect.

# F3. True 2D vector map
Make it a true vector map like most web mapping apps, but keep it as optimized as possible

# F4. Performance Optimizations
Self-explanatory, optimize wherever while keeping visual fidelity. Can come last, but better to make efficient renderers/frontends and code wherever possible.


## Recommended first milestone — "Map identity"

**A1 + A2 + B1 + B3: labels, offline namebases, lakes, and named geographic
features.** Together these turn the simulation viewer into something that
reads as a *map*, and none of them require new simulation machinery — lakes
fall out of the existing depression fill, and labels/names consume data the
pipeline already computes or can derive in one pass.

A detailed implementation design exists for this milestone. Summary:

- **Phasing** (~2–3 weeks total): (1) namebase engine + civ naming hooks;
  (2) hydrology refactor + `Lake` entity + lake biomes/colors; (3) feature
  detection stage + geographic naming + Inspector integration; (4) label
  model + 2D/3D renderers + declutter + export labels. Phases 1 and 2 are
  independent and can land as separate PRs.
- **Data model**: `Lake`, `GeoFeature`, `RiverInfo` types; `Cell.lakeId`;
  `WorldParams.nameStyle`; `LAKE`/`SALT_LAKE` biomes; JSON save schema bumps
  to v1.5 with a names-only `geoNames` container (geometry is never
  serialized — it regenerates deterministically, matching the existing
  save philosophy in `utils/export.ts`).
- **Determinism rules**: every new random draw uses a fresh RNG side-stream
  (`new RNG(seed + '_<purpose>')`) so existing seeds keep their terrain and
  civ geometry byte-identical; geographic names key off the terrain seed,
  civ names off `civSeed`.
- **Top risks**: label decluttering complexity (scoped to greedy
  priority + bbox rejection; curved text deferred), label paths crossing
  Dymaxion net seams (straight anchored labels only there), and 200k-cell
  performance (O(n) typed-array clustering, labels drawn in the cheap
  composite pass, allocation-free 3D culling).

## Suggested order after that

1. **A3 + A4** (map styles, hillshading/contours) — presentation payoff
   compounds with labels.
2. **A5 + E1/E2** (geodesic ruler, SVG/GeoJSON export) — small, high-value,
   and differentiating.
3. **C1 → C3** (cultures → religions → routes) — worldbuilding depth, in
   dependency order.
4. **D1 → D2** (seasons, currents) — the sphere-native simulation showcase.
5. **D4** (submaps) and the rest as demand dictates.
