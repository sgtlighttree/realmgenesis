# Civilization & Worldbuilding Layers

Everything above the terrain: cultures, religions, factions, provinces, towns,
routes, user markers, and AI lore. All of it (except markers) is
**deterministic from the seeds and regenerated, not serialized** — a save file
carries only names and edits (see [export.md](export.md)).

Two seeds drive it: terrain-derived layers (cultures, geography names) key off
`seed`; faction geometry keys off `civSeed`. Re-rolling civs never renames
mountains.

The whole stack is re-runnable without regenerating terrain — the hook exposes
"update civs" and "update provinces" that call `recalculateCivs` /
`recalculateProvinces` directly.

## Order of operations

`recalculateCivs` (also the terminal call of `generateWorld`) runs the chain:

```
recalculateCultures  →  faction expansion  →  recalculateProvinces
                                                 └── recalculateReligions
```

Routes are computed separately and lazily (below).

## Cultures (C1)

`recalculateCultures` in `utils/worldGen.ts`. `numCultures` (default 4) home
cells are seeded, then culture regions spread by **terrain affinity** (forest /
desert / coastal profiles biased per culture) across land cells. Each
`CultureData` carries a `nameStyle` (a namebase from A2). Cultures **predate**
factions: a faction inherits the naming style of the culture its capital sits
in, which is how multi-ethnic states get consistent regional naming.
`cell.cultureId` is set on land cells only.

## Factions

Faction expansion (inside `recalculateCivs`):
1. Place `numFactions` capitals with a scale-independent minimum separation
   (`capitalSpacing`² · 4 / numFactions).
2. Grow territory by **Dijkstra** over terrain costs: `waterCrossingCost`
   multiplies water edges, mountains/deserts add penalties, `borderRoughness`
   injects noise for irregular borders.
3. `civSizeVariance` draws a per-faction size factor; cheaper movement wins the
   competitive frontier race, so larger factors → larger factions.
4. Water within `territorialWaters` graph-distance of land is claimed as
   territorial waters.

`cell.regionId` is the faction id (`POLITICAL_ERASER_ID = -1` clears it).

## Provinces & towns

`recalculateProvinces`: subdivides each faction into provinces sized by
`provinceSize`, places towns, and assigns population by biome suitability
(fertile biomes higher, deserts/tundra lower). Province 0 holds the capital
town. Sets `cell.provinceId`, `isCapital`, `isTown`, `population`.

## Religions (C2)

`recalculateReligions` runs at the end of `recalculateProvinces` (so holy towns
exist). Two kinds:
- **Folk** — one per culture, covering that culture's land by default.
- **Organized** — spread outward from a holy city (a town cell) along the graph
  within a bounded budget, converting folk cells they reach; unreached land
  stays folk.

`cell.religionId` on land cells; `ReligionData.cultureId` is set for folk,
`holyCellId` for organized.

## Routes (C3)

`computeRoutes` in `utils/routes.ts`. A*/Dijkstra land **roads** and sea
**searoutes** between towns over the terrain-cost model. **Not part of the core
pipeline** — computed lazily in `useWorldEngine` when the routes toggle is on,
because it's O(towns · A*). `recalculateProvinces` clears `world.routes` so they
recompute against fresh towns. `RouteData` is never serialized.

## Markers (C4)

User-placed points of interest (`MarkerData`: dungeon / ruin / battlefield /
portal / poi, with a free-text `note`). Unlike everything else here, markers
**are persisted** and are anchored to a sphere `position`, not a `cellId`, so
they survive regeneration (a marker may end up over water after terrain
changes — acceptable). Managed via the marker-mode handlers in the hook.

## Editing (C5)

`utils/civEdit.ts` provides in-place faction/province operations: merge
factions, rename province/town, relocate capital. Faction split is deferred.

## AI lore (optional)

`services/gemini.ts` `generateWorldLore(world)` calls Google Gemini with a
prompt whose depth is set by `loreLevel` (1 names / 2 + provinces & towns / 3 +
backstories). It **mutates `world.civData` in place** — the caller must
`setWorld({ ...world })` to re-render (see [invariants.md](invariants.md) §7).
Lore is an optional enhancer: offline namebases (A2) already name every world
without a key. The key is ephemeral (invariant §19).
