# C3 — Roads & Trade Routes (design)

**Date:** 2026-07-24
**Status:** Approved for planning
**Roadmap item:** C3 (last pre-D6 additive feature)

## Goal

Add land roads and sea trade routes between towns, derived deterministically
from the existing civilization + terrain data, rendered like rivers across all
three views (3D globe, 2D Mercator, 2D Dymaxion) and PNG/SVG export, behind a
single toggle. Purely additive: no change to terrain, civ, culture, or religion
geometry for any existing seed.

## Non-goals (deferred, on purpose)

- **Route → town-importance feedback loop.** ROADMAP wants route connectivity to
  raise town population/importance, but that makes C3 non-additive (it feeds
  back into civ generation). Deferred. `RouteData.fromCellId/toCellId` are stored
  now so it becomes a small future step.
- New tunable `WorldParams`. Route density is governed by module constants for
  this feature; promoting them to sliders is future work (and would belong to F1
  UI reorganization, not here).

## Data model (`types.ts`)

```ts
export interface RouteData {
  path: Point[];               // unit-sphere polyline, same shape rivers use
  kind: 'road' | 'searoute';
  fromCellId: number;          // town endpoints (cheap; unlocks deferred feedback)
  toCellId: number;
}

// WorldData:
routes?: RouteData[];
```

Routes are **derived, never serialized** — identical philosophy to `rivers` and
`lakes`. They regenerate deterministically from seed + civData, so **the JSON
save schema does not change**. This is a deliberate constraint: adding C3 must
not touch the save format.

## Generation — `utils/routes.ts` (new)

Entry point: `computeRoutes(world: WorldData, params: WorldParams): RouteData[]`.

### Shared cost extraction (prerequisite)

Extract the per-step traversal cost currently inline in `recalculateCivs`
(`worldGen.ts:1458-1470`) into an exported `cellTraversalCost(from, to, params)`
so roads and civ expansion share one source of truth. The land cost keeps the
existing terms: base `landCost`, `slope * 20`, biome multipliers (HOT_DESERT ×2,
ICE_CAP ×4, VOLCANIC ×5). Water is impassable for the land cost.

> Note: the civ pass also multiplies cost by `(1 + civRng.next() * borderRoughness)`
> and by a per-faction `costMult`. Those are civ-expansion-specific and must
> **not** be part of the shared deterministic `cellTraversalCost` used for roads
> (roads have no per-faction cost scaling and draw no RNG). Extract only the
> terrain terms; the civ pass keeps applying its own multipliers on top.

### Pass A — Roads (land network, a forest)

The road network is a **forest**: one tree per land-reachable town cluster. This
is the fix for the great-circle-MST-strands-mainland-towns bug (see Risks).

1. **Land-connected components.** BFS/union-find over the land cell graph
   (cells with `height >= seaLevel` and not lake cells) to label each land cell
   with a component id. O(n), reuses the existing `neighbors` adjacency.
2. **Group towns** by `(factionId, landComponentId)`. Each group is the set of
   same-faction towns that are mutually land-reachable.
3. **MST per group** using great-circle distance between town cells as the edge
   weight (avoids O(towns²) A\* runs). Edge selection is deterministic: sort
   candidate edges by `(weight, minCellId, maxCellId)` so equal-weight ties break
   by cell id.
4. **A\* per MST edge** (great-circle heuristic × min step cost) over the land
   cost field to produce the actual road polyline. Because both endpoints are in
   the same land component, a land path is guaranteed to exist. Store as
   `kind: 'road'`.
5. **Capital trunk roads:** for each pair of bordering factions (adjacency read
   from existing `regionId` neighbor relationships), if both capitals are in the
   same land component, add an A\* road capital↔capital. Cross-water capital
   pairs are left to the sea layer.

### Pass B — Sea routes (maritime web)

Uses a **sea cost field**: water cells cheap (uniform), land impassable except
the port endpoint cells.

1. **Ports** = coastal towns (reuse the `coastalIds` set precomputed at
   `worldGen.ts:1005`). **Major** ports = top fraction by population (module
   constant, e.g. top 40%, with a small floor so tiny worlds still get routes).
2. For each major port, take its **nearest ~3 major ports** by great-circle,
   collect as undirected candidate edges, dedup. (User chose the dense "all major
   ports" web.)
3. **A\*** each candidate over the sea cost field → `kind: 'searoute'` polyline.
   Drop any edge with no sea path within the expansion cap.

### Determinism

- No RNG is required for MST or A\*. A tie-break stream `new RNG(civSeed +
  '_routes')` is available if any future tie-break needs it, but the primary
  determinism guarantees are: (a) MST edge sort breaks ties by cell id;
  (b) A\* relies on deterministic neighbor iteration order, exactly as the
  existing civ Dijkstra does (already guarded by determinism tests).
- Because routes derive only from towns + terrain, re-rolling civs regenerates
  routes and **terrain-only param changes leave civ geometry byte-identical**.

### Wiring

`computeRoutes` is called at the **tail of `recalculateProvinces`** (towns must
exist first). Verified: `generateWorld` → `recalculateCivs` → `return
recalculateProvinces(...)`, so fresh generation *and* standalone civ recalcs both
reach it. `world.routes` is assigned there.

## Rendering (mirrors rivers exactly)

- **3D (`WorldViewer.tsx`):** new `RouteLines` component, keyed on `world.routes`
  (stable across paint strokes), batched into `LineSegments`. Two materials:
  roads solid tan; sea routes **dashed** teal via `LineDashedMaterial` +
  `computeLineDistances()` (the one bit of new render machinery). Rendered at a
  radius just above the surface, like rivers. Matching disposal effect
  (per the invariant: every `useMemo` geometry has a disposal effect).
- **2D (`Map2D.tsx`):** drawn in the composite pass on both Mercator and
  Dymaxion, roads solid / sea dashed, reusing the river antimeridian-jump
  polyline handling.
- **Export (`utils/export.ts`):** mirrored in both raster export paths like
  rivers/contours; included in SVG export (E1) as styled polyline groups.
- **Toggle:** a single `showRoutes` boolean threaded App → Controls ("Map
  Overlays" group) → WorldViewer / Map2D / export. **One** toggle, not two:
  the pickup notes direct that C-tier UI be kept minimal because F1 will
  reorganize it. Default off.

## Tests (`tests/routes.test.ts`)

1. **Determinism:** same seed → identical `routes` (deep-equal on paths + kinds).
2. **Purity:** `computeRoutes` never mutates `height`, `biome`, `regionId`,
   `provinceId`, or `civData`.
3. **Forest invariant (scoped to MST edges only):** considering only the
   `kind: 'road'` edges produced by the per-group MST (i.e. *excluding* capital
   trunk roads), any two same-faction towns in the same land component are
   connected, and that MST subgraph is acyclic (edge count == towns − number of
   `(faction, component)` groups). Capital trunk roads are inter-faction edges
   that may add inter-group cycles and are therefore excluded from this test —
   they are covered by a separate assertion that every trunk endpoint is a
   capital in a bordering faction sharing the same land component.
4. **Sea bridging:** towns in different land components of the same faction are
   *not* required to be road-connected (they rely on sea routes) — the negative
   of test 3, so the two can't contradict.

Extend `paramLiveness` only if a new `WorldParams` key is introduced — it is not,
so no liveness change is required.

## Risks & mitigations

- **[Resolved in design] Great-circle MST stranding land towns.** A per-faction
  great-circle MST could pick A–island–B, then A\* drops both water edges,
  disconnecting mainland towns A and B and contradicting the connectivity test.
  Fixed by the land-component forest: MST is built per land component, so every
  MST edge is guaranteed land-traversable and the invariant is both true and
  testable.
- **A\* count at the 200k-cell cap.** MST keeps land A\* runs to ~(towns − trees);
  sea to ~3× major ports. Mitigated by the great-circle heuristic (directed
  search, most cells never expanded), a per-edge expansion cap, and skipping
  edges whose great-circle length already exceeds a max-route-length constant.
- **[Tuning knob, non-blocking] Coastal sea clutter.** "Nearest 3 major ports"
  can pick a port a few cells down the same coast, drawing tiny sea hops that
  parallel the coast road. If it reads as clutter, dedup sea candidates against
  already-road-connected pairs or impose a minimum crossing distance. The user
  chose the dense web, so this is post-hoc tuning, not a design change.

## Files touched

| File | Change |
|---|---|
| `types.ts` | `RouteData` interface; `WorldData.routes?`. |
| `utils/routes.ts` | **NEW** — `computeRoutes`, land-component BFS, per-component MST, A\*, sea pass. |
| `utils/worldGen.ts` | Extract `cellTraversalCost`; call `computeRoutes` at tail of `recalculateProvinces`. |
| `components/WorldViewer.tsx` | `RouteLines` component + disposal; `showRoutes` prop. |
| `components/Map2D.tsx` | Draw routes on Mercator + Dymaxion; `showRoutes` prop. |
| `components/Controls.tsx` | `showRoutes` checkbox in Map Overlays. |
| `App.tsx` | `showRoutes` state, prop-drilled. |
| `utils/export.ts` | Routes in raster + SVG export. |
| `tests/routes.test.ts` | **NEW** — determinism, purity, forest invariant, sea bridging. |

## Gates

typecheck 0, lint 0 errors at the **30-warning ratchet (add none)**, all tests
green, build OK. Browser-verify routes on all three views + a PNG export before
commit. Do not push (standing rule).
