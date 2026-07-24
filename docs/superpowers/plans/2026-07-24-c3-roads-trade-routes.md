# C3 — Roads & Trade Routes Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking. **Execution is inline/self by maintainer directive — do NOT delegate to subagents.**

**Goal:** Add land roads (a per-land-component forest of town MSTs plus capital trunks) and a dense sea-route web between major ports, derived deterministically from existing civ+terrain data and rendered across all three views + export behind one toggle.

**Architecture:** A new `utils/routes.ts` computes routes from the finished `world` at the tail of `recalculateProvinces`. Path cost logic is extracted from the civ Dijkstra into a shared `utils/pathfinding.ts` (also home to the moved `MinHeap` and a new `aStar`). Routes are derived, never serialized — same philosophy as `rivers`. Rendering mirrors `RiverLines` exactly, roads solid / sea dashed.

**Tech Stack:** TypeScript (strict), React 19, Three.js/R3F, Canvas2D, Vitest.

## Global Constraints

- **Relative imports only** — no `@/` alias (copy verbatim from CLAUDE.md).
- **Gates:** `npm run typecheck` = 0 errors; `npm run lint` = 0 errors at the **30-warning ratchet, add none**; `npm test` all green; `npm run build` OK.
- **Determinism:** existing seeds must keep byte-identical terrain + civ geometry. New random draws (none required here) would use a fresh `new RNG(seed + '_<purpose>')` side-stream.
- **No save-schema change** — routes are derived and not serialized.
- **Rendering is verified manually in the browser** (per CLAUDE.md); engine logic is Vitest-tested.
- **Do NOT push.** Commit locally in focused chunks.
- **Do NOT touch** `CLAUDE.md` / `ROADMAP.md`; `HANDOFF.md` is updated only in the final task.

---

### Task 1: Shared pathfinding module + cost extraction

Extract the reusable primitives so roads and the civ pass share one cost source, with **byte-identical** civ output (guarded by the existing determinism suite).

**Files:**
- Create: `utils/pathfinding.ts`
- Modify: `utils/worldGen.ts` (move `MinHeap` out; replace inline civ step cost)
- Modify: `types.ts` (add `RouteData`, `WorldData.routes`)
- Test: `tests/pathfinding.test.ts`

**Interfaces:**
- Consumes: `Cell`, `WorldParams`, `BiomeType` from `types.ts`.
- Produces:
  - `class MinHeap<T>` (constructor `(scoreFunction: (t: T) => number)`, `push`, `pop`, `size`) — moved verbatim from `worldGen.ts`.
  - `isWaterCell(cell: Cell, seaLevel: number): boolean`
  - `landTerrainStepCost(from: Cell, to: Cell): number`
  - `aStar(cells: Cell[], startId: number, goalId: number, stepCost: (fromId: number, toId: number) => number, heuristic: (id: number) => number, maxExpansions: number): number[] | null`
  - `types.ts`: `interface RouteData { path: Point[]; kind: 'road' | 'searoute'; fromCellId: number; toCellId: number; }` and `WorldData.routes?: RouteData[]`.

- [ ] **Step 1: Write the failing test** (`tests/pathfinding.test.ts`)

```ts
import { describe, it, expect } from 'vitest';
import { isWaterCell, landTerrainStepCost, aStar, MinHeap } from '../utils/pathfinding';
import { BiomeType, type Cell } from '../types';

const cell = (over: Partial<Cell>): Cell =>
  ({ id: 0, center: { x: 0, y: 0, z: 1 }, height: 0.5, neighbors: [], biome: BiomeType.GRASSLAND, ...over }) as Cell;

describe('pathfinding primitives', () => {
  it('isWaterCell: below sea level or a lake biome is water', () => {
    expect(isWaterCell(cell({ height: 0.1 }), 0.3)).toBe(true);
    expect(isWaterCell(cell({ height: 0.5, biome: BiomeType.LAKE }), 0.3)).toBe(true);
    expect(isWaterCell(cell({ height: 0.5, biome: BiomeType.SALT_LAKE }), 0.3)).toBe(true);
    expect(isWaterCell(cell({ height: 0.5 }), 0.3)).toBe(false);
  });

  it('landTerrainStepCost: base 1 + slope*20, biome multipliers on the target', () => {
    const flat = cell({ height: 0.5 });
    expect(landTerrainStepCost(flat, cell({ height: 0.5 }))).toBeCloseTo(1);
    expect(landTerrainStepCost(flat, cell({ height: 0.55 }))).toBeCloseTo(1 + 0.05 * 20);
    expect(landTerrainStepCost(flat, cell({ height: 0.5, biome: BiomeType.HOT_DESERT }))).toBeCloseTo(2);
    expect(landTerrainStepCost(flat, cell({ height: 0.5, biome: BiomeType.ICE_CAP }))).toBeCloseTo(4);
    expect(landTerrainStepCost(flat, cell({ height: 0.5, biome: BiomeType.VOLCANIC }))).toBeCloseTo(5);
  });

  it('aStar: finds a path on a simple line graph and returns null past the expansion cap', () => {
    // 0 - 1 - 2 - 3 line
    const cells: Cell[] = [0, 1, 2, 3].map(i =>
      cell({ id: i, height: 0.5, neighbors: [i - 1, i + 1].filter(n => n >= 0 && n <= 3) }));
    const path = aStar(cells, 0, 3, () => 1, () => 0, 100);
    expect(path).toEqual([0, 1, 2, 3]);
    expect(aStar(cells, 0, 3, () => 1, () => 0, 1)).toBeNull();
  });
});
```

- [ ] **Step 2: Run to verify it fails**

Run: `npx vitest run tests/pathfinding.test.ts`
Expected: FAIL — cannot resolve `../utils/pathfinding`.

- [ ] **Step 3: Create `utils/pathfinding.ts`**

Copy the existing `MinHeap` class body verbatim from `utils/worldGen.ts` (lines ~10-60) into this file and `export` it, then add the helpers:

```ts
import { BiomeType, type Cell, type WorldParams } from '../types';

export class MinHeap<T> {
  // ... paste the exact existing MinHeap implementation from worldGen.ts, unchanged, with `export` ...
}

export function isWaterCell(cell: Cell, seaLevel: number): boolean {
  return cell.height < seaLevel || cell.biome === BiomeType.LAKE || cell.biome === BiomeType.SALT_LAKE;
}

// Terrain-only step cost, extracted from recalculateCivs. Base land cost 1,
// biome multipliers on the destination, slope penalty. NO water override and
// NO civ-specific multipliers (borderRoughness / costMult) — callers add those.
export function landTerrainStepCost(from: Cell, to: Cell): number {
  let c = 1;
  if (to.biome === BiomeType.ICE_CAP) c *= 4;
  if (to.biome === BiomeType.HOT_DESERT) c *= 2;
  if (to.biome === BiomeType.VOLCANIC) c *= 5;
  c += Math.abs(to.height - from.height) * 20;
  return c;
}

// Directed search. stepCost returns Infinity for impassable edges. heuristic is
// a per-cell lower-ish bound on remaining cost (accelerator; routes need only be
// good and deterministic, not provably shortest). Returns cell-id path or null.
export function aStar(
  cells: Cell[],
  startId: number,
  goalId: number,
  stepCost: (fromId: number, toId: number) => number,
  heuristic: (id: number) => number,
  maxExpansions: number,
): number[] | null {
  const open = new MinHeap<{ id: number; f: number }>(x => x.f);
  const g = new Map<number, number>();
  const cameFrom = new Map<number, number>();
  g.set(startId, 0);
  open.push({ id: startId, f: heuristic(startId) });
  const closed = new Set<number>();
  let expansions = 0;

  while (open.size() > 0) {
    const { id } = open.pop()!;
    if (id === goalId) {
      const path: number[] = [id];
      let cur = id;
      while (cameFrom.has(cur)) { cur = cameFrom.get(cur)!; path.push(cur); }
      return path.reverse();
    }
    if (closed.has(id)) continue;
    closed.add(id);
    if (++expansions > maxExpansions) return null;
    // Deterministic neighbor order (neighbors arrays are already stable).
    for (const nId of cells[id].neighbors) {
      const step = stepCost(id, nId);
      if (!isFinite(step)) continue;
      const tentative = g.get(id)! + step;
      if (!g.has(nId) || tentative < g.get(nId)!) {
        g.set(nId, tentative);
        cameFrom.set(nId, id);
        open.push({ id: nId, f: tentative + heuristic(nId) });
      }
    }
  }
  return null;
}
```

- [ ] **Step 4: Refactor `worldGen.ts` to consume the shared module**

Remove the local `MinHeap` class definition. Add at the top with the other util imports:

```ts
import { MinHeap, landTerrainStepCost } from './pathfinding';
```

In `recalculateCivs`, replace the inline biome/slope block (`worldGen.ts:1458-1466`):

```ts
let moveCost = landCost;
if (nCell.biome === BiomeType.ICE_CAP) moveCost *= 4;
if (nCell.biome === BiomeType.HOT_DESERT) moveCost *= 2;
if (nCell.biome === BiomeType.VOLCANIC) moveCost *= 5;
const slope = Math.abs(nCell.height - currCell.height);
moveCost += slope * 20;
const isWater = nCell.height < params.seaLevel || isLakeCell(nCell);
if (isWater) moveCost = waterCost;
```

with:

```ts
let moveCost = landTerrainStepCost(currCell, nCell);
const isWater = nCell.height < params.seaLevel || isLakeCell(nCell);
if (isWater) moveCost = waterCost;
```

Leave the following `moveCost *= (1 + civRng.next() * params.borderRoughness)` and `moveCost *= costMult[region]` lines untouched. (The `landCost` const may now be unused — delete it if lint flags it.)

Add to `types.ts`:

```ts
export interface RouteData {
  path: Point[];
  kind: 'road' | 'searoute';
  fromCellId: number;
  toCellId: number;
}
```

and inside `WorldData`, next to `rivers?`:

```ts
routes?: RouteData[];
```

- [ ] **Step 5: Run the new test + the full suite**

Run: `npx vitest run tests/pathfinding.test.ts`
Expected: PASS.
Run: `npm test`
Expected: ALL PASS — the generation-determinism tests confirm the civ refactor is byte-identical.

- [ ] **Step 6: typecheck + lint**

Run: `npm run typecheck && npm run lint`
Expected: 0 errors; warnings ≤ 30.

- [ ] **Step 7: Commit**

```bash
git add utils/pathfinding.ts utils/worldGen.ts types.ts tests/pathfinding.test.ts
git commit -m "Extract shared pathfinding + RouteData type"
```

---

### Task 2: Roads pass — land-component forest + capital trunks

**Files:**
- Create: `utils/routes.ts`
- Modify: `utils/worldGen.ts` (call `computeRoutes` at tail of `recalculateProvinces`)
- Test: `tests/routes.test.ts`

**Interfaces:**
- Consumes: `MinHeap`, `isWaterCell`, `landTerrainStepCost`, `aStar` from `./pathfinding`; `RouteData`, `WorldData`, `WorldParams`, `Cell`, `Point` from `../types`.
- Produces: `computeRoutes(world: WorldData, params: WorldParams): RouteData[]` (roads only this task; sea added in Task 3). Internal (not exported): `landComponentIds`, town gathering, Kruskal MST.

- [ ] **Step 1: Write the failing tests** (`tests/routes.test.ts`)

```ts
import { describe, it, expect } from 'vitest';
import { generateWorld } from '../utils/worldGen';
import { DEFAULT_PARAMS } from '../types';

async function build(seed = 'route-test') {
  const params = { ...DEFAULT_PARAMS, seed, civSeed: seed };
  return generateWorld(params);
}

describe('computeRoutes — roads', () => {
  it('is deterministic: same seed → identical routes', async () => {
    const a = await build();
    const b = await build();
    expect(JSON.stringify(a.routes)).toEqual(JSON.stringify(b.routes));
  });

  it('produces at least some road routes with town endpoints', async () => {
    const w = await build();
    const roads = (w.routes ?? []).filter(r => r.kind === 'road');
    expect(roads.length).toBeGreaterThan(0);
    for (const r of roads) {
      expect(r.path.length).toBeGreaterThanOrEqual(2);
      expect(r.fromCellId).toBeGreaterThanOrEqual(0);
      expect(r.toCellId).toBeGreaterThanOrEqual(0);
    }
  });

  it('purity: computing routes does not mutate terrain or civ fields', async () => {
    const w = await build();
    const before = w.cells.map(c => `${c.height}|${c.biome}|${c.regionId}|${c.provinceId}`).join(',');
    // recompute (idempotent) — routes are pure over the finished world
    const { computeRoutes } = await import('../utils/routes');
    computeRoutes(w, { ...DEFAULT_PARAMS, seed: 'route-test', civSeed: 'route-test' });
    const after = w.cells.map(c => `${c.height}|${c.biome}|${c.regionId}|${c.provinceId}`).join(',');
    expect(after).toEqual(before);
  });
});
```

> If `DEFAULT_PARAMS` is not the exact export name, use the constant the other tests import (grep `tests/generation.test.ts` for the params factory and match it).

- [ ] **Step 2: Run to verify it fails**

Run: `npx vitest run tests/routes.test.ts`
Expected: FAIL — `w.routes` is undefined / `../utils/routes` unresolved.

- [ ] **Step 3: Create `utils/routes.ts` (roads)**

```ts
import { MinHeap, isWaterCell, landTerrainStepCost, aStar } from './pathfinding';
import type { Cell, Point, RouteData, WorldData, WorldParams } from '../types';

const MAX_ROAD_EXPANSIONS = 20000;   // per-edge A* cap
const MAX_ROAD_ANGLE = 1.2;          // skip MST edges longer than this (radians)

function angle(a: Point, b: Point): number {
  const d = a.x * b.x + a.y * b.y + a.z * b.z;
  return Math.acos(Math.max(-1, Math.min(1, d)));
}

// BFS label every land cell with a component id (-1 for water).
function landComponentIds(world: WorldData, params: WorldParams): number[] {
  const comp = new Array<number>(world.cells.length).fill(-1);
  let next = 0;
  for (const c of world.cells) {
    if (comp[c.id] !== -1 || isWaterCell(c, params.seaLevel)) continue;
    const id = next++;
    const stack = [c.id];
    comp[c.id] = id;
    while (stack.length) {
      const cur = stack.pop()!;
      for (const nId of world.cells[cur].neighbors) {
        const n = world.cells[nId];
        if (comp[nId] === -1 && !isWaterCell(n, params.seaLevel)) { comp[nId] = id; stack.push(nId); }
      }
    }
  }
  return comp;
}

interface Town { cellId: number; factionId: number; population: number; isCapital: boolean; }

function gatherTowns(world: WorldData): Town[] {
  const towns: Town[] = [];
  const factions = world.civData?.factions ?? [];
  factions.forEach((f, fi) => {
    for (const prov of f.provinces) {
      for (const t of prov.towns) {
        towns.push({ cellId: t.cellId, factionId: fi, population: t.population, isCapital: t.isCapital });
      }
    }
  });
  return towns;
}

// Kruskal MST over the given town indices, weighted by great-circle distance.
// Deterministic: edges sorted by (weight, minCellId, maxCellId).
function mstEdges(towns: Town[], idx: number[], cells: Cell[]): Array<[number, number]> {
  const edges: Array<{ a: number; b: number; w: number }> = [];
  for (let i = 0; i < idx.length; i++) {
    for (let j = i + 1; j < idx.length; j++) {
      const a = towns[idx[i]].cellId, b = towns[idx[j]].cellId;
      edges.push({ a, b, w: angle(cells[a].center, cells[b].center) });
    }
  }
  edges.sort((e1, e2) =>
    e1.w - e2.w || Math.min(e1.a, e1.b) - Math.min(e2.a, e2.b) || Math.max(e1.a, e1.b) - Math.max(e2.a, e2.b));
  const parent = new Map<number, number>();
  const find = (x: number): number => { while (parent.get(x)! !== x) { parent.set(x, parent.get(parent.get(x)!)!); x = parent.get(x)!; } return x; };
  for (const e of edges) { parent.set(e.a, e.a); parent.set(e.b, e.b); }
  const out: Array<[number, number]> = [];
  for (const e of edges) {
    const ra = find(e.a), rb = find(e.b);
    if (ra !== rb) { parent.set(ra, rb); out.push([e.a, e.b]); }
  }
  return out;
}

function landHeuristic(cells: Cell[], goalId: number, stepsPerRadian: number): (id: number) => number {
  const goal = cells[goalId].center;
  return (id: number) => angle(cells[id].center, goal) * stepsPerRadian;
}

export function computeRoutes(world: WorldData, params: WorldParams): RouteData[] {
  const cells = world.cells;
  const comp = landComponentIds(world, params);
  const towns = gatherTowns(world);
  const routes: RouteData[] = [];

  // Very rough lower bound on cells crossed per radian, for the A* heuristic.
  const stepsPerRadian = Math.max(1, Math.sqrt(cells.length / (4 * Math.PI)));

  const landStep = (fromId: number, toId: number): number =>
    isWaterCell(cells[toId], params.seaLevel) ? Infinity : landTerrainStepCost(cells[fromId], cells[toId]);

  // Roads: one MST per (faction, land-component) group.
  const groups = new Map<string, number[]>();
  towns.forEach((t, i) => {
    const key = `${t.factionId}:${comp[t.cellId]}`;
    if (comp[t.cellId] < 0) return; // town on water (shouldn't happen) — skip
    (groups.get(key) ?? groups.set(key, []).get(key)!).push(i);
  });

  for (const idx of groups.values()) {
    if (idx.length < 2) continue;
    for (const [a, b] of mstEdges(towns, idx, cells)) {
      if (angle(cells[a].center, cells[b].center) > MAX_ROAD_ANGLE) continue;
      const path = aStar(cells, a, b, landStep, landHeuristic(cells, b, stepsPerRadian), MAX_ROAD_EXPANSIONS);
      if (path && path.length >= 2) {
        routes.push({ path: path.map(id => cells[id].center), kind: 'road', fromCellId: a, toCellId: b });
      }
    }
  }

  // Capital trunk roads between bordering factions in the same land component.
  const capitals = towns.filter(t => t.isCapital);
  for (let i = 0; i < capitals.length; i++) {
    for (let j = i + 1; j < capitals.length; j++) {
      const A = capitals[i], B = capitals[j];
      if (A.factionId === B.factionId) continue;
      if (comp[A.cellId] < 0 || comp[A.cellId] !== comp[B.cellId]) continue;
      if (angle(cells[A.cellId].center, cells[B.cellId].center) > MAX_ROAD_ANGLE) continue;
      const path = aStar(cells, A.cellId, B.cellId, landStep, landHeuristic(cells, B.cellId, stepsPerRadian), MAX_ROAD_EXPANSIONS);
      if (path && path.length >= 2) {
        routes.push({ path: path.map(id => cells[id].center), kind: 'road', fromCellId: A.cellId, toCellId: B.cellId });
      }
    }
  }

  return routes;
}
```

> The `groups.get(key) ?? groups.set(...).get(key)!` idiom is terse; if lint objects, expand to an explicit `if (!groups.has(key)) groups.set(key, [])`.

- [ ] **Step 4: Wire into the pipeline** (`utils/worldGen.ts`)

Add import near the top:

```ts
import { computeRoutes } from './routes';
```

At the very end of `recalculateProvinces`, immediately before its `return world;`:

```ts
world.routes = computeRoutes(world, params);
```

- [ ] **Step 5: Run tests**

Run: `npx vitest run tests/routes.test.ts`
Expected: PASS (determinism, non-empty roads, purity).
Run: `npm test`
Expected: ALL PASS (no existing signature drift).

- [ ] **Step 6: typecheck + lint**

Run: `npm run typecheck && npm run lint`
Expected: 0 errors; warnings ≤ 30.

- [ ] **Step 7: Commit**

```bash
git add utils/routes.ts utils/worldGen.ts tests/routes.test.ts
git commit -m "Add roads pass: per-land-component town MST forest + trunks"
```

---

### Task 3: Sea routes pass — dense major-port web

**Files:**
- Modify: `utils/routes.ts`
- Test: `tests/routes.test.ts` (extend)

**Interfaces:**
- Consumes: everything from Task 2 plus `world.civData` town populations.
- Produces: `computeRoutes` now also appends `kind: 'searoute'` entries. No new exports.

- [ ] **Step 1: Write the failing test** (append to `tests/routes.test.ts`)

```ts
describe('computeRoutes — sea', () => {
  it('produces sea routes between coastal towns, dashed-kind, distinct from roads', async () => {
    const w = await build('sea-test');
    const sea = (w.routes ?? []).filter(r => r.kind === 'searoute');
    // Not every world has multiple landmasses; assert structure when present.
    for (const r of sea) {
      expect(r.path.length).toBeGreaterThanOrEqual(2);
      expect(r.fromCellId).not.toEqual(r.toCellId);
    }
    // Default-ish worlds have oceans + multiple coasts; expect at least one.
    expect(sea.length).toBeGreaterThan(0);
  });
});
```

> If the chosen seed yields a single landmass (no sea routes), pick a seed in the test that has oceans separating coasts — verify by logging `w.lakes`/ocean coverage during Step 2, then pin that seed.

- [ ] **Step 2: Run to verify it fails**

Run: `npx vitest run tests/routes.test.ts -t sea`
Expected: FAIL — `sea.length` is 0.

- [ ] **Step 3: Implement the sea pass** (in `utils/routes.ts`)

Add constants near the top:

```ts
const MAX_SEA_EXPANSIONS = 40000;
const PORT_TOP_FRACTION = 0.4;   // "major" ports = top 40% by population
const SEA_NEIGHBORS = 3;         // each major port links to nearest ~3 ports
const SEA_STEP = 1;              // uniform water step cost
```

Add a coastal-town detector and the sea pass. Insert before the final `return routes;` in `computeRoutes`:

```ts
  // --- Sea routes: dense web over major ports (coastal towns). ---
  const coastal = new Set<number>();
  for (const c of cells) {
    if (isWaterCell(c, params.seaLevel)) continue;
    for (const nId of c.neighbors) {
      if (isWaterCell(cells[nId], params.seaLevel)) { coastal.add(c.id); break; }
    }
  }
  const ports = towns.filter(t => coastal.has(t.cellId));
  if (ports.length >= 2) {
    const sorted = [...ports].sort((a, b) => b.population - a.population || a.cellId - b.cellId);
    const majorCount = Math.max(2, Math.ceil(sorted.length * PORT_TOP_FRACTION));
    const major = sorted.slice(0, majorCount);

    const seaStep = (fromId: number, toId: number): number => {
      const to = cells[toId];
      if (isWaterCell(to, params.seaLevel)) return SEA_STEP;
      return majorSet.has(toId) ? SEA_STEP : Infinity; // land only enterable at a port endpoint
    };
    const majorSet = new Set(major.map(p => p.cellId));

    const seen = new Set<string>();
    for (const p of major) {
      const nearest = [...major]
        .filter(q => q.cellId !== p.cellId)
        .sort((a, b) => angle(cells[p.cellId].center, cells[a.cellId].center) - angle(cells[p.cellId].center, cells[b.cellId].center) || a.cellId - b.cellId)
        .slice(0, SEA_NEIGHBORS);
      for (const q of nearest) {
        const key = `${Math.min(p.cellId, q.cellId)}-${Math.max(p.cellId, q.cellId)}`;
        if (seen.has(key)) continue;
        seen.add(key);
        const path = aStar(cells, p.cellId, q.cellId, seaStep, landHeuristic(cells, q.cellId, stepsPerRadian), MAX_SEA_EXPANSIONS);
        if (path && path.length >= 2) {
          routes.push({ path: path.map(id => cells[id].center), kind: 'searoute', fromCellId: p.cellId, toCellId: q.cellId });
        }
      }
    }
  }
```

> Move the `const majorSet` declaration **above** `seaStep` (hoisting-safe ordering) — it's referenced inside the closure. Fix during implementation so it isn't a temporal-dead-zone bug.

- [ ] **Step 4: Run tests**

Run: `npx vitest run tests/routes.test.ts`
Expected: PASS (roads + sea).
Run: `npm test`
Expected: ALL PASS.

- [ ] **Step 5: typecheck + lint**

Run: `npm run typecheck && npm run lint`
Expected: 0 errors; warnings ≤ 30.

- [ ] **Step 6: Commit**

```bash
git add utils/routes.ts tests/routes.test.ts
git commit -m "Add sea-route pass: dense major-port web"
```

---

### Task 4: `showRoutes` toggle plumbing (App + Controls)

Wire the toggle through state before rendering consumes it. Default **off**.

**Files:**
- Modify: `App.tsx`
- Modify: `components/Controls.tsx`

**Interfaces:**
- Produces: `showRoutes: boolean` / `setShowRivers`-style setter prop-drilled to `WorldViewer` and `Map2D`, and a checkbox in the Controls "Map Overlays" group.

- [ ] **Step 1: Add state in `App.tsx`**

Next to `const [showRivers, setShowRivers] = useState(true);` (line ~122):

```ts
const [showRoutes, setShowRoutes] = useState(false);
```

Pass `showRoutes={showRoutes} setShowRoutes={setShowRoutes}` to `<Controls .../>` (mirror the `showRivers` prop site, line ~648) and `showRoutes={showRoutes}` to `<WorldViewer .../>` and `<Map2D .../>` (lines ~685, ~712).

- [ ] **Step 2: Add the prop + checkbox in `components/Controls.tsx`**

Add to the props interface next to `showRivers: boolean;`:

```ts
showRoutes: boolean;
setShowRoutes: (v: boolean) => void;
```

Destructure `showRoutes, setShowRoutes` alongside `showRivers`. Add a checkbox row mirroring the rivers one (around line 524), using a route-ish lucide icon (e.g. `Route` or `Milestone` from `lucide-react`, importing it with the others):

```tsx
<label className="flex items-center gap-2 cursor-pointer">
  <Route size={12} className={showRoutes ? "text-amber-400" : "text-gray-600"} />
  <span>Roads & Routes</span>
  <input type="checkbox" checked={showRoutes} onChange={e => setShowRoutes(e.target.checked)} className="ml-auto" />
</label>
```

> Match the exact class names / structure of the neighboring `showRivers` label so styling is consistent; the snippet above is the shape, not necessarily the exact markup.

- [ ] **Step 3: typecheck + lint + build**

Run: `npm run typecheck && npm run lint && npm run build`
Expected: 0 errors; warnings ≤ 30; build OK. (No rendering yet — toggle is inert, that's expected.)

- [ ] **Step 4: Commit**

```bash
git add App.tsx components/Controls.tsx
git commit -m "Add showRoutes toggle state + Map Overlays checkbox"
```

---

### Task 5: 3D rendering — `RouteLines`

**Files:**
- Modify: `components/WorldViewer.tsx`

**Interfaces:**
- Consumes: `world.routes`, `showRoutes` prop.
- Produces: a `RouteLines` component mirroring `RiverLines`, rendered inside the scene where `<RiverLines .../>` is (line ~1052); `showRoutes` threaded into the `WorldViewer` prop type + inner component (mirror `showRivers`, lines ~758/773/1079/1163).

- [ ] **Step 1: Add `RouteLines` near `RiverLines`** (`components/WorldViewer.tsx`)

```tsx
const RouteLines: React.FC<{ world: WorldData; visible: boolean }> = ({ world, visible }) => {
  const routes = world.routes;
  const road = useMemo(() => buildRouteGeometry(routes, 'road', visible), [routes, visible]);
  const sea = useMemo(() => buildRouteGeometry(routes, 'searoute', visible), [routes, visible]);
  useEffect(() => () => { road?.dispose(); }, [road]);
  useEffect(() => () => { sea?.dispose(); }, [sea]);
  if (!visible) return null;
  return (
    <>
      {road && (
        <LineSegments geometry={road}>
          <LineBasicMaterial color="#c8a25a" opacity={0.85} transparent linewidth={1.5} />
        </LineSegments>
      )}
      {sea && (
        <LineSegments geometry={sea}>
          <LineDashedMaterial color="#5eb8c8" opacity={0.85} transparent dashSize={0.02} gapSize={0.012} />
        </LineSegments>
      )}
    </>
  );
};
```

with a module-scope helper (mirrors the `RiverLines` batching but lifts routes slightly off the surface and computes line distances for the dashed material):

```tsx
function buildRouteGeometry(
  routes: RouteData[] | undefined,
  kind: 'road' | 'searoute',
  visible: boolean,
): THREE.BufferGeometry | null {
  if (!routes || !visible) return null;
  const positions: number[] = [];
  const LIFT = 1.008; // just above surface; rivers sit at r≈1.0
  for (const r of routes) {
    if (r.kind !== kind || r.path.length < 2) continue;
    const vectors = r.path.map(p => new THREE.Vector3(p.x, p.y, p.z).multiplyScalar(LIFT));
    const curve = new THREE.CatmullRomCurve3(vectors);
    const pts = curve.getPoints(Math.min(60, vectors.length * 4));
    for (let i = 0; i < pts.length - 1; i++) {
      positions.push(pts[i].x, pts[i].y, pts[i].z, pts[i + 1].x, pts[i + 1].y, pts[i + 1].z);
    }
  }
  if (positions.length === 0) return null;
  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
  if (kind === 'searoute') geo.computeLineDistances(); // required for LineDashedMaterial
  return geo;
}
```

Ensure `LineDashedMaterial` and `RouteData` are imported (add to the existing `three` / `types` imports; `LineSegments` and `LineBasicMaterial` are already used by `RiverLines`).

- [ ] **Step 2: Thread `showRoutes` + render**

- Add `showRoutes: boolean` to both the inner scene component props and the outer `WorldViewer` props (mirror every `showRivers` occurrence at lines ~758, ~773, ~1079, ~1163), default `false`.
- Render `<RouteLines world={world} visible={showRoutes} />` immediately after `<RiverLines world={world} visible={showRivers} />` (line ~1052).

- [ ] **Step 3: typecheck + lint + build**

Run: `npm run typecheck && npm run lint && npm run build`
Expected: 0 errors; warnings ≤ 30; build OK.

- [ ] **Step 4: Browser-verify** (manual, per CLAUDE.md)

Start the dev server if not already running (`npm run dev`, port 3000). Generate a world, toggle **Roads & Routes** on. Confirm: solid tan roads connect towns within landmasses; dashed teal lines cross water between coasts; lines sit just above the surface; toggling off removes them cleanly; regenerating keeps them consistent. Note any z-fighting or clutter.

- [ ] **Step 5: Commit**

```bash
git add components/WorldViewer.tsx
git commit -m "Render roads & sea routes on the 3D globe"
```

---

### Task 6: 2D rendering — Map2D (Mercator + Dymaxion)

**Files:**
- Modify: `components/Map2D.tsx`

**Interfaces:**
- Consumes: `world.routes`, `showRoutes` prop.
- Produces: routes drawn in the composite pass of both projections, reusing the river antimeridian-jump handling.

- [ ] **Step 1: Add the `showRoutes` prop**

Add `showRoutes?: boolean;` to the `Map2D` props interface next to `showRivers?`, destructure with `showRoutes = false`, and add `showRoutes` to the redraw dependency arrays wherever `showRivers` appears (line ~725).

- [ ] **Step 2: Draw routes after rivers in the Mercator path**

Immediately after the Mercator river-drawing block (`components/Map2D.tsx:406-`), add a routes block that reuses the exact same equirectangular projection + antimeridian-jump logic, styled per kind:

```ts
if (showRoutes && world.routes) {
  for (const r of world.routes) {
    ctx.beginPath();
    ctx.strokeStyle = r.kind === 'road' ? '#b8863a' : '#3f97a8';
    ctx.lineWidth = r.kind === 'road' ? 1.4 : 1.2;
    ctx.setLineDash(r.kind === 'searoute' ? [5, 4] : []);
    // ...project each path point with the same helper the river loop uses,
    //    breaking the polyline on antimeridian jumps exactly as rivers do...
    ctx.stroke();
  }
  ctx.setLineDash([]);
}
```

Copy the point-projection + jump-break inner loop verbatim from the adjacent river block so behavior is identical.

- [ ] **Step 3: Draw routes after rivers in the Dymaxion path**

Repeat in the Dymaxion composite pass (after the river block at `components/Map2D.tsx:594-`), reusing that path's Dymaxion projector and its seam-break handling, same per-kind styling as Step 2.

- [ ] **Step 4: typecheck + lint + build**

Run: `npm run typecheck && npm run lint && npm run build`
Expected: 0 errors; warnings ≤ 30; build OK.

- [ ] **Step 5: Browser-verify** (manual)

Toggle Roads & Routes on both 2D Mercator and 2D Dymaxion. Confirm roads solid / sea dashed, no stray segments across the antimeridian or Dymaxion net seams, and that lines align with towns.

- [ ] **Step 6: Commit**

```bash
git add components/Map2D.tsx
git commit -m "Render roads & sea routes on 2D Mercator + Dymaxion"
```

---

### Task 7: Export (raster + SVG)

**Files:**
- Modify: `utils/export.ts`

**Interfaces:**
- Consumes: `world.routes`, the live `showRoutes` visibility (WYSIWYG with on-screen toggles).
- Produces: routes drawn in the raster export paths and emitted as an SVG layer group.

- [ ] **Step 1: Thread `showRoutes` into the export signature**

Find how `showRivers` reaches `exportMap` / `exportDymaxionRaster` / the SVG exporter and add a parallel `showRoutes` parameter (default false), passed from the caller's live toggle exactly like rivers/contours.

- [ ] **Step 2: Raster routes**

In each raster export path, after the river-drawing step, draw `world.routes` with the same mirror-corrected projection the rivers use, roads solid / sea dashed (`setLineDash([5,4])`), matching Task 6 colors.

- [ ] **Step 3: SVG routes layer**

In the SVG exporter, emit a `<g id="routes">` group with `<polyline>`/`<path>` elements per route, `stroke` per kind and `stroke-dasharray="5 4"` on sea routes, honoring the same geo-group mirror/counter-mirror convention the existing SVG rivers/borders use.

- [ ] **Step 4: typecheck + lint + build**

Run: `npm run typecheck && npm run lint && npm run build`
Expected: 0 errors; warnings ≤ 30; build OK.

- [ ] **Step 5: Verify export output** (manual)

Export a 4K PNG with routes on — confirm they appear correctly placed. Export SVG and confirm the `routes` group is present and, if convenient, that `xmllint --noout` accepts it.

- [ ] **Step 6: Commit**

```bash
git add utils/export.ts
git commit -m "Include roads & sea routes in PNG + SVG export"
```

---

### Task 8: Full gate pass, browser regression, HANDOFF

**Files:**
- Modify: `HANDOFF.md`

- [ ] **Step 1: Run all four gates**

Run: `npm run typecheck && npm run lint && npm test && npm run build`
Expected: typecheck 0, lint 0 errors / ≤30 warnings, all tests green, build OK. Paste the real output — no success claim without it.

- [ ] **Step 2: Browser regression** (manual)

Generate at least two seeds; verify routes across 3D, 2D Mercator, 2D Dymaxion, toggle on/off, regenerate, and a PNG export — plus a quick check that rivers/labels/borders still render (no regressions from shared-code changes).

- [ ] **Step 3: Update `HANDOFF.md`**

Add a dated Session entry: C3 shipped, the land-component forest decision (and the great-circle-MST bug it fixed), the `pathfinding.ts` extraction, the deferred route→importance feedback loop, and the sea-clutter tuning knob as a known follow-up. Mark C3 done in the pickup notes.

- [ ] **Step 4: Commit**

```bash
git add HANDOFF.md
git commit -m "Record C3 completion in HANDOFF"
```

---

## Self-Review

**Spec coverage:** RouteData/routes (T1) ✓; shared cost extraction (T1) ✓; land-component forest roads + trunks (T2) ✓; dense major-port sea web (T3) ✓; determinism + purity + forest invariant tests (T2/T3) ✓; single toggle (T4) ✓; 3D dashed/solid (T5) ✓; 2D both projections (T6) ✓; raster + SVG export (T7) ✓; no save-schema change (never serialized) ✓; deferred importance feedback (noted, T8) ✓.

> **Gap noted for execution:** the spec's Test 3 "forest invariant (MST edges only, acyclic, edge count == towns − group count)" and Task 3's "sea bridging" negative test are described in the spec but only partially encoded in T2/T3 test code (which assert non-empty + determinism + purity). When executing T2, **add the explicit forest-invariant assertion**: group towns by `(faction, component)`, count `kind:'road'` MST edges excluding trunks, and assert connectivity per group + acyclicity. This was left as a strengthening step rather than fully inlined because it needs the internal grouping exposed for test; expose a small `__roadDebug` export or assert via graph reconstruction from `fromCellId/toCellId`.

**Placeholder scan:** no TBD/TODO; every code step carries real code. Two spots flagged inline for the implementer to resolve mechanically (the `majorSet` hoist in T3; the exact Controls markup in T4) — these are correctness notes, not placeholders.

**Type consistency:** `RouteData { path, kind, fromCellId, toCellId }` used identically across T1 (def), T2/T3 (produce), T5/T6/T7 (consume). `computeRoutes(world, params)` signature stable. `landTerrainStepCost(from, to)` (no params) consistent between pathfinding.ts, worldGen refactor, and routes.ts. `isWaterCell(cell, seaLevel)` consistent. `aStar(cells, start, goal, stepCost, heuristic, maxExpansions)` consistent across pathfinding.ts and both route passes.
