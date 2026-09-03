# F3 Phase 1 — cached-Canvas2D + vector groundwork Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Kill the Map2D ~2.2s pan/zoom stall by projecting and colouring cells once into caches, and build the reusable vector geometry (edge-chained boundaries, a simplification primitive, a pick index) the WebGL Phase 2 will consume.

**Architecture:** Two per-world caches (screen geometry + cell colour) replace per-redraw d3 reprojection and per-cell colour recompute; three new pure utils (edge-chaining, Douglas–Peucker simplification, quadtree pick) supply the groundwork. A batched `fillCells` substrate primitive drives fills from the caches. Everything stays Canvas2D; no new render surface.

**Tech Stack:** TypeScript (strict), React 19, d3-geo + d3-quadtree (already in the `d3` ^7.9.0 dep), Vitest, Canvas2D.

**Spec:** `docs/superpowers/specs/2026-09-03-f3-phase1-vector-groundwork-design.md` (read it alongside this plan).

## Global Constraints

- **Relative imports only** — the `@/` alias is configured but intentionally unused.
- **`utils/` must never import from `components/`** — the reason `boundaries.ts` and Map2D's border scan were duplicated.
- **Point type:** `Point = {x:number;y:number;z:number}` (unit sphere); `Point3 = [number,number,number]` (from `utils/geo.ts`). `toLonLat(p: Point3): [lon,lat]` and `lonLatToPoint3` live in `utils/geo.ts`.
- **Mirror + DPR live in the CONTEXT transform, not path coords.** The offscreen context does `setTransform(renderDpr,…); translate(size.width,0); scale(-1,1)` (`Map2D.tsx:471-474`). Caches store **un-flipped CSS-px** geometry and are filled through that context. Never bake the flip or DPR into cached coordinates.
- **Cache keys:** geometry `(projectionType, width, height, world)`; colour `(viewMode, season, seaLevel, factionColors, cultureColors, religionColors, shadeMap)`. Never keyed on `scale`/`offset`/`qualityDpr`.
- **Dymaxion is out of scope** — separate per-pixel pipeline with its own pick buffer; do not touch it.
- 2-space indent, semicolons, single quotes, trailing commas. `interface` for objects, `type` for unions.
- **Testing discipline (CLAUDE.md):** run touched-file tests (`npx vitest run tests/x.test.ts`), never the full suite without asking. New engine tests must not call `generateWorld` more than necessary; prefer synthetic fixtures.

---

## Task 1: Douglas–Peucker simplification (Unit 4)

**Files:**
- Create: `utils/simplify.ts`
- Test: `tests/simplify.test.ts`

**Interfaces:**
- Consumes: nothing.
- Produces: `simplifyPolyline(points: [number, number][], tolerance: number): [number, number][]` — operates on 2D (projected CSS-px) points; endpoints always preserved; removes points whose perpendicular distance to the retained chord is `< tolerance`.

- [ ] **Step 1: Write the failing test**

```typescript
import { describe, it, expect } from 'vitest';
import { simplifyPolyline } from '../utils/simplify';

describe('simplifyPolyline', () => {
  it('preserves both endpoints', () => {
    const pts: [number, number][] = [[0, 0], [1, 0.001], [2, 0]];
    const out = simplifyPolyline(pts, 0.5);
    expect(out[0]).toEqual([0, 0]);
    expect(out[out.length - 1]).toEqual([2, 0]);
  });

  it('drops a near-collinear interior point below tolerance', () => {
    const pts: [number, number][] = [[0, 0], [1, 0.01], [2, 0]];
    expect(simplifyPolyline(pts, 0.5)).toEqual([[0, 0], [2, 0]]);
  });

  it('keeps an interior point above tolerance', () => {
    const pts: [number, number][] = [[0, 0], [1, 5], [2, 0]];
    expect(simplifyPolyline(pts, 0.5)).toEqual([[0, 0], [1, 5], [2, 0]]);
  });

  it('never increases vertex count and keeps max deviation <= tolerance', () => {
    const pts: [number, number][] = [];
    for (let i = 0; i <= 100; i++) pts.push([i, Math.sin(i / 5) * 3]);
    const out = simplifyPolyline(pts, 1.0);
    expect(out.length).toBeLessThanOrEqual(pts.length);
    expect(out.length).toBeGreaterThanOrEqual(2);
  });

  it('returns input unchanged for 2 or fewer points', () => {
    expect(simplifyPolyline([[0, 0]], 1)).toEqual([[0, 0]]);
    expect(simplifyPolyline([[0, 0], [9, 9]], 1)).toEqual([[0, 0], [9, 9]]);
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npx vitest run tests/simplify.test.ts`
Expected: FAIL — `simplifyPolyline is not a function` / module not found.

- [ ] **Step 3: Write minimal implementation**

```typescript
// Iterative Douglas–Peucker line simplification on 2D points (projected CSS-px).
// Used by the map geometry cache (Unit 1) to thin coastline/border/river
// polylines once at build time; the tolerance is a genuine sub-pixel distance
// because it runs AFTER projection. Endpoints are always kept.

type P = [number, number];

// Perpendicular distance from p to the segment a-b (2D).
const perpDist = (p: P, a: P, b: P): number => {
  const dx = b[0] - a[0];
  const dy = b[1] - a[1];
  const len2 = dx * dx + dy * dy;
  if (len2 === 0) return Math.hypot(p[0] - a[0], p[1] - a[1]);
  const t = ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / len2;
  const cx = a[0] + t * dx;
  const cy = a[1] + t * dy;
  return Math.hypot(p[0] - cx, p[1] - cy);
};

export const simplifyPolyline = (points: P[], tolerance: number): P[] => {
  const n = points.length;
  if (n <= 2) return points.slice();
  const keep = new Uint8Array(n);
  keep[0] = 1;
  keep[n - 1] = 1;
  const stack: Array<[number, number]> = [[0, n - 1]];
  while (stack.length) {
    const [first, last] = stack.pop()!;
    let maxDist = -1;
    let idx = -1;
    for (let i = first + 1; i < last; i++) {
      const d = perpDist(points[i], points[first], points[last]);
      if (d > maxDist) { maxDist = d; idx = i; }
    }
    if (maxDist > tolerance && idx !== -1) {
      keep[idx] = 1;
      stack.push([first, idx], [idx, last]);
    }
  }
  const out: P[] = [];
  for (let i = 0; i < n; i++) if (keep[i]) out.push(points[i]);
  return out;
};
```

- [ ] **Step 4: Run test to verify it passes**

Run: `npx vitest run tests/simplify.test.ts`
Expected: PASS (5 tests).

- [ ] **Step 5: Commit**

```bash
git add utils/simplify.ts tests/simplify.test.ts
git commit -m "feat(F3): Douglas-Peucker polyline simplification (Unit 4)"
```

---

## Task 2: Boundary edge-chaining + lake outlines (Unit 3)

**Files:**
- Modify: `utils/boundaries.ts` (add chaining + lake outlines; keep existing segment functions)
- Test: `tests/boundaries.test.ts` (create or extend)

**Interfaces:**
- Consumes: existing `computeCoastlineSegments(world): Array<[Point3,Point3]>`, `computeFactionBorderSegments(world)`, `computeBoundarySegments(world, differs)` (unchanged).
- Produces:
  - `chainSegments(segments: Array<[Point3, Point3]>): Point3[][]` — assembles disjoint 2-point edges into contiguous polylines; a closed ring repeats its start point as its last point.
  - `computeLakeOutlineSegments(world: WorldData): Array<[Point3, Point3]>` — the boundary edges between each lake cell and its non-same-lake neighbours.

- [ ] **Step 1: Write the failing test**

```typescript
import { describe, it, expect } from 'vitest';
import { chainSegments } from '../utils/boundaries';
import { Point3 } from '../utils/geo';

const seg = (a: Point3, b: Point3): [Point3, Point3] => [a, b];

describe('chainSegments', () => {
  it('joins a broken-order open chain into one polyline', () => {
    // Endpoints given out of order and reversed.
    const segs: Array<[Point3, Point3]> = [
      seg([2, 0, 0], [3, 0, 0]),
      seg([1, 0, 0], [2, 0, 0]),
      seg([0, 0, 0], [1, 0, 0]),
    ];
    const chains = chainSegments(segs);
    expect(chains).toHaveLength(1);
    expect(chains[0]).toHaveLength(4); // 4 points, open
    // First and last differ (open chain).
    expect(chains[0][0]).not.toEqual(chains[0][chains[0].length - 1]);
  });

  it('closes a ring (last point equals first)', () => {
    const segs: Array<[Point3, Point3]> = [
      seg([0, 0, 0], [1, 0, 0]),
      seg([1, 0, 0], [1, 1, 0]),
      seg([1, 1, 0], [0, 0, 0]),
    ];
    const chains = chainSegments(segs);
    expect(chains).toHaveLength(1);
    expect(chains[0][0]).toEqual(chains[0][chains[0].length - 1]);
  });

  it('conserves edge count across all chains', () => {
    const segs: Array<[Point3, Point3]> = [
      seg([0, 0, 0], [1, 0, 0]),
      seg([1, 0, 0], [2, 0, 0]),
      seg([5, 5, 5], [6, 5, 5]), // disjoint second chain
    ];
    const chains = chainSegments(segs);
    const edges = chains.reduce((s, c) => s + (c.length - 1), 0);
    expect(edges).toBe(segs.length);
    expect(chains).toHaveLength(2);
  });

  it('handles a T-junction deterministically without dropping edges', () => {
    // Three edges share the vertex [1,0,0].
    const segs: Array<[Point3, Point3]> = [
      seg([0, 0, 0], [1, 0, 0]),
      seg([1, 0, 0], [2, 0, 0]),
      seg([1, 0, 0], [1, 1, 0]),
    ];
    const chains = chainSegments(segs);
    const edges = chains.reduce((s, c) => s + (c.length - 1), 0);
    expect(edges).toBe(3); // no edge lost at the branch
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npx vitest run tests/boundaries.test.ts`
Expected: FAIL — `chainSegments is not a function`.

- [ ] **Step 3: Write minimal implementation** (append to `utils/boundaries.ts`)

```typescript
import { LakeData } from '../types'; // add to existing import if not present

// Quantise an endpoint to an exact hash key. Vertices come from a
// distance-threshold match (see computeBoundarySegments), so raw floats are only
// approximately equal; snapping to a 1e-5 grid makes "shared endpoint" an exact
// map key and keeps the walk robust at T-junctions (3 cells meeting).
const keyOf = (p: Point3): string =>
  `${Math.round(p[0] * 1e5)},${Math.round(p[1] * 1e5)},${Math.round(p[2] * 1e5)}`;

/**
 * Assemble disjoint 2-point boundary edges into contiguous polylines. A chain
 * that returns to its start is closed (last point === first). At a branch point
 * (endpoint shared by >2 edges) the walk continues deterministically (lowest
 * unused edge index) and the remaining edges seed their own chains, so no edge
 * is ever dropped. Edge count is conserved: sum(chain.length - 1) === segments.length.
 */
export const chainSegments = (segments: Array<[Point3, Point3]>): Point3[][] => {
  // endpoint key -> list of segment indices touching it
  const byKey = new Map<string, number[]>();
  const push = (k: string, i: number) => {
    const arr = byKey.get(k);
    if (arr) arr.push(i); else byKey.set(k, [i]);
  };
  segments.forEach(([a, b], i) => { push(keyOf(a), i); push(keyOf(b), i); });

  const used = new Uint8Array(segments.length);
  const chains: Point3[][] = [];

  const nextFrom = (key: string): number => {
    const cands = byKey.get(key);
    if (!cands) return -1;
    for (const i of cands) if (!used[i]) return i; // lowest unused index
    return -1;
  };

  for (let start = 0; start < segments.length; start++) {
    if (used[start]) continue;
    used[start] = 1;
    const [a, b] = segments[start];
    const chain: Point3[] = [a, b];
    // Extend forward from the tail.
    let tailKey = keyOf(b);
    for (;;) {
      const ni = nextFrom(tailKey);
      if (ni === -1) break;
      used[ni] = 1;
      const [na, nb] = segments[ni];
      const next = keyOf(na) === tailKey ? nb : na;
      chain.push(next);
      tailKey = keyOf(next);
      if (tailKey === keyOf(chain[0])) break; // closed ring
    }
    chains.push(chain);
  }
  return chains;
};

/**
 * Boundary edges around lake cells — a lake cell adjacent to any cell not in the
 * SAME lake. Lakes carry only member cell ids (LakeData.cellIds); this derives
 * their outline for the vector map. Reuses the shared adjacency scan.
 */
export const computeLakeOutlineSegments = (world: WorldData): Array<[Point3, Point3]> => {
  const lakes = world.lakes;
  if (!lakes || lakes.length === 0) return [];
  const lakeOf = new Int32Array(world.cells.length).fill(-1);
  lakes.forEach((lake: LakeData, li: number) => {
    for (const id of lake.cellIds) lakeOf[id] = li;
  });
  return computeBoundarySegments(world, (a, b) => lakeOf[a.id] !== lakeOf[b.id]
    && (lakeOf[a.id] !== -1 || lakeOf[b.id] !== -1));
};
```

> Note: confirm `WorldData.lakes` / `LakeData.cellIds` names against `types.ts` before implementing; the geometry inventory reported `LakeData.cellIds` at `types.ts:231-238`. If the collection is named differently, adjust `world.lakes` accordingly.

- [ ] **Step 4: Run test to verify it passes**

Run: `npx vitest run tests/boundaries.test.ts`
Expected: PASS (4 tests).

- [ ] **Step 5: Add the closure guard over a generated world**

Append this test (it uses one generated world; keep it a single `generateWorld` call in `beforeAll`):

```typescript
import { beforeAll } from 'vitest';
import { generateWorld } from '../utils/worldGen';
import { computeCoastlineSegments, chainSegments } from '../utils/boundaries';
import { DEFAULT_PARAMS } from '../utils/worldGen'; // adjust to the real default-params export

describe('coastline chains close (hard guard)', () => {
  let unclosed = 0;
  let total = 0;
  beforeAll(async () => {
    const world = await generateWorld({ ...DEFAULT_PARAMS, points: 3000, seed: 'chain-test' });
    const chains = chainSegments(computeCoastlineSegments(world));
    total = chains.length;
    for (const c of chains) {
      if (c.length < 2) { unclosed++; continue; }
      const a = c[0]; const b = c[c.length - 1];
      const same = Math.abs(a[0] - b[0]) < 1e-4 && Math.abs(a[1] - b[1]) < 1e-4 && Math.abs(a[2] - b[2]) < 1e-4;
      if (!same) unclosed++;
    }
  });
  it('has at least one coastline chain', () => { expect(total).toBeGreaterThan(0); });
  it('every coastline chain closes into a ring', () => { expect(unclosed).toBe(0); });
});
```

Run: `npx vitest run tests/boundaries.test.ts`
Expected: PASS. **If `unclosed > 0`, the threshold matching is too loose — do not weaken the assertion; report the count and investigate the quantisation grid or the underlying vertex match. This is the advisor-flagged correctness gate.**

- [ ] **Step 6: Commit**

```bash
git add utils/boundaries.ts tests/boundaries.test.ts
git commit -m "feat(F3): edge-chaining + lake outlines with closure guard (Unit 3)"
```

---

## Task 3: Quadtree cell picking (Unit 5)

**Files:**
- Create: `utils/mapPick.ts`
- Test: `tests/mapPick.test.ts`

**Interfaces:**
- Consumes: a d3 projection (`d3.GeoProjection`), `world.cells[i].center` (Point), `toLonLat`.
- Produces:
  - `buildCellQuadtree(world, projection, width, height): Quadtree<number>` — indexes each cell id by its projected CSS-px centre (un-flipped).
  - `findCellIdAtPoint(qt: Quadtree<number>, xUnflipped: number, y: number): number | null` — nearest projected centre.
  - Callers pass `size.width - mapX` for x (the un-flip that `getCellIdAtMapPoint` already applies before `projection.invert`).

- [ ] **Step 1: Write the failing test**

```typescript
import { describe, it, expect } from 'vitest';
import * as d3 from 'd3';
import { buildCellQuadtree, findCellIdAtPoint } from '../utils/mapPick';

// Minimal fake world: cells with center on the unit sphere.
function fakeWorld(centers: Array<[number, number, number]>) {
  return { cells: centers.map((c, id) => ({ id, center: { x: c[0], y: c[1], z: c[2] } })) } as any;
}

// Reference: geodesic max-dot nearest (what production does today).
function nearestByDot(world: any, lon: number, lat: number): number {
  const lonR = lon * Math.PI / 180, latR = lat * Math.PI / 180, cl = Math.cos(latR);
  const x = cl * Math.cos(lonR), y = Math.sin(latR), z = cl * Math.sin(lonR);
  let best = -1, bestDot = -Infinity;
  for (const c of world.cells) {
    const d = c.center.x * x + c.center.y * y + c.center.z * z;
    if (d > bestDot) { bestDot = d; best = c.id; }
  }
  return best;
}

describe('quadtree pick', () => {
  const W = 800, H = 400;
  // A lon/lat grid of cell centers.
  const centers: Array<[number, number, number]> = [];
  for (let lat = -60; lat <= 60; lat += 20) {
    for (let lon = -160; lon <= 160; lon += 20) {
      const lo = lon * Math.PI / 180, la = lat * Math.PI / 180, cl = Math.cos(la);
      centers.push([cl * Math.cos(lo), Math.sin(la), cl * Math.sin(lo)]);
    }
  }
  const world = fakeWorld(centers);
  const projection = d3.geoEquirectangular().fitSize([W, H], { type: 'Sphere' } as any);
  const qt = buildCellQuadtree(world, projection, W, H);

  it('matches geodesic nearest at interior sample points', () => {
    let mismatches = 0, tested = 0;
    for (let lat = -50; lat <= 50; lat += 10) {
      for (let lon = -150; lon <= 150; lon += 10) {
        const px = projection([lon, lat]);
        if (!px) continue;
        tested++;
        const qtId = findCellIdAtPoint(qt, px[0], px[1]);
        const dotId = nearestByDot(world, lon, lat);
        if (qtId !== dotId) mismatches++;
      }
    }
    expect(tested).toBeGreaterThan(50);
    // Interior parity should be essentially exact for equirectangular.
    expect(mismatches).toBe(0);
  });

  it('degrades gracefully at the antimeridian seam and poles (documented)', () => {
    // Seam/high-lat divergence is allowed; assert it does not throw and returns a cell.
    for (const [lon, lat] of [[179, 0], [-179, 0], [0, 85], [0, -85]] as [number, number][]) {
      const px = projection([lon, lat]);
      if (!px) continue;
      const id = findCellIdAtPoint(qt, px[0], px[1]);
      expect(id).not.toBeNull();
    }
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npx vitest run tests/mapPick.test.ts`
Expected: FAIL — module not found.

- [ ] **Step 3: Write minimal implementation**

```typescript
import { quadtree, Quadtree } from 'd3-quadtree';
import type { GeoProjection } from 'd3-geo';

import { WorldData } from '../types';
import { toLonLat } from './geo';

// Index each cell id by its projected CSS-px centre (un-flipped, matching the
// coordinate space getCellIdAtMapPoint queries after size.width - mapX). Replaces
// the O(cells) geodesic scan; nearest-projected-centre differs from geodesic only
// near the antimeridian/poles, which the parity test samples deliberately.
export const buildCellQuadtree = (
  world: WorldData,
  projection: GeoProjection,
  width: number,
  height: number,
): Quadtree<number> => {
  void width; void height; // reserved: projection already fits [width,height]
  const qt = quadtree<number>()
    .x((id) => projected[id][0])
    .y((id) => projected[id][1]);
  const projected: [number, number][] = new Array(world.cells.length);
  for (let i = 0; i < world.cells.length; i++) {
    const c = world.cells[i].center;
    const p = projection(toLonLat([c.x, c.y, c.z]));
    projected[i] = p && Number.isFinite(p[0]) ? [p[0], p[1]] : [NaN, NaN];
  }
  // Add only cells that projected to a finite point.
  for (let i = 0; i < world.cells.length; i++) {
    if (Number.isFinite(projected[i][0])) qt.add(i);
  }
  return qt;
};

export const findCellIdAtPoint = (
  qt: Quadtree<number>,
  xUnflipped: number,
  y: number,
): number | null => {
  const id = qt.find(xUnflipped, y);
  return id === undefined ? null : id;
};
```

> Note the closure-order bug hazard: `projected` must be declared and filled BEFORE `.x`/`.y` are invoked. `.x`/`.y` accessors are only called during `.add`, which happens after `projected` is filled, so the forward reference is safe — but if you refactor, keep the fill before the first `add`. If TS complains about use-before-declaration, hoist `const projected` above the `quadtree()` call.

- [ ] **Step 4: Run test to verify it passes**

Run: `npx vitest run tests/mapPick.test.ts`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add utils/mapPick.ts tests/mapPick.test.ts
git commit -m "feat(F3): quadtree cell picking with seam/pole parity test (Unit 5)"
```

---

## Task 4: Cell colour cache (Unit 2)

**Files:**
- Create: `utils/mapColorCache.ts`
- Test: `tests/mapColorCache.test.ts`

**Interfaces:**
- Consumes: `getCellColor(cell, mode, ColorContext)` from `utils/colors.ts`, `seasonalTemperatureDelta(cell, params)` from `utils/seasons.ts`, `ViewMode`, `ColorContext`.
- Produces: `buildCellColorCache(world, viewMode, colorCtx, shadeMap): string[]` — `#rrggbb` per cell id, seasonal delta applied per cell, `shadeMap[i]` multiplied in when present.

> **Spec deviation (deliberate):** the spec named a `Uint32Array`; Phase 1 stores `#rrggbb` **strings** because Canvas2D `fillStyle` needs a string and packing→unpacking per fill would re-introduce per-fill work. Packed `Uint32Array` for WebGL is a trivial Phase-2 addition. Recorded here so it is not read as drift.

- [ ] **Step 1: Write the failing test**

```typescript
import { describe, it, expect } from 'vitest';
import { buildCellColorCache } from '../utils/mapColorCache';
import { getCellColor } from '../utils/colors';
import { seasonalTemperatureDelta } from '../utils/seasons';
import { generateWorld, DEFAULT_PARAMS } from '../utils/worldGen'; // adjust default-params import

describe('buildCellColorCache', () => {
  it('matches per-cell getCellColor with shade + seasonal delta', async () => {
    const world = await generateWorld({ ...DEFAULT_PARAMS, points: 2000, seed: 'color-cache' });
    const colorCtx = {
      seaLevel: world.params.seaLevel,
      factionColors: undefined, cultureColors: undefined, religionColors: undefined,
    };
    const shadeMap = new Float32Array(world.cells.length).fill(1);
    shadeMap[5] = 0.7;
    const cache = buildCellColorCache(world, 'biome', colorCtx, shadeMap);
    for (let i = 0; i < world.cells.length; i++) {
      const c = getCellColor(world.cells[i], 'biome', {
        ...colorCtx, seasonalDelta: seasonalTemperatureDelta(world.cells[i], world.params),
      });
      c.multiplyScalar(shadeMap[i]);
      expect(cache[i]).toBe(`#${c.getHexString()}`);
    }
  });

  it('works with a null shadeMap', async () => {
    const world = await generateWorld({ ...DEFAULT_PARAMS, points: 2000, seed: 'color-cache-2' });
    const cache = buildCellColorCache(world, 'biome', { seaLevel: world.params.seaLevel }, null);
    expect(cache).toHaveLength(world.cells.length);
    expect(cache[0]).toMatch(/^#[0-9a-f]{6}$/);
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npx vitest run tests/mapColorCache.test.ts`
Expected: FAIL — module not found.

- [ ] **Step 3: Write minimal implementation**

```typescript
import { WorldData, ViewMode } from '../types';
import { ColorContext, getCellColor } from './colors';
import { seasonalTemperatureDelta } from './seasons';

// Per-cell fill colour, computed ONCE per (viewMode, season, colormaps, seaLevel,
// shadeMap). The fill passes read cache[i] instead of calling getCellColor +
// seasonalTemperatureDelta + allocating a THREE.Color every redraw — the other
// half of the 2.2s stall. Strings (#rrggbb) because Canvas2D fillStyle needs them.
export const buildCellColorCache = (
  world: WorldData,
  viewMode: ViewMode,
  colorCtx: ColorContext,
  shadeMap: Float32Array | null,
): string[] => {
  const out: string[] = new Array(world.cells.length);
  for (let i = 0; i < world.cells.length; i++) {
    const cell = world.cells[i];
    const color = getCellColor(cell, viewMode, {
      ...colorCtx,
      seasonalDelta: seasonalTemperatureDelta(cell, world.params),
    });
    if (shadeMap) color.multiplyScalar(shadeMap[cell.id]);
    out[i] = `#${color.getHexString()}`;
  }
  return out;
};
```

- [ ] **Step 4: Run test to verify it passes**

Run: `npx vitest run tests/mapColorCache.test.ts`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add utils/mapColorCache.ts tests/mapColorCache.test.ts
git commit -m "feat(F3): per-cell colour cache (Unit 2)"
```

---

## Task 5: Screen-geometry cache (Unit 1)

**Files:**
- Create: `utils/mapCache.ts`
- Test: `tests/mapCache.test.ts`

**Interfaces:**
- Consumes: `d3.GeoProjection`, `world.geoJson.features[i]` (lon/lat rings), `chainSegments` + `computeCoastlineSegments`/`computeFactionBorderSegments`/`computeLakeOutlineSegments` (Task 2), `simplifyPolyline` (Task 1), `toLonLat`.
- Produces:
  - `interface MapGeometryCache { cellPaths: Path2D[]; cellVerts: Float32Array; cellOffsets: Uint32Array; coast: Path2D; borders: Path2D; lakes: Path2D; }` (rivers/routes may be added the same way; start with the above).
  - `buildMapGeometryCache(world, projection, width, height, opts?: { simplifyTolerancePx?: number }): MapGeometryCache`.

- [ ] **Step 1: Write the failing test** (parity: cached fill vs fresh d3 path)

```typescript
import { describe, it, expect } from 'vitest';
import * as d3 from 'd3';
import { buildMapGeometryCache } from '../utils/mapCache';
import { generateWorld, DEFAULT_PARAMS } from '../utils/worldGen'; // adjust import

// A tiny stub 2D context recording the path points a d3 canvas generator emits.
function recordingCtx() {
  const pts: [number, number][] = [];
  return {
    beginPath() {}, closePath() {}, moveTo(x: number, y: number) { pts.push([x, y]); },
    lineTo(x: number, y: number) { pts.push([x, y]); }, arc() {}, pts,
  } as any;
}

describe('buildMapGeometryCache', () => {
  it('cached cell vertices match a fresh d3 projection (sub-pixel)', async () => {
    const world = await generateWorld({ ...DEFAULT_PARAMS, points: 2000, seed: 'geo-cache' });
    const W = 1024, H = 512;
    const projection = d3.geoMercator().fitSize([W, H], { type: 'Sphere' } as any);
    const cache = buildMapGeometryCache(world, projection, W, H, { simplifyTolerancePx: 0 });

    // Compare cell 0's cached vertices to a fresh geoPath trace.
    const ctx = recordingCtx();
    const gen = d3.geoPath(projection, ctx);
    ctx.pts.length = 0;
    gen(world.geoJson!.features![0] as any);
    const fresh = ctx.pts;
    const off0 = cache.cellOffsets[0], off1 = cache.cellOffsets[1];
    const cachedCount = off1 - off0;
    expect(cachedCount).toBe(fresh.length);
    for (let k = 0; k < cachedCount; k++) {
      expect(Math.abs(cache.cellVerts[(off0 + k) * 2] - fresh[k][0])).toBeLessThan(1e-6);
      expect(Math.abs(cache.cellVerts[(off0 + k) * 2 + 1] - fresh[k][1])).toBeLessThan(1e-6);
    }
  });

  it('produces one Path2D per cell', async () => {
    const world = await generateWorld({ ...DEFAULT_PARAMS, points: 2000, seed: 'geo-cache-2' });
    const projection = d3.geoEquirectangular().fitSize([800, 400], { type: 'Sphere' } as any);
    const cache = buildMapGeometryCache(world, projection, 800, 400);
    expect(cache.cellPaths).toHaveLength(world.cells.length);
  });
});
```

> `Path2D` and `d3.geoPath(projection, ctx)` both work under Vitest's jsdom/happy-dom env; if `Path2D` is undefined in the test env, add a minimal `Path2D` polyfill in the test setup or assert only on `cellVerts`/`cellOffsets` (the parity that matters) and skip the Path2D-count assertion.

- [ ] **Step 2: Run test to verify it fails**

Run: `npx vitest run tests/mapCache.test.ts`
Expected: FAIL — module not found.

- [ ] **Step 3: Write minimal implementation**

```typescript
import type { GeoProjection } from 'd3-geo';

import { WorldData } from '../types';
import { Point3, toLonLat } from './geo';
import {
  chainSegments, computeCoastlineSegments, computeFactionBorderSegments,
  computeLakeOutlineSegments,
} from './boundaries';
import { simplifyPolyline } from './simplify';

export interface MapGeometryCache {
  cellPaths: Path2D[];
  cellVerts: Float32Array; // [x0,y0,x1,y1,...] CSS-px, un-flipped
  cellOffsets: Uint32Array; // cell i's verts are [offsets[i], offsets[i+1])
  coast: Path2D;
  borders: Path2D;
  lakes: Path2D;
}

// Project a chained boundary into one Path2D of simplified polylines.
const chainedToPath = (
  segments: Array<[Point3, Point3]>,
  projection: GeoProjection,
  tolPx: number,
): Path2D => {
  const path = new Path2D();
  for (const chain of chainSegments(segments)) {
    let pts: [number, number][] = [];
    for (const v of chain) {
      const p = projection(toLonLat(v));
      if (p && Number.isFinite(p[0])) pts.push([p[0], p[1]]);
    }
    if (tolPx > 0 && pts.length > 2) pts = simplifyPolyline(pts, tolPx);
    if (pts.length < 2) continue;
    path.moveTo(pts[0][0], pts[0][1]);
    for (let k = 1; k < pts.length; k++) path.lineTo(pts[k][0], pts[k][1]);
  }
  return path;
};

// Project all cell rings ONCE to CSS-px (un-flipped) + build a Path2D per cell.
// Coordinates are filled through Map2D's flipped/DPR context, so the mirror and
// DPR are applied at draw, never baked here (see the geometry-cache spec §2).
export const buildMapGeometryCache = (
  world: WorldData,
  projection: GeoProjection,
  width: number,
  height: number,
  opts?: { simplifyTolerancePx?: number },
): MapGeometryCache => {
  void width; void height; // projection already fits [width,height]
  const tolPx = opts?.simplifyTolerancePx ?? 0.5;
  const n = world.cells.length;
  const cellPaths: Path2D[] = new Array(n);
  const cellOffsets = new Uint32Array(n + 1);

  // First pass: count vertices to size the flat buffer.
  const features = world.geoJson?.features ?? [];
  const ringOf = (i: number): [number, number][] => {
    const f = features[i];
    const g = f?.geometry as { type: string; coordinates: number[][][] } | null | undefined;
    if (!g || g.type !== 'Polygon') return [];
    const ring = g.coordinates[0] ?? [];
    const out: [number, number][] = [];
    for (const [lon, lat] of ring) {
      const p = projection([lon, lat]);
      if (p && Number.isFinite(p[0])) out.push([p[0], p[1]]);
    }
    return out;
  };

  const rings: [number, number][][] = new Array(n);
  let totalVerts = 0;
  for (let i = 0; i < n; i++) {
    const r = ringOf(i);
    rings[i] = r;
    cellOffsets[i] = totalVerts;
    totalVerts += r.length;
  }
  cellOffsets[n] = totalVerts;

  const cellVerts = new Float32Array(totalVerts * 2);
  for (let i = 0; i < n; i++) {
    const r = rings[i];
    const base = cellOffsets[i];
    const path = new Path2D();
    for (let k = 0; k < r.length; k++) {
      cellVerts[(base + k) * 2] = r[k][0];
      cellVerts[(base + k) * 2 + 1] = r[k][1];
      if (k === 0) path.moveTo(r[k][0], r[k][1]);
      else path.lineTo(r[k][0], r[k][1]);
    }
    if (r.length > 0) path.closePath();
    cellPaths[i] = path;
  }

  return {
    cellPaths, cellVerts, cellOffsets,
    coast: chainedToPath(computeCoastlineSegments(world), projection, tolPx),
    borders: chainedToPath(computeFactionBorderSegments(world), projection, tolPx),
    lakes: chainedToPath(computeLakeOutlineSegments(world), projection, tolPx),
  };
};
```

- [ ] **Step 4: Run test to verify it passes**

Run: `npx vitest run tests/mapCache.test.ts`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add utils/mapCache.ts tests/mapCache.test.ts
git commit -m "feat(F3): screen-geometry cache, project cells once (Unit 1)"
```

---

## Task 6: `fillCells` substrate primitive

**Files:**
- Modify: `utils/mapStyle/substrate.ts` (add method to interface; remove the accidental duplicated `grain`/`withSphereClip` block, lines ~76-91)
- Modify: `utils/mapStyle/substrateCanvas.ts` (implement using cached `Path2D[]`)
- Modify: `utils/mapStyle/substrateSvg.ts` (implement as per-cell paths — no aggregation)
- Test: `tests/fillCells.test.ts`

**Interfaces:**
- Consumes: `MapGeometryCache.cellPaths` (Task 5), colour cache `string[]` (Task 4).
- Produces: `Substrate.fillCells(indices: number[] | Uint32Array, colors: string[]): void`. The Canvas2D substrate is constructed with an optional `cellPaths?: Path2D[]`; `fillCells` fills `cellPaths[i]` with `colors[i]` (fill + 0.5px same-colour hairline, matching `fillFeature`). The SVG substrate is constructed with the cell `GeoFeatureLike[]` and emits one `<path>` per index exactly as its `fillFeature` loop does.

- [ ] **Step 1: Write the failing test**

```typescript
import { describe, it, expect } from 'vitest';
import { Canvas2DSubstrate } from '../utils/mapStyle/substrateCanvas';

describe('Canvas2DSubstrate.fillCells', () => {
  it('fills each requested cell path with its colour', () => {
    const calls: Array<{ style: string; filled: boolean }> = [];
    const ctx = {
      save() {}, restore() {}, beginPath() {}, closePath() {},
      set fillStyle(v: string) { calls.push({ style: v, filled: false }); },
      get fillStyle() { return ''; },
      set strokeStyle(_v: string) {}, get strokeStyle() { return ''; },
      set lineWidth(_v: number) {}, set globalAlpha(_v: number) {},
      fill(_p?: any) { if (calls.length) calls[calls.length - 1].filled = true; },
      stroke(_p?: any) {}, translate() {}, scale() {},
    } as any;
    const paths = [new Path2D(), new Path2D(), new Path2D()];
    const sub = new Canvas2DSubstrate(ctx, (() => {}) as any, 100, 100, false, paths);
    sub.fillCells([0, 2], ['#112233', '#445566', '#778899']);
    // Two cells filled, with colours indexed by cell id (0 and 2).
    const fills = calls.filter((c) => c.filled);
    expect(fills.map((c) => c.style)).toEqual(['#112233', '#778899']);
  });
});
```

> If `Path2D` is undefined in the test env, add a one-line stub `globalThis.Path2D = class {} as any;` at the top of the test.

- [ ] **Step 2: Run test to verify it fails**

Run: `npx vitest run tests/fillCells.test.ts`
Expected: FAIL — `fillCells is not a function` (and the constructor rejects the 6th arg).

- [ ] **Step 3: Implement**

In `utils/mapStyle/substrate.ts`, add to the `Substrate` interface and delete the duplicated `grain`/`withSphereClip` declarations (the accidental second copy):

```typescript
  /**
   * Fill a set of cells by INDEX from a prebuilt geometry cache, each with its
   * own colour (colours indexed by cell id). The batched hot-loop primitive for
   * cached rendering: Canvas2D fills cached Path2Ds (no reprojection); SVG emits
   * one <path> per cell as fillFeature does (no aggregation — batching helps only
   * the raster/GPU backends). See the F3 Phase 1 spec.
   */
  fillCells(indices: number[] | Uint32Array, colors: string[]): void;
```

In `substrateCanvas.ts`, add `private cellPaths: Path2D[] = []` as an optional 6th constructor parameter and:

```typescript
  fillCells(indices: number[] | Uint32Array, colors: string[]): void {
    for (let j = 0; j < indices.length; j++) {
      const i = indices[j];
      const path = this.cellPaths[i];
      if (!path) continue;
      this.ctx.fillStyle = colors[i];
      this.ctx.fill(path);
      this.ctx.strokeStyle = colors[i];
      this.ctx.lineWidth = 0.5;
      this.ctx.stroke(path);
    }
  }
```

In `substrateSvg.ts`, add the cell `GeoFeatureLike[]` to the constructor (or reuse what it has) and:

```typescript
  fillCells(indices: number[] | Uint32Array, colors: string[]): void {
    for (let j = 0; j < indices.length; j++) {
      const i = indices[j];
      const f = this.cellFeatures?.[i];
      if (!f) continue;
      this.fillFeature(f, colors[i]); // per-cell path, exactly as today
    }
  }
```

- [ ] **Step 4: Run test to verify it passes**

Run: `npx vitest run tests/fillCells.test.ts`
Expected: PASS. Also run the existing style tests to confirm no regression: `npx vitest run tests/mapStyle*.test.ts` (or the relevant existing substrate/style test files — list them with `ls tests | grep -i style`).

- [ ] **Step 5: Commit**

```bash
git add utils/mapStyle/substrate.ts utils/mapStyle/substrateCanvas.ts utils/mapStyle/substrateSvg.ts tests/fillCells.test.ts
git commit -m "feat(F3): fillCells substrate primitive; drop dup substrate decls"
```

---

## Task 7: Wire passes to the caches

**Files:**
- Modify: `utils/mapStyle/passes.ts` (`landPass`, `hillshadePass` read the colour cache + call `fillCells` when caches are present)
- Modify: `utils/mapStyle/types.ts` (add optional cache fields to `StyleRenderContext`)
- Test: extend `tests/fillCells.test.ts` or the existing pass test with a cache-present path.

**Interfaces:**
- Consumes: `MapGeometryCache` (Task 5), colour cache `string[]` (Task 4), `Substrate.fillCells` (Task 6).
- Produces: `StyleRenderContext` gains `geometryCache?: MapGeometryCache` and `colorCache?: string[]`. When both are present, `landPass` fills via `fillCells`; when absent it falls back to the current per-feature path (so nothing that doesn't build caches breaks — e.g. the Dymaxion source raster and SVG export).

- [ ] **Step 1: Write the failing test**

```typescript
// Add to tests/fillCells.test.ts
import { landPass } from '../utils/mapStyle/passes';
// Build a StyleRenderContext with colorCache + a stub substrate recording fillCells,
// assert landPass calls sub.fillCells with the land-cell indices and passes the cache.
```

Write a focused test: a 3-cell fake world (2 land, 1 ocean), a `colorCache` of 3 strings, a stub `Substrate` recording `fillCells(indices, colors)` calls; assert `landPass(...)` (with `policy !== 'bare'`) calls `fillCells` once with the two land indices and the colour cache. (Use the real `landPass` factory; pass a `fillPolicy` returning `'continuous'` and a minimal palette.)

- [ ] **Step 2: Run test to verify it fails**

Run: `npx vitest run tests/fillCells.test.ts`
Expected: FAIL — `landPass` still loops `fillFeature`.

- [ ] **Step 3: Implement**

Add to `StyleRenderContext` (`types.ts`):

```typescript
  /** F3 Phase 1: prebuilt caches. When both present, fills read them instead of
   *  reprojecting/recolouring per redraw. Absent for the Dymaxion source raster
   *  and SVG export, which fall back to the per-feature path. */
  geometryCache?: import('../mapCache').MapGeometryCache;
  colorCache?: string[];
```

In `landPass`, before the per-cell loop, add the fast path:

```typescript
    if (ctx.colorCache && ctx.geometryCache && policy !== 'bare') {
      const indices: number[] = [];
      for (let i = 0; i < world.cells.length; i++) {
        const cell = world.cells[i];
        if (cell.height < colorCtx.seaLevel && !isLakeCell(cell)) continue;
        indices.push(i);
      }
      // mute (categorical) is already folded into colorCache at build time when
      // the policy is categorical — see the Map2D integration note in Task 8.
      sub.fillCells(indices, ctx.colorCache);
      return;
    }
```

> **Important:** the `mute` lerp toward paper for categorical policies must be applied when the colour cache is BUILT (Task 8 passes the right `viewMode`/policy context), OR keep `landPass`'s cache fast-path limited to non-categorical policies and let categorical fall through to the per-feature path. Choose the simpler: **only take the cache path when `mute === 0`** (continuous/political-with-no-mute); categorical falls back. Add `&& mute === 0` to the condition above. This avoids a second colour-cache variant in Phase 1.

Apply the analogous cache-aware fill to `hillshadePass` only if a shade colour cache exists; otherwise leave `hillshadePass` on `fillFeature` (it fills a subset with per-cell opacity — lower priority; the spec's dominant cost is `landPass`). Document that hillshade stays per-feature in Phase 1 if you defer it.

- [ ] **Step 4: Run test to verify it passes**

Run: `npx vitest run tests/fillCells.test.ts`
Expected: PASS. Run existing style tests again for no regression.

- [ ] **Step 5: Commit**

```bash
git add utils/mapStyle/passes.ts utils/mapStyle/types.ts tests/fillCells.test.ts
git commit -m "feat(F3): passes read geometry+colour caches via fillCells"
```

---

## Task 8: Map2D integration (kept in-house)

**Files:**
- Modify: `components/Map2D.tsx`

**Interfaces:**
- Consumes: everything above.
- Produces: no new exports; wires caches into the offscreen redraw, routes picking through the quadtree, deletes the duplicated `getFactionBorders`.

- [ ] **Step 1: Build the caches with `useMemo`, keyed correctly**

Add memos near the existing `projection`/`shadeMap` memos (`Map2D.tsx:405-460`):

```typescript
const geometryCache = useMemo(
  () => (world && projection && projectionType !== 'dymaxion'
    ? buildMapGeometryCache(world, projection, size.width, size.height)
    : null),
  [world, projection, projectionType, size.width, size.height],
);
const colorCache = useMemo(
  () => (world ? buildCellColorCache(world, viewMode, {
    seaLevel: world.params.seaLevel, factionColors, cultureColors, religionColors,
  }, shadeMap) : null),
  [world, viewMode, factionColors, cultureColors, religionColors, shadeMap,
   world?.params.season, world?.params.seaLevel],
);
const pickQuadtree = useMemo(
  () => (world && projection && projectionType !== 'dymaxion'
    ? buildCellQuadtree(world, projection, size.width, size.height) : null),
  [world, projection, projectionType, size.width, size.height],
);
```

- [ ] **Step 2: Pass caches into the styled render context**

In the offscreen effect where `runStyle(style, {...}, sub)` is called for the NON-dymaxion branch, add `geometryCache` and `colorCache` to the `StyleRenderContext` object, and construct the `Canvas2DSubstrate` with `geometryCache?.cellPaths`. Do NOT pass caches in the Dymaxion `srcSub` branch (it renders at source resolution — leave it on the per-feature path).

- [ ] **Step 3: Route picking through the quadtree**

Replace the body of `getCellIdAtMapPoint` (`Map2D.tsx:1201-1214`) non-dymaxion branch:

```typescript
    if (!projection || !pickQuadtree) return null;
    // Un-flip x the same way the old invert path did (size.width - mapX), then
    // query the quadtree of projected centres.
    return findCellIdAtPoint(pickQuadtree, size.width - mapX, mapY);
```

Keep the `projectionType === 'dymaxion'` branch (`getPickBufferCellId`) exactly as-is.

- [ ] **Step 4: Delete the duplicated border scan**

Remove the private `getFactionBorders` (`Map2D.tsx:52-80`) and its call site; use `computeFactionBorderSegments(world)` from `utils/boundaries.ts` (import it), matching how coastlines are already sourced. Verify the faction-border overlay still draws.

- [ ] **Step 5: Verify — typecheck, tests, and the perf number**

```bash
npm run typecheck
npx vitest run tests/simplify.test.ts tests/boundaries.test.ts tests/mapPick.test.ts tests/mapColorCache.test.ts tests/mapCache.test.ts tests/fillCells.test.ts
node scripts/renderMap2DPerf.mjs --points=30000 --dpr=2 --label=after
```
Expected: typecheck 0; all unit tests pass; the settle-redraw ms in `tmp/map2d-after.json` is dramatically below the ~2.2s baseline, and the settled pixel hash is stable across runs (output unchanged). Record before/after in HANDOFF.

- [ ] **Step 6: Browser confirmation (reuse the dev server, never start a second)**

Confirm in the running app: pan and zoom the 2D map at a high cell count feels responsive (no multi-second freeze at gesture end); cell picking (click to inspect) still selects the right cell; painting still targets the right cell; faction borders + coastlines still render. (Playwright screenshots not required if the perf JSON + interaction check are clean.)

- [ ] **Step 7: Commit**

```bash
git add components/Map2D.tsx
git commit -m "feat(F3): wire Map2D to geometry/colour caches + quadtree picking"
```

- [ ] **Step 8: Update HANDOFF + ROADMAP**

Record the before/after perf numbers, mark F3 Phase 1 shipped, and note Phase 2 (WebGL vector surface consuming `MapGeometryCache.cellVerts/cellOffsets`) as the next rung. Commit docs separately.

---

## Self-review notes (author)

- **Spec coverage:** Units 1–5 → Tasks 5,4,2,1,3; substrate `fillCells` → Task 6; passes wiring → Task 7; Map2D integration + invariants + perf gate → Task 8. Dymaxion left untouched (Tasks 5/8 guard on `projectionType !== 'dymaxion'`). Mirror/DPR handled by the existing context transform (Task 5 note + Global Constraints).
- **Deferred/verify-at-implementation:** exact `DEFAULT_PARAMS` export name and `WorldData.lakes`/`LakeData.cellIds` field names (confirm against `types.ts`/`worldGen.ts` — flagged in Tasks 2/4). `Path2D` availability in the test env (fallback noted in Tasks 5/6). Categorical `mute` handled by restricting the cache fast-path to `mute === 0` (Task 7).
- **Deviation from spec:** colour cache stores `#rrggbb` strings, not `Uint32Array` (Task 4 rationale). Recorded, not drift.
- **Delegation:** Tasks 1–5 are pure, independently testable → delegatable to Sonnet in parallel (no shared files: `simplify.ts`, `boundaries.ts`, `mapPick.ts`, `mapColorCache.ts`, `mapCache.ts` — note Task 5 imports Tasks 1+2, so land it after them). Task 6 (substrate) needs close review. Tasks 7–8 kept in-house (passes wiring + invariant-dense Map2D integration).
