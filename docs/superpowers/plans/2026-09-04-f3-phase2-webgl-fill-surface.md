# F3 Phase 2 — WebGL Fill Surface Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the Default-style, non-Dymaxion 2D map a true vector surface — crisp cell fills, coastlines, and labels at any zoom — by replacing the offscreen-blit magnification with a Three.js GPU fill layer plus a per-frame Canvas2D vector overlay.

**Architecture:** Two synchronized layers driven by the same `(scale, offset, dpr, mirror)`: a Three.js `OrthographicCamera` surface drawing earcut-triangulated cells with per-vertex colour in one draw call, under a transparent Canvas2D overlay that redraws linework and labels each frame through the same composed transform. WebGL init failure or context loss falls back to the existing blit renderer; Parchment, Dymaxion, and the 3D globe keep the current path.

**Tech Stack:** React 19, Three.js 0.182.0, `earcut` (new dep), d3-geo, TypeScript strict, Vitest.

**Spec:** `docs/superpowers/specs/2026-09-04-f3-phase2-webgl-fill-surface-design.md`

## Global Constraints

- **Relative imports only** — the `@/` alias is configured but intentionally unused.
- **Lint ratchet is 30/30 warnings, zero errors, zero headroom.** The next `any` anywhere breaks CI. Type all new Three/GL code from the first commit; no `as any`.
- **`getCellColor` takes a `ColorContext` object**, with required `ctx.seaLevel = world.params.seaLevel` (never hardcoded); categorical modes need their colour maps present.
- **Reuse geometry across paint strokes** — `world.cells` identity is the structural key; never replace the `cells` array outside full regeneration. Paint mutates colours only.
- **Every GPU/`useMemo` resource has a matching disposal effect** (WorldViewer pattern).
- **The mirror trap** (`docs/map-styles.md` §"The mirror trap"): 2D paths flip horizontally; polygons/strokes land correct through the flipped context, text/glyphs are compensated inside their own draw functions. Reproduce the existing composed transform exactly — do not add or remove a flip.
- **Test discipline** (`CLAUDE.md`): default to touched-file runs (`npx vitest run tests/foo.test.ts`); full suite only before merge, with Matt's speed choice; never alongside browser automation; never two vitest invocations at once. Do not add a test that generates more than one world per case.
- **Determinism / no pipeline change** — Phase 2 touches only rendering. `generateWorld` output must stay byte-identical.
- **Commit in small scoped chunks.** End commit messages with the Co-Authored-By and Claude-Session trailers used on the spec commit.

---

### Task 1: Subpath-boundary index in `mapCache`

`cellVerts` concatenates every sub-ring's points with no boundary marker, so it cannot be tessellated as-is (its own header comment says a consumer must add its own index first). Emit per-cell sub-ring offsets from the `perCell` data already built during capture.

**Files:**
- Modify: `utils/mapCache.ts` (the `MapGeometryCache` interface + `buildMapGeometryCache`)
- Test: `tests/mapCache.test.ts` (extend the existing first `it` — reuse its generated world; no new `generateWorld`)

**Interfaces:**
- Produces: `MapGeometryCache.cellSubOffsets: Uint32Array` — flat concatenation, per cell, of each sub-ring's **start index into `cellVerts` (vertex units, not float units)**, terminated by that cell's `cellOffsets[i+1]`. Indexed by a companion `cellSubStart: Uint32Array` of length `n+1` where cell `i`'s sub-ring starts occupy `cellSubOffsets[cellSubStart[i] .. cellSubStart[i+1])`. A cell with zero sub-rings occupies zero entries. Consumed by Task 2 (`tessellate`).

- [ ] **Step 1: Write the failing test** — add to the existing parity `it` block in `tests/mapCache.test.ts`, inside the per-cell loop, after the `cachedCount` assertions. It reuses `ctx.subpaths` (already recorded per cell):

```ts
// --- Task 1: sub-ring boundary index parity ---
// cellSubOffsets/cellSubStart must reproduce the SAME subpath split d3 emits.
const subStart = cache.cellSubStart[i];
const subEnd = cache.cellSubStart[i + 1];
const subCount = subEnd - subStart;
expect(subCount).toBe(ctx.subpaths.length); // one entry per emitted sub-ring
let acc = off0;
for (let s = 0; s < subCount; s++) {
  // each recorded sub-ring's start vertex index, in cellVerts vertex units
  expect(cache.cellSubOffsets[subStart + s]).toBe(acc);
  acc += ctx.subpaths[s].length;
}
expect(acc).toBe(off1); // sub-rings tile the whole cell range exactly
```

- [ ] **Step 2: Run it to verify it fails**

Run: `npx vitest run tests/mapCache.test.ts`
Expected: FAIL — `cache.cellSubStart` / `cache.cellSubOffsets` are `undefined`.

- [ ] **Step 3: Implement in `utils/mapCache.ts`** — add the two fields to the interface and populate them in the existing per-cell build loop (the one that already walks `perCell[i]` sub-paths building `cellPaths`). Add near the `cellVerts`/`cellOffsets` declarations:

```ts
// In MapGeometryCache:
cellSubStart: Uint32Array;   // length n+1: cell i's sub-ring entries are [cellSubStart[i], cellSubStart[i+1])
cellSubOffsets: Uint32Array; // per sub-ring: start vertex index into cellVerts (vertex units)
```

Build them alongside the existing loops. First pass (the one computing `cellOffsets`) also counts sub-rings:

```ts
const cellSubStart = new Uint32Array(n + 1);
let totalSubs = 0;
for (let i = 0; i < n; i++) {
  cellSubStart[i] = totalSubs;
  totalSubs += perCell[i].length;
}
cellSubStart[n] = totalSubs;
const cellSubOffsets = new Uint32Array(totalSubs);
```

In the second loop (the one filling `cellVerts` with `cursor`), record each sub-ring's start before writing its points:

```ts
let subCursor = cellSubStart[i];
for (const sp of subpaths) {
  cellSubOffsets[subCursor++] = cursor; // cursor is the vertex index at sub-ring start
  sp.points.forEach(([x, y], k) => { /* existing body unchanged */ });
  if (sp.closed) path.closePath();
}
```

Add both to the returned object.

- [ ] **Step 4: Run to verify it passes**

Run: `npx vitest run tests/mapCache.test.ts`
Expected: PASS (both parity assertions and the multi-subpath guard).

- [ ] **Step 5: Typecheck + lint, then commit**

Run: `npm run typecheck && npm run lint`

```bash
git add utils/mapCache.ts tests/mapCache.test.ts
git commit -m "feat(F3): emit sub-ring boundary index from mapCache for tessellation"
```

---

### Task 2: `utils/tessellate.ts` — pure earcut triangulation

**Files:**
- Create: `utils/tessellate.ts`
- Create: `tests/tessellate.test.ts`
- Modify: `package.json` (add `earcut` + `@types/earcut`)

**Interfaces:**
- Consumes: `cellVerts: Float32Array`, `cellOffsets: Uint32Array`, `cellSubStart: Uint32Array`, `cellSubOffsets: Uint32Array` (Task 1).
- Produces:
  ```ts
  export interface CellTessellation {
    positions: Float32Array;   // [x,y, x,y, ...] triangle-list, CSS-px, un-flipped
    cellTriRange: Uint32Array; // length n+1: cell i's positions occupy vertex indices [cellTriRange[i], cellTriRange[i+1])
  }
  export const tessellateCells = (
    cellVerts: Float32Array, cellOffsets: Uint32Array,
    cellSubStart: Uint32Array, cellSubOffsets: Uint32Array,
    cellCount: number,
  ): CellTessellation
  ```
  Consumed by Tasks 3, 4, 7.

- [ ] **Step 1: Add the dependency**

Run: `npm install earcut && npm install -D @types/earcut`

- [ ] **Step 2: Write the failing test** — `tests/tessellate.test.ts`. Uses hand-built buffers (no world generation):

```ts
import { describe, it, expect } from 'vitest';
import { tessellateCells } from '../utils/tessellate';

// One square cell as a single sub-ring: (0,0)-(10,0)-(10,10)-(0,10).
function squareCell() {
  const verts = new Float32Array([0,0, 10,0, 10,10, 0,10]);
  const cellOffsets = new Uint32Array([0, 4]);
  const cellSubStart = new Uint32Array([0, 1]);
  const cellSubOffsets = new Uint32Array([0]);
  return { verts, cellOffsets, cellSubStart, cellSubOffsets };
}

describe('tessellateCells', () => {
  it('triangulates a convex quad into two triangles covering its area', () => {
    const { verts, cellOffsets, cellSubStart, cellSubOffsets } = squareCell();
    const t = tessellateCells(verts, cellOffsets, cellSubStart, cellSubOffsets, 1);
    const triVerts = t.cellTriRange[1] - t.cellTriRange[0];
    expect(triVerts).toBe(6); // 2 triangles
    // area of emitted triangles ~= 100
    let area = 0;
    for (let i = 0; i < t.positions.length; i += 6) {
      const [ax,ay,bx,by,cx,cy] = t.positions.slice(i, i+6);
      area += Math.abs((bx-ax)*(cy-ay) - (cx-ax)*(by-ay)) / 2;
    }
    expect(area).toBeCloseTo(100, 5);
  });

  it('triangulates a NON-convex (L-shaped) ring without triangles outside the polygon', () => {
    // L-shape: earcut must not fan across the concavity.
    const verts = new Float32Array([0,0, 10,0, 10,4, 4,4, 4,10, 0,10]);
    const cellOffsets = new Uint32Array([0, 6]);
    const cellSubStart = new Uint32Array([0, 1]);
    const cellSubOffsets = new Uint32Array([0]);
    const t = tessellateCells(verts, cellOffsets, cellSubStart, cellSubOffsets, 1);
    let area = 0;
    for (let i = 0; i < t.positions.length; i += 6) {
      const [ax,ay,bx,by,cx,cy] = t.positions.slice(i, i+6);
      area += Math.abs((bx-ax)*(cy-ay) - (cx-ax)*(by-ay)) / 2;
    }
    expect(area).toBeCloseTo(76, 4); // L area = 100 - 24; a bad fan would overshoot
  });

  it('emits an empty range for a cell with zero sub-rings (polar/degenerate)', () => {
    const verts = new Float32Array([]);
    const cellOffsets = new Uint32Array([0, 0]);
    const cellSubStart = new Uint32Array([0, 0]);
    const cellSubOffsets = new Uint32Array([]);
    const t = tessellateCells(verts, cellOffsets, cellSubStart, cellSubOffsets, 1);
    expect(t.cellTriRange[1] - t.cellTriRange[0]).toBe(0);
  });

  it('triangulates a cell with two disjoint sub-rings (antimeridian split) independently', () => {
    // two separate squares as two sub-rings of one cell
    const verts = new Float32Array([0,0, 2,0, 2,2, 0,2,  10,0, 12,0, 12,2, 10,2]);
    const cellOffsets = new Uint32Array([0, 8]);
    const cellSubStart = new Uint32Array([0, 2]);
    const cellSubOffsets = new Uint32Array([0, 4]);
    const t = tessellateCells(verts, cellOffsets, cellSubStart, cellSubOffsets, 1);
    expect(t.cellTriRange[1] - t.cellTriRange[0]).toBe(12); // 4 triangles total
  });
});
```

- [ ] **Step 3: Run to verify it fails**

Run: `npx vitest run tests/tessellate.test.ts`
Expected: FAIL — module not found.

- [ ] **Step 4: Implement `utils/tessellate.ts`**

```ts
import earcut from 'earcut';

export interface CellTessellation {
  positions: Float32Array;
  cellTriRange: Uint32Array;
}

// Earcut each cell's sub-rings independently and concatenate the resulting
// triangle-list positions. Cells are convex on the sphere but NOT after
// projection (Mercator polar stretch) or antimeridian clipping (a sub-ring
// bounded partly by the map edge), so a centroid fan would emit triangles
// outside the polygon exactly at cell corners — earcut avoids that. A cell
// whose clip produced no sub-rings (a pole under a cylindrical projection)
// contributes zero triangles.
export const tessellateCells = (
  cellVerts: Float32Array,
  cellOffsets: Uint32Array,
  cellSubStart: Uint32Array,
  cellSubOffsets: Uint32Array,
  cellCount: number,
): CellTessellation => {
  const out: number[] = [];
  const cellTriRange = new Uint32Array(cellCount + 1);
  let vertCount = 0;
  for (let i = 0; i < cellCount; i++) {
    cellTriRange[i] = vertCount;
    const subLo = cellSubStart[i];
    const subHi = cellSubStart[i + 1];
    for (let s = subLo; s < subHi; s++) {
      const ringStart = cellSubOffsets[s];             // vertex index
      const ringEnd = s + 1 < subHi ? cellSubOffsets[s + 1] : cellOffsets[i + 1];
      const count = ringEnd - ringStart;
      if (count < 3) continue;                          // degenerate sub-ring
      const coords: number[] = new Array(count * 2);
      for (let k = 0; k < count; k++) {
        coords[k * 2] = cellVerts[(ringStart + k) * 2];
        coords[k * 2 + 1] = cellVerts[(ringStart + k) * 2 + 1];
      }
      const tris = earcut(coords, undefined, 2);        // no holes; indices into `coords`
      for (const idx of tris) {
        out.push(coords[idx * 2], coords[idx * 2 + 1]);
        vertCount++;
      }
    }
  }
  cellTriRange[cellCount] = vertCount;
  return { positions: new Float32Array(out), cellTriRange };
};
```

- [ ] **Step 5: Run to verify it passes**

Run: `npx vitest run tests/tessellate.test.ts`
Expected: PASS (all four cases).

- [ ] **Step 6: Typecheck + lint, commit**

Run: `npm run typecheck && npm run lint`

```bash
git add utils/tessellate.ts tests/tessellate.test.ts package.json package-lock.json
git commit -m "feat(F3): earcut cell tessellation (utils/tessellate)"
```

---

### Task 3: `utils/mapFillColorBuffer.ts` — per-vertex colours

**Files:**
- Create: `utils/mapFillColorBuffer.ts`
- Create: `tests/mapFillColorBuffer.test.ts`

**Interfaces:**
- Consumes: the existing per-cell colour cache `string[]` of `#rrggbb` (from `buildCellColorCache`, `utils/mapColorCache.ts`) and `cellTriRange` (Task 2).
- Produces:
  ```ts
  export const buildFillColorBuffer = (
    perCellHex: string[], cellTriRange: Uint32Array, out?: Float32Array,
  ): Float32Array  // length = lastTriRange*3, RGB in 0..1, one triple per position vertex
  export const writeCellColors = (
    buf: Float32Array, perCellHex: string[], cellTriRange: Uint32Array, cellIds: Iterable<number>,
  ): void  // rewrites only the named cells' vertex ranges in place
  ```
  Consumed by Tasks 4, 7.

- [ ] **Step 1: Write the failing test** — `tests/mapFillColorBuffer.test.ts`:

```ts
import { describe, it, expect } from 'vitest';
import { buildFillColorBuffer, writeCellColors } from '../utils/mapFillColorBuffer';

describe('mapFillColorBuffer', () => {
  // cell 0 -> 3 verts (1 tri), cell 1 -> 6 verts (2 tris)
  const cellTriRange = new Uint32Array([0, 3, 9]);

  it('fills each cell vertex range with the cell colour as normalized RGB', () => {
    const hex = ['#ff0000', '#0080ff'];
    const buf = buildFillColorBuffer(hex, cellTriRange);
    expect(buf.length).toBe(9 * 3);
    // cell 0 vertices -> red
    expect(buf[0]).toBeCloseTo(1); expect(buf[1]).toBeCloseTo(0); expect(buf[2]).toBeCloseTo(0);
    // cell 1 first vertex (index 3) -> (0, 128/255, 1)
    expect(buf[3 * 3]).toBeCloseTo(0);
    expect(buf[3 * 3 + 1]).toBeCloseTo(128 / 255, 3);
    expect(buf[3 * 3 + 2]).toBeCloseTo(1);
  });

  it('writeCellColors rewrites only the named cell, leaving others untouched', () => {
    const hex = ['#ff0000', '#0080ff'];
    const buf = buildFillColorBuffer(hex, cellTriRange);
    const hex2 = ['#00ff00', '#0080ff'];
    writeCellColors(buf, hex2, cellTriRange, [0]);
    expect(buf[0]).toBeCloseTo(0); expect(buf[1]).toBeCloseTo(1); expect(buf[2]).toBeCloseTo(0);
    // cell 1 unchanged
    expect(buf[3 * 3 + 2]).toBeCloseTo(1);
  });
});
```

- [ ] **Step 2: Run to verify it fails**

Run: `npx vitest run tests/mapFillColorBuffer.test.ts`
Expected: FAIL — module not found.

- [ ] **Step 3: Implement `utils/mapFillColorBuffer.ts`**

```ts
// Expand the per-cell #rrggbb colour cache (utils/mapColorCache) into a
// per-vertex RGB buffer aligned to the tessellation's positions, so the GPU
// fill surface reads one colour per triangle vertex. Parsing the shared hex
// cache keeps a single source of truth for cell colour (no second getCellColor
// path that could drift from the Canvas2D fill).
const parseHex = (hex: string): [number, number, number] => {
  const n = parseInt(hex.slice(1), 16);
  return [((n >> 16) & 255) / 255, ((n >> 8) & 255) / 255, (n & 255) / 255];
};

export const buildFillColorBuffer = (
  perCellHex: string[],
  cellTriRange: Uint32Array,
  out?: Float32Array,
): Float32Array => {
  const total = cellTriRange[cellTriRange.length - 1];
  const buf = out && out.length === total * 3 ? out : new Float32Array(total * 3);
  for (let i = 0; i < perCellHex.length; i++) {
    const [r, g, b] = parseHex(perCellHex[i]);
    for (let v = cellTriRange[i]; v < cellTriRange[i + 1]; v++) {
      buf[v * 3] = r; buf[v * 3 + 1] = g; buf[v * 3 + 2] = b;
    }
  }
  return buf;
};

export const writeCellColors = (
  buf: Float32Array,
  perCellHex: string[],
  cellTriRange: Uint32Array,
  cellIds: Iterable<number>,
): void => {
  for (const i of cellIds) {
    const [r, g, b] = parseHex(perCellHex[i]);
    for (let v = cellTriRange[i]; v < cellTriRange[i + 1]; v++) {
      buf[v * 3] = r; buf[v * 3 + 1] = g; buf[v * 3 + 2] = b;
    }
  }
};
```

- [ ] **Step 4: Run to verify it passes**

Run: `npx vitest run tests/mapFillColorBuffer.test.ts`
Expected: PASS.

- [ ] **Step 5: Typecheck + lint, commit**

Run: `npm run typecheck && npm run lint`

```bash
git add utils/mapFillColorBuffer.ts tests/mapFillColorBuffer.test.ts
git commit -m "feat(F3): per-vertex fill colour buffer from the cell colour cache"
```

---

### Task 4: `GLFillSurface` — static GPU fill layer

Renders the tessellated fills with per-vertex colour under an orthographic camera. No pan/zoom yet — a fixed identity+mirror transform proving fills appear correct and right-way-round. Interaction arrives in Task 6.

**Files:**
- Create: `components/map2d/GLFillSurface.tsx`
- Verification: browser (rendering is verified in the browser per repo convention). No world-generating unit test.

**Interfaces:**
- Consumes: `positions` (Task 2), `colors` (Task 3), `size {width,height}`, `dpr`.
- Produces:
  ```ts
  export interface GLFillHandle {
    setTransform: (scale: number, offsetX: number, offsetY: number) => void; // used in Task 6
    setColors: (colors: Float32Array, range?: { start: number; count: number }) => void; // Task 7
    redraw: () => void;
    isContextLost: () => boolean;
  }
  export interface GLFillSurfaceProps {
    positions: Float32Array; colors: Float32Array;
    width: number; height: number; dpr: number;
    onContextLost?: () => void;
    handleRef?: React.MutableRefObject<GLFillHandle | null>;
  }
  export const GLFillSurface: React.FC<GLFillSurfaceProps>;
  ```

- [ ] **Step 1: Implement the component.** Key points, all load-bearing:
  - `WebGLRenderer({ canvas, antialias: true, alpha: true })` — `alpha:true` so the overlay/host background shows through where there are no fills; `antialias:true` for MSAA (Three's default is false → jaggy).
  - `OrthographicCamera(left=0, right=width, top=0, bottom=height, near=-1, far=1)` in CSS-px, so positions (CSS-px) map 1:1 at scale 1. Renderer `setPixelRatio(dpr)` and `setSize(width, height, false)`.
  - **Mirror:** the current 2D path flips X (`translate(width,0); scale(-1,1)`). Reproduce it as a mesh/camera transform: `mesh.scale.x = -1; mesh.position.x = width;` (or equivalently bake into the projection). Verify in-browser that coastlines match the blit path and are NOT left-right reversed.
  - `BufferGeometry` with `position` = `Float32BufferAttribute(positions, 2)` **plus a z of 0** — use a 2-component position attribute and an ortho camera in the z=0 plane (set `gl_Position.z` via the material; simplest is `new THREE.BufferAttribute(positions, 2)` with `MeshBasicMaterial({ vertexColors: true })` and a 2D→3D expansion, OR store positions as 3-component with z=0). Choose 3-component `[x,y,0]` to avoid a custom shader: expand once at upload.
  - `color` = `Float32BufferAttribute(colors, 3)`; `MeshBasicMaterial({ vertexColors: true })`.
  - `setTransform` updates `mesh.scale`/`mesh.position` (pan/zoom over the mirror) and calls `redraw`.
  - `setColors(colors, range?)`: assign to the color attribute's array (or copy the range), set `needsUpdate = true`; if `range`, call the Three r182 partial-upload API. **Pull the exact partial-update call (`addUpdateRange`/`clearUpdateRanges` vs legacy `updateRange`) from Context7 for `three@0.182.0` before writing this line — do not guess.**
  - `webglcontextlost` listener → `event.preventDefault()`, set an internal flag, call `onContextLost`. `isContextLost()` returns it.
  - **Disposal effect:** on unmount / prop-geometry change, `geometry.dispose()`, `material.dispose()`, `renderer.dispose()`, and remove listeners (WorldViewer invariant).
  - Type everything; no `any`. Use `THREE.` types for renderer/camera/mesh refs.

- [ ] **Step 2: Browser verification** — mount it behind a temporary dev toggle (or via Task 6 wiring later) at Default style + Mercator. Confirm: fills render in the right colours, coastline outline matches the blit path's shape, and text/east-west orientation is correct (compare a known coastline). Because this is pre-interaction, verify at scale 1 only.

- [ ] **Step 3: Typecheck + lint, commit**

Run: `npm run typecheck && npm run lint`

```bash
git add components/map2d/GLFillSurface.tsx
git commit -m "feat(F3): static Three.js GPU fill surface for the 2D map"
```

---

### Task 5: `VectorOverlay` — per-frame crisp linework + labels

A transparent Canvas2D layer over the GL canvas, redrawn each frame through the exact composed transform the current offscreen uses, so all existing draw functions run unchanged and land crisp. Adds zoom-settle boundary LOD.

**Files:**
- Create: `components/map2d/VectorOverlay.tsx`
- Reference (do not duplicate — import/reuse): the existing draw helpers used in `Map2D.tsx` (`drawMapLabels` from `utils/labels`, `drawRiverPaths`, route/contour/grid/current/marker passes) and `chainedToPath`-style boundary rebuild from `utils/mapCache` / `utils/boundaries`.
- Verification: browser.

**Interfaces:**
- Consumes: `size`, `dpr`, live `(scale, offset)`, the world + overlay toggles, and the cached boundary `Path2D`s (`coast`/`borders`/`lakes` from `MapGeometryCache`), plus a `boundaryTolPx` that the host lowers on settle.
- Produces:
  ```ts
  export interface VectorOverlayHandle { redraw: () => void; }
  export const VectorOverlay: React.FC<VectorOverlayProps>;
  ```

- [ ] **Step 1: Implement the overlay.** The composed transform, applied at the top of every `redraw` (this is the same order the offscreen uses, plus pan/zoom):

```ts
const ctx = canvas.getContext('2d')!;
ctx.setTransform(1, 0, 0, 1, 0, 0);
ctx.clearRect(0, 0, canvas.width, canvas.height);
// pan/zoom (display) ∘ dpr, matching the current blit's displayDpr*scale form:
ctx.setTransform(dpr * scale, 0, 0, dpr * scale, dpr * offsetX, dpr * offsetY);
// mirror (same as the offscreen build): flip X about width
ctx.translate(size.width, 0);
ctx.scale(-1, 1);
// now call the EXISTING draw functions unchanged — they already compensate
// text/glyphs for this flip (see docs/map-styles.md §mirror trap):
strokeBoundary(ctx, cache.coast, ...);      // coast/borders/lakes Path2D
drawRiverPaths(ctx, world.rivers, projection, 1.5 / (scale), overlayInk);
// contours, routes, grid, currents, markers — same calls as Map2D today
drawMapLabels(ctx, ...);                     // crisp text
```

  - Line widths that were previously divided by `qualityDpr` now divide by `scale` (the vector transform already carries dpr), so strokes keep constant screen weight as you zoom. Verify weights look right in-browser.
  - **Boundary LOD:** the `coast/borders/lakes` paths in the cache were simplified at `tol = 0.5` at base fit. Accept a `boundaryPaths` prop that the host swaps on zoom-settle for versions rebuilt at `tol = 0.5 / scale` (host rebuild in Task 6). Comment clearly: *this is settle-refined GEOMETRY (vector-tile LOD), not settle-refined pixels — the overlay is resolution-independent at every zoom; only the simplification budget tracks zoom.*

- [ ] **Step 2: Browser verification** — coastlines, borders, rivers, routes, labels all crisp at 1× and 6×; text forward-facing; stroke weights constant across zoom; nothing haloed off the fills.

- [ ] **Step 3: Typecheck + lint, commit**

Run: `npm run typecheck && npm run lint`

```bash
git add components/map2d/VectorOverlay.tsx
git commit -m "feat(F3): per-frame Canvas2D vector overlay with zoom-settle boundary LOD"
```

---

### Task 6: Wire Map2D — layer selection, transform drive, fallback

**Files:**
- Modify: `components/Map2D.tsx` (add the GL+overlay path; keep the blit path intact as fallback/out-of-scope route)
- Create: `utils/map2dRenderPath.ts` (pure selection predicate) + `tests/map2dRenderPath.test.ts`

**Interfaces:**
- Consumes: Tasks 2–5 outputs and existing Map2D state (`scale`, `offset`, `size`, `projection`, `projectionType`, `style`, `geometryCache`, `colorCache`, `world`).
- Produces: `export const shouldUseGLFill = (o: { styleId: string; projectionType: string; webglAvailable: boolean; contextLost: boolean }) => boolean`.

- [ ] **Step 1: Write the failing test** — `tests/map2dRenderPath.test.ts`:

```ts
import { describe, it, expect } from 'vitest';
import { shouldUseGLFill } from '../utils/map2dRenderPath';

describe('shouldUseGLFill', () => {
  const base = { styleId: 'default', projectionType: 'mercator', webglAvailable: true, contextLost: false };
  it('uses GL for default style on a non-dymaxion projection with webgl', () => {
    expect(shouldUseGLFill(base)).toBe(true);
  });
  it('falls back for parchment', () => {
    expect(shouldUseGLFill({ ...base, styleId: 'parchment' })).toBe(false);
  });
  it('falls back for dymaxion', () => {
    expect(shouldUseGLFill({ ...base, projectionType: 'dymaxion' })).toBe(false);
  });
  it('falls back when webgl unavailable', () => {
    expect(shouldUseGLFill({ ...base, webglAvailable: false })).toBe(false);
  });
  it('falls back after context loss', () => {
    expect(shouldUseGLFill({ ...base, contextLost: true })).toBe(false);
  });
});
```

- [ ] **Step 2: Run to verify it fails**

Run: `npx vitest run tests/map2dRenderPath.test.ts`
Expected: FAIL — module not found.

- [ ] **Step 3: Implement `utils/map2dRenderPath.ts`**

```ts
export const shouldUseGLFill = (o: {
  styleId: string; projectionType: string; webglAvailable: boolean; contextLost: boolean;
}): boolean =>
  o.styleId === 'default' &&
  o.projectionType !== 'dymaxion' &&
  o.webglAvailable &&
  !o.contextLost;
```

- [ ] **Step 4: Run to verify it passes**

Run: `npx vitest run tests/map2dRenderPath.test.ts`
Expected: PASS.

- [ ] **Step 5: Wire into `Map2D.tsx`.** Behind `shouldUseGLFill`:
  - Compute `webglAvailable` once (try creating a probe context; memoize). Track `contextLost` state; `GLFillSurface.onContextLost` sets it → selection flips to blit on next render.
  - `tessellation = useMemo(() => tessellateCells(cache.cellVerts, cache.cellOffsets, cache.cellSubStart, cache.cellSubOffsets, world.cells.length), [cache])`.
  - `fillColors = useMemo(() => buildFillColorBuffer(colorCache, tessellation.cellTriRange), [colorCache, tessellation])`.
  - Render `<GLFillSurface positions=… colors=… width height dpr … handleRef=glRef onContextLost=…/>` and `<VectorOverlay … handleRef=ovRef/>` **instead of** the offscreen `<canvas>` blit when `shouldUseGLFill`. Otherwise render the existing canvas/blit untouched.
  - **DPR on the GL path is fixed** `dpr = Math.min(2, devicePixelRatio || 1)` — the `qualityDpr`/`MAX_SHARP_*` ramp is bypassed here (resolution is zoom-independent now). Leave that machinery in place for the blit path.
  - Drive interaction: in the existing wheel/drag handlers, when on the GL path, call `glRef.setTransform(scale, offset.x, offset.y)` + `glRef.redraw()` and `ovRef.redraw()` each frame instead of blitting. No geometry rebuild.
  - **Boundary LOD:** on zoom-settle (reuse the existing `settleTimer`), rebuild `coast/borders/lakes` at `tol = 0.5 / scale` (call the same `chainedToPath` path builder `mapCache` uses, exported for reuse) and hand them to `VectorOverlay`.
  - **Two canvases agree to the pixel:** size both from the same `size`×`dpr`; position them in the same stacking container.

- [ ] **Step 6: Browser verification (the big one)** — at Default + Mercator: crisp fills/lines/labels at 1× and 6×; pan/zoom smooth; **switch to Parchment → blit path (unchanged look); switch to Dymaxion → unchanged; simulate context loss (`renderer.forceContextLoss()` via devtools) → falls back to blit without a blank frame**. Pick still selects the right cell at zoom (mapPick unaffected).

- [ ] **Step 7: Typecheck + lint, commit**

Run: `npm run typecheck && npm run lint`

```bash
git add components/Map2D.tsx utils/map2dRenderPath.ts tests/map2dRenderPath.test.ts utils/mapCache.ts
git commit -m "feat(F3): wire GL fill + vector overlay into Map2D with blit fallback"
```

---

### Task 7: Recolour + paint partial-update path

**Files:**
- Modify: `components/Map2D.tsx` (recolour effects), `components/map2d/GLFillSurface.tsx` (already exposes `setColors`)
- Verification: browser + reuse the pure tests from Task 3.

**Interfaces:** consumes `writeCellColors`, `buildFillColorBuffer`, `GLFillHandle.setColors`.

- [ ] **Step 1: View-mode / season / seaLevel / hillshade recolour.** When `colorCache` changes (its memo already keys on viewMode/season/seaLevel/colormaps/shadeMap), rebuild the fill colour buffer into the SAME `Float32Array` (`buildFillColorBuffer(colorCache, cellTriRange, existingBuffer)`) and call `glRef.setColors(buffer)` (full upload) + `redraw`. Topology (positions) untouched.

- [ ] **Step 2: Paint-stroke partial update.** In the edit/paint path, after a stroke mutates cells in place, call `writeCellColors(buffer, freshHexForCells, cellTriRange, changedCellIds)` then `glRef.setColors(buffer, { start: cellTriRange[minCell]*3, count: … })` using the affected contiguous range where possible, else per-cell ranges. Do NOT re-upload the whole attribute per stroke. **Confirm the r182 partial-upload API via Context7.**

- [ ] **Step 3: Browser verification** — switch biomes→political→height: instant, correct recolour, no geometry flicker. Paint a stroke: only painted cells change; undo restores; performance stays smooth.

- [ ] **Step 4: Typecheck + lint, commit**

Run: `npm run typecheck && npm run lint`

```bash
git add components/Map2D.tsx components/map2d/GLFillSurface.tsx
git commit -m "feat(F3): GPU recolour on view-mode change + paint partial buffer update"
```

---

### Task 8: Perf pass at 200k

**Files:** likely `components/Map2D.tsx`, possibly `workers/worldGen.worker.ts` / `utils/worldGenClient.ts` if tessellation moves to the worker. Create `scripts/renderMap2DPerf.mjs` additions only if needed.

- [ ] **Step 1: Measure** build time (tessellate + colour buffer + upload) and steady-state pan/zoom frame time at `points=200000`, Default + Mercator, dpr 2, using the existing `scripts/renderMap2DPerf.mjs` harness (extend it to the GL path). Record numbers in HANDOFF.
- [ ] **Step 2: Decide worker offload.** If main-thread tessellation of 200k stalls interaction (> ~1 frame of jank on generate), move `tessellateCells` into the generation worker (it is pure and its outputs are transferable `Float32Array`/`Uint32Array` — transfer, don't copy). Keep it on the main thread if build time is negligible. Document the decision + measured threshold in HANDOFF (n and evidence).
- [ ] **Step 3: Overlay cadence guard (only if needed).** If the per-frame overlay drops frames at 200k during active drag, redraw it at a reduced cadence during drag and sharp on settle (GL fills stay live). Do not add this unless measurement shows it is needed — YAGNI.
- [ ] **Step 4: Commit** whatever changes the measurement justified (may be none beyond the HANDOFF note).

```bash
git add -p   # only the files the perf work touched
git commit -m "perf(F3): 200k tessellation/upload measured; <worker offload decision>"
```

---

### Task 9: Verification, docs, and merge prep

**Files:** `ROADMAP.md`, `HANDOFF.md`, `docs/map-styles.md` (note the two-layer Default path), `CLAUDE.md` invariants if a new one is load-bearing.

- [ ] **Step 1: Full browser sweep** — every Default view mode at 1× and 6×, both non-Dymaxion projections; parchment/dymaxion/3D untouched; context-loss fallback; pick at zoom; PNG/SVG export **unchanged** (exports do not use the GL path — confirm they still route through the existing renderer).
- [ ] **Step 2: Docs.** ROADMAP F3 Phase 2 → ✅ DONE with the measured sharpness result; HANDOFF session entry (decisions + measured perf, at honest confidence); a `docs/map-styles.md` note that Default 2D now renders as GL fills + Canvas2D overlay while Parchment stays on the blit path.
- [ ] **Step 3: Gates.** `npm run typecheck && npm run lint` (must stay ≤30 warnings, 0 errors). Then **ask Matt** before the full suite, offering all three speeds (full/capped/single) — never run it alongside browser automation.
- [ ] **Step 4: Final commit + advisor review** of the accumulated diff before declaring done.

```bash
git add ROADMAP.md HANDOFF.md docs/map-styles.md
git commit -m "docs(F3): mark Phase 2 done; record two-layer Default render path + perf"
```

---

## Self-Review

**Spec coverage:** GL fills (T4), overlay+labels (T5), subpath index (T1), earcut (T2), per-vertex colour (T3), boundary LOD (T5/T6), fallback on init+context-loss (T6), recolour+paint partial update (T7), 200k perf/worker decision (T8), mirror correctness (T4/T5 browser gates), antialias (T4), two-canvas pixel agreement (T6), disposal (T4), pick unchanged (T6), exports unchanged (T9), lint ratchet (every task's gate). All spec sections map to a task.

**Placeholder scan:** no "TBD"/"handle edge cases"/"write tests for the above" — pure units carry real test + impl code; component tasks are honestly marked browser-verified per repo convention with concrete implementation notes and the two Context7-verification points (Three partial-upload API) called out rather than guessed.

**Type consistency:** `cellSubStart`/`cellSubOffsets` (T1) → consumed verbatim by `tessellateCells` (T2) → `cellTriRange`/`positions` consumed by `buildFillColorBuffer`/`writeCellColors` (T3), `GLFillSurface` (T4), Map2D wiring (T6/T7). `GLFillHandle.setColors(colors, range?)` defined T4, used T7. `shouldUseGLFill` signature defined and consumed in T6. Names consistent across tasks.
