# F2 draped graticule — implementation plan

Spec: `docs/superpowers/specs/2026-08-15-f2-drape-graticule-design.md`
Branch: `f2-drape-graticule`

## Task 1 — `nearestCellWalk` helper + unit test (self-contained)

**File:** `utils/nearestCell.ts` (new)

```ts
import { Cell } from '../types';

// Greedy hill-climb over dot(dir, cell.center) across the symmetric Voronoi
// neighbor graph → the display cell nearest to a unit direction. Seed `startId`
// with the previous sample for O(1)-ish walks along a polyline; brute-force the
// very first sample. Cap defends against a degenerate graph (dot is strictly
// increasing on a real Voronoi graph, so the cap never binds).
export function nearestCellWalk(
  cells: Cell[], dx: number, dy: number, dz: number, startId: number,
): number {
  let cur = startId >= 0 && startId < cells.length ? startId : 0;
  let curDot = cells[cur].center.x * dx + cells[cur].center.y * dy + cells[cur].center.z * dz;
  for (let iter = 0; iter < cells.length; iter++) {
    let best = cur, bestDot = curDot;
    for (const n of cells[cur].neighbors) {
      const c = cells[n].center;
      const d = c.x * dx + c.y * dy + c.z * dz;
      if (d > bestDot) { bestDot = d; best = n; }
    }
    if (best === cur) break;
    cur = best; curDot = bestDot;
  }
  return cur;
}

// One-time seed for the first sample of a draw (O(n)).
export function nearestCellBrute(cells: Cell[], dx: number, dy: number, dz: number): number {
  let best = 0, bestDot = -Infinity;
  for (let i = 0; i < cells.length; i++) {
    const c = cells[i].center;
    const d = c.x * dx + c.y * dy + c.z * dz;
    if (d > bestDot) { bestDot = d; best = i; }
  }
  return best;
}
```

**Test** (`tests/nearestCell.test.ts`): build a tiny cell array (e.g. 6 cells at
known directions with symmetric `neighbors`), assert `nearestCellWalk` returns the
true argmax-dot from several `startId`s (including the far side), and that it
equals `nearestCellBrute` for random directions on a slightly larger synthetic
graph. Keep it pure — construct minimal `Cell`-shaped objects (only `center` and
`neighbors` are read).

**Gate:** `npx vitest run tests/nearestCell.test.ts`, typecheck.

→ **Delegate to OpenCode first** (self-contained; 2-min wall-clock cap/attempt,
3 attempts). On repeated failure → Sonnet.

## Task 2 — drape the graticule (integration; in-house or Sonnet, NOT OpenCode)

**File:** `components/overlays/tenants.ts` `drawGraticuleTenant`

- Import `nearestCellWalk`, `nearestCellBrute`.
- Keep the `smooth` path exactly as-is (R = 1).
- Raised path (`!smooth`): replace the fixed `R = 1 + seaLevel*0.05` with a
  per-sample draped radius. Maintain one `let startId` across the whole function;
  seed it with `nearestCellBrute` on the first projected sample. For each sample
  direction `(x,y,z)` (already unit — the lat/lon points are on the unit sphere):
  `startId = nearestCellWalk(cells, x, y, z, startId);`
  `const h = Math.max(cells[startId].height, world.params.seaLevel);`
  `const r = 1 + h * 0.05;` then pass `(x*r, y*r, z*r)` to `project`.
- `cells = world.cells`.

Note the parallels/meridians build points as `(cr*cos, cy, cr*sin)` etc. — those
are unit directions when R=1; keep them unit for the walk, apply `r` only when
handing to `project`. Simplest: compute the unit dir, walk, then scale.

**Gate:** typecheck, lint, build.

## Task 3 — verify + finish (in-house)

- Confirm the redraw gate covers paint (ScreenOverlay `[nCells, world]` resets
  lastKey) — read, don't add unless broken.
- In-browser (raised globe, grid on, smooth toggled OFF): coastline meets terrain;
  grid ripples over mountains; ocean limb unchanged; 0 console errors; perf fine.
- Full gates: typecheck 0, lint 0 errors, `npm test`, build.
- HANDOFF Session 18 entry; note the "drape lands ⇒ Matt may drop grid→smooth
  coupling" option. Commit; leave branch for Matt to merge (per his call).
