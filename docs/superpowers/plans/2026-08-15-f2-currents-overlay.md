# F2 Screen-Space Overlay + Currents Visualization — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Persist the D2 ocean-current field on `WorldData` and draw it as a toggleable screen-space overlay on the globe and 2D map, on a new generic `ScreenOverlay` layer that also absorbs the graticule (proving the abstraction).

**Architecture:** The velocity field + SST anomaly (already computed in the worker, currently discarded) are kept on `WorldData.currents` (optional, export-excluded). A new `ScreenOverlay` component inside the R3F `<Canvas>` projects visible cells to a sibling 2D canvas using an analytic horizon test (no depth-buffer readback), and dispatches to tenant draw callbacks. v1 tenants: currents (arrows + warm/cold tint) and the graticule (migrated off its 3D `LatLongGrid`). The 2D `Map2D` draws currents in its existing composite pass.

**Tech Stack:** React 19, three.js / react-three-fiber, Canvas2D, TypeScript (strict), Vitest, Web Worker generation.

**Spec:** `docs/superpowers/specs/2026-08-15-f2-currents-overlay-design.md`

## Global Constraints

- **Relative imports only** — the `@/` alias is configured but intentionally unused.
- **Determinism is a test invariant** — no RNG in this feature; same seed → byte-identical `currents` arrays.
- **`currents` is optional on `WorldData`** — present iff `params.currentStrength > 0`; absent (undefined) otherwise. This preserves D2's byte-identical `currentStrength === 0` escape hatch.
- **`currents` is never JSON-serialized** — strip before save in `export.ts`; regenerates on load.
- **Rendering is verified manually in the browser** (project convention) — pure logic is unit-tested with Vitest; React/canvas glue gets structured in-browser (Playwright) verification steps.
- **Horizon test epsilon: ε = 0.08** (absorbs the ≤1.05× relief grazing edge).
- **Gates before any task's final commit:** `npm run typecheck` (0), `npm run lint` (≤30 warnings ratchet), relevant `npm test` green. Watch the documented M1 parallel-load flake — re-run a suspected failure in isolation before believing it.
- **Delegation:** Tasks marked **[DELEGABLE]** may be sent to a Claude subagent (Sonnet) or `opencode-go/deepseek-v4-flash` against this plan. **DeepSeek output gets a full review pass, not a spot-check** (different model + harness) — read the entire diff, run typecheck+lint+tests locally, and confirm behavior before committing.

---

## File Structure

- **`types.ts`** (modify) — add `OceanCurrentData` interface + optional `currents?` on `WorldData` (near line 263).
- **`utils/worldGen.ts`** (modify) — stop discarding `currentField`/`sstAnomaly` (computed ~542–543); attach to `world` at assembly (line 626).
- **`utils/export.ts`** (modify) — strip `currents` before serialize on save.
- **`utils/screenProject.ts`** (create) — pure screen-projection + analytic horizon-test helpers (no React, no THREE component types). Unit-tested.
- **`components/overlays/ScreenOverlay.tsx`** (create) — the R3F-Canvas-child screen-space overlay layer + tenant dispatch.
- **`components/overlays/tenants.ts`** (create) — pure tenant draw functions: `drawCurrentsTenant`, `drawGraticuleTenant`.
- **`components/WorldViewer.tsx`** (modify) — mount `ScreenOverlay`; remove the 3D `LatLongGrid` graticule (keep `TiltAxisLine`); thread `showCurrents`.
- **`components/Map2D.tsx`** (modify) — draw currents in the composite pass.
- **`hooks/useWorldEngine.ts`** (modify) — `showCurrents` state (mirror `showContours`, lines 78–82 + return object line 655).
- **`components/Controls.tsx`** (modify) — "Currents" toggle (mirror the Contours/Hillshade toggle).
- **`tests/currentsPersistence.test.ts`** (create) — determinism + presence/absence.
- **`tests/screenProject.test.ts`** (create) — horizon test + projection math.

---

## Task 1: Persist the current field on `WorldData` — [DELEGABLE]

**Files:**
- Modify: `types.ts:263` (WorldData) + new interface
- Modify: `utils/worldGen.ts:539-543, 626`
- Modify: `utils/export.ts` (save serialization path)
- Test: `tests/currentsPersistence.test.ts`

**Interfaces:**
- Produces: `interface OceanCurrentData { vx: Float32Array; vy: Float32Array; vz: Float32Array; sst: Float32Array }` and `WorldData.currents?: OceanCurrentData`.

- [ ] **Step 1: Write the failing test**

```ts
// tests/currentsPersistence.test.ts
import { describe, it, expect } from 'vitest';
import { generateWorld } from '../utils/worldGen';
import { makeParams } from './helpers';

describe('WorldData.currents persistence', () => {
  it('present and deterministic when currentStrength > 0', async () => {
    const p = makeParams({ seed: 'currents-persist', currentStrength: 1.0 });
    const a = await generateWorld(p);
    const b = await generateWorld(p);
    expect(a.currents).toBeDefined();
    expect(a.currents!.vx.length).toBe(a.cells.length);
    expect(a.currents!.sst.length).toBe(a.cells.length);
    // byte-identical reruns
    expect(Array.from(a.currents!.vx)).toEqual(Array.from(b.currents!.vx));
    expect(Array.from(a.currents!.sst)).toEqual(Array.from(b.currents!.sst));
  });

  it('absent when currentStrength === 0 (escape hatch)', async () => {
    const p = makeParams({ seed: 'currents-persist', currentStrength: 0 });
    const w = await generateWorld(p);
    expect(w.currents).toBeUndefined();
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npx vitest run tests/currentsPersistence.test.ts`
Expected: FAIL — `a.currents` is `undefined` (field not yet attached).

- [ ] **Step 3: Add the type** — in `types.ts`, above `interface WorldData` (line 263):

```ts
export interface OceanCurrentData {
  vx: Float32Array;   // per-cell tangential velocity (unit-sphere frame)
  vy: Float32Array;
  vz: Float32Array;
  sst: Float32Array;  // per-cell SST anomaly, °C (ocean cells; 0 elsewhere)
}
```
Then add to `WorldData`: `currents?: OceanCurrentData;`

- [ ] **Step 4: Attach in worldGen** — `utils/worldGen.ts`. The block at ~539–543 computes `currentField` / `sstAnomaly` inside the `currentStrength > 0` guard. At the world assembly (line 626), attach them:

```ts
const world: WorldData = { cells, params, geoJson: polygons, rivers, lakes };
if (currentField && sstAnomaly) {
  world.currents = {
    vx: currentField.vx, vy: currentField.vy, vz: currentField.vz, sst: sstAnomaly,
  };
}
```
(When `currentStrength === 0`, both stay `null` → `currents` omitted.)

- [ ] **Step 5: Exclude from save** — `utils/export.ts`, in the world→JSON serialization: ensure the serialized object never includes `currents` (it is not part of the save schema; only `params`+`civData`+`markers` are persisted, geometry regenerates). If the save builds an explicit object, do not add `currents`. If it spreads `world`, destructure it out: `const { currents: _drop, ...saveable } = world;` and serialize `saveable`. Verify no `Float32Array` leaks into the JSON.

- [ ] **Step 6: Run tests + gates**

Run: `npx vitest run tests/currentsPersistence.test.ts tests/paramLiveness.test.ts && npm run typecheck`
Expected: PASS (currents present/deterministic/absent-at-0); paramLiveness still green (no new case needed — `currentStrength` is already exercised via temperature/moisture).

- [ ] **Step 7: Commit**

```bash
git add types.ts utils/worldGen.ts utils/export.ts tests/currentsPersistence.test.ts
git commit -m "feat(F2): persist ocean-current field on WorldData (optional, export-excluded)"
```

---

## Task 2: Pure screen-projection + horizon-test helper

**Files:**
- Create: `utils/screenProject.ts`
- Test: `tests/screenProject.test.ts`

**Interfaces:**
- Produces:
  - `isVisible(px, py, pz, camx, camy, camz, eps=0.08): boolean` — analytic horizon test. The point `p` (unit-sphere cell center) is visible iff the vector from `p` to the camera has a positive-enough dot with the surface normal (= `p` itself on the unit sphere): `dot(normalize(cam − p), p) > eps`.
  - `ProjectedCells { x: Float32Array; y: Float32Array; visible: Uint8Array; n: number }` (screen-pixel coords; `visible[i]` = 1/0).

- [ ] **Step 1: Write the failing test**

```ts
// tests/screenProject.test.ts
import { describe, it, expect } from 'vitest';
import { isVisible } from '../utils/screenProject';

describe('analytic horizon test', () => {
  const cam = [0, 0, 3]; // camera on +Z
  it('front-facing point (near +Z) is visible', () => {
    expect(isVisible(0, 0, 1, cam[0], cam[1], cam[2])).toBe(true);
  });
  it('back-facing point (near −Z) is culled', () => {
    expect(isVisible(0, 0, -1, cam[0], cam[1], cam[2])).toBe(false);
  });
  it('grazing point at the horizon (near +X) is culled within eps', () => {
    // p=(1,0,0): cam−p=(-1,0,3) normalized dot p = -1/sqrt(10) < 0.08 → false
    expect(isVisible(1, 0, 0, cam[0], cam[1], cam[2])).toBe(false);
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npx vitest run tests/screenProject.test.ts`
Expected: FAIL — module/function not found.

- [ ] **Step 3: Implement**

```ts
// utils/screenProject.ts
export interface ProjectedCells {
  x: Float32Array; y: Float32Array; visible: Uint8Array; n: number;
}

export function isVisible(
  px: number, py: number, pz: number,
  camx: number, camy: number, camz: number,
  eps = 0.08,
): boolean {
  let dx = camx - px, dy = camy - py, dz = camz - pz;
  const len = Math.hypot(dx, dy, dz);
  if (len < 1e-9) return true;
  dx /= len; dy /= len; dz /= len;
  // surface normal at a unit-sphere point IS the point itself
  return dx * px + dy * py + dz * pz > eps;
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `npx vitest run tests/screenProject.test.ts`
Expected: PASS (3/3).

- [ ] **Step 5: Commit**

```bash
git add utils/screenProject.ts tests/screenProject.test.ts
git commit -m "feat(F2): pure screen-projection + analytic horizon-test helper"
```

---

## Task 3: `ScreenOverlay` layer (R3F-Canvas child)

**Files:**
- Create: `components/overlays/ScreenOverlay.tsx`
- Create: `components/overlays/tenants.ts` (stubs this task; filled in Tasks 4–5)
- Modify: `components/WorldViewer.tsx` (mount inside `<Canvas>`, near where `LatLongGrid` renders ~1114)

**Interfaces:**
- Consumes: `isVisible`, `ProjectedCells` (Task 2); `WorldData.currents` (Task 1).
- Produces:
  - `interface OverlayTenant { id: string; visible: boolean; draw(ctx: CanvasRenderingContext2D, proj: ProjectedCells, world: WorldData): void; }`
  - `<ScreenOverlay world={WorldData} tenants={OverlayTenant[]} />` — mounted inside `<Canvas>`.

**Design notes (no unit test — in-browser verified per convention):**
- Renders a sibling `<canvas>` absolutely positioned over the WebGL canvas via a portal to the R3F container, `pointer-events: none`, sized `clientWidth*DPR × clientHeight*DPR`, CSS-scaled back down.
- Uses `useThree()` for `camera`, `gl`, `size`; a `useFrame` callback that:
  1. early-returns unless the camera matrix OR `world`/tenant identity changed since last draw (store last `camera.matrixWorld.elements` hash / a version ref);
  2. fills `ProjectedCells` once: for each cell, `isVisible(center, cameraPos)`; if visible, project `center` via `new THREE.Vector3(cx,cy,cz).project(camera)` → NDC → pixel `(x+1)/2*W`, `(1-(y+1)/2)*H`;
  3. clears the 2D ctx and calls each visible tenant's `draw(ctx, proj, world)`.
- Cleanup effect: remove the canvas node on unmount. (No THREE geometry is created, so the "every useMemo geometry needs a disposal effect" invariant is N/A here — a deliberate simplification.)

- [ ] **Step 1: Create `tenants.ts` stubs**

```ts
// components/overlays/tenants.ts
import { WorldData } from '../../types';
import { ProjectedCells } from '../../utils/screenProject';

export function drawCurrentsTenant(ctx: CanvasRenderingContext2D, proj: ProjectedCells, world: WorldData): void {
  // filled in Task 4
}
export function drawGraticuleTenant(ctx: CanvasRenderingContext2D, proj: ProjectedCells, world: WorldData): void {
  // filled in Task 5
}
```

- [ ] **Step 2: Implement `ScreenOverlay.tsx`** per the design notes above. Export `OverlayTenant`. The projection loop must be allocation-light (reuse a scratch `Vector3`, not one-per-cell).

- [ ] **Step 3: Mount in `WorldViewer`** — inside `<Canvas>`, add `{world && <ScreenOverlay world={world} tenants={overlayTenants} />}` where `overlayTenants` is a `useMemo` array (empty for now; Tasks 4–5 populate based on `showCurrents`/`showGrid`).

- [ ] **Step 4: Typecheck + lint**

Run: `npm run typecheck && npm run lint`
Expected: 0 errors.

- [ ] **Step 5: In-browser smoke** (dev server on :3000, reuse if running): the app renders, globe rotates, 0 console errors, no visual change yet (tenants empty). Confirm the overlay `<canvas>` exists in the DOM over the WebGL canvas.

- [ ] **Step 6: Commit**

```bash
git add components/overlays/ScreenOverlay.tsx components/overlays/tenants.ts components/WorldViewer.tsx
git commit -m "feat(F2): ScreenOverlay screen-space layer + tenant dispatch (no tenants yet)"
```

---

## Task 4: Currents tenant (arrows + SST tint) on the globe

**Files:**
- Modify: `components/overlays/tenants.ts` (`drawCurrentsTenant`)
- Modify: `components/WorldViewer.tsx` (include the currents tenant when `showCurrents && world.currents`)

**Interfaces:**
- Consumes: `world.currents` (`vx,vy,vz,sst`), `ProjectedCells`.

**Implementation:**
- For each visible **ocean** cell (`world.currents` present): project both the cell center AND a second point offset along the velocity `(vx,vy,vz)` by a small arc, take the screen-space delta as the arrow direction; length ∝ speed (clamp). Skip cells with speed below a threshold. Subsample for legibility — first pass: draw a cell only if `speed > SPEED_MIN` (start `0.05`) and stride every Nth by id; tune in-browser.
- Color: lerp cold `#3b6fb0` ↔ warm `#c0442e` by `sst` (clamp `sst` to ~[−6, 6] °C). Arrowheads optional in v1 (a short line + a 2px dot at the head is fine).

- [ ] **Step 1: Implement `drawCurrentsTenant`** per above (guard `if (!world.currents) return;`).
- [ ] **Step 2: Wire the tenant** in `WorldViewer` — add `{ id: 'currents', visible: !!showCurrents && !!world.currents, draw: drawCurrentsTenant }` to `overlayTenants`.
- [ ] **Step 3: Typecheck + lint** — `npm run typecheck && npm run lint` → 0.
- [ ] **Step 4: In-browser verify** (Playwright, 5k cells, seed `realmgenesis`, `currentStrength` default 1.0): toggle Currents on → arrows appear over oceans on the **near** hemisphere only (far side occluded); warm/cold tint visible; rotating the globe re-culls correctly; toggle off → arrows gone; 0 console errors. Capture a screenshot to the session scratchpad.
- [ ] **Step 5: Commit**

```bash
git add components/overlays/tenants.ts components/WorldViewer.tsx
git commit -m "feat(F2): currents tenant — arrows + warm/cold SST tint on globe"
```

---

## Task 5: Migrate the graticule onto `ScreenOverlay`

**Files:**
- Modify: `components/overlays/tenants.ts` (`drawGraticuleTenant`)
- Modify: `components/WorldViewer.tsx` — remove the 3D `LatLongGrid` (line 1114 render + the `LatLongGrid` component ~1282); **keep `TiltAxisLine`** (1115); add the graticule tenant.

**Implementation:**
- Draw 10° meridians and parallels as 2D polylines: sample each line at N points on the unit sphere, project each via the same `isVisible` + projection, and `moveTo/lineTo` across runs of visible points (break the path where a point is culled, so lines don't cross the globe silhouette). Thin white, ~0.35 alpha (match the old grid's look).

- [ ] **Step 1: Implement `drawGraticuleTenant`** per above.
- [ ] **Step 2: Remove `LatLongGrid`** 3D component + its render site; add `{ id: 'graticule', visible: !!showGrid, draw: drawGraticuleTenant }` to `overlayTenants`. Delete the now-unused `LatLongGrid` geometry/disposal.
- [ ] **Step 3: Typecheck + lint** → 0.
- [ ] **Step 4: In-browser verify:** toggle Grid on → graticule reads correctly, occludes on the far hemisphere (an improvement over the old always-visible 3D grid), rotates correctly; `TiltAxisLine` still renders; 0 console errors.
- [ ] **Step 5: Commit**

```bash
git add components/overlays/tenants.ts components/WorldViewer.tsx
git commit -m "feat(F2): migrate graticule off 3D LatLongGrid onto ScreenOverlay"
```

---

## Task 6: Currents on the 2D map — [DELEGABLE]

**Files:**
- Modify: `components/Map2D.tsx` (composite pass)

**Implementation:**
- In the composite/overlay pass, when `showCurrents && world.currents`, draw arrows over ocean cells using `Map2D`'s existing sphere→screen projection (Mercator/Dymaxion) — the same one used for labels. Same tint/speed rules as Task 4. This does NOT use `ScreenOverlay` (the 2D map is already screen-space and has its own projection).

- [ ] **Step 1: Locate the composite pass** and the existing per-cell projection helper in `Map2D.tsx`; add a `drawCurrents2D` sub-pass gated on `showCurrents && world.currents`.
- [ ] **Step 2: Implement** the arrow draw (reuse tint/speed constants from `tenants.ts` — export them there to keep a single source).
- [ ] **Step 3: Typecheck + lint** → 0.
- [ ] **Step 4: In-browser verify:** switch to 2D (Mercator), toggle Currents → field draws over oceans, warm/cold tint; Dymaxion path also draws; 0 console errors.
- [ ] **Step 5 (if delegated to DeepSeek): FULL review** — read the entire diff, run `npm run typecheck && npm run lint && npx vitest run`, confirm no unrelated edits, verify in-browser before committing.
- [ ] **Step 6: Commit**

```bash
git add components/Map2D.tsx components/overlays/tenants.ts
git commit -m "feat(F2): draw currents field on the 2D map composite pass"
```

---

## Task 7: `showCurrents` state + Controls toggle

**Files:**
- Modify: `hooks/useWorldEngine.ts` (state + return object)
- Modify: `components/Controls.tsx` (toggle UI)
- Modify: `components/WorldViewer.tsx` + `ShellApp`/prop chain (thread `showCurrents` prop, mirroring `showContours`)

**Implementation — mirror `showContours` exactly:**

- [ ] **Step 1:** In `useWorldEngine.ts`, after line 82 (`const [showContours, setShowContours] = useState(false);`): `const [showCurrents, setShowCurrents] = useState(false);`. Add `showCurrents, setShowCurrents` to the returned object (line ~655).
- [ ] **Step 2:** Thread `showCurrents` through the prop chain to `WorldViewer` and `Map2D` (follow every place `showContours` is passed — same call sites).
- [ ] **Step 3:** In `Controls.tsx`, add a "Currents" toggle beside Contours/Hillshade. Disable it (greyed + hint "requires ocean currents on") when `world?.currents` is absent (i.e. `currentStrength === 0`).
- [ ] **Step 4: Typecheck + lint** → 0.
- [ ] **Step 5: In-browser verify:** toggle drives both globe + 2D; greyed when `currentStrength = 0`; enabled at default; 0 console errors.
- [ ] **Step 6: Commit**

```bash
git add hooks/useWorldEngine.ts components/Controls.tsx components/WorldViewer.tsx components/Map2D.tsx
git commit -m "feat(F2): showCurrents overlay toggle + state wiring"
```

---

## Final verification (before offering the branch for merge)

- [ ] Full gate: `npm run typecheck` (0), `npm run lint` (≤30), `npm test` (green — re-run any suspected M1 flake in isolation), `npm run build` (OK; note worker chunk size — the new code is render-side, so the worker chunk should NOT grow beyond the D2 baseline ~86.8KB).
- [ ] In-browser end-to-end (Playwright, 0 console errors): globe currents + graticule occlude correctly and rotate; 2D currents draw (Mercator + Dymaxion); toggle greyed at `currentStrength=0`; save→load drops `currents` and re-renders.
- [ ] Update `HANDOFF.md` (Session 16 entry) and flip ROADMAP F2 → 🟡 PARTIAL with the shipped-vs-queued split.

---

## Self-review notes (author)

- **Spec coverage:** §2.1 persistence+traps → Task 1; §2.2 horizon test → Task 2; §3.2 ScreenOverlay → Task 3; §3.3 tenants → Tasks 4–5; §3.4 Map2D → Task 6; §3.5 state/UI → Task 7. All covered.
- **Type consistency:** `OceanCurrentData`/`WorldData.currents` (Task 1) consumed identically in Tasks 3/4/6; `OverlayTenant`/`ProjectedCells`/`isVisible` names consistent across Tasks 2–5.
- **Testing altitude:** real TDD on the two pure units (persistence determinism, horizon math); render tasks use the repo's manual-in-browser convention with explicit verification steps — not a placeholder, the deliberate house rule.
