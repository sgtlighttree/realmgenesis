# D2 Ocean Currents Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a steady-state, deterministic ocean-current field (wind-driven, Coriolis-deflected, continent-blocked) whose sea-surface-temperature anomaly moderates coastal temperature and modulates evaporation into the existing 8-pass moisture transport.

**Architecture:** New pure worker-safe module `utils/currents.ts` runs inside the `worldGen.ts` climate block (after wind vectors, before moisture + temperature). Two fixed-pass relaxation solvers (velocity, then heat advection) — identical control-flow shape to the existing 8-pass moisture loop, so byte-identical reruns hold with no RNG. A `currentStrength` param (default 1.0, default-on) gates it; `0` short-circuits the whole stage → byte-identical to pre-D2.

**Tech Stack:** TypeScript (strict), Vitest, React (Controls slider). No new deps. Runs in the generation Web Worker.

**Spec:** `docs/superpowers/specs/2026-08-15-d2-ocean-currents-design.md` (read it alongside this plan).

## Global Constraints

- **Relative imports only** (`@/` alias is intentionally unused).
- **`utils/currents.ts` imports `types` + local math ONLY** — never `three`/`d3`. Verify the worker chunk size is unchanged after wiring (currents must not leak THREE/d3 into the worker bundle).
- **Determinism is the load-bearing test invariant.** Both solvers: `cells.forEach` + preallocated `Float32Array` scratch, fixed pass counts, fixed accumulation order, **no RNG**. Same construction as the 8-pass moisture loop.
- **`currentStrength === 0` is an early return of the whole stage**, never a multiply-by-zero. At 0 the moisture seed stays the literal `1.0` and the temperature loop adds nothing → byte-identical to pre-D2.
- **Never loosen an existing assertion to accommodate a re-baseline** — replace the seed/constant, per S9/S12/S14 precedent.
- **Milestone-D determinism break is authorized by Matt.** `currentStrength` is default-on.
- Gates that must pass before merge: `npm run typecheck` (0), `npm run lint` (0 errors / ≤29 warnings), `npm test` (full suite; honor the S10/S11 M1 parallel-load paramLiveness flake — re-run a failing paramLiveness in isolation before believing it), `npm run build`.

---

## File Structure

- **Create `utils/currents.ts`** — the two pure solvers + coupling constants. One responsibility: turn wind + geometry into a current velocity field and an SST anomaly.
- **Create `tests/currents.test.ts`** — pure-function unit tests.
- **Modify `types.ts`** — add `currentStrength` to `WorldParams`.
- **Modify `utils/defaultParams.ts`** — `currentStrength: 1.0` (single source; `makeParams` inherits it).
- **Modify `utils/worldGen.ts`** — call the solvers in the climate block; add temperature + moisture coupling with the 0-short-circuit.
- **Modify `utils/export.ts`** — `validateWorldParams` bound + `withParamDefaults` back-compat.
- **Modify `components/Controls.tsx`** — Climate-tab slider + regen dep list.
- **Modify `services/gemini.ts`** — lore prompt line.
- **Modify `tests/paramLiveness.test.ts`** — `currentStrength` perturbation.
- **Modify docs** — `ROADMAP.md`, `HANDOFF.md`, `docs/ENGINEERING-NOTES.md`, `docs/params-reference.md`.

---

## Task 1: Pure current solvers — `utils/currents.ts`

**Files:**
- Create: `utils/currents.ts`
- Test: `tests/currents.test.ts`

**Interfaces:**
- Consumes: `Cell` from `../types`, existing wind vectors as `Vec3[]` (`{x,y,z}`), `annualMeanLatTemp` from `./seasons`.
- Produces:
  - `export interface OceanCurrentField { vx: Float32Array; vy: Float32Array; vz: Float32Array; }`
  - `export function computeOceanCurrents(cells: Cell[], windVectors: { x: number; y: number; z: number }[], seaLevel: number, currentStrength: number): OceanCurrentField`
  - `export function computeSstAnomaly(cells: Cell[], field: OceanCurrentField, params: WorldParams, seaLevel: number): Float32Array`
  - Exported coupling constants `export const COAST_K = 0.6;` and `export const EVAP_K = 0.02;` (consumed by worldGen; tuned in Task 5).

- [ ] **Step 1: Write the failing tests**

```ts
// tests/currents.test.ts
import { describe, it, expect } from 'vitest';
import { generateWorld } from '../utils/worldGen';
import { computeOceanCurrents, computeSstAnomaly, OceanCurrentField } from '../utils/currents';
import { makeParams } from './helpers';
import type { Cell } from '../types';

// Build a lightweight world once for structural fixtures.
const buildCells = async (): Promise<Cell[]> =>
  (await generateWorld(makeParams({ seed: 'currents_test' }))).cells;

// A wind field that pushes every cell 'eastward' (tangent), for a deterministic probe.
const eastwardWind = (cells: Cell[]) =>
  cells.map(c => {
    const len = Math.hypot(c.center.x, c.center.z) || 1;
    return { x: -c.center.z / len, y: 0, z: c.center.x / len };
  });

describe('ocean currents', () => {
  it('is deterministic: same inputs -> identical velocity and anomaly', async () => {
    const cells = await buildCells();
    const wind = eastwardWind(cells);
    const p = makeParams({ seed: 'currents_test' });
    const a = computeOceanCurrents(cells, wind, p.seaLevel, 1.0);
    const b = computeOceanCurrents(cells, wind, p.seaLevel, 1.0);
    expect(Array.from(a.vx)).toEqual(Array.from(b.vx));
    expect(Array.from(a.vy)).toEqual(Array.from(b.vy));
    expect(Array.from(a.vz)).toEqual(Array.from(b.vz));
    const sa = computeSstAnomaly(cells, a, p, p.seaLevel);
    const sb = computeSstAnomaly(cells, b, p, p.seaLevel);
    expect(Array.from(sa)).toEqual(Array.from(sb));
  });

  it('leaves land cells with zero velocity', async () => {
    const cells = await buildCells();
    const p = makeParams({ seed: 'currents_test' });
    const field = computeOceanCurrents(cells, eastwardWind(cells), p.seaLevel, 1.0);
    cells.forEach((c, i) => {
      if (c.height >= p.seaLevel) {
        expect(field.vx[i]).toBe(0);
        expect(field.vy[i]).toBe(0);
        expect(field.vz[i]).toBe(0);
      }
    });
  });

  it('boundary tangency: coastal ocean flow has no strong into-land component', async () => {
    const cells = await buildCells();
    const p = makeParams({ seed: 'currents_test' });
    const field = computeOceanCurrents(cells, eastwardWind(cells), p.seaLevel, 1.0);
    let checked = 0;
    cells.forEach((c, i) => {
      if (c.height >= p.seaLevel) return;
      for (const n of c.neighbors) {
        if (cells[n].height < p.seaLevel) continue; // only land neighbours
        const dx = cells[n].center.x - c.center.x;
        const dy = cells[n].center.y - c.center.y;
        const dz = cells[n].center.z - c.center.z;
        const len = Math.hypot(dx, dy, dz) || 1;
        const into = (field.vx[i] * dx + field.vy[i] * dy + field.vz[i] * dz) / len;
        const speed = Math.hypot(field.vx[i], field.vy[i], field.vz[i]);
        // into-land component must be small relative to speed (tangency enforced)
        if (speed > 1e-6) expect(into).toBeLessThan(speed * 0.5 + 1e-6);
        checked++;
      }
    });
    expect(checked).toBeGreaterThan(0); // fixture actually has coastlines
  });

  it('poleward current warms high latitudes: SST anomaly correlates with poleward flow', async () => {
    const cells = await buildCells();
    const p = makeParams({ seed: 'currents_test' });
    const field = computeOceanCurrents(cells, eastwardWind(cells), p.seaLevel, 1.0);
    const anomaly = computeSstAnomaly(cells, field, p, p.seaLevel);
    // Sign check: at a high-latitude ocean cell whose current flows poleward
    // (|y| increasing), anomaly should be >= 0; equatorward flow <= 0. Aggregate.
    let polewardWarm = 0, polewardCold = 0;
    cells.forEach((c, i) => {
      if (c.height >= p.seaLevel || Math.abs(c.center.y) < 0.3) return;
      const poleward = Math.sign(c.center.y) * field.vy[i]; // >0 means flowing toward its pole
      if (poleward > 1e-4) { anomaly[i] >= 0 ? polewardWarm++ : polewardCold++; }
    });
    expect(polewardWarm).toBeGreaterThanOrEqual(polewardCold);
  });
});
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `npx vitest run tests/currents.test.ts`
Expected: FAIL — `computeOceanCurrents` / `computeSstAnomaly` not exported.

- [ ] **Step 3: Implement `utils/currents.ts`**

```ts
import { Cell, WorldParams } from '../types';
import { annualMeanLatTemp } from './seasons';

// D2 ocean currents — a fixed-pass, fixed-order relaxation (no Poisson solve,
// no RNG) so byte-identical reruns hold. See docs/superpowers/specs/2026-08-15-d2-ocean-currents-design.md.

export interface OceanCurrentField {
  vx: Float32Array;
  vy: Float32Array;
  vz: Float32Array;
}

// Tuned in the render-verify step (Task 5). Keep as named constants.
const DRAG = 0.6;        // wind-stress -> initial current speed
const CORIOLIS_K = 0.9;  // deflection magnitude scale
const MIX = 0.3;         // advective neighbour blend per pass
const N_VEL = 12;        // velocity relaxation passes
const MAX_SPEED = 1.0;   // speed ceiling
const N_HEAT = 12;       // heat-advection passes
const HEAT_MIX = 0.5;    // upstream blend per pass

// Coupling constants consumed by worldGen (centralized here for tuning).
export const COAST_K = 0.6; // land coastal-moderation weight
export const EVAP_K = 0.02; // evaporation boost per +°C of warm anomaly

const projectTangent = (vx: number, vy: number, vz: number, nx: number, ny: number, nz: number) => {
  const d = vx * nx + vy * ny + vz * nz;
  return { x: vx - d * nx, y: vy - d * ny, z: vz - d * nz };
};

export function computeOceanCurrents(
  cells: Cell[],
  windVectors: { x: number; y: number; z: number }[],
  seaLevel: number,
  currentStrength: number,
): OceanCurrentField {
  const n = cells.length;
  const vx = new Float32Array(n);
  const vy = new Float32Array(n);
  const vz = new Float32Array(n);
  const isOcean = (i: number) => cells[i].height < seaLevel;

  // Init from wind stress (tangent), ocean only.
  for (let i = 0; i < n; i++) {
    if (!isOcean(i)) continue;
    const c = cells[i].center;
    const t = projectTangent(windVectors[i].x, windVectors[i].y, windVectors[i].z, c.x, c.y, c.z);
    vx[i] = t.x * DRAG; vy[i] = t.y * DRAG; vz[i] = t.z * DRAG;
  }

  const sx = new Float32Array(n), sy = new Float32Array(n), sz = new Float32Array(n);
  for (let pass = 0; pass < N_VEL; pass++) {
    for (let i = 0; i < n; i++) {
      if (!isOcean(i)) { sx[i] = 0; sy[i] = 0; sz[i] = 0; continue; }
      const c = cells[i].center;
      // (a) Coriolis: rotate about local normal by theta ~ currentStrength*sin(lat).
      // sin(lat) == c.y (unit sphere). Rodrigues about normal k=c: for a tangent v,
      // v_rot = v*cosT + (k x v)*sinT (v is tangent so (k·v)k term vanishes).
      const lat = c.y; // = sin(latitude)
      const theta = CORIOLIS_K * currentStrength * lat;
      const cosT = Math.cos(theta), sinT = Math.sin(theta);
      let ax = vx[i], ay = vy[i], az = vz[i];
      const kx = c.x, ky = c.y, kz = c.z; // normal
      const crx = ky * az - kz * ay;
      const cry = kz * ax - kx * az;
      const crz = kx * ay - ky * ax;
      ax = ax * cosT + crx * sinT;
      ay = ay * cosT + cry * sinT;
      az = az * cosT + crz * sinT;
      // (b) advective smoothing toward mean ocean-neighbour velocity.
      let mx = 0, my = 0, mz = 0, cnt = 0;
      for (const nb of cells[i].neighbors) {
        if (!isOcean(nb)) continue;
        mx += vx[nb]; my += vy[nb]; mz += vz[nb]; cnt++;
      }
      if (cnt > 0) {
        ax = ax * (1 - MIX) + (mx / cnt) * MIX;
        ay = ay * (1 - MIX) + (my / cnt) * MIX;
        az = az * (1 - MIX) + (mz / cnt) * MIX;
      }
      // (c) boundary tangency: remove into-land component for each land neighbour.
      for (const nb of cells[i].neighbors) {
        if (isOcean(nb)) continue;
        let dx = cells[nb].center.x - c.x, dy = cells[nb].center.y - c.y, dz = cells[nb].center.z - c.z;
        // project edge dir tangent then normalize
        const dp = dx * c.x + dy * c.y + dz * c.z;
        dx -= dp * c.x; dy -= dp * c.y; dz -= dp * c.z;
        const dl = Math.hypot(dx, dy, dz);
        if (dl < 1e-9) continue;
        dx /= dl; dy /= dl; dz /= dl;
        const into = ax * dx + ay * dy + az * dz;
        if (into > 0) { ax -= into * dx; ay -= into * dy; az -= into * dz; }
      }
      // reproject to tangent + clamp speed
      const t = projectTangent(ax, ay, az, c.x, c.y, c.z);
      const sp = Math.hypot(t.x, t.y, t.z);
      const scale = sp > MAX_SPEED ? MAX_SPEED / sp : 1;
      sx[i] = t.x * scale; sy[i] = t.y * scale; sz[i] = t.z * scale;
    }
    vx.set(sx); vy.set(sy); vz.set(sz);
  }
  return { vx, vy, vz };
}

export function computeSstAnomaly(
  cells: Cell[],
  field: OceanCurrentField,
  params: WorldParams,
  seaLevel: number,
): Float32Array {
  const n = cells.length;
  const isOcean = (i: number) => cells[i].height < seaLevel;
  const base = new Float32Array(n);
  const T = new Float32Array(n);
  for (let i = 0; i < n; i++) {
    if (!isOcean(i)) continue;
    const phi = Math.asin(Math.max(-1, Math.min(1, cells[i].center.y)));
    base[i] = annualMeanLatTemp(phi, params);
    T[i] = base[i];
  }
  const scratch = new Float32Array(n);
  for (let pass = 0; pass < N_HEAT; pass++) {
    for (let i = 0; i < n; i++) {
      if (!isOcean(i)) { scratch[i] = 0; continue; }
      const c = cells[i].center;
      let incoming = 0, cnt = 0;
      for (const nb of cells[i].neighbors) {
        if (!isOcean(nb)) continue;
        // edgeDir from nb toward i (tangent, unnormalized ok for sign)
        const dx = c.x - cells[nb].center.x;
        const dy = c.y - cells[nb].center.y;
        const dz = c.z - cells[nb].center.z;
        const dot = dx * field.vx[nb] + dy * field.vy[nb] + dz * field.vz[nb];
        if (dot > 0) { incoming += T[nb]; cnt++; } // nb flows toward i => upstream
      }
      if (cnt === 0) scratch[i] = T[i] * 0.95 + base[i] * 0.05;
      else scratch[i] = T[i] * (1 - HEAT_MIX) + (incoming / cnt) * HEAT_MIX;
    }
    T.set(scratch);
  }
  const anomaly = new Float32Array(n);
  for (let i = 0; i < n; i++) if (isOcean(i)) anomaly[i] = T[i] - base[i];
  return anomaly;
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `npx vitest run tests/currents.test.ts`
Expected: PASS (4 tests). If the poleward-sign test is noisy on this seed, widen to `toBeGreaterThanOrEqual` aggregate (already is) — do not weaken determinism/tangency tests.

- [ ] **Step 5: Typecheck + commit**

```bash
npm run typecheck
git add utils/currents.ts tests/currents.test.ts
git commit -m "feat(D2): pure ocean-current + SST-anomaly solvers"
```

---

## Task 2: Add `currentStrength` param (type + default + validation)

**Files:**
- Modify: `types.ts` (WorldParams, Climate section near `starClass`)
- Modify: `utils/defaultParams.ts`
- Modify: `utils/export.ts:18` (`withParamDefaults`) and `:480` region (`validateWorldParams` bounds)

**Interfaces:**
- Produces: `WorldParams.currentStrength: number` available to worldGen (Task 3) and Controls (Task 4).

- [ ] **Step 1: Add to the type.** In `types.ts`, in the `// Climate` block after the `starClass` line:

```ts
  currentStrength: number; // D2: ocean-current intensity (0-2). 0 = disabled/byte-identical; scales Coriolis, coastal temperature moderation, and warm-current evaporation.
```

- [ ] **Step 2: Add the default.** In `utils/defaultParams.ts`, after `starClass: 'G',`:

```ts
  currentStrength: 1.0,
```

- [ ] **Step 3: Back-compat + validation.** In `utils/export.ts` `withParamDefaults` (after the `starClass:` line ~25):

```ts
  // D2: pre-D2 saves lack currentStrength → default to 1.0 (default-on).
  currentStrength: typeof params.currentStrength === 'number' && isFinite(params.currentStrength) ? params.currentStrength : 1.0,
```

And in the `validateWorldParams` bounds object (near `season: [0, 1],`):

```ts
        currentStrength: [0, 2],
```

- [ ] **Step 4: Verify typecheck passes (both param lists complete via single-source).**

Run: `npm run typecheck`
Expected: 0 errors. (`makeParams` inherits `currentStrength` free via the S15 spread — no test file edit needed here.)

- [ ] **Step 5: Commit**

```bash
git add types.ts utils/defaultParams.ts utils/export.ts
git commit -m "feat(D2): currentStrength param — type, default, validation, back-compat"
```

---

## Task 3: Wire currents into the climate block (temperature + moisture coupling)

**Files:**
- Modify: `utils/worldGen.ts:514-589` (the climate block: wind vectors → moisture 8-pass → temperature loop)
- Test: reuse `tests/worldGen.test.ts` determinism + add a no-op assertion (below)

**Interfaces:**
- Consumes: `computeOceanCurrents`, `computeSstAnomaly`, `COAST_K`, `EVAP_K` from `./currents`.
- Produces: `cell.temperature` and `cell.moisture` now reflect currents when `currentStrength > 0`.

- [ ] **Step 1: Write the failing no-op test.** Add to `tests/worldGen.test.ts`:

```ts
it('D2: currentStrength=0 is byte-identical to the pre-D2 climate (no-op escape hatch)', async () => {
  const off = await generateWorld(makeParams({ currentStrength: 0 }));
  const on = await generateWorld(makeParams({ currentStrength: 1.0 }));
  // Off must match a run that never touched currents; On must differ (coasts move).
  const sigOff = off.cells.map(c => `${c.temperature.toFixed(6)}|${c.moisture.toFixed(6)}`).join(';');
  const off2 = await generateWorld(makeParams({ currentStrength: 0 }));
  const sigOff2 = off2.cells.map(c => `${c.temperature.toFixed(6)}|${c.moisture.toFixed(6)}`).join(';');
  const sigOn = on.cells.map(c => `${c.temperature.toFixed(6)}|${c.moisture.toFixed(6)}`).join(';');
  expect(sigOff).toBe(sigOff2);   // deterministic + stage skipped
  expect(sigOn).not.toBe(sigOff); // default-on actually changes climate
});
```

- [ ] **Step 2: Run to verify it fails.**

Run: `npx vitest run tests/worldGen.test.ts -t "D2"`
Expected: FAIL — `sigOn` currently equals `sigOff` (currents not wired yet).

- [ ] **Step 3: Add the import** at the top of `utils/worldGen.ts` (with the other `./` imports):

```ts
import { computeOceanCurrents, computeSstAnomaly, COAST_K, EVAP_K, OceanCurrentField } from './currents';
```

- [ ] **Step 4: Compute currents after `windVectors`** (right after the `const windVectors = ...` block ends, ~line 531). Insert:

```ts
  // D2: ocean currents. currentStrength === 0 short-circuits the whole stage
  // (moisture seed stays literal 1.0, temperature adds nothing) -> byte-identical to pre-D2.
  const currentStrength = params.currentStrength ?? 1.0;
  let currentField: OceanCurrentField | null = null;
  let sstAnomaly: Float32Array | null = null;
  if (currentStrength > 0) {
    currentField = computeOceanCurrents(cells, windVectors, params.seaLevel, currentStrength);
    sstAnomaly = computeSstAnomaly(cells, currentField, params, params.seaLevel);
  }
```

- [ ] **Step 5: Modify the moisture seed** (the `cells.forEach` at ~533 and the ocean branch inside the 8-pass at ~543). Change the initial seed:

```ts
  cells.forEach((c, i) => {
    if (c.height < params.seaLevel) {
      c.moisture = sstAnomaly ? 1.0 + EVAP_K * Math.max(0, sstAnomaly[i]) : 1.0;
    } else {
      c.moisture = 0.1 * params.rainfallMultiplier;
    }
  });
```

And inside the 8-pass loop, the ocean pin (`if (c.height < params.seaLevel) { newMoisture[i] = 1.0; return; }`) becomes:

```ts
      if (c.height < params.seaLevel) {
        newMoisture[i] = sstAnomaly ? 1.0 + EVAP_K * Math.max(0, sstAnomaly[i]) : 1.0;
        return;
      }
```

(When `sstAnomaly` is null the expression is the literal `1.0` — identical to today.)

- [ ] **Step 6: Modify the temperature loop** (~575). After `let temp = applyStarClass(...)` and before the elevation lapse, add the current term:

```ts
      // D2: ocean cells take their own SST anomaly; land cells take a 1-ring
      // coastal-moderation blend of adjacent ocean anomalies. Null field => no-op.
      if (sstAnomaly) {
        if (c.height < params.seaLevel) {
          temp += sstAnomaly[cells.indexOf(c)] * currentStrength;
        } else {
          let sum = 0, cnt = 0;
          for (const nb of c.neighbors) if (cells[nb].height < params.seaLevel) { sum += sstAnomaly[nb]; cnt++; }
          if (cnt > 0) temp += COAST_K * currentStrength * (sum / cnt);
        }
      }
```

**NOTE:** `cells.indexOf(c)` is O(n) — the surrounding loop is `cells.forEach(c => ...)` without an index. Change the loop signature to `cells.forEach((c, ci) => ...)` and use `sstAnomaly[ci]` instead of `cells.indexOf(c)`. Verify the forEach at ~575 exposes the index; if not, add it.

- [ ] **Step 7: Run the no-op + determinism tests.**

Run: `npx vitest run tests/worldGen.test.ts`
Expected: PASS — `sigOff === sigOff2`, `sigOn !== sigOff`, existing determinism test still green.

- [ ] **Step 8: Typecheck + commit**

```bash
npm run typecheck && npm run lint
git add utils/worldGen.ts tests/worldGen.test.ts
git commit -m "feat(D2): couple currents into temperature + moisture (0 = short-circuit no-op)"
```

---

## Task 4: UI slider + lore + param-liveness

**Files:**
- Modify: `components/Controls.tsx` (Climate tab ~944-960; regen dep list ~214)
- Modify: `services/gemini.ts:49` region
- Modify: `tests/paramLiveness.test.ts:42` region (`TERRAIN_PERTURBATIONS`)

**Interfaces:**
- Consumes: `params.currentStrength` (Task 2), live engine effect (Task 3).

- [ ] **Step 1: Add the param-liveness case.** In `tests/paramLiveness.test.ts` `TERRAIN_PERTURBATIONS`, after the `starClass` line:

```ts
  currentStrength: { currentStrength: 0 }, // D2: disabling currents changes the climate signature (default is 1.0)
```

- [ ] **Step 2: Run it — expect PASS** (Task 3 made currents live, so `currentStrength: 0` differs from the default-on baseline).

Run: `npx vitest run tests/paramLiveness.test.ts -t "terrain"`
Expected: PASS. If it flakes with a 120s timeout under load (S10/S11), re-run in isolation before believing failure.

- [ ] **Step 3: Add the Controls slider.** In `components/Controls.tsx`, in the Climate tab after the star-class `Select` (~960), mirror the season slider's structure:

```tsx
                <div>
                  <label>Ocean Currents (Coastal Climate)</label>
                  <input
                    type="range" min={0} max={2} step={0.05}
                    value={params.currentStrength ?? 1.0}
                    onChange={(e) => { handleChange('currentStrength', parseFloat(e.target.value)); }}
                  />
                  <span>{(params.currentStrength ?? 1.0) === 0 ? 'Off' : `${(params.currentStrength ?? 1.0).toFixed(2)}×`}</span>
                </div>
```

(Match the exact class names / wrapper markup of the adjacent sliders — copy the season slider's wrapper.)

- [ ] **Step 4: Add to the auto-update regen dep list.** In `components/Controls.tsx` ~214 (the array containing `params.starClass,`), add:

```tsx
      params.currentStrength,
```

- [ ] **Step 5: Surface in lore.** In `services/gemini.ts` ~49, after the star-class line:

```ts
    - Ocean currents ${(world.params.currentStrength ?? 1.0) === 0 ? 'are absent (no coastal moderation)' : 'moderate coastal climate'}.
```

- [ ] **Step 6: Gates + commit**

```bash
npm run typecheck && npm run lint && npx vitest run tests/paramLiveness.test.ts
git add components/Controls.tsx services/gemini.ts tests/paramLiveness.test.ts
git commit -m "feat(D2): currentStrength slider, lore surfacing, param-liveness case"
```

---

## Task 5: Render-verify, tune constants, re-verify the D1 invariant

**Files:**
- Modify (tuning only): `utils/currents.ts` constants (`DRAG`/`CORIOLIS_K`/`MIX`/`COAST_K`/`EVAP_K`) if the render calls for it.

This task is **self / browser work — not delegatable.** No new automated test; it produces evidence recorded in HANDOFF.

- [ ] **Step 1: Start the dev server** (reuse if already running on :3000; only start one you own).

Run: `npm run dev` (background)

- [ ] **Step 2: Render `realmgenesis` at default** (`currentStrength = 1.0`), satellite + temperature + biome views, 2D Mercator. Confirm: gyre-like circulation reads in the current-affected temperature (warm tongues reaching poleward on west sides of oceans / cool on east coasts), no NaN/black cells, 0 console errors. Screenshot to the session scratchpad. Tune constants if gyres are invisible (raise `CORIOLIS_K` / boundary weight) or overwhelming (lower `DRAG`/`COAST_K`); re-run Tasks 1/3 tests after any constant change.

- [ ] **Step 3: D1 escape-hatch proof.** Set `currentStrength = 0`, 5k cells, Mercator, biome view, season 0.5. Confirm the pixel checksum is **909197 exactly** (the S14 D1 baseline). This proves the escape hatch end-to-end. If it is NOT 909197, the short-circuit is leaking — STOP and fix Task 3 before proceeding.

- [ ] **Step 4: D1 new-baseline proof.** Set `currentStrength = 1.0`. Capture the new neutral (season 0.5) biome-view checksum `B`. Drive the season slider 0.5 → 0.15 → 0.5 and confirm it returns to `B` **exactly** (neutral == canonical + return-to-exact against the new baseline).

- [ ] **Step 5: Record evidence + commit any tuning.**

```bash
git add utils/currents.ts   # only if constants changed
git commit -m "tune(D2): current solver constants from render-verify"   # skip if no change
```

Record in HANDOFF (Task 7): the two checksums, gyre observations, and final constants.

---

## Task 6: Re-baseline the climate-sensitive seed fixtures  [DELEGATABLE]

**Files:**
- Modify: `tests/lakes.test.ts` (`SALT_SEED` / `FRESH_SEED`), `tests/routes.test.ts` (`SEA_SEED`), and any pinned expectations in `tests/biomes.test.ts` / `tests/features.test.ts` / `tests/religions.test.ts` that shift.

**Delegation brief (tier `sonnet-low`, run this task in a subagent):**

> Temperature AND moisture now shift for every coastal world (D2 currents, default-on). Re-baseline the pinned-seed climate fixtures WITHOUT loosening any assertion — replace the seed constant / expected count only (S9/S12/S14 precedent).
>
> **Method (mandatory, from S14):** any seed-scanning script must write results to a JSON file via `fs.writeFileSync` and read that back. **Do NOT trust piped vitest `console.log`** — S14 documents a piped scan emitting zero output. **Run every test/scan synchronously in the foreground, per-file, never backgrounded** (S12: a subagent cannot be resumed by a background notification).
>
> Steps: (1) `npx vitest run tests/lakes.test.ts` — if the salt-endorheic / fresh assertions fail, scan seeds `s1..s200` for a single-salt-endorheic-lake seed and a fresh-lake seed exactly as the prior fixtures required (see the file's existing structure), update `SALT_SEED` / `FRESH_SEED` + cell counts. (2) Same for `tests/routes.test.ts` `SEA_SEED`. (3) Run `tests/biomes.test.ts`, `tests/features.test.ts`, `tests/religions.test.ts`; re-baseline any pinned seed/count that shifted. (4) Report the full `npx vitest run` results per file and the exact constant changes made.

- [ ] **Step 1: Dispatch the subagent** with the brief above.
- [ ] **Step 2: Review the diff** — confirm only seed constants / expected counts changed, no assertion was weakened.
- [ ] **Step 3: Run the full suite yourself** to confirm green (honor the paramLiveness isolation rule).

```bash
npm test
```

- [ ] **Step 4: Commit**

```bash
git add tests/lakes.test.ts tests/routes.test.ts tests/biomes.test.ts tests/features.test.ts tests/religions.test.ts
git commit -m "test(D2): re-baseline climate-sensitive seed fixtures for ocean currents"
```

---

## Task 7: Docs — ROADMAP, ENGINEERING-NOTES seam, params-reference, HANDOFF

**Files:**
- Modify: `ROADMAP.md` (D2 → ✅ DONE or 🟡 PARTIAL with note), `docs/ENGINEERING-NOTES.md` (monsoon composability seam), `docs/params-reference.md` (`currentStrength` row), `HANDOFF.md` (Session 15 continuation entry with the two D1 checksums + final constants).

- [ ] **Step 1: ROADMAP D2** — flip the status line + add a one-line shipped note (fixed-pass relaxation, temp+moisture, `currentStrength` default-on, monsoon seam deferred).
- [ ] **Step 2: ENGINEERING-NOTES** — add a "D2 → deferred seasonal-moisture seam" entry: the current-driven evaporation term is the ANNUAL ocean-moisture baseline; the future monsoon solver layers a per-season overlay ON TOP of it and must not overwrite the ocean seed (introduces no per-season field, so D1/D3's free O(n) biome-at-season recompute is preserved).
- [ ] **Step 3: params-reference** — add the `currentStrength` row (range 0–2, default 1.0, default-on, 0 = byte-identical).
- [ ] **Step 4: HANDOFF** — extend the Session 15 entry: D2 shipped, the two D1 checksums (909197 at 0; new baseline B at default), final tuned constants, re-baselined fixtures, NOT pushed.
- [ ] **Step 5: Commit**

```bash
git add ROADMAP.md docs/ENGINEERING-NOTES.md docs/params-reference.md HANDOFF.md
git commit -m "docs(D2): ROADMAP, ENGINEERING-NOTES monsoon seam, params-reference, HANDOFF"
```

---

## Self-Review (completed by author)

**Spec coverage:** §3 placement → Task 3; §4 solver → Task 1; §5 anomaly → Task 1; §6 temp coupling → Task 3; §7 moisture coupling + seam → Task 3 (code) + Task 7 (seam doc); §8 determinism/escape-hatch/D1 re-verify → Task 3 (no-op test) + Task 5 (browser); §9 param surface → Tasks 2+4; §10 testing → Tasks 1,3,4,6; §11 delegation → Task 6; §12 risks → covered across tasks. No gaps.

**Placeholder scan:** constants carry concrete starting values (tuned in Task 5, explicitly flagged — not TBD). All code steps show code. No "similar to Task N."

**Type consistency:** `OceanCurrentField {vx,vy,vz}`, `computeOceanCurrents(cells, windVectors, seaLevel, currentStrength)`, `computeSstAnomaly(cells, field, params, seaLevel)`, `COAST_K`, `EVAP_K`, `currentStrength` used identically across Tasks 1/3/4. The Task-3 index fix (`forEach((c, ci))`) is called out explicitly so `sstAnomaly[ci]` is valid.
