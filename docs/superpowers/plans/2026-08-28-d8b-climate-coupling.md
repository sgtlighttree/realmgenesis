# D8b — Datum Climate Coupling (`physicalClimate`) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Feed the D8a metre datum back into generation so three tuned magic constants (lapse rate, orographic rain shadow, snow line) become grounded physics, behind a default-on `physicalClimate` flag whose off-state is byte-identical to current `main`.

**Architecture:** One new boolean `WorldParams` key gates two edit sites in `utils/worldGen.ts` — the per-cell lapse-rate term and the 8-pass orographic moisture term. Both grounded branches call `elevationMetres` (already in `utils/datum.ts`); new physical constants live beside it. `determineBiome` is untouched — the snow line moves for free through temperature. The moisture retune is folded into the orographic rewrite and tuned empirically against a measured target.

**Tech Stack:** TypeScript (strict), Vite, Vitest, React 19. Engine runs in a Web Worker; generation is pure functions over a cell graph.

**Spec:** `docs/superpowers/specs/2026-08-23-d8b-climate-coupling-design.md`

## Global Constraints

- **`physicalClimate` default is `true`.** Off reproduces current `main` (post-curve, post-grassland) byte-identical — NOT ancient main.
- **Byte-identical baseline is CURRENT `main`.** The grassland biome and the `frac^2` curve are already in `main`; the hatch reproduces the engine as it stands after S25b.
- **Escape hatch is byte-identical BY CONSTRUCTION.** Each grounded site keeps the exact old expression verbatim in its `else` branch — `elevation * 60`, and `heightDiff > 0.02 → *1.5` / `< -0.02 → *0.2`. Do not "clean up" the else branch; a rewrite risks drift.
- **Datum pick is fixed: physical lapse 6.5 °C/km + datum 9000 + the `frac^2` curve, all three.** Do not reopen this — the curve dissolved the former pick-two tradeoff (ICE_CAP stays ~7%).
- **`maxElevationM` becomes a GENERATION param (Matt's call, 2026-08-28).** Moving it requires a regenerate, like Mountain Height. Lapse/orographic read `params.maxElevationM`. Its live display-only sync in `hooks/useWorldEngine.ts` is REMOVED for `maxElevationM` (season stays). This closes the determinism hole where a live datum would not match the world it baked. It enters `TERRAIN_PERTURBATIONS`.
- **New physical constants live in `utils/datum.ts`** (single source beside the datum they depend on), imported by `worldGen.ts`. Relative imports only (`./datum`).
- **Seed breakage is authorized** (ROADMAP:34) and Matt accepts the civ-layout movement default-ON causes. This is not a regression.
- **`determineBiome` is NOT modified.** The snow-line change is emergent through temperature.
- Code style: 2-space indent, semicolons, single quotes, trailing commas. Verify with `npm run lint` (zero errors, ≤30 warnings) and `npm run typecheck` (strict, zero errors).

---

## File Structure

- `utils/datum.ts` — **modify.** Add physical constants: `LAPSE_RATE_C_PER_KM`, `OROG_WINDWARD_PER_KM`, `OROG_WINDWARD_CAP`, `OROG_LEEWARD_PER_KM`, `OROG_LEEWARD_FLOOR`.
- `utils/worldGen.ts` — **modify.** Import `elevationMetres` + the constants; gate the lapse term (~L621-622) and the orographic term (~L582-584) on `params.physicalClimate`.
- `types.ts` — **modify.** Add `physicalClimate: boolean` to `WorldParams` (near `season`, ~L134).
- `utils/defaultParams.ts` — **modify.** `physicalClimate: true`.
- `utils/export.ts` — **modify.** `withParamDefaults` defaults missing key to `true`; `validateWorldParams` gets a boolean type-check line (NOT a `numericBounds` entry).
- `hooks/useWorldEngine.ts` — **modify.** Remove `maxElevationM` from the live display-only sync effect (~L138-153); leave `season` there. It is now a generation param (regenerate to apply).
- `components/Controls.tsx` — **modify.** A toggle in the Climate tab with warning copy.
- `tests/datum.test.ts` — **modify.** Assert the new constants exist with expected values.
- `tests/paramLiveness.test.ts` — **modify.** Move `maxElevationM` into `TERRAIN_PERTURBATIONS`; add `physicalClimate: false`; drop `maxElevationM` from the display-only allowlist comment.
- `tests/climateCoupling.test.ts` — **create.** Unit tests: hatch-off keeps old lapse/orographic math; hatch-on changes temperature at altitude.
- `scripts/captureBaseline.mjs` — **used, not modified** for the manual byte-identical gate (see Task 7). No committed fixture — the script's own header forbids it (V8 `Math.*` drift).

---

## Task 1: Add physical constants to `utils/datum.ts`

**Files:**
- Modify: `utils/datum.ts` (after `HYPSOMETRIC_EXPONENT`, ~L41)
- Test: `tests/datum.test.ts`

**Interfaces:**
- Consumes: nothing.
- Produces: `LAPSE_RATE_C_PER_KM: number = 6.5`, `OROG_WINDWARD_PER_KM: number = 0.5`, `OROG_WINDWARD_CAP: number = 2.5`, `OROG_LEEWARD_PER_KM: number = 0.3`, `OROG_LEEWARD_FLOOR: number = 0.5`. These are STARTING values for the orographic pair (Task 5 tunes them); lapse is fixed at 6.5.

- [ ] **Step 1: Write the failing test**

Add to `tests/datum.test.ts`:

```ts
import {
  LAPSE_RATE_C_PER_KM,
  OROG_WINDWARD_PER_KM,
  OROG_WINDWARD_CAP,
  OROG_LEEWARD_PER_KM,
  OROG_LEEWARD_FLOOR,
} from '../utils/datum';

describe('D8b physical constants', () => {
  it('uses the standard environmental lapse rate', () => {
    expect(LAPSE_RATE_C_PER_KM).toBe(6.5);
  });
  it('exposes bounded orographic constants', () => {
    expect(OROG_WINDWARD_PER_KM).toBeGreaterThan(0);
    expect(OROG_WINDWARD_CAP).toBeGreaterThan(1);
    expect(OROG_LEEWARD_PER_KM).toBeGreaterThan(0);
    expect(OROG_LEEWARD_FLOOR).toBeGreaterThan(0);
    expect(OROG_LEEWARD_FLOOR).toBeLessThan(1);
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npx vitest run tests/datum.test.ts -t "D8b physical constants"`
Expected: FAIL — the constants are not exported.

- [ ] **Step 3: Add the constants**

In `utils/datum.ts`, after the `HYPSOMETRIC_EXPONENT` export (~L41):

```ts
/**
 * D8b — physical climate coupling constants (gated by `physicalClimate`).
 *
 * `LAPSE_RATE_C_PER_KM` is the standard environmental lapse rate: air cools
 * 6.5 °C for every kilometre of altitude. It replaces the invented `* 60`
 * multiplier on normalized height in worldGen's temperature finalize.
 *
 * The OROG_* constants scale the orographic (rain-shadow) moisture term by the
 * real barrier metres crossed. Windward ascent boosts rain (capped so coasts do
 * not saturate); leeward descent dries it (floored so a single edge cannot lose
 * everything). STARTING values — tuned empirically against the land-moisture
 * target in the D8b spec §6 (scripts/queryWorld.mjs climate/biomes).
 */
export const LAPSE_RATE_C_PER_KM = 6.5;
export const OROG_WINDWARD_PER_KM = 0.5;
export const OROG_WINDWARD_CAP = 2.5;
export const OROG_LEEWARD_PER_KM = 0.3;
export const OROG_LEEWARD_FLOOR = 0.5;
```

- [ ] **Step 4: Run test to verify it passes**

Run: `npx vitest run tests/datum.test.ts`
Expected: PASS (all datum tests, new block included).

- [ ] **Step 5: Commit**

```bash
git add utils/datum.ts tests/datum.test.ts
git commit -m "feat(D8b): add physical lapse + orographic constants to datum"
```

---

## Task 2: Add the `physicalClimate` param end-to-end (type, default, validator, defaults)

**Files:**
- Modify: `types.ts` (~L134, near `season`)
- Modify: `utils/defaultParams.ts` (~L42)
- Modify: `utils/export.ts` — `withParamDefaults` (~L19-31), `validateWorldParams` (~L509-515)
- Test: `tests/export.test.ts` (or the existing import-validation test file)

**Interfaces:**
- Consumes: nothing.
- Produces: `WorldParams.physicalClimate: boolean`, defaulting to `true` everywhere. `worldGen.ts` (Tasks 3-4) reads `params.physicalClimate`.

**Verified passthrough (no extra task needed):**
- `tests/helpers.ts` `makeParams` spreads `DEFAULT_PARAMS`, so adding the key there (Step 5) flows into every test with no helper edit.
- The worker passes the whole params object: `worldGenClient.ts:104` posts `{ type: 'generate', params }` and `worldGen.worker.ts:25` calls `generateWorld(e.data.params, ...)`. `physicalClimate` reaches the browser engine — no worker change required.

- [ ] **Step 1: Find the import-validation test file**

Run: `grep -rln "validateWorldParams\|withParamDefaults" tests/`
Read the matching file. Add the new tests there. If none exists, create `tests/paramDefaults.test.ts` importing from `../utils/export`.

- [ ] **Step 2: Write the failing tests**

```ts
import { validateWorldParams, withParamDefaults } from '../utils/export';
import { makeParams } from './helpers';

describe('physicalClimate param', () => {
  it('defaults a missing physicalClimate to true (old saves get grounded climate)', () => {
    const p = makeParams();
    delete (p as Record<string, unknown>).physicalClimate;
    expect(withParamDefaults(p).physicalClimate).toBe(true);
  });
  it('preserves an explicit false', () => {
    const p = makeParams({ physicalClimate: false });
    expect(withParamDefaults(p).physicalClimate).toBe(false);
  });
  it('rejects a non-boolean physicalClimate', () => {
    const bad = { ...makeParams(), physicalClimate: 'yes' } as unknown;
    expect(validateWorldParams(bad)).toBe(false);
  });
  it('accepts a boolean physicalClimate', () => {
    expect(validateWorldParams({ ...makeParams(), physicalClimate: true })).toBe(true);
  });
});
```

- [ ] **Step 3: Run test to verify it fails**

Run: `npx vitest run tests/<the file>.ts -t "physicalClimate param"`
Expected: FAIL — the type/default/validator do not exist yet.

- [ ] **Step 4: Add the type**

In `types.ts`, after the `season` line (~L134):

```ts
  physicalClimate: boolean; // D8b: grounded lapse rate + orographic rain shadow from the metre datum. Default true; false = byte-identical old formulas. Changes generated climate AND civ layout for existing seeds.
```

- [ ] **Step 5: Add the default**

In `utils/defaultParams.ts`, near `season: 0.5` (~L42):

```ts
  physicalClimate: true,
```

- [ ] **Step 6: Add the withParamDefaults fallback**

In `utils/export.ts` `withParamDefaults` (after the `maxElevationM` line, ~L30):

```ts
  // D8b: pre-D8b saves lack physicalClimate → default to grounded (true). Old
  // saves get grounded climate on load — the accepted civ-move; reproducible
  // with false.
  physicalClimate: typeof params.physicalClimate === 'boolean' ? params.physicalClimate : true,
```

- [ ] **Step 7: Add the validator boolean check**

In `utils/export.ts` `validateWorldParams`, alongside the other type-checks (after the `starClass` line, ~L515):

```ts
    if ('physicalClimate' in p && typeof p.physicalClimate !== 'boolean') return false;
```

Do NOT add `physicalClimate` to `numericBounds` — it is a boolean.

- [ ] **Step 8: Run tests to verify they pass**

Run: `npx vitest run tests/<the file>.ts -t "physicalClimate param"`
Expected: PASS. Then `npm run typecheck` — expected zero errors.

- [ ] **Step 9: Commit**

```bash
git add types.ts utils/defaultParams.ts utils/export.ts tests/<the file>.ts
git commit -m "feat(D8b): add physicalClimate param (type, default, validator, save-default)"
```

---

## Task 3: Site 1 — grounded lapse rate (gated)

**Files:**
- Modify: `utils/worldGen.ts` — import block (~L11) and the lapse term (~L621-622)
- Create: `tests/climateCoupling.test.ts`

**Interfaces:**
- Consumes: `params.physicalClimate` (Task 2), `elevationMetres` + `LAPSE_RATE_C_PER_KM` (Task 1).
- Produces: temperature that cools by `(elevationMetres/1000) * 6.5` on land when the hatch is on; the exact `elevation * 60` when off.

- [ ] **Step 1: Write the failing test**

Create `tests/climateCoupling.test.ts`:

```ts
import { describe, it, expect } from 'vitest';
import { generateWorld } from '../utils/worldGen';
import { makeParams, terrainSignature } from './helpers';

describe('D8b lapse-rate coupling', () => {
  it('hatch off vs on produces a different terrain signature (lapse is live)', async () => {
    const off = await generateWorld(makeParams({ physicalClimate: false }));
    const on = await generateWorld(makeParams({ physicalClimate: true }));
    expect(terrainSignature(on)).not.toBe(terrainSignature(off));
  }, 120000);

  it('a taller datum cools high peaks further when the hatch is on', async () => {
    // Larger maxElevationM => more metres per height unit => colder peaks =>
    // a different terrain signature. Proves lapse reads the datum, not raw height.
    const lo = await generateWorld(makeParams({ physicalClimate: true, maxElevationM: 4500 }));
    const hi = await generateWorld(makeParams({ physicalClimate: true, maxElevationM: 18000 }));
    expect(terrainSignature(hi)).not.toBe(terrainSignature(lo));
  }, 120000);
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npx vitest run tests/climateCoupling.test.ts -t "lapse-rate"`
Expected: FAIL — worldGen ignores `physicalClimate`, so both signatures match.

- [ ] **Step 3: Add the import**

In `utils/worldGen.ts`, extend the imports (~L11, after the currents import):

```ts
import { elevationMetres, LAPSE_RATE_C_PER_KM, OROG_WINDWARD_PER_KM, OROG_WINDWARD_CAP, OROG_LEEWARD_PER_KM, OROG_LEEWARD_FLOOR } from './datum';
```

(The OROG_* imports are unused until Task 4 — that is fine within this branch; if lint flags unused imports as an error mid-task, add only `elevationMetres` and `LAPSE_RATE_C_PER_KM` here and the rest in Task 4. Check `npm run lint` after Step 4.)

- [ ] **Step 4: Gate the lapse term**

Replace `utils/worldGen.ts:621-622`:

```ts
      const elevation = Math.max(0, c.height - params.seaLevel);
      temp -= elevation * 60;
```

with:

```ts
      if (params.physicalClimate) {
        const elevM = Math.max(0, elevationMetres(c.height, params.seaLevel, params.maxElevationM));
        temp -= (elevM / 1000) * LAPSE_RATE_C_PER_KM;
      } else {
        const elevation = Math.max(0, c.height - params.seaLevel);
        temp -= elevation * 60; // exact old formula — byte-identical hatch
      }
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `npx vitest run tests/climateCoupling.test.ts -t "lapse-rate"`
Expected: PASS (both cases).
Then: `npm run lint` — expected zero errors (OROG_* may warn as unused until Task 4; a warning is allowed, an error is not — if it errors, apply the Step-3 note).

- [ ] **Step 6: Commit**

```bash
git add utils/worldGen.ts tests/climateCoupling.test.ts
git commit -m "feat(D8b): ground the lapse rate on the metre datum (gated)"
```

---

## Task 4: Site 2 — grounded orographic precipitation (gated)

**Files:**
- Modify: `utils/worldGen.ts` — the orographic term (~L582-584, inside the `if (dot > 0)` windward branch)
- Test: `tests/climateCoupling.test.ts`

**Interfaces:**
- Consumes: `params.physicalClimate`, `elevationMetres`, `OROG_WINDWARD_PER_KM`, `OROG_WINDWARD_CAP`, `OROG_LEEWARD_PER_KM`, `OROG_LEEWARD_FLOOR` (Task 1).
- Produces: moisture whose rain-shadow scaling reads real barrier metres when the hatch is on; the exact `1.5 / 0.2` steps when off.

- [ ] **Step 1: Write the failing test**

Add to `tests/climateCoupling.test.ts`:

```ts
describe('D8b orographic coupling', () => {
  it('changes the land-moisture signature vs the old normalized threshold', async () => {
    // Isolate moisture: same seed, only the hatch flips. terrainSignature
    // includes moisture, so a change here proves the orographic branch is live.
    // (Lapse also moves under the hatch — this test asserts the pair is live,
    // Task 5 measures the moisture distribution specifically.)
    const off = await generateWorld(makeParams({ physicalClimate: false }));
    const on = await generateWorld(makeParams({ physicalClimate: true }));
    const sea = makeParams().seaLevel;
    const moist = (w: typeof on) =>
      w.cells.filter(c => c.height >= sea).map(c => c.moisture.toFixed(4)).join(',');
    expect(moist(on)).not.toBe(moist(off));
  }, 120000);
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npx vitest run tests/climateCoupling.test.ts -t "orographic"`
Expected: FAIL — after Task 3 the moisture is genuinely identical between hatch states. The lapse term runs AFTER the 8-pass moisture loop (`worldGen.ts:621` vs `:563-595`), so it cannot change moisture; only gating the orographic term will. The test fails cleanly until Step 3.

- [ ] **Step 3: Gate the orographic term**

Replace `utils/worldGen.ts:582-584`:

```ts
                 const heightDiff = c.height - n.height;
                 if (heightDiff > 0.02) carry *= 1.5;
                 else if (heightDiff < -0.02) carry *= 0.2;
```

with:

```ts
                 if (params.physicalClimate) {
                   // Land-side above-sea metres only (ocean neighbour → 0), so a
                   // coastal barrier is measured from sea level and we never cross
                   // the seaLevel slope kink inside the loop (spec §5, advisor item 3).
                   const cM = Math.max(0, elevationMetres(c.height, params.seaLevel, params.maxElevationM));
                   const nM = Math.max(0, elevationMetres(n.height, params.seaLevel, params.maxElevationM));
                   const dKm = (cM - nM) / 1000;
                   if (dKm > 0) carry *= Math.min(OROG_WINDWARD_CAP, 1 + OROG_WINDWARD_PER_KM * dKm);
                   else if (dKm < 0) carry *= Math.max(OROG_LEEWARD_FLOOR, 1 + OROG_LEEWARD_PER_KM * dKm);
                 } else {
                   const heightDiff = c.height - n.height;
                   if (heightDiff > 0.02) carry *= 1.5;
                   else if (heightDiff < -0.02) carry *= 0.2; // exact old formula — byte-identical hatch
                 }
```

Note `OROG_LEEWARD_PER_KM` is positive and `dKm` is negative in the leeward branch, so `1 + 0.3 * dKm` decreases toward the floor. This matches the spec §5 sign convention.

- [ ] **Step 4: Run tests to verify they pass**

Run: `npx vitest run tests/climateCoupling.test.ts`
Expected: PASS (lapse + orographic). Then `npm run lint` — expected zero errors (OROG_* now used).

- [ ] **Step 5: Commit**

```bash
git add utils/worldGen.ts tests/climateCoupling.test.ts
git commit -m "feat(D8b): ground orographic rain shadow on barrier metres (gated)"
```

---

## Task 5: Moisture retune — tune OROG constants to the measured target

This task has NO clean closed form (spec §6). Its deliverable is tuned constant values that move the land-moisture distribution toward the target, verified by measurement, not a pass/fail unit test.

**Files:**
- Modify: `utils/datum.ts` — the `OROG_*` values ONLY.
- Measure with: `scripts/queryWorld.mjs`

**The `* 0.98` land decay (`worldGen.ts:592`) is OFF-LIMITS here.** It is ungated shared code, so changing it would break the byte-identical hatch (the off-state would no longer match old `main`). Tuning lives entirely in the `OROG_*` pair — this keeps the "exactly two gated sites" invariant. If the `OROG_*` pair alone cannot reach the target, STOP and raise it with Matt (gating the decay as a third site is a scope change, not a silent addition), do not edit the decay.

**Interfaces:**
- Consumes: the gated orographic branch (Task 4).
- Produces: final `OROG_*` values. No signature change — same names, new numbers.

**Target (spec §6):** land moisture `< 0.15` share drops from the current ~42% toward Earth-like ~33% arid+semiarid, WITHOUT flipping the map to all-forest. Grassland and temperate forest should grow; steppe falls to a realistic minority.

- [ ] **Step 1: Record the BEFORE numbers**

Run, capturing output to a scratch note:

```bash
node scripts/queryWorld.mjs climate --seed=realmgenesis --points=20000
node scripts/queryWorld.mjs biomes  --seed=realmgenesis --points=20000
```

Repeat for two more seeds (e.g. `--seed=alpha`, `--seed=delta`). Record the `< 0.15` land-moisture share and the biome-share table for each.

- [ ] **Step 2: Adjust constants and re-measure**

Iterate on `utils/datum.ts` `OROG_*` in the priority order from spec §6:
1. Soften the leeward floor first (`OROG_LEEWARD_FLOOR`, `OROG_LEEWARD_PER_KM`) — the old `*0.2` lost 80%/step and stranded interiors.
2. Keep `OROG_WINDWARD_CAP` bounded so coasts do not saturate.
3. If the pair alone cannot reach the target, STOP and raise it with Matt (see the off-limits note above) — do NOT touch the `* 0.98` decay.

After each change re-run the three-seed `climate` + `biomes` measurement.

- [ ] **Step 3: Confirm the target is met without over-wetting**

Stop when, averaged across the three seeds: `< 0.15` land-moisture share is near ~33% (accept ~30–36%), grassland/forest shares rose, steppe fell, and no seed collapsed to a single dominant biome. Record the final numbers.

- [ ] **Step 4: Re-run the coupling unit tests**

Run: `npx vitest run tests/climateCoupling.test.ts`
Expected: PASS (the retune must not break the live-ness assertions).

- [ ] **Step 5: Commit**

```bash
git add utils/datum.ts
git commit -m "tune(D8b): retune orographic constants to Earth-like inland moisture"
```

---

## Task 6: sync-effect removal, paramLiveness, Controls toggle

**Files:**
- Modify: `hooks/useWorldEngine.ts` (the D1/D8a sync effect ~L132-153)
- Modify: `tests/paramLiveness.test.ts` (allowlist comment ~L13, `TERRAIN_PERTURBATIONS` ~L21-46)
- Modify: `components/Controls.tsx` (Climate tab, near the Ocean Currents toggle ~L1012-1023)

**Interfaces:**
- Consumes: `physicalClimate` live in generation (Tasks 3-4).
- Produces: liveness coverage for `maxElevationM` (now a generation input) and `physicalClimate`; a UI toggle bound via `handleChange('physicalClimate', ...)`; `maxElevationM` no longer live-synced (regenerate to apply).

- [ ] **Step 1: Remove `maxElevationM` from the live sync effect**

`maxElevationM` is now a generation param (Matt's call). Leaving it in the live sync would let the slider rescale readouts without regenerating — the determinism hole (`world.params.maxElevationM` would not match the datum that baked the world). Remove it from the effect; keep `season`.

Replace `hooks/useWorldEngine.ts:132-153`:

```ts
  // D1/D8a: season and maxElevationM are render-only params — they must take
  // effect without regenerating. The viewers and Inspector read world.params (a
  // generation snapshot), so push the live values into world.params here, keeping
  // world.cells identity so WorldMesh geometry is reused and only colors/readouts
  // recompute (paint-stroke pattern). maxElevationM only rescales displayed
  // metres; it touches no terrain.
  useEffect(() => {
    setWorld(prev =>
      prev &&
      (prev.params.season !== params.season ||
        prev.params.maxElevationM !== params.maxElevationM)
        ? {
            ...prev,
            params: {
              ...prev.params,
              season: params.season,
              maxElevationM: params.maxElevationM,
            },
          }
        : prev,
    );
  }, [params.season, params.maxElevationM]);
```

with:

```ts
  // D1: season is render-only — it takes effect without regenerating. The
  // viewers and Inspector read world.params (a generation snapshot), so push the
  // live value into world.params here, keeping world.cells identity so WorldMesh
  // geometry is reused and only colors recompute (paint-stroke pattern).
  // D8b: maxElevationM is NO LONGER synced here — it became a generation param
  // (it drives lapse + orographic), so it must regenerate to take effect. A live
  // push would make world.params.maxElevationM disagree with the datum that baked
  // the world (a determinism hole on save/reload).
  useEffect(() => {
    setWorld(prev =>
      prev && prev.params.season !== params.season
        ? { ...prev, params: { ...prev.params, season: params.season } }
        : prev,
    );
  }, [params.season]);
```

- [ ] **Step 2: Typecheck**

Run: `npm run typecheck`
Expected: zero errors.

- [ ] **Step 3: Write the paramLiveness change**

In `tests/paramLiveness.test.ts`:

- Remove the `maxElevationM` line from the display-only allowlist comment (~L13). Add a note that under D8b default-ON it drives lapse + orographic.
- Add to `TERRAIN_PERTURBATIONS` (~L21-46):

```ts
  maxElevationM: { maxElevationM: 4500 }, // D8b: datum drives lapse + orographic when physicalClimate is on (default)
  physicalClimate: { physicalClimate: false }, // D8b: toggling off changes the climate signature (default is true)
```

- [ ] **Step 4: Run the terrain-liveness test**

Run: `npx vitest run tests/paramLiveness.test.ts -t "terrain params change"`
Expected: PASS — both new perturbations change the terrain signature. If `maxElevationM` does NOT change it, `physicalClimate` default is not reaching the engine — re-check Task 2/3.

Note: this file is the LOAD CANARY. If it times out with zero assertion failures, run `uptime` before suspecting the code; re-run isolated.

- [ ] **Step 5: Add the Controls toggle**

In `components/Controls.tsx`, in the Climate tab after the Ocean Currents block (~L1023), following the checkbox pattern at L564-576:

```tsx
              <div className="flex items-center justify-between text-xs text-ink-muted pt-2">
                <div className="flex flex-col">
                  <label>Physical Climate (grounded)</label>
                  <span className="text-[10px] text-ink-faint">
                    Real lapse rate + rain shadow from the elevation datum. Off = classic formulas.
                    Changes climate and civ layout for existing worlds.
                  </span>
                </div>
                <input
                  type="checkbox"
                  checked={params.physicalClimate ?? true}
                  onChange={(e) => { handleChange('physicalClimate', e.target.checked); }}
                  className="bg-surface-hover"
                />
              </div>
```

- [ ] **Step 6: Typecheck + lint**

Run: `npm run typecheck && npm run lint`
Expected: zero errors, ≤30 warnings.

- [ ] **Step 7: Commit**

```bash
git add hooks/useWorldEngine.ts tests/paramLiveness.test.ts components/Controls.tsx
git commit -m "feat(D8b): maxElevationM/physicalClimate become generation params + Climate toggle"
```

---

## Task 7: Verification — byte-identical hatch, full gate, docs

**Files:**
- Uses: `scripts/captureBaseline.mjs` (byte-identical gate)
- Modify: `HANDOFF.md`, `ROADMAP.md` (D8b status)

**Interfaces:**
- Consumes: everything above.
- Produces: recorded evidence that the hatch is byte-identical and all gates pass.

- [ ] **Step 1: Prove the escape hatch is byte-identical (worktree before/after)**

The escape hatch is byte-identical BY CONSTRUCTION (verbatim else branches). This step confirms it empirically. A committed fixture is NOT used — `captureBaseline.mjs`'s header explains why (V8 `Math.*` drift makes a committed bit-exact baseline dishonest). `digestWorld` (`tests/helpers/worldDigest.ts`) is the right instrument: it hashes per-cell fields with exact IEEE-754 bits (not `toFixed`), so a match is a true bit-exact match, not a coarse-histogram match.

Capture the BEFORE from a pristine pre-D8b worktree, then compare the AFTER with the hatch forced OFF:

```bash
# 1. Capture BEFORE from a pristine pre-D8b worktree (the SHA of main before Task 1).
mkdir -p tmp
git worktree add --detach /tmp/d8b-before <pre-D8b-SHA>
ln -s "$(pwd)/node_modules" /tmp/d8b-before/node_modules
( cd /tmp/d8b-before && node scripts/captureBaseline.mjs before && cp tmp/baseline-before.json "$OLDPWD/tmp/" )
```

The `before` worlds use `makeParams()`, which pre-D8b has no `physicalClimate` key (→ the old path). To reproduce that path AFTER, generate with the hatch OFF. This project is ESM (`package.json` may set `"type": "module"`), so `require` is not available — write a `.mjs` scratch script, not `node -e`:

Create `tmp/d8b-hatch-check.mjs`:

```js
import { createServer } from 'vite';
import { readFileSync } from 'node:fs';

const s = await createServer({ server: { middlewareMode: true } });
const { generateWorld } = await s.ssrLoadModule('/utils/worldGen.ts');
const { digestWorld } = await s.ssrLoadModule('/tests/helpers/worldDigest.ts');
const { makeParams } = await s.ssrLoadModule('/tests/helpers.ts');

const before = JSON.parse(readFileSync('tmp/baseline-before.json', 'utf8')).digest;
const out = {};
for (const [name, p] of [
  ['small', makeParams({ physicalClimate: false })],
  ['medium', makeParams({ points: 5000, erosionIterations: 3, physicalClimate: false })],
]) {
  const d = digestWorld(await generateWorld(p));
  for (const [k, v] of Object.entries(d)) out[`${name}.${k}`] = v;
}
await s.close();

const diffs = Object.keys(before).filter(k => before[k] !== out[k]);
console.log(diffs.length ? `DIFFERS: ${diffs.join(', ')}` : 'IDENTICAL — hatch reproduces pre-D8b main');
process.exit(diffs.length ? 1 : 0);
```

Run: `node tmp/d8b-hatch-check.mjs`

Expected: `IDENTICAL`. If it differs, an else branch drifted from the old expression — fix the verbatim math before proceeding.

Clean up: `git worktree remove /tmp/d8b-before` and `rm tmp/d8b-hatch-check.mjs` (scratch script — do not commit it, per the "clean up scratchpad scripts" rule).

- [ ] **Step 2: Run the full gate**

Do NOT run with a browser open (M1 cannot carry both — HANDOFF trap). Redirect the test log to a file, do not pipe through `tail`.

```bash
npm run typecheck
npm run lint
npm test > /tmp/d8b-test.log 2>&1; echo "exit: $?"
npm run build
```

Expected: typecheck 0, lint 0 errors / ≤30 warnings, tests green, build OK. If `paramLiveness` alone times out, run `uptime`; re-run it isolated before suspecting code.

- [ ] **Step 3: Browser verification**

Use the S25b headless Playwright recipe (HANDOFF "Playwright recipe"). Reuse Matt's `:3000` server — never kill it. Generate a default world (hatch on), confirm: no console errors, cold high peaks, wetter interiors, more grassland/forest, less steppe. Then toggle Physical Climate off and regenerate — confirm the classic climate returns. Kill the chromium tree after: `pkill -f chrome-headless-shell`.

- [ ] **Step 4: Update HANDOFF and ROADMAP**

- `ROADMAP.md` D8b: `⬜ TODO` → `✅ DONE`, with the final measured biome/moisture shift and the byte-identical-hatch confirmation.
- `HANDOFF.md`: new session entry — what shipped, the tuned `OROG_*` values, the measured before/after moisture numbers, the byte-identical result, gate state. Note the deferred items (air-temperature orographic factor, volcanic decoupling) stayed out of scope.

- [ ] **Step 5: Commit**

```bash
git add HANDOFF.md ROADMAP.md
git commit -m "docs(D8b): mark D8b done, record measured shift + byte-identical hatch"
```

---

## Self-Review

**Spec coverage:**
- §2 escape hatch → Task 2 (param) + verbatim else branches (Tasks 3-4) + Task 7 Step 1 (proof). ✓
- §4 lapse rate → Task 3. ✓
- §5 orographic → Task 4 (land-side metres, ocean→0, cap/floor). ✓
- §6 moisture retune → Task 5 (empirical, measured target). ✓
- §7 snow line → emergent, no task needed (determineBiome untouched — asserted in Global Constraints). ✓
- §8 plumbing → Task 2 (types/defaults/export) + Task 6 (sync-effect removal, paramLiveness, Controls). ✓
- Blocker resolved (advisor + Matt, 2026-08-28): the D8a live `maxElevationM` sync collides with it becoming a generation input → Matt chose "generation param (regen)"; Task 6 Step 1 removes it from the sync effect. ✓
- §9 testing → Task 1/2/3/4 unit tests, Task 5 empirical gate, Task 7 byte-identical + full gate + browser. ✓
- §11 out of scope → recorded in Task 7 Step 4. ✓

**Placeholder scan:** every code step shows the exact edit; the one intentionally open variable (OROG_* values) is scoped to Task 5 with a measured target and a starting value, not a "TODO". ✓

**Type consistency:** `physicalClimate: boolean` used identically in types.ts, defaultParams, export validator/defaults, worldGen guards, paramLiveness, Controls. `elevationMetres(height, seaLevel, maxElevationM)` matches the datum.ts signature. Constant names identical between Task 1 (definition) and Tasks 3-4 (use). ✓

**Known risk carried forward:** the ~21-27% biome-change figure is a FLOOR (spec §4 note) — Task 7 re-measures rather than asserting it. The moisture tune is the only open variable; the hatch guarantees a bad tune is never forced on anyone.
