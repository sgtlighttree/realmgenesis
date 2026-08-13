# D7 Polish — Plate-Shape Realism Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make tectonic plates read as tectonically shaped rather than Voronoi-blobby, via three additive heuristic-tier changes, with no rebuild of the plate-assignment substrate.

**Architecture:** All engine work is in `utils/tectonicsV3.ts`. Item 1 (band seeding) grows each plate's Dijkstra source from one cell into a connected chain along the plate's Euler-pole velocity — elongated seed, elongated plate, connectivity preserved by construction. Item 3 (age perturbation) breaks GDH1's clean age bands with fresh-stream noise. Item 2 (transform fracture) is conditional — built only if 1+3 don't kill the Voronoi read — and rides the existing multi-timestep loop, feeding step *k*'s transform classification into step *k+1*'s cost field.

**Tech Stack:** TypeScript (strict), React 19, Vitest. Deterministic seeded RNG (`utils/rng.ts`), spherical vector helpers (`utils/spherical.ts`).

## Global Constraints

- **Relative imports only** — `@/` alias is intentionally unused.
- **Determinism:** every new random draw uses a fresh RNG side-stream (`new RNG(seed + '_<purpose>')`). **No new draw may touch `plateRng` or `microRng` order.** Same seed twice → byte-identical `height` + `plateId`.
- **0-exclave invariant:** every plate must remain exactly one connected component over the cell neighbor graph. Verify with the connectivity probe (Task 1).
- **Gates (all must pass before merge):** `npm run typecheck` (0 errors), `npm run lint` (0 errors / ≤30 warnings), `npm test`, `npm run build`.
- **Param liveness:** any new `WorldParams` key must influence generated output (`tests/paramLiveness.test.ts`).
- **Known flake:** `paramLiveness` can *time out* under 23-file parallel load on M1 (Session 10) — a timeout is the flake; an assertion failure is real. Re-run in isolation (`npx vitest run tests/paramLiveness.test.ts`) to disambiguate.
- **Re-baseline is expected:** `plateElongation` default 0.4 changes plate shapes for all seeds; absolute-value seed baselines (lakes) re-scan. Relative signature-change tests (paramLiveness, determinism) keep working.

---

### Task 1: Band / chain seeding (`plateElongation`)

The dominant fix. Adds the param, threads it into `assignPlatesDijkstra`, grows each plate's single source cell into a connected chain along its velocity direction, and proves connectivity + liveness.

**Files:**
- Modify: `types.ts` (add `plateElongation` to `WorldParams`)
- Modify: `hooks/useWorldEngine.ts` (default `plateElongation: 0.4`)
- Modify: `tests/helpers.ts` (default `plateElongation: 0.4`)
- Modify: `utils/tectonicsV3.ts` — `assignPlatesDijkstra` (296) signature + seed loop; both call sites (567, 624); thread `plateElongation` from `simulateTectonics` params
- Modify: `components/Controls.tsx` (Geo tab slider)
- Modify: `tests/paramLiveness.test.ts` (add `plateElongation` case)
- Create: `tests/plateConnectivity.test.ts` (0-exclave probe)

**Interfaces:**
- Consumes: `PlateState` (`{ id, eulerPole: {axis: Point, rate: number}, dominantCrustType, seedPosition }`), `cross3`, `scale3`, `normalizeVec`, `dot3`, `magnitude` from `./spherical`.
- Produces: `assignPlatesDijkstra(macroPoints, macroNeighbors, edgeCosts, rotatedSeeds, plateSpeeds, activeMask, plateIds, plates, elongation)` — two new trailing params `plates: PlateState[]`, `elongation: number`.

- [ ] **Step 1: Write the failing connectivity test**

Create `tests/plateConnectivity.test.ts`:

```ts
import { describe, it, expect } from 'vitest';
import { generateWorld } from '../utils/worldGen';
import { makeParams } from './helpers';

// Every plateId must form exactly one connected component over cell.neighbors.
function strayComponents(cells: { plateId: number; neighbors: number[] }[]): number {
  const seen = new Uint8Array(cells.length);
  const compsByPlate = new Map<number, number>();
  for (let i = 0; i < cells.length; i++) {
    if (seen[i]) continue;
    const pid = cells[i].plateId;
    const stack = [i];
    seen[i] = 1;
    while (stack.length) {
      const c = stack.pop()!;
      for (const n of cells[c].neighbors) {
        if (!seen[n] && cells[n].plateId === pid) { seen[n] = 1; stack.push(n); }
      }
    }
    compsByPlate.set(pid, (compsByPlate.get(pid) ?? 0) + 1);
  }
  let stray = 0;
  for (const count of compsByPlate.values()) stray += count - 1;
  return stray;
}

describe('band seeding preserves plate connectivity', () => {
  it('has 0 stray plate components at elongation 0.4', async () => {
    const world = await generateWorld(makeParams({ seed: 'realmgenesis', plateElongation: 0.4 }));
    expect(strayComponents(world.cells as never)).toBe(0);
  }, 120000);
});
```

Then write the genuine TDD-**RED** gate in `tests/paramLiveness.test.ts` (this is the assertion that honestly fails pre-implementation — before Steps 4–5, `plateElongation` is ignored, so 0.0 and 1.0 produce identical worlds). Add inside the V3-params test:
```ts
    expect(terrainSignature(await generateWorld(makeParams({ plateElongation: 0.0 }))))
      .not.toBe(terrainSignature(await generateWorld(makeParams({ plateElongation: 1.0 }))));
```

- [ ] **Step 2: Run both — connectivity GREEN (guard), liveness RED**

Run: `npx vitest run tests/plateConnectivity.test.ts tests/paramLiveness.test.ts`
Expected: `plateConnectivity` **PASSES** — it is a guard that must hold on today's code and after Step 4 (vitest does not typecheck, so the extra `plateElongation` key is inert at runtime; the world generates and connectivity already holds). `paramLiveness` **FAILS** on the new `.not.toBe` — identical worlds until the algorithm lands. That failure is the real RED gate; do not "fix" it by touching the type or the test. (If `paramLiveness` times out instead of asserting, re-run it in isolation — Session-10 flake.)

- [ ] **Step 3: Add the param to the type + defaults**

`types.ts` — in `WorldParams`, near `microplateIntensity`:
```ts
  plateElongation: number; // 0–1: seed-chain length → plate elongation (0 = round blobs)
```
`hooks/useWorldEngine.ts` — in the defaults object, near `microplateIntensity`:
```ts
  plateElongation: 0.4,
```
`tests/helpers.ts` — in the default params, near `microplateIntensity`:
```ts
  plateElongation: 0.4,
```

- [ ] **Step 4: Grow the seed chain in `assignPlatesDijkstra`**

`utils/tectonicsV3.ts` — extend the signature (296) with two trailing params:
```ts
function assignPlatesDijkstra(
  macroPoints: Point[],
  macroNeighbors: number[][],
  edgeCosts: number[][],
  rotatedSeeds: Point[],
  plateSpeeds: number[],
  activeMask: (j: number) => boolean,
  plateIds: Int32Array,
  plates: PlateState[],
  elongation: number,
): void {
```

Replace the single-source seed loop (315–325) with a chain seeder. **A `claimed` set shared across the whole seed loop is load-bearing** (see the exclave note below):
```ts
  const chainLen = 1 + Math.round(elongation * 4);
  // Cells already used as a dist-0 source by an EARLIER plate's chain. Disjoint
  // source sets across plates are what keep connectivity by construction: two
  // plates must never both seed the same cell, or one plate's chain gets severed
  // into an exclave. A dist-0 pop always beats any positive-dist front, so every
  // unclaimed source settles for its own plate.
  const claimed = new Uint8Array(n);
  for (let j = 0; j < rotatedSeeds.length; j++) {
    if (!activeMask(j)) continue;
    // nearest UNCLAIMED macro cell to the (rotated) seed
    let best = -1, bestD = Infinity;
    const s = rotatedSeeds[j];
    for (let i = 0; i < n; i++) {
      if (claimed[i]) continue;
      const d = chordDistance(macroPoints[i], s);
      if (d < bestD) { bestD = d; best = i; }
    }
    if (best < 0) continue;

    // chain axis = plate velocity direction at the seed (tangent to sphere)
    const vel = scale3(cross3(plates[j].eulerPole.axis, s), plates[j].eulerPole.rate);
    const axis = magnitude(vel) > 1e-9 ? normalizeVec(vel) : null;

    // walk a connected chain of UNCLAIMED cells forward along the axis; all become dist-0 sources
    claimed[best] = 1;
    heap.push({ cell: best, plate: j, dist: 0 });
    let tip = best;
    let chainSize = 1;
    while (axis && chainSize < chainLen) {
      let nextCell = -1, nextDot = 0; // require forward progress (dot > 0)
      for (const nId of macroNeighbors[tip]) {
        if (claimed[nId]) continue;
        const dir = normalizeVec({
          x: macroPoints[nId].x - macroPoints[tip].x,
          y: macroPoints[nId].y - macroPoints[tip].y,
          z: macroPoints[nId].z - macroPoints[tip].z,
        });
        const d = dot3(dir, axis);
        if (d > nextDot || (d === nextDot && nextCell >= 0 && nId < nextCell)) { nextDot = d; nextCell = nId; }
      }
      if (nextCell < 0) break;
      claimed[nextCell] = 1;
      heap.push({ cell: nextCell, plate: j, dist: 0 });
      tip = nextCell;
      chainSize++;
    }
  }
```
> **Why `claimed` (fable-advisor catch):** chains are 3–5 cells and proto-plate/microplate seeds are dense and boundary-adjacent, so two plates' chains can otherwise both push the same cell B at dist 0. The earlier plate wins the pop (plate-index tiebreak, line 312), severing the later plate's chain into two independently-grown regions with a foreign plate between them → **exclave**. It is seed-dependent, so a single-seed test can pass while other seeds violate the invariant. Sharing `claimed` across the seed loop makes the per-plate source sets disjoint, restoring connectivity by construction.
>
> Note: `n`, `heap`, `chordDistance` are already in scope. The chain walk is pure geometry — **no RNG draw** — so `plateRng` order is untouched and `elongation === 0` (`chainLen === 1`) reproduces today's single-source seeding byte-identically (the single nearest cell is never pre-claimed since each plate's nearest differs; if a rare tie pre-claims it, the plate falls to its next-nearest — still deterministic).

- [ ] **Step 5: Update both call sites + thread the param**

In `simulateTectonics`, read the param once near the other param reads (~548):
```ts
  const plateElongation = params.plateElongation ?? 0.4;
```
Call site 1 (567–570):
```ts
  assignPlatesDijkstra(
    macroPoints, macroNeighbors, edgeCosts, rotatedSeeds, plateSpeeds,
    () => true, plateIds, plates, plateElongation,
  );
```
Call site 2 (624–627):
```ts
    assignPlatesDijkstra(
      macroPoints, macroNeighbors, edgeCosts, rotatedSeeds, plateSpeeds,
      j => postMergeCounts[j] > 0, plateIds, plates, plateElongation,
    );
```

- [ ] **Step 6: Run the connectivity test — expect PASS**

Run: `npx vitest run tests/plateConnectivity.test.ts`
Expected: PASS (0 stray components). If it fails with stray > 0, a chain crossed into a disconnected region — check the `chain.has(nId)` guard and that all chain cells share `plate: j`.

- [ ] **Step 7: Confirm the liveness gate is now GREEN**

The `plateElongation` `.not.toBe` assertion was written in Step 1. Re-run it — now that Steps 4–5 make the param live, it must PASS:

Run: `npx vitest run tests/paramLiveness.test.ts`
Expected: PASS (elongation 0.0 vs 1.0 now produce different terrain signatures). This RED→GREEN transition is the proof the param is wired.

- [ ] **Step 8: Add the Controls slider**

`components/Controls.tsx` — in the Geo/Advanced tab, next to the Session-10 Microplates slider, add a `plateElongation` slider (0–1, step 0.05) following the exact pattern of the adjacent slider (label "Plate Elongation", `accent-brand-soft`, functional `setParams(prev => ({ ...prev, plateElongation: v }))`). Copy the neighboring slider block verbatim and change the key, label, and range.

- [ ] **Step 9: Re-baseline lakes seeds**

Run: `npx vitest run tests/lakes.test.ts`
If it fails (absolute lake cells shifted by the new terrain), re-scan the seeds exactly as Session 9/10 did: read the new lake output the test logs and update the expected cell counts/ids in `tests/lakes.test.ts`. **Do not** loosen assertions — update the expected values to the new deterministic output.

- [ ] **Step 10: Run all gates**

Run: `npm run typecheck && npm run lint && npx vitest run tests/plateConnectivity.test.ts tests/paramLiveness.test.ts tests/lakes.test.ts && npm run build`
Expected: typecheck 0, lint 0 errors, the three test files pass, build OK. (Full `npm test` in Step 11.)

- [ ] **Step 11: Full suite + isolation re-run if flake**

Run: `npm test`
Expected: all pass. If `paramLiveness` throws a **timeout**, re-run `npx vitest run tests/paramLiveness.test.ts` — isolation pass confirms the Session-10 flake, not a regression.

- [ ] **Step 12: Commit**

```bash
git add types.ts hooks/useWorldEngine.ts tests/helpers.ts utils/tectonicsV3.ts components/Controls.tsx tests/paramLiveness.test.ts tests/plateConnectivity.test.ts tests/lakes.test.ts
git commit -m "d7 polish: band/chain plate seeding (plateElongation)"
```

---

### Task 2: Seafloor age-band perturbation

Breaks GDH1's artificially clean, symmetric age bands. Cosmetic, bathymetry-only, fresh side-stream.

**Files:**
- Modify: `utils/tectonicsV3.ts` — `computeSeafloorAge` (the `age[i] = Math.min(...)` block, ~435–440) + its `RNG`/`SimplexNoise` wiring
- Modify: `tests/paramLiveness.test.ts` (only if a signature re-baseline is needed — see Step 4)

**Interfaces:**
- Consumes: `RNG` (`../rng`), `SimplexNoise` (`../rng`), `macroPoints: Point[]`, the plumbed `params.seed`.
- Produces: no new exported symbol; perturbs `age` in place.

- [ ] **Step 1: Thread the seed into `computeSeafloorAge`**

Confirmed: the current signature ends `..., rotatedSeeds: Point[], spreadRate: number)` — no seed. Add a trailing `seed: string` param:
```ts
function computeSeafloorAge(
  macroPoints: Point[],
  macroNeighbors: number[][],
  plateIds: Int32Array,
  crustTypes: Uint8Array,
  plates: PlateState[],
  rotatedSeeds: Point[],
  spreadRate: number,
  seed: string,
): Float32Array {
```
And update the call site (~707): `computeSeafloorAge(macroPoints, macroNeighbors, plateIds, crust.crustTypes, plates, rotatedSeeds, spreadRate, params.seed)`.

- [ ] **Step 2: Add the perturbation**

At the top of `computeSeafloorAge`, build a fresh simplex:
```ts
  const ageNoise = new SimplexNoise(new RNG(seed + '_agenoise_v3'));
```
Replace the final age assignment (currently `age[i] = Math.min(MAX_SEAFLOOR_AGE, d / rate);`):
```ts
    const raw = d / rate;
    const p = macroPoints[i];
    const perturbed = raw * (1 + 0.1 * ageNoise.noise3D(p.x * 2, p.y * 2, p.z * 2));
    age[i] = Math.max(0, Math.min(MAX_SEAFLOOR_AGE, perturbed));
```

- [ ] **Step 3: Determinism sanity — same seed twice is identical**

Run: `npx vitest run tests/worldGen.test.ts`
Expected: PASS — the existing determinism test compares two in-process runs; the fresh stream is deterministic so it stays green.

- [ ] **Step 4: Re-baseline if needed + gates**

Run: `npm run typecheck && npm run lint && npx vitest run tests/lakes.test.ts && npm run build`
If `lakes` shifts again (ocean-height change can nudge coastal lakes), re-scan as in Task 1 Step 9. Then `npm test` for the full suite.

- [ ] **Step 5: Commit**

```bash
git add utils/tectonicsV3.ts tests/lakes.test.ts
git commit -m "d7 polish: perturb seafloor age bands (_agenoise_v3)"
```

---

### CHECKPOINT (manual, with maintainer)

**Before Task 3:** render seed `realmgenesis` in Plates + Height views (2D Mercator) and compare against the pre-Task-1 baseline screenshots in the session scratchpad. Judge: does band seeding alone kill the Voronoi read?
- **If yes (likely — it's the dominant term):** Task 3 is optional cheap polish. Decide with the maintainer whether to proceed.
- **If Task 3's expected jaggedness would be indistinguishable from high `boundaryRoughness`:** **drop Task 3** (advisor fallback (d)) and close D7. Do not chase it.

Proceed to Task 3 only on an explicit go.

---

### Task 3 (CONDITIONAL): Transform-edge fracture

Rides the existing timestep loop: feed step *k*'s transform classification into step *k+1*'s cost field. **Recompute the effective cost field as a SET each step — never accumulate**, or boundaries run away.

**Files:**
- Modify: `utils/tectonicsV3.ts` — `simulateTectonics` timestep loop (609–onwards): add a per-edge transform flag captured in 4c, and a per-step effective-cost rebuild consumed by 4b's `assignPlatesDijkstra`
- Create: `tests/transformFracture.test.ts` (connectivity still holds)

**Interfaces:**
- Consumes: the static `edgeCosts: number[][]` (553), `macroNeighbors`, `RNG` (`_fracture_v3`).
- Produces: `effectiveEdgeCosts: number[][]` rebuilt each step; passed to `assignPlatesDijkstra` in place of the static `edgeCosts` inside the loop.

- [ ] **Step 1: Write the failing connectivity test**

Create `tests/transformFracture.test.ts` mirroring `tests/plateConnectivity.test.ts` (import the same `strayComponents` helper — extract it to `tests/helpers.ts` and import in both), asserting 0 stray components at `microplateIntensity: 0.35, plateElongation: 0.4` (defaults). This guards that time-varying costs don't break region-growth connectivity.

- [ ] **Step 2: Run it — expect PASS already (guard test)**

Run: `npx vitest run tests/transformFracture.test.ts`
Expected: PASS on current code (region-growth is connected). This test exists to FAIL if Step 3 regresses connectivity.

- [ ] **Step 3: Build the per-edge fracture multiplier + transform flags**

Before the loop (after `edgeCosts` at 556):
```ts
  const fractureRng = new RNG(params.seed + '_fracture_v3');
  // Per-directed-edge fracture multiplier, keyed by a stable per-edge draw.
  const fractureMul: number[][] = macroNeighbors.map(row => row.map(() => 1 + fractureRng.next() * 1.5));
  // Transform flag for each directed edge, from the PREVIOUS step's classification.
  const isTransform: Uint8Array[] = macroNeighbors.map(row => new Uint8Array(row.length));
  // Effective (time-varying) cost field, rebuilt as a SET each step.
  const effectiveEdgeCosts: number[][] = macroNeighbors.map(row => row.map(() => 0));
```

- [ ] **Step 4: Rebuild effective costs at the top of each step, feed to 4b**

At the start of the loop body (before 4b's `assignPlatesDijkstra`), rebuild from static costs — never from the previous effective value:
```ts
    for (let i = 0; i < numMacro; i++) {
      const row = effectiveEdgeCosts[i];
      for (let t = 0; t < row.length; t++) {
        row[t] = edgeCosts[i][t] * (isTransform[i][t] ? fractureMul[i][t] : 1);
      }
    }
```
Change 4b's call to pass `effectiveEdgeCosts` instead of `edgeCosts`.

- [ ] **Step 5: Set transform flags in 4c — indexed loop + per-cell reset**

Two edits to 4c, both load-bearing:

**(a) Convert 4c's neighbor loop from `for...of` to indexed**, so a directed-edge index `t` exists. The loop currently reads `for (const nId of macroNeighbors[i]) {` (~637). Change to:
```ts
      const nRow = macroNeighbors[i];
      for (let t = 0; t < nRow.length; t++) {
        const nId = nRow[t];
        if (plateIds[nId] === plateIds[i]) continue;
        // ... existing body unchanged (relV, edgeNormal, vn, vnMag, vtMag, classification) ...
```

**(b) Reset the whole row before the neighbor loop, then set only cross-plate transform edges.** Because 4c `continue`s on same-plate edges (638), a flag left set from a prior step would never be cleared when an edge stops being a transform — the accumulation failure the task header warns against. At the **top of each cell `i`'s 4c body**, before the neighbor loop:
```ts
      isTransform[i].fill(0);
```
Then, inside the loop after the existing `vn` (650) and `vtMag` (651) are computed, reuse those locals (do NOT recompute — one transform definition, consistent with `injectMicroplates`' shear at line 474):
```ts
        // transform = tangential relative velocity dominates the normal component
        isTransform[i][t] = vtMag > vnMag ? 1 : 0;
```
> `vnMag` (= `Math.abs(vn)`) and `vtMag` are already computed at 649–651. The flag is thus **rebuilt as a set** each step: cleared by `fill(0)`, then set to 1 only on current transform edges. Edges that stop being transforms revert to 0 next step; same-plate edges (skipped by `continue`) stay 0.

- [ ] **Step 6: Run connectivity + determinism**

Run: `npx vitest run tests/transformFracture.test.ts tests/worldGen.test.ts`
Expected: both PASS. Stray > 0 means the effective-cost swing broke a front — verify costs are always ≥ the static base (fractureMul ≥ 1) so no negative/zero-cost shortcut appears.

- [ ] **Step 7: Re-baseline + full gates**

Run: `npm run typecheck && npm run lint && npx vitest run tests/lakes.test.ts && npm run build && npm test`
Re-scan lakes if shifted. Isolation re-run on any paramLiveness timeout.

- [ ] **Step 8: Commit**

```bash
git add utils/tectonicsV3.ts tests/helpers.ts tests/transformFracture.test.ts tests/plateConnectivity.test.ts tests/lakes.test.ts
git commit -m "d7 polish: transform-edge fracture via time-varying cost field (_fracture_v3)"
```

---

## Post-implementation

- Update `HANDOFF.md` (Session entry) and `ROADMAP.md` (D7 tag → done, Cortial recorded as deliberate NO-GO).
- Re-run the 0-exclave probe across `realmgenesis` / `route-test` / `abcxyz` as a final invariant check.
- Screenshot the final result for the maintainer.

## Self-review notes

- **Spec coverage:** Item 1 → Task 1; Item 3 → Task 2; Item 2 → Task 3 (conditional, with the checkpoint gate the spec mandates). NO-GO decision + re-baseline + determinism rules → Global Constraints.
- **Placeholders:** none — every code step carries real code; the one deliberately-manual step (Controls slider, Task 1 Step 8) points at an exact existing pattern to copy rather than inventing UI.
- **Type consistency:** `assignPlatesDijkstra` gains `plates: PlateState[], elongation: number` (used identically at both call sites); `strayComponents` is defined once (extracted to `tests/helpers.ts` in Task 3 Step 1, used by both connectivity tests); `_agenoise_v3` / `_fracture_v3` streams are distinct and never reuse `plateRng`/`microRng`.
