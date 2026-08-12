# D6 Stage 1 — Worker Migration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Move the entire world-generation pipeline off the main thread into a Web Worker with **zero algorithm change**, proven by bit-exact comparison of main-thread output against worker output.

**Architecture:** `generateWorld` and its transitive imports become worker-safe (no `three`, no `d3`, no DOM). A new `utils/worldTransfer.ts` flattens `WorldData` — an array-of-objects with variable-length `vertices`/`neighbors` — into structure-of-arrays typed buffers that cross the worker boundary as zero-copy Transferables, and rehydrates them into the identical `Cell[]` on the main thread. `utils/worldGenClient.ts` presents the *same call signature* as `generateWorld`, so `useWorldEngine` changes by one import. Abort is by `worker.terminate()`, which lets every `setTimeout(0)` yield be deleted.

**Tech Stack:** TypeScript (strict), Vite 6 (`?worker` import), Web Workers + Transferable ArrayBuffers, Vitest, Playwright (via MCP) for the browser-only gate.

---

## Global Constraints

- **Relative imports only** — no `@/` alias, even though it is configured.
- **Gates, all four, after every task:** `npm run typecheck` = 0 errors; `npm run lint` = 0 errors at the **29-warning ratchet, add none**; `npm test` all green; `npm run build` OK.
- **Zero algorithm change.** This stage must not alter a single generated value. Any step that changes output is a bug in that step, not a design decision to be made here.
- **No new `WorldParams` keys.** `tests/paramLiveness.test.ts` is untouched by this stage.
- **No save-schema change.** `utils/export.ts` is not modified.
- **2-space indent, semicolons, single quotes, trailing commas.** `interface` for objects, `type` for unions.
- **Do NOT push.** Commit locally in focused, single-topic chunks.
- **Do NOT touch** `CLAUDE.md` / `ROADMAP.md`. `HANDOFF.md` is updated in Task 7 only.
- **The dev server on `:3000` may be Matt's.** Do not kill it. Start your own on another port for verification.

---

## Why this stage exists, and what it is not

Stage 1 ships **no visible improvement**. It is sequenced first for exactly one reason: a byte-identical output gate is only possible *before* Stage 2 changes the algorithm on purpose. Once V3 lands, "did the migration preserve behaviour?" becomes an unanswerable question forever.

The second reason is scope discovery. `Cell` is an array-of-objects with ragged `vertices` and `neighbors`; structured-cloning 200k of them is the single biggest risk in the whole D6 project (spec §6, Stage 1). Discovering the flattening contract during Stage 2 — while values are also legitimately changing — would make every bug ambiguous.

**Not in scope:** any change to crust/plate decoupling, Euler poles, seam elimination, coarse simulation, or sub-cell detail. Those are Stage 2 and Stage 3. Do not "improve" the algorithm while you are in there.

---

## Instruments, and the limits of each

This plan uses three different correctness instruments. Each catches something the others cannot, and **each has a stated blind spot**. Do not collapse them into one.

| Instrument | Catches | Blind to |
|---|---|---|
| **Existing determinism suite** (`tests/worldGen.test.ts`) | run-to-run nondeterminism, in-process | any change applied consistently to both runs; 4 fields only, rounded to 6 decimals |
| **`worldDigest` + same-session baseline** (Task 1) | a refactor silently changing values on this machine, this engine, across all fields at full float precision | cross-engine ULP drift; cannot compare main thread to worker |
| **Browser main-vs-worker compare** (Task 5) | serialization loss across the real worker boundary — dropped fields, `Float64`→`Float32`, `undefined`→`0`, `-0`, key-set changes | nothing about whether the *algorithm* is right; only that both paths agree |

**Why no committed golden fixture.** The obvious design is a committed pre-migration `WorldData` snapshot that CI compares against. It was rejected: `Math.sin`/`cos`/`pow` are implementation-defined in ECMAScript, so a bit-exact fixture drifts by a last-ULP across V8 versions and turns into a CI flake that gets `toBeCloseTo`'d into uselessness — at which point it no longer catches a `Float32` downcast, the one thing it was built for. Instead the baseline is captured **into a gitignored temp file on the same engine in the same session** immediately before each risky change. That is the Session 6e method (334 computed styles captured before, compared after) and it has no drift term at all.

**Record this if you are tempted to add the committed fixture back:** the weakness of the temp-file baseline is that it is a discipline, not a gate — an implementer who skips the capture gets no failure. Task 2 and Task 6 therefore make the capture an explicit numbered step, and the digest helper writes a timestamp + git SHA into the baseline file so a stale one is visible.

---

## File structure

**Create:**

| File | Responsibility |
|---|---|
| `utils/palette.ts` | The three color arrays + precomputed folk colors. **Leaf module, zero imports beyond `types.ts`.** Exists to keep `three` out of the worker bundle. |
| `utils/vec.ts` | `normalizeVec` and friends. Leaf module. Exists to keep `d3` out of the worker bundle. |
| `utils/worldTransfer.ts` | `serializeWorld` / `deserializeWorld` — the SoA transfer contract. The load-bearing file of this stage. |
| `workers/worldGen.worker.ts` | Worker entry: receives params, runs `generateWorld`, streams progress/log messages, posts the serialized world with transferables. |
| `utils/worldGenClient.ts` | Main-thread client. Same signature as `generateWorld`. Owns worker lifecycle and terminate-based abort. |
| `tests/helpers/worldDigest.ts` | Full-surface, exact-float digest of a `WorldData`. The measuring instrument. |
| `tests/worldDigest.test.ts` | **Mutation tests proving the instrument fails when it should.** |
| `tests/worldTransfer.test.ts` | In-process round-trip: serialize → deserialize → digest-identical. |
| `dev/goldenCompare.html` | Dev-only harness page: runs main-thread and worker paths in one browser session and prints per-field PASS/FAIL. Not shipped, not in CI. |
| `scripts/captureBaseline.mjs` | Writes a digest baseline to `tmp/` for before/after refactor comparison. |

**Modify:**

| File | Change |
|---|---|
| `utils/colors.ts` | Re-export the palette from `utils/palette.ts`; keep `darkenForFolk` for UI use. No behaviour change. |
| `utils/features.ts` | Import `normalizeVec` from `utils/vec.ts` instead of `utils/geo.ts`. |
| `utils/geo.ts` | Re-export `normalizeVec` from `utils/vec.ts`. |
| `utils/worldGen.ts` | Import palette from `utils/palette.ts`; use `FOLK_COLORS`; delete the `setTimeout(0)` yields (Task 6 only). |
| `hooks/useWorldEngine.ts` | Two `generateWorld` call sites → `generateWorldInWorker`. |
| `vite.config.ts` | `worker: { format: 'es' }`. |
| `.gitignore` | Add `tmp/`. |

---

### Task 1: The digest instrument, and proof that it works

The single most important task in this plan. Everything downstream trusts this file, so it gets mutation-tested first: a digest that silently passes is worse than no digest.

**Files:**
- Create: `tests/helpers/worldDigest.ts`
- Create: `tests/worldDigest.test.ts`

**Interfaces:**
- Consumes: `WorldData`, `Cell`, `BiomeType` from `../types` (via `../../types` from inside `tests/helpers/`).
- Produces:
  - `type WorldDigest = Record<string, string>` — field path → hex digest.
  - `digestWorld(world: WorldData): WorldDigest`
  - `diffDigests(a: WorldDigest, b: WorldDigest): string[]` — field paths that differ or are present in only one side.
  - `CELL_FIELDS: readonly string[]` — the canonical, explicit `Cell` field list. **Not `Object.keys`** — a fixed list is what makes "this field vanished" detectable.

- [ ] **Step 1: Write the instrument**

Create `tests/helpers/worldDigest.ts`:

```ts
import { createHash } from 'node:crypto';
import { WorldData, Cell } from '../../types';

export type WorldDigest = Record<string, string>;

// Canonical Cell field list. Deliberately explicit rather than Object.keys:
// a field silently disappearing across the worker boundary must show up as a
// digest change, and Object.keys would simply stop emitting it on both sides.
export const CELL_FIELDS = [
  'id', 'height', 'plateId', 'temperature', 'moisture', 'biome', 'flux',
  'regionId', 'provinceId', 'isCapital', 'isTown', 'population',
  'cultureId', 'religionId',
] as const;

// Exact IEEE-754 bits, not toFixed. A Float64 -> Float32 downcast of 0.5123456
// survives toFixed(6) but changes these bits, which is precisely the class of
// bug this stage exists to catch.
const f64 = new DataView(new ArrayBuffer(8));
const bits = (n: number): string => {
  f64.setFloat64(0, n);
  return f64.getBigUint64(0).toString(16).padStart(16, '0');
};

// undefined, null, -0 and NaN each get a distinct, stable encoding.
const enc = (v: unknown): string => {
  if (v === undefined) return 'u';
  if (v === null) return 'n';
  if (typeof v === 'number') return Object.is(v, -0) ? '-0' : bits(v);
  if (typeof v === 'boolean') return v ? 'T' : 'F';
  if (typeof v === 'string') return 's' + v;
  if (Array.isArray(v)) return '[' + v.map(enc).join(',') + ']';
  if (typeof v === 'object') {
    const o = v as Record<string, unknown>;
    return '{' + Object.keys(o).sort().map(k => k + ':' + enc(o[k])).join(',') + '}';
  }
  return 'x' + String(v);
};

const hash = (s: string): string => createHash('sha256').update(s).digest('hex').slice(0, 32);

export const digestWorld = (world: WorldData): WorldDigest => {
  const d: WorldDigest = {};
  d['cellCount'] = hash(String(world.cells.length));

  for (const field of CELL_FIELDS) {
    d[`cell.${field}`] = hash(
      world.cells.map(c => enc((c as unknown as Record<string, unknown>)[field])).join(';'),
    );
  }
  d['cell.center'] = hash(world.cells.map(c => enc(c.center)).join(';'));
  d['cell.vertices'] = hash(world.cells.map(c => enc(c.vertices)).join(';'));
  d['cell.neighbors'] = hash(world.cells.map(c => enc(c.neighbors)).join(';'));

  d['geoJson'] = hash(enc(world.geoJson as unknown));
  d['params'] = hash(enc(world.params as unknown));
  d['rivers'] = hash(enc(world.rivers as unknown));
  d['lakes'] = hash(enc(world.lakes as unknown));
  d['features'] = hash(enc(world.features as unknown));
  d['cultures'] = hash(enc(world.cultures as unknown));
  d['religions'] = hash(enc(world.religions as unknown));
  d['civData'] = hash(enc(world.civData as unknown));
  d['markers'] = hash(enc(world.markers as unknown));
  d['routes'] = hash(enc(world.routes as unknown));
  return d;
};

export const diffDigests = (a: WorldDigest, b: WorldDigest): string[] => {
  const keys = new Set([...Object.keys(a), ...Object.keys(b)]);
  return [...keys].filter(k => a[k] !== b[k]).sort();
};
```

- [ ] **Step 2: Write the mutation tests**

Create `tests/worldDigest.test.ts`. These prove the instrument **fails** when it should — a digest nobody has tried to fool is not evidence.

```ts
import { describe, it, expect } from 'vitest';
import { generateWorld } from '../utils/worldGen';
import { digestWorld, diffDigests } from './helpers/worldDigest';
import { makeParams } from './helpers';
import { WorldData, Cell } from '../types';

const clone = (w: WorldData): WorldData =>
  ({ ...w, cells: w.cells.map(c => ({ ...c, center: { ...c.center }, vertices: c.vertices.map(v => ({ ...v })), neighbors: [...c.neighbors] })) });

describe('worldDigest catches what the existing determinism suite cannot', () => {
  it('reports no differences for an unmodified copy', async () => {
    const w = await generateWorld(makeParams());
    expect(diffDigests(digestWorld(w), digestWorld(clone(w)))).toEqual([]);
  }, 30000);

  it('catches a Float64 -> Float32 downcast of height', async () => {
    const w = await generateWorld(makeParams());
    const m = clone(w);
    m.cells.forEach(c => { c.height = Math.fround(c.height); });
    expect(diffDigests(digestWorld(w), digestWorld(m))).toContain('cell.height');
  }, 30000);

  it('catches a 1e-12 perturbation that toFixed(6) hides', async () => {
    const w = await generateWorld(makeParams());
    const m = clone(w);
    m.cells[0].temperature += 1e-12;
    expect(diffDigests(digestWorld(w), digestWorld(m))).toContain('cell.temperature');
  }, 30000);

  it('catches undefined collapsing to 0', async () => {
    const w = await generateWorld(makeParams());
    const m = clone(w);
    m.cells.forEach(c => { if (c.regionId === undefined) (c as Cell).regionId = 0; });
    expect(diffDigests(digestWorld(w), digestWorld(m))).toContain('cell.regionId');
  }, 30000);

  it('catches a dropped neighbor link', async () => {
    const w = await generateWorld(makeParams());
    const m = clone(w);
    m.cells[0].neighbors.pop();
    expect(diffDigests(digestWorld(w), digestWorld(m))).toContain('cell.neighbors');
  }, 30000);

  it('catches a changed geoJson coordinate', async () => {
    const w = await generateWorld(makeParams());
    const m = clone(w);
    m.geoJson = JSON.parse(JSON.stringify(w.geoJson));
    const g = m.geoJson.features.find(f => f.geometry)!;
    g.geometry!.coordinates[0][0][0] += 1e-9;
    expect(diffDigests(digestWorld(w), digestWorld(m))).toContain('geoJson');
  }, 30000);

  it('catches a lost vertices ring', async () => {
    const w = await generateWorld(makeParams());
    const m = clone(w);
    m.cells[5].vertices = [];
    expect(diffDigests(digestWorld(w), digestWorld(m))).toContain('cell.vertices');
  }, 30000);
});
```

- [ ] **Step 3: Run the tests and confirm they pass**

Run: `npx vitest run tests/worldDigest.test.ts`
Expected: 7 passed. If the "unmodified copy" case reports differences, the digest is reading something non-deterministic — fix that before continuing; if it passes but a mutation case does *not* appear in the diff, the instrument has a hole and nothing downstream can be trusted.

- [ ] **Step 4: Write the baseline capture script**

Create `scripts/captureBaseline.mjs`:

```js
// Captures a worldDigest baseline into tmp/ for before/after refactor comparison.
// Deliberately NOT committed as a fixture: Math.sin/cos/pow are
// implementation-defined, so a committed bit-exact baseline drifts across V8
// versions. Same engine, same session, same machine is the only honest use.
// Usage: node scripts/captureBaseline.mjs before
//        node scripts/captureBaseline.mjs after && node scripts/captureBaseline.mjs compare
import { mkdirSync, writeFileSync, readFileSync } from 'node:fs';
import { execSync } from 'node:child_process';
import { createServer } from 'vite';

const label = process.argv[2];
if (!['before', 'after', 'compare'].includes(label)) {
  console.error('usage: node scripts/captureBaseline.mjs before|after|compare');
  process.exit(1);
}

if (label === 'compare') {
  const a = JSON.parse(readFileSync('tmp/baseline-before.json', 'utf8'));
  const b = JSON.parse(readFileSync('tmp/baseline-after.json', 'utf8'));
  console.log(`before: ${a.capturedAt} @ ${a.gitSha}`);
  console.log(`after:  ${b.capturedAt} @ ${b.gitSha}`);
  const keys = [...new Set([...Object.keys(a.digest), ...Object.keys(b.digest)])].sort();
  const diffs = keys.filter(k => a.digest[k] !== b.digest[k]);
  if (diffs.length === 0) { console.log('IDENTICAL — 0 fields differ'); process.exit(0); }
  console.error(`DIFFERS in ${diffs.length} field(s):\n  ${diffs.join('\n  ')}`);
  process.exit(1);
}

const server = await createServer({ server: { middlewareMode: true } });
const { generateWorld } = await server.ssrLoadModule('/utils/worldGen.ts');
const { digestWorld } = await server.ssrLoadModule('/tests/helpers/worldDigest.ts');
const { makeParams } = await server.ssrLoadModule('/tests/helpers.ts');

const worlds = {};
for (const [name, params] of [['small', makeParams()], ['medium', makeParams({ points: 5000, erosionIterations: 3 })]]) {
  const w = await generateWorld(params);
  const d = digestWorld(w);
  for (const [k, v] of Object.entries(d)) worlds[`${name}.${k}`] = v;
}
await server.close();

mkdirSync('tmp', { recursive: true });
writeFileSync(`tmp/baseline-${label}.json`, JSON.stringify({
  capturedAt: new Date().toISOString(),
  gitSha: execSync('git rev-parse HEAD').toString().trim(),
  node: process.version,
  digest: worlds,
}, null, 2));
console.log(`wrote tmp/baseline-${label}.json (${Object.keys(worlds).length} fields)`);
```

- [ ] **Step 5: Gitignore `tmp/` and smoke the script**

Append `tmp/` to `.gitignore`, then run:

```bash
node scripts/captureBaseline.mjs before && node scripts/captureBaseline.mjs after && node scripts/captureBaseline.mjs compare
```

Expected: `IDENTICAL — 0 fields differ`, exit 0. (Both captures ran the same unmodified code, so anything else means `generateWorld` is nondeterministic — stop and investigate, that is a pre-existing bug this stage would otherwise inherit.)

- [ ] **Step 6: Run all four gates**

```bash
npm run typecheck && npm run lint && npm test && npm run build
```

Expected: 0 type errors; 0 lint errors at ≤29 warnings; all tests green (138 existing + 7 new = 145); build OK.

- [ ] **Step 7: Commit**

```bash
git add tests/helpers/worldDigest.ts tests/worldDigest.test.ts scripts/captureBaseline.mjs .gitignore
git commit -m "Add full-surface world digest instrument

The existing determinism suite compares two in-process runs over four
per-cell fields at toFixed(6). That passes a Float32 downcast, a dropped
flux field, and an undefined-to-0 collapse — the exact failure modes a
worker migration introduces.

worldDigest hashes every Cell field plus geoJson, rivers, lakes,
features, cultures, religions and civData at exact IEEE-754 bits, over an
explicit field list so a vanishing field is detectable. Mutation tests
prove it fails on each of those cases rather than assuming it would."
```

---

### Task 2: Make the generation closure worker-safe

`utils/worldGen.ts` transitively imports `three` (via `utils/colors.ts`) and all of `d3` (via `utils/features.ts` → `utils/geo.ts`). Neither breaks in a worker, but both land in the worker bundle and are re-parsed on every worker spawn — and Task 6 spawns one per generation. This task removes both with pure moves.

**Files:**
- Create: `utils/palette.ts`, `utils/vec.ts`
- Modify: `utils/colors.ts`, `utils/geo.ts`, `utils/features.ts`, `utils/worldGen.ts:4`, `utils/worldGen.ts:1135`
- Test: `tests/palette.test.ts`

**Interfaces:**
- Produces:
  - `utils/palette.ts`: `FACTION_COLORS: string[]`, `CULTURE_COLORS: string[]`, `RELIGION_COLORS: string[]`, `FOLK_COLORS: string[]` (index-aligned to `CULTURE_COLORS`).
  - `utils/vec.ts`: `normalizeVec(p: Point): Point` — moved verbatim from `utils/geo.ts`.
- Consumes: nothing from Task 1 at runtime; uses `scripts/captureBaseline.mjs` as its gate.

**The `darkenForFolk` problem, and why `FOLK_COLORS` is a frozen table.** `worldGen.ts:1135` is the only call site, and it is always called on a `CULTURE_COLORS` entry. `darkenForFolk` uses `THREE.Color`, whose `setHex` applies an sRGB→working-colorspace conversion under `THREE.ColorManagement`. Hand-porting that to plain math is exactly the kind of "obviously equivalent" rewrite that silently shifts a value. So: precompute the eight outputs, freeze them into `palette.ts`, and pin them with a test that recomputes via THREE and compares. Byte-identical by construction, and the test fails loudly if a future `three` upgrade changes the conversion.

- [ ] **Step 1: Capture the before-baseline**

```bash
node scripts/captureBaseline.mjs before
```

Expected: `wrote tmp/baseline-before.json`. **Do not skip this.** It is the only evidence this task changed nothing.

- [ ] **Step 2: Print the folk color table**

```bash
node -e "
const { CULTURE_COLORS, darkenForFolk } = await import('./utils/colors.ts');
console.log(JSON.stringify(CULTURE_COLORS.map(darkenForFolk), null, 2));
" --experimental-strip-types 2>/dev/null || npx tsx -e "
import { CULTURE_COLORS, darkenForFolk } from './utils/colors';
console.log(JSON.stringify(CULTURE_COLORS.map(darkenForFolk), null, 2));
"
```

Copy the printed array verbatim into `FOLK_COLORS` in the next step. Do not retype the hex values by hand.

- [ ] **Step 3: Create `utils/palette.ts`**

Move the three arrays out of `utils/colors.ts` verbatim (cut, do not retype) and add the table:

```ts
// Leaf palette module. Imports NOTHING but types, on purpose: utils/worldGen.ts
// pulls this in, and worldGen runs inside a Web Worker where three.js is dead
// weight re-parsed on every spawn. Keep this file dependency-free.

export const FACTION_COLORS = [ /* moved verbatim from utils/colors.ts */ ];
export const CULTURE_COLORS = [ /* moved verbatim from utils/colors.ts */ ];
export const RELIGION_COLORS = [ /* moved verbatim from utils/colors.ts */ ];

// Precomputed darkenForFolk(CULTURE_COLORS[i]), index-aligned.
//
// Generated once from the THREE implementation rather than ported, because
// THREE.Color applies an sRGB -> working-colorspace conversion in setHex that a
// hand-written HSL port would not reproduce exactly. tests/palette.test.ts
// recomputes these via THREE and fails if a three.js upgrade shifts them.
export const FOLK_COLORS = [ /* paste the array printed in Step 2 */ ];
```

- [ ] **Step 4: Create `utils/vec.ts` and rewire the importers**

Cut `normalizeVec` out of `utils/geo.ts` into `utils/vec.ts`:

```ts
import { Point } from '../types';

// Leaf vector math. Split out of utils/geo.ts so utils/features.ts — and
// therefore the worker bundle — does not pull in all of d3.
export const normalizeVec = (p: Point): Point => { /* moved verbatim from utils/geo.ts */ };
```

Then:
- `utils/geo.ts`: add `export { normalizeVec } from './vec';` so existing importers are untouched.
- `utils/features.ts`: change `import { normalizeVec } from './geo';` → `from './vec';`
- `utils/colors.ts`: add `export { FACTION_COLORS, CULTURE_COLORS, RELIGION_COLORS } from './palette';` so existing importers are untouched. Keep `darkenForFolk` here.
- `utils/worldGen.ts:4`: `import { FACTION_COLORS, CULTURE_COLORS, RELIGION_COLORS, darkenForFolk } from './colors';` → `import { FACTION_COLORS, CULTURE_COLORS, RELIGION_COLORS, FOLK_COLORS } from './palette';`
- `utils/worldGen.ts:1135`: `color: darkenForFolk(culture.color),` → `color: FOLK_COLORS[culture.id % FOLK_COLORS.length],`

**Verify the index assumption before making that last edit.** Read the surrounding block and confirm `culture.color` is `CULTURE_COLORS[culture.id % CULTURE_COLORS.length]`. If cultures draw their color any other way (a shuffle, a filtered subset), `FOLK_COLORS[culture.id % …]` is wrong — in that case keep the lookup keyed on the color string: `FOLK_BY_CULTURE_COLOR.get(culture.color)!`, a `Map` built in `palette.ts` from the same generated table. The baseline compare in Step 7 catches it either way, but knowing which shape you need saves a cycle.

- [ ] **Step 5: Pin the folk table against THREE**

Create `tests/palette.test.ts`:

```ts
import { describe, it, expect } from 'vitest';
import { CULTURE_COLORS, FOLK_COLORS } from '../utils/palette';
import { darkenForFolk } from '../utils/colors';

describe('palette', () => {
  it('FOLK_COLORS matches the THREE-based darkenForFolk exactly', () => {
    expect(FOLK_COLORS).toEqual(CULTURE_COLORS.map(darkenForFolk));
  });

  it('FOLK_COLORS is index-aligned to CULTURE_COLORS', () => {
    expect(FOLK_COLORS.length).toBe(CULTURE_COLORS.length);
  });
});
```

- [ ] **Step 6: Prove the worker closure is free of `three` and `d3`**

```bash
npx madge --json utils/worldGen.ts 2>/dev/null | grep -E '"(three|d3[^"]*)"' && echo "STILL COUPLED" || echo "clean"
```

If `madge` is unavailable, use the direct check instead:

```bash
node -e "
const seen = new Set(); const bad = [];
const fs = require('fs'); const path = require('path');
(function walk(f) {
  if (seen.has(f) || !fs.existsSync(f)) return; seen.add(f);
  for (const m of fs.readFileSync(f, 'utf8').matchAll(/from '([^']+)'/g)) {
    const s = m[1];
    if (!s.startsWith('.')) { if (/^(three|d3)/.test(s)) bad.push(f + ' -> ' + s); continue; }
    const base = path.join(path.dirname(f), s);
    walk(fs.existsSync(base + '.ts') ? base + '.ts' : base + '.tsx');
  }
})('utils/worldGen.ts');
console.log(bad.length ? 'STILL COUPLED:\n' + bad.join('\n') : 'clean');
"
```

Expected: `clean`. `d3-geo-voronoi` is still imported and that is correct — it is the Voronoi builder and must run in the worker. The check above only rejects bare `three` and `d3`.

- [ ] **Step 7: Prove nothing changed**

```bash
node scripts/captureBaseline.mjs after && node scripts/captureBaseline.mjs compare
```

Expected: `IDENTICAL — 0 fields differ`. **If any field differs, this task has a bug — revert and find it.** There is no legitimate reason for a module move to change output.

- [ ] **Step 8: Run all four gates**

```bash
npm run typecheck && npm run lint && npm test && npm run build
```

Expected: 0 errors, ≤29 warnings, 147 tests green, build OK.

- [ ] **Step 9: Commit**

```bash
git add utils/palette.ts utils/vec.ts utils/colors.ts utils/geo.ts utils/features.ts utils/worldGen.ts tests/palette.test.ts
git commit -m "Extract palette and vector leaves out of the worker closure

utils/worldGen.ts transitively imported three (via colors.ts) and all of
d3 (via features.ts -> geo.ts). Neither breaks in a worker, but both get
re-parsed on every spawn, and the next commit spawns one per generation.

darkenForFolk is precomputed into a frozen FOLK_COLORS table rather than
ported, because THREE.Color.setHex applies an sRGB conversion a hand
port would not reproduce; tests/palette.test.ts recomputes it via THREE
and fails if a three upgrade shifts it.

Digest baseline identical across all fields at both 300 and 5000 cells."
```

---

### Task 3: The SoA transfer contract

The load-bearing file. `Cell[]` is an array-of-objects with ragged `vertices` and `neighbors`; structured-cloning 200k of those is the risk the D6 spec names as the biggest in the project. This flattens to typed arrays that cross as zero-copy Transferables.

**Files:**
- Create: `utils/worldTransfer.ts`
- Create: `tests/worldTransfer.test.ts`

**Interfaces:**
- Consumes: `digestWorld`, `diffDigests` from `tests/helpers/worldDigest` (test only).
- Produces:
  - `interface WorldPayload` — the wire format (all fields below).
  - `serializeWorld(world: WorldData): { payload: WorldPayload; transfer: ArrayBuffer[] }`
  - `deserializeWorld(payload: WorldPayload): WorldData`

**Three decisions, with rationale, so they are not re-litigated:**

1. **The main thread rehydrates to `Cell[]` (AoS).** Every consumer — `utils/colors.ts`, `Map2D`, `WorldViewer`, `paintUtils`, `civEdit`, `export` — reads `cell.height`, `cell.neighbors`. Making them read SoA is a rendering-layer rewrite that belongs to F4, not here. Rehydration costs one allocation pass; structured clone costs that *plus* serialization on both sides, so SoA transfer is strictly cheaper even with the shim. **Say out loud that the shim is temporary:** the long-term shape is SoA end-to-end, and this file is where that migration would start.
2. **Optionals need presence bits.** `regionId`, `provinceId`, `population`, `cultureId`, `religionId`, `flux`, `isCapital`, `isTown` are all `?:`, and live code tests `=== undefined`. A sentinel like `-1` is not enough — `undefined` and `-1` must round-trip distinctly, and the Task 1 mutation test for `undefined → 0` exists to enforce that.
3. **Roster-scale data is structured-cloned as-is.** `civData`, `cultures`, `religions`, and the non-`cellIds` parts of `lakes`/`features` number in the hundreds, not the hundred-thousands. Flattening them would be work with no payoff. Only per-cell-scale arrays get the SoA treatment.

- [ ] **Step 1: Write the round-trip test first**

Create `tests/worldTransfer.test.ts`:

```ts
import { describe, it, expect } from 'vitest';
import { generateWorld } from '../utils/worldGen';
import { serializeWorld, deserializeWorld } from '../utils/worldTransfer';
import { digestWorld, diffDigests } from './helpers/worldDigest';
import { makeParams } from './helpers';
import { BiomeType } from '../types';

describe('worldTransfer round trip', () => {
  it('is digest-identical across serialize -> deserialize', async () => {
    const w = await generateWorld(makeParams());
    const { payload } = serializeWorld(w);
    expect(diffDigests(digestWorld(w), digestWorld(deserializeWorld(payload)))).toEqual([]);
  }, 30000);

  it('preserves undefined optionals as undefined, not 0', async () => {
    const w = await generateWorld(makeParams());
    const back = deserializeWorld(serializeWorld(w).payload);
    const water = w.cells.findIndex(c => c.regionId === undefined);
    expect(water).toBeGreaterThanOrEqual(0);
    expect(back.cells[water].regionId).toBeUndefined();
    expect(back.cells[water]).not.toHaveProperty('regionId', 0);
  }, 30000);

  it('preserves a cell with null geoJson geometry', async () => {
    const w = await generateWorld(makeParams());
    w.geoJson.features[0].geometry = null;
    const back = deserializeWorld(serializeWorld(w).payload);
    expect(back.geoJson.features[0].geometry).toBeNull();
  }, 30000);

  it('preserves an empty vertices ring', async () => {
    const w = await generateWorld(makeParams());
    w.cells[3].vertices = [];
    const back = deserializeWorld(serializeWorld(w).payload);
    expect(back.cells[3].vertices).toEqual([]);
  }, 30000);

  it('preserves every biome value through the Uint8 encoding', async () => {
    const w = await generateWorld(makeParams());
    const all = Object.values(BiomeType);
    w.cells.forEach((c, i) => { c.biome = all[i % all.length]; });
    const back = deserializeWorld(serializeWorld(w).payload);
    expect(back.cells.map(c => c.biome)).toEqual(w.cells.map(c => c.biome));
  }, 30000);

  it('lists every buffer it hands out as transferable', async () => {
    const w = await generateWorld(makeParams());
    const { payload, transfer } = serializeWorld(w);
    const inPayload = new Set<ArrayBuffer>();
    for (const v of Object.values(payload as Record<string, unknown>)) {
      if (ArrayBuffer.isView(v)) inPayload.add(v.buffer as ArrayBuffer);
    }
    expect(new Set(transfer)).toEqual(inPayload);
  }, 30000);
});
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `npx vitest run tests/worldTransfer.test.ts`
Expected: FAIL — `Failed to resolve import "../utils/worldTransfer"`.

- [ ] **Step 3: Write the transfer module**

Create `utils/worldTransfer.ts`:

```ts
import { BiomeType, Cell, CivData, CultureData, GeoJsonCollection, LakeData, GeoFeature, MarkerData, Point, ReligionData, RouteData, WorldData, WorldParams } from '../types';

// Stable biome <-> byte mapping. Order is the declaration order of the enum and
// must never be reordered — it is a wire format, not a display list.
const BIOME_LIST = Object.values(BiomeType) as BiomeType[];
const BIOME_INDEX = new Map<BiomeType, number>(BIOME_LIST.map((b, i) => [b, i]));

// Presence bits. Every optional Cell field needs one: live code tests
// `c.regionId === undefined`, so undefined and a sentinel must round-trip
// distinctly. A shared bitfield beats one array per optional.
const P_FLUX = 1, P_REGION = 2, P_PROVINCE = 4, P_POPULATION = 8,
      P_CULTURE = 16, P_RELIGION = 32, P_CAPITAL = 64, P_TOWN = 128;
// Value bits for the two booleans (presence says "is set", this says "is true").
const V_CAPITAL = 1, V_TOWN = 2;

export interface Ragged { offsets: Uint32Array; data: Float64Array }
export interface RaggedI32 { offsets: Uint32Array; data: Int32Array }

export interface WorldPayload {
  cellCount: number;
  // per-cell scalars
  height: Float64Array;
  temperature: Float64Array;
  moisture: Float64Array;
  flux: Float64Array;
  population: Float64Array;
  plateId: Int32Array;
  regionId: Int32Array;
  provinceId: Int32Array;
  cultureId: Int32Array;
  religionId: Int32Array;
  biome: Uint8Array;
  presence: Uint8Array;
  bools: Uint8Array;
  center: Float64Array;        // 3 per cell
  // ragged per-cell
  vertOffsets: Uint32Array; vertData: Float64Array;   // 3 per vertex
  nbrOffsets: Uint32Array;  nbrData: Int32Array;
  ringOffsets: Uint32Array; ringData: Float64Array;   // 2 per coord (lon, lat)
  geomPresent: Uint8Array;                            // geoJson geometry !== null
  // geoJson feature properties. VERIFIED against d3-geo-voronoi v2.1: every
  // feature carries { site: [lon,lat], sitecoordinates: [lon,lat],
  // neighbours: number[] }. It is per-CELL-scale, not roster-scale, so it gets
  // the SoA treatment like everything else at 200k. `neighbours` is d3's own
  // adjacency and is NOT the same array as Cell.neighbors (which is built from
  // links(), deduped, and differently ordered) — transfer both, never alias.
  propsPresent: Uint8Array;
  geoSite: Float64Array;                              // 2 per cell
  geoSiteCoords: Float64Array;                        // 2 per cell
  geoNbrOffsets: Uint32Array; geoNbrData: Int32Array;
  // rivers: array of paths of Points
  riverOffsets: Uint32Array; riverData: Float64Array; // 3 per point
  hasRivers: boolean;
  // ragged cellIds for lakes and features
  lakeIdOffsets: Uint32Array; lakeIdData: Int32Array;
  featIdOffsets: Uint32Array; featIdData: Int32Array;
  // roster-scale: plain structured clone
  params: WorldParams;
  civData?: CivData;
  cultures?: CultureData[];
  religions?: ReligionData[];
  markers?: MarkerData[];
  routes?: RouteData[];
  lakesMeta?: Omit<LakeData, 'cellIds'>[];
  featuresMeta?: Omit<GeoFeature, 'cellIds'>[];
}

const packRagged = <T>(rows: T[][], stride: number, write: (t: T, out: Float64Array, at: number) => void): Ragged => {
  const offsets = new Uint32Array(rows.length + 1);
  let total = 0;
  for (let i = 0; i < rows.length; i++) { offsets[i] = total; total += rows[i].length; }
  offsets[rows.length] = total;
  const data = new Float64Array(total * stride);
  let at = 0;
  for (const row of rows) for (const item of row) { write(item, data, at); at += stride; }
  return { offsets, data };
};

const packRaggedI32 = (rows: number[][]): RaggedI32 => {
  const offsets = new Uint32Array(rows.length + 1);
  let total = 0;
  for (let i = 0; i < rows.length; i++) { offsets[i] = total; total += rows[i].length; }
  offsets[rows.length] = total;
  const data = new Int32Array(total);
  let at = 0;
  for (const row of rows) for (const n of row) data[at++] = n;
  return { offsets, data };
};

const writePoint = (p: Point, out: Float64Array, at: number) => { out[at] = p.x; out[at + 1] = p.y; out[at + 2] = p.z; };
const readPoint = (d: Float64Array, at: number): Point => ({ x: d[at], y: d[at + 1], z: d[at + 2] });

export const serializeWorld = (world: WorldData): { payload: WorldPayload; transfer: ArrayBuffer[] } => {
  const cells = world.cells;
  const n = cells.length;

  const height = new Float64Array(n), temperature = new Float64Array(n), moisture = new Float64Array(n);
  const flux = new Float64Array(n), population = new Float64Array(n);
  const plateId = new Int32Array(n), regionId = new Int32Array(n), provinceId = new Int32Array(n);
  const cultureId = new Int32Array(n), religionId = new Int32Array(n);
  const biome = new Uint8Array(n), presence = new Uint8Array(n), bools = new Uint8Array(n);
  const center = new Float64Array(n * 3);

  for (let i = 0; i < n; i++) {
    const c = cells[i];
    height[i] = c.height; temperature[i] = c.temperature; moisture[i] = c.moisture;
    plateId[i] = c.plateId;
    biome[i] = BIOME_INDEX.get(c.biome)!;
    writePoint(c.center, center, i * 3);
    let p = 0, b = 0;
    if (c.flux !== undefined) { flux[i] = c.flux; p |= P_FLUX; }
    if (c.regionId !== undefined) { regionId[i] = c.regionId; p |= P_REGION; }
    if (c.provinceId !== undefined) { provinceId[i] = c.provinceId; p |= P_PROVINCE; }
    if (c.population !== undefined) { population[i] = c.population; p |= P_POPULATION; }
    if (c.cultureId !== undefined) { cultureId[i] = c.cultureId; p |= P_CULTURE; }
    if (c.religionId !== undefined) { religionId[i] = c.religionId; p |= P_RELIGION; }
    if (c.isCapital !== undefined) { p |= P_CAPITAL; if (c.isCapital) b |= V_CAPITAL; }
    if (c.isTown !== undefined) { p |= P_TOWN; if (c.isTown) b |= V_TOWN; }
    presence[i] = p; bools[i] = b;
  }

  const verts = packRagged(cells.map(c => c.vertices), 3, writePoint);
  const nbrs = packRaggedI32(cells.map(c => c.neighbors));

  const geomPresent = new Uint8Array(n);
  const rings = world.geoJson.features.map((f, i) => {
    if (!f.geometry) return [];
    geomPresent[i] = 1;
    return f.geometry.coordinates[0];
  });
  const ringPacked = packRagged(rings as number[][][], 2, (coord, out, at) => { out[at] = coord[0]; out[at + 1] = coord[1]; });

  const propsPresent = new Uint8Array(n);
  const geoSite = new Float64Array(n * 2);
  const geoSiteCoords = new Float64Array(n * 2);
  const geoNbrRows: number[][] = world.geoJson.features.map((f, i) => {
    const props = f.properties as { site?: number[]; sitecoordinates?: number[]; neighbours?: number[] } | null;
    if (!props) return [];
    propsPresent[i] = 1;
    if (props.site) { geoSite[i * 2] = props.site[0]; geoSite[i * 2 + 1] = props.site[1]; }
    if (props.sitecoordinates) { geoSiteCoords[i * 2] = props.sitecoordinates[0]; geoSiteCoords[i * 2 + 1] = props.sitecoordinates[1]; }
    return props.neighbours ?? [];
  });
  const geoNbrs = packRaggedI32(geoNbrRows);

  const riverPacked = packRagged(world.rivers ?? [], 3, writePoint);
  const lakeIds = packRaggedI32((world.lakes ?? []).map(l => l.cellIds));
  const featIds = packRaggedI32((world.features ?? []).map(f => f.cellIds));

  const payload: WorldPayload = {
    cellCount: n,
    height, temperature, moisture, flux, population,
    plateId, regionId, provinceId, cultureId, religionId,
    biome, presence, bools, center,
    vertOffsets: verts.offsets, vertData: verts.data,
    nbrOffsets: nbrs.offsets, nbrData: nbrs.data,
    ringOffsets: ringPacked.offsets, ringData: ringPacked.data, geomPresent,
    propsPresent, geoSite, geoSiteCoords,
    geoNbrOffsets: geoNbrs.offsets, geoNbrData: geoNbrs.data,
    riverOffsets: riverPacked.offsets, riverData: riverPacked.data,
    hasRivers: world.rivers !== undefined,
    lakeIdOffsets: lakeIds.offsets, lakeIdData: lakeIds.data,
    featIdOffsets: featIds.offsets, featIdData: featIds.data,
    params: world.params,
    civData: world.civData,
    cultures: world.cultures,
    religions: world.religions,
    markers: world.markers,
    routes: world.routes,
    lakesMeta: world.lakes?.map(({ cellIds, ...rest }) => rest),
    featuresMeta: world.features?.map(({ cellIds, ...rest }) => rest),
  };

  const transfer: ArrayBuffer[] = [];
  for (const v of Object.values(payload as Record<string, unknown>)) {
    if (ArrayBuffer.isView(v)) transfer.push(v.buffer as ArrayBuffer);
  }
  return { payload, transfer };
};

export const deserializeWorld = (p: WorldPayload): WorldData => {
  const n = p.cellCount;
  const cells: Cell[] = new Array(n);

  for (let i = 0; i < n; i++) {
    const vs = p.vertOffsets[i], ve = p.vertOffsets[i + 1];
    const vertices: Point[] = new Array(ve - vs);
    for (let k = vs; k < ve; k++) vertices[k - vs] = readPoint(p.vertData, k * 3);

    const ns = p.nbrOffsets[i], ne = p.nbrOffsets[i + 1];
    const neighbors: number[] = new Array(ne - ns);
    for (let k = ns; k < ne; k++) neighbors[k - ns] = p.nbrData[k];

    const pr = p.presence[i], bo = p.bools[i];
    const c: Cell = {
      id: i,
      center: readPoint(p.center, i * 3),
      vertices,
      neighbors,
      height: p.height[i],
      plateId: p.plateId[i],
      temperature: p.temperature[i],
      moisture: p.moisture[i],
      biome: BIOME_LIST[p.biome[i]],
    };
    if (pr & P_FLUX) c.flux = p.flux[i];
    if (pr & P_REGION) c.regionId = p.regionId[i];
    if (pr & P_PROVINCE) c.provinceId = p.provinceId[i];
    if (pr & P_POPULATION) c.population = p.population[i];
    if (pr & P_CULTURE) c.cultureId = p.cultureId[i];
    if (pr & P_RELIGION) c.religionId = p.religionId[i];
    if (pr & P_CAPITAL) c.isCapital = (bo & V_CAPITAL) !== 0;
    if (pr & P_TOWN) c.isTown = (bo & V_TOWN) !== 0;
    cells[i] = c;
  }

  const unpackIds = (offsets: Uint32Array, data: Int32Array, i: number): number[] => {
    const s = offsets[i], e = offsets[i + 1];
    const out: number[] = new Array(e - s);
    for (let k = s; k < e; k++) out[k - s] = data[k];
    return out;
  };

  const geoJson: GeoJsonCollection = {
    type: 'FeatureCollection',
    features: new Array(n).fill(null).map((_, i) => {
      let geometry: { type: 'Polygon'; coordinates: number[][][] } | null = null;
      if (p.geomPresent[i]) {
        const s = p.ringOffsets[i], e = p.ringOffsets[i + 1];
        const ring: number[][] = new Array(e - s);
        for (let k = s; k < e; k++) ring[k - s] = [p.ringData[k * 2], p.ringData[k * 2 + 1]];
        geometry = { type: 'Polygon', coordinates: [ring] };
      }
      const properties = p.propsPresent[i]
        ? {
            site: [p.geoSite[i * 2], p.geoSite[i * 2 + 1]],
            sitecoordinates: [p.geoSiteCoords[i * 2], p.geoSiteCoords[i * 2 + 1]],
            neighbours: unpackIds(p.geoNbrOffsets, p.geoNbrData, i),
          }
        : null;
      return { type: 'Feature' as const, geometry, properties };
    }),
  };

  const rivers: Point[][] | undefined = p.hasRivers
    ? new Array(p.riverOffsets.length - 1).fill(null).map((_, i) => {
        const s = p.riverOffsets[i], e = p.riverOffsets[i + 1];
        const path: Point[] = new Array(e - s);
        for (let k = s; k < e; k++) path[k - s] = readPoint(p.riverData, k * 3);
        return path;
      })
    : undefined;

  const lakes: LakeData[] | undefined = p.lakesMeta?.map((m, i) => ({ ...m, cellIds: unpackIds(p.lakeIdOffsets, p.lakeIdData, i) }));
  const features: GeoFeature[] | undefined = p.featuresMeta?.map((m, i) => ({ ...m, cellIds: unpackIds(p.featIdOffsets, p.featIdData, i) }));

  const world: WorldData = { cells, params: p.params, geoJson };
  if (p.civData !== undefined) world.civData = p.civData;
  if (rivers !== undefined) world.rivers = rivers;
  if (lakes !== undefined) world.lakes = lakes;
  if (features !== undefined) world.features = features;
  if (p.cultures !== undefined) world.cultures = p.cultures;
  if (p.religions !== undefined) world.religions = p.religions;
  if (p.markers !== undefined) world.markers = p.markers;
  if (p.routes !== undefined) world.routes = p.routes;
  return world;
};
```

**Note on `properties` — already verified, do not re-derive.** An early draft of this plan reconstructed `properties: null`. That is wrong: `d3-geo-voronoi` v2.1's `polygons()` attaches `{ site: [lon, lat], sitecoordinates: [lon, lat], neighbours: number[] }` to **every** feature. Probed directly on 2026-07-27, output verbatim:

```json
{ "site": [-180, -80], "sitecoordinates": [-180, -80],
  "neighbours": [19, 10, 1, 28, 31, 22, 13, 16, 37, 25, 34] }
```

Two consequences the code above already handles: (a) `neighbours` is **per-cell-scale** ragged data, so it cannot go in the roster-scale structured clone — at 200k cells that is 200k objects, exactly what this file exists to avoid; (b) `properties.neighbours` is d3's own adjacency and is **not** the same array as `Cell.neighbors`, which `worldGen.ts` builds separately from `voronoi.links()`, dedupes through a `Set`, and orders differently. Transfer both. Never alias one to the other to "save space".

- [ ] **Step 4: Run the tests to verify they pass**

Run: `npx vitest run tests/worldTransfer.test.ts`
Expected: 6 passed. If "digest-identical" fails, `diffDigests` names the exact field — fix that field's encoding, do not loosen the test.

- [ ] **Step 5: Run all four gates**

```bash
npm run typecheck && npm run lint && npm test && npm run build
```

Expected: 0 errors, ≤29 warnings, 153 tests green, build OK.

- [ ] **Step 6: Commit**

```bash
git add utils/worldTransfer.ts tests/worldTransfer.test.ts
git commit -m "Add SoA world transfer contract

Cell is an array-of-objects with ragged vertices and neighbors;
structured-cloning 200k of them across the worker boundary is the
biggest risk the D6 spec names. This flattens per-cell data to typed
arrays transferred zero-copy, and rehydrates to the identical Cell[] on
the main thread.

Optional fields carry presence bits rather than sentinels, because live
code tests === undefined and undefined must not collapse to 0 or -1.
Roster-scale data (civData, cultures, religions) is structured-cloned
as-is — hundreds of objects, not hundreds of thousands.

The AoS rehydration is a deliberate temporary shim: every renderer reads
cell.height today, and making them read SoA is F4 work, not this stage."
```

---

### Task 4: The worker and its client

**Files:**
- Create: `workers/worldGen.worker.ts`, `utils/worldGenClient.ts`
- Modify: `vite.config.ts`

**Interfaces:**
- Consumes: `serializeWorld`, `deserializeWorld`, `WorldPayload` from `utils/worldTransfer`.
- Produces: `generateWorldInWorker(params: WorldParams, onLog?: (msg: string) => void, signal?: AbortSignal, onProgress?: (stage: number, total: number) => void): Promise<WorldData>` — **identical signature to `generateWorld`**, so the call sites in Task 6 change by one identifier.

**Abort is `terminate()`, not a message. This is the decision that makes the yields deletable.** A worker running a synchronous loop cannot receive `postMessage` — message events are macrotasks, so an abort message would only be seen at a yield, which is exactly what this stage removes. `SharedArrayBuffer` + `Atomics` would work but requires COOP/COEP headers, a deployment constraint on Netlify for zero benefit here. So: **one worker per generation, terminated on abort.** Worker spawn is ~1–5ms against a multi-second generation, and Task 2 already stripped `three` and `d3` from the bundle it has to parse.

**Consequence to accept:** `AbortController` semantics change subtly. Today an aborted generation *throws* `"Generation Cancelled"` from inside `generateWorld`. With `terminate()` the worker simply stops, and the client rejects with the same message from the main-thread side. `useWorldEngine` only inspects `e.message === "Generation Cancelled"`, so behaviour is preserved — but the error no longer originates in `checkAbort`.

- [ ] **Step 1: Enable ES-module workers in Vite**

In `vite.config.ts`, add `worker: { format: 'es' }` to the returned config object, as a sibling of `plugins`:

```ts
      plugins: [react()],
      worker: {
        // ES-module worker: the generation closure uses ESM imports
        // (d3-geo-voronoi et al.) and the classic-worker format cannot.
        format: 'es',
      },
```

- [ ] **Step 2: Write the worker entry**

Create `workers/worldGen.worker.ts`:

```ts
import { generateWorld } from '../utils/worldGen';
import { serializeWorld, WorldPayload } from '../utils/worldTransfer';
import { WorldParams } from '../types';

export type WorkerRequest = { type: 'generate'; params: WorldParams };
export type WorkerResponse =
  | { type: 'log'; msg: string }
  | { type: 'progress'; stage: number; total: number }
  | { type: 'done'; payload: WorldPayload }
  | { type: 'error'; message: string };

// This project's tsconfig includes the DOM lib, so bare `self` types as Window
// and `self.onmessage` will not accept a MessageEvent<WorkerRequest> handler.
// Adding `/// <reference lib="webworker" />` collides with DOM in the same
// program — the classic time-sink. One local alias avoids both.
const ctx = self as unknown as Worker;

// Abort is by terminate() from the client, never by message: a synchronous
// generation loop cannot drain the message queue, and draining it is exactly
// what deleting the setTimeout(0) yields gives up. See the client for why.
ctx.onmessage = async (e: MessageEvent<WorkerRequest>) => {
  if (e.data.type !== 'generate') return;
  try {
    const world = await generateWorld(
      e.data.params,
      msg => ctx.postMessage({ type: 'log', msg } satisfies WorkerResponse),
      undefined,
      (stage, total) => ctx.postMessage({ type: 'progress', stage, total } satisfies WorkerResponse),
    );
    const { payload, transfer } = serializeWorld(world);
    ctx.postMessage({ type: 'done', payload } satisfies WorkerResponse, transfer);
  } catch (err) {
    ctx.postMessage({ type: 'error', message: (err as Error).message } satisfies WorkerResponse);
  }
};
```

- [ ] **Step 3: Write the client**

Create `utils/worldGenClient.ts`:

```ts
import { WorldData, WorldParams } from '../types';
import { deserializeWorld } from './worldTransfer';
import type { WorkerRequest, WorkerResponse } from '../workers/worldGen.worker';
import WorldGenWorker from '../workers/worldGen.worker?worker';

// Drop-in replacement for generateWorld with the identical signature, so the
// call sites in useWorldEngine change by one identifier.
//
// One worker per generation, terminated on abort. A long-lived worker would
// need a message-based abort, and a synchronous generation loop cannot receive
// messages — which is the whole reason the setTimeout(0) yields can go away.
// Spawn cost is ~1-5ms against a multi-second generation.
export const generateWorldInWorker = (
  params: WorldParams,
  onLog?: (msg: string) => void,
  signal?: AbortSignal,
  onProgress?: (stage: number, total: number) => void,
): Promise<WorldData> =>
  new Promise<WorldData>((resolve, reject) => {
    if (signal?.aborted) { reject(new Error('Generation Cancelled')); return; }

    const worker = new WorldGenWorker();
    let settled = false;
    const finish = (fn: () => void) => {
      if (settled) return;
      settled = true;
      signal?.removeEventListener('abort', onAbort);
      worker.terminate();
      fn();
    };
    const onAbort = () => finish(() => reject(new Error('Generation Cancelled')));
    signal?.addEventListener('abort', onAbort);

    worker.onmessage = (e: MessageEvent<WorkerResponse>) => {
      const m = e.data;
      if (m.type === 'log') onLog?.(m.msg);
      else if (m.type === 'progress') onProgress?.(m.stage, m.total);
      else if (m.type === 'done') finish(() => resolve(deserializeWorld(m.payload)));
      else if (m.type === 'error') finish(() => reject(new Error(m.message)));
    };
    worker.onerror = (e: ErrorEvent) => finish(() => reject(new Error(e.message || 'Worker failed')));

    worker.postMessage({ type: 'generate', params } satisfies WorkerRequest);
  });
```

- [ ] **Step 4: Add the worker type declaration if `?worker` does not resolve**

If `npm run typecheck` reports `Cannot find module '../workers/worldGen.worker?worker'`, create `vite-env.d.ts` at the repo root (or append to the existing one):

```ts
/// <reference types="vite/client" />
```

That pulls in Vite's own declarations for `?worker`. Do not hand-write a module shim.

- [ ] **Step 5: Run all four gates**

```bash
npm run typecheck && npm run lint && npm test && npm run build
```

Expected: 0 errors, ≤29 warnings, 153 tests green, build OK. **Check the build output for a separate worker chunk** — you should see a `worldGen.worker-<hash>.js` asset. If the worker code got inlined into the main bundle instead, `format: 'es'` did not apply.

- [ ] **Step 6: Commit**

```bash
git add workers/worldGen.worker.ts utils/worldGenClient.ts vite.config.ts vite-env.d.ts
git commit -m "Add generation worker and drop-in client

generateWorldInWorker has the identical signature to generateWorld, so
switching the call sites is a one-identifier change.

Abort is worker.terminate(), not a message. A synchronous generation
loop cannot drain the message queue, so a message-based abort would only
be seen at a yield — and deleting those yields is the point of the
migration. SharedArrayBuffer + Atomics would work but needs COOP/COEP
headers on Netlify for no benefit. One worker per generation instead;
spawn is ~1-5ms against a multi-second run.

Not yet wired: useWorldEngine still calls generateWorld directly."
```

---

### Task 5: The gate — main thread vs worker, in one browser session

This is the gate the whole stage is sequenced around. It runs in a real browser because that is the only place a real `Worker`, a real structured clone, and real Transferables exist. **It is not in CI**, and that limit is stated on the page itself.

**Files:**
- Create: `dev/goldenCompare.html`

**Interfaces:**
- Consumes: `generateWorld` (`utils/worldGen`), `generateWorldInWorker` (`utils/worldGenClient`).
- Produces: nothing importable. A page that prints PASS/FAIL and the differing field names.

- [ ] **Step 1: Write the harness page**

Create `dev/goldenCompare.html`. Note the digest is re-implemented inline against `crypto.subtle` — `tests/helpers/worldDigest.ts` uses `node:crypto` and cannot load in a browser. **This is a feature, not duplication to be DRY'd away:** an independent second implementation of the instrument is exactly the "instrument that does not share the migration's assumptions" rule from Session 6e.

```html
<meta charset="utf-8">
<title>D6 Stage 1 — main vs worker</title>
<style>body{font:13px ui-monospace,monospace;background:#0a0a0a;color:#e5e5e5;padding:2rem;line-height:1.6}
.pass{color:#4ade80}.fail{color:#f87171}h1{font-size:15px}</style>
<h1>D6 Stage 1 — main thread vs worker</h1>
<p>Not a CI gate: needs a real Worker, real structured clone, real Transferables.
Run it by hand before wiring the client into useWorldEngine, and again after.</p>
<pre id="out">running…</pre>
<script type="module">
import { generateWorld } from '/utils/worldGen.ts';
import { generateWorldInWorker } from '/utils/worldGenClient.ts';

const out = document.getElementById('out');
const log = (s, cls = '') => { out.innerHTML += `\n<span class="${cls}">${s}</span>`; };

// Deliberately a SECOND implementation of the digest, not an import of
// tests/helpers/worldDigest.ts (which needs node:crypto). An instrument that
// shares the migration's assumptions cannot detect the migration's mistakes.
const dv = new DataView(new ArrayBuffer(8));
const enc = v => {
  if (v === undefined) return 'u';
  if (v === null) return 'n';
  if (typeof v === 'number') { if (Object.is(v, -0)) return '-0'; dv.setFloat64(0, v); return dv.getBigUint64(0).toString(16); }
  if (typeof v === 'boolean') return v ? 'T' : 'F';
  if (typeof v === 'string') return 's' + v;
  if (Array.isArray(v)) return '[' + v.map(enc).join(',') + ']';
  if (typeof v === 'object') return '{' + Object.keys(v).sort().map(k => k + ':' + enc(v[k])).join(',') + '}';
  return 'x' + String(v);
};
const sha = async s => {
  const b = await crypto.subtle.digest('SHA-256', new TextEncoder().encode(s));
  return [...new Uint8Array(b)].map(x => x.toString(16).padStart(2, '0')).join('').slice(0, 32);
};
const FIELDS = ['id','height','plateId','temperature','moisture','biome','flux','regionId','provinceId','isCapital','isTown','population','cultureId','religionId'];
const digest = async w => {
  const d = { cellCount: String(w.cells.length) };
  for (const f of FIELDS) d['cell.' + f] = await sha(w.cells.map(c => enc(c[f])).join(';'));
  d['cell.center'] = await sha(w.cells.map(c => enc(c.center)).join(';'));
  d['cell.vertices'] = await sha(w.cells.map(c => enc(c.vertices)).join(';'));
  d['cell.neighbors'] = await sha(w.cells.map(c => enc(c.neighbors)).join(';'));
  for (const k of ['geoJson','params','rivers','lakes','features','cultures','religions','civData','markers','routes'])
    d[k] = await sha(enc(w[k]));
  return d;
};

const params = {
  mapName: 'golden', points: 5000, seed: 'golden_seed', planetRadius: 6371, axialTilt: 23.5,
  landStyle: 'Continents', cellJitter: 0.5, noiseScale: 0.4, ridgeBlend: 0.1,
  mountainHeight: 1.0, oceanDepth: 1.0, maskType: 'None', warpStrength: 0.5,
  plateInfluence: 0.5, erosionIterations: 3, plates: 8, seaLevel: 0.55, roughness: 0.5,
  detailLevel: 3, baseTemperature: 30, poleTemperature: -30, rainfallMultiplier: 1.0,
  moistureTransport: 0.5, temperatureVariance: 5, numFactions: 4, civSeed: 'golden_civs',
  borderRoughness: 0.2, civSizeVariance: 0.5, waterCrossingCost: 0.8, territorialWaters: 0.15,
  capitalSpacing: 0.5, provinceSize: 0.5, numCultures: 4, nameStyle: 'fantasy', loreLevel: 1,
};

out.textContent = 'generating on the main thread…';
const t0 = performance.now();
const a = await generateWorld(params);
const tMain = performance.now() - t0;

out.textContent = 'generating in the worker…';
const t1 = performance.now();
const b = await generateWorldInWorker(params);
const tWorker = performance.now() - t1;

out.textContent = 'digesting…';
const [da, db] = [await digest(a), await digest(b)];
const keys = [...new Set([...Object.keys(da), ...Object.keys(db)])].sort();
const diffs = keys.filter(k => da[k] !== db[k]);

out.textContent = '';
log(`cells: ${a.cells.length}   main: ${tMain.toFixed(0)}ms   worker: ${tWorker.toFixed(0)}ms`);
log(`fields compared: ${keys.length}`);
if (diffs.length === 0) log('PASS — worker output is bit-identical to the main thread', 'pass');
else { log(`FAIL — ${diffs.length} field(s) differ:`, 'fail'); diffs.forEach(k => log('  ' + k, 'fail')); }
</script>
```

- [ ] **Step 2: Run it**

```bash
npm run dev -- --port 4180
```

Then open `http://localhost:4180/dev/goldenCompare.html`. **Do not use port 3000 — that server may be Matt's.**

Expected on the page: `PASS — worker output is bit-identical to the main thread`, over ~28 fields at 5000 cells.

- [ ] **Step 3: If it fails, read the field names before touching anything**

The failing field names localize the bug precisely:
- `cell.flux` / `cell.population` alone → a presence-bit case in `worldTransfer.ts`.
- `cell.regionId` + `cell.provinceId` + `cell.cultureId` together → the whole optional-int path.
- `geoJson` alone → the `properties` reconstruction warned about in Task 3 Step 3.
- `cell.vertices` or `cell.neighbors` → a ragged offset off by one.
- *Everything* differs → not serialization. Check whether the worker is somehow seeing different params, or `Math` differs — but it should not; same engine, same thread pool.

- [ ] **Step 4: Widen the identity check to 20000 cells**

Edit `points: 5000` → `points: 20000` and reload. Expected: still `PASS`.

**Stop there for identity. Do not run the digest at 200k.** Both digest implementations build one giant string per field before hashing — `w.cells.map(c => enc(c.vertices)).join(';')` at 200k cells is roughly 200000 × 7 vertices × 3 floats × 17 hex chars ≈ 70MB in a single string, and `enc(world.geoJson)` is worse. That OOMs or stalls the tab for minutes, and the failure reads as *"the transfer can't handle 200k"* when the truth is *"the instrument can't"*. Misdiagnosing that would be the expensive mistake of this task.

- [ ] **Step 5: Measure timing and survival at 200000 — no digest**

Add a URL-driven fast path to the page so the cap can be measured without hashing. At the top of the script:

```js
const CELLS = Number(new URLSearchParams(location.search).get('points') ?? 5000);
const TIMING_ONLY = new URLSearchParams(location.search).has('timing');
```

Use `points: CELLS` in the params object, and guard the digest block:

```js
if (TIMING_ONLY) {
  log(`cells: ${a.cells.length}   main: ${tMain.toFixed(0)}ms   worker: ${tWorker.toFixed(0)}ms`);
  log(`worker cells rehydrated: ${b.cells.length}, geoJson features: ${b.geoJson.features.length}`);
  log('TIMING ONLY — identity not checked at this size (see plan Task 5 Step 4)', 'fail');
} else { /* the existing digest + diff block */ }
```

Then open `http://localhost:4180/dev/goldenCompare.html?points=200000&timing`.

Record: main-thread ms, worker ms, and whether the run completed at all. This is the first real measurement of transfer cost at the cap and Stage 2 will want it.

**If the tab dies at 200k**, that is a finding, not a reason to skip: record the failure mode (OOM during rehydration vs `postMessage` failure vs a stall that never returns) in `HANDOFF.md` in Task 7, because it directly bounds what Stage 2 can assume about display resolution.

- [ ] **Step 6: Commit**

```bash
git add dev/goldenCompare.html
git commit -m "Add main-vs-worker bit-identity harness

Runs both paths in one browser session and diffs a full-surface digest.
Not in CI: a real Worker, real structured clone and real Transferables
only exist in a browser, and a committed golden fixture was rejected
because Math.sin/cos/pow are implementation-defined and would drift a
last-ULP across V8 versions until someone loosened it into uselessness.

The digest here is a deliberate second implementation rather than an
import of tests/helpers/worldDigest.ts — an instrument that shares the
migration's assumptions cannot catch the migration's mistakes.

Identity verified at 5000 and 20000 cells; deliberately NOT at 200000,
where the digest's own per-field string concatenation is a ~70MB string
that OOMs the tab and would misread as a transfer failure. The cap is
measured for timing and survival only.

Measured: 5000 cells main <N>ms / worker <N>ms, 0 of <N> fields differ.
200000 cells main <N>ms / worker <N>ms, <completed | failed how>."
```

---

### Task 6: Wire it up and delete the yields

**Files:**
- Modify: `hooks/useWorldEngine.ts:187`, `hooks/useWorldEngine.ts:228`, `hooks/useWorldEngine.ts:3`
- Modify: `utils/worldGen.ts` — remove `await new Promise(r => setTimeout(r, 0))` at lines 118, 156, 189, 221, 305, 343, 559, 731 and `setTimeout(r, 10)` at 520; remove the now-dead `checkAbort` calls that followed them

- [ ] **Step 1: Capture the before-baseline**

```bash
node scripts/captureBaseline.mjs before
```

- [ ] **Step 2: Switch the two call sites**

In `hooks/useWorldEngine.ts`:
- Line 3: `import { generateWorld, recalculateCivs, recalculateProvinces } from '../utils/worldGen';` → `import { recalculateCivs, recalculateProvinces } from '../utils/worldGen';` plus a new line `import { generateWorldInWorker } from '../utils/worldGenClient';`
- Line 187: `await generateWorld(p, …)` → `await generateWorldInWorker(p, …)`
- Line 228: `await generateWorld(newParams, …)` → `await generateWorldInWorker(newParams, …)`

Nothing else in the hook changes. `recalculateCivs` and `recalculateProvinces` stay on the main thread — they run on an already-materialized world in response to slider changes, and moving them is a separate decision with its own transfer cost.

**Watch for this, and note it in `HANDOFF.md` if it bites.** Nothing under `tests/` imports `utils/worldGenClient.ts` today, so the `?worker` import suffix has never been through Vitest's transform. This step makes `useWorldEngine.ts` import it — so **the first test that ever imports the hook will fail on an unresolvable `?worker` module**, and the error will look nothing like a worker problem. No such test exists right now (the suite covers the pure engine only), so this is a tripwire for a future session, not a blocker here. The fix when it happens is `@vitest/web-worker` or a `vi.mock` of the client, not deleting the import.

- [ ] **Step 3: Delete the yields**

In `utils/worldGen.ts`, remove every `await new Promise(r => setTimeout(r, 0));` (and the one `setTimeout(r, 10)`), and the `checkAbort(signal);` line that immediately follows each. Keep the `signal` parameters and the `checkAbort` helper: `generateWorld` is still directly callable (tests, the compare harness), and `tests/worldGen.test.ts` asserts that an already-aborted signal throws before any work happens.

**Verify that specific test still passes** — it aborts *before* calling, so the entry-point `checkAbort` must survive. If no `checkAbort` remains at the top of `generateWorld`, add one there explicitly rather than restoring a yield.

- [ ] **Step 4: Prove the algorithm still produces identical values**

```bash
node scripts/captureBaseline.mjs after && node scripts/captureBaseline.mjs compare
```

Expected: `IDENTICAL — 0 fields differ`. Removing a `setTimeout` cannot change output — it consumes no RNG — so any difference here means a `checkAbort` deletion took a line of real logic with it.

- [ ] **Step 5: Re-run the browser gate**

Restart the dev server and reload `http://localhost:4180/dev/goldenCompare.html`.
Expected: still `PASS`, and the worker timing should now be **lower** than in Task 5 — the yields are gone. Note the delta.

- [ ] **Step 6: Interactive smoke on the real app**

The app itself has no automated coverage of these paths (HANDOFF, Session 5 advisor note 3). Run each by hand at `http://localhost:4180/?shell=1`:

- [ ] Generate a world — progress bar advances and reaches 100%, log lines stream in, globe renders.
- [ ] **Drag a slider mid-generation to trigger a regenerate** — the first run must abort cleanly and the second must complete. This is the path `abortControllerRef` guards, and it is now `terminate()`-based; it has never had coverage.
- [ ] Paint a terrain stroke, then undo. Confirms the rehydrated `cells` array identity works with the paint model (`WorldMesh` geometry is keyed on `world.cells` identity — CLAUDE.md invariant).
- [ ] Inspect a cell — Inspector shows live height/biome/temperature/moisture, not `undefined`.
- [ ] Switch to 2D Mercator and Dymaxion — both render. Confirms the rehydrated `geoJson` is what d3 expects.
- [ ] Save to a slot, reload the page, load it back. Exercises the `handleLoadWorld` worker call site.
- [ ] **Confirm the UI stays responsive during generation** — drag the globe while it generates. This is the actual user-facing payoff of the whole stage; if it still stutters, something is still running on the main thread.

- [ ] **Step 7: Run all four gates**

```bash
npm run typecheck && npm run lint && npm test && npm run build
```

Expected: 0 errors, ≤29 warnings, 153 tests green, build OK.

- [ ] **Step 8: Commit**

```bash
git add hooks/useWorldEngine.ts utils/worldGen.ts
git commit -m "Run generation in the worker and delete the yield staging

useWorldEngine's two generateWorld call sites become
generateWorldInWorker; the signatures are identical so nothing else in
the hook moves. recalculateCivs/recalculateProvinces stay on the main
thread — they run on an already-materialized world and moving them is a
separate call with its own transfer cost.

The nine setTimeout(0) yields existed only to keep the main thread
barely responsive between stages. Off-thread they cost time and buy
nothing, so they are gone along with the checkAbort calls that rode on
them; the entry-point checkAbort stays because generateWorld is still
directly callable and the abort test asserts it.

Digest baseline identical; browser harness still bit-identical; the
interactive smoke covered generate, abort-mid-generate, paint+undo,
inspect, both 2D projections, and save/load."
```

---

### Task 7: Record what was learned

**Files:**
- Modify: `HANDOFF.md`

- [ ] **Step 1: Write the session entry**

Add a new session section at the top of `HANDOFF.md` (below Matt's notes, above the current pickup block), and **replace the "NEW-THREAD PICKUP" block** so the next thread is pointed at Stage 2 rather than at this plan.

The entry must carry, at minimum:

- **The instrument limits table** from this plan's preamble — which of the three instruments catches what, and the blind spot of each. This is the part most likely to be re-derived wrongly.
- **Why there is no committed golden fixture** (implementation-defined `Math.sin`/`cos`/`pow`, cross-engine ULP drift, the loosening death spiral). Someone *will* propose adding one.
- **Why abort is `terminate()`** and not a message — with the one-line reason (a synchronous loop cannot drain the message queue), because "just post an abort message" is the obvious-looking thing to try.
- **The measured numbers**: main-thread vs worker generation time at 5000 and 200000 cells, and whether 200k survived at all.
- **The AoS rehydration shim** is temporary and its removal is F4 work, not Stage 2's.
- Anything that failed or surprised, at the confidence level you actually have — say "n=1, unconfirmed" out loud where that is the truth.

Do **not** narrate progress ("moved the palette, then wrote the worker"). That is what `git log` is for. The test is: *does the next session need this to avoid repeating a mistake or re-litigating a decision?*

- [ ] **Step 2: Point the pickup at Stage 2**

The next block of work is D6 Stage 2 (the V3 terrain model). It has four unresolved questions in §9 of `docs/superpowers/specs/2026-07-26-d6-terrain-v3-design.md` that Stage 1 does not answer, and it starts with `brainstorming` or `writing-plans` against §3–§5 of that spec — not with code. Say so explicitly, and repeat the §5.1 refuted-hypothesis warning (accumulating uplift over timesteps is **not** the seam fix and makes the wall thinner and taller) so a fresh thread meets it without opening the spec.

- [ ] **Step 3: Commit**

```bash
git add HANDOFF.md
git commit -m "Record Session 7: D6 Stage 1 worker migration"
```

---

## Self-review notes

**Spec coverage.** This plan implements exactly D6 spec §6 Stage 1 and nothing else. Stage 1's stated scope — "move the existing pipeline into a Web Worker with no algorithm change", gated on byte-identical output, with the Voronoi build moving in too and returning flattened typed arrays — is covered by Tasks 2–6. The spec's §7 determinism rules are honored trivially: this stage adds no RNG draws and no `WorldParams` keys, so `tests/paramLiveness.test.ts` is untouched. §6's claim that "the determinism suite already tests exactly this" is **contradicted by this plan on purpose** — the suite compares two in-process runs over four rounded fields, which is why Tasks 1 and 5 exist. That correction should propagate into `HANDOFF.md` in Task 7.

**Deliberately deferred to Stage 2+:** everything in spec §3, §4, §5 (crust/plate decoupling, Euler poles, coarse simulation, seam elimination), all of §9's open questions, and Stage 3's continuous height field.

**Known soft spots in this plan, stated rather than hidden:**

- **Task 5 is not a CI gate.** It requires a human to open a page. The mitigation is that it is the *only* step that can prove main-thread/worker equivalence at all, and it is run twice (Tasks 5 and 6).
- **Task 2 Step 4's `FOLK_COLORS[culture.id % …]` indexing is an assumption** with a stated fallback and a gate that catches it either way.
- **The 200k-cell transfer is unmeasured**, and Task 5 Step 5 measures timing and survival only — bit-identity at the cap is out of reach because the digest itself does not scale there. If the cap fails, that is a finding that constrains Stage 2's display resolution, not a reason to abandon the stage: the worker still wins at every count below it.

**Resolved during review, recorded so it is not re-derived:** an early draft reconstructed `geoJson.features[].properties` as `null`. Probing `d3-geo-voronoi` v2.1 directly showed it always carries `{ site, sitecoordinates, neighbours }`, and that `neighbours` is per-cell-scale — so it needed the SoA treatment, not the roster-scale clone. Task 3 now handles it, and the note there carries the verbatim probe output.
