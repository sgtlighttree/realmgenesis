# Testing

The strategy in one line: **test the seeded pure engine exhaustively; verify the
canvas/WebGL layer by hand.** The generation core is DOM-free and deterministic,
so tests are cheap and stable there; R3F/Canvas testing has poor cost/benefit,
so rendering is checked in the browser.

## What's covered

Vitest suite in `tests/` (25 files, run `npm test`; single file
`npx vitest run tests/<name>.test.ts`). It covers the pure engine: RNG,
biome classification, generation determinism, param liveness, paint utils,
lakes, cultures, religions, routes, markers, features, measurement, namegen,
pathfinding, plate connectivity, worldGen client/transfer, vector export, and
import validation.

`simulationResolution` in `tests/helpers.ts` is **10000** — production value, so
every test file runs a full 10k-macro V3 simulation. `vite.config.ts` sets
`testTimeout: 120000` because V3 generation is ~9 s vs the old ~0.3 s. Consequence:
`npm test` takes ~3 minutes. If you need it faster, lower
`simulationResolution` in the helper — but see the flake note before assuming a
timeout is a real failure.

## The param-liveness contract

`tests/paramLiveness.test.ts` fails if any tunable `WorldParams` key stops
influencing generated output. **Extend it whenever you add a param.** Two
patterns matter:

- **Terrain params** must change the `terrainSignature`. Add the key to the V3
  perturbation map with a value ≠ its default.
- **A param can read "dead" at the default test world but be live elsewhere.**
  `provinceSize` needs enough faction density to subdivide; `capitalSpacing`
  only binds when capitals are dense enough for the min-separation to reject a
  candidate. These get **dedicated binding-density cases** (higher point count /
  faction count), not the generic loop. When a new param looks dead, check
  whether the default world simply can't exercise it before assuming a bug.

This test is the safety net that catches a slider wired to nothing (a real past
defect: `civSizeVariance`/`detailLevel` were once live UI bound to unread
params).

## Determinism: three instruments, three blind spots — do not collapse them

The Web-Worker migration is guarded by three separate checks. Each is blind to
something the others catch; keeping all three is deliberate.

| Instrument | Catches | Blind to |
|-----------|---------|----------|
| `tests/worldGen.test.ts` | run-to-run nondeterminism, in-process | anything applied to both runs; a few fields at 6 decimals |
| `tests/helpers/worldDigest.ts` + `scripts/captureBaseline.mjs` | a refactor changing values on this engine, all fields at exact IEEE-754 bits | cross-engine ULP drift; can't compare main thread to worker |
| `dev/goldenCompare.html` | serialization loss across the real worker boundary | whether the algorithm is *right* — only that both paths agree |

`tests/worldGen.test.ts` alone would green-light a broken migration: it compares
two in-process runs of the *same code* at `toFixed(6)`, which passes a
`Float64`→`Float32` downcast, a dropped field, and an `undefined`→`0` collapse.
That's why the digest (exact bits, all fields) and the golden-compare (real
`postMessage`) exist alongside it.

### Why there is NO committed golden fixture

Someone will propose adding a checked-in bit-exact baseline. Don't.
`Math.sin`/`cos`/`pow` are **implementation-defined** in ECMAScript, so a
committed fixture drifts a last-ULP across V8 versions → becomes a CI flake →
gets `toBeCloseTo`'d → and at that point no longer catches the `Float64`→`Float32`
downcast it existed for. Baselines are captured **same-engine, same-session**
into gitignored `tmp/` instead — no drift term.

**Trap in `captureBaseline.mjs`:** it stamps `gitSha` = HEAD, which does *not*
capture working-tree state. A `before` captured *after* editing looks identical
to a correctly-sequenced one and proves nothing. Capture `before` from a pristine
`git worktree add --detach <pre-change-sha>` (symlink `node_modules` in).

## Known flake: paramLiveness timeout under M1 parallel load

`npm test` runs 25 files in parallel, each launching a 10k-macro V3 sim. On the
M1 this occasionally throws **one** `paramLiveness` *timeout* (120 s) — not a
dead-param assertion. It passes in isolation
(`npx vitest run tests/paramLiveness.test.ts`) and did not recur on re-run. It is
**not a real failure and not a dead param.** If CI flakes here, the fix is to
lower `simulationResolution` in `tests/helpers.ts` or cap Vitest concurrency —
not to weaken an assertion.

## Rendering is verified by hand

There are no component/WebGL tests — an explicit trade-off (canvas/R3F testing is
brittle and expensive here). Render behavior is checked in the browser. Two
method traps worth knowing:

- **Synthetic `MouseEvent`s do not drive `Map2D` painting** — only real events
  do. Don't conclude paint is broken from a synthetic-canvas pixel-readback;
  read the app's own undo counter instead.
- **Playwright's YAML aria-snapshot elides the name of a button whose text sits
  in a nested `<div>`** — a snapshot-formatting artifact, not a real a11y defect.
  Confirm against the DOM (`getByRole('button', {name})`) before chasing one.
