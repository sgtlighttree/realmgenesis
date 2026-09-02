# Performance workflow (F4)

How optimisation work is done in this repo. F4 is not a milestone that finishes —
it runs alongside the rest of the roadmap — so it needs a *method* rather than a
plan, and the method is what this document is.

Written after the session-28 attempt, where four agents were dispatched in
parallel to "investigate and optimise" four subsystems and produced **zero
optimisations between them**. What survived was their instruments. That outcome
is the argument for everything below.

---

## The one rule

**No optimisation lands without a before number, an after number, and a
correctness number.**

The third is the one people skip, and it is the one that matters. A faster
routine that changes the output is not a fast routine, it is a bug with good
timing. Every harness here therefore emits a correctness signal alongside the
timing:

| Layer | Timing signal | Correctness signal |
|---|---|---|
| Generation | stage ms | **byte-identical world digest**, or index-for-index agreement with the routine being replaced |
| 3D globe | ms/frame per tenant | **max screen error (px)** and **visibility flips** vs the baseline maths |
| 2D map | ms per redraw / per gesture frame | **`settled.hash`** of the rendered pixels |
| React | render counts per interaction | the app still behaves — exercised, not assumed |

A frame-time graph cannot show you that a projection moved a point half a pixel
or flipped one cell across the horizon. Ask for the correctness number first,
then the speed.

---

## The instruments

They already exist. Writing them is the expensive part of perf work; do not
rebuild them.

### `scripts/renderGlobeBench.mjs` → `scripts/perf/globeBench.ts`

Per-tenant `ScreenOverlay` cost, geometry-refill cost, and candidate projector
variants.

```bash
node scripts/renderGlobeBench.mjs --points=30000 --frames=120 --out=tmp/before.txt
# make the change
node scripts/renderGlobeBench.mjs --points=30000 --frames=120 --out=tmp/after.txt
diff tmp/before.txt tmp/after.txt
```

Tenants are drawn against a **null 2D context**, so each number is the JS cost —
projection, culling, polyline building — with rasterisation removed. That is the
part an optimisation can actually move; the raster cost belongs to the browser.

It is a **script and not a test**, deliberately. It generates a 30k world and can
measure for minutes; `npm test` must never pick it up. It arrived as
`tests/f4bench.test.ts` and was converted for exactly that reason.

### `scripts/renderMap2DPerf.mjs` → `scripts/preview/map2dPerf.html`

Mounts the real `Map2D` and scripts a drag and a wheel gesture against it.

```bash
node scripts/renderMap2DPerf.mjs --points=30000 --style=parchment --dpr=2 --label=before --out=tmp
```

**`--dpr` matters and its failure is silent.** At `deviceScaleFactor` 1 the
settled quality DPR equals `INTERACTION_DPR`, so `setQualityDpr` bails, the
offscreen is never redrawn on mousedown, and the gesture redraws this script
exists to measure **do not happen at all**. It will happily report numbers for
work it never did. 2 is what the machine actually runs at.

Emits `settled.hash` — a hash of the settled frame's pixels — beside the
timings, so one run carries both the measurement and the proof nothing moved.

### `scripts/perf/renderProbe.ts`

Render counts with **prop attribution**. Add `useRenderProbe('Foo', { propA })`
to the components under study, measure, then **revert the call sites**. It lives
under `scripts/` and not `hooks/` precisely so that committing those call sites
feels wrong; nothing imports it, so it never bundles.

"`Controls` rendered 40 times during a slider drag" is not actionable.
"…and `factionColors` changed identity on 40 of them" names the fix.

### Determinism digest, for generation work

`tests/helpers/worldDigest.ts` hashes exact IEEE-754 bits. Use it as the gate on
any generation change claiming to preserve output.

**Known gap:** `CELL_FIELDS` omits `crustType`, `crustThickness` and
`upliftAccum`. Any change to the coarse→fine tectonic projection is therefore
*invisible* to the standard digest — a match proves nothing about that stage.
Adding them is not a one-liner; see HANDOFF S28 for why (the worker payload's
`presence` bitmask is a full `Uint8Array`). Until then, hash those fields
yourself when touching that stage, as `e890463` did.

---

## The order of operations

1. **Measure before forming a hypothesis.** The session-28 globe bench found the
   graticule costing 1.282 ms/frame against 0.05–0.33 for every other tenant —
   and 0.818 ms of that was the *walk*, not the projection. Neither fact was in
   anybody's guess list.
2. **Attack the measured hot spot, not the ugly code.** They are rarely the same
   place.
3. **Take the correctness number in the same run.** Not afterwards, when the
   change already feels finished.
4. **Record what you tried and reverted, with the number that killed it.** A
   measured dead end is a real result and stops the next person repeating it.
   The repo's HANDOFF discipline covers this; perf work generates more refuted
   hypotheses than most work, so it matters more here.

---

## Delegating perf work

The session-28 failure is worth stating precisely, because the obvious lesson is
the wrong one.

**What happened.** Four `opus-high` agents were dispatched simultaneously, each
told to investigate and optimise one subsystem, each starting from a cold
context. All four hit a hard spend ceiling **during orientation** — while reading
the repo, before producing anything. Zero commits. One had written a genuinely
good bucket-grid nearest-neighbour and never ran a single verification on it; it
was later verified by hand and landed as `e890463` at 54–118× on its stage.

**What that does not show.** It is not evidence about model throughput or about
parallelism in general. Nobody got past reading.

**What it does show.** Four cold contexts each paying full repo-orientation cost
*before* producing anything is the wrong shape whenever a hard ceiling exists
anywhere in the run. Orientation is pure overhead multiplied by the fan-out.

**The shape that works instead:**

- **Measure first, centrally.** One session runs the harnesses and produces a
  ranked list of hot spots with numbers.
- **Then delegate against a written plan**, one hot spot per agent, with the
  measurement already in the brief. An agent handed "the graticule tenant costs
  1.282 ms/frame, 0.818 ms of which is the walk; here is the harness command;
  reduce it without changing `max screen error` or `visibility flips`" needs no
  orientation phase at all.
- **Never delegate "go investigate."** That is the orientation-cost trap by
  another name.
- **Verify centrally.** The dispatching session runs the gate, in a clean tree.
  An agent's self-reported "green" is not evidence.
- **Cap local compute in the brief.** This machine is a fanless M1 Air that is
  often also running After Effects. `VITEST_WORKERS=1`, targeted test files
  only, never the full suite, reuse the dev server on port 3000 and never start
  a second Vite. See the run policy in `CLAUDE.md`.

**One trap specific to agent worktrees:** they are created under
`.claude/worktrees/`, *inside* the repo, so any gate that globs the tree sees N
copies of the codebase. This blew the lint ratchet to 146 warnings against a cap
of 30 with four agents running. `.claude/**` is now ignored in
`eslint.config.js` and `.gitignore`, but check any new gate you add.

---

## Open leads, with numbers

Measured, not guessed. None of these are started.

| Lead | Evidence |
|---|---|
| **Graticule tenant dominates overlay cost** | 1.282 ms/frame vs 0.05–0.33 for every other tenant, at 6k cells. 0.818 ms is the walk with a no-op projector, so the win is in the traversal, not the projection. |
| **Colour math dominates geometry refill** | 1.662 ms of a 2.234 ms full refill at 6k. Splitting position and colour updates only pays when colours are stable — the colour computation itself is the cost, so a split alone will disappoint. |
| **Fused MVP projector** | 0.642 → 0.460 ms/frame for 6k cells. Candidate is implemented in the bench; it needs its `max screen error` / `visibility flips` numbers before anyone trusts it. |
| **`computeShadeMap` is refill-scale** | 2.160 ms at 6k, comparable to a full geometry refill. Check whether it recomputes when it need not. |

Re-measure at 30k+ before acting: every number above is 6k, and these costs do
not all scale the same way.
