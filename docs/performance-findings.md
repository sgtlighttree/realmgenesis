# Performance findings (F4) — measured baseline

Companion to [performance-workflow.md](performance-workflow.md), which holds the
*method*. This file holds the *numbers*: a measured baseline and a ranked
implementation queue, so the next session can start optimising instead of
starting by finding out what is slow.

**Provenance.** Four agents were dispatched at F4 in session 28 and all four
died on a spend ceiling before reporting. Their instruments survived and are now
in the tree; everything below was measured by running those instruments, not
taken from the agents' claims. Where an agent left a *lead* rather than a
measurement, it is labelled as such.

**Confidence.** One seed, one machine (8-core M1 Air), **with After Effects
running**. Absolute milliseconds are therefore pessimistic and machine-specific.
**Ratios between rows of the same table are the trustworthy part** — they were
measured back-to-back in one process. Re-measure before and after any change
rather than diffing against these numbers months later.

---

## The headline

Two findings dominate everything else:

1. **A Map2D pan or zoom at 30k cells costs ~2.2 seconds of blocking main-thread
   redraw, at each end of the gesture.** This is the worst number in the
   project.
2. **The globe's screen projection can be made ~11x faster with output that is
   exact to 4.55e-13 px and zero visibility flips** — and the work is already
   written and validated inside the bench. It is an implementation task, not a
   research task.

---

## 1. 3D globe — `node scripts/renderGlobeBench.mjs --points=30000 --frames=120`

### Baseline, 30k cells / 179,988 tris

| Item | ms/frame |
|---|---|
| project all 30000 cells | **3.680** |
| graticule tenant | **2.375** |
| currents tenant | 1.431 |
| contours tenant | 1.406 |
| rivers tenant | 0.768 |
| borders tenant | 0.401 |
| labels tenant | 0.063 |
| dymaxion tenant | 0.057 |
| routes tenant | 0.003 |
| matrixKey gate | 0.002 |

Overlays alone are **~10.2 ms/frame** (projection + tenants) against a 16.7 ms
budget at 60 fps. That is the whole problem in one line.

### The validated ladder — do this first

> **Rung 1 SHIPPED (S29, `9bc20fe`):** the `hypot`→`sqrt` swap is now in
> `utils/screenProject.ts`. The per-tenant speedups below are realised.
>
> **Rung 2 SHIPPED (S29, `2342249`):** the staged typed-array + fused-MVP
> projector is ported into the app. It lives as pure functions
> (`stageCellPoints`/`projectStaged`/`projectLocalPoint`) in
> `utils/screenProject.ts`, wired into `ScreenOverlay`; guarded by
> `tests/screenProjectStaged.test.ts` (parity vs the naive path, 0 flips). Bench
> 1.373 → 0.242 ms/frame at 20k. **The globe projection ladder is complete** —
> the remaining globe levers below are the graticule drape walk and the
> geometry-refill colour cache, not projection.

The globe agent's last words before it died were *"typed-array staging and
hypot→sqrt"*. Both are **already implemented and correctness-checked** in
`scripts/perf/globeBench.ts`. Three steps, each measured:

| Step | project 30k cells | vs baseline |
|---|---|---|
| baseline | 3.680 ms | — |
| + fused MVP (hoist matrix elements out of the loop) | 3.134 ms | 1.17x |
| + staged typed arrays + `sqrt` | **0.320 ms** | **11.5x** |

`stage()` costs 0.738 ms and runs **per world, not per frame**.

**Correctness, already measured:** `max screen error 4.55e-13 px over 8797 pts,
0 visibility flips` — for both the fused and the sqrt variants. That is
floating-point noise, not a behaviour change.

**The root cause is one call.** `utils/screenProject.ts:25,29` uses
`Math.hypot`, which V8 implements with overflow-safe scaling and is far slower
than the naive form. Replacing it with `Math.sqrt(dx*dx + dy*dy + dz*dz)`
cascades into every tenant, because they all project through `isVisible`:

| Tenant | baseline | with `sqrt` | speedup |
|---|---|---|---|
| contours | 1.406 | **0.403** | 3.5x |
| borders | 0.401 | **0.100** | 4.0x |
| rivers | 0.768 | **0.158** | 4.9x |
| graticule (smooth) | 0.639 | **0.180** | 3.6x |
| currents | 1.431 | 1.017 | 1.4x |

**Start with the `hypot`→`sqrt` swap.** It is a two-line change in one file,
it is the largest single lever, and its correctness number already exists.

### The graticule is a separate problem

It is the most expensive tenant, and the cost is **not** projection:

| Variant | ms/frame |
|---|---|
| graticule, full | 2.375 |
| …with a no-op projector (pure walk cost) | **1.872** (79%) |
| …smooth mode, no drape walk | 0.639 |
| …smooth + sqrt | 0.180 |

79% of the cost is the drape walk — subdividing each graticule line to follow
terrain relief. Smooth mode skips it and is 3.7x cheaper. The lever here is
**adaptive subdivision** (fewer samples where relief is flat), not a faster
projector. Untried.

### Geometry refill — the obvious idea is the wrong one

| Item | ms/frame |
|---|---|
| full refill (179,988 tris) | 11.692 |
| **colour math only** | **10.140** |
| seasonalDelta only | 0.498 |
| `computeShadeMap` | **12.768** |

**Do not split position and colour updates.** That is the intuitive fix and the
numbers kill it: colour math is 87% of the refill, so a split that still
recomputes colour saves ~13%. The real lever is **not recomputing colour at
all** when `viewMode`, season and the colour maps are unchanged — i.e. cache the
colour buffer and refill positions only.

**`computeShadeMap` costs more than an entire geometry refill.** Check whether
it is memoised against `(cells, seaLevel)` before anything else; if it runs per
frame or per refill, that is a larger win than everything above it in this
table.

---

## 2. 2D map — `node scripts/renderMap2DPerf.mjs --points=30000 --dpr=2`

**This is the worst result in the project.** At 30k cells, parchment, DPR 2
(canvas 2396x1198):

| Phase | p50 | max |
|---|---|---|
| pan — mousedown redraw | **2192.8 ms** | 2192.8 ms |
| pan — moves (blit only) | 0.7 ms | 6.0 ms |
| pan — settle redraw | **2186.4 ms** | 2186.4 ms |
| zoom — first wheel | 1.0 ms | **2194.4 ms** |
| zoom — settle | 5.3 ms | 2173.7 ms |

Grabbing the map costs **~2.2 s of blocked main thread**, and letting go costs
another ~2.2 s. Dragging itself is free (0.7 ms — it blits), so the cost is
entirely the two DPR-switch redraws bracketing the gesture.

It also appears to scale badly: the same harness at 8k measured ~216 ms for the
equivalent redraw. 3.75x the cells for ~10x the time — **apparent superlinearity,
n=1 at each size, worth confirming before relying on it.**

Nobody has looked at where those 2.2 s go. That is the first job: the harness
reports the gesture, not the breakdown, so instrument inside the redraw
(glyph placement? `computeShadeMap` again? per-cell d3-geo path generation?)
before choosing a fix.

**Trap, already paid for:** run this at `--dpr=2`. At DPR 1 the settled quality
DPR equals `INTERACTION_DPR`, `setQualityDpr` bails, the offscreen is never
redrawn on mousedown, and **the redraws this harness exists to measure never
happen** — it reports fast numbers for work it did not do.

Every run also emits `settled.hash`, a hash of the settled frame's pixels. Use
it as the "nothing changed" evidence: baseline at 30k/parchment is `e3edcbbd`.

---

## 3. Generation — largely done

**Landed (`e890463`): the coarse→fine nearest-macro lookup.** Was a brute-force
argmin over `simulationResolution` (default 10000) macro points per display
cell — 2e8 chord distances at 20k cells, 2e9 at 200k. Now a bucket grid:

| display cells | brute | grid | speedup |
|---|---|---|---|
| 3,000 | 83 ms | 6 ms | 13.6x |
| 20,000 | 387 ms | 7 ms | 53.6x |
| 80,000 | 2520 ms | 21 ms | 117.9x |

Verified index-for-index against the original loop (0 mismatches, including a
zero-jitter lattice — the maximum-tie case) and by byte-identical world digest.
Guarded by `tests/macroPointGrid.test.ts`.

**Unverified lead, from the agent's own comment:** it claimed that lookup was
"60% of a 80k-point world". That figure was never substantiated and should be
treated as a guess. Whole-world timings at 20k moved 5915→2995 ms and
5294→3886 ms, but those are n=1 and noisy.

**Watch out when measuring generation:** `CELL_FIELDS` in
`tests/helpers/worldDigest.ts` omits `crustType`, `crustThickness` and
`upliftAccum`, so the standard digest **cannot see** changes to the coarse→fine
projection stage. Hash those fields by hand there. Widening it means widening
the worker payload's `presence` bitmask past a full `Uint8Array` — see HANDOFF
S28.

---

## 4. React re-renders — not measured

The React agent built `scripts/perf/renderProbe.ts` and died before taking a
single measurement, so **there is no baseline here at all**. Nothing in this
section is a finding; it is a starting point.

Its instrumentation technique is worth copying, though, because it is
non-invasive: rather than editing each component, it wrapped them at the call
site in `ShellApp`.

```tsx
const probed = (name, C) => (props) => { useRenderProbe(name, props); return <C {...props} />; };
const ProbedWorldViewer = probed('WorldViewer', WorldViewer);
```

Every inserted line carried a `/*F4PROBE*/` marker so the whole instrumentation
is greppable and removable in one pass. Do that.

It chose to probe `ShellApp`, `WorldViewer`, `Map2D`, `Controls`,
`ViewControls`, `Inspector`, `MiniMap`, `Select`, `CheckboxMenu` and
`EditToolbar` — a reasonable target list inherited free.

`window.__rc.dump()` returns counts **and which prop changed identity**, which
is the part that makes a result actionable. Measure four interactions: slider
drag, overlay toggle, view-mode change, map pan.

**Architectural constraint:** all state lives in `useWorldEngine` and is
prop-drilled by design. Memoisation, callback stability and component splitting
are in scope; introducing Context or a state library is **not** — that is an
architecture change, not an optimisation.

---

## Ranked queue

Ordered by (measured win) / (risk x effort). The first three all have their
correctness numbers already.

| # | Change | Expected | Confidence |
|---|---|---|---|
| 1 | `Math.hypot` → `Math.sqrt` in `utils/screenProject.ts` | 3.5–4.9x on most tenants | **Measured, 0 visibility flips** |
| 2 | Check `computeShadeMap` memoisation | up to 12.768 ms/frame | Measured cost; memoisation status unchecked |
| 3 | Fused MVP + staged typed arrays in the projector | 3.680 → 0.320 ms/frame | **Measured, written, validated** |
| 4 | Find where Map2D's 2.2 s gesture redraw goes | unknown, but it is the biggest number here | Cost measured, cause unknown |
| 5 | Cache the colour buffer; refill positions only | up to 10.140 ms/frame | Measured; needs an invalidation rule |
| 6 | Adaptive graticule subdivision | up to 1.872 ms/frame | Measured cost, approach untried |
| 7 | Baseline the React re-render counts | unknown | No data at all |

**Nothing here is byte-identical-optional.** Every item is a pure optimisation:
if output changes, the change is wrong. Take the correctness number in the same
run as the timing, not afterwards.
