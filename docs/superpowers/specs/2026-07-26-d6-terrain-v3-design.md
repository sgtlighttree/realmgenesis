# D6 — Terrain Generation V3: design

**Status:** design approved by Matt, 2026-07-26. **No implementation plan yet, deliberately.**
Matt asked for brainstorm output only, for a future session to pick up. The
`writing-plans` step was consciously skipped, not forgotten.

**Companion research** (read before implementing):
- [`docs/research/2026-07-25-realistic-terrain-generation.md`](../../research/2026-07-25-realistic-terrain-generation.md) — prior-art survey (Antigravity/Gemini)
- [`docs/research/2026-07-25-tectonics-adversarial-pass.md`](../../research/2026-07-25-tectonics-adversarial-pass.md) — adversarial red-team (Claude sonnet-medium). **This one invalidated part of the first design.**

---

## 1. Why V3 exists

Three defects, each verified by reading `utils/worldGen.ts`, not inferred.

### 1.1 Continents ARE plates — the defect Matt actually reported

> "the biggest issue is how painfully obvious the plates are. A lot of the time
> there is a hard divide between land and water with the plate acting as the
> boundary. In some cases, continents just look like big islands defined by the
> plates underneath."

The cause is exact. Each plate is flipped wholesale to land or ocean:

```js
const isLand = pRng.next() < landChance;
plateHeights[i] = isLand ? (landLevel + pRng.next() * 0.3) : (oceanLevel + pRng.next() * 0.3);
```

and that value reaches the cell smoothed over **one ring of neighbours only**:

```js
c.neighbors.forEach(n => { baseSum += plateHeights[cells[n].plateId]; bCount++; });
const avgBase = baseSum / bCount;
let height = avgBase * influence + structuralNoise * (1.2 - influence);
```

So elevation carries a piecewise-constant, per-plate offset, and the coastline —
a level set of that field — traces plate edges by construction. At
`plateInfluence = 1.0` continents and plates are identical.

The independent research pass named this, unprompted, as the classic failure:
"if a plate is entirely continental, continents end up looking exactly like the
Voronoi polygons of the plates themselves."

### 1.2 The faceted seam along plate boundaries

Stress is computed only on boundary cells, then box-blurred over the neighbour
graph, and uplift is applied through a linear ramp behind a hard threshold:

```js
const edgeProx = Math.max(0, 1.0 - distToEdge[c.id] * 0.5);
if (edgeProx > 0) {
    if (stress > 0.05) { const mtnHeight = stress * edgeProx * 1.5; ... }
```

Three compounding causes: a step function on a discrete lattice *is* a faceted
edge; `distToEdge` is a `+0.1`-per-iteration diffusion approximation rather than
a real geodesic distance; and the falloff spans about one ring.

### 1.3 `detailLevel` cannot produce sub-cell detail

Height is sampled once per cell, at `c.center`. Raising the FBM octave count
aliases higher frequencies into a single value per cell. Sub-cell detail
necessarily requires evaluating height at more points than there are cells.

### 1.4 Plate motion is not a rigid motion

```js
const plateDrift = plateVectors.map(() => ({ x: rng()-0.5, y: rng()-0.5, z: rng()-0.5 }));
```

A random linear vector is not a rigid motion on a sphere — it has no axis, so a
plate's "drift" is inconsistent across its own extent. Since relative velocity
drives *all* boundary classification, the whole tectonic signal is built on a
malformed primitive. Real plates rotate about an **Euler pole**, and velocity
varies with angular distance from it. (Found by the research pass, not by us.)

---

## 2. Decisions, with rationale

| Decision | Rationale |
|---|---|
| **Full rebuild, not patches** | Matt: "if we can tear down the flawed foundation for a new one, why not?" The defects are in the model, not the code quality. |
| **Crust and plates are independent fields** | The direct fix for §1.1. A plate must be able to carry continent *and* ocean. |
| **Bounded kinematic simulation** | Gives asymmetric, history-bearing boundaries. Full geological history (Wilson cycle) was offered and **declined** as a research project in its own right. |
| **Euler poles** | §1.4. Also produces natural variation in boundary intensity along a single margin. |
| **Generation moves to a Web Worker** | See §6.0 — the justification changed during design; recorded honestly. |
| **Sub-cell detail is RENDERING ONLY** | Cells stay the single authority on height, so biomes, rivers, civs, the Inspector and the paint/undo model all keep agreeing. Avoids a town rendering offshore. |
| **Simulate coarse, project once** | No implementation found runs multi-step CPU simulation at 200k cells. |

---

## 3. Architecture

Two **independent** per-cell fields. This is the core of V3.

| Field | Source | Governs |
|---|---|---|
| `crustType` (continental/oceanic), `crustThickness` | low-frequency noise on its own RNG stream | the shape of continents |
| `plateId`, Euler pole, angular rate | plate seeds | tectonics **only** |

**Plates never contribute base elevation.** They deform crust; they do not decide
where land is.

### 3.1 Coupling — Matt's mountains concern

Fully independent fields mean mountain belts track a different random field than
coastlines, so ranges cut across continents arbitrarily. Matt asked whether
seeding can be *coupled* instead. It can, cheaply, by reusing machinery the seam
fix needs anyway.

The boundary-membership test requires a domain warp regardless (§5). Let one
weighted term of that warp be driven by the **gradient of the crust field**:

```
warp(p) = noiseWarp(p) + marginCoupling * ∇crust(p)
```

Boundaries then bend toward continental margins where the coupling term is
strong, producing Andes-style active margins, without margins *determining*
boundaries.

**Read the GRADIENT, not the field value.** Coupling to the value would make
boundaries hug or avoid continents wholesale, which is a weaker restatement of
"continents are plates" — the defect V3 exists to remove.

**Bound the knob; do not expose 0–1 as a free slider.** At `marginCoupling → 1`
every boundary hugs a margin: every continent gets a range on all sides, and
there are no mid-ocean ridges and no intra-continental rifts. That is as wrong as
0, in the opposite direction. **Default 0.3–0.4.**

**Stated limitation (do not re-derive this):** this is *geometric correlation,
not causal history*. It makes belts tend to run along margins; it does not make
the margin exist *because* the crust broke there. That is the Wilson cycle, it
requires the full-history model Matt declined, and it is the honest reason
coupling is a knob rather than a mechanism.

### 3.2 Crust is never advected

**Do not move the crust field.** Resampling it at back-rotated positions each
timestep interpolates repeatedly, and after ~30 steps continents diffuse into
mush and coastlines go soft.

Instead: keep crust fixed and **inverse-rotate each cell's position into each
plate's rest frame**, testing membership by nearest rotated seed. Nothing is
interpolated, so nothing blurs.

This is what five independent implementations converge on — either exactly-rotated
Lagrangian markers with nearest-neighbour reassignment (Blender/Cortial reimpl.,
Clustered Convection, GPlates tracer methods), or not moving crust at all (Red
Blob Games). Fancier interpolation (BFECC, MacCormack) does *not* solve it; those
still interpolate every step.

**Corollary:** any scalar accumulated on top (uplift, thickness) must stay
attached to the crust identity being tracked, never resampled independently, or
diffusion re-enters through the back door.

---

## 4. The simulation model

```
seed crust field           (own RNG stream, independent of plates)
seed plates                (Euler pole + angular rate each)

for t in 0..N:                                   # N ≈ 20–40, at COARSE resolution
    for each cell: inverse-rotate into plate rest frames,
                   assign membership by nearest rotated seed
                   (warped per §3.1 — this is also the coupling channel)
    for each boundary pair:
        vn = relative velocity ⟂ boundary        # converge / diverge
        vt = relative velocity ∥ boundary        # shear
        continental + continental, vn<0  -> collision belt
        oceanic + continental,     vn<0  -> trench on the OCEANIC side,
                                            arc offset onto the CONTINENTAL side
        oceanic + oceanic,         vn<0  -> trench + island arc (older subducts)
        vn > 0                           -> rift, then new oceanic crust
        |vt| dominant                    -> transform: offset, little relief
    accumulate uplift

project coarse result -> full cell mesh          # ONE nearest-neighbour pass
compose height = f(crustThickness, upliftAccum, isostasy) + structural noise
normalize; then existing erosion / climate / biomes / rivers, unchanged
```

### 4.1 Resolution decoupling

Run the loop at **5k–20k macro-cells regardless of the display cell count**, then
project once onto the full mesh, adding structural noise and isostasy at full
resolution in that final pass.

Precedent: World Orogen simulates ~20k regions and projects to 2.56M display
cells; tectonics.js caps near 10k cells for interactive JS simulation. Nothing
found runs 20–40 CPU steps at 200k.

This generalises the project's existing principle: *simulation resolution, like
sub-cell detail, is decoupled from display resolution.*

### 4.2 Cases that must be handled explicitly

- **Triple junctions.** Pairwise classification is undefined where three plates
  meet. Without a tie-break rule, expect visible notches or spikes at every plate
  corner. These points are geometrically prominent.
- **Passive margins.** Because crust is seeded independently, continental/oceanic
  transitions will sit in plate *interiors*, far from any classified boundary.
  Earth's Atlantic coast is this case. If the compose step assumes all relief
  comes from boundary uplift, these render as unmotivated cliffs. They need a
  gentle isostatic shelf slope.
- **Euler pole placement.** A pole inside its own plate's territory gives
  near-zero velocity near the pole and fast motion far from it, so uplift
  intensity varies arbitrarily along one margin. Constrain placement.
- **Isostasy must saturate.** A linear thickness→elevation relation lets a
  long-lived collision belt grow without bound over 30 steps. Real crustal
  thickening saturates. Use a clamped or asymptotic response.
- **Do not double-count.** Uplift is the *kinematic forcing*; isostasy is the
  *static buoyancy response* to the resulting thickness. One collision must not
  contribute two independent upward terms.

---

## 5. Seam elimination

### 5.1 REFUTED HYPOTHESIS — read this before re-proposing it

**We believed:** accumulating uplift over 20–40 timesteps would replace the
one-ring ramp with a wide, irregular belt, and that this alone would remove the
faceted ridge. This was the headline seam fix in the first draft of this design.

**It is wrong, and can make the artifact worse.** If plate membership is a
nearest-seed test and per-step rotation is small relative to cell spacing, **the
same cell-graph edge is re-selected as "the boundary" on most steps**. Uplift
piles onto that one edge, producing a *taller, thinner* wall precisely aligned to
the Voronoi cut — a sharper version of the very bug it was meant to fix.

**Refuted by** the adversarial research pass, 2026-07-25. Kept here per the
HANDOFF discipline: a wrong idea that looked right is worth more than silence,
because the next session will otherwise re-derive it.

### 5.2 What actually fixes it

1. **Break cell-graph alignment explicitly.** Perturb the boundary membership
   test with per-cell noise / domain warp so the effective boundary wanders
   across a *band* of cells instead of snapping to one graph edge. (This is the
   same warp that carries the §3.1 coupling term.)
2. **Let the coarse→fine projection set the visual boundary width**, not the
   coarse simulation's cell graph.
3. **Replace the hard `stress > 0.05` threshold with a smooth falloff.** A step
   function on a discrete lattice is a faceted edge by definition.
4. **Use a real geodesic distance field** from the boundary, not the
   `+0.1`-per-iteration diffusion approximation.

Uplift accumulation is still worth keeping for *asymmetry and history* — it is
simply not the seam fix.

---

## 6. Staging

Three stages, each independently shippable and independently judgeable.

### 6.0 The worker's justification changed — do not restate the original one

The Web Worker was proposed, and agreed with Matt, on the grounds that it
**unlocks iterative simulation** that the main thread could not afford.

That reasoning is now weaker. Once the simulation runs at coarse resolution
(§4.1), 20–40 timesteps over 5k–20k macro-cells is cheap enough that it was
never really the blocker.

The worker still earns its place, but for a different reason: **the existing
pipeline is what stutters** — the 200k-cell Voronoi build, erosion, climate and
rivers, currently kept barely responsive by `setTimeout(0)` yields between
stages. Moving those off the main thread is the actual win, and it removes the
yield-based staging entirely.

Recorded because a future session reading only the conclusion would repeat the
original, now-inaccurate, justification.

### Stage 1 — Worker migration, algorithm unchanged
Move the existing pipeline into a Web Worker with **no algorithm change**.

**Gate: byte-identical output for a fixed seed.** The determinism suite already
tests exactly this, so the infrastructure change gets free proof of correctness
*before* any algorithm change makes output comparison meaningless.

This follows the DECOMPOSE mode in `CLAUDE.md` and mirrors the `useWorldEngine`
extraction (verbatim move, frozen return block, compiler carries the proof).

**In scope, and the biggest risk in the whole project:** `Cell` is an
array-of-objects with variable-length `vertices` and `neighbors`. Structured-
cloning 200k of those across the worker boundary will be slow enough to matter.
The Voronoi build should move into the worker too, returning flattened typed
arrays (SoA + offset buffers). **This is stage 1 scope — discovering it during
stage 2 would be the expensive way to learn it.**

### Stage 2 — The V3 terrain model
Crust/plate decoupling, Euler poles, bounded kinematic simulation at coarse
resolution, coarse→fine projection, the §5.2 seam fixes, and the §4.2 cases.

**Gate: Matt's eye.** There is no automated test for "does this look tectonic."

### Stage 3 — Continuous height field
`h(p)` evaluable at any point, sampled per-vertex or into a texture for sub-cell
rendering detail. 3D noise on Cartesian position — **never** 2D noise via
lat/lon UV, which pinches at the poles and seams.

Cells remain the authority on height. Stage 3 is a rendering concern.

---

## 7. Determinism and testing

- Every new random draw takes a **fresh RNG side-stream** (`seed + '_crust'`,
  `seed + '_plates_v3'`, …), per the existing project convention.
- The determinism suite must stay green: same seed → same world.
- `tests/paramLiveness.test.ts` must be extended for every new `WorldParams` key
  (`marginCoupling`, timestep count, simulation resolution). That suite fails if
  a param stops influencing output — extend it, do not exempt it.
- Stage 1's gate is byte-identical output. Stages 2 and 3 have no such gate by
  nature; they change values deliberately.

---

## 8. Accepted consequences

- **Every existing seed produces a different world.** Terrain regenerates
  deterministically from the seed, so saved worlds will look different after V3.
  Accepted: D6 breaks *values*, not *interfaces*; derived layers (civs, cultures,
  routes) regenerate by design.
- **`plateInfluence` changes meaning** — from "how much plates dictate elevation"
  to "how strongly tectonics deform crust." Old saved values will read
  differently. Consider renaming rather than silently repurposing.
- **Belts are not causally tied to coastlines.** Per §3.1, coupling is geometric.
  Mountain ranges will sometimes cut across continents in ways Earth's do not.
  This is the accepted price of fixing "continents are plates," and no number of
  timesteps removes it.

---

## 9. Open questions for the plan phase

1. Does V3 replace V2 outright, or ship behind a flag for side-by-side comparison
   during development? (Leaning: temporary flag during stage 2, removed at the
   end.)
2. Should erosion move to edge-length-weighted diffusion? The research recommends
   it; the current Priority-Flood + flux + thermal model already avoids the
   particle-based approach the research warns against, so this may be a
   non-issue. Verify before changing.
3. Is Lloyd's relaxation of the Fibonacci seed points worth it? Fibonacci is
   already near-uniform, so probably not — but it is the standard mitigation for
   lattice artifacts and worth one measurement.
4. Exact `N` timesteps and coarse resolution: pick empirically against Matt's eye.

---

## 10. Sources and confidence

**Well-sourced** (five independent implementations agree): crust should not be
advected by interpolated resampling; simulate coarse and project once; 3D noise
rather than UV mapping; worker plus typed-array transfer.

**Inference, not confirmed:** every claim specific to **Cortial et al. 2019**.
Both research passes failed to reach the primary text — Wiley is paywalled and
the HAL and LIRIS mirrors were unreachable across several attempts. The
cross-implementation pattern above does **not** depend on Cortial and stands on
its own. Do not treat "Cortial does X" in either research doc as established.

**Struck:** the first research pass cites *"Erleben, K. et al., Lattice-aligned
artifact mitigation via Lloyd's Relaxation and Voronoi smoothing."* Two
independent checks failed to find it; treat it as fabricated. The underlying
advice about Lloyd relaxation is sound, but cite Lloyd 1982 or Red Blob Games
instead.

Correct citation: Y. Cortial, A. Peytavie, E. Galin, E. Guérin, "Procedural
Tectonic Planets," *Computer Graphics Forum* 38(2), Eurographics 2019,
DOI [10.1111/cgf.13614](https://onlinelibrary.wiley.com/doi/abs/10.1111/cgf.13614).
