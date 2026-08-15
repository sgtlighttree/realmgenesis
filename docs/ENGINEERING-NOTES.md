# Engineering Notes

Decisions with their rationale, hypotheses that were tried and refuted, and
levers that were designed but deliberately shelved. This is the "why" that
evaporates fastest — recorded so a later session doesn't re-derive a dead end or
re-litigate a settled call.

HANDOFF.md is where these start life (at the confidence level of the moment);
they graduate here once they've stopped changing. Session numbers refer to
`HANDOFF.md` entries.

---

## Refuted hypotheses

Keep these visible. A wrong idea that looked right is useful — it stops the next
session re-deriving it.

### Uplift accumulation over 20–40 timesteps (REFUTED)
The V3 design's headline seam-fix was to accumulate tectonic uplift over many
small timesteps. **Wrong.** With small per-step rotation, the same cell-graph
edge is re-selected as the boundary every step, so uplift piles onto one edge
and produces a *taller, thinner* wall exactly on the Voronoi cut — the seam it
was meant to remove. Recorded in the V3 design spec §5.1.

### Microplates as post-hoc peeling (REFUTED)
Peeling cells off existing plates *after* assignment desyncs the
`plates[plateIds[i]]` lookup and can reintroduce exclaves. The shipped approach
is **seed injection**: append a new `PlateState`, plant one seed cell, and let
the existing per-timestep Dijkstra region-growth grow it into a connected
region — reusing connectivity-by-construction. Verified 0 exclaves across seeds,
plate count 8→10 (Session 10).

### Divergent rift-lowering applied to oceanic crust (REFUTED)
The old code lowered **all** divergent-boundary cells via `upliftAccum`. For
oceanic ridges that fought the GDH1 age→depth model, which is now the sole
source of oceanic-divergence elevation. Rift-lowering is restricted to
continental–continental divergence (`crustA===1 && crustB===1`) in
`utils/tectonicsV3.ts`; oceanic divergence elevation is GDH1's job (Session 10).

### Band/chain seeding (`plateElongation`) as a macro-shape lever (REFUTED as significant)
`plateElongation` grows each plate's Dijkstra source into a velocity-aligned
chain, but the chains cap at ~3 cells (0.4) / ~5 cells (1.0) against plates that
Dijkstra grows to hundreds of cells — visually **near-inert** at the
macro-silhouette level (renders at 0.4 and 1.0 look identical to baseline). The
real de-blob levers turned out to be `plateJitter` (size/position variety) and
`boundaryRoughness` (jagged boundaries), both range-extended to 0–3 and
defaulted to 1.5 (Session 12). `plateElongation` is kept at default 0.4 because
it's cheap and the shelved anisotropy lever below reuses the param.

---

## Shelved, not abandoned — the next D7 levers

The heuristic-tier plate model is considered **done** (ROADMAP D7 ✅). These are
recorded for a future revisit, in rough order of value. None is scheduled.

### 1. Anisotropic Dijkstra growth cost — the real shape lever
Elongate the *whole* plate along its motion (not just the seed) by making the
Dijkstra relax step cheaper along the plate's velocity direction. In the inner
loop, per popped cell:

```
cost × (1 + k·(1 − |dot(edgeDir, v̂)|))
  v̂ = normalize(cross3(plates[plate].eulerPole.axis, macroPoints[cell]))
  k  = plateElongation · ~0.8   // clamp so max stretch ≈ 2–3:1
```

Deterministic, connectivity-preserving (still region-growth), computed once per
popped cell. Clamp `k` — over-strong reads as fake cigars. This is the "real"
version of what band-seeding failed to do. Design detail is in the Task-1b brief
under the D7-polish plan. Deferred because `plateJitter` + `boundaryRoughness`
already delivered the de-blob.

### 2. Transform-edge fracture
Feed step *k*'s transform-boundary classification into step *k+1*'s cost field
(recompute as a **set** each step, never accumulate — see the refuted
accumulation hypothesis). Shelved: `boundaryRoughness` at 1.5 already provides
the jaggedness; this is marginal on top.

### 3. Cortial boundary-curve rebuild — deliberate NO-GO
The "properly grounded" model: plates as a topological graph of simulated
boundary curves + terranes, not seed-grown regions. **Both the `advisor` tool
and the `fable-advisor` subagent independently ruled NO-GO** (Session 12), and a
render confirmed it:

- It replaces the entire plate-assignment substrate — killing
  `assignPlatesDijkstra`, `computeEdgeCosts`, `mergeSmallPlates`,
  `injectMicroplates`, and the **0-exclave-by-construction invariant**
  (Session 9). Curve intersection + retriangulation on a sphere is weeks of
  geometric-robustness edge cases.
- The commissioned research report's **own** engineering recommendation for a
  ~10k-cell browser budget is the heuristic tier — which is already shipped.

The path stays open but is not recommended for a browser budget. Method is in
the research report cited in HANDOFF Session 10/12.

---

## Deferred — Milestone D climate depth

### Seasonal wind + moisture (monsoons) — deferred from D1
D1 (seasonal cycle) shipped with **temperature + ice/biome edges** moving through
the year, but **wind and moisture stay annual**. The richer variant: rerun the
wind-band classification and the 8-pass moisture transport solver per season so
wet/dry seasons and monsoon reversals appear (subsolar point drags the ITCZ, wind
bands migrate, rain shifts hemisphere).

Deferred deliberately, two reasons:
1. **Cost.** The 8-pass moisture solver is the expensive part of the climate pass.
   D1's temperature excursion is a closed-form per-cell formula
   (`utils/seasons.ts`), recomputed live in the render layer for free. Seasonal
   moisture cannot be — it needs the solver, so it would be either a per-season
   regeneration or a set of precomputed seasonal moisture snapshots.
2. **Architecture.** D1's seasonal biome recompute is a free O(n) pure function
   because moisture is annual: `determineBiome(height, T(s), moistureAnnual, sl)`.
   If moisture also varied by season, biome-at-season would depend on a per-season
   moisture field, breaking the "recompute biome live in the color path" seam.

If revisited: precompute N seasonal moisture fields at generation (N≈4–12),
store per-cell, and interpolate — or gate it behind an explicit "simulate seasonal
moisture" toggle so the default stays cheap. Pairs naturally with D2 (ocean
currents), which also feeds the moisture model.

**D2 composability seam (read this before building seasonal moisture).** D2 (ocean
currents, shipped S15) injects a warm-current **evaporation** term into the ocean
moisture *seed* (`1.0 + EVAP_K·max(0, sstAnomaly)` in the `worldGen` 8-pass, from
`utils/currents.ts`). This is the **annual** ocean-moisture baseline and introduces
**no per-season field**, so D1/D3's free O(n) biome-at-season recompute is preserved.
When seasonal moisture lands, it must **layer its per-season overlay on top of** this
current-modified annual baseline — do **not** overwrite the ocean seed back to a bare
`1.0`, or warm-current coasts silently lose their extra rainfall. The current field
itself stays annual/steady-state (no seasonal gyre reversal) at this tier.

### Sea-level coupling for ice — deferred from D3
D3 (sea-ice) ships as a render overlay: cells below the seawater freeze point
render as ice, no change to height or coastline. The "grounded" extension —
glacial/sea ice mass lowering global sea level (more ice → lower `seaLevel` →
exposed shelf) — is deliberately deferred. It changes the **coastline**, so it is
a generation-stage / regeneration concern, not a render overlay, and it would ripple
into hydrology, biomes, and civs. If revisited, it belongs alongside a seasonal-
or climate-driven `seaLevel` adjustment computed at generation, not in `getCellColor`.
Glacial land cover beyond the existing `ICE_CAP` biome + snow is also deferred —
revisit only if polar lowlands read wrong once sea-ice is in.

### D5 planetary params — day length, gravity, moons deferred
D5 shipped **host star class only** (`starClass`, a live temperature hook). The
other three listed params were deliberately NOT built, because none has a
principled mechanical hook in the current model, and each would ship a permanent
save-schema field:
- **Day length** — no diurnal cycle exists; its only real hook is Coriolis
  strength, which is D2 (ocean currents / wind) territory. Revisit with D2.
- **Gravity** — there is **no principled mapping** from *g* to the normalized
  `[0,1]` height field. Wiring it into the Stage-9b remap would be a fudge factor
  that duplicates `mountainHeight` — a knob that lies. Do **not** add it as a
  relief scaler. A real hook would need a physical erosion/isostasy model.
- **Moons** — no tide model, so it would be a pure lore string. Revisit only if a
  tide simulation is ever built.
The guiding rule (same as the D3 sea-ice-temp decision): a param must do something
real and non-duplicative, or it doesn't ship — a hollow param is worse than none
because it's permanent once in the save schema.

---

## Standing decisions (settled, with rationale)

### The V3 terrain model is the only path
V3 (independent crust fields + Euler-pole kinematics + Dijkstra plate growth)
replaced the V2 "crust-is-plates" height model. The `V3_ENABLED` flag and all
V2 dead code (`assignPlates` argmin, `enforceConnectivity`, `randomVector`) were
removed at Session 11. `plateInfluence` was renamed `tectonicStrength` (Session
8); old saved values silently die — an accepted consequence.

### Plate connectivity is guaranteed at the macro level, not the fine mesh
Multi-source Dijkstra region-growth over the macro neighbor graph makes each
macro plate exactly one connected component by construction (a shared `claimed`
set keeps per-plate seed sets disjoint). **But** the macro→fine
nearest-macro-cell downsample can pinch a thin macro plate into disconnected
*display* strays. This is pre-existing and fires even at `plateElongation` 0.
`tests/plateConnectivity.test.ts` guards it for seed `realmgenesis` specifically
— it is **not** a general invariant (its own comment says so). Fixing the
downsample (BFS cleanup or a better projection) is a real open task.

### Seafloor bathymetry from crust age (GDH1)
Oceanic floor with a valid crust age follows GDH1 subsidence (Stein & Stein
1992) instead of flat isostasy: deeper away from spreading ridges. Depth is
mapped into the existing oceanic height *band* (ridge ≈ −0.5, old floor ≈ −0.85),
**not** meters into the raw field — otherwise the global min shifts and
normalization rescales the land fraction. Fed in before normalization so
`seaLevel`, the Stage-9b remap, climate, and erosion are untouched (Session 10).

### seafloorDepth is a linear datum; oceanDepth is a contrast curve
The two ocean-floor knobs are complementary, not redundant. `oceanDepth` is a
power curve (reshapes trench-vs-shelf contrast); `seafloorDepth` is a linear
multiplier on depth-below-`seaLevel` (shifts mean water depth, preserves
relative shape). Both hold the coastline fixed and live in the same Stage-9b
block. `seafloorDepth` replaced the former `seafloorDetail` texture knob, whose
two internal effects (abyssal-hill amplitude, GDH1 noise-damping) were baked at
their 0.5 default (2026-08-14). Full rationale:
`docs/superpowers/specs/2026-08-14-seafloor-depth-datum-design.md`.

### Projection noise is damped over deep ocean
`projectTectonicsToDisplay` blends structural noise at weight `1.2 −
tectonicStrength`; over deep-ocean cells it's further scaled down (baked
0.675 factor) so the GDH1 age→depth gradient isn't washed out, plus abyssal-hill
noise at a fixed low amplitude (Session 10 fable-advisor correction).
