# D2 — Ocean Currents (semi-simulated gyres → temperature + moisture)

**Status:** design approved 2026-08-15 (Session 15). Roadmap item **D2**.
**Author:** Claude (Opus 4.8) + Matt, via superpowers:brainstorming.
**Advisor:** three consults (fidelity/coupling framing, determinism constraint,
this pre-spec sharpening pass).

---

## 1. Goal & scope

Give the climate model its first ocean-driven signal. Today temperature is
purely latitude + star-class + elevation lapse + noise, and ocean cells are inert
(moisture pinned to `1.0`, no thermal effect). D2 adds a **steady-state ocean
current field** — wind-driven, Coriolis-deflected, blocked by continents into
emergent gyres and western-boundary currents — that carries a **sea-surface
temperature (SST) anomaly**. That anomaly:

- **moderates coastal temperature** (a warm poleward current warms downwind
  coasts — the Gulf-Stream / NW-Europe effect; a cold equatorward current chills
  west coasts), and
- **modulates evaporation** (warm-current coasts get wetter, cold-current coasts
  drier) by seeding the existing 8-pass moisture transport.

**Fidelity tier (Matt's call): semi-simulated gyres**, realized as a **fixed-pass,
fixed-order relaxation** — never a divergence-free / Poisson solver. Determinism
(byte-identical reruns) is the project's load-bearing test invariant; a solver
with convergence thresholds or variable iteration counts breaks it *and* the
browser budget. The relaxation forms *emergent* gyres without enforcing strict
incompressibility. That is the one honest limitation: gyres are approximate, not
PDE-exact. It is cosmetic — the coastal heat/moisture anomalies that actually feed
climate are correct whether or not a gyre closes perfectly.

**Coupling (Matt's call): temperature AND moisture.** This pulls forward part of
the deferred D1 seasonal-moisture work; §7 specifies the clean seam so the future
monsoon solver builds on D2 rather than rewriting it.

**Currents are annual / steady-state** — no seasonal reversal (monsoon-style
current flips stay with the deferred seasonal-moisture work). Consistent with D1's
annual-climate scope.

### Non-goals (YAGNI)

- A currents **view-mode / arrow overlay** (ROADMAP F2 territory). Climate
  coupling lands first; visualization is a clean follow-up.
- Seasonal current reversal / monsoon gyres (deferred with D1 seasonal-moisture).
- A second tuning param for evaporation strength — one `currentStrength` knob
  scales both couplings.

---

## 2. Why relaxation, not Sverdrup (rejected alternative)

The Sverdrup streamfunction approach (wind-stress curl → integrate to a
streamfunction → velocity = its perpendicular gradient) is more principled and
gives automatic western intensification. **Rejected** because its westward
integration and ocean-basin identification have no natural expression on an
irregular geodesic cell graph — there is no axis-aligned "west" and basins are not
axis-aligned regions. That is exactly the geometric-robustness edge-case territory
the advisor flagged; the "principled" advantage largely evaporates once the
integral is approximated on a messy graph. The relaxation is deterministic by
construction, cheap, robust, and forms gyres + boundary currents from the real
drivers (wind bands + continent geometry + Coriolis).

A full divergence-free (Poisson-projection) solver was rejected outright: it breaks
determinism and the browser budget (advisor's hard NO-GO tier).

---

## 3. Pipeline placement

New pure module **`utils/currents.ts`** (worker-safe — imports `types` and math
helpers only; no `three` / `d3`, matching `seasons.ts` / `planetary.ts`).

It runs **inside the existing climate block** in `utils/worldGen.ts`
(`onLog("Calculating Climate (Wind & Rain)...")`), in this order:

1. `windVectors` computed (existing).
2. **NEW: `currentStrength === 0` → skip everything below and keep today's exact
   behavior** (early return of the stage, §6).
3. **NEW: `computeOceanCurrents(...)`** → per-ocean-cell velocity field.
4. **NEW: `computeSstAnomaly(...)`** → per-ocean-cell SST anomaly (heat advection).
5. Moisture 8-pass (existing) — **modified** ocean-cell seed (§5).
6. Temperature loop (existing) — **modified** to add the current moderation term
   (§4).

Coastline is final by this point (post Stage-9b height remap), so ocean/land
classification (`height < seaLevel`) is stable.

---

## 4. The solver — `computeOceanCurrents(cells, windVectors, params)`

Returns a velocity per cell as a tangent vector `{x,y,z}` (or three parallel
`Float32Array`s), zero on land.

**Initialization.** For each ocean cell (`height < seaLevel`):
`v = DRAG · projectTangent(windVectors[i], normal_i)`, where `normal_i` is the cell
centre unit vector and `DRAG` is a small constant (~0.6). Land cells: `v = 0`.

**Fixed `N_VEL ≈ 12` passes**, each a `cells.forEach` writing into a
`Float32Array` scratch triple, then a copy-back — identical control-flow shape to
the existing 8-pass moisture loop, so iteration and accumulation order are fixed
and the result is deterministic with **no RNG** (currents are wind-derived):

For each ocean cell `i`:

- **(a) Coriolis deflection.** Rotate `v_i` about `normal_i` by
  `θ = CORIOLIS_K · currentStrength · sin(lat_i)` — right-handed (clockwise seen
  from outside) in the northern hemisphere (`lat > 0`), left-handed in the south
  (the sign of `sin(lat)` handles this). Rodrigues rotation about `normal_i`.
  **At `lat ≈ 0`, `θ → 0`: equatorial flow is pure wind-driven zonal flow with no
  deflection. This is physically correct (equatorial currents are wind-driven and
  zonal) and must NOT be "fixed" — it is specified behavior, not a degenerate bug.**
- **(b) Advective smoothing.** Blend toward the mean velocity of ocean neighbours:
  `v_i ← (1−MIX)·v_i + MIX·mean(v_n : n ∈ neighbors_i, ocean)`. Skips land
  neighbours. If a cell has zero ocean neighbours (isolated water pixel), keep
  `v_i` unchanged.
- **(c) Boundary tangency (no-normal-flow).** For each **land** neighbour `n`,
  remove the component of `v_i` pointing toward `n`: with `d = normalize(center_n −
  center_i)` projected tangent, if `dot(v_i, d) > 0` then `v_i ← v_i − dot(v_i,d)·d`.
  This makes currents run *along* coasts and produces western-boundary
  intensification from wind + continent geometry.

After the passes, clamp speed to a `MAX_SPEED` ceiling (guards against runaway
accumulation). Re-project to tangent each pass to keep `v` on the sphere surface.

Constants (`DRAG`, `CORIOLIS_K`, `MIX`, `N_VEL`, `MAX_SPEED`) are module-level
named constants, tuned during the render-verify step (§8), documented inline.

---

## 5. Heat advection → SST anomaly — `computeSstAnomaly(cells, velocity, params)`

Returns a per-cell `Float32Array sstAnomaly` (0 on land).

**Seed.** Each ocean cell starts at its own latitude temperature
`base_i = annualMeanLatTemp(lat_i, params)` (the temperature it would have with no
current — reuse `seasons.ts annualMeanLatTemp`; star-class scaling is applied
uniformly and cancels in the anomaly, so seed with the un-scaled latitude curve).
A working field `T_i` is initialized to `base_i`.

**Fixed `N_HEAT ≈ 12` passes** — same forEach + scratch shape. For each ocean cell
`i`, **accumulate over ALL upstream neighbours** (mirroring the 8-pass moisture
rule exactly — NOT a single-argmax upstream pick, which would reintroduce a
tie-break problem à la `assignPlatesDijkstra`):

```
incoming = 0; count = 0
for n in neighbors_i where n is ocean:
    edgeDir = normalize(center_i − center_n)   // from n toward i
    if dot(edgeDir, velocity_n) > 0:           // current at n flows toward i → n is upstream
        incoming += T_n; count++
if count === 0:                                 // convergence zone / no upstream: relax toward own base
    T'_i = T_i·0.95 + base_i·0.05
else:
    T'_i = T_i·(1−HEAT_MIX) + (incoming/count)·HEAT_MIX
```

**Anomaly.** After the passes: `sstAnomaly_i = T_i − base_i`. Poleward-flowing
(warm) currents carry low-latitude heat to high latitudes → positive high-lat
anomaly; equatorward currents → negative anomaly.

---

## 6. Coupling into temperature (`worldGen.ts` temperature loop)

Add, before the elevation lapse and noise, a current term:

- **Ocean cells** (`height < seaLevel`): `temp += sstAnomaly_i · currentStrength`.
  (So D3 sea-ice responds — a warm current holds back polar ice; a cold current
  freezes a coast that latitude alone would leave open.)
- **Land cells**: a **coastal moderation** term = proximity-weighted blend of
  adjacent ocean anomalies. Concretely, for a land cell, average `sstAnomaly` over
  its ocean neighbours (1-ring); if it has none, term is 0 (inland cells are
  unaffected). `temp += COAST_K · currentStrength · meanOceanNeighborAnomaly`.
  A 1-ring reach keeps it O(n) and physically honest (direct coastal contact);
  deeper inland penetration is a YAGNI knob.

**D1 preservation.** The current anomaly is baked into the stored annual-mean
`cell.temperature`. The D1 seasonal excursion (`seasons.ts
seasonalTemperatureDelta`, computed from *geometric* latitude, added at render) is
orthogonal and unchanged: at season 0.5 the excursion is 0 regardless of the
baseline, so the **neutral == canonical** and **return-to-exact** invariants still
hold — but against a **new baseline world** (coastal temps now differ). See §8 for
the re-verification protocol.

---

## 7. Coupling into moisture (`worldGen.ts` 8-pass seed) + monsoon seam

Today ocean cells are pinned to `moisture = 1.0` both at init and every pass. D2
replaces the pin with an **evaporation-boosted seed**:

```
ocean cell: moisture = 1.0 + EVAP_K · max(0, sstAnomaly_i)         // warm water evaporates more
```

(Only positive anomalies boost; cold currents don't suppress evaporation below the
`1.0` baseline — they cool the air, handled by the temperature term, and their
drying effect emerges downwind from carrying less-elevated moisture.) The existing
8-pass upwind transport then carries this elevated moisture onto downwind coasts →
warm-current coasts wetter. `EVAP_K` tuned in §8.

**Composability seam with deferred D1 seasonal-moisture / monsoons (advisor's
flag).** This term modifies the **annual** ocean-moisture baseline. It introduces
**no per-season moisture field**, so D1/D3's free O(n) biome-at-season recompute
(`getCellColor` display-biome fast path uses a single stored annual moisture) is
**preserved**. When the deferred monsoon solver lands, it must **layer its
per-season overlay on top of** this current-modified annual baseline, not overwrite
the ocean seed. This is documented here and will be noted in `ENGINEERING-NOTES.md`
so the future work inherits the seam cleanly.

---

## 8. Determinism, the `currentStrength = 0` escape hatch, and D1 re-verification

**`currentStrength = 0` is an early return of the entire stage** — NOT a
multiply-by-zero. When `currentStrength === 0`, `computeOceanCurrents` /
`computeSstAnomaly` are not called, the moisture seed stays the literal `1.0`
assignment it is today, and the temperature loop adds nothing. This makes the
world **byte-identical to pre-D2** and provable by the existing `worldGen`
determinism test, exactly as D5's `starClass === 'G'` short-circuit did. (A
`1.0 + evap·0` seed happens to equal `1.0` in IEEE-754, but a zero-angle Rodrigues
rotation accumulated through floats is not guaranteed bit-identical — so the
short-circuit is the correct mechanism, not reliance on ×0.)

**Determinism of the active path.** Both solvers are `cells.forEach` + preallocated
`Float32Array` scratch, fixed pass counts, fixed accumulation order, **no RNG** —
the same construction the 8-pass moisture loop already relies on for byte-identical
reruns. A `tests/currents.test.ts` case asserts run-to-run identity of the velocity
+ anomaly fields.

**D1 checksum re-verification protocol (in this order):**

1. At **`currentStrength = 0`**, re-run the exact D1 browser check (5k cells,
   Mercator, biome view, season 0.5) and confirm the pixel checksum is
   **909197 exactly**. This is the end-to-end proof the escape hatch is real —
   stronger than the unit test. (Do NOT pre-declare 909197 "obsolete" — prove it
   survives at 0 first.)
2. At **default `currentStrength = 1.0`**, capture the **new** neutral baseline
   checksum, then drive the season slider off-neutral and back and confirm it
   returns to that new baseline **exactly** (neutral == canonical + return-to-exact
   against the fresh baseline).

---

## 9. Parameter surface

New `WorldParams.currentStrength: number` — range **0–2**, **default 1.0**,
**default-on**. (Matt authorized milestone-D determinism breaks; S14 precedent.)
`0` = disabled/byte-identical escape hatch. One knob scales Coriolis magnitude,
temperature coupling, and evaporation coupling.

Wiring (all required for param-liveness + old-save back-compat):

- `types.ts` — add to `WorldParams` (Climate section) with an inline-comment doc.
- `utils/defaultParams.ts` — `currentStrength: 1.0` (single source; `makeParams`
  inherits it free via the S15 spread).
- `components/Controls.tsx` — slider in the **Climate** tab; add to the
  auto-update regen dependency list (it's a generation param).
- `utils/export.ts` — `validateWorldParams` bound `[0, 2]`; `withParamDefaults`
  defaults missing `currentStrength` to `1.0` (old-save back-compat, pattern of
  `season`/`starClass`).
- `tests/paramLiveness.test.ts` — add a `currentStrength` perturbation (e.g.
  `0 vs 2`) proving it changes the terrain/climate signature.
- `services/gemini.ts` — surface in the lore prompt (pattern of `starClass`).

---

## 10. Testing & re-baseline

**New `tests/currents.test.ts` (pure-function):**
1. **Determinism** — same seed → identical velocity + `sstAnomaly`.
2. **Poleward → warm** — a constructed/seed config where a poleward current
   produces a positive high-latitude anomaly; equatorward → negative.
3. **`currentStrength = 0` no-op** — full `generateWorld` byte-identical to a
   pre-D2 baseline digest (or asserted via the determinism harness).
4. **Boundary tangency** — an ocean cell adjacent to land has ~no velocity
   component pointing into the land neighbour after the solve.

**Re-baseline sweep (explicit plan tasks — do NOT loosen assertions, replace the
seed/constant per S9/S12/S14 precedent).** Temperature *and* moisture both move, so
expect several fixtures to need rescan:
- `tests/lakes.test.ts` — `SALT_SEED` (s7) / `FRESH_SEED` (lakeworld): rescan for a
  salt-endorheic + a fresh replacement. **Delegatable** (see §11).
- `tests/routes.test.ts` — `SEA_SEED` (islands): confirm or replace.
- `tests/biomes.test.ts`, `tests/features.test.ts`, `tests/religions.test.ts` —
  run, re-baseline any pinned-seed expectations that shift.
- `tests/paramLiveness.test.ts` — new `currentStrength` case (above).

**Gates:** typecheck 0, lint 0 errors / ≤29 warnings, full vitest suite green
(honor the S10/S11 M1 parallel-load paramLiveness flake — re-run in isolation
before believing a failure), build OK, worker chunk unchanged (currents.ts must
not pull THREE/d3 into the worker — verify chunk size).

**Performance.** Runs in the worker. `N_VEL + N_HEAT ≈ 24` passes over ocean cells
only (~40–60% of cells); at default 5k this is <~100ms, dwarfed by the tectonic
sim. O(n·passes) scales fine to 200k. Measure and record.

---

## 11. Delegation plan

- **Self (not delegated):** `utils/currents.ts` (the solver — its correctness *is*
  the work) and the `worldGen.ts` climate-block wiring (a serial chain through one
  file). Param wiring is also self (small, spread across interacting files).
- **Delegatable — the seed re-baseline rescan.** Mechanical, high-volume,
  judgment-free once the predicate is specified. Brief must carry the **S14
  method**: the scan writes results to a JSON scratchpad via `fs.writeFileSync` and
  reads that back — **never trust piped vitest console** (S14 documents a piped scan
  returning zero output). And the **S12 rule**: run **synchronous foreground,
  per-file, never backgrounded** — a subagent cannot be resumed by a background
  notification. Tier: `sonnet-low` (mechanical brief-execution). Sequential if it
  shares any test file with other work.

---

## 12. Risk register

| Risk | Mitigation |
|---|---|
| Gyres don't visibly form / look wrong | Render-verify on `realmgenesis` before finalizing constants (§8); tune `DRAG`/`CORIOLIS_K`/`MIX`. If emergent gyres fail, fall back to a stronger boundary-tangency weight — still fixed-pass. |
| Determinism break in the active path | forEach + Float32Array scratch, fixed passes, no RNG; `currents.test.ts` determinism case. |
| Escape hatch not truly byte-identical | Early-return short-circuit (§8), proven by determinism test + the `currentStrength=0 → 909197` browser check. |
| D1 invariant silently broken | Two-step re-verification protocol (§8), browser checksum, run at 0 first. |
| Worker bloat (THREE/d3 leak) | `currents.ts` imports types + math only; verify worker chunk size unchanged. |
| Re-baseline larger than expected | Blast radius enumerated (§10); rescan delegated with S14/S12-safe brief. |
| Monsoon work forced to rewrite moisture seed | Composability seam documented (§7) + ENGINEERING-NOTES. |
