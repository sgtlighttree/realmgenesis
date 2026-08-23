# D8b — Datum simulation coupling (`physicalClimate`)

**Status:** design approved (S25b, 2026-08-23), implementation pending.
**Author:** agent + Matt (decisions), advisor-reviewed.
**Depends on:** D8a datum + the `frac^2` hypsometric curve (both shipped S25/S25b,
`utils/datum.ts`). The GRASSLAND biome (shipped S25b) receives cells that grow
wetter under the moisture retune here.

---

## 1. Purpose

Three tuned magic constants in generation are physical quantities wearing
normalized clothes. D8a gave heights a real datum (metres); D8b feeds that datum
back into the simulation so the constants become grounded physics:

1. **Lapse rate** — `temp -= elevation * 60` (`worldGen.ts:622`). The `60` is an
   invented cooling rate. In metres it should be ~6.5 °C/km.
2. **Orographic precipitation** — the 8-pass moisture transport keys rain shadow
   off a normalized `heightDiff` threshold (`worldGen.ts:582-584`). Rain-shadow
   strength physically depends on barrier height in metres.
3. **Snow line** — currently emergent from temperature only; a grounded lapse rate
   makes it vary with altitude AND latitude for free (no `determineBiome` change).

The moisture retune also fixes a measured realism defect (see §6): 42–45 % of land
sits at moisture < 0.15 (Earth arid+semiarid ≈ 33 %), which is the real cause of
steppe/dry-biome dominance — not elevation.

---

## 2. The escape hatch (non-negotiable)

A single boolean param **`physicalClimate: boolean`**, default **`true`**.

- `true` (default) → grounded physics (this spec).
- `false` → the EXACT current formulas, byte-identical to the S25b engine.

Naming: `physicalClimate`, NOT `datumCoupling` — every other hatch in this project
is named for the effect it controls (`currentStrength`, `starClass`,
`seafloorDepth`), not the mechanism.

**Default-ON is a deliberate, authorized choice (Matt, S25b).** It changes output
for existing seeds — including **civ layout** (biome → settlement suitability →
`recalculateCivs` → factions/towns/populations), not just climate tint. ROADMAP:34
authorizes seed breakage; Matt explicitly accepted the civ-layer movement.

**Byte-identical baseline is CURRENT `main` (post-curve, post-grassland), NOT
ancient main.** The grassland biome and the hypsometric curve were separate
deliberate changes already in `main`. `physicalClimate=false` reproduces the
engine as it stands after S25b, nothing older.

---

## 3. Architecture

D8b touches exactly **two code sites**, both in `utils/worldGen.ts`, both gated by
`params.physicalClimate`. `determineBiome` is **NOT modified** — the snow-line
change is emergent through temperature.

```
generateWorld
 ├─ 8-pass moisture transport loop     ← SITE 2 (orographic), gated
 └─ per-cell climate finalize
      ├─ temp -= lapse                  ← SITE 1 (lapse rate), gated
      └─ c.biome = determineBiome(...)  ← UNCHANGED (snow line emergent)
```

New physical constants live in `utils/datum.ts` alongside the datum they depend on
(single source of truth), imported by `worldGen.ts`.

---

## 4. Site 1 — Lapse rate

**Current (`worldGen.ts:621-622`):**
```ts
const elevation = Math.max(0, c.height - params.seaLevel);
temp -= elevation * 60;
```

**Grounded:**
```ts
if (params.physicalClimate) {
  const elevM = Math.max(0, elevationMetres(c.height, params.seaLevel, params.maxElevationM));
  temp -= (elevM / 1000) * LAPSE_RATE_C_PER_KM;   // LAPSE_RATE_C_PER_KM = 6.5
} else {
  const elevation = Math.max(0, c.height - params.seaLevel);
  temp -= elevation * 60;                          // exact old formula
}
```

`elevationMetres` already applies the `frac^2` curve, so most land stays low and
only genuine peaks cool hard. **This is why the accurate datum pick works** —
measured (curve + 6.5 + datum 9000): ICE_CAP 7.0→6.7 %, biome change ~21–27 %,
vs 11→28 % / 15 % new-ice under the OLD linear datum. The curve dissolved the
former pick-two tradeoff; take physical lapse 6.5 + datum 9000 + curve, all three.

`LAPSE_RATE_C_PER_KM = 6.5` — the standard environmental lapse rate. New constant
in `utils/datum.ts`.

**Note:** the `~21–27 %` figure is a FLOOR. It was measured lapse-only (post-hoc
temperature substitution). The full engine also changes moisture (Site 2), which
`determineBiome` reads, so the real shift is larger. Re-measure with the tool after
implementation; do not treat 21–27 % as the prediction.

---

## 5. Site 2 — Orographic precipitation

**Current (`worldGen.ts:582-584`), inside the `if (dot > 0)` windward branch:**
```ts
const heightDiff = c.height - n.height;
if (heightDiff > 0.02) carry *= 1.5;       // windward (air rising) — more rain
else if (heightDiff < -0.02) carry *= 0.2; // leeward (descending) — rain shadow
```

**Grounded:** scale by the real barrier metres crossed, land-side only.

```ts
if (params.physicalClimate) {
  // Land-side above-sea metres only, to avoid the seaLevel slope kink. An ocean
  // neighbour contributes 0 (elevationMetres clamps its land term at coast).
  const cM = Math.max(0, elevationMetres(c.height, params.seaLevel, params.maxElevationM));
  const nM = Math.max(0, elevationMetres(n.height, params.seaLevel, params.maxElevationM));
  const dKm = (cM - nM) / 1000;
  if (dKm > 0)      carry *= 1 + OROG_WINDWARD_PER_KM * dKm;                  // capped
  else if (dKm < 0) carry *= Math.max(OROG_LEEWARD_FLOOR, 1 + OROG_LEEWARD_PER_KM * dKm);
} else {
  const heightDiff = c.height - n.height;
  if (heightDiff > 0.02) carry *= 1.5;
  else if (heightDiff < -0.02) carry *= 0.2;   // exact old formula
}
```

**Advisor blocker resolved (item 3):** use `max(0, elevationMetres(...))` per cell
(land-side, above-sea) rather than an `elevationMetres` DIFFERENCE that would cross
the seaLevel kink where the land slope (`maxElevationM/(1−sea)`) meets the ocean
slope (`MAX_DEPTH_M/sea`). Ocean neighbours read 0, so a coastal barrier is
measured from sea level, which is correct.

**Constants (tuned to a TARGET, not prescribed):**
- `OROG_WINDWARD_PER_KM`, `OROG_LEEWARD_PER_KM`, `OROG_LEEWARD_FLOOR` (and a
  windward cap) live in `utils/datum.ts`.
- **Do NOT guess these.** Tune them against `scripts/queryWorld.mjs climate` until
  the land moisture distribution matches the target in §6. The current `1.5 / 0.2`
  pair produced 42–45 % of land at moisture < 0.15; the retune must bring that
  toward Earth's ~33 % arid+semiarid without over-wetting coasts.

**Air-temperature factor: DEFERRED.** ROADMAP notes rain-shadow strength also
depends on air temperature, but `temp` is finalized AFTER the moisture loop
(Site 1 runs later). Using it needs a reorder or a pre-pass temperature estimate.
Barrier height is the dominant term; ship it first, note the temp factor as a
follow-on that requires the reorder.

---

## 6. Moisture retune target (folded into Site 2)

**Measured problem** (`queryWorld.mjs climate`, 3 seeds): median land moisture ≈
0.22, bimodal — coasts soak (p75 ≈ 0.78), interiors starve (p25 ≈ 0.086), 42–45 %
of land < 0.15. The 8-pass transport under-delivers inland; the harsh `*0.2`
leeward step and cumulative `*0.98` land decay strand interiors.

**Target:** land moisture < 0.15 share drops from ~42 % toward Earth-like ~33 %
arid+semiarid, WITHOUT flipping the map to all-forest. Grassland (added S25b) and
temperate forest should grow as interiors wet; steppe falls to a realistic minority.

**Levers (in priority order):**
1. Soften the leeward floor (`*0.2` loses 80 %/step — too harsh for a single edge).
2. Consider the land decay (`worldGen.ts:592` `*= 0.98`) if interiors still starve.
3. Windward boost stays bounded so coasts don't saturate.

**Verify by measuring, not eyeballing** — `queryWorld.mjs climate` + `biomes`
before/after, 3+ seeds. This is the one part of D8b with no clean closed form; it
is empirical tuning against a measured target.

---

## 7. Snow line — emergent, no code

`determineBiome` (`worldGen.ts:354`) decides ICE_CAP/TUNDRA from temperature. With
grounded lapse (Site 1), temperature already varies with altitude and latitude, so
the snow line becomes an altitude-and-latitude surface for free. **No
`determineBiome` change.** When `physicalClimate=false`, temperature uses the old
lapse, so biomes are unchanged — the hatch covers the snow line automatically.

**Known consequence:** colder peaks nearly eliminate VOLCANIC (`landH>0.85 &&
temp>-5`) — measured VOL → 0 %. A snow-capped volcano is realistic, so this is
acceptable. `landH>0.85` is already datum-relative. Decoupling volcanic from
temperature (it is a rock type, not a climate) is a SEPARATE biome decision, out of
D8b scope — note it, do not fix it here.

---

## 8. Plumbing

- **`types.ts`:** add `physicalClimate: boolean` to `WorldParams` (near `season`).
- **`utils/defaultParams.ts`:** `physicalClimate: true`.
- **`utils/datum.ts`:** add `LAPSE_RATE_C_PER_KM = 6.5`, `OROG_*` constants.
- **`utils/export.ts`:**
  - `withParamDefaults`: default missing `physicalClimate` to `true` (old saves get
    grounded climate on load — the accepted civ-move; reproducible with `false`).
  - `validateWorldParams`: a BOOLEAN type-check line alongside `mapName`/`seed`,
    NOT an entry in `numericBounds` (advisor item 5).
- **`components/Controls.tsx`:** a toggle in a climate/presentation section (Clim
  tab is natural). Warning text: grounded physics changes generated climate AND
  civ layout for existing worlds; off reproduces current output.
- **`tests/paramLiveness.test.ts`:**
  - `maxElevationM` moves from the display-only allowlist INTO a live perturbation
    — default-ON makes it a genuine generation input (it drives lapse + orographic).
    Add it to `TERRAIN_PERTURBATIONS`.
  - Add `physicalClimate: false` to `TERRAIN_PERTURBATIONS` (toggling it must change
    the terrain signature).

---

## 9. Testing

1. **`tests/datum.test.ts`** — unit-test the new constants exist and the grounded
   lapse/orographic helper formulas (if extracted) compute expected metres.
2. **Byte-identical hatch test** — generate with `physicalClimate:false` and assert
   the terrain+biome signature equals the S25b baseline. Capture the baseline from
   current `main` BEFORE implementing (worktree at the pre-D8b SHA, or
   `scripts/captureBaseline.mjs before`). This is the D2/D5 discipline.
3. **paramLiveness** — as §8 (maxElevationM live, physicalClimate live).
4. **Empirical moisture check** — `queryWorld.mjs climate`/`biomes`, 3+ seeds,
   before/after; record the numbers in HANDOFF. Not a pass/fail unit test — a tuning
   gate with a target (§6).
5. **Browser** — grounded world renders, no console errors; cold high peaks,
   wetter interiors, more grassland/forest, less steppe. Use the S25b headless
   Playwright recipe (HANDOFF); reuse Matt's :3000; kill chromium after.

---

## 10. Risks

- **Moisture tuning is the only open variable.** Bad constants destabilize biomes
  and civs across 8 coupled passes. Mitigate by tuning against the measured target,
  not by intuition; keep the hatch so a bad tune is never forced on anyone.
- **The 21–27 % biome-change floor will grow** once moisture moves (§4 note). Set
  expectations: existing seeds will look meaningfully different. That is authorized.
- **Worker chunk changes** — these are generation-stage edits, so `worldGen.worker`
  rebuilds. Expected, not a regression (unlike the render-only S16–S24 work).

---

## 11. Out of scope (decide, don't drift)

- Air-temperature factor in orographic (needs a loop reorder) — §5.
- Volcanic decoupling from temperature — §7.
- D5 gravity (isostasy cap on peak height) and D3 sea-level coupling — D8b unblocks
  them per ROADMAP, but they are their own features.
- Lore-prompt elevation text (deferred from D8a).
