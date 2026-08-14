# D1 — Seasonal Cycle (design)

**Date:** 2026-08-14 · **Milestone:** ROADMAP D (planet-scale simulation) · **Status:** approved, implementing

A season slider that shifts climate through the orbital year. `axialTilt` currently
bakes a static, permanent offset; D1 reinterprets it as the amplitude of a seasonal
excursion and lets temperature (and derived ice/biome edges) move as the subsolar
point migrates between the hemispheres.

## Decisions (settled in brainstorm + advisor review)

1. **Compute = live formula, not precomputed snapshots.** The seasonal excursion is
   closed-form in latitude and params; it is derived in the render/color layer from the
   slider, not baked.
2. **Scope = temperature + ice extent only.** Wind and the 8-pass moisture solver stay
   annual (no per-season rerun). The richer "wind/moisture shift → monsoons" variant is
   deferred and recorded in `docs/ENGINEERING-NOTES.md`.
3. **Accept the determinism break.** Reinterpreting `axialTilt` changes existing seeds'
   baseline climate. The maintainer explicitly accepted this while milestone D is in
   flight — do not preserve byte-compat with pre-D1 seeds.
4. **Option B — biomes shift for display, civs frozen.** Canonical `cell.biome` stays the
   annual-mean biome (drives civs, export, labels). The *displayed* biome at a non-neutral
   season is re-derived live. Civ borders never move with the season.

## Model

No new per-cell field is stored. The temperature formula is closed-form in
latitude/height/params, and the season-independent terms (`−elevation·60`, the simplex
noise) are additive, so the seasonal excursion is a pure function of latitude and params:

```
φ        = asin(clamp(cell.center.y, -1, 1))        // geometric latitude, radians
Tlat(ψ)  = baseTemperature·(1 − r²) + poleTemperature·r²,   r = |ψ| / (π/2)
δ(s)     = axialTilt_rad · sin(2π·s)                 // subsolar declination, s ∈ [0,1)
Tmean_lat(φ) = ⟨ Tlat(φ − δ(s)) ⟩ over one orbit     // 12–24 sample average
ΔT(φ, s) = Tlat(φ − δ(s)) − Tmean_lat(φ)             // seasonal excursion, mean 0 by construction
```

- **Stored (generation):** `cell.temperature = Tmean_lat(φ) − elevation·60 + noise·temperatureVariance`.
  This is the orbit-averaged annual mean. It keeps `axialTilt` affecting generated
  output — `Tlat` is quadratic, so by Jensen's inequality the orbital average shifts with
  tilt (well clear of float noise at the paramLiveness 60° perturbation). This is the fix
  for the "`axialTilt` goes dead → paramLiveness fails" trap (the `roughness`/S9 failure mode).
- **Displayed (render):** shown temperature = `cell.temperature + ΔT(φ, s)`.
- **No** `worldTransfer.ts` change, **no** serialization/save-schema bump, **no** `_season`
  RNG side-stream — none are needed because nothing new is persisted or randomly drawn.

**Wind coherence fix.** The wind-band block currently uses the *tilted-axis* latitude.
Move it to geometric latitude (`asin(cell.center.y)`) so winds and temperature share one
axis. Winds remain annual/static per the scope decision — this is purely a correctness fix,
one line.

## New parameter

`season: number` — orbital position, `0 ≤ season < 1`. **`0.5` is the neutral point**
(δ = `axialTilt·sin(π)` = 0 → ΔT ≡ 0 → shown = annual mean). Chosen so the default world
looks like the annual mean. Default `0.5`. Not a terrain param — it must NOT enter
`generateWorld` output (it is render-only), so it is **excluded** from paramLiveness terrain
perturbations and does not trigger regeneration.

- `types.ts`: add to `WorldParams`.
- `hooks/useWorldEngine.ts`: `DEFAULT_PARAMS.season = 0.5`.
- `tests/helpers.ts`: add `season: 0.5`.
- `utils/export.ts` `validateWorldParams`: bound `season: [0, 1]`.

## Render threading

`getCellColor` reads `cell.temperature` at two points: the snow threshold
(`colors.ts:122`, `> 20`) and the temperature-viewmode normalization (`colors.ts:164`).
Both must use the seasonal value. Approach: add a trailing **optional `seasonalDelta?: number`**
that `getCellColor` adds to `cell.temperature` before those two reads (default 0 → neutral,
byte-identical to today). A trailing optional (over a 7th positional) is chosen because
`getCellColor` is already 6-positional and `DymaxionPreview2D` calls it with only 3 args.

Each render path computes the per-cell delta once via a shared pure helper
`seasonalTemperatureDelta(cell, params)` (returns 0 when `season === 0.5`) and passes it.
**9 call sites:** `Map2D.tsx` (×2), `WorldViewer.tsx`, `MiniMap.tsx`, `DymaxionPreview2D.tsx`,
`utils/export.ts` (×2), `utils/exportVector.ts`, `utils/exportGLB.ts`.

- **Displayed biome:** in `Map2D`/`WorldViewer` biome + satellite modes, the color path
  derives the seasonal biome from the seasonal temperature for display only. `cell.biome`
  is never mutated.
- **Snow caps** advance/retreat with season for free via the threshold read.

## Divergence decisions (view vs canonical data)

- **Inspector** (`Inspector.tsx:264,283`): when `season !== 0.5`, show both the annual
  `cell.temperature`/`cell.biome` and the seasonal values. Prevents "map says X, panel says Y".
- **Export** (PNG/SVG/GeoJSON/GLB): render **as-displayed** at the current season, so an
  export always matches the screen. Neutral (0.5) by default.
- **Paint** (`paintUtils.ts:121`): stays **canonical** — `determineBiome` keeps reading
  `cell.temperature`. Never season-aware. Stated here so it is not "fixed" later.

## UI

Season slider in `Controls.tsx` (Climate/Geo group): `0–1`, neutral at `0.5`, labelled by
season quarter (e.g. 0.5 = "equinox baseline", 0.75/0.25 = solstices) rather than raw float.
Moving it recolors only — it must **not** schedule regeneration (it is not in the
generation-affecting param set of the auto-update effect).

## Testing

- `tests/seasons.test.ts` (new): `seasonalTemperatureDelta` — (a) orbital mean ≈ 0 across a
  full `s` sweep; (b) hemisphere sign flip (north warms while south cools at a given `s`);
  (c) `axialTilt = 0` ⇒ ΔT ≡ 0 for all `s`; (d) neutral `season = 0.5` ⇒ ΔT = 0.
- `tests/paramLiveness.test.ts`: `axialTilt` must STILL be live after the orbit-average
  change — run `npx vitest run tests/paramLiveness.test.ts` and confirm the terrain
  signature responds to the 60° perturbation. `season` is NOT added to terrain perturbations.
- Full suite + typecheck + lint + build gates.

## Decomposition

1. **Engine** (serial, first): orbit-averaged `Tmean`, untilt wind block, `seasonalTemperatureDelta`
   helper, `season` param plumbing. Judgment-heavy — done directly, not delegated.
2. **Render threading** (after 1): `getCellColor` seasonal arg + 9 call sites + displayed biome.
3. **UI** (parallel with 2): season slider + `useWorldEngine`/Controls wiring.
4. **Tests** (after 1): `seasons.test.ts` + paramLiveness confirm.
5. **ENGINEERING-NOTES**: deferred wind/moisture-monsoon variant entry (independent).

Small feature: implemented mostly directly with the test suite as the gate, rather than
fanned out to subagents (round-trip cost would exceed the savings on ~9 one-line edits).

## Deferred (to ENGINEERING-NOTES)

- **Seasonal wind + moisture (monsoons).** Rerun the wind-band + 8-pass moisture solver per
  season so wet/dry seasons and monsoon reversals appear. Expensive (the solver is the costly
  part of climate) and it would make biome-at-season depend on a per-season moisture field
  rather than the annual one — no longer a free O(n) recompute. Deferred deliberately.
