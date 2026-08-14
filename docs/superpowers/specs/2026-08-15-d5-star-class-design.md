# D5 — Planetary parameters: host star class (design)

**Date:** 2026-08-15 · **Milestone:** ROADMAP D · **Status:** implemented (partial — star class only)

## Scope decision: ship `starClass` only, defer the rest

D5 as listed is "day length, gravity, star class, moons — mostly lore/export
metadata at first." Building all four now would ship **hollow params into the save
schema permanently**. Only `starClass` has a real, cheap, *live* mechanical hook:

- **Star class** → scales global insolation → temperature → (via D1) biomes and
  sea-ice all follow for free. One enum, one multiplication, immediately live in
  `paramLiveness` with no allowlist entry. **Ship it.**
- **Day length** — no diurnal cycle in the model; the only real hook (Coriolis
  strength) is D2 territory. **Defer.**
- **Gravity** — there is no principled map from *g* to the normalized `[0,1]`
  height field. Wiring it into Stage-9b would be a fudge factor that duplicates
  `mountainHeight` — a param that *lies*. **Defer** (do not invent the mapping).
- **Moons** — no tide model. Pure string. **Defer.**

So ROADMAP D5 → 🟡 PARTIAL, exactly matching its own prediction ("mostly metadata at
first, but D1 and tides give them mechanical hooks over time"). Deferred set recorded
in ENGINEERING-NOTES.

## Model

`starClass: 'O'|'B'|'A'|'F'|'G'|'K'|'M'`, default `'G'` (Sun-like). Applied in the
climate pass (`utils/planetary.ts applyStarClass`), scaling the latitude temperature
**in Kelvin** (Stefan-Boltzmann style, T ∝ L^¼) before the elevation lapse:

```
temp = applyStarClass(annualMeanLatTemp(φ, params), starClass)
     = (Tlat + 273.15) · f(class) − 273.15
```

Kelvin, **not** a multiply on °C: a plain °C multiply would warp the *negative* pole
temperature the wrong way (a dimmer star would spuriously *warm* the poles). In
Kelvin a dimmer star cools the whole world and drives the poles further negative.

Factors `f` are a **stylized insolation range** (0.93 … 1.07), not literal stellar
luminosity — a real O star is ~10^5× solar and would vaporize any planet outside a
razor-thin habitable zone. The range is kept narrow so biomes stay varied at both
extremes. **`G = 1.0` is an exact identity** (short-circuited), so default worlds are
byte-identical to pre-D5.

## Verification

- **Determinism / liveness:** `paramLiveness` gains `starClass: {starClass: 'M'}`
  (changes the terrain signature → live); default G keeps the baseline. `worldGen`
  determinism + full suite green.
- **Biome-variety census** (the advisor's blocking check, 300-cell default seed):
  - G: 11 distinct land biomes, temp −26.9…32.1 (= pre-D5 baseline, confirms no-op).
  - O: 8 biomes, −9.8…53.1 (hot/tropical world — rainforest/desert/savanna).
  - M: 6 biomes, −44…11 (frozen world — tundra/ice-cap/boreal).
  Both extremes stay varied (no monoculture), so the range needs no narrowing.

## Surfacing

- **Controls** (Climate tab): a `Select` for host star class (labels like
  "K — orange dwarf"). It's a generation param → added to the auto-update
  regeneration dep list.
- **Lore** (`services/gemini.ts`): the prompt now states "Orbits a `<class>`-class star."
- **Save/back-compat:** `validateWorldParams` accepts the 7 classes; `withParamDefaults`
  defaults a missing `starClass` to `'G'` (pre-D5 saves).

## Deferred (ENGINEERING-NOTES)

Day length, gravity, moons — no principled mechanical hook yet. Revisit when D2
(currents/Coriolis) gives day length a hook, and if a tide model ever gives moons one.
Gravity specifically must **not** be wired as a relief fudge factor.
