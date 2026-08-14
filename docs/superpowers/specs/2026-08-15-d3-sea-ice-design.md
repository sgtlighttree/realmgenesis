# D3 — Sea-ice & glaciers (design)

**Date:** 2026-08-15 · **Milestone:** ROADMAP D · **Status:** implemented

Seasonal sea-ice on the ocean, building directly on D1's seasonal temperature.

## Scope (deliberately narrow)

- **Sea-ice only.** Open-ocean cells whose *seasonal* temperature is below the
  seawater freezing point render as ice. Land ice is **left to the existing**
  `ICE_CAP` biome + snow-on-elevation logic — no second whitening path (a third
  overlay interacting with both muddies the result).
- **Render overlay, not a stored biome.** `cell.biome` stays `OCEAN`/`DEEP_OCEAN`;
  ice is a color decision in `getCellColor`. No new biome, no civ/navigation
  impact, no generation change. Mirrors D1's "display shifts, canonical frozen".
- **No parameter.** Seawater freezes at −2°C — a physical constant
  (`SEAWATER_FREEZE_C` in `utils/seasons.ts`), not a creative knob. The iciness
  levers already exist and are live/regenerating: `poleTemperature`,
  `baseTemperature`, `axialTilt`. (A freeze-point slider would only let the user
  lie about physics, and every render-only param costs a paramLiveness exclusion +
  defaults + validation + save-schema entry.)

## Census (the go/no-go check, run before writing color code)

`generateWorld(makeParams())` at defaults (`poleTemperature −30`, `tempVariance 5`),
237 ocean cells: min −26.9°C, max 32.1°C, **35 cells (14.8%) below −2°C**, 38 below
0°C. So sea-ice covers ~15% polar ocean at defaults — a visible win, neither absent
nor swamping. The −2..0°C edge band is ~3 cells; `temperatureVariance` (±5°C) makes
the edge mildly ragged, which is physically fine (real sea-ice edges are ragged) and
finer at higher cell counts.

## Implementation

In `getCellColor`, before the view-mode switch:

```
if ((mode === 'satellite' || mode === 'biome')     // physical views only
    && cell.height < seaLevel                       // open water
    && biome !== LAKE && biome !== SALT_LAKE         // lakes read as lakes
    && seasonalTemp < SEAWATER_FREEZE_C) {           // seasonal, from D1
  iciness = clamp01((FREEZE − seasonalTemp) / 15);   // colder → whiter
  color = lerp(#bcd4e6 /*edge*/, #ffffff /*core*/, iciness);
  return;
}
```

- **`satellite` + `biome` only.** Explicitly **not** `height`/`height_bw`
  (elevation data views — ice would corrupt the read) or `temperature` (already
  encodes the datum as color). Do not "fix" this later.
- `seasonalTemp = cell.temperature + seasonalDelta`, so ice extent is annual-mean
  at neutral season and migrates hemisphere-to-hemisphere as the slider moves.

## Testing

- Pure-function test (`tests/seasons.test.ts`): a mid-high latitude exists whose
  shown temperature is above freeze in its hemisphere's summer and below in winter
  → the ice edge migrates. Plus a north-hemisphere winter-colder-than-summer check.
- In-browser (Playwright checksum method, 5k cells): at neutral, both poles render
  full sea-ice (0 blue ocean in the polar bands), equator stays blue. At season
  0.25 (N summer) the north cap melts to open water (ice 17.3k→0.5k, ocean
  0→16.8k) while the south stays frozen. 0 console errors.

## Deferred (ENGINEERING-NOTES)

- **Sea-level coupling** (more ice → lower sea level): changes the coastline, so it
  is a generation-stage / regeneration concern, not a render overlay. Out of scope
  for the render-only tier.
- **Glacial land cover beyond ICE_CAP/snow:** revisit only if polar *lowlands* read
  wrong once sea-ice lands — a separate, better-informed change.
