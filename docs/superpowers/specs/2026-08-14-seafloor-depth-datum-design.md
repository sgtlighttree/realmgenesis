# Seafloor Depth datum — design

**Date:** 2026-08-14
**Status:** approved (Matt, inline)
**Scope:** repurpose the "Seafloor Detail" slider into a linear ocean-floor depth datum.

## Goal

Matt's HANDOFF note: *"Seafloor Detail slider function more like a sea level
controller."* Clarified to: a knob that **raises/lowers the whole ocean floor**
(mean water depth) **without moving the coastline** `seaLevel`.

`oceanDepth` already exists and holds the coastline fixed while reshaping
bathymetry — but it is a **power curve** (contrast: trenches disproportionately
deeper than shelves). The new knob is a **linear scale** (datum: whole floor up
or down, relative shape preserved). The two are complementary and both kept.

## The knob

New param **`seafloorDepth`**, range **0.3–2.0**, default **1.0**
(byte-identical to current behavior at default). Applied to every below-sea-level
cell in Stage 9b (`utils/worldGen.ts`), folded into the existing
`mountainHeight`/`oceanDepth` remap block, **after** final normalization:

```
water cell (h < seaLevel):
  t      = (seaLevel - h) / seaLevel          // normalized depth 0..1
  shaped = pow(t, 1/oceanDepth)               // existing power-curve (∈ [0,1])
  h'     = seaLevel - seaLevel * min(1, shaped * seafloorDepth)
```

- `seafloorDepth < 1` → floor rises toward `seaLevel` → shallower seas, less mean depth.
- `seafloorDepth > 1` → floor sinks → deep abyssal everywhere; deepest cells clamp at `h = 0`.
- `seafloorDepth == 1` → `min(1, shaped)` is a no-op (`shaped ∈ [0,1]`), so output is
  **byte-identical** to the current remap. Coastline (`t == 0 → h = seaLevel`) never moves.

The block's guard becomes `if (mh !== 1.0 || od !== 1.0 || sd !== 1.0)`.

## Retiring `seafloorDetail`

The old `seafloorDetail` param (0–1, default 0.5) did two internal jobs in
`utils/tectonicsV3.ts`: abyssal-hill noise amplitude (`* seafloorDetail * 0.06`)
and GDH1 noise-damping (`1 − 0.65·seafloorDetail`, the fable-advisor Session-10
correction that stops structural noise washing out the age→depth gradient).

**Both are baked at their current default (0.5)** as constants. Every default
world stays visually identical; the GDH1 protection is preserved. The param is
removed from `WorldParams`, defaults, and the slider.

- Precedent: `plateInfluence`→`tectonicStrength` (Session 8, old key silently dead — accepted).
- Old saves carrying `seafloorDetail` ignore it; missing `seafloorDepth` falls back to 1.0.

**Rejected alternative:** keep `seafloorDetail` as a hidden detail knob and add
`seafloorDepth` as a fourth ocean slider. Rejected — Matt flagged the Seafloor
Detail slider *specifically* as the thing to repurpose, and four ocean-shape
sliders (seaLevel / oceanDepth / seafloorDepth / seafloorDetail) is knob-soup.

## Edit sites

- `types.ts` — `seafloorDetail` → `seafloorDepth` (comment updated).
- `hooks/useWorldEngine.ts`, `tests/helpers.ts` — default `1.0`.
- `components/Controls.tsx` — slider relabel "Seafloor Depth", range 0.3–2.0, default 1.0.
- `utils/worldGen.ts` — Stage 9b remap: add `sd` guard + water-cell scale.
- `utils/tectonicsV3.ts` — bake the two `seafloorDetail` uses at constant 0.5.
- `tests/paramLiveness.test.ts` — swap the `seafloorDetail` live-case for
  `seafloorDepth: 0.5` (live because the default world has ocean at `seaLevel 0.55`).

## Not touched

`seaLevel`, `oceanDepth`, `mountainHeight`, GDH1/tectonicsV3 macro heights,
climate, erosion, biomes, rivers, civ.

## Verification

Four gates (typecheck / lint / test / build). paramLiveness proves the new param
is live and the retired one is gone. Byte-identical-at-default is guaranteed by
the `min(1, shaped)` no-op and the 0.5 bake, not just asserted.
