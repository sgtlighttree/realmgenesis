# Parchment art direction — labels, ink, furniture, projections

**Branch** `a3-map-style` (pushed, not merged). Follows S27e, which fixed the
globe texture's `flipY` inversion. Decisions below are Matt's, taken 2026-08-29.

## Decisions

- **Type:** `@fontsource/cinzel` for faction and province names (engraved Roman
  caps), `@fontsource/im-fell-english` italic for water and geographic names.
  Follows the cartographic rule — roman for land, italic for water — executed
  rustically. Rejected: all-IM-Fell (most rustic, worst legibility at 5k labels)
  and EB Garamond (most legible, reads as "a serif font" not as art direction).
- **Scope:** ink pass **and** map furniture.
- **Furniture is 2D ONLY.** Matt: "Parchment on 3D is more of a novelty more
  than anything." A compass rose and a neatline on a sphere are meaningless, and
  every furniture element is one more place 2D and 3D can diverge.
- **Also add Equirectangular and Winkel Tripel as 2D view projections.**

## W1 — Type

1. `npm i @fontsource/cinzel @fontsource/im-fell-english` (both 5.3.0).
2. Import only the faces used. Cinzel 400/700 latin; IM Fell English 400 italic.
3. `utils/labels.ts:262` is the ONE place a label font is set
   (`Inter, ui-sans-serif...`). It becomes style-aware, fed from the palette.
4. `LABEL_CONFIG` currently hardcodes `fill` (`#dbeafe`, `#fde68a`). Those move
   to the palette too, or parchment labels stay pale blue.

**THE TRAP, and it is the one this repo keeps falling into.** Canvas2D does not
wait for webfonts. If `ctx.font` names a family that has not finished loading,
the browser **silently falls back to sans-serif and draws anyway** — the map
looks exactly as it does today, and nothing errors. So:

- Await `document.fonts.load('700 14px Cinzel')` and the italic face explicitly.
- Re-render once they resolve; a map drawn before the fonts land keeps the
  fallback until something else invalidates it.
- Gate it in one hook, not per call site.
- Verify by screenshot, not by "the import is there".

**SVG export is a separate problem.** `exportSVG` emits `font-family` by name;
the file only renders correctly on a machine that HAS the font. Either embed the
woff2 as a base64 `@font-face` in `<defs>`, or accept the limitation and say so
in the export UI. Decide before shipping — do not leave it silent.

## W2 — Ink pass (the S27d finding)

`StylePalette` gains overlay entries: river, border, road, searoute, graticule,
contour, labelInk, labelHalo, marker, ruler.

One accessor — `overlayInk(mapStyleId)` — consumed by every site. Today the river
colour alone is hardcoded in **four files and five call sites**:

| File | Sites |
|---|---|
| `components/overlays/tenants.ts` | rivers, borders, roads, sea routes, graticule, contours, ruler, dymaxion cage |
| `components/Map2D.tsx` | rivers at `:480` AND `:717`, `drawFactionBorders` at `:76` |
| `utils/labels.ts` | label fill and halo |
| `utils/export.ts` `:353`, `utils/exportVector.ts` `:138` | borders, rivers |

**Collapse the duplicates while threading the style.** If the five river sites
survive this pass, the next style diverges in exactly the same way — that is the
whole lesson of the finding.

## W3 — Map furniture (2D only)

Compass rose, ruled neatline, edge vignette. These are canvas furniture, not
per-cell passes, so they need an explicit `showFurniture?: boolean` on
`StyleRenderContext`.

**Do NOT key off `wrapsHorizontally`** just because the globe bake happens to be
the only caller that sets it. Two unrelated meanings on one flag is how the next
reader gets it wrong. Set `showFurniture` from Map2D and the export paths; leave
it unset in `bakeStyleTexture`.

Furniture must sit OUTSIDE the horizontal mirror the 2D paths apply
(`translate(w,0); scale(-1,1)`), or the compass rose reads backwards — the same
trap glyphs already documented.

## W4 — Equirectangular + Winkel Tripel on screen

The export path ALREADY supports equirectangular, mercator, winkeltripel,
robinson, mollweide and orthographic (`utils/exportVector.ts` `buildProjection`,
via `d3-geo-projection`, already a dependency). Only the on-screen `Map2D` is
limited to Mercator plus Dymaxion.

1. Extract `buildProjection` into `utils/projections.ts` so screen and export
   share one factory. Two registries would drift.
2. `DisplayMode` (`types.ts:313`) gains `'equirectangular' | 'winkeltripel'`.
3. `Map2D`'s projection memo uses the factory. Pan and zoom are projection
   agnostic.
4. Toolbar gains the two buttons; the footer hint stops hardcoding
   "2D Mercator" (`Map2D.tsx:1244`).

**Verified up front:** `geoWinkel3().invert()` exists and round-trips
(37.5, −12.25 → exact), so cell picking reuses the Mercator path unchanged. No
pick buffer needed — that is a Dymaxion-only requirement.

## Order

W4 is independent; W1 → W2 → W3 are sequential. Ship W1 first — it is the most
visible and it carries the font-loading trap that must be proved by screenshot.

## Verification

`scripts/renderMapPreview.mjs --style=parchment --png` for the flat map, and the
running app through Playwright for the fonts, since the fallback is invisible to
any headless SVG check. Per S27e: plausibility is not a check. Compare against
ground truth, and look at the actual pixels.
