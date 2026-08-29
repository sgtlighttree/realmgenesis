# Map styles

_A3. Settled context — how the style layer works and how to add one._

A **style** changes how the 2D map and its exports are *drawn*. It is a separate
axis from `ViewMode`: `viewMode` decides what the map shows, a style decides how
it looks. Parchment + political is a valid, and intended, combination.

The 3D globe **is** styled, via a baked texture (see below). It was originally
out of scope — it paints per-cell vertex colours, so hatching, stroke weights and
paper grain have nowhere to live, and per-cell colour alone can only ever produce
a beige globe rather than a drawn map.

## The problem the architecture solves

Parchment has to render on two different surfaces:

- **Canvas2D** — `Map2D` (screen and the Dymaxion source buffer) and the PNG export
- **SVG string** — the vector export

The failure mode is writing the look twice and watching the two drift. Every
decision below serves avoiding that.

## Substrate

`utils/mapStyle/substrate.ts` defines one narrow drawing interface with two
implementations, `Canvas2DSubstrate` and `SvgSubstrate`. Style passes are written
once, against the interface.

The interface is deliberately narrow: a pass may only do what **both** surfaces
express natively. Anything raster-only would silently degrade the SVG export.

Three details are corrections, not preferences — do not "simplify" them back:

| Detail | Why |
|---|---|
| `fillFeature(feature, fill, opacity?)` takes opacity as a **parameter** | Canvas2D accepts `rgba()`; **SVG 1.1 has no `rgba()` colour syntax** and renders it inconsistently in Illustrator and Inkscape. Colour strings stay opaque. |
| Opacity applies to the seam **stroke** as well as the fill | The hairline stroke closes the antialiasing seam between adjacent Voronoi polygons. Without `stroke-opacity`, a 5%-opacity hillshade drew a full-strength outline on every cell and the map rendered as a dark honeycomb. |
| `hatchFeatures(features[], spec)` is **plural** | A world has 13,000–17,000 ocean cells and Map2D re-renders on every viewport change. Per-feature hatching meant a clip plus a full-canvas line sweep each. One composite region, one sweep. |

## Passes

A style is a palette plus an ordered list of passes (`utils/mapStyle/passes.ts`).
Parchment draws:

1. `paperPass` — paper tone plus grain
2. `oceanFillPass` — flat sea tone on ocean cells
3. `landPass` — obeys the fill policy below
4. `hillshadePass` — land only, tinted from the palette
5. `oceanHatchPass` — one composite ocean region, **after** the shading
6. `coastlinePass` — a single strong ink line
7. `glyphPass` — relief and vegetation marks

`MapStyle.passes` is also the single test for "does this style draw anything".
Render paths check `passes.length`, **never** the style id, so the invariant
lives in one place. The `default` style has an empty pass list and is a visual
no-op.

## Fill policy — the ViewMode collision

"Bare parchment, glyphs carry the terrain" collides with twelve view modes. The
rule is explicit in `styleParchment.ts`:

| Group | Modes | Treatment |
|---|---|---|
| Glyph-carried | `satellite`, `biome`, `height_bw` | Bare paper |
| Categorical | `political`, `province`, `culture`, `religion`, `plates` | Muted fill on paper |
| Continuous ramp | `height`, `temperature`, `moisture`, `population` | Keep the mode's own fill; glyphs suppressed |

The third row matters: a ramp's entire information content **is** its fill.
Rendering it on bare paper renders nothing. Glyphs are suppressed there too —
they would fight the ramp.

## Glyph placement

`placeGlyphs(cells, projection, widthPx, opts)` is a **pure function** with no
canvas, DOM or SVG contact. That is structural, not incidental: glyph size and
collision are screen-space decisions, and Mercator vs Mollweide distort area
wildly at high latitude. Inside a per-cell draw callback this logic would be
written once per substrate and would drift.

- Candidates come from elevation band and biome; relief outranks vegetation, so
  a forested mountain reads as a mountain.
- Thinning is greedy against a uniform grid, in **pixels**, so density is
  constant across output sizes.
- Size and spacing both scale by `widthPx / 1024`, so an 8192px export is the
  same map at higher resolution, not a denser one. Always pass the **output**
  width.

Shapes are procedural SVG path strings (`glyphPaths.ts`). One string serves both
substrates: SVG embeds it, Canvas2D feeds it to `Path2D`, which accepts SVG path
syntax. That is what guarantees both surfaces draw the identical shape.

## The mirror trap

**Every** 2D render path here flips horizontally — Map2D's main render, Map2D's
Dymaxion source buffer, both PNG export paths, and `exportSVG`. Polygons and
strokes come from the same mirrored path generator and land correctly. **Glyphs
do not**: they are drawn from their own coordinates and come out backwards.

Substrates take a `mirrored` flag and handle it in `drawGlyph` alone, by applying
the same flip again — two mirrors compose to the identity — and flipping the
glyph's x so it lands back over its own cell. Passes stay free of surface
geometry.

Do **not** hoist glyphs into a separate unmirrored group. They would be flipped
twice.

## The globe

`bakeStyleTexture` renders the real 2D style to an equirectangular canvas, and
the cell mesh samples it through UVs. The mesh keeps its displacement, so relief
still reads; the look is identical to the 2D map because it IS the 2D map.

Three things are load-bearing:

- **The bake is NOT mirrored.** The 2D screen and export paths flip horizontally
  for their own reasons, but UVs come straight from longitude, so a flip here
  puts the world back to front.
- **The antimeridian seam.** A cell straddling lon ±180 has vertices at u ≈ 0.99
  and u ≈ 0.01; interpolated as-is, that triangle samples the entire texture
  backwards and draws a bright smear down the seam. `buildGlobeUVs` pushes the
  low values past 1.0 instead, which is correct only because the texture uses
  `RepeatWrapping`.
- **The material is UNLIT and its JSX `key` matters.** Unlit because the baked
  texture already contains the hillshade pass — a lit material shades it twice
  and warm paper renders as dark grey-brown. The `key` because two branches are
  both `meshStandardMaterial`: without distinct keys React patches props onto one
  material instance, and three.js compiles its shader from the material's feature
  set, so adding `map` to a material that never had one does nothing until the
  program is rebuilt. The globe silently kept its vertex colours.

Default bake is 2048×1024 (~8MB) — chosen for a 16GB M1 with tight thermals.

## Adding a style

1. Implement `MapStyle` in `utils/mapStyle/style<Name>.ts`: an id, a name, a
   `StylePalette`, a `fillPolicy`, and a `passes` array.
2. Add the id to `MapStyleId` in `types.ts` and register it in `index.ts`.
3. Reuse the passes in `passes.ts` where they fit; they are parameterised by
   palette.

Nothing else needs to change. `Map2D`, both PNG paths and the SVG export all
drive the registry.

## Where state lives

`mapStyleId` sits beside `viewMode` in `hooks/useWorldEngine.ts` and is
prop-drilled through `ShellApp`, per the project's no-Context rule.

It is **not** a `WorldParam`, deliberately: it never influences generation, so
`tests/paramLiveness.test.ts` would fail it, and it must be changeable without
regenerating the world.

## Testing note

The unit tests cover every piece in isolation, and all of them passed while the
map was visibly broken in two ways — coastlines collapsed to a point, and a
full-opacity seam stroke turned the map into a honeycomb. Both were found by
rendering a real SVG and looking at it.

When changing a pass or a substrate, render one and look at it. The recipe:
export an SVG to a file, open it in headless Chromium, screenshot it.
