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
the cell mesh samples it. The mesh keeps its displacement, so relief still reads;
the look is identical to the 2D map because it IS the 2D map.

Three things are load-bearing:

- **The texture coordinate is computed in the FRAGMENT SHADER, not per vertex.**
  `createStyledGlobeMaterial` injects `u = atan2(z, x)/2π + 0.5`,
  `v = 0.5 − asin(y)/π` into `MeshBasicMaterial`. An equirectangular `u` IS a
  longitude, and longitude is not linear across a triangle on a sphere, so a
  per-vertex `uv` attribute is interpolated wrongly by construction. The old
  `buildGlobeUVs` needed two guards to survive that — an antimeridian seam wrap
  and a polar collapse — and the collapse traded polar streaks for a spiral
  rosette at the pole. Per-fragment coordinates need neither guard.
  `customProgramCacheKey` is required, or three.js may return a cached program
  compiled without the injection.
- **The direction comes from `sphereDir`, not from `position`.**
  `buildGlobeDirs` emits the UNDISPLACED unit sphere direction per vertex, laid
  out to match the position buffer. `position` carries per-cell height, so
  normalizing it would make two neighbours at different heights sample different
  content either side of their shared edge — a UV jog on every cell boundary.
- **The texture has NO MIPMAPS** (`createStyleTexture`). `atan` jumps by a full
  turn across the antimeridian, so the screen-space derivative explodes on that
  one line and mip selection falls to the coarsest level — a blurred band down
  the seam. Linear filtering uses no derivative. The globe sits near 1:1 against
  a 2048-wide texture, so there is little to minify. `RepeatWrapping` on `u` lets
  the seam texels blend; `ClampToEdge` on `v` is right because the poles ARE the
  edge.
- **The bake is NOT mirrored.** The 2D screen and export paths flip horizontally
  for their own reasons, but the globe's coordinate comes straight from
  longitude, so a flip here puts the world back to front.
- **A wrapped surface needs a phase-aligned hatch.** The globe joins the bake's
  left and right edges, and a hatch is drawn in output PIXELS, so unless its
  pattern repeats a whole number of times across the width the phase differs
  either side of the join and the diagonals jog down the antimeridian as a thin
  vertical line. `StyleRenderContext.wrapsHorizontally` tells `oceanHatchPass` to
  snap the spacing; the adjustment is a fraction of a percent. Screenshot at
  **lon 180** to check it — a view at lon 0 is 180 degrees from the seam and
  shows nothing.
- **The pole is a CONTENT problem, not a geometry one.** Equirectangular content
  genuinely converges at the pole; no coordinate scheme changes that. What made
  it read as a defect was the ocean hatch — a fixed-frequency pattern wound
  round the singularity into a rosette. `oceanHatchPass` fades the hatch out
  over a latitude band so there is no pattern left to swirl; the pole degrades to
  plain paper instead. Coastlines and fills pinch, which is geographically
  honest.
- **Line weights scale by `lineScale`.** The globe shows one hemisphere filling
  the viewport, so weights tuned for a whole-world flat map read as heavy blobs.
  The bake halves them. Hatch spacing scales with output width for the same
  reason glyph size does — a fixed pixel spacing gives an 8192 export hair-fine
  hatching and the globe corduroy.
- **The material is UNLIT.** The baked texture already contains the hillshade
  pass, so a lit material shades it twice and warm paper renders as dark
  grey-brown. It is also built imperatively and attached with `<primitive>`,
  because `onBeforeCompile` is not a prop React can safely patch: where two
  branches render the same element type, React reconciles them into one material
  instance, and three.js compiles its shader from a material's feature set — so
  adding `map` to a material that never had one does nothing until the program is
  rebuilt. That is why the remaining declared branches still carry distinct
  `key`s.

**Looking at it.** `scripts/renderGlobePreview.mjs` screenshots the real styled
globe from fixed camera latitudes (90° is straight down on the north pole).
`renderMapPreview.mjs` renders the flat SVG and cannot see globe defects at all —
the polar rosette survived three fix attempts partly because nothing was looking
at the globe.

Default bake is 2048×1024 (~8MB) — chosen for a 16GB M1 with tight thermals.

## Adding a style

1. Implement `MapStyle` in `utils/mapStyle/style<Name>.ts`: an id, a name, a
   `StylePalette`, a `fillPolicy`, and a `passes` array.
2. Add the id to `MapStyleId` in `types.ts` and register it in `index.ts`.
3. Reuse the passes in `passes.ts` where they fit; they are parameterised by
   palette.

Nothing else needs to change. `Map2D`, both PNG paths and the SVG export all
drive the registry.

### The desk

`palette.desk` is not painted by a pass. `paintDesk` (`desk.ts`) runs in Map2D's
display effect, in SCREEN space, for **every** style — `default` included — and
never in an export. Three layers: the flat `desk` colour, an optional tiled
`deskTexture` image, and an optional `deskEdge` vignette.

Two things to know before setting them:

- **The desk is seen in the letterbox margin down each side, not in the corners.**
  A viewport is wide, so the mid-edge sits at ~0.87 of the corner radius. A
  vignette that reaches full strength there paints the margin near-black, which
  is the raw black the desk exists to replace. The stops are late and shallow on
  purpose.
- **`deskTexture` must load through `useDeskTexture`.** An undecoded image draws
  nothing and reports nothing — the same silent failure as the webfonts. It must
  also tile: it is drawn as a `repeat` pattern in device pixels.

One thing the bake supplies that a flat render does not: `polarHatchFadeDeg` and
`wrapsHorizontally` on the render context. A new style that hatches the ocean
inherits both — a polar fade and a seam-aligned spacing — without asking for
them. That is right for a sphere texture, but it is not obvious from the pass.

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
