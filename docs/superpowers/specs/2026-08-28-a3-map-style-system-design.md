# A3 — Map style system (parchment)

_Design spec. Session 27, 2026-08-28._

## 1. Goal

Make RealmGenesis output **pinboard-worthy rather than diagnostic**: a parchment
style that renders the world as a drawn map — bare paper land, hand-drawn relief
glyphs, coastline strokes, hatched ocean, paper grain.

Scope decisions taken with Matt this session:

| Question | Answer |
|---|---|
| How many styles? | **One preset — parchment — done properly.** Other preset ideas are logged in §9, not built. |
| Which surfaces? | **Map2D and PNG/SVG export.** The 3D globe stays diagnostic and is untouched. |
| Relief representation? | **Glyphs over a soft hillshade.** |
| Land fill? | **Flat parchment. Glyphs carry the terrain information.** |

The globe is excluded deliberately. It paints per-cell vertex colours, so
hatching, stroke weights and paper grain have no natural home there — they would
need a baked equirectangular texture or custom shaders, roughly doubling the work
and risking M1 performance for a view that is used for inspection, not output.

## 2. The problem this architecture solves

Parchment must render on two different substrates:

- **Canvas2D** — `Map2D.tsx` (screen + Dymaxion source buffer), `export.ts` (PNG)
- **SVG string** — `exportVector.ts`

The failure mode is writing parchment twice and watching the two drift. Every
design decision below serves avoiding that.

## 3. Architecture

### 3.1 Substrate adapter

One interface, two implementations. Style passes are written **once** against it.

```ts
interface Substrate {
  fillPath(d: PathSpec, fill: string): void;
  strokePath(d: PathSpec, stroke: string, width: number): void;
  fillHatch(d: PathSpec, hatch: HatchSpec): void;   // pattern fill
  grain(spec: GrainSpec): void;                      // paper texture, full-bleed
  glyph(g: PlacedGlyph, ink: string): void;
  text(t: TextSpec): void;
}
```

- `Canvas2DSubstrate` wraps a `CanvasRenderingContext2D` + `d3.GeoPath`.
- `SvgSubstrate` accumulates element strings and emits `<defs>` for patterns and
  filters.

SVG loses nothing: hatching is `<pattern>`, and paper grain is
`<filter><feTurbulence>`. Both are native.

### 3.2 Style as an ordered pass list

A style is data plus an ordered list of passes:

```ts
interface MapStyle {
  id: 'default' | 'parchment';
  name: string;
  palette: StylePalette;
  passes: StylePass[];
}
type StylePass = (ctx: StyleRenderContext, sub: Substrate) => void;
```

Parchment's pass order:

1. `paperPass` — full-bleed base tone + grain
2. `oceanPass` — flat sea tone, hatch bands parallel to coast
3. `landPass` — bare parchment fill (or categorical fill, see §4)
4. `hillshadePass` — existing `computeShadeMap`, multiplied in at low opacity
5. `coastlinePass` — heavy ink stroke, plus a lighter offset "swash" line
6. `glyphPass` — draws the placed-glyph list from §3.3
7. `bordersPass` / `riversPass` / `routesPass` — restyled existing segment data
8. `labelsPass` — existing label data, parchment typography

Passes 7 and 8 reuse the segment/label computations that `exportVector.ts` and
`Map2D` already have. The style changes their ink, weight and dash — not their
geometry.

### 3.3 Glyph placement is a separate, substrate-independent stage

**This is the load-bearing split.** Glyph size and collision are inherently
screen-space decisions, and Mercator vs Mollweide distort area wildly at high
latitude. If placement lived inside a per-cell draw callback it would be written
twice — the exact failure this architecture exists to prevent.

Placement is a pure function:

```ts
placeGlyphs(
  cells: Cell[],
  projection: d3.GeoProjection,
  widthPx: number,
  opts: GlyphOptions,
): PlacedGlyph[]   // { x, y, kind, scale, seedRot }
```

Verified feasible: **all four render paths already build a `d3.GeoProjection`
and a `d3.geoPath`** (`Map2D.tsx:336`, `Map2D.tsx:362`, `export.ts`,
`exportVector.ts`). `projection([lon, lat])` yields pixel coordinates in every
one of them, so both substrates simply draw the returned list.

Placement rules:

- **Candidate selection** — a cell earns a glyph from its biome and elevation
  band: mountain above a height threshold, hill below it, forest/conifer from
  forest biomes, dune from `HOT_DESERT`, marsh from wetland-ish cells.
- **Thinning** — greedy pass over candidates sorted by prominence (elevation for
  relief, arbitrary-but-seeded for vegetation), rejecting any candidate within
  `minSpacingPx` of an accepted one. Spacing is in **pixels**, so density stays
  constant across output sizes and projections.
- **Scale** — glyph size derives from `widthPx`, not from cell size, so a 2048px
  and an 8192px export look like the same map at different resolutions.
- **Variation** — `seedRot` and a shape index come from a hash of `cell.id` plus
  `world.params.seed`, so glyphs vary but are deterministic and stable across
  re-renders.

Glyphs are **procedural vector paths**, not sprites or a font. That keeps SVG
export genuinely vector, avoids shipping assets, and lets each glyph vary by
seed.

### 3.4 Where the code lives

New directory `utils/mapStyle/`:

```
utils/mapStyle/
  types.ts        Substrate, MapStyle, StylePass, PlacedGlyph, palettes
  substrateCanvas.ts
  substrateSvg.ts
  placeGlyphs.ts  the pure placement stage
  glyphPaths.ts   procedural bezier shapes per glyph kind
  styleDefault.ts the existing look, ported as a style
  styleParchment.ts
  index.ts        registry + lookup
```

`Map2D.tsx` is already 1,155 lines and `export.ts` 672. None of the new code goes
in them; they gain a call into the pass runner.

## 4. Parchment × ViewMode — the explicit rule

"Flat parchment, glyphs carry the info" collides with the twelve existing view
modes, and the collision must be settled in the spec rather than discovered when
a mode renders blank.

`ViewMode` splits three ways:

| Group | Modes | Parchment treatment |
|---|---|---|
| **Glyph-carried** | `satellite`, `biome`, `height_bw` | **Bare paper.** Terrain is read from glyphs, coastline and hillshade. |
| **Categorical fill** | `political`, `province`, `culture`, `religion`, `plates` | **Bare paper + the categorical fill on top**, muted toward the paper and given a drawn edge. This is the money shot. |
| **Continuous ramp** | `height`, `temperature`, `moisture`, `population` | **Keep the existing ramp fill.** Parchment supplies paper, coastline, hillshade, glyph suppression and typography around it. |

Rationale for the third row: the entire information content of a ramp mode *is*
the fill. Rendering it on bare paper renders nothing. Glyphs are suppressed there
too, because they would fight the ramp.

Style is **orthogonal to `ViewMode`** — a separate axis, per ROADMAP A3's wording
("a style layer over the existing view modes"). It is not a thirteenth view mode.

## 5. The `ColorContext` refactor — its own commit, first

`getCellColor` currently takes **7 positional arguments** and is called from
**8 sites** (`Map2D` ×2, `WorldViewer`, `MiniMap`, `DymaxionPreview2D`,
`DymaxionNetPreview`, `export.ts` ×2, `exportVector.ts`, `exportGLB.ts`). Adding
style makes it 8+.

Collapse the tail into one object:

```ts
getCellColor(cell: Cell, mode: ViewMode, ctx: ColorContext): THREE.Color
// ColorContext = { seaLevel, factionColors?, cultureColors?, religionColors?, seasonalDelta?, style? }
```

This is **not** mechanical, despite looking it:

- Two CLAUDE.md invariants are keyed to the *positional* form and must be
  rewritten in the same commit, or the next agent reads a doc describing a
  signature that no longer exists:
  - "`seaLevel` must be passed to `getCellColor` as the third argument"
  - "`factionColors` map required for political rendering"
- `tests/` must be checked for callers whose expectations encode argument order.

It ships as a **separate commit before any style code**, with zero behaviour
change, verified by the full suite. Mixing it with the feature would give a
rendering regression two candidate causes.

## 5b. Where style state lives — NOT in `WorldParams`

Style is a **render** choice, not a generation input. It must not go into
`WorldParams`:

- `WorldParams` changes are what `tests/paramLiveness.test.ts` polices — it fails
  if a key stops influencing *generated output*. A style key would never
  influence generation, so it would either fail that test or need an exemption.
- Style must be changeable without regenerating the world, like `viewMode`.
- `export.ts`'s `withParamDefaults` and the import validator would both need
  entries for a value that has nothing to do with the world.

So `mapStyleId` lives beside `viewMode` in the shell's render state
(`hooks/useWorldEngine.ts`, prop-drilled by `ShellApp.tsx`, per the project's
no-Context rule), and is passed into the render paths the same way `viewMode` is.

**UI:** a Style selector in `Controls.tsx`, adjacent to the view-mode control,
listing the registry's styles by `name`. The export panel takes the same value so
a PNG or SVG exports in the style shown on screen.

## 6. Data flow

```
world.cells ─┬─> getCellColor(cell, mode, ColorContext) ──> fill colour
             ├─> computeShadeMap(cells, seaLevel) ────────> hillshade factor
             ├─> coastline / border / river segments ─────> stroked paths
             └─> placeGlyphs(cells, projection, widthPx) ─> PlacedGlyph[]
                                    │
                        MapStyle.passes (ordered)
                                    │
                     ┌──────────────┴──────────────┐
              Canvas2DSubstrate              SvgSubstrate
             (Map2D, PNG export)             (SVG export)
```

## 7. Known gotchas

1. **Map2D mirrors its Dymaxion source buffer** (`translate(w,0); scale(-1,1)` at
   `Map2D.tsx:367`) — CLAUDE.md requires the pick buffer to mirror the visible
   rasterisation. Glyphs and text drawn inside that transform come out
   **flipped**. Glyph and label passes must counter-transform, or run after the
   mirror is restored. Verify visually on the Dymaxion projection specifically.
2. **`world.geoJson` features are the fill geometry**; glyph placement uses
   `cell.center` through the projection instead. They must agree — a glyph must
   land inside its own cell's polygon.
3. **Antimeridian** — Map2D's route pass already handles the jump. Glyphs near
   lon ±180 must not be drawn twice or smeared; clip by pixel bounds after
   projection.
4. **D10 changed `shading.ts`** — water now hillshades. The parchment ocean pass
   must decide whether it wants that (probably yes, at low opacity, under the
   hatching).

## 8. Testing

The existing suite verifies the pure engine; rendering is checked manually in the
browser. This feature follows that split:

- **Unit-testable, and will be tested:**
  - `placeGlyphs` — determinism for a fixed seed; density scales with `widthPx`;
    no two accepted glyphs closer than `minSpacingPx`; every glyph lands inside
    its source cell's projected polygon.
  - `getCellColor` with `ColorContext` — the existing colour tests, ported.
  - Style registry lookup and the ViewMode rule table from §4.
- **Manual, in-browser:** the actual look, on Mercator and Dymaxion, plus a PNG
  and an SVG export opened and inspected.
- **Gates unchanged:** typecheck 0, lint 0 errors / ≤30 warnings, full suite.

## 9. Preset backlog — logged, NOT built

Matt asked that other preset ideas be recorded so they are not lost. None of
these are in scope for A3.

- **Blueprint** — cyan paper, white line work, technical annotation, grid always on.
- **Ink & wash** — monochrome linework with a single wash colour for elevation.
- **Antique nautical** — rhumb lines, compass roses, sea monsters, heavy graticule.
- **Satellite-clean** — the current satellite look with better typography and no
  chrome; an export-grade version of what exists.
- **Risk/boardgame** — flat saturated territory fills, thick black borders,
  province numbers.
- **Topographic** — contour-led (A4 contours already exist), minimal fill,
  survey-map typography.

## 10. Open question for Matt

**Which branch does A3 start from?** D10 is complete but unmerged on
`d10-seafloor-relief`, and A3 touches `shading.ts` and both export paths — files
D10 just changed. Branching from `main` guarantees a merge conflict in those
files; branching from D10 inherits work Matt has not yet eyeballed in the
browser. Merging D10 first is the clean option.
