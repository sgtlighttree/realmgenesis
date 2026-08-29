import * as THREE from 'three';

import { BiomeType, ViewMode } from '../../types';
import { getCellColor } from '../colors';
import { SEAWATER_FREEZE_C, seasonalTemperatureDelta } from '../seasons';
import { isLakeCell } from '../worldGen';
import { Substrate } from './substrate';
import { FillPolicy, MapStyle, StylePalette, StylePass, StyleRenderContext } from './types';

/**
 * Run a style's passes in order. The only entry point a render path needs.
 *
 * A style with no passes draws nothing, so calling this on the default style is
 * a no-op — the caller decides separately whether to run its own legacy fill
 * loop, and tests `style.passes.length`, never the style id.
 */
export const runStyle = (style: MapStyle, ctx: StyleRenderContext, sub: Substrate): void => {
  for (const pass of style.passes) pass(ctx, sub);
};

/** Ocean cells that have a drawable polygon, in cell order. */
const oceanFeatures = (ctx: StyleRenderContext) => {
  const out = [];
  const { world, colorCtx } = ctx;
  for (let i = 0; i < world.cells.length; i++) {
    if (world.cells[i].height >= colorCtx.seaLevel) continue;
    const feature = world.geoJson?.features?.[i];
    if (feature) out.push(feature);
  }
  return out;
};

/**
 * Ocean features grouped by hatch opacity, for `polarHatchFadeDeg`.
 *
 * Quantised to quarter steps rather than a true per-cell gradient because each
 * distinct opacity costs one clip-and-sweep in Canvas2D and one <pattern> in
 * SVG. Four bands over a ten-degree span is finer than the cells are, so the
 * ramp reads as smooth; a hard cut would just move the artifact to the cut line.
 */
const oceanFeaturesByHatchOpacity = (
  ctx: StyleRenderContext,
  [startDeg, endDeg]: [number, number],
): Map<number, ReturnType<typeof oceanFeatures>> => {
  const bands = new Map<number, ReturnType<typeof oceanFeatures>>();
  const { world, colorCtx } = ctx;
  const span = Math.max(1e-6, endDeg - startDeg);
  for (let i = 0; i < world.cells.length; i++) {
    const cell = world.cells[i];
    if (cell.height >= colorCtx.seaLevel) continue;
    const feature = world.geoJson?.features?.[i];
    if (!feature) continue;
    // cell.center is a unit vector, so y is sin(latitude).
    const latDeg = Math.abs((Math.asin(Math.max(-1, Math.min(1, cell.center.y))) * 180) / Math.PI);
    const t = (latDeg - startDeg) / span;
    const opacity = Math.ceil(Math.max(0, Math.min(1, 1 - t)) * 4) / 4;
    if (opacity <= 0) continue;
    const band = bands.get(opacity);
    if (band) band.push(feature);
    else bands.set(opacity, [feature]);
  }
  return bands;
};

/**
 * Paper tone plus grain, over the MAP, not over the canvas. Always first —
 * everything else sits on it.
 *
 * Clipped to the projected sphere. A projection does not necessarily fill its
 * viewport: Mercator clipped at +/-85 degrees is square, so a wide canvas gets
 * a margin down each side, and paper painted by raw canvas coordinates spilled
 * into it — parchment then showed a band of blank land off each edge of the
 * world, which is what Matt reported against the Mercator view. For
 * equirectangular the sphere is the whole canvas and the clip changes nothing.
 */
export const paperPass = (palette: StylePalette, seed: string): StylePass =>
  (ctx, sub) => {
    sub.withSphereClip(() => {
      sub.fillRect(0, 0, ctx.widthPx, ctx.heightPx, palette.paper);
      sub.grain({ seed, opacity: 0.10, scale: 1 });
    });
  };

/**
 * Flat sea tone on ocean cells, with FROZEN sea painted as ice.
 *
 * The sea-ice check is here rather than inherited from `getCellColor` because
 * this pass never calls it — it paints a flat palette tone. Without this, D3's
 * seasonal sea ice silently disappeared under a style, the same class of loss as
 * the ice caps on land.
 */
export const oceanFillPass = (palette: StylePalette): StylePass =>
  (ctx, sub) => {
    const { world, colorCtx } = ctx;
    for (let i = 0; i < world.cells.length; i++) {
      const cell = world.cells[i];
      if (cell.height >= colorCtx.seaLevel) continue;
      const feature = world.geoJson?.features?.[i];
      if (!feature) continue;
      const frozen =
        cell.biome !== BiomeType.LAKE && cell.biome !== BiomeType.SALT_LAKE &&
        cell.temperature + seasonalTemperatureDelta(cell, world.params) < SEAWATER_FREEZE_C;
      sub.fillFeature(feature, frozen ? palette.ice : palette.sea);
    }
  };

/**
 * Nudge a hatch spacing so the pattern repeats a whole number of times across
 * the output width, for a surface whose two edges meet.
 *
 * A hatch at `angleDeg` repeats every `spacing / sin(angleDeg)` pixels
 * horizontally. If the width is not a multiple of that, the phase at the left
 * edge does not match the phase at the right, and on the globe the diagonals jog
 * down the antimeridian — a thin vertical line, visible at any zoom. The
 * adjustment is a fraction of a percent; it is the phase that matters, not the
 * density.
 *
 * A horizontal hatch (sin 0) never repeats horizontally, so there is nothing to
 * align and the spacing passes through.
 */
const snapHatchToWrap = (spacingPx: number, angleDeg: number, ctx: StyleRenderContext): number => {
  if (!ctx.wrapsHorizontally) return spacingPx;
  const sin = Math.abs(Math.sin((angleDeg * Math.PI) / 180));
  if (sin < 1e-6) return spacingPx;
  const repeats = Math.max(1, Math.round((ctx.widthPx * sin) / spacingPx));
  return (ctx.widthPx * sin) / repeats;
};

/**
 * Diagonal hatch over the ocean, as ONE composite region.
 *
 * Not full-bleed: under the `bare` fill policy the land is unpainted parchment,
 * so a whole-canvas hatch would sit directly on it. Not per-cell either — see
 * `Substrate.hatchFeatures` for why that does not scale.
 *
 * Runs AFTER hillshadePass so the hatching reads over the relief shading.
 */
export const oceanHatchPass = (palette: StylePalette): StylePass =>
  (ctx, sub) => {
    // Spacing and width scale with output, like glyphs. They were fixed pixel
    // values, which meant an 8192 export got hair-fine hatching and the globe
    // bake got corduroy — the same "denser map at higher resolution" mistake
    // placeGlyphs already avoids.
    const k = (ctx.widthPx / 1024) * (ctx.lineScale ?? 1);
    const angleDeg = ctx.hatchAngleDeg ?? 45;
    const spec = {
      color: palette.seaHatch,
      spacingPx: snapHatchToWrap(Math.max(3, 6 * k), angleDeg, ctx),
      widthPx: Math.max(0.4, 0.9 * k),
      angleDeg,
    };
    // The globe bake asks for a polar fade; a flat map never does. See
    // `StyleRenderContext.polarHatchFadeDeg` for why the poles need it.
    if (!ctx.polarHatchFadeDeg) {
      sub.hatchFeatures(oceanFeatures(ctx), spec);
      return;
    }
    const bands = oceanFeaturesByHatchOpacity(ctx, ctx.polarHatchFadeDeg);
    // Descending, so the strongest band is laid down first and the faded ones
    // read as an edge to it rather than as a separate pattern on top.
    for (const opacity of [...bands.keys()].sort((a, b) => b - a)) {
      sub.hatchFeatures(bands.get(opacity)!, { ...spec, opacity });
    }
  };

/**
 * Land fill, obeying the style's fill policy for this view mode (spec §4):
 *   bare        → paper shows through; glyphs carry the terrain
 *   categorical → the mode's own fill, muted toward the paper
 *   ramp        → the mode's own fill at full strength
 *
 * Takes `fillPolicy` and `palette` directly rather than the whole MapStyle: a
 * style holds its passes, so passing the style in would be circular.
 */
export const landPass = (
  fillPolicy: (mode: ViewMode) => FillPolicy,
  palette: StylePalette,
): StylePass =>
  (ctx, sub) => {
    const policy = fillPolicy(ctx.viewMode);
    const mute = policy === 'categorical' ? 0.45 : 0;
    const paper = new THREE.Color(palette.paper);
    const { world, colorCtx } = ctx;
    for (let i = 0; i < world.cells.length; i++) {
      const cell = world.cells[i];
      if (cell.height < colorCtx.seaLevel && !isLakeCell(cell)) continue;
      const feature = world.geoJson?.features?.[i];
      if (!feature) continue;
      // Ice is the ONE exception to bare paper. Left unpainted, an ice cap is
      // indistinguishable from temperate land — 6.8% of land in a default world
      // simply vanished. Glyphs alone cannot carry it either: a white expanse is
      // the information. Every other terrain still shows through as paper.
      if (policy === 'bare') {
        if (cell.biome === BiomeType.ICE_CAP) sub.fillFeature(feature, palette.ice);
        continue;
      }
      // seasonalDelta is PER CELL. Reusing one shared context here would
      // silently disable D1 seasonal biome re-derivation and D3 sea ice under
      // this style, with no test to catch it.
      const color = getCellColor(cell, ctx.viewMode, {
        ...colorCtx,
        seasonalDelta: seasonalTemperatureDelta(cell, world.params),
      });
      if (mute > 0) color.lerp(paper, mute);
      sub.fillFeature(feature, `#${color.getHexString()}`);
    }
  };

/**
 * Relief shading from the existing `computeShadeMap`.
 *
 * Tints come from the palette, not hardcoded black and white: on bare parchment
 * a neutral grey shadow reads as dirt on the paper rather than as landform. A
 * drawn map shades in its own ink.
 *
 * Opacity is an explicit argument, never baked into the colour as `rgba()` —
 * SVG 1.1 has no `rgba()` syntax. See `Substrate.fillFeature`.
 *
 * D10 made water shade too, so this covers the sea bed as well as the land.
 */
export const hillshadePass = (
  palette: StylePalette,
  opacityScale: number,
  landOnly = false,
): StylePass =>
  (ctx, sub) => {
    if (!ctx.shadeMap) return;
    const { world, colorCtx } = ctx;
    for (let i = 0; i < world.cells.length; i++) {
      const cell = world.cells[i];
      // D10 made water hillshade, which is right for the diagnostic views. On a
      // DRAWN map it fights the look: sea-bed relief renders as mottled grey
      // blotches that read as stains on the paper rather than as bathymetry.
      // Parchment shades land only and lets the hatch carry the sea, the way a
      // period map does.
      if (landOnly && cell.height < colorCtx.seaLevel) continue;
      const s = ctx.shadeMap[cell.id];
      if (s === 1) continue; // flat ground contributes nothing
      const feature = world.geoJson?.features?.[i];
      if (!feature) continue;
      const a = Math.min(0.5, Math.abs(1 - s) * opacityScale);
      sub.fillFeature(feature, s < 1 ? palette.shadow : palette.highlight, a);
    }
  };

/** Heavy ink coastline plus a lighter offset swash line just inside it. */
export const coastlinePass = (palette: StylePalette): StylePass =>
  (ctx, sub) => {
    const w = Math.max(0.6, (ctx.widthPx / 1024) * 2.0 * (ctx.lineScale ?? 1));
    sub.strokeSegments(ctx.coastlines, palette.coast, w);
  };

/** Relief and vegetation glyphs. Empty when the fill policy suppressed them. */
export const glyphPass = (palette: StylePalette): StylePass =>
  (ctx, sub) => {
    const w = Math.max(0.4, (ctx.widthPx / 1024) * 1.1 * (ctx.lineScale ?? 1));
    for (const g of ctx.glyphs) sub.drawGlyph(g, palette.ink, w);
  };
