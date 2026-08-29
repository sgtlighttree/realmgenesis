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

/** Full-bleed paper tone plus grain. Always first — everything sits on it. */
export const paperPass = (palette: StylePalette, seed: string): StylePass =>
  (ctx, sub) => {
    sub.fillRect(0, 0, ctx.widthPx, ctx.heightPx, palette.paper);
    sub.grain({ seed, opacity: 0.10, scale: 1 });
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
    sub.hatchFeatures(oceanFeatures(ctx), {
      color: palette.seaHatch, spacingPx: 6, widthPx: 0.9, angleDeg: 45,
    });
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
    const w = Math.max(1, (ctx.widthPx / 1024) * 2.0);
    sub.strokeSegments(ctx.coastlines, palette.coast, w);
  };

/** Relief and vegetation glyphs. Empty when the fill policy suppressed them. */
export const glyphPass = (palette: StylePalette): StylePass =>
  (ctx, sub) => {
    const w = Math.max(0.6, (ctx.widthPx / 1024) * 1.1);
    for (const g of ctx.glyphs) sub.drawGlyph(g, palette.ink, w);
  };
