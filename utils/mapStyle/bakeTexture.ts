import * as d3 from 'd3';

import { WorldData, ViewMode } from '../../types';
import { buildFactionColorMap, buildCultureColorMap, buildReligionColorMap } from '../colors';
import { computeCoastlineSegments } from '../boundaries';
import { computeShadeMap } from '../shading';
import { runStyle } from './passes';
import { placeGlyphs } from './placeGlyphs';
import { Canvas2DSubstrate } from './substrateCanvas';
import { MapStyle } from './types';

/**
 * Default bake size. 2048x1024 is ~8MB as an RGBA texture — deliberate on a
 * 16GB M1 with tight thermals. 4096 is sharper under close zoom and four times
 * the memory; make it a choice, not a default.
 */
export const BAKE_WIDTH = 2048;

/**
 * Render a map style to an equirectangular canvas, for use as a globe texture.
 *
 * The 3D globe paints each cell as a flat-coloured triangle fan, so there is no
 * surface for hatching, paper grain or glyphs to live on — per-cell colour can
 * only ever produce a beige globe, not a drawn map. Baking the real 2D style and
 * wrapping it keeps the actual look, and the cell mesh keeps its displacement so
 * relief still reads.
 *
 * NOT mirrored: the 2D screen and export paths flip horizontally for their own
 * reasons, but the globe's texture coordinate is derived straight from
 * longitude, so a flip here would put the world back to front.
 *
 * Returns null where there is no DOM, or if the style draws nothing.
 */
export const bakeStyleTexture = (
  world: WorldData,
  viewMode: ViewMode,
  style: MapStyle,
  showHillshade: boolean,
  width: number = BAKE_WIDTH,
): HTMLCanvasElement | null => {
  if (style.passes.length === 0) return null;
  if (typeof document === 'undefined') return null;

  const height = Math.round(width / 2);
  const canvas = document.createElement('canvas');
  canvas.width = width;
  canvas.height = height;
  const ctx = canvas.getContext('2d');
  if (!ctx) return null;

  const projection = d3.geoEquirectangular()
    .fitSize([width, height], { type: 'Sphere' } as d3.GeoPermissibleObjects);
  const pathGenerator = d3.geoPath(projection, ctx);

  /**
   * The globe shows one HEMISPHERE filling the viewport, so roughly half the
   * texture is magnified to full screen. Line weights, hatch density and glyph
   * sizes tuned for a flat map — where the whole world is in view — come out
   * heavy and blobby here, swamping the coastline shape. Everything scales down
   * by this factor; glyphs also get more spacing so they read as marks rather
   * than as scattered dirt.
   */
  const lineScale = 0.5;

  const glyphs = style.fillPolicy(viewMode) === 'ramp'
    ? []
    : placeGlyphs(world.cells, projection, width, {
      seaLevel: world.params.seaLevel,
      seed: world.params.seed,
      sizePx: 16 * lineScale * 1.4,
      minSpacingPx: 22 * 1.6,
    });

  const sub = new Canvas2DSubstrate(
    ctx, pathGenerator as unknown as (o: unknown) => unknown, width, height, false,
  );

  runStyle(style, {
    world,
    viewMode,
    widthPx: width,
    heightPx: height,
    glyphs,
    lineScale,
    // Fade the ocean hatch out towards the poles. Values chosen by looking:
    // `scripts/renderGlobePreview.mjs --views=90` puts the camera straight down
    // on the north pole, and the tight part of the swirl sat inside ~78 deg.
    polarHatchFadeDeg: [66, 82],
    shadeMap: showHillshade ? computeShadeMap(world.cells, world.params.seaLevel) : null,
    coastlines: computeCoastlineSegments(world),
    colorCtx: {
      seaLevel: world.params.seaLevel,
      factionColors: buildFactionColorMap(world.civData),
      cultureColors: buildCultureColorMap(world.cultures),
      religionColors: buildReligionColorMap(world.religions),
    },
  }, sub);

  return canvas;
};

/**
 * Per-vertex UNIT SPHERE DIRECTIONS for the globe's triangle-fan cell mesh, laid
 * out to match the position buffer exactly: 3 vertices (centre, v1, v2) per
 * triangle.
 *
 * This is NOT a UV buffer, and the difference is the whole point. The globe
 * material derives its texture coordinate in the FRAGMENT shader from the
 * interpolated direction, so no `u` or `v` is ever interpolated across a
 * triangle. See `WorldViewer`'s styled material for why.
 *
 * Directions are taken BEFORE displacement (`displayRadius`) and before the
 * cell-overhang inflation. Sampling the displaced surface instead would make
 * the texture coordinate depend on cell height, so two neighbours at different
 * heights would sample different content either side of their shared edge — a
 * visible jog along every cell boundary.
 */
export const buildGlobeDirs = (world: WorldData): Float32Array => {
  let triCount = 0;
  for (const cell of world.cells) triCount += cell.vertices.length;
  const dirs = new Float32Array(triCount * 9);

  const put = (o: number, p: { x: number; y: number; z: number }) => {
    const len = Math.hypot(p.x, p.y, p.z) || 1;
    dirs[o] = p.x / len;
    dirs[o + 1] = p.y / len;
    dirs[o + 2] = p.z / len;
  };

  let o = 0;
  for (const cell of world.cells) {
    for (let i = 0; i < cell.vertices.length; i++) {
      put(o, cell.center);
      put(o + 3, cell.vertices[i]);
      put(o + 6, cell.vertices[(i + 1) % cell.vertices.length]);
      o += 9;
    }
  }
  return dirs;
};
