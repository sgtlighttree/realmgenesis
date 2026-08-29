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
 * reasons, but a texture is sampled through UVs computed straight from
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
 * Per-vertex UVs for the globe's triangle-fan cell mesh, laid out to match the
 * position buffer exactly: 3 vertices (centre, v1, v2) per triangle.
 *
 * **Antimeridian seam.** A cell straddling lon +/-180 has vertices at u ~ 0.99
 * and u ~ 0.01. Interpolated naively, that triangle samples the ENTIRE texture
 * backwards and draws a bright smear down the seam. Where a triangle's u range
 * exceeds half the texture, the low values are pushed past 1.0 instead; with
 * RepeatWrapping they sample correctly and the triangle stays narrow.
 */
export const buildGlobeUVs = (world: WorldData): Float32Array => {
  let triCount = 0;
  for (const cell of world.cells) triCount += cell.vertices.length;
  const uv = new Float32Array(triCount * 6);

  // u from longitude, v from latitude — the same mapping d3.geoEquirectangular
  // uses for the bake, so texel and cell agree.
  const uOf = (x: number, z: number) => (Math.atan2(z, x) / (2 * Math.PI)) + 0.5;
  const vOf = (y: number) => 0.5 - (Math.asin(Math.max(-1, Math.min(1, y))) / Math.PI);

  let o = 0;
  for (const cell of world.cells) {
    const c = cell.center;
    const cu = uOf(c.x, c.z);
    const cv = vOf(c.y);
    for (let i = 0; i < cell.vertices.length; i++) {
      const v1 = cell.vertices[i];
      const v2 = cell.vertices[(i + 1) % cell.vertices.length];
      const us = [cu, uOf(v1.x, v1.z), uOf(v2.x, v2.z)];
      const vs = [cv, vOf(v1.y), vOf(v2.y)];

      // Wrap each vertex to the equivalent u NEAREST the cell centre — the
      // shortest way round the cylinder.
      //
      // An earlier version pushed every u below 0.5 past 1.0 instead. That is
      // wrong whenever the three vertices do not straddle 0.5 together: a cell
      // at u = [0.1, 0.6, 0.9] came out [1.1, 0.6, 0.9], still spanning half the
      // texture. Measured on a 5k world, 76 triangles kept a span up to 0.747 —
      // each one smeared three quarters of the map across a single cell, which
      // showed as dark streaks scattered over open ocean.
      for (let k = 1; k < 3; k++) {
        if (us[k] - us[0] > 0.5) us[k] -= 1;
        else if (us[0] - us[k] > 0.5) us[k] += 1;
      }

      // POLES. Equirectangular has a genuine singularity there: a cell touching
      // a pole spans a wide longitude range legitimately, and no UV assignment
      // can avoid stretching it. Measured on a 5k world, 76 triangles — all at
      // |lat| >= 84 degrees — still spanned up to 0.63 of the texture after the
      // wrap above, each smearing most of the map across one cell as dark
      // streaks near the poles.
      //
      // Collapsing them to the centre's longitude samples a narrow strip
      // instead. That is not a fudge: near a pole every longitude maps to
      // near-identical content, so the strip is representative, and a thin
      // correct band beats a wide wrong one.
      if (Math.max(...us) - Math.min(...us) > 0.08) {
        us[1] = us[0];
        us[2] = us[0];
      }

      for (let k = 0; k < 3; k++) {
        uv[o] = us[k];
        uv[o + 1] = vs[k];
        o += 2;
      }
    }
  }
  return uv;
};
