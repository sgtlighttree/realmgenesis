import { describe, it, expect } from 'vitest';
import * as d3 from 'd3';

import { buildMapGeometryCache } from '../utils/mapCache';
import { generateWorld } from '../utils/worldGen';
import { makeParams } from './helpers';

// This suite runs under plain Node (no jsdom/happy-dom — see vite.config.ts),
// so the global `Path2D` the DOM provides in a browser does not exist here.
// mapCache.ts calls `new Path2D()` while building cellPaths/coast/borders/
// lakes, so a minimal stub is required just to let buildMapGeometryCache run
// at all; it doesn't need to record anything since these tests only assert on
// cellVerts/cellOffsets (the parity that matters) and cellPaths.length.
const g = globalThis as { Path2D?: typeof Path2D };
if (typeof g.Path2D === 'undefined') {
  g.Path2D = class Path2D {
    moveTo() {}
    lineTo() {}
    closePath() {}
  } as unknown as typeof Path2D;
}

// A tiny stub 2D context recording the path points a d3 canvas generator emits,
// split by moveTo into subpaths — a real generated world always has some cells
// straddling the antimeridian (split into 2+ closed pieces by d3's clip
// stream) and, for cylindrical projections like Mercator, polar cells that
// legitimately emit nothing at all. Comparing only a flat point count (as a
// naive recorder would) can't tell "matches the real trace" apart from "both
// happened to emit the same number of points for the wrong reason" on those
// cells, so this records per-subpath structure and lengths too.
function recordingCtx() {
  const subpaths: [number, number][][] = [];
  let current: [number, number][] | null = null;
  return {
    beginPath() {},
    closePath() {},
    moveTo(x: number, y: number) { current = [[x, y]]; subpaths.push(current); },
    lineTo(x: number, y: number) { current?.push([x, y]); },
    arc() {},
    subpaths,
    get pts() { return subpaths.flat(); },
  } as unknown as CanvasRenderingContext2D & { subpaths: [number, number][][]; pts: [number, number][] };
}

describe('buildMapGeometryCache', () => {
  it('cached cell vertices match a fresh d3 projection for every cell (sub-pixel at float32 precision)', async () => {
    const world = await generateWorld(makeParams({ points: 2000, seed: 'geo-cache' }));
    const W = 1024, H = 512;
    const projection = d3.geoMercator().fitSize([W, H], { type: 'Sphere' as const });
    const cache = buildMapGeometryCache(world, projection, W, H, { simplifyTolerancePx: 0 });

    const ctx = recordingCtx();
    const gen = d3.geoPath(projection, ctx);

    let exercisedNonEmptyCell = false;
    let exercisedMultiSubpathCell = false;

    for (let i = 0; i < world.cells.length; i++) {
      ctx.subpaths.length = 0;
      gen(world.geoJson.features[i] as d3.GeoPermissibleObjects);
      const fresh = ctx.pts;
      const off0 = cache.cellOffsets[i], off1 = cache.cellOffsets[i + 1];
      const cachedCount = off1 - off0;

      // cellOffsets bookkeeping must agree with the fresh trace's point count
      // for EVERY cell, including polar cells that emit zero points — an
      // off-by-one in the offset math would only show up on a cell whose
      // neighbor(s) are non-vacuous, which single-cell checks can't catch.
      expect(cachedCount).toBe(fresh.length);
      if (fresh.length > 0) exercisedNonEmptyCell = true;
      if (ctx.subpaths.length > 1) exercisedMultiSubpathCell = true;

      // Float32Array can't hold d3's float64 coordinates exactly (ULP at
      // x~900 is ~6e-5) — compare against the value AS STORED at float32
      // precision, which is the honest contract for this buffer.
      for (let k = 0; k < cachedCount; k++) {
        expect(cache.cellVerts[(off0 + k) * 2]).toBe(Math.fround(fresh[k][0]));
        expect(cache.cellVerts[(off0 + k) * 2 + 1]).toBe(Math.fround(fresh[k][1]));
      }
    }

    // Guard against a vacuously-passing suite: this seed/point-count must
    // actually exercise both a normal cell and an antimeridian-split
    // (multi-subpath) cell, or the loop above proves nothing.
    expect(exercisedNonEmptyCell).toBe(true);
    expect(exercisedMultiSubpathCell).toBe(true);
  });

  it('produces one Path2D per cell', async () => {
    const world = await generateWorld(makeParams({ points: 2000, seed: 'geo-cache-2' }));
    const projection = d3.geoEquirectangular().fitSize([800, 400], { type: 'Sphere' as const });
    const cache = buildMapGeometryCache(world, projection, 800, 400);
    expect(cache.cellPaths).toHaveLength(world.cells.length);
  });
});
