import { describe, it, expect } from 'vitest';
import * as d3 from 'd3';

import { generateWorld } from '../utils/worldGen';
import { makeParams } from './helpers';
import { placeGlyphs } from '../utils/mapStyle/placeGlyphs';

const WIDTH = 1024;

const projectionFor = (width: number) =>
  d3.geoEquirectangular().fitSize([width, width / 2], { type: 'Sphere' } as d3.GeoPermissibleObjects);

describe('placeGlyphs', () => {
  it('is deterministic for the same world and width', async () => {
    const world = await generateWorld(makeParams({ points: 2000 }));
    const opts = { seaLevel: world.params.seaLevel, seed: world.params.seed };
    const a = placeGlyphs(world.cells, projectionFor(WIDTH), WIDTH, opts);
    const b = placeGlyphs(world.cells, projectionFor(WIDTH), WIDTH, opts);
    expect(a).toEqual(b);
    expect(a.length).toBeGreaterThan(0);
  }, 30000);

  it('never places two glyphs closer than minSpacingPx', async () => {
    const world = await generateWorld(makeParams({ points: 2000 }));
    const minSpacingPx = 24;
    const placed = placeGlyphs(world.cells, projectionFor(WIDTH), WIDTH, {
      seaLevel: world.params.seaLevel, seed: world.params.seed, minSpacingPx,
    });
    for (let i = 0; i < placed.length; i++) {
      for (let j = i + 1; j < placed.length; j++) {
        const d = Math.hypot(placed[i].x - placed[j].x, placed[i].y - placed[j].y);
        expect(d).toBeGreaterThanOrEqual(minSpacingPx);
      }
    }
  }, 30000);

  it('places no glyph on a water cell', async () => {
    const world = await generateWorld(makeParams({ points: 2000 }));
    const placed = placeGlyphs(world.cells, projectionFor(WIDTH), WIDTH, {
      seaLevel: world.params.seaLevel, seed: world.params.seed,
    });
    for (const g of placed) {
      expect(world.cells[g.cellId].height).toBeGreaterThanOrEqual(world.params.seaLevel);
    }
  }, 30000);

  it('scales glyph size with output width, keeping density roughly constant', async () => {
    const world = await generateWorld(makeParams({ points: 2000 }));
    const opts = { seaLevel: world.params.seaLevel, seed: world.params.seed };
    const small = placeGlyphs(world.cells, projectionFor(512), 512, opts);
    const large = placeGlyphs(world.cells, projectionFor(2048), 2048, opts);
    // Size is proportional to width.
    expect(large[0].scale).toBeCloseTo(small[0].scale * 4, 5);
    // Density (glyphs per unit area) stays within a factor of two, because
    // spacing scales with width too — a big export is the same map, not a
    // denser one.
    const ratio = large.length / Math.max(1, small.length);
    expect(ratio).toBeGreaterThan(0.5);
    expect(ratio).toBeLessThan(2.0);
  }, 30000);

  it('emits only known glyph kinds', async () => {
    const world = await generateWorld(makeParams({ points: 2000 }));
    const placed = placeGlyphs(world.cells, projectionFor(WIDTH), WIDTH, {
      seaLevel: world.params.seaLevel, seed: world.params.seed,
    });
    const kinds = new Set(placed.map(g => g.kind));
    for (const k of kinds) {
      expect(['mountain', 'hill', 'forest', 'conifer', 'dune', 'marsh', 'ice']).toContain(k);
    }
  }, 30000);
});
