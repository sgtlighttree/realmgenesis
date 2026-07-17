import { describe, it, expect, beforeAll } from 'vitest';
import { generateWorld } from '../utils/worldGen';
import { getCellsInRadius, applyTerrainStroke, applyFlattenStroke, applyPoliticalStroke, applyBiomeStroke, refreshBiomes } from '../utils/paintUtils';
import { BiomeType, WorldData } from '../types';
import { makeParams } from './helpers';

let world: WorldData;

beforeAll(async () => {
  world = await generateWorld(makeParams());
}, 30000);

describe('getCellsInRadius', () => {
  it('returns the center at ring 0 and unique cells only', () => {
    const center = world.cells[0];
    const brush = getCellsInRadius(center, 2, world.cells);
    expect(brush[0].cell.id).toBe(center.id);
    expect(brush[0].ring).toBe(0);
    const ids = brush.map(b => b.cell.id);
    expect(new Set(ids).size).toBe(ids.length);
    brush.forEach(b => expect(b.ring).toBeLessThanOrEqual(2));
  });

  it('grows with brush size', () => {
    const center = world.cells[0];
    const small = getCellsInRadius(center, 0, world.cells);
    const large = getCellsInRadius(center, 3, world.cells);
    expect(small.length).toBe(1);
    expect(large.length).toBeGreaterThan(small.length);
  });
});

describe('stroke functions', () => {
  it('terrain strokes clamp heights to [0, 1]', () => {
    const center = world.cells[10];
    const brush = getCellsInRadius(center, 2, world.cells);

    brush.forEach(({ cell }) => { cell.height = 0.99; });
    for (let i = 0; i < 50; i++) applyTerrainStroke(brush, 2, 'raise', 'freeform', world.cells, 1.0);
    brush.forEach(({ cell }) => {
      expect(cell.height).toBeLessThanOrEqual(1);
    });

    for (let i = 0; i < 200; i++) applyTerrainStroke(brush, 2, 'lower', 'freeform', world.cells, 1.0);
    brush.forEach(({ cell }) => {
      expect(cell.height).toBeGreaterThanOrEqual(0);
    });
  });

  it('flatten pulls heights toward the sampled target', () => {
    const center = world.cells[20];
    const brush = getCellsInRadius(center, 1, world.cells);
    brush.forEach(({ cell }) => { cell.height = 0.9; });
    const before = center.height;
    applyFlattenStroke(brush, 1, 0.2, 1.0);
    expect(center.height).toBeLessThan(before);
    expect(center.height).toBeGreaterThanOrEqual(0.2);
  });

  it('political stroke writes and clears region/province ownership', () => {
    const center = world.cells[30];
    const brush = getCellsInRadius(center, 1, world.cells);

    applyPoliticalStroke(brush, 2, 1);
    brush.forEach(({ cell }) => {
      expect(cell.regionId).toBe(2);
      expect(cell.provinceId).toBe(1);
    });

    applyPoliticalStroke(brush, undefined, undefined); // eraser
    brush.forEach(({ cell }) => {
      expect(cell.regionId).toBeUndefined();
      expect(cell.provinceId).toBeUndefined();
    });
  });

  it('biome stroke paints and refreshBiomes restores physical classification', () => {
    const center = world.cells[40];
    const brush = getCellsInRadius(center, 1, world.cells);

    applyBiomeStroke(brush, BiomeType.HOT_DESERT);
    brush.forEach(({ cell }) => expect(cell.biome).toBe(BiomeType.HOT_DESERT));

    brush.forEach(({ cell }) => { cell.height = 0.1; }); // deep below sea level
    refreshBiomes(brush.map(b => b.cell), world.params.seaLevel);
    brush.forEach(({ cell }) => {
      expect([BiomeType.OCEAN, BiomeType.DEEP_OCEAN]).toContain(cell.biome);
    });
  });
});
