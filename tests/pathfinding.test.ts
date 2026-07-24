import { describe, it, expect } from 'vitest';
import { isWaterCell, landTerrainStepCost, aStar } from '../utils/pathfinding';
import { BiomeType, type Cell } from '../types';

const cell = (over: Partial<Cell>): Cell =>
  ({ id: 0, center: { x: 0, y: 0, z: 1 }, height: 0.5, neighbors: [], biome: BiomeType.STEPPE, ...over }) as Cell;

describe('pathfinding primitives', () => {
  it('isWaterCell: below sea level or a lake biome is water', () => {
    expect(isWaterCell(cell({ height: 0.1 }), 0.3)).toBe(true);
    expect(isWaterCell(cell({ height: 0.5, biome: BiomeType.LAKE }), 0.3)).toBe(true);
    expect(isWaterCell(cell({ height: 0.5, biome: BiomeType.SALT_LAKE }), 0.3)).toBe(true);
    expect(isWaterCell(cell({ height: 0.5 }), 0.3)).toBe(false);
  });

  it('landTerrainStepCost: base 1 + slope*20, biome multipliers on the target', () => {
    const flat = cell({ height: 0.5 });
    expect(landTerrainStepCost(flat, cell({ height: 0.5 }))).toBeCloseTo(1);
    expect(landTerrainStepCost(flat, cell({ height: 0.55 }))).toBeCloseTo(1 + 0.05 * 20);
    expect(landTerrainStepCost(flat, cell({ height: 0.5, biome: BiomeType.HOT_DESERT }))).toBeCloseTo(2);
    expect(landTerrainStepCost(flat, cell({ height: 0.5, biome: BiomeType.ICE_CAP }))).toBeCloseTo(4);
    expect(landTerrainStepCost(flat, cell({ height: 0.5, biome: BiomeType.VOLCANIC }))).toBeCloseTo(5);
  });

  it('aStar: finds a path on a simple line graph and returns null past the expansion cap', () => {
    // 0 - 1 - 2 - 3 line
    const cells: Cell[] = [0, 1, 2, 3].map(i =>
      cell({ id: i, height: 0.5, neighbors: [i - 1, i + 1].filter(n => n >= 0 && n <= 3) }));
    const path = aStar(cells, 0, 3, () => 1, () => 0, 100);
    expect(path).toEqual([0, 1, 2, 3]);
    expect(aStar(cells, 0, 3, () => 1, () => 0, 1)).toBeNull();
  });
});
