import { describe, it, expect } from 'vitest';
import { settlementScore } from '../utils/settlements';
import { BiomeType, Cell, WorldData, WorldParams } from '../types';

const cell = (over: Partial<Cell> = {}): Cell => ({
  height: 0.3, biome: BiomeType.GRASSLAND, flux: 0, neighbors: [], center: { x: 1, y: 0, z: 0 }, ...over,
} as Cell);
const params = { seaLevel: 0.2 } as WorldParams;
const world = { cells: [] } as unknown as WorldData;

describe('settlementScore', () => {
  it('is zero for cells at or below sea level', () => {
    expect(settlementScore(cell({ height: 0.1 }), world, params)).toBe(0);
  });
  it('ranks Mediterranean above hot desert', () => {
    const med = settlementScore(cell({ biome: BiomeType.MEDITERRANEAN }), world, params);
    const des = settlementScore(cell({ biome: BiomeType.HOT_DESERT }), world, params);
    expect(med).toBeGreaterThan(des);
  });
  it('never returns negative for high peaks', () => {
    expect(settlementScore(cell({ height: 0.95, biome: BiomeType.VOLCANIC }), world, params)).toBeGreaterThanOrEqual(0);
  });
});
