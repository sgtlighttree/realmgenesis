import { describe, it, expect } from 'vitest';
import { determineBiome } from '../utils/worldGen';
import { BiomeType } from '../types';

// Mirrors the classification table in ARCHITECTURE.md. seaLevel = 0.55.
const SL = 0.55;

describe('determineBiome', () => {
  it('classifies water by depth', () => {
    expect(determineBiome(0.1, 20, 0.5, SL)).toBe(BiomeType.DEEP_OCEAN); // < SL * 0.6
    expect(determineBiome(0.5, 20, 0.5, SL)).toBe(BiomeType.OCEAN);
  });

  it('classifies coastal fringe as beach only when warm', () => {
    const coastal = SL + 0.01 * (1 - SL); // landH = 0.01 < 0.02
    expect(determineBiome(coastal, 20, 0.5, SL)).toBe(BiomeType.BEACH);
    expect(determineBiome(coastal, 10, 0.5, SL)).not.toBe(BiomeType.BEACH);
  });

  it('classifies extreme elevations as volcanic unless frozen', () => {
    const high = SL + 0.9 * (1 - SL); // landH = 0.9 > 0.85
    expect(determineBiome(high, 10, 0.5, SL)).toBe(BiomeType.VOLCANIC);
    expect(determineBiome(high, -20, 0.5, SL)).toBe(BiomeType.ICE_CAP); // temp < -5 skips volcanic
  });

  it('classifies polar bands by temperature', () => {
    const land = SL + 0.3 * (1 - SL);
    expect(determineBiome(land, -15, 0.5, SL)).toBe(BiomeType.ICE_CAP);
    expect(determineBiome(land, -5, 0.5, SL)).toBe(BiomeType.TUNDRA);
  });

  it('classifies arid bands (moisture < 0.15) by temperature', () => {
    const land = SL + 0.3 * (1 - SL);
    expect(determineBiome(land, 30, 0.1, SL)).toBe(BiomeType.HOT_DESERT);
    expect(determineBiome(land, 15, 0.1, SL)).toBe(BiomeType.STEPPE);
    expect(determineBiome(land, 5, 0.1, SL)).toBe(BiomeType.COLD_DESERT);
  });

  it('classifies semi-arid bands (moisture < 0.4) by temperature', () => {
    const land = SL + 0.3 * (1 - SL);
    expect(determineBiome(land, 30, 0.3, SL)).toBe(BiomeType.TROPICAL_SAVANNA);
    expect(determineBiome(land, 15, 0.3, SL)).toBe(BiomeType.MEDITERRANEAN);
    expect(determineBiome(land, 5, 0.3, SL)).toBe(BiomeType.STEPPE);
  });

  it('classifies humid bands (moisture >= 0.4) by temperature', () => {
    const land = SL + 0.3 * (1 - SL);
    expect(determineBiome(land, 30, 0.7, SL)).toBe(BiomeType.TROPICAL_RAINFOREST);
    expect(determineBiome(land, 20, 0.7, SL)).toBe(BiomeType.TEMPERATE_RAINFOREST);
    expect(determineBiome(land, 10, 0.7, SL)).toBe(BiomeType.TEMPERATE_FOREST);
    expect(determineBiome(land, 2, 0.7, SL)).toBe(BiomeType.BOREAL_FOREST);
  });
});
