import { describe, it, expect, beforeAll } from 'vitest';
import { generateWorld, isLakeCell } from '../utils/worldGen';
import { placeSettlements } from '../utils/settlements';
import { makeParams } from './helpers';
import { WorldData } from '../types';

// Task 3: per-faction quota-then-score settlement placement with tiers.
// placeSettlements is not yet wired into the pipeline (a later task wires it
// into recalculateCivs), so these tests generate a world once, then call it
// manually against that world's civData/regionId/capitalId state.
describe('placeSettlements', () => {
  let world: WorldData;

  beforeAll(async () => {
    world = await generateWorld(makeParams({ points: 3000 }));
  }, 60000);

  it('gives every faction with land at least one settlement, and total scales with settlementDensity', () => {
    placeSettlements(world, world.params);
    const settlements = world.civData!.settlements;
    expect(settlements.length).toBeGreaterThan(0);

    for (const faction of world.civData!.factions) {
      const land = world.cells.filter(c => c.regionId === faction.id && c.height >= world.params.seaLevel && !isLakeCell(c));
      if (land.length === 0) continue;
      const factionSettlements = settlements.filter(s => s.factionId === faction.id);
      expect(factionSettlements.length).toBeGreaterThanOrEqual(1);
    }

    // Total scales with settlementDensity: a higher density should not produce
    // fewer settlements than the default.
    const lowWorld = structuredClone(world);
    const lowParams = { ...world.params, settlementDensity: 0.3 };
    placeSettlements(lowWorld, lowParams);
    const lowCount = lowWorld.civData!.settlements.length;

    const highWorld = structuredClone(world);
    const highParams = { ...world.params, settlementDensity: 3.0 };
    placeSettlements(highWorld, highParams);
    const highCount = highWorld.civData!.settlements.length;

    expect(highCount).toBeGreaterThan(lowCount);
  }, 60000);

  it('assigns exactly one capital-tier settlement per faction, at the faction capitalId', () => {
    const w = structuredClone(world);
    placeSettlements(w, w.params);
    const settlements = w.civData!.settlements;

    for (const faction of w.civData!.factions) {
      const land = w.cells.filter(c => c.regionId === faction.id && c.height >= w.params.seaLevel && !isLakeCell(c));
      if (land.length === 0) continue;
      const capitals = settlements.filter(s => s.factionId === faction.id && s.tier === 'capital');
      expect(capitals).toHaveLength(1);
      expect(capitals[0].cellId).toBe(faction.capitalId);
      expect(capitals[0].isCapital).toBe(true);
    }
  }, 60000);

  it('names every settlement — never empty or null', () => {
    const w = structuredClone(world);
    placeSettlements(w, w.params);
    const settlements = w.civData!.settlements;
    expect(settlements.length).toBeGreaterThan(0);
    for (const s of settlements) {
      expect(typeof s.name).toBe('string');
      expect(s.name.length).toBeGreaterThan(0);
    }
  }, 60000);

  it('is deterministic for a fixed civSeed and varies when civSeed changes', () => {
    const w1 = structuredClone(world);
    placeSettlements(w1, w1.params);
    const w2 = structuredClone(world);
    placeSettlements(w2, w2.params);

    const sig = (w: WorldData) => w.civData!.settlements
      .map(s => `${s.cellId}|${s.tier}|${s.name}|${s.factionId}|${s.isPort}`)
      .join(';');
    expect(sig(w1)).toBe(sig(w2));

    const w3 = structuredClone(world);
    const differentParams = { ...w3.params, civSeed: w3.params.civSeed + '_alt' };
    placeSettlements(w3, differentParams);
    expect(sig(w3)).not.toBe(sig(w1));
  }, 60000);

  it('spreads settlements across biomes — the single best biome does not dominate', () => {
    const w = structuredClone(world);
    placeSettlements(w, w.params);
    const settlements = w.civData!.settlements;

    let checked = 0;
    for (const faction of w.civData!.factions) {
      const factionSettlements = settlements.filter(s => s.factionId === faction.id);
      if (factionSettlements.length < 5) continue; // too small a sample to judge distribution

      const biomeCounts = new Map<string, number>();
      for (const s of factionSettlements) {
        const biome = w.cells[s.cellId].biome;
        biomeCounts.set(biome, (biomeCounts.get(biome) || 0) + 1);
      }
      if (biomeCounts.size < 3) continue; // guard only applies when the faction spans >=3 biomes

      checked++;
      const topShare = Math.max(...biomeCounts.values()) / factionSettlements.length;
      expect(topShare).toBeLessThan(0.6);
    }
    // Prove the guard above actually ran against a real faction — without this,
    // an unlucky world where no faction qualifies would pass having asserted
    // nothing about distribution at all.
    expect(checked).toBeGreaterThan(0);
  }, 60000);
});
