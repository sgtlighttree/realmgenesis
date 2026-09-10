import { describe, it, expect, beforeAll } from 'vitest';
import { generateWorld, isLakeCell } from '../utils/worldGen';
import { placeSettlements, settlementScore } from '../utils/settlements';
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

  // S36: the city tier was DEAD below ~8 settlements per faction. The capital is
  // forced into `accepted` at index 0, so the first real candidate sat at i=1 and
  // cleared `i < accepted.length * 0.12` only when a faction held >8.3 settlements.
  // Measured: zero cities at 5k points. Rank fractions are now taken over the
  // non-capital candidates, so the top-ranked non-capital is always a city.
  it('gives every faction with 2+ settlements at least one city-tier settlement', () => {
    const w = structuredClone(world);
    placeSettlements(w, w.params);
    const settlements = w.civData!.settlements;

    let checked = 0;
    for (const faction of w.civData!.factions) {
      const mine = settlements.filter(s => s.factionId === faction.id);
      if (mine.length < 2) continue;
      checked++;
      expect(mine.filter(s => s.tier === 'city').length).toBeGreaterThanOrEqual(1);
    }
    // Non-vacuous: this world must actually contain a multi-settlement faction.
    expect(checked).toBeGreaterThan(0);
  }, 60000);

  // S36 review finding 1b. The browser-save load path (Controls.handleLoadBrowser
  // -> onLoadWorld) passes stored params through without withParamDefaults, so
  // placeSettlements must defend at its own read site. Without the `?? 1.0`, quota
  // is NaN and NaN defeats both loop guards, accepting EVERY ranked land cell.
  it('survives a missing settlementDensity instead of settling every land cell', () => {
    const w = structuredClone(world);
    const p = { ...w.params };
    delete (p as unknown as Record<string, unknown>).settlementDensity;
    placeSettlements(w, p);
    const settlements = w.civData!.settlements;
    const land = w.cells.filter(c => c.height >= w.params.seaLevel && !isLakeCell(c));

    expect(settlements.length).toBeGreaterThan(0);
    expect(Number.isFinite(settlements.length)).toBe(true);
    // The failure mode is "every land cell": stay far below it, and match what
    // an explicit density of 1.0 produces.
    expect(settlements.length).toBeLessThan(land.length * 0.25);
    const explicit = structuredClone(world);
    placeSettlements(explicit, { ...w.params, settlementDensity: 1.0 });
    expect(settlements.length).toBe(explicit.civData!.settlements.length);
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

    // What this guards (spec §Risks 1): quota-derived SPACING is what stops the
    // greedy rank-walk striping the best-scoring biome band. So the test compares
    // placement against the thing spacing exists to beat — the same number of
    // cells taken straight off the score ranking — rather than an absolute cap.
    //
    // It used to assert `topShare < 0.6`. That only holds while settlements are
    // sparse: the area-based quota (1b998f6) made them dense, and a faction whose
    // land is 46% one biome then correctly puts ~62% of its settlements there.
    // That is tracking the land, not a biome dominating it.
    //
    // Distance is total variation between the settlement and land biome mixes,
    // weighted by settlement count across factions. Measured over 6 seeds:
    // placement 0.13–0.20, score-only 0.44–0.65; with spacing disabled placement
    // collapses to >= 0.95 of score-only. The bound is the midpoint.
    const mix = (biomes: string[]): Map<string, number> => {
      const m = new Map<string, number>();
      for (const b of biomes) m.set(b, (m.get(b) || 0) + 1 / biomes.length);
      return m;
    };
    const tvd = (a: Map<string, number>, b: Map<string, number>): number => {
      let s = 0;
      for (const k of new Set([...a.keys(), ...b.keys()])) s += Math.abs((a.get(k) || 0) - (b.get(k) || 0));
      return s / 2;
    };

    let checked = 0;
    let placedDist = 0;
    let scoreOnlyDist = 0;
    for (const faction of w.civData!.factions) {
      const factionSettlements = settlements.filter(s => s.factionId === faction.id);
      if (factionSettlements.length < 5) continue; // too small a sample to judge distribution
      const placed = factionSettlements.map(s => w.cells[s.cellId].biome);
      if (new Set(placed).size < 3) continue; // guard only applies when the faction spans >=3 biomes

      const land = w.cells.filter(c =>
        c.regionId === faction.id && c.height >= w.params.seaLevel && !isLakeCell(c));
      const scoreOnly = [...land]
        .sort((a, b) => settlementScore(b, w, w.params) - settlementScore(a, w, w.params))
        .slice(0, placed.length)
        .map(c => c.biome);
      const landMix = mix(land.map(c => c.biome));

      checked++;
      placedDist += tvd(mix(placed), landMix) * placed.length;
      scoreOnlyDist += tvd(mix(scoreOnly), landMix) * placed.length;
    }
    expect(placedDist).toBeLessThan(scoreOnlyDist / 2);
    // Prove the guard above actually ran against a real faction — without this,
    // an unlucky world where no faction qualifies would pass having asserted
    // nothing about distribution at all.
    expect(checked).toBeGreaterThan(0);
  }, 60000);
});
