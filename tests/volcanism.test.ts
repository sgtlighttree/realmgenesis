import { describe, expect, it } from 'vitest';

import { BiomeType, Cell } from '../types';
import { generateWorld } from '../utils/worldGen';
import { makeParams } from './helpers';

// VOLCANIC was extinct before S28: `determineBiome` gated it on
// `landH > 0.85 && temp > -5`, and D8b's 6.5 C/km lapse rate made those two
// conditions mutually exclusive — landH 0.85 is ~6.5 km, i.e. ~42 C of cooling.
// Nothing caught it because every biome test used synthetic inputs, so the
// classifier passed while no generated world contained the biome.
//
// The lesson generalises, and it is why these are presence tests over a REAL
// world rather than more unit tests: a biome can be unreachable while the
// function that returns it is perfectly correct.
const PTS = 12000;
const share = (cells: Cell[], seaLevel: number): number => {
  const land = cells.filter(c => c.height >= seaLevel);
  return land.filter(c => c.biome === BiomeType.VOLCANIC).length / land.length;
};

describe('volcanism', () => {
  it('puts volcanic ground on a default world', async () => {
    const w = await generateWorld(makeParams({ seed: 'volc-a', points: PTS }));
    expect(share(w.cells, w.params.seaLevel)).toBeGreaterThan(0);
  }, 120000);

  it('produces none at all when the param is zero', async () => {
    const w = await generateWorld(makeParams({ seed: 'volc-a', points: PTS, volcanism: 0 }));
    expect(share(w.cells, w.params.seaLevel)).toBe(0);
  }, 120000);

  it('gives more volcanic ground as the param rises', async () => {
    const low = await generateWorld(makeParams({ seed: 'volc-a', points: PTS, volcanism: 0.5 }));
    const high = await generateWorld(makeParams({ seed: 'volc-a', points: PTS, volcanism: 2 }));
    expect(share(high.cells, high.params.seaLevel))
      .toBeGreaterThan(share(low.cells, low.params.seaLevel));
  }, 240000);

  it('stays a minority land cover at the default setting', async () => {
    // Calibrated to a median of ~1% of land across 5 seeds at 20k, ranging
    // 0.5%-3%. The bound is deliberately loose: the threshold is ABSOLUTE, so a
    // tectonically active world SHOULD carry more volcanic ground than a quiet
    // one, and pinning a narrow band would be pinning the seed rather than the
    // behaviour. This guards the failure that matters — volcanic rock taking
    // over the map.
    const w = await generateWorld(makeParams({ seed: 'volc-b', points: PTS }));
    expect(share(w.cells, w.params.seaLevel)).toBeLessThan(0.15);
  }, 120000);

  it('only ever claims land cells that carry high uplift', async () => {
    const w = await generateWorld(makeParams({ seed: 'volc-a', points: PTS }));
    const sea = w.params.seaLevel;
    const volcanic = w.cells.filter(c => c.biome === BiomeType.VOLCANIC);
    expect(volcanic.length).toBeGreaterThan(0);
    for (const c of volcanic) {
      expect(c.height).toBeGreaterThanOrEqual(sea);
      // The threshold at volcanism = 1. Asserting the CONTRACT rather than the
      // absence of lakes: "no volcanic cell is also a lake" is vacuous, since a
      // cell holds exactly one biome — the lake pass simply overwrites it, and
      // a test cannot distinguish that from the lake never overlapping.
      expect(c.upliftAccum ?? 0).toBeGreaterThanOrEqual(2.5);
    }
  }, 120000);
});
