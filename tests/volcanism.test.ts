import { describe, expect, it } from 'vitest';

import { BiomeType } from '../types';
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
const landCells = (cells: { height: number }[], seaLevel: number) =>
  cells.filter(c => c.height >= seaLevel);
const share = (cells: { height: number; biome: BiomeType }[], seaLevel: number) => {
  const land = landCells(cells, seaLevel);
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

  it('never claims ocean or lake cells', async () => {
    const w = await generateWorld(makeParams({ seed: 'volc-a', points: PTS }));
    const sea = w.params.seaLevel;
    // The pass runs before the lake pass precisely so standing water wins.
    for (const c of w.cells) {
      if (c.biome === BiomeType.VOLCANIC) expect(c.height).toBeGreaterThanOrEqual(sea);
    }
    const lakeIsLake = w.cells.every(
      c => !(c.biome === BiomeType.LAKE || c.biome === BiomeType.SALT_LAKE)
        || c.biome !== BiomeType.VOLCANIC,
    );
    expect(lakeIsLake).toBe(true);
  }, 120000);
});
