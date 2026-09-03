import { describe, it, expect } from 'vitest';
import { buildCellColorCache } from '../utils/mapColorCache';
import { getCellColor } from '../utils/colors';
import { seasonalTemperatureDelta } from '../utils/seasons';
import { generateWorld } from '../utils/worldGen';
import { makeParams } from './helpers';

describe('buildCellColorCache', () => {
  it('matches per-cell getCellColor with shade + seasonal delta', async () => {
    const world = await generateWorld(makeParams({ points: 2000, seed: 'color-cache' }));
    const colorCtx = {
      seaLevel: world.params.seaLevel,
      factionColors: undefined, cultureColors: undefined, religionColors: undefined,
    };
    const shadeMap = new Float32Array(world.cells.length).fill(1);
    shadeMap[5] = 0.7;
    const cache = buildCellColorCache(world, 'biome', colorCtx, shadeMap);
    for (let i = 0; i < world.cells.length; i++) {
      const c = getCellColor(world.cells[i], 'biome', {
        ...colorCtx, seasonalDelta: seasonalTemperatureDelta(world.cells[i], world.params),
      });
      c.multiplyScalar(shadeMap[i]);
      expect(cache[i]).toBe(`#${c.getHexString()}`);
    }
  });

  it('works with a null shadeMap', async () => {
    const world = await generateWorld(makeParams({ points: 2000, seed: 'color-cache-2' }));
    const cache = buildCellColorCache(world, 'biome', { seaLevel: world.params.seaLevel }, null);
    expect(cache).toHaveLength(world.cells.length);
    expect(cache[0]).toMatch(/^#[0-9a-f]{6}$/);
  });
});
