import { describe, it, expect } from 'vitest';
import { generateWorld } from '../utils/worldGen';
import { makeParams, terrainSignature } from './helpers';

describe('D8b lapse-rate coupling', () => {
  it('hatch off vs on produces a different terrain signature (lapse is live)', async () => {
    const off = await generateWorld(makeParams({ physicalClimate: false }));
    const on = await generateWorld(makeParams({ physicalClimate: true }));
    expect(terrainSignature(on)).not.toBe(terrainSignature(off));
  }, 120000);

  it('a taller datum cools high peaks further when the hatch is on', async () => {
    // Larger maxElevationM => more metres per height unit => colder peaks =>
    // a different terrain signature. Proves lapse reads the datum, not raw height.
    const lo = await generateWorld(makeParams({ physicalClimate: true, maxElevationM: 4500 }));
    const hi = await generateWorld(makeParams({ physicalClimate: true, maxElevationM: 18000 }));
    expect(terrainSignature(hi)).not.toBe(terrainSignature(lo));
  }, 120000);
});
