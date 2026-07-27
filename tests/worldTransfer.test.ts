import { describe, it, expect } from 'vitest';
import { generateWorld } from '../utils/worldGen';
import { serializeWorld, deserializeWorld } from '../utils/worldTransfer';
import { digestWorld, diffDigests } from './helpers/worldDigest';
import { makeParams } from './helpers';
import { BiomeType } from '../types';

describe('worldTransfer round trip', () => {
  it('is digest-identical across serialize -> deserialize', async () => {
    const w = await generateWorld(makeParams());
    const { payload } = serializeWorld(w);
    expect(diffDigests(digestWorld(w), digestWorld(deserializeWorld(payload)))).toEqual([]);
  }, 30000);

  it('preserves undefined optionals as undefined, not 0', async () => {
    const w = await generateWorld(makeParams());
    const back = deserializeWorld(serializeWorld(w).payload);
    const water = w.cells.findIndex(c => c.regionId === undefined);
    expect(water).toBeGreaterThanOrEqual(0);
    expect(back.cells[water].regionId).toBeUndefined();
    expect(back.cells[water]).not.toHaveProperty('regionId', 0);
  }, 30000);

  it('preserves a cell with null geoJson geometry', async () => {
    const w = await generateWorld(makeParams());
    w.geoJson.features[0].geometry = null;
    const back = deserializeWorld(serializeWorld(w).payload);
    expect(back.geoJson.features[0].geometry).toBeNull();
  }, 30000);

  it('preserves an empty vertices ring', async () => {
    const w = await generateWorld(makeParams());
    w.cells[3].vertices = [];
    const back = deserializeWorld(serializeWorld(w).payload);
    expect(back.cells[3].vertices).toEqual([]);
  }, 30000);

  it('preserves every biome value through the Uint8 encoding', async () => {
    const w = await generateWorld(makeParams());
    const all = Object.values(BiomeType);
    w.cells.forEach((c, i) => { c.biome = all[i % all.length]; });
    const back = deserializeWorld(serializeWorld(w).payload);
    expect(back.cells.map(c => c.biome)).toEqual(w.cells.map(c => c.biome));
  }, 30000);

  it('lists every buffer it hands out as transferable', async () => {
    const w = await generateWorld(makeParams());
    const { payload, transfer } = serializeWorld(w);
    const inPayload = new Set<ArrayBuffer>();
    for (const v of Object.values(payload as unknown as Record<string, unknown>)) {
      if (ArrayBuffer.isView(v)) inPayload.add(v.buffer as ArrayBuffer);
    }
    expect(new Set(transfer)).toEqual(inPayload);
  }, 30000);
});
