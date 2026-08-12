import { describe, it, expect } from 'vitest';
import { generateWorld } from '../utils/worldGen';
import { serializeWorld, deserializeWorld } from '../utils/worldTransfer';
import { digestWorld, diffDigests } from './helpers/worldDigest';
import { makeParams } from './helpers';
import { BiomeType } from '../types';

// Independent instrument for the "every buffer is listed as transferable"
// check. Deliberately does NOT reuse the shallow Object.values walk that
// utils/worldTransfer.ts uses to build `transfer` — it descends into plain
// objects and arrays so it can catch a buffer nested inside a wrapper (e.g.
// a future payload field of shape `{ tileOffsets: Ragged }` instead of two
// flat top-level keys), which a shallow walk would silently skip on both
// sides at once.
const collectBuffers = (
  v: unknown,
  out: Set<ArrayBuffer> = new Set(),
  seen: Set<object> = new Set(),
): Set<ArrayBuffer> => {
  if (v === null || typeof v !== 'object') return out;
  if (ArrayBuffer.isView(v)) { out.add((v as ArrayBufferView).buffer as ArrayBuffer); return out; }
  if (seen.has(v)) return out;
  seen.add(v);
  if (Array.isArray(v)) {
    for (const item of v) collectBuffers(item, out, seen);
  } else {
    for (const val of Object.values(v as Record<string, unknown>)) collectBuffers(val, out, seen);
  }
  return out;
};

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
    const reachable = collectBuffers(payload);
    expect(new Set(transfer)).toEqual(reachable);
  }, 30000);

  it('the recursive collector catches a buffer nested one level inside a plain object', () => {
    // Proves the previous test is not tautological. A payload field shaped
    // like `{ wrapper: { inner: Float64Array } }` hides its buffer from a
    // shallow Object.values walk (the same walk both the production
    // `transfer` builder and the old test used) but must not hide from the
    // recursive collector used above.
    const shaped = { wrapper: { inner: new Float64Array([1, 2, 3]) } };

    const shallow = new Set<ArrayBuffer>();
    for (const v of Object.values(shaped)) {
      if (ArrayBuffer.isView(v)) shallow.add((v as ArrayBufferView).buffer as ArrayBuffer);
    }
    expect(shallow.size).toBe(0); // the naive walk misses it entirely

    const found = collectBuffers(shaped);
    expect(found.has(shaped.wrapper.inner.buffer)).toBe(true); // the recursive one does not
  });

  it('throws on encode when a cell has a biome value unknown to the wire format', async () => {
    const w = await generateWorld(makeParams());
    w.cells[0].biome = 'Not A Real Biome' as BiomeType;
    expect(() => serializeWorld(w)).toThrow(/unknown biome/i);
  }, 30000);

  it('throws on decode when a biome byte is out of range for the wire format', async () => {
    const w = await generateWorld(makeParams());
    const { payload } = serializeWorld(w);
    payload.biome[0] = 255;
    expect(() => deserializeWorld(payload)).toThrow(/out-of-range biome/i);
  }, 30000);
});
