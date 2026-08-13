import { describe, it, expect } from 'vitest';
import { generateWorld } from '../utils/worldGen';
import { BiomeType, WorldData } from '../types';
import { makeParams } from './helpers';

const isLakeBiome = (b: BiomeType) => b === BiomeType.LAKE || b === BiomeType.SALT_LAKE;

// Deterministic seeds discovered by scanning at 300 points under the V3 terrain
// model. Rescanned for D7 plateElongation (band/chain plate seeding) —
// 'k2' no longer produces a lake under the shifted terrain; 'basin' does
// and matches the same shape 'k2' used to (1-cell salt endorheic).
//   'basin'     -> exactly one salt lake (1 cell, hot + arid endorheic basin)
//   'lakeworld' -> exactly one fresh lake (2 cells, cool basin) — unchanged
// The default 'test_seed' produces no depressions, which is why the existing
// determinism/regression signatures stay byte-identical.
const SALT_SEED = 'basin';
const FRESH_SEED = 'lakeworld';

describe('lakes', () => {
  it('produces the expected lakes for the salt-lake seed', async () => {
    const world = await generateWorld(makeParams({ seed: SALT_SEED }));
    expect(world.lakes).toBeDefined();
    expect(world.lakes!.length).toBe(1);
    const lake = world.lakes![0];
    expect(lake.cellIds.length).toBe(1);
    expect(lake.isSalt).toBe(true);
    expect(lake.isEndorheic).toBe(true);
    // Every lake cell carries the matching hydrology biome.
    lake.cellIds.forEach(id => {
      expect(world.cells[id].biome).toBe(BiomeType.SALT_LAKE);
    });
  }, 30000);

  it('classifies a cool basin as a fresh lake', async () => {
    const world = await generateWorld(makeParams({ seed: FRESH_SEED }));
    expect(world.lakes!.length).toBeGreaterThan(0);
    const fresh = world.lakes!.find(l => !l.isSalt);
    expect(fresh).toBeDefined();
    fresh!.cellIds.forEach(id => {
      expect(world.cells[id].biome).toBe(BiomeType.LAKE);
    });
  }, 30000);

  it('is deterministic from the terrain seed', async () => {
    const sig = (w: WorldData) =>
      (w.lakes ?? [])
        .map(l => `${l.cellIds.join(',')}|${l.surfaceLevel.toFixed(6)}|${l.outflowCellId}|${l.isSalt}|${l.isEndorheic}`)
        .join(';');
    const a = await generateWorld(makeParams({ seed: SALT_SEED }));
    const b = await generateWorld(makeParams({ seed: SALT_SEED }));
    expect(sig(a)).toBe(sig(b));
    expect(sig(a).length).toBeGreaterThan(0);
  }, 30000);

  it('every lake cell sits above sea level but below its lake surface, and lakes are contiguous', async () => {
    const world = await generateWorld(makeParams({ seed: SALT_SEED }));
    const seaLevel = world.params.seaLevel;

    world.lakes!.forEach(lake => {
      expect(lake.cellIds.length).toBeGreaterThan(0);

      lake.cellIds.forEach(id => {
        const c = world.cells[id];
        expect(isLakeBiome(c.biome)).toBe(true);
        // A depression: land elevation, flooded above its own terrain.
        expect(c.height).toBeGreaterThanOrEqual(seaLevel);
        expect(c.height).toBeLessThan(lake.surfaceLevel);
      });

      // Contiguity: the cell set is connected via the neighbor graph.
      const set = new Set(lake.cellIds);
      const seen = new Set<number>([lake.cellIds[0]]);
      const stack = [lake.cellIds[0]];
      while (stack.length) {
        const id = stack.pop()!;
        for (const nId of world.cells[id].neighbors) {
          if (set.has(nId) && !seen.has(nId)) {
            seen.add(nId);
            stack.push(nId);
          }
        }
      }
      expect(seen.size).toBe(lake.cellIds.length);

      // Outflow, when present, is a cell just outside the lake.
      if (lake.outflowCellId !== null) {
        expect(set.has(lake.outflowCellId)).toBe(false);
      }
    });
  }, 30000);

  it('places no capital, town, or population on lake cells', async () => {
    const world = await generateWorld(makeParams({ seed: SALT_SEED }));
    world.cells.forEach(c => {
      if (isLakeBiome(c.biome)) {
        expect(c.isCapital).toBeFalsy();
        expect(c.isTown).toBeFalsy();
        expect(c.population ?? 0).toBe(0);
      }
    });
  }, 30000);

  it('never starts a river in a lake nor threads one through a lake interior', async () => {
    const world = await generateWorld(makeParams({ seed: SALT_SEED }));
    expect(world.rivers!.length).toBeGreaterThan(0); // guards against a vacuous test

    // Map river render points back to cell centers so we can flag lake cells.
    // River points are scaled off the (unit) cell center by a radius factor, so
    // normalizing recovers the exact center key.
    const lakeCenters = new Set<string>();
    world.lakes!.forEach(l => l.cellIds.forEach(id => {
      const c = world.cells[id].center;
      lakeCenters.add(`${c.x.toFixed(4)},${c.y.toFixed(4)},${c.z.toFixed(4)}`);
    }));
    const key = (p: { x: number; y: number; z: number }) => {
      const r = Math.hypot(p.x, p.y, p.z);
      return `${(p.x / r).toFixed(4)},${(p.y / r).toFixed(4)},${(p.z / r).toFixed(4)}`;
    };

    world.rivers!.forEach(path => {
      // A river may only touch a lake as its final vertex (the shore). Its
      // source and every interior vertex must be dry land — lakes are open
      // water, not river channel. (Lake cells are excluded from source
      // candidates and the trace is cut on lake entry.)
      path.forEach((p, i) => {
        if (i < path.length - 1) {
          expect(lakeCenters.has(key(p))).toBe(false);
        }
      });
    });
  }, 30000);
});
