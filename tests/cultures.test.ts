import { describe, it, expect } from 'vitest';
import { generateWorld, isLakeCell } from '../utils/worldGen';
import { RNG } from '../utils/rng';
import { createNameGenerator } from '../utils/namegen';
import { makeParams, cultureSignature } from './helpers';

// C1: the culture layer. Cultures are seeded and expanded on their own RNG
// side-stream (civSeed + '_cultures') before faction seeding runs, so these
// tests exercise: determinism, full land coverage, home-cell spacing, the
// naming handoff from culture -> faction/province/town, and numCultures
// liveness. See tests/paramLiveness.test.ts for the geometry-isolation
// assertion (numCultures must not perturb civ/terrain signatures).

describe('culture layer (C1)', () => {
  it('is deterministic — same civSeed yields identical cultureId assignment and culture roster', async () => {
    const params = makeParams({ civSeed: 'cultures_determinism' });
    const a = await generateWorld(params);
    const b = await generateWorld(params);

    expect(a.cultures).toBeDefined();
    expect(a.cultures!.length).toBeGreaterThan(0);
    expect(cultureSignature(a)).toBe(cultureSignature(b));
  }, 30000);

  it('assigns a cultureId to every land cell and none to water cells', async () => {
    const world = await generateWorld(makeParams());
    const seaLevel = world.params.seaLevel;

    for (const cell of world.cells) {
      const isWater = cell.height < seaLevel || isLakeCell(cell);
      if (isWater) {
        expect(cell.cultureId, `water cell ${cell.id} should have no cultureId`).toBeUndefined();
      } else {
        expect(cell.cultureId, `land cell ${cell.id} is missing a cultureId`).toBeDefined();
      }
    }
  }, 30000);

  it('spaces culture home cells apart — no two cultures share a home cell', async () => {
    const world = await generateWorld(makeParams({ numCultures: 6, points: 800 }));
    const homeIds = world.cultures!.map(c => c.homeCellId);
    expect(new Set(homeIds).size).toBe(homeIds.length);
  }, 30000);

  it('numCultures changes how many culture regions are seeded', async () => {
    const world4 = await generateWorld(makeParams({ numCultures: 4 }));
    const world8 = await generateWorld(makeParams({ numCultures: 8 }));
    expect(world4.cultures!.length).toBe(4);
    expect(world8.cultures!.length).toBe(8);
  }, 30000);

  // Naming integration: with a single culture, every land cell (and thus
  // every faction capital) belongs to culture 0, which always speaks
  // params.nameStyle. Faction names are then generated in capital-acceptance
  // order from the '_facnames_<style>' stream — replayable independently.
  it('faction names are drawn from their capital culture\'s naming style', async () => {
    const params = makeParams({ numCultures: 1, nameStyle: 'norse', civSeed: 'naming_integration_civs' });
    const world = await generateWorld(params);
    expect(world.civData?.factions.length).toBeGreaterThan(0);

    for (const cell of world.cells) {
      const seaLevel = world.params.seaLevel;
      if (cell.height >= seaLevel && !isLakeCell(cell)) {
        expect(cell.cultureId).toBe(0);
      }
    }
    expect(world.cultures![0].nameStyle).toBe('norse');

    const rng = new RNG(params.civSeed + '_facnames_norse');
    const gen = createNameGenerator('norse', () => rng.next());
    const expectedNames = world.civData!.factions.map(() => gen.faction());
    expect(world.civData!.factions.map(f => f.name)).toEqual(expectedNames);
  }, 30000);
});
