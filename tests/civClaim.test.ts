import { describe, it, expect } from 'vitest';
import { generateWorld, isLakeCell } from '../utils/worldGen';
import { makeParams } from './helpers';

// Regression guards for the 2026-08-22 territorial-waters fix. Before it,
// recalculateCivs charged one ocean step at `waterCrossingCost * 50` (40 at
// defaults) against a budget of `territorialWaters * 50` (7.5), so no water cell
// was ever claimable from any capital at any distance — and because lakes were
// folded into the same `isWater` test, every lake was both unclaimable AND an
// impassable wall that expansion routed around.
describe('civ expansion claims water correctly', () => {
  it('claims lake cells — an enclosed lake belongs to the realm around it', async () => {
    const world = await generateWorld(makeParams({ points: 3000 }));
    const lakes = world.cells.filter(isLakeCell);

    // A guard on the guard: no lakes would make this pass trivially.
    expect(lakes.length).toBeGreaterThan(0);
    const unclaimed = lakes.filter((c) => c.regionId === undefined);
    expect(unclaimed).toHaveLength(0);
  }, 60000);

  it('claims territorial waters off the coast, but not the open ocean', async () => {
    const world = await generateWorld(makeParams({ points: 3000 }));
    const sea = world.params.seaLevel;
    const ocean = world.cells.filter((c) => c.height < sea && !isLakeCell(c));

    expect(ocean.length).toBeGreaterThan(0);
    const claimed = ocean.filter((c) => c.regionId !== undefined).length;
    const share = claimed / ocean.length;

    // Some coastal water is claimed...
    expect(claimed).toBeGreaterThan(0);
    // ...but territorial waters are a coastal fringe, not the whole sea. The
    // step allowance at default territorialWaters (0.15 -> 2 ocean steps) keeps
    // this well under half; the loose bound is deliberate, since the exact share
    // depends on how convoluted the coastline of a given seed is.
    expect(share).toBeLessThan(0.5);
  }, 60000);

  it('leaves no unclaimed land on a generated world', async () => {
    // The old `cost > 200` frontier cap claimed a cell and then refused to
    // expand from it, so expansion died partway across a landmass — worst on
    // mountains, where landTerrainStepCost burns the budget fastest.
    // Genuinely isolated islands beyond the water-step allowance stay unclaimed
    // and SHOULD; hence a small tolerance rather than zero.
    const world = await generateWorld(makeParams({ points: 3000 }));
    const sea = world.params.seaLevel;
    const land = world.cells.filter((c) => c.height >= sea && !isLakeCell(c));

    expect(land.length).toBeGreaterThan(0);
    const unclaimed = land.filter((c) => c.regionId === undefined).length;
    expect(unclaimed / land.length).toBeLessThan(0.02);
  }, 60000);

  it('territorialWaters = 0 claims no ocean at all', async () => {
    // The slider floor must be expressible: 0 steps means no territorial waters.
    const world = await generateWorld(makeParams({ points: 3000, territorialWaters: 0.01 }));
    const sea = world.params.seaLevel;
    const ocean = world.cells.filter((c) => c.height < sea && !isLakeCell(c));

    expect(ocean.filter((c) => c.regionId !== undefined)).toHaveLength(0);
  }, 60000);
});
