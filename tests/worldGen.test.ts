import { describe, it, expect } from 'vitest';
import { generateWorld } from '../utils/worldGen';
import { BiomeType } from '../types';
import { makeParams, terrainSignature, civSignature } from './helpers';

describe('generateWorld', () => {
  it('is fully deterministic for the same seed and params', async () => {
    const params = makeParams();
    const a = await generateWorld(params);
    const b = await generateWorld(params);
    expect(terrainSignature(a)).toBe(terrainSignature(b));
    expect(civSignature(a)).toBe(civSignature(b));
  }, 30000);

  it('produces a structurally sound world', async () => {
    const world = await generateWorld(makeParams());

    expect(world.cells.length).toBe(300);
    expect(world.geoJson.features.length).toBe(world.cells.length);

    world.cells.forEach((cell, i) => {
      expect(cell.id).toBe(i);
      expect(cell.height).toBeGreaterThanOrEqual(0);
      expect(cell.height).toBeLessThanOrEqual(1);
      expect(cell.moisture).toBeGreaterThanOrEqual(0);
      expect(cell.moisture).toBeLessThanOrEqual(1);
      expect(cell.neighbors.length).toBeGreaterThan(0);
      // Neighbor links are symmetric
      cell.neighbors.forEach(nId => {
        expect(world.cells[nId].neighbors).toContain(cell.id);
      });
      // Center is on the unit sphere
      const r = Math.hypot(cell.center.x, cell.center.y, cell.center.z);
      expect(r).toBeCloseTo(1, 6);
    });

    // Water/land classification matches heights
    world.cells.forEach(cell => {
      const isWaterBiome = cell.biome === BiomeType.OCEAN || cell.biome === BiomeType.DEEP_OCEAN;
      expect(isWaterBiome).toBe(cell.height < world.params.seaLevel);
    });
  }, 30000);

  it('generates civilizations with capitals on suitable land', async () => {
    const world = await generateWorld(makeParams());
    expect(world.civData).toBeDefined();
    const factions = world.civData!.factions;
    expect(factions.length).toBeGreaterThan(0);
    expect(factions.length).toBeLessThanOrEqual(4);

    factions.forEach(f => {
      const cap = world.cells[f.capitalId];
      expect(cap.isCapital).toBe(true);
      expect(cap.regionId).toBe(f.id);
      expect(cap.height).toBeGreaterThanOrEqual(world.params.seaLevel);
      expect(f.provinces.length).toBeGreaterThan(0);
      expect(f.totalPopulation).toBeGreaterThanOrEqual(0);
    });

    // No cell may carry a negative population (regression: high-elevation
    // suitability used to go negative)
    world.cells.forEach(cell => {
      expect(cell.population ?? 0).toBeGreaterThanOrEqual(0);
    });
  }, 30000);

  it('supports zero erosion iterations', async () => {
    const world = await generateWorld(makeParams({ erosionIterations: 0 }));
    expect(world.cells.length).toBe(300);
  }, 30000);

  it('reports progress reaching exactly the declared total', async () => {
    let last = 0;
    let total = 0;
    await generateWorld(makeParams(), undefined, undefined, (stage, t) => { last = stage; total = t; });
    expect(last).toBe(total);

    // The erosion tick must fire even when erosion is skipped
    last = 0;
    await generateWorld(makeParams({ erosionIterations: 0 }), undefined, undefined, (stage, t) => { last = stage; total = t; });
    expect(last).toBe(total);
  }, 60000);

  it('throws "Generation Cancelled" when aborted', async () => {
    const controller = new AbortController();
    controller.abort();
    await expect(generateWorld(makeParams(), undefined, controller.signal)).rejects.toThrow('Generation Cancelled');
  }, 30000);

  it('D2: currentStrength=0 is a byte-identical no-op; default-on changes climate', async () => {
    const climateSig = (w: Awaited<ReturnType<typeof generateWorld>>) =>
      w.cells.map(c => `${c.temperature.toFixed(6)}|${c.moisture.toFixed(6)}`).join(';');
    const off = climateSig(await generateWorld(makeParams({ currentStrength: 0 })));
    const off2 = climateSig(await generateWorld(makeParams({ currentStrength: 0 })));
    const on = climateSig(await generateWorld(makeParams({ currentStrength: 1.0 })));
    expect(off).toBe(off2);    // stage skipped + deterministic
    expect(on).not.toBe(off);  // default-on actually moves temperature + moisture
  }, 60000);
});
