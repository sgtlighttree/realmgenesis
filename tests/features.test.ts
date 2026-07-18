import { describe, it, expect } from 'vitest';
import { generateWorld, recalculateCivs } from '../utils/worldGen';
import { GeoFeature, BiomeType, WorldData } from '../types';
import { makeParams } from './helpers';

// Same 'lakeworld' seed the B1 lake suite pins: it yields exactly one salt lake,
// so the lakes-map-1:1 assertion has something concrete to bite on.
const LAKE_SEED = 'lakeworld';

const FOREST_BIOMES = new Set<BiomeType>([
  BiomeType.TROPICAL_RAINFOREST,
  BiomeType.TEMPERATE_RAINFOREST,
  BiomeType.TEMPERATE_FOREST,
  BiomeType.BOREAL_FOREST,
]);
const DESERT_BIOMES = new Set<BiomeType>([BiomeType.HOT_DESERT, BiomeType.COLD_DESERT]);
const isLake = (b: BiomeType) => b === BiomeType.LAKE || b === BiomeType.SALT_LAKE;

const featureSig = (features: GeoFeature[]): string =>
  features
    .map(f => `${f.kind}|${f.name}|${f.cellIds.join(',')}|` +
      `${f.anchor.x.toFixed(6)},${f.anchor.y.toFixed(6)},${f.anchor.z.toFixed(6)}|${f.size}`)
    .join(';');

describe('geographic features (B3)', () => {
  it('is deterministic — same seed yields identical features and names', async () => {
    const a = await generateWorld(makeParams({ seed: LAKE_SEED }));
    const b = await generateWorld(makeParams({ seed: LAKE_SEED }));
    expect(a.features).toBeDefined();
    expect(a.features!.length).toBeGreaterThan(0);
    expect(featureSig(a.features!)).toBe(featureSig(b.features!));
  }, 30000);

  it('classifies every feature consistently with its cells', async () => {
    const world = await generateWorld(makeParams({ seed: LAKE_SEED }));
    const seaLevel = world.params.seaLevel;
    const highT = seaLevel + 0.55 * (1 - seaLevel);
    const cells = world.cells;

    for (const f of world.features!) {
      expect(f.cellIds.length).toBeGreaterThan(0);
      switch (f.kind) {
        case 'range':
          for (const id of f.cellIds) {
            expect(cells[id].height).toBeGreaterThan(highT);
            expect(isLake(cells[id].biome)).toBe(false);
          }
          break;
        case 'desert':
          for (const id of f.cellIds) expect(DESERT_BIOMES.has(cells[id].biome)).toBe(true);
          break;
        case 'forest':
          for (const id of f.cellIds) expect(FOREST_BIOMES.has(cells[id].biome)).toBe(true);
          break;
        case 'sea':
        case 'ocean':
          for (const id of f.cellIds) expect(cells[id].height).toBeLessThan(seaLevel);
          break;
        case 'island':
          expect(f.cellIds.length).toBeLessThanOrEqual(25);
          for (const id of f.cellIds) {
            expect(cells[id].height).toBeGreaterThanOrEqual(seaLevel);
            expect(isLake(cells[id].biome)).toBe(false);
          }
          break;
        case 'lake':
          for (const id of f.cellIds) expect(isLake(cells[id].biome)).toBe(true);
          break;
      }
    }
  }, 30000);

  it('maps lake features 1:1 onto world.lakes', async () => {
    const world = await generateWorld(makeParams({ seed: LAKE_SEED }));
    const lakeFeatures = world.features!.filter(f => f.kind === 'lake');
    expect(lakeFeatures.length).toBe(world.lakes!.length);
    expect(lakeFeatures.length).toBeGreaterThan(0);
    // Each lake feature reuses a lake's exact cell set.
    const lakeSets = world.lakes!.map(l => l.cellIds.join(','));
    for (const f of lakeFeatures) {
      expect(lakeSets).toContain(f.cellIds.join(','));
    }
  }, 30000);

  it('names every feature with a kind-appropriate pattern', async () => {
    const world = await generateWorld(makeParams({ seed: LAKE_SEED }));
    for (const f of world.features!) {
      expect(f.name.length).toBeGreaterThan(0);
      switch (f.kind) {
        case 'range':
          expect(/(Mountains|Range)$/.test(f.name)).toBe(true);
          break;
        case 'desert':
          expect(/(Desert|Wastes)$/.test(f.name)).toBe(true);
          break;
        case 'forest':
          expect(/(Forest|wood)$/.test(f.name)).toBe(true);
          break;
        case 'sea':
          expect(/(^Sea of |Sea$)/.test(f.name)).toBe(true);
          break;
        case 'ocean':
          expect(/Ocean$/.test(f.name)).toBe(true);
          break;
        case 'island':
          expect(/(Isle|Island)$/.test(f.name)).toBe(true);
          break;
        case 'lake':
          expect(/(^Lake |Lake$|Flats$)/.test(f.name)).toBe(true);
          break;
      }
    }
  }, 30000);

  it('leaves features untouched when civs are re-rolled with a new civSeed', async () => {
    const world: WorldData = await generateWorld(makeParams({ seed: LAKE_SEED, civSeed: 'civs_a' }));
    const before = featureSig(world.features!);
    // Re-roll only the civ layer with a different civSeed — mountains must not
    // rename (ROADMAP requirement) and no feature may shift.
    recalculateCivs(world, { ...world.params, civSeed: 'civs_b' });
    expect(world.features).toBeDefined();
    expect(featureSig(world.features!)).toBe(before);
  }, 30000);
});
