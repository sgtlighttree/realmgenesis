import { describe, it, expect } from 'vitest';
import { buildMapConfig, validateCivData, validateWorldParams } from '../utils/export';
import { BiomeType, Cell, GeoJsonCollection, WorldData, WorldParams, Settlement, SettlementTier } from '../types';

// Hand-built fixture: the save payload only needs cells + civData, so we avoid
// paying for a real generateWorld() run (CPU-constrained machine).
const makeCell = (id: number, overrides: Partial<Cell> = {}): Cell => ({
  id,
  center: { x: id, y: 0, z: 0 },
  vertices: [],
  neighbors: [],
  height: 0.6,
  plateId: 0,
  temperature: 15,
  moisture: 0.5,
  biome: BiomeType.GRASSLAND,
  regionId: 0,
  provinceId: 0,
  ...overrides,
});

const makeSettlement = (id: number, cellId: number, tier: SettlementTier, overrides: Partial<Settlement> = {}): Settlement => ({
  id,
  cellId,
  name: `Settlement ${id}`,
  population: 100,
  tier,
  factionId: 0,
  isCapital: tier === 'capital',
  isPort: false,
  ...overrides,
});

const makeWorld = (towns: Settlement[]): WorldData => {
  const maxCellId = Math.max(...towns.map(t => t.cellId));
  const cells: Cell[] = [];
  for (let id = 0; id <= maxCellId; id++) cells[id] = makeCell(id);
  return {
    cells,
    params: { seaLevel: 0.3 } as WorldParams,
    geoJson: { type: 'FeatureCollection', features: [] } as GeoJsonCollection,
    civData: {
      factions: [
        {
          id: 0,
          name: 'Test Faction',
          color: '#fff',
          capitalId: towns.find(t => t.tier === 'capital')?.cellId ?? towns[0].cellId,
          provinces: [
            {
              id: 0,
              name: 'Test Province',
              towns,
              totalPopulation: towns.reduce((s, t) => s + t.population, 0),
            },
          ],
          totalPopulation: towns.reduce((s, t) => s + t.population, 0),
        },
      ],
      settlements: towns,
    },
  } as WorldData;
};

describe('buildMapConfig', () => {
  const towns = [makeSettlement(0, 0, 'capital'), makeSettlement(1, 1, 'town')];
  const world = makeWorld(towns);
  const params = { seed: 'abc' } as WorldParams;

  it('does not serialize settlements (they regenerate from the seed)', () => {
    const payload = buildMapConfig(params, world);
    const json = JSON.parse(JSON.stringify(payload));
    expect(json.civData.settlements).toBeUndefined();
  });

  it('bumps the schema version to 1.5', () => {
    const payload = buildMapConfig(params, world);
    expect(payload.version).toBe('1.5');
  });

  it('keeps the nested town names intact after the strip', () => {
    const payload = buildMapConfig(params, world);
    const json = JSON.parse(JSON.stringify(payload));
    expect(json.civData.factions[0].provinces[0].towns[0].name).toBe(towns[0].name);
  });

  it('produces a civData that still passes validateCivData', () => {
    const payload = buildMapConfig(params, world);
    const json = JSON.parse(JSON.stringify(payload));
    expect(validateCivData(json.civData)).toBe(true);
  });
});

describe('validateWorldParams settlementDensity bound', () => {
  it('accepts settlementDensity within [0, 2]', () => {
    expect(validateWorldParams({ settlementDensity: 1.0 })).toBe(true);
  });

  it('rejects settlementDensity out of bounds', () => {
    expect(validateWorldParams({ settlementDensity: 5 })).toBe(false);
  });
});
