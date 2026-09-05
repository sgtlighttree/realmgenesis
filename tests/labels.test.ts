import { describe, it, expect } from 'vitest';
import { collectLabels, LABEL_CONFIG, ZOOM_THRESHOLDS } from '../utils/labels';
import { BiomeType, Cell, GeoJsonCollection, WorldData, WorldParams, Settlement, SettlementTier } from '../types';

// Hand-built fixture: label emission only needs cells + civData, so we avoid
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
  };
};

describe('labels — tier kinds', () => {
  it('emits kind "city" for a city-tier settlement', () => {
    const town = makeSettlement(1, 1, 'city');
    const world = makeWorld([town]);
    const labels = collectLabels(world);
    const label = labels.find(l => l.name === town.name);
    expect(label).toBeDefined();
    expect(label!.kind).toBe('city');
  });

  it('emits kind "capital" for a capital-tier settlement', () => {
    const town = makeSettlement(1, 1, 'capital');
    const world = makeWorld([town]);
    const labels = collectLabels(world);
    const label = labels.find(l => l.name === town.name);
    expect(label).toBeDefined();
    expect(label!.kind).toBe('capital');
  });

  it('emits kind "village" for a village-tier settlement', () => {
    const town = makeSettlement(1, 1, 'village');
    const world = makeWorld([town]);
    const labels = collectLabels(world);
    const label = labels.find(l => l.name === town.name);
    expect(label).toBeDefined();
    expect(label!.kind).toBe('village');
  });

  it('emits kind "hamlet" for a hamlet-tier settlement', () => {
    const town = makeSettlement(1, 1, 'hamlet');
    const world = makeWorld([town]);
    const labels = collectLabels(world);
    const label = labels.find(l => l.name === town.name);
    expect(label).toBeDefined();
    expect(label!.kind).toBe('hamlet');
  });
});

describe('labels — ZOOM_THRESHOLDS ordering', () => {
  it('village zooms in later than town', () => {
    expect(ZOOM_THRESHOLDS.village).toBeGreaterThan(ZOOM_THRESHOLDS.town);
  });

  it('hamlet zooms in later than village', () => {
    expect(ZOOM_THRESHOLDS.hamlet).toBeGreaterThan(ZOOM_THRESHOLDS.village);
  });
});

describe('labels — LABEL_CONFIG visibility keys', () => {
  it('maps each tier kind to its own visibility toggle', () => {
    expect(LABEL_CONFIG.city.visibilityKey).toBe('cities');
    expect(LABEL_CONFIG.village.visibilityKey).toBe('villages');
    expect(LABEL_CONFIG.hamlet.visibilityKey).toBe('hamlets');
  });
});
