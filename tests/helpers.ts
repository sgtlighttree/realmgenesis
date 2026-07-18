import { WorldData, WorldParams } from '../types';

// Small, fast baseline for engine tests. Mirrors App.tsx DEFAULT_PARAMS with
// a reduced point count and one erosion pass so a full generation stays fast.
export const makeParams = (overrides: Partial<WorldParams> = {}): WorldParams => ({
  mapName: 'test',
  points: 300,
  planetRadius: 6371,
  axialTilt: 23.5,
  plates: 8,
  seaLevel: 0.55,
  roughness: 0.5,
  detailLevel: 3,
  landStyle: 'Continents',
  cellJitter: 0.5,
  noiseScale: 0.4,
  ridgeBlend: 0.1,
  mountainHeight: 1.0,
  oceanDepth: 1.0,
  maskType: 'None',
  warpStrength: 0.5,
  plateInfluence: 0.5,
  erosionIterations: 1,
  baseTemperature: 30,
  poleTemperature: -30,
  rainfallMultiplier: 1.0,
  moistureTransport: 0.5,
  temperatureVariance: 5,
  numFactions: 4,
  civSeed: 'test_civs',
  borderRoughness: 0.2,
  civSizeVariance: 0.5,
  waterCrossingCost: 0.8,
  territorialWaters: 0.15,
  capitalSpacing: 0.5,
  provinceSize: 0.5,
  nameStyle: 'fantasy',
  loreLevel: 1,
  seed: 'test_seed',
  ...overrides,
});

// Compact signature of everything the terrain pipeline produces per cell.
export const terrainSignature = (world: WorldData): string =>
  world.cells
    .map(c => `${c.height.toFixed(6)}|${c.biome}|${c.temperature.toFixed(3)}|${c.moisture.toFixed(6)}`)
    .join(';');

// Compact signature of the civilization layer.
export const civSignature = (world: WorldData): string =>
  world.cells
    .map(c => `${c.regionId ?? '-'}|${c.provinceId ?? '-'}|${c.population ?? 0}`)
    .join(';');

// Compact signature of every generated name (factions, provinces, towns).
export const nameSignature = (world: WorldData): string =>
  (world.civData?.factions ?? [])
    .map(f => `${f.name}[${f.provinces.map(p => `${p.name}:${p.towns.map(t => t.name).join(',')}`).join('|')}]`)
    .join(';');
