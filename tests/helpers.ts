import { WorldData, WorldParams } from '../types';
import { DEFAULT_PARAMS } from '../utils/defaultParams';

// Small, fast baseline for engine tests. Spreads the app's canonical
// DEFAULT_PARAMS (single source of truth — divergence is now structurally
// impossible) and overrides only the handful of keys that keep a full
// generation fast: a reduced point count, fewer plates/factions, one erosion
// pass, and stable test seeds.
export const makeParams = (overrides: Partial<WorldParams> = {}): WorldParams => ({
  ...DEFAULT_PARAMS,
  mapName: 'test',
  points: 300,
  plates: 8,
  erosionIterations: 1,
  numFactions: 4,
  civSeed: 'test_civs',
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

// Compact signature of the culture layer (C1) — per-cell cultureId, plus the
// culture roster's own names/styles so home-cell/naming changes show up too.
export const cultureSignature = (world: WorldData): string =>
  world.cells.map(c => `${c.cultureId ?? '-'}`).join(';') +
  '|' +
  (world.cultures ?? []).map(c => `${c.name}:${c.nameStyle}:${c.homeCellId}`).join(',');

// Compact signature of the religion layer (C2) — per-cell religionId, plus
// the religion roster's own identity fields so naming/holy-city changes
// show up too.
export const religionSignature = (world: WorldData): string =>
  world.cells.map(c => `${c.religionId ?? '-'}`).join(';') +
  '|' +
  (world.religions ?? []).map(r => `${r.name}:${r.kind}:${r.cultureId ?? '-'}:${r.holyCellId ?? '-'}`).join(',');
