// Auto-detected & named geographic features (feature B3).
//
// Pure, read-only pass over the cell neighbor graph: BFS clusters cells into
// mountain ranges, deserts, forests, oceans/seas, islands, and (reused from B1)
// lakes, then names each deterministically from the TERRAIN seed. Detection
// never mutates cells — scratch state lives in local typed arrays only — so
// terrain/civ signatures stay byte-identical. All clustering is O(n): each
// predicate does one linear sweep with a flood fill over neighbors.

import { WorldData, GeoFeature, GeoFeatureKind, BiomeType, Point } from '../types';
import { normalizeVec } from './vec';
import { RNG } from './rng';
import { createNameGenerator, NameStyle } from './namegen';

// Minimum member counts per kind. Below these a cluster is noise, not a place.
const MIN_RANGE = 5;
const MIN_DESERT = 8;
const MIN_FOREST = 10;
const MIN_WATER = 15; // seas/oceans
const MAX_ISLAND = 25; // land components larger than this are continents, not islands

// A range is "well above" sea level: 55% of the way from the shoreline to the
// summit (height 1). Tuned to catch cordilleras without flagging gentle uplands.
const rangeThreshold = (seaLevel: number): number => seaLevel + 0.55 * (1 - seaLevel);

// Lakes sit at land height (>= seaLevel), so the ocean-water predicate excludes
// them automatically. Inlined here (rather than importing worldGen.isLakeCell)
// to avoid a features <-> worldGen import cycle.
const isLake = (biome: BiomeType): boolean =>
  biome === BiomeType.LAKE || biome === BiomeType.SALT_LAKE;

const FOREST_BIOMES = new Set<BiomeType>([
  BiomeType.TROPICAL_RAINFOREST,
  BiomeType.TEMPERATE_RAINFOREST,
  BiomeType.TEMPERATE_FOREST,
  BiomeType.BOREAL_FOREST,
]);

const DESERT_BIOMES = new Set<BiomeType>([
  BiomeType.HOT_DESERT,
  BiomeType.COLD_DESERT,
]);

export const detectFeatures = (world: WorldData): GeoFeature[] => {
  const { cells, params, lakes } = world;
  const seaLevel = params.seaLevel;
  const n = cells.length;

  // Flood-fill the neighbor graph into maximal connected components of every
  // cell satisfying `predicate`. Read-only: visited state is a local array.
  const components = (predicate: (id: number) => boolean): number[][] => {
    const visited = new Uint8Array(n);
    const comps: number[][] = [];
    for (let i = 0; i < n; i++) {
      if (visited[i] || !predicate(i)) continue;
      const comp: number[] = [];
      const stack: number[] = [i];
      visited[i] = 1;
      while (stack.length) {
        const id = stack.pop() as number;
        comp.push(id);
        for (const nb of cells[id].neighbors) {
          if (!visited[nb] && predicate(nb)) {
            visited[nb] = 1;
            stack.push(nb);
          }
        }
      }
      comps.push(comp);
    }
    return comps;
  };

  const anchorOf = (ids: number[]): Point => {
    let x = 0, y = 0, z = 0;
    for (const id of ids) {
      const c = cells[id].center;
      x += c.x; y += c.y; z += c.z;
    }
    const [nx, ny, nz] = normalizeVec([x, y, z]);
    return { x: nx, y: ny, z: nz };
  };

  const features: GeoFeature[] = [];
  const push = (kind: GeoFeatureKind, ids: number[]): void => {
    features.push({ id: features.length, kind, name: '', cellIds: ids, anchor: anchorOf(ids), size: ids.length });
  };

  // Ranges: contiguous high-elevation land (lakes excluded — they are water).
  const highT = rangeThreshold(seaLevel);
  for (const comp of components(i => !isLake(cells[i].biome) && cells[i].height > highT)) {
    if (comp.length >= MIN_RANGE) push('range', comp);
  }

  // Deserts and forests: contiguous biome families.
  for (const comp of components(i => DESERT_BIOMES.has(cells[i].biome))) {
    if (comp.length >= MIN_DESERT) push('desert', comp);
  }
  for (const comp of components(i => FOREST_BIOMES.has(cells[i].biome))) {
    if (comp.length >= MIN_FOREST) push('forest', comp);
  }

  // Water bodies: connected sub-sea-level cells. The largest 1-2 components are
  // oceans; the rest (>= MIN_WATER) are seas. Second-largest only counts as an
  // ocean if it is at least half the largest — otherwise a single dominant
  // ocean shouldn't drag a modest basin into the same class.
  const water = components(i => cells[i].height < seaLevel)
    .filter(c => c.length >= MIN_WATER)
    .sort((a, b) => b.length - a.length);
  water.forEach((comp, idx) => {
    const isOcean = idx === 0 || (idx === 1 && comp.length >= water[0].length * 0.5);
    push(isOcean ? 'ocean' : 'sea', comp);
  });

  // Islands: small connected landmasses (lakes excluded). Larger components are
  // continents and get no feature.
  for (const comp of components(i => cells[i].height >= seaLevel && !isLake(cells[i].biome))) {
    if (comp.length >= 1 && comp.length <= MAX_ISLAND) push('island', comp);
  }

  // Lakes: reuse the B1 entities 1:1 rather than re-detecting.
  for (const lake of lakes ?? []) {
    push('lake', lake.cellIds);
  }

  nameFeatures(features, params.seed, params.nameStyle ?? 'fantasy', lakes);
  return features;
};

// Deterministic naming from the TERRAIN seed (params.seed), NOT civSeed: re-
// rolling civilizations must never rename a mountain range. One RNG stream feeds
// both the base-name generator and the template-variant picks, so the whole pass
// is reproducible given the seed and a fixed feature order.
const nameFeatures = (
  features: GeoFeature[],
  seed: string,
  nameStyle: NameStyle,
  lakes: WorldData['lakes'],
): void => {
  const rng = new RNG(seed + '_geonames');
  const nameGen = createNameGenerator(nameStyle, () => rng.next());
  const pick = <T,>(variants: T[]): T => variants[Math.floor(rng.next() * variants.length)];

  const saltLakeCells = new Set<number>();
  for (const lake of lakes ?? []) {
    if (lake.isSalt) lake.cellIds.forEach(id => saltLakeCells.add(id));
  }

  for (const f of features) {
    const base = nameGen.town();
    switch (f.kind) {
      case 'range':
        f.name = pick([`${base} Mountains`, `${base} Range`]);
        break;
      case 'desert':
        f.name = pick([`${base} Desert`, `${base} Wastes`]);
        break;
      case 'forest':
        // 'Xwood' compounds the base directly, e.g. Eldoriawood.
        f.name = pick([`${base} Forest`, `${base}wood`]);
        break;
      case 'sea':
        f.name = pick([`Sea of ${base}`, `${base} Sea`]);
        break;
      case 'ocean':
        f.name = `${base} Ocean`;
        break;
      case 'island':
        f.name = pick([`${base} Isle`, `${base} Island`]);
        break;
      case 'lake': {
        // Salt basins occasionally read as 'Flats'; fresh lakes never do.
        const isSalt = f.cellIds.length > 0 && saltLakeCells.has(f.cellIds[0]);
        f.name = isSalt
          ? pick([`Lake ${base}`, `${base} Lake`, `${base} Flats`])
          : pick([`Lake ${base}`, `${base} Lake`]);
        break;
      }
    }
  }
};
