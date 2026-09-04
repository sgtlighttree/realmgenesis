import { Cell, BiomeType, WorldData, WorldParams } from '../types';
import { RNG } from './rng';
import { createNameGenerator, NameGenerator, NameStyle } from './namegen';

// Lakes are land-height depressions, not ocean. Civ suitability, drainage, and
// coast checks must treat them as water even though height >= seaLevel.
export function isLakeCell(cell: Cell): boolean {
  return cell.biome === BiomeType.LAKE || cell.biome === BiomeType.SALT_LAKE;
}

// Looks up the naming style a cell's culture speaks, falling back to the
// world's base nameStyle for cells with no culture assignment yet (e.g. an
// oceanic world where recalculateCultures found zero land candidates).
export function getCultureNameStyle(world: WorldData, params: WorldParams, cellId: number): NameStyle {
  const cell = world.cells[cellId];
  const culture = cell?.cultureId !== undefined ? world.cultures?.find(c => c.id === cell.cultureId) : undefined;
  return culture?.nameStyle ?? (params.nameStyle ?? 'fantasy');
}

// One NameGenerator per NameStyle, lazily created on first use and cached for
// the life of the caller. Each style gets its own RNG side-stream (seedPrefix
// + style) so faction/province naming stays fully deterministic regardless of
// which cultures happen to touch which capitals/towns.
export function makeStyleNameGenCache(seedPrefix: string): (style: NameStyle) => NameGenerator {
  const cache = new Map<NameStyle, NameGenerator>();
  return (style: NameStyle) => {
    let gen = cache.get(style);
    if (!gen) {
      const rng = new RNG(seedPrefix + style);
      gen = createNameGenerator(style, () => rng.next());
      cache.set(style, gen);
    }
    return gen;
  };
}
