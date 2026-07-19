import { WorldData } from '../types';

// Surgical editors for the Civ tab's "editor completeness" operations (C5):
// merge, rename, capital relocation. Unlike recalculateProvinces (worldGen.ts),
// these never regenerate names or boundaries — they mutate the passed world
// in place and let the caller shallow-copy + setWorld, matching the paint/edit
// convention used throughout App.tsx.

/**
 * Merges srcId's provinces, towns, and cells into dstId. src is left with no
 * provinces and removed from civData.factions entirely.
 *
 * Returns false (no-op) when srcId === dstId or either faction is missing.
 */
export function mergeFactions(world: WorldData, srcId: number, dstId: number): boolean {
  if (!world.civData || srcId === dstId) return false;
  const src = world.civData.factions.find(f => f.id === srcId);
  const dst = world.civData.factions.find(f => f.id === dstId);
  if (!src || !dst) return false;

  // Local province ids are per-faction and (per Inspector.tsx's
  // faction.provinces[cell.provinceId] array-index lookup) must stay
  // contiguous 0..n-1 within that array. Appending src's provinces onto
  // dst's array and renumbering them to dst's next free ids preserves that
  // invariant. Build the FULL old-id -> new-id map before touching any
  // cell — src and dst both start their own province ids at 0, so
  // remapping cells mid-build would collide/double-apply.
  const provinceIdMap = new Map<number, number>();
  src.provinces.forEach(province => {
    const oldId = province.id;
    const newId = dst.provinces.length;
    province.id = newId;
    dst.provinces.push(province);
    provinceIdMap.set(oldId, newId);
  });
  src.provinces = [];

  world.cells.forEach(cell => {
    if (cell.regionId !== srcId) return;
    cell.regionId = dstId;
    if (cell.provinceId !== undefined) {
      cell.provinceId = provinceIdMap.get(cell.provinceId) ?? cell.provinceId;
    }
  });

  // isCapital is a dual flag (Cell + TownData) — demote src's capital on
  // both. dst's capital is untouched.
  const srcCapitalTown = dst.provinces
    .flatMap(p => p.towns)
    .find(t => t.cellId === src.capitalId);
  if (srcCapitalTown) srcCapitalTown.isCapital = false;
  if (world.cells[src.capitalId]) world.cells[src.capitalId].isCapital = false;

  dst.totalPopulation += src.totalPopulation;
  world.civData.factions = world.civData.factions.filter(f => f.id !== srcId);
  return true;
}

/** Renames a province. Rejects (returns false) empty/whitespace-only names. */
export function renameProvince(world: WorldData, factionId: number, provinceId: number, name: string): boolean {
  const trimmed = name.trim();
  if (!trimmed) return false;
  const faction = world.civData?.factions.find(f => f.id === factionId);
  const province = faction?.provinces.find(p => p.id === provinceId);
  if (!province) return false;
  province.name = trimmed;
  return true;
}

/** Renames a town. Rejects (returns false) empty/whitespace-only names. */
export function renameTown(
  world: WorldData,
  factionId: number,
  provinceId: number,
  cellId: number,
  name: string,
): boolean {
  const trimmed = name.trim();
  if (!trimmed) return false;
  const faction = world.civData?.factions.find(f => f.id === factionId);
  const province = faction?.provinces.find(p => p.id === provinceId);
  const town = province?.towns.find(t => t.cellId === cellId);
  if (!town) return false;
  town.name = trimmed;
  return true;
}

/**
 * Promotes an existing town of the faction to capital. Returns false if
 * townCellId isn't a town of that faction, or is already the capital.
 */
export function relocateCapital(world: WorldData, factionId: number, townCellId: number): boolean {
  const faction = world.civData?.factions.find(f => f.id === factionId);
  if (!faction) return false;
  if (townCellId === faction.capitalId) return false;

  const allTowns = faction.provinces.flatMap(p => p.towns);
  const newCapitalTown = allTowns.find(t => t.cellId === townCellId);
  if (!newCapitalTown) return false;

  const oldCapitalTown = allTowns.find(t => t.cellId === faction.capitalId);
  if (oldCapitalTown) oldCapitalTown.isCapital = false;
  if (world.cells[faction.capitalId]) world.cells[faction.capitalId].isCapital = false;

  newCapitalTown.isCapital = true;
  if (world.cells[townCellId]) world.cells[townCellId].isCapital = true;
  faction.capitalId = townCellId;
  return true;
}
