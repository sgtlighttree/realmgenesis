import { describe, it, expect } from 'vitest';
import { generateWorld } from '../utils/worldGen';
import { makeParams } from './helpers';
import { mergeFactions, renameProvince, renameTown, relocateCapital } from '../utils/civEdit';
import { WorldData } from '../types';

// Reproduces App.tsx's recalculatePoliticalTotals inline (rather than
// importing App.tsx) — matches the "caller recomputes totals after a
// civEdit mutation" convention the handlers in App.tsx follow.
function recalcTotals(world: WorldData): void {
  if (!world.civData) return;
  const factionById = new Map(world.civData.factions.map(f => [f.id, f]));
  world.civData.factions.forEach(faction => {
    faction.totalPopulation = 0;
    faction.provinces.forEach(p => { p.totalPopulation = 0; });
  });
  world.cells.forEach(cell => {
    if (cell.regionId === undefined || cell.provinceId === undefined) return;
    const faction = factionById.get(cell.regionId);
    const province = faction?.provinces.find(p => p.id === cell.provinceId);
    if (!faction || !province) return;
    province.totalPopulation += cell.population || 0;
  });
  world.civData.factions.forEach(faction => {
    faction.totalPopulation = faction.provinces.reduce((sum, p) => sum + p.totalPopulation, 0);
  });
}

// Denser provinces than the default helper params so at least one faction
// reliably ends up with several provinces — needed to exercise local
// province-id remapping during a merge.
const civParams = () => makeParams({ points: 1200, provinceSize: 0.1, numFactions: 3, landStyle: 'Continents' });

describe('civEdit', () => {
  describe('mergeFactions', () => {
    it('merges provinces with remapped local ids and rewrites every affected cell', async () => {
      const world = await generateWorld(civParams());
      const factions = world.civData!.factions;

      const sortedByProvinceCount = [...factions].sort((a, b) => b.provinces.length - a.provinces.length);
      const dstFaction = sortedByProvinceCount[0];
      const srcFaction = sortedByProvinceCount.find(f => f.id !== dstFaction.id && f.provinces.length >= 2);
      // Precondition: the test needs a source faction with >=2 provinces to
      // actually exercise id remapping (a single-province merge can't tell
      // "remapped" from "coincidentally unchanged"). If params drift and this
      // stops holding, fail loudly instead of silently under-testing.
      expect(srcFaction, 'expected a faction with >=2 provinces to merge away').toBeDefined();
      const src = srcFaction!;
      const dst = dstFaction;

      const dstProvinceCountBefore = dst.provinces.length;
      const srcProvinceNamesInOrder = src.provinces.map(p => p.name);
      const srcCapitalId = src.capitalId;
      const dstCapitalId = dst.capitalId;

      // Snapshot every src cell's provinceId before the merge so we can
      // verify the exact remap afterwards.
      const srcCellSnapshots = world.cells
        .filter(c => c.regionId === src.id)
        .map(c => ({ id: c.id, oldProvinceId: c.provinceId }));
      expect(srcCellSnapshots.length).toBeGreaterThan(0);

      const ok = mergeFactions(world, src.id, dst.id);
      expect(ok).toBe(true);

      // src is gone
      expect(world.civData!.factions.find(f => f.id === src.id)).toBeUndefined();

      const mergedDst = world.civData!.factions.find(f => f.id === dst.id)!;
      expect(mergedDst).toBeDefined();
      expect(mergedDst.provinces.length).toBe(dstProvinceCountBefore + srcProvinceNamesInOrder.length);

      // Local province ids stay contiguous 0..n-1 matching array index
      // (Inspector.tsx indexes faction.provinces[cell.provinceId] directly).
      mergedDst.provinces.forEach((p, idx) => expect(p.id).toBe(idx));

      // The appended tail holds src's provinces, in their original order,
      // renumbered starting at dst's old province count.
      const appendedNames = mergedDst.provinces.slice(dstProvinceCountBefore).map(p => p.name);
      expect(appendedNames).toEqual(srcProvinceNamesInOrder);

      // Every src cell now belongs to dst, with its provinceId shifted by
      // exactly dst's pre-merge province count.
      srcCellSnapshots.forEach(({ id, oldProvinceId }) => {
        const cell = world.cells[id];
        expect(cell.regionId).toBe(dst.id);
        if (oldProvinceId !== undefined) {
          expect(cell.provinceId).toBe(oldProvinceId + dstProvinceCountBefore);
        }
      });

      // src's capital is demoted (both the cell flag and the TownData copy
      // that now lives inside dst's provinces); dst's capital is untouched.
      expect(world.cells[srcCapitalId].isCapital).toBe(false);
      const demotedTown = mergedDst.provinces.flatMap(p => p.towns).find(t => t.cellId === srcCapitalId);
      expect(demotedTown?.isCapital).toBe(false);

      expect(mergedDst.capitalId).toBe(dstCapitalId);
      expect(world.cells[dstCapitalId].isCapital).toBe(true);
    }, 30000);

    it('no-ops when src and dst are the same id', async () => {
      const world = await generateWorld(civParams());
      const factionId = world.civData!.factions[0].id;
      const factionsBefore = world.civData!.factions;
      const ok = mergeFactions(world, factionId, factionId);
      expect(ok).toBe(false);
      expect(world.civData!.factions).toBe(factionsBefore);
    }, 30000);

    it('no-ops when either faction id does not exist', async () => {
      const world = await generateWorld(civParams());
      const realId = world.civData!.factions[0].id;
      const badId = -999;
      expect(mergeFactions(world, badId, realId)).toBe(false);
      expect(mergeFactions(world, realId, badId)).toBe(false);
      expect(mergeFactions(world, badId, badId)).toBe(false);
    }, 30000);

    it('lets the caller recompute totals after a merge (recalculatePoliticalTotals pattern)', async () => {
      const world = await generateWorld(civParams());
      const factions = world.civData!.factions;
      const sortedByProvinceCount = [...factions].sort((a, b) => b.provinces.length - a.provinces.length);
      const dst = sortedByProvinceCount[0];
      const src = sortedByProvinceCount.find(f => f.id !== dst.id)!;

      expect(mergeFactions(world, src.id, dst.id)).toBe(true);
      recalcTotals(world);

      world.civData!.factions.forEach(faction => {
        const expectedTotal = world.cells
          .filter(c => c.regionId === faction.id)
          .reduce((sum, c) => sum + (c.population || 0), 0);
        expect(faction.totalPopulation).toBe(expectedTotal);
      });
    }, 30000);
  });

  describe('renameProvince', () => {
    it('renames an existing province', async () => {
      const world = await generateWorld(civParams());
      const faction = world.civData!.factions[0];
      const province = faction.provinces[0];
      const ok = renameProvince(world, faction.id, province.id, '  New Province Name  ');
      expect(ok).toBe(true);
      expect(province.name).toBe('New Province Name');
    }, 30000);

    it('rejects an empty or whitespace-only name', async () => {
      const world = await generateWorld(civParams());
      const faction = world.civData!.factions[0];
      const province = faction.provinces[0];
      const originalName = province.name;
      expect(renameProvince(world, faction.id, province.id, '   ')).toBe(false);
      expect(renameProvince(world, faction.id, province.id, '')).toBe(false);
      expect(province.name).toBe(originalName);
    }, 30000);

    it('rejects unknown faction or province ids', async () => {
      const world = await generateWorld(civParams());
      const faction = world.civData!.factions[0];
      expect(renameProvince(world, -999, faction.provinces[0].id, 'X')).toBe(false);
      expect(renameProvince(world, faction.id, -999, 'X')).toBe(false);
    }, 30000);
  });

  describe('renameTown', () => {
    it('renames an existing town', async () => {
      const world = await generateWorld(civParams());
      const faction = world.civData!.factions[0];
      const province = faction.provinces[0];
      const town = province.towns[0];
      const ok = renameTown(world, faction.id, province.id, town.cellId, '  New Town  ');
      expect(ok).toBe(true);
      expect(town.name).toBe('New Town');
    }, 30000);

    it('rejects an empty or whitespace-only name', async () => {
      const world = await generateWorld(civParams());
      const faction = world.civData!.factions[0];
      const province = faction.provinces[0];
      const town = province.towns[0];
      const originalName = town.name;
      expect(renameTown(world, faction.id, province.id, town.cellId, '   ')).toBe(false);
      expect(town.name).toBe(originalName);
    }, 30000);

    it('rejects an unknown cellId within an otherwise valid province', async () => {
      const world = await generateWorld(civParams());
      const faction = world.civData!.factions[0];
      const province = faction.provinces[0];
      expect(renameTown(world, faction.id, province.id, -999, 'X')).toBe(false);
    }, 30000);
  });

  describe('relocateCapital', () => {
    it('moves the capital flag pair and updates faction.capitalId', async () => {
      const world = await generateWorld(civParams());
      const faction = world.civData!.factions.find(f => f.provinces.length >= 2)!;
      expect(faction, 'expected a faction with a non-capital town to relocate to').toBeDefined();

      const oldCapitalId = faction.capitalId;
      // Any province other than the capital's (province 0) holds a plain,
      // non-capital town in this generator (one town per province).
      const targetTown = faction.provinces.find(p => p.id !== 0)!.towns[0];
      expect(targetTown.isCapital).toBe(false);

      const ok = relocateCapital(world, faction.id, targetTown.cellId);
      expect(ok).toBe(true);

      expect(faction.capitalId).toBe(targetTown.cellId);
      expect(world.cells[targetTown.cellId].isCapital).toBe(true);
      expect(targetTown.isCapital).toBe(true);

      expect(world.cells[oldCapitalId].isCapital).toBe(false);
      const oldTown = faction.provinces.flatMap(p => p.towns).find(t => t.cellId === oldCapitalId);
      expect(oldTown?.isCapital).toBe(false);
    }, 30000);

    it('rejects a target that is not a town of that faction', async () => {
      const world = await generateWorld(civParams());
      const faction = world.civData!.factions[0];
      const nonTownCell = world.cells.find(c => c.regionId === faction.id && !c.isTown);
      expect(nonTownCell, 'expected at least one non-town cell in the faction').toBeDefined();
      const capitalBefore = faction.capitalId;
      expect(relocateCapital(world, faction.id, nonTownCell!.id)).toBe(false);
      expect(faction.capitalId).toBe(capitalBefore);
    }, 30000);

    it('rejects a target that is already the capital', async () => {
      const world = await generateWorld(civParams());
      const faction = world.civData!.factions[0];
      expect(relocateCapital(world, faction.id, faction.capitalId)).toBe(false);
    }, 30000);

    it('rejects an unknown faction id', async () => {
      const world = await generateWorld(civParams());
      const faction = world.civData!.factions[0];
      const someTown = faction.provinces[0].towns[0];
      expect(relocateCapital(world, -999, someTown.cellId)).toBe(false);
    }, 30000);
  });
});
