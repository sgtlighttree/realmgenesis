import { describe, it, expect, beforeAll } from 'vitest';
import { generateWorld, isLakeCell, recalculateProvinces } from '../utils/worldGen';
import { settlementScore } from '../utils/settlements';
import { RNG } from '../utils/rng';
import { WorldData } from '../types';
import { makeParams } from './helpers';

// Task-4 pipeline wiring: placeSettlements runs inside recalculateCivs, the tier
// multiplier lands on cell.population (mirrored to Settlement.population), and
// provinces are seeded from the settlement set with towns[0] the province seat.
//
// The world is sized so the assertions BITE: at the default 300-cell/4-faction
// baseline every faction has quota 1 and every province one town, which would
// make the seat/identity checks vacuous. points:1000/numFactions:3/density:4
// yields multi-settlement, multi-province factions.
describe('settlement pipeline', () => {
  let world: WorldData;

  beforeAll(async () => {
    world = await generateWorld(makeParams({ points: 1000, numFactions: 3, settlementDensity: 4 }));
  }, 60000);

  it('produces a non-trivial settlement/province structure (guards vacuity)', () => {
    const civ = world.civData!;
    // more settlements than factions: the capital-only degenerate case is escaped
    expect(civ.settlements.length).toBeGreaterThan(civ.factions.length);
    // at least one province holds 2+ towns, so the seat-ordering check is real
    expect(civ.factions.some(f => f.provinces.some(p => p.towns.length >= 2))).toBe(true);
  });

  it('faction.totalPopulation equals Σ cell.population over its cells', () => {
    world.civData!.factions.forEach(f => {
      const sum = world.cells
        .filter(c => c.regionId === f.id)
        .reduce((s, c) => s + (c.population || 0), 0);
      expect(f.totalPopulation).toBe(sum);
    });
  });

  it('multiplied population stays within a bounded ratio of the base seeding', () => {
    // Replay the provRng '_prov' stream EXACTLY as recalculateProvinces does:
    // no .next() in the water branch. baseSum is the un-multiplied Σ; the live
    // total adds the tier multipliers on the small settlement subset.
    const rng = new RNG(world.params.civSeed + '_prov');
    let baseSum = 0;
    world.cells.forEach(c => {
      if (c.height < world.params.seaLevel || isLakeCell(c)) return;
      const suit = settlementScore(c, world, world.params);
      baseSum += Math.floor(suit * 10000 * (0.8 + rng.next() * 0.4));
    });
    const total = world.cells.reduce((s, c) => s + (c.population || 0), 0);
    const ratio = total / baseSum;
    // > 1: the multiplier must actually have fired (guards an empty byCell map).
    // < 3: settlements are a small fraction of land cells, so even with the
    // capital's ×8 the lift over base is bounded well under 3× (measured ~1.1–1.2).
    expect(baseSum).toBeGreaterThan(0);
    expect(ratio).toBeGreaterThan(1);
    expect(ratio).toBeLessThan(3);
  });

  it('towns[0] is the province seat (highest population in the province)', () => {
    world.civData!.factions.forEach(f => f.provinces.forEach(p => {
      expect(p.towns.length).toBeGreaterThan(0);
      const maxPop = Math.max(...p.towns.map(t => t.population));
      expect(p.towns[0].population).toBe(maxPop);
    }));
  });

  it('every province town is the SAME object as a CivData.settlements entry', () => {
    const set = new Set(world.civData!.settlements);
    world.civData!.factions.forEach(f => f.provinces.forEach(p =>
      p.towns.forEach(t => expect(set.has(t)).toBe(true))));
  });

  // S36 review findings 3 + 10. Political painting (applyPoliticalStroke) rewrites
  // cell.regionId WITHOUT updating Settlement.factionId, so after a repaint a
  // settlement can be listed under faction A while its cell belongs to B. The
  // attach loop used to read the cell's provinceId and push into A's province
  // array on a bounds check alone — a valid index into A's array, wrong faction.
  // Repainting ALL of a faction's settlements away also used to leave it with land
  // cells and zero provinces.
  it('keeps settlements with their cell\'s faction after a political repaint', () => {
    const w = structuredClone(world);
    const civ = w.civData!;
    const [a, b] = civ.factions;
    expect(b).toBeDefined();

    // Move every one of A's settlement cells to B, the way the political brush
    // does: cell.regionId only, Settlement.factionId deliberately untouched.
    const aSettlements = civ.settlements.filter(s => s.factionId === a.id);
    expect(aSettlements.length).toBeGreaterThan(0);
    for (const s of aSettlements) w.cells[s.cellId].regionId = b.id;

    recalculateProvinces(w, w.params);

    // No settlement may sit in a province of a faction its cell does not belong to.
    for (const f of civ.factions) {
      for (const prov of f.provinces) {
        for (const t of prov.towns) {
          expect(w.cells[t.cellId].regionId).toBe(f.id);
        }
      }
    }

    // Finding 10: a faction that still holds land must still hold a province.
    for (const f of civ.factions) {
      const land = w.cells.filter(c => c.regionId === f.id
        && c.height >= w.params.seaLevel && !isLakeCell(c));
      if (land.length === 0) continue;
      expect(f.provinces.length).toBeGreaterThanOrEqual(1);
    }
  }, 60000);

  // S36 regression. recalculateProvinces runs standalone from the province-size
  // slider (useWorldEngine.handleUpdateProvinces) but had no provinceId reset of
  // its own — only recalculateCivs did. Its Dijkstra refuses to relax into a cell
  // whose provinceId is already set, so a second standalone call died at the seeds
  // and each province accumulated only its own seed cell: measured ~87% collapse
  // in faction totals (6.87M -> 888k) while Sigma cell.population stayed 7.00M.
  // Pre-existing on main. Uses a private clone so the shared `world` above is not
  // repartitioned out from under the other cases.
  it('repartitions correctly when run standalone and repeatedly', () => {
    const w = structuredClone(world);
    const totals = () => w.civData!.factions.map(f => f.totalPopulation);
    const cellSum = () => w.cells.reduce((s, c) => s + (c.population || 0), 0);

    const baseline = totals();
    const baseSum = cellSum();
    expect(baseline.reduce((a, b) => a + b, 0)).toBeGreaterThan(baseSum * 0.5);

    recalculateProvinces(w, w.params);
    expect(totals()).toEqual(baseline);
    expect(cellSum()).toEqual(baseSum);

    // Idempotent across repeats, and across a province-count change and back.
    recalculateProvinces(w, w.params);
    expect(totals()).toEqual(baseline);

    recalculateProvinces(w, { ...w.params, provinceSize: 0.9 });
    const wide = totals();
    expect(wide.reduce((a, b) => a + b, 0)).toBeGreaterThan(baseSum * 0.5);

    // No cell may point at a province its faction no longer has.
    for (const c of w.cells) {
      if (c.regionId === undefined || c.provinceId === undefined) continue;
      const f = w.civData!.factions.find(x => x.id === c.regionId);
      if (!f) continue;
      expect(c.provinceId).toBeLessThan(f.provinces.length);
    }

    recalculateProvinces(w, w.params);
    expect(totals()).toEqual(baseline);
  }, 60000);
});
