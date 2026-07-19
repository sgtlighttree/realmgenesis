import { describe, it, expect } from 'vitest';
import { generateWorld, isLakeCell } from '../utils/worldGen';
import { makeParams, religionSignature } from './helpers';

// C2: the religion layer. One 'folk' religion per culture (covers every
// land cell by default) plus a handful of 'organized' religions spreading
// from holy cities (top-population towns) that convert folk cells within a
// limited budget. Religions are computed on their own RNG side-stream
// (civSeed + '_religions') at the end of recalculateProvinces, after
// cultures/factions/provinces/towns are all in place. See
// tests/paramLiveness.test.ts, tests/cultures.test.ts and
// tests/worldGen.test.ts for the geometry-isolation assertions that must
// keep passing UNMODIFIED alongside these.

describe('religion layer (C2)', () => {
  it('is deterministic — same civSeed yields identical religionId assignment and religion roster', async () => {
    const params = makeParams({ civSeed: 'religions_determinism', points: 800 });
    const a = await generateWorld(params);
    const b = await generateWorld(params);

    expect(a.religions).toBeDefined();
    expect(a.religions!.length).toBeGreaterThan(0);
    expect(religionSignature(a)).toBe(religionSignature(b));
  }, 30000);

  it('assigns a religionId to every land cell and none to water cells', async () => {
    const world = await generateWorld(makeParams({ points: 800 }));
    const seaLevel = world.params.seaLevel;

    for (const cell of world.cells) {
      const isWater = cell.height < seaLevel || isLakeCell(cell);
      if (isWater) {
        expect(cell.religionId, `water cell ${cell.id} should have no religionId`).toBeUndefined();
      } else {
        expect(cell.religionId, `land cell ${cell.id} is missing a religionId`).toBeDefined();
      }
    }
  }, 30000);

  it('creates exactly one folk religion per culture, matching cultureId', async () => {
    const world = await generateWorld(makeParams({ numCultures: 6, points: 800 }));
    const folkReligions = world.religions!.filter(r => r.kind === 'folk');

    expect(folkReligions.length).toBe(world.cultures!.length);
    // Every culture has exactly one folk religion whose cultureId points back to it.
    const cultureIds = world.cultures!.map(c => c.id).sort((a, b) => a - b);
    const folkCultureIds = folkReligions.map(r => r.cultureId).sort((a, b) => (a ?? 0) - (b ?? 0));
    expect(folkCultureIds).toEqual(cultureIds);
    folkReligions.forEach(r => { expect(r.holyCellId).toBeNull(); });
  }, 30000);

  it('honors the organized-religion count formula: max(1, floor(numCultures / 2))', async () => {
    const world4 = await generateWorld(makeParams({ numCultures: 4, points: 800 }));
    const world8 = await generateWorld(makeParams({ numCultures: 8, points: 1200 }));
    const world1 = await generateWorld(makeParams({ numCultures: 1, points: 800 }));

    const organizedCount = (w: typeof world4) => w.religions!.filter(r => r.kind === 'organized').length;
    expect(organizedCount(world4)).toBe(2);
    expect(organizedCount(world8)).toBe(4);
    expect(organizedCount(world1)).toBe(1); // max(1, floor(1/2)) === 1
  }, 60000);

  it('seeds organized religions at town holy cities with minimum spacing', async () => {
    const world = await generateWorld(makeParams({ numCultures: 4, points: 800 }));
    const organizedReligions = world.religions!.filter(r => r.kind === 'organized');
    expect(organizedReligions.length).toBeGreaterThan(0);

    for (const r of organizedReligions) {
      expect(r.holyCellId).not.toBeNull();
      const holyCell = world.cells[r.holyCellId!];
      expect(holyCell.isTown).toBe(true);
      expect(r.cultureId).toBeNull();
    }

    // Minimum chord spacing between holy cities (mirrors capital/culture-home
    // spacing) — with 800 points and only 2 holy cities needed, the map has
    // ample well-spaced towns, so the no-spacing fallback should not fire.
    const HOLY_CITY_SPACING = 0.5;
    const minChordSq = HOLY_CITY_SPACING ** 2 * (4 / organizedReligions.length);
    for (let i = 0; i < organizedReligions.length; i++) {
      for (let j = i + 1; j < organizedReligions.length; j++) {
        const a = world.cells[organizedReligions[i].holyCellId!];
        const b = world.cells[organizedReligions[j].holyCellId!];
        const d = (a.center.x - b.center.x) ** 2 + (a.center.y - b.center.y) ** 2 + (a.center.z - b.center.z) ** 2;
        expect(d).toBeGreaterThanOrEqual(minChordSq);
      }
    }
  }, 30000);

  it('every organized religion claims cells including its own holy city', async () => {
    const world = await generateWorld(makeParams({ numCultures: 4, points: 800 }));
    const organizedReligions = world.religions!.filter(r => r.kind === 'organized');

    for (const r of organizedReligions) {
      const claimedCells = world.cells.filter(c => c.religionId === r.id);
      expect(claimedCells.length, `organized religion "${r.name}" claimed no cells`).toBeGreaterThan(0);
      expect(world.cells[r.holyCellId!].religionId, `holy city of "${r.name}" should belong to it`).toBe(r.id);
    }
  }, 30000);

  it('leaves a meaningful fraction of land as folk faith beyond organized budgets', async () => {
    // Large-ish fixture so the per-religion cell-count budget (~half an
    // equal land share) has room to plateau well short of full coverage.
    const world = await generateWorld(makeParams({ numCultures: 4, points: 1500 }));
    const seaLevel = world.params.seaLevel;
    const landCells = world.cells.filter(c => c.height >= seaLevel && !isLakeCell(c));
    const folkReligionIds = new Set(world.religions!.filter(r => r.kind === 'folk').map(r => r.id));
    const folkCells = landCells.filter(c => c.religionId !== undefined && folkReligionIds.has(c.religionId));

    expect(folkCells.length, 'expected some land to remain folk faith').toBeGreaterThan(0);
    // Budgets are sized so organized faiths cover roughly half the land at
    // most — folk should retain a non-trivial share, not be nearly wiped out.
    expect(folkCells.length / landCells.length).toBeGreaterThan(0.15);
  }, 30000);

  it('recalculateReligions never perturbs terrain or civ geometry', async () => {
    // Religions run on their own civSeed + '_religions' stream and never
    // touch civRng/provRng/cultureRng or any cell field but religionId, so
    // changing numCultures (which reshapes the religion layer through the
    // culture/organized-count formula, and legitimately reshapes cultures
    // too — see tests/cultures.test.ts) must not perturb terrain or civ
    // (regionId/provinceId/population) geometry at all. The dedicated
    // liveness assertions already cover this in tests/paramLiveness.test.ts
    // (run unmodified); this is a direct spot-check alongside them.
    const { terrainSignature, civSignature } = await import('./helpers');
    const baseline = await generateWorld(makeParams({ numCultures: 4 }));
    const changed = await generateWorld(makeParams({ numCultures: 8 }));

    expect(religionSignature(changed)).not.toBe(religionSignature(baseline));
    expect(civSignature(changed)).toBe(civSignature(baseline));
    expect(terrainSignature(changed)).toBe(terrainSignature(baseline));
  }, 30000);
});
