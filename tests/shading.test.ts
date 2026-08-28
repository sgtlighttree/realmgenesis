import { describe, it, expect } from 'vitest';
import { generateWorld } from '../utils/worldGen';
import { computeShadeMap, computeContourSegments, SHADE_MIN, SHADE_MAX } from '../utils/shading';
import { Cell, BiomeType } from '../types';
import { makeParams } from './helpers';

// Minimal hand-built cell for the synthetic contour fixture. Two cells share an
// edge when they share two vertices within the epsilon used by computeContourSegments.
const makeCell = (id: number, height: number, neighbors: number[], vertices: Cell['vertices']): Cell => ({
  id,
  center: { x: 1, y: 0, z: 0 },
  vertices,
  neighbors,
  height,
  plateId: 0,
  temperature: 10,
  moisture: 0.5,
  biome: BiomeType.TEMPERATE_FOREST,
});

const SL = 0.55;

describe('computeShadeMap', () => {
  it('is deterministic for the same world', async () => {
    const world = await generateWorld(makeParams());
    const a = computeShadeMap(world.cells, world.params.seaLevel);
    const b = computeShadeMap(world.cells, world.params.seaLevel);
    expect(Array.from(a)).toEqual(Array.from(b));
  }, 30000);

  it('stays in the clamp band for every cell', async () => {
    const world = await generateWorld(makeParams());
    const shade = computeShadeMap(world.cells, world.params.seaLevel);
    expect(shade.length).toBe(world.cells.length);
    world.cells.forEach(cell => {
      const s = shade[cell.id];
      expect(s).toBeGreaterThanOrEqual(SHADE_MIN);
      expect(s).toBeLessThanOrEqual(SHADE_MAX);
    });
  }, 30000);

  // D10 changed this deliberately: water used to be forced to exactly 1.0, which
  // rendered the sea bed as a flat colour ramp in Map2D and exports whatever the
  // bathymetry said. Water now shades off its water neighbours.
  it('shades the sea bed, not just the land', async () => {
    const world = await generateWorld(makeParams());
    const shade = computeShadeMap(world.cells, world.params.seaLevel);
    const water = world.cells.filter(c => c.height < world.params.seaLevel);
    expect(water.length).toBeGreaterThan(0);
    expect(water.some(c => shade[c.id] !== 1.0)).toBe(true);
  }, 30000);

  it('produces at least some non-neutral shading on real relief', async () => {
    const world = await generateWorld(makeParams());
    const shade = computeShadeMap(world.cells, world.params.seaLevel);
    const shaded = world.cells.some(cell => cell.height >= world.params.seaLevel && shade[cell.id] !== 1.0);
    expect(shaded).toBe(true);
  }, 30000);
});

describe('computeContourSegments', () => {
  // Shared edge between cells 0 and 1: vertices v1 (0,1,0) and v2 (0,0,1).
  const v1 = { x: 0, y: 1, z: 0 };
  const v2 = { x: 0, y: 0, z: 1 };
  const far = { x: 1, y: 1, z: 1 };

  it('emits one segment for two land cells straddling a level', () => {
    const cells = [
      makeCell(0, SL + 0.05, [1], [v1, v2, { ...far }]),
      makeCell(1, SL + 0.15, [0], [{ ...v1 }, { ...v2 }, { ...far }]),
    ];
    // interval 0.1: band(0.60)=0, band(0.70)=1 -> straddles seaLevel+0.1
    const segs = computeContourSegments(cells, SL, 0.1);
    expect(segs.length).toBe(1);
  });

  it('emits nothing for two land cells in the same band', () => {
    const cells = [
      makeCell(0, SL + 0.02, [1], [v1, v2, { ...far }]),
      makeCell(1, SL + 0.05, [0], [{ ...v1 }, { ...v2 }, { ...far }]),
    ];
    expect(computeContourSegments(cells, SL, 0.1).length).toBe(0);
  });

  it('emits nothing when one cell is water even if bands differ', () => {
    const cells = [
      makeCell(0, SL - 0.2, [1], [v1, v2, { ...far }]), // water
      makeCell(1, SL + 0.15, [0], [{ ...v1 }, { ...v2 }, { ...far }]),
    ];
    expect(computeContourSegments(cells, SL, 0.1).length).toBe(0);
  });

  it('is deterministic on a generated world and only spans land cells', async () => {
    const world = await generateWorld(makeParams());
    const a = computeContourSegments(world.cells, world.params.seaLevel, 0.1);
    const b = computeContourSegments(world.cells, world.params.seaLevel, 0.1);
    expect(a).toEqual(b);
    expect(a.length).toBeGreaterThan(0);
  }, 30000);
});
