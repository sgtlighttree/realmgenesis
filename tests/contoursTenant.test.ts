import { describe, it, expect } from 'vitest';
import { drawContoursTenant, CONTOUR_LIFT } from '../components/overlays/tenants';
import { computeContourSegments, ContourSegment, CONTOUR_INDEX_EVERY } from '../utils/shading';
import { Cell, BiomeType } from '../types';
import { makeFakeCtx, makeProjector } from './helpers/overlayCanvas';

const SL = 0.55;
const v1 = { x: 0, y: 1, z: 0 };
const v2 = { x: 0, y: 0, z: 1 };
const far = { x: 1, y: 1, z: 1 };

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

// Two neighboring land cells sharing edge v1-v2, at the given heights.
const pair = (hA: number, hB: number): Cell[] => [
  makeCell(0, hA, [1], [v1, v2, { ...far }]),
  makeCell(1, hB, [0], [{ ...v1 }, { ...v2 }, { ...far }]),
];

const seg = (height: number, index: boolean): ContourSegment => ({
  a: [0, 1, 0], b: [0, 0, 1], height, index,
});

describe('contour segment metadata', () => {
  it('carries the height of the TALLER of the two cells', () => {
    // The globe renders a cell boundary as a vertical step; riding the taller
    // cell's radius crowns that step instead of cutting into it.
    const segs = computeContourSegments(pair(SL + 0.05, SL + 0.15), SL, 0.1);
    expect(segs).toHaveLength(1);
    expect(segs[0].height).toBeCloseTo(SL + 0.15, 6);
  });

  it('carries the taller height regardless of cell order', () => {
    const segs = computeContourSegments(pair(SL + 0.15, SL + 0.05), SL, 0.1);
    expect(segs[0].height).toBeCloseTo(SL + 0.15, 6);
  });

  it('flags an index contour on every CONTOUR_INDEX_EVERY-th level', () => {
    // bands 0 and 1 straddle level 1; bands 1 and 2 straddle level 2.
    const lvl1 = computeContourSegments(pair(SL + 0.05, SL + 0.15), SL, 0.1);
    const lvl2 = computeContourSegments(pair(SL + 0.15, SL + 0.25), SL, 0.1);
    expect(lvl1[0].index).toBe(1 % CONTOUR_INDEX_EVERY === 0);
    expect(lvl2[0].index).toBe(2 % CONTOUR_INDEX_EVERY === 0);
    // With the shipped setting the two must differ — that IS the alternation.
    expect(lvl1[0].index).not.toBe(lvl2[0].index);
  });
});

describe('contours overlay tenant', () => {
  it('projects each segment at its own cell radius plus the lift', () => {
    const segs = [seg(0, false), seg(1.0, false)];
    const p = makeProjector();
    drawContoursTenant(makeFakeCtx(), segs, p.project, false);

    // two endpoints per segment
    expect(p.radii).toHaveLength(4);
    expect(p.radii[0]).toBeCloseTo(1 + CONTOUR_LIFT, 6);
    expect(p.radii[1]).toBeCloseTo(1 + CONTOUR_LIFT, 6);
    expect(p.radii[2]).toBeCloseTo(1.05 + CONTOUR_LIFT, 6);
    expect(p.radii[3]).toBeCloseTo(1.05 + CONTOUR_LIFT, 6);
  });

  it('collapses to the unit sphere plus lift when the globe is smooth', () => {
    const p = makeProjector();
    drawContoursTenant(makeFakeCtx(), [seg(1.0, false)], p.project, true);

    for (const r of p.radii) expect(r).toBeCloseTo(1 + CONTOUR_LIFT, 6);
  });

  it('drops a segment when either endpoint is past the horizon', () => {
    // b = (0,0,1) is culled; the segment must not be drawn half-way.
    const p = makeProjector((_x, _y, z) => z <= 0);
    const ctx = makeFakeCtx();
    drawContoursTenant(ctx, [seg(0.8, false)], p.project, false);

    expect(ctx.ops.filter((o) => o.op === 'moveTo')).toHaveLength(0);
    expect(ctx.ops.filter((o) => o.op === 'lineTo')).toHaveLength(0);
  });

  it('draws index contours in a second pass, over the intermediates', () => {
    const ctx = makeFakeCtx();
    drawContoursTenant(ctx, [seg(0.8, true), seg(0.8, false)], makeProjector().project, false);

    // Two strokes: intermediates first, index second.
    const strokes = ctx.ops.filter((o) => o.op === 'stroke');
    expect(strokes).toHaveLength(2);
    const firstStroke = ctx.ops.findIndex((o) => o.op === 'stroke');
    const moves = ctx.ops.map((o, i) => ({ ...o, i })).filter((o) => o.op === 'moveTo');
    expect(moves).toHaveLength(2);
    expect(moves[0].i).toBeLessThan(firstStroke); // intermediate in pass 1
    expect(moves[1].i).toBeGreaterThan(firstStroke); // index in pass 2
  });

  it('is a no-op with no segments', () => {
    const ctx = makeFakeCtx();
    drawContoursTenant(ctx, [], makeProjector().project, false);
    expect(ctx.ops).toHaveLength(0);
  });
});
