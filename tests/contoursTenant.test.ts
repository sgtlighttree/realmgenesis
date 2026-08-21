import { describe, it, expect } from 'vitest';
import { drawContoursTenant, CONTOUR_LIFT } from '../components/overlays/tenants';
import {
  computeContourSegments, ContourSegment, CONTOUR_INDEX_EVERY,
  contourInterval, CONTOUR_TARGET_LEVELS, contourLabel,
} from '../utils/shading';
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

const seg = (height: number, index: boolean, elevation = SL + 0.1): ContourSegment => ({
  a: [0, 1, 0], b: [0, 0, 1], height, index, elevation,
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
    // Asserts the RULE across levels rather than one constant's alternation:
    // the previous form encoded CONTOUR_INDEX_EVERY === 2 and broke when the
    // adaptive interval let the standard 5th be restored.
    const interval = 0.02;
    for (let level = 1; level <= 2 * CONTOUR_INDEX_EVERY; level++) {
      // Two cells straddling exactly this level: band level-1 and band level.
      const below = SL + (level - 1) * interval + interval / 2;
      const above = SL + level * interval + interval / 2;
      const segs = computeContourSegments(pair(below, above), SL, interval);
      expect(segs).toHaveLength(1);
      expect(segs[0].index).toBe(level % CONTOUR_INDEX_EVERY === 0);
    }
  });

  it('carries the elevation of the contour LINE, not either cell height', () => {
    const interval = 0.02;
    const segs = computeContourSegments(pair(SL + 0.01, SL + 0.03), SL, interval);
    // bands 0 and 1 straddle level 1 -> the line sits at seaLevel + 1*interval.
    expect(segs[0].elevation).toBeCloseTo(SL + interval, 6);
    expect(segs[0].elevation).not.toBeCloseTo(segs[0].height, 6);
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

describe('adaptive contour interval', () => {
  const world = (heights: number[]): Cell[] =>
    heights.map((h, i) => makeCell(i, h, [], [v1, v2, { ...far }]));

  it('yields roughly CONTOUR_TARGET_LEVELS bands regardless of relief', () => {
    // A flat world and an alpine one must both get a usable number of lines.
    // The shipped fixed 0.1 gave four possible levels on a default world
    // (seaLevel 0.55), which is what made contours read as blobby outlines.
    for (const max of [0.6, 0.75, 1.0]) {
      const interval = contourInterval(world([SL, max]), SL);
      const levels = (max - SL) / interval;
      expect(levels).toBeGreaterThan(CONTOUR_TARGET_LEVELS / 3);
      expect(levels).toBeLessThan(CONTOUR_TARGET_LEVELS * 3);
    }
  });

  it('beats the old fixed 0.1 on a default-relief world', () => {
    const cells = world([SL, 0.75]);
    const adaptive = (0.75 - SL) / contourInterval(cells, SL);
    const fixed = (0.75 - SL) / 0.1;
    expect(fixed).toBeLessThan(3);           // the bug: two lines of terrain
    expect(adaptive).toBeGreaterThan(fixed); // the fix
  });

  it('returns 0 when there is no land above sea level', () => {
    expect(contourInterval(world([0.1, 0.3, SL]), SL)).toBe(0);
    // and computeContourSegments treats 0 as "no contours" rather than dividing
    expect(computeContourSegments(world([0.1, 0.3]), SL, 0)).toHaveLength(0);
  });
});

describe('contour elevation labels', () => {
  it('reads out the line elevation as a percentage', () => {
    expect(contourLabel(0.68)).toBe('68%');
    expect(contourLabel(SL)).toBe('55%');
  });

  it('thins labels so they do not stack, and caps the total', () => {
    // 200 index segments all projecting near the same spot: the min-gap rule
    // must collapse them to a handful, not draw 200 overlapping readouts.
    const segs: ContourSegment[] = [];
    for (let i = 0; i < 200; i++) segs.push(seg(0.8, true, SL + 0.1));
    const ctx = makeFakeCtx();
    drawContoursTenant(ctx, segs, makeProjector().project, true);

    const drawn = ctx.ops.filter(o => o.op === 'fillText');
    expect(drawn.length).toBeGreaterThan(0);
    expect(drawn.length).toBeLessThan(5); // all coincident -> essentially one
    expect(drawn[0].text).toBe('65%');
  });

  it('labels only index contours', () => {
    const ctx = makeFakeCtx();
    drawContoursTenant(ctx, [seg(0.8, false, SL + 0.1)], makeProjector().project, true);
    expect(ctx.ops.filter(o => o.op === 'fillText')).toHaveLength(0);
  });
});
