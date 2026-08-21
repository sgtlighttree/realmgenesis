import { describe, it, expect } from 'vitest';
import { computeBorderSegments } from '../utils/borders';
import { drawBordersTenant, BORDER_LIFT } from '../components/overlays/tenants';
import { Cell, Point } from '../types';
import { makeFakeCtx, makeProjector } from './helpers/overlayCanvas';

const v = (x: number, y: number, z: number): Point => ({ x, y, z });

// Shared corners of the edge between cell 0 and cell 1.
const SHARED_A = v(0, 1, 0);
const SHARED_B = v(0, -1, 0);

// Two adjacent cells that share exactly one edge (SHARED_A..SHARED_B), plus one
// private corner each. `makeRingWorld` cannot be reused here: its cells carry no
// `vertices`, and the shared-edge scan is the whole point of the extraction.
function makePair(hA: number, hB: number, rA?: number, rB?: number): Cell[] {
  return [
    {
      id: 0, height: hA, regionId: rA, neighbors: [1],
      center: v(1, 0, 0), vertices: [SHARED_A, SHARED_B, v(1, 0, 0)],
    },
    {
      id: 1, height: hB, regionId: rB, neighbors: [0],
      center: v(-1, 0, 0), vertices: [SHARED_A, SHARED_B, v(-1, 0, 0)],
    },
  ] as unknown as Cell[];
}

describe('computeBorderSegments', () => {
  it('emits one segment for an edge between two different regions', () => {
    const segs = computeBorderSegments(makePair(0.6, 0.7, 1, 2));
    expect(segs).toHaveLength(1);
    expect(segs[0].a).toEqual(SHARED_A);
    expect(segs[0].b).toEqual(SHARED_B);
  });

  it('emits nothing when both cells belong to the same region', () => {
    expect(computeBorderSegments(makePair(0.6, 0.7, 1, 1))).toHaveLength(0);
  });

  it('emits nothing when neither cell is claimed', () => {
    // `undefined !== undefined` is false, so unclaimed ocean never borders itself.
    expect(computeBorderSegments(makePair(0.2, 0.2))).toHaveLength(0);
  });

  it('emits a segment between a claimed cell and unclaimed territory', () => {
    expect(computeBorderSegments(makePair(0.6, 0.2, 1, undefined))).toHaveLength(1);
  });

  it('processes each neighbor pair once, not twice', () => {
    // Both cells list the other as a neighbor; the id guard must dedupe.
    expect(computeBorderSegments(makePair(0.6, 0.7, 1, 2))).toHaveLength(1);
  });

  it('carries the height of the TALLER of the two cells', () => {
    expect(computeBorderSegments(makePair(0.2, 0.9, 1, 2))[0].height).toBe(0.9);
    expect(computeBorderSegments(makePair(0.9, 0.2, 1, 2))[0].height).toBe(0.9);
  });
});

describe('borders overlay tenant', () => {
  it('projects on the unit sphere plus lift when the globe is smooth', () => {
    const segs = computeBorderSegments(makePair(1.0, 0.4, 1, 2));
    const p = makeProjector();
    drawBordersTenant(makeFakeCtx(), segs, p.project, true);

    expect(p.radii).toHaveLength(2); // both ends of the one segment
    for (const r of p.radii) expect(r).toBeCloseTo(1 + BORDER_LIFT, 6);
  });

  it('rides the taller cell radius on the raised globe', () => {
    // The mesh draws a cell boundary as a vertical step; the border must crown
    // that step, so the radius comes from the peak, not the valley or a mean.
    const segs = computeBorderSegments(makePair(0, 1.0, 1, 2));
    const p = makeProjector();
    drawBordersTenant(makeFakeCtx(), segs, p.project, false);

    expect(p.radii).toHaveLength(2);
    for (const r of p.radii) expect(r).toBeCloseTo(1.05 + BORDER_LIFT, 6);
  });

  it('uses the RAW height, never clamped to sea level', () => {
    // The S18 graticule bug, guarded for borders: ocean cells render at their
    // true seafloor radius, so clamping to seaLevel would float a coastal border.
    const segs = computeBorderSegments(makePair(0, 0, 1, undefined));
    const p = makeProjector();
    drawBordersTenant(makeFakeCtx(), segs, p.project, false);

    for (const r of p.radii) expect(r).toBeCloseTo(1 + BORDER_LIFT, 6);
  });

  it('drops a segment entirely when either end is past the horizon', () => {
    // One cell edge is far too short to be worth clipping against the limb, so
    // a partly hidden segment is skipped rather than drawn as a stub.
    const segs = computeBorderSegments(makePair(0.6, 0.7, 1, 2));
    const ctx = makeFakeCtx();
    const p = makeProjector((_x, y) => y > 0); // SHARED_B (y = -1) is culled
    drawBordersTenant(ctx, segs, p.project, true);

    expect(ctx.ops.filter((o) => o.op === 'moveTo')).toHaveLength(0);
    expect(ctx.ops.filter((o) => o.op === 'lineTo')).toHaveLength(0);
  });

  it('is a no-op with no segments', () => {
    const ctx = makeFakeCtx();
    drawBordersTenant(ctx, [], makeProjector().project, false);
    expect(ctx.ops).toHaveLength(0);
  });
});
