import { describe, it, expect } from 'vitest';
import { drawRulerTenant, RULER_RADIUS_RAISED, RULER_RADIUS_SMOOTH } from '../components/overlays/tenants';
import { makeFakeCtx, makeProjector } from './helpers/overlayCanvas';
import { Point } from '../types';

// A short great-circle arc of unit vectors in the XZ plane (equator).
function arc(n: number, turnFrac = 0.25): Point[] {
  const pts: Point[] = [];
  for (let i = 0; i < n; i++) {
    const a = (i / (n - 1)) * turnFrac * Math.PI * 2;
    pts.push({ x: Math.cos(a), y: 0, z: Math.sin(a) });
  }
  return pts;
}

describe('ruler overlay tenant', () => {
  it('draws the arc as one unbroken subpath plus two endpoint dots when all visible', () => {
    const ctx = makeFakeCtx();
    const points = arc(6);
    drawRulerTenant(ctx, points, makeProjector().project, false);

    const moves = ctx.ops.filter((o) => o.op === 'moveTo');
    const lines = ctx.ops.filter((o) => o.op === 'lineTo');
    const arcs = ctx.ops.filter((o) => o.op === 'arc');
    expect(moves).toHaveLength(1);
    expect(lines).toHaveLength(5); // 6 points -> 5 segments
    expect(arcs).toHaveLength(2);  // two endpoint dots
  });

  it('breaks the polyline at the horizon instead of drawing a chord across the globe', () => {
    // Cull the far half (z < 0). Arc runs front -> back -> front, so it must
    // start a new subpath when it re-emerges.
    const ctx = makeFakeCtx();
    // full ring so the arc dips behind and returns
    const pts: Point[] = [];
    const order = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1]; // start front, cross back, return
    for (const k of order) {
      const a = (k / 12) * Math.PI * 2;
      pts.push({ x: Math.cos(a), y: 0, z: Math.sin(a) });
    }
    const p = makeProjector((_x, _y, z) => z >= 0);
    drawRulerTenant(ctx, pts, p.project, false);
    const moves = ctx.ops.filter((o) => o.op === 'moveTo');
    expect(moves.length).toBeGreaterThanOrEqual(2); // no chord across the far side
  });

  it('projects every arc point at the fixed raised radius, never draped', () => {
    const p = makeProjector();
    drawRulerTenant(makeFakeCtx(), arc(6), p.project, false);
    // All radii (arc samples + endpoint dots) equal the fixed raised radius.
    for (const r of p.radii) expect(r).toBeCloseTo(RULER_RADIUS_RAISED, 6);
  });

  it('uses the smooth radius on the smooth globe', () => {
    const p = makeProjector();
    drawRulerTenant(makeFakeCtx(), arc(6), p.project, true);
    for (const r of p.radii) expect(r).toBeCloseTo(RULER_RADIUS_SMOOTH, 6);
  });

  it('draws nothing for a degenerate (single-point) ruler', () => {
    const ctx = makeFakeCtx();
    drawRulerTenant(ctx, arc(1), makeProjector().project, false);
    expect(ctx.ops).toHaveLength(0);
  });
});
