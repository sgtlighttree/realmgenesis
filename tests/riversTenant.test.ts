import { describe, it, expect } from 'vitest';
import { drawRiversTenant, RIVER_LIFT } from '../components/overlays/tenants';
import { computeRiverPolylines } from '../utils/riverPaths';
import { Point } from '../types';
import { makeFakeCtx, makeProjector } from './helpers/overlayCanvas';

// Builds one river point the way generation does. `getRenderPoint`
// (utils/worldGen.ts ~L310) pre-scales every river point to
// r = 1 + cell.height·0.05 + 0.005, so a river point is NOT a unit direction —
// which is the whole trap drawRiversTenant exists to respect.
//
// The `x, y, z` argument here is a DIRECTION and is normalized before scaling,
// so the returned point's |P| really is that baked radius whatever the caller
// passes. Without that the helper silently lies for any non-unit direction, and
// a radius assertion written against it would be checking the caller's own
// arithmetic rather than the tenant's.
function bakedPoint(height: number, x = 1, y = 0, z = 0): Point {
  const r = 1 + height * 0.05 + 0.005;
  const len = Math.hypot(x, y, z) || 1;
  const k = r / len;
  return { x: x * k, y: y * k, z: z * k };
}

describe('rivers overlay tenant', () => {
  it('projects on the unit sphere plus lift when the globe is smooth', () => {
    const path = [bakedPoint(0.2, 1, 0, 0), bakedPoint(0.8, 0, 1, 0)];
    const p = makeProjector();
    drawRiversTenant(makeFakeCtx(), [path], p.project, true);

    expect(p.radii).toHaveLength(2);
    for (const r of p.radii) expect(r).toBeCloseTo(1 + RIVER_LIFT, 6);
  });

  it('projects at each point\'s own baked radius on the raised globe, radii differ', () => {
    const path = [bakedPoint(0.1, 1, 0, 0), bakedPoint(0.9, 0, 1, 0)];
    const p = makeProjector();
    drawRiversTenant(makeFakeCtx(), [path], p.project, false);

    expect(p.radii).toHaveLength(2);
    expect(p.radii[0]).toBeCloseTo(1 + 0.1 * 0.05 + 0.005, 6);
    expect(p.radii[1]).toBeCloseTo(1 + 0.9 * 0.05 + 0.005, 6);
    expect(p.radii[0]).not.toBeCloseTo(p.radii[1], 3);
  });

  it('does not double-scale a baked point on the raised globe (the TRAP guard)', () => {
    const baked = bakedPoint(0.5, 1, 0, 0); // r = 1 + 0.5*0.05 + 0.005
    const p = makeProjector();
    drawRiversTenant(makeFakeCtx(), [[baked, bakedPoint(0.5, 0, 1, 0)]], p.project, false);

    const expectedR = 1 + 0.5 * 0.05 + 0.005;
    expect(p.radii[0]).toBeCloseTo(expectedR, 6);
  });

  it('breaks the polyline at the horizon instead of drawing a chord across the globe', () => {
    // Visibility keyed on the SIGN of z, the same discriminator the routes
    // suite uses. A threshold on a scaled coordinate is too easy to get wrong:
    // an earlier draft of this test culled only one of the two points it
    // claimed to cull, and still passed for the wrong reason.
    const path = [
      bakedPoint(0.2, 1, 0, 0.5),   // z > 0, visible
      bakedPoint(0.2, 0, 1, -0.5),  // z < 0, culled
      bakedPoint(0.2, -1, 0, -0.5), // z < 0, culled
      bakedPoint(0.2, 0, -1, 0.5),  // z > 0, visible
    ];
    const ctx = makeFakeCtx();
    const p = makeProjector((_x, _y, z) => z >= 0);
    drawRiversTenant(ctx, [path], p.project, false);

    const moves = ctx.ops.filter((o) => o.op === 'moveTo');
    expect(moves).toHaveLength(2); // one per visible run, no chord between them

    const drawOps = ctx.ops.filter((o) => o.op === 'moveTo' || o.op === 'lineTo');
    expect(drawOps[0].op).toBe('moveTo');
  });

  it('is a no-op on empty input', () => {
    const ctx = makeFakeCtx();
    drawRiversTenant(ctx, [], makeProjector().project, false);
    expect(ctx.ops).toHaveLength(0);
  });

  it('is a no-op on a single-point path', () => {
    const ctx = makeFakeCtx();
    drawRiversTenant(ctx, [[bakedPoint(0.4)]], makeProjector().project, false);
    expect(ctx.ops.filter((o) => o.op === 'moveTo' || o.op === 'lineTo')).toHaveLength(0);
  });
});

describe('computeRiverPolylines', () => {
  it('drops paths shorter than two points', () => {
    const rivers: Point[][] = [[bakedPoint(0.3)], []];
    expect(computeRiverPolylines(rivers)).toHaveLength(0);
  });

  it('preserves the baked radius of its inputs (does not normalize)', () => {
    // Both endpoints share a baked radius far from 1 (height 0.9 -> r ~1.05).
    // If computeRiverPolylines normalized control points to unit before
    // smoothing (as `RiverLines` did for smoothGlobe), every output point would
    // collapse toward r=1; it must not, since the offset now moves to the
    // tenant, applied only at draw time.
    const a = bakedPoint(0.9, 1, 0, 0);
    const b = bakedPoint(0.9, 0.999, 0.0447, 0); // ~2.6 deg away, same radius
    const [smoothed] = computeRiverPolylines([[a, b]]);

    expect(smoothed.length).toBeGreaterThan(0);
    const expectedR = 1 + 0.9 * 0.05 + 0.005; // 1.05
    for (const p of smoothed) {
      const r = Math.hypot(p.x, p.y, p.z);
      expect(r).toBeCloseTo(expectedR, 2);
    }
  });
});
