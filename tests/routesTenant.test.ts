import { describe, it, expect } from 'vitest';
import { drawRoutesTenant, ROUTE_LIFT } from '../components/overlays/tenants';
import { ProjectedCells } from '../utils/screenProject';
import { WorldData, RouteData } from '../types';
import { makeFakeCtx, makeProjector, makeRingWorld } from './helpers/overlayCanvas';

const emptyProj: ProjectedCells = {
  x: new Float32Array(0), y: new Float32Array(0), visible: new Uint8Array(0), n: 0,
};

// A route along the equatorial ring of makeRingWorld, cells `ids` in order.
function ringRoute(world: WorldData, ids: number[], kind: RouteData['kind']): RouteData {
  return {
    path: ids.map((i) => world.cells[i].center),
    kind,
    fromCellId: ids[0],
    toCellId: ids[ids.length - 1],
  };
}

function withRoutes(world: WorldData, routes: RouteData[]): WorldData {
  return { ...world, routes } as WorldData;
}

describe('routes overlay tenant', () => {
  it('draws a road as one unbroken subpath when every point is visible', () => {
    const base = makeRingWorld(16, [1.0], 0.5);
    const world = withRoutes(base, [ringRoute(base, [0, 1, 2, 3], 'road')]);
    const ctx = makeFakeCtx();
    drawRoutesTenant(ctx, emptyProj, world, makeProjector().project, true);

    const moves = ctx.ops.filter((o) => o.op === 'moveTo');
    const lines = ctx.ops.filter((o) => o.op === 'lineTo');
    expect(moves).toHaveLength(1);
    expect(lines).toHaveLength(3);
  });

  it('breaks the polyline at the horizon instead of drawing a chord across the globe', () => {
    // Cull the far half of the ring (z < 0). A route that runs through the far
    // side must start a NEW subpath when it re-emerges, never lineTo across.
    // Ring cell i sits at angle i/16 turns, so z >= 0 for ids 0..8 and z < 0 for
    // 9..15. This path runs visible (4..8) -> culled (9..15) -> visible (0..3),
    // so the far-side crossing must produce a second subpath.
    const base = makeRingWorld(16, [1.0], 0.5);
    const ids = [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3];
    const world = withRoutes(base, [ringRoute(base, ids, 'road')]);
    const ctx = makeFakeCtx();
    const p = makeProjector((_x, _y, z) => z >= 0);
    drawRoutesTenant(ctx, emptyProj, world, p.project, true);

    const moves = ctx.ops.filter((o) => o.op === 'moveTo');
    expect(moves).toHaveLength(2); // one per visible run, no chord between them

    // No lineTo may follow a culled point directly: every lineTo must be
    // preceded by a moveTo or another lineTo in the same subpath.
    const drawOps = ctx.ops.filter((o) => o.op === 'moveTo' || o.op === 'lineTo');
    expect(drawOps[0].op).toBe('moveTo');
  });

  it('dashes sea routes and resets the shared line dash before returning', () => {
    const base = makeRingWorld(16, [0], 0.5);
    const world = withRoutes(base, [ringRoute(base, [0, 1, 2], 'searoute')]);
    const ctx = makeFakeCtx();
    drawRoutesTenant(ctx, emptyProj, world, makeProjector().project, true);

    // A dash pattern was applied at some point...
    expect(ctx.dashHistory.some((d) => d.length > 0)).toBe(true);
    // ...and the context was handed back to the next tenant undashed.
    expect(ctx.lineDash).toEqual([]);
  });

  it('skips degenerate routes with fewer than two points', () => {
    const base = makeRingWorld(16, [1.0], 0.5);
    const world = withRoutes(base, [ringRoute(base, [0], 'road')]);
    const ctx = makeFakeCtx();
    drawRoutesTenant(ctx, emptyProj, world, makeProjector().project, true);

    expect(ctx.ops.filter((o) => o.op === 'moveTo')).toHaveLength(0);
  });

  it('is a no-op when the world has no routes', () => {
    const world = makeRingWorld(16, [1.0], 0.5);
    const ctx = makeFakeCtx();
    drawRoutesTenant(ctx, emptyProj, world, makeProjector().project, false);

    expect(ctx.ops).toHaveLength(0);
  });

  it('projects on the unit sphere plus lift when the globe is smooth', () => {
    const base = makeRingWorld(16, [1.0], 0.5);
    const world = withRoutes(base, [ringRoute(base, [0, 1, 2], 'road')]);
    const p = makeProjector();
    drawRoutesTenant(makeFakeCtx(), emptyProj, world, p.project, true);

    expect(p.radii.length).toBe(3);
    for (const r of p.radii) expect(r).toBeCloseTo(1 + ROUTE_LIFT, 6);
  });

  it('projects at the sea-level radius plus lift on the raised globe (phase 1, flat)', () => {
    const base = makeRingWorld(16, [1.0], 0.5);
    const world = withRoutes(base, [ringRoute(base, [0, 1, 2], 'road')]);
    const p = makeProjector();
    drawRoutesTenant(makeFakeCtx(), emptyProj, world, p.project, false);

    for (const r of p.radii) expect(r).toBeCloseTo(1 + 0.5 * 0.05 + ROUTE_LIFT, 6);
  });
});
