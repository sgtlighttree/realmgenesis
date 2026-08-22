import { describe, it, expect } from 'vitest';
import * as THREE from 'three';

import { drawDymaxionTenant, DYMAXION_CAGE_RADIUS } from '../components/overlays/tenants';
import { cageEdges, ICOSA_VERTS, ICOSA_EDGES } from '../utils/dymaxionCage';
import { makeFakeCtx, makeProjector } from './helpers/overlayCanvas';
import { Point, DymaxionSettings } from '../types';

const SAMPLES = 16; // must match DYMAXION_SAMPLES in the tenant

function settings(lon: number, lat: number, roll: number): DymaxionSettings {
  return { lon, lat, roll, showOverlay: true, layout: 'blender', mode: 'overlay' };
}

describe('dymaxion cage geometry', () => {
  it('has 12 unit vertices and 30 edges', () => {
    expect(ICOSA_VERTS).toHaveLength(12);
    expect(ICOSA_EDGES).toHaveLength(30);
    for (const v of ICOSA_VERTS) {
      expect(Math.hypot(v.x, v.y, v.z)).toBeCloseTo(1, 6);
    }
  });

  it('cageEdges returns 30 unit-length edge pairs', () => {
    const edges = cageEdges(settings(37, -12, 5));
    expect(edges).toHaveLength(30);
    for (const [a, b] of edges) {
      expect(Math.hypot(a.x, a.y, a.z)).toBeCloseTo(1, 6);
      expect(Math.hypot(b.x, b.y, b.z)).toBeCloseTo(1, 6);
    }
  });

  it('identity settings leave the base orientation untouched', () => {
    const edges = cageEdges(settings(0, 0, 0));
    const [a, b] = edges[0]; // ICOSA_EDGES[0] = [1, 2]
    expect(a.x).toBeCloseTo(ICOSA_VERTS[1].x, 6);
    expect(a.z).toBeCloseTo(ICOSA_VERTS[1].z, 6);
    expect(b.y).toBeCloseTo(ICOSA_VERTS[2].y, 6);
  });

  it('rotates by the same YXZ euler the old 3D cage used', () => {
    const s = settings(37, -12, 5);
    const euler = new THREE.Euler(
      THREE.MathUtils.degToRad(s.lat),
      -THREE.MathUtils.degToRad(s.lon),
      THREE.MathUtils.degToRad(s.roll),
      'YXZ',
    );
    const expected = new THREE.Vector3(ICOSA_VERTS[0].x, ICOSA_VERTS[0].y, ICOSA_VERTS[0].z).applyEuler(euler);
    // ICOSA_EDGES[3] = [0, 2], so vertex 0 is the first endpoint.
    const got = cageEdges(s)[3][0];
    expect(got.x).toBeCloseTo(expected.x, 6);
    expect(got.y).toBeCloseTo(expected.y, 6);
    expect(got.z).toBeCloseTo(expected.z, 6);
  });
});

describe('dymaxion cage overlay tenant', () => {
  const oneEdge: [Point, Point][] = [[{ x: 1, y: 0, z: 0 }, { x: 0, y: 1, z: 0 }]];

  it('draws a fully visible edge as one sampled subpath', () => {
    const ctx = makeFakeCtx();
    drawDymaxionTenant(ctx, oneEdge, makeProjector().project, false);
    const moves = ctx.ops.filter((o) => o.op === 'moveTo');
    const lines = ctx.ops.filter((o) => o.op === 'lineTo');
    expect(moves).toHaveLength(1);
    expect(lines).toHaveLength(SAMPLES); // SAMPLES+1 points -> SAMPLES segments
  });

  it('draws nothing for an empty cage', () => {
    const ctx = makeFakeCtx();
    drawDymaxionTenant(ctx, [], makeProjector().project, false);
    expect(ctx.ops.filter((o) => o.op === 'moveTo')).toHaveLength(0);
  });

  it('projects every sample at the fixed cage radius, never draped, in both globe modes', () => {
    const raised = makeProjector();
    drawDymaxionTenant(makeFakeCtx(), oneEdge, raised.project, false);
    const smooth = makeProjector();
    drawDymaxionTenant(makeFakeCtx(), oneEdge, smooth.project, true);
    // Fixed radius: raised and smooth must be byte-identical (not draped).
    expect(smooth.radii).toEqual(raised.radii);
    // Endpoints sit at exactly the cage radius; the chord dips below it.
    expect(Math.max(...raised.radii)).toBeCloseTo(DYMAXION_CAGE_RADIUS, 6);
    for (const r of raised.radii) expect(r).toBeLessThanOrEqual(DYMAXION_CAGE_RADIUS + 1e-9);
  });

  it('breaks at the horizon instead of drawing a chord across the far side', () => {
    // Cull the near band so the edge is visible, then hidden, then visible:
    // a straight chord from +z to -z, kept only where |z| is large.
    const ctx = makeFakeCtx();
    const edge: [Point, Point][] = [[{ x: 0, y: 0, z: 1 }, { x: 0, y: 0, z: -1 }]];
    const p = makeProjector((_x, _y, z) => Math.abs(z) > 0.5 * DYMAXION_CAGE_RADIUS);
    drawDymaxionTenant(ctx, edge, p.project, false);
    const moves = ctx.ops.filter((o) => o.op === 'moveTo');
    expect(moves.length).toBeGreaterThanOrEqual(2); // two visible runs, no chord across
  });

  it('drops an edge whose samples are all behind the horizon', () => {
    const ctx = makeFakeCtx();
    drawDymaxionTenant(ctx, oneEdge, makeProjector(() => false).project, false);
    expect(ctx.ops.filter((o) => o.op === 'moveTo')).toHaveLength(0);
    expect(ctx.ops.filter((o) => o.op === 'lineTo')).toHaveLength(0);
  });
});
