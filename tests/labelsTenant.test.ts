import { describe, it, expect } from 'vitest';
import { drawLabelsTenant, LABEL_LIFT } from '../components/overlays/tenants';
import { computeLabelAnchors, LabelAnchor } from '../utils/labelAnchors';
import { nearestCellBrute } from '../utils/nearestCell';
import { MapLabel } from '../utils/labels';
import { Cell, DEFAULT_LABEL_VISIBILITY } from '../types';
import { makeFakeCtx, makeProjector } from './helpers/overlayCanvas';

// `marker` is used as the default test kind: its globe-LOD gate
// (passesGlobeLod in tenants.ts) always passes regardless of camDist, so
// these tests don't have to thread a specific camera distance through every
// case — only the horizon/declutter/RAW-height/scope behaviour is under test.
// Within every kind's globe-LOD cutoff used below (town <= 2, capital <= 4),
// so no test result is silently decided by LOD filtering instead of the
// behaviour actually under test.
const CAM_DIST = 1.5;
const VIS = DEFAULT_LABEL_VISIBILITY; // markers: true

const label = (overrides: Partial<MapLabel> = {}): MapLabel => ({
  kind: 'marker',
  name: 'Test Label',
  position: { x: 1, y: 0, z: 0 },
  priority: 1,
  ...overrides,
});

const anchor = (height: number, overrides: Partial<MapLabel> = {}): LabelAnchor => ({
  label: label(overrides),
  height,
});

describe('labels overlay tenant', () => {
  it('projects on the unit sphere plus lift when the globe is smooth, regardless of height', () => {
    const anchors = [
      anchor(0.2, { position: { x: 1, y: 0, z: 0 } }),
      anchor(0.9, { name: 'Other', position: { x: 0, y: 1, z: 0 } }),
    ];
    const p = makeProjector();
    drawLabelsTenant(makeFakeCtx(), anchors, p.project, true, CAM_DIST, VIS);

    expect(p.radii).toHaveLength(2);
    for (const r of p.radii) expect(r).toBeCloseTo(1 + LABEL_LIFT, 6);
  });

  it('projects at 1 + h*0.05 + LABEL_LIFT on the raised globe; different heights give different radii', () => {
    const anchors = [
      anchor(0.6, { position: { x: 1, y: 0, z: 0 } }),
      anchor(0.2, { name: 'Other', position: { x: 0, y: 1, z: 0 } }),
    ];
    const p = makeProjector();
    drawLabelsTenant(makeFakeCtx(), anchors, p.project, false, CAM_DIST, VIS);

    expect(p.radii).toHaveLength(2);
    expect(p.radii[0]).toBeCloseTo(1 + 0.6 * 0.05 + LABEL_LIFT, 6);
    expect(p.radii[1]).toBeCloseTo(1 + 0.2 * 0.05 + LABEL_LIFT, 6);
    expect(p.radii[0]).not.toBeCloseTo(p.radii[1], 3);
  });

  it('uses the RAW height, never clamped to sea level, for an ocean-height anchor', () => {
    // Height 0 (ocean-floor) must NOT be lifted to any sea-level radius — the
    // S18 graticule bug, guarded here for labels (see drawGraticuleTenant).
    const anchors = [anchor(0, { kind: 'ocean', position: { x: 1, y: 0, z: 0 } })];
    const p = makeProjector();
    drawLabelsTenant(makeFakeCtx(), anchors, p.project, false, CAM_DIST, VIS);

    expect(p.radii).toHaveLength(1);
    expect(p.radii[0]).toBeCloseTo(1 + LABEL_LIFT, 6);
  });

  it('drops a label entirely when it falls past the horizon', () => {
    const anchors = [
      anchor(0.3, { name: 'Visible', position: { x: 1, y: 0, z: 0 } }),
      anchor(0.3, { name: 'Culled', position: { x: 0, y: 1, z: 0 } }),
    ];
    const ctx = makeFakeCtx();
    // Culls exactly the y-direction label (x === 0); the x-direction label
    // (x === 1*r > 0) survives. Verified against both labels individually
    // below, not just the aggregate count (S21 rivers-suite trap, plan §6).
    const p = makeProjector((x) => x > 0);
    drawLabelsTenant(ctx, anchors, p.project, false, CAM_DIST, VIS);

    const texts = ctx.ops.filter((o) => o.op === 'fillText').map((o) => o.text);
    expect(texts).toEqual(['Visible']);
  });

  it('declutters: two labels at the same screen point emit only the higher-priority fillText', () => {
    // Same position -> identical projected rect -> the greedy declutter in
    // drawMapLabels keeps whichever is FIRST in the array, so the higher
    // priority (lower `priority` number, per collectLabels' ascending sort)
    // must come first here.
    const anchors = [
      anchor(0.3, { name: 'Capital', kind: 'capital', priority: 1, position: { x: 1, y: 0, z: 0 } }),
      anchor(0.3, { name: 'Town', kind: 'town', priority: 3, position: { x: 1, y: 0, z: 0 } }),
    ];
    const ctx = makeFakeCtx();
    const p = makeProjector();
    drawLabelsTenant(ctx, anchors, p.project, false, CAM_DIST, VIS);

    const texts = ctx.ops.filter((o) => o.op === 'fillText').map((o) => o.text);
    expect(texts).toEqual(['Capital']);
  });

  it('never emits a faction label — faction names stay 3D (plan §1 scope guard)', () => {
    const anchors = [
      anchor(0.5, { name: 'Empire of Test', kind: 'faction', priority: 0, position: { x: 1, y: 0, z: 0 } }),
      anchor(0.5, { name: 'Some Marker', kind: 'marker', priority: 3, position: { x: 0, y: 1, z: 0 } }),
    ];
    const ctx = makeFakeCtx();
    const p = makeProjector();
    drawLabelsTenant(ctx, anchors, p.project, false, CAM_DIST, VIS);

    const texts = ctx.ops.filter((o) => o.op === 'fillText').map((o) => o.text);
    expect(texts).toEqual(['Some Marker']);
  });

  it('is a no-op with no anchors', () => {
    const ctx = makeFakeCtx();
    drawLabelsTenant(ctx, [], makeProjector().project, false, CAM_DIST, VIS);
    expect(ctx.ops).toHaveLength(0);
  });
});

describe('computeLabelAnchors', () => {
  // Four cells spread around the sphere so nearest-cell lookup is unambiguous,
  // each with a distinct height. `neighbors` links them in a ring so
  // nearestCellWalk can hill-climb between any pair.
  function makeCells(): Cell[] {
    const dirs: [number, number, number][] = [
      [1, 0, 0], [0, 1, 0], [-1, 0, 0], [0, -1, 0],
    ];
    const heights = [0.1, 0.4, 0.7, 0.95];
    return dirs.map((d, i) => ({
      id: i,
      center: { x: d[0], y: d[1], z: d[2] },
      height: heights[i],
      neighbors: [(i + 3) % 4, (i + 1) % 4],
    })) as unknown as Cell[];
  }

  it('returns one anchor per input label, in input order, with the nearest cell\'s height', () => {
    const cells = makeCells();
    const labels: MapLabel[] = [
      label({ name: 'A', position: { x: 0, y: 1, z: 0 } }), // nearest cell 1 (h 0.4)
      label({ name: 'B', position: { x: 1, y: 0, z: 0 } }), // nearest cell 0 (h 0.1)
      label({ name: 'C', position: { x: -1, y: 0, z: 0 } }), // nearest cell 2 (h 0.7)
      label({ name: 'D', position: { x: 0, y: -1, z: 0 } }), // nearest cell 3 (h 0.95)
    ];

    const anchors = computeLabelAnchors(labels, cells);

    expect(anchors).toHaveLength(4);
    anchors.forEach((a, i) => {
      expect(a.label).toBe(labels[i]); // input order preserved, same object
      const expectedId = nearestCellBrute(cells, labels[i].position.x, labels[i].position.y, labels[i].position.z);
      expect(a.height).toBe(cells[expectedId].height);
    });
    expect(anchors.map((a) => a.height)).toEqual([0.4, 0.1, 0.7, 0.95]);
  });

  it('returns an empty array for no labels', () => {
    expect(computeLabelAnchors([], makeCells())).toHaveLength(0);
  });

  it('handles an empty cell set without throwing (height 0)', () => {
    const labels = [label({ position: { x: 1, y: 0, z: 0 } })];
    const anchors = computeLabelAnchors(labels, []);
    expect(anchors).toEqual([{ label: labels[0], height: 0 }]);
  });
});
