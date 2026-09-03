import { describe, it, expect, beforeAll } from 'vitest';
import { chainSegments, computeCoastlineSegments } from '../utils/boundaries';
import { Point3 } from '../utils/geo';
import { generateWorld } from '../utils/worldGen';
import { makeParams } from './helpers';

const seg = (a: Point3, b: Point3): [Point3, Point3] => [a, b];

describe('chainSegments', () => {
  it('joins a broken-order open chain into one polyline', () => {
    // Endpoints given out of order and reversed.
    const segs: Array<[Point3, Point3]> = [
      seg([2, 0, 0], [3, 0, 0]),
      seg([1, 0, 0], [2, 0, 0]),
      seg([0, 0, 0], [1, 0, 0]),
    ];
    const chains = chainSegments(segs);
    expect(chains).toHaveLength(1);
    expect(chains[0]).toHaveLength(4); // 4 points, open
    // First and last differ (open chain).
    expect(chains[0][0]).not.toEqual(chains[0][chains[0].length - 1]);
  });

  it('closes a ring (last point equals first)', () => {
    const segs: Array<[Point3, Point3]> = [
      seg([0, 0, 0], [1, 0, 0]),
      seg([1, 0, 0], [1, 1, 0]),
      seg([1, 1, 0], [0, 0, 0]),
    ];
    const chains = chainSegments(segs);
    expect(chains).toHaveLength(1);
    expect(chains[0][0]).toEqual(chains[0][chains[0].length - 1]);
  });

  it('conserves edge count across all chains', () => {
    const segs: Array<[Point3, Point3]> = [
      seg([0, 0, 0], [1, 0, 0]),
      seg([1, 0, 0], [2, 0, 0]),
      seg([5, 5, 5], [6, 5, 5]), // disjoint second chain
    ];
    const chains = chainSegments(segs);
    const edges = chains.reduce((s, c) => s + (c.length - 1), 0);
    expect(edges).toBe(segs.length);
    expect(chains).toHaveLength(2);
  });

  it('handles a T-junction deterministically without dropping edges', () => {
    // Three edges share the vertex [1,0,0].
    const segs: Array<[Point3, Point3]> = [
      seg([0, 0, 0], [1, 0, 0]),
      seg([1, 0, 0], [2, 0, 0]),
      seg([1, 0, 0], [1, 1, 0]),
    ];
    const chains = chainSegments(segs);
    const edges = chains.reduce((s, c) => s + (c.length - 1), 0);
    expect(edges).toBe(3); // no edge lost at the branch
  });
});

describe('coastline chains close (hard guard)', () => {
  let unclosed = 0;
  let total = 0;
  let segmentCount = 0;
  let edgeCount = 0;
  beforeAll(async () => {
    const world = await generateWorld(makeParams({ points: 3000, seed: 'chain-test' }));
    const segments = computeCoastlineSegments(world);
    segmentCount = segments.length;
    const chains = chainSegments(segments);
    total = chains.length;
    edgeCount = chains.reduce((s, c) => s + (c.length - 1), 0);
    for (const c of chains) {
      if (c.length < 2) { unclosed++; continue; }
      const a = c[0]; const b = c[c.length - 1];
      const same = Math.abs(a[0] - b[0]) < 1e-4 && Math.abs(a[1] - b[1]) < 1e-4 && Math.abs(a[2] - b[2]) < 1e-4;
      if (!same) unclosed++;
    }
  });
  it('has at least one coastline chain', () => { expect(total).toBeGreaterThan(0); });
  it('conserves edge count at scale (no edge dropped or double-counted by the welder)', () => {
    expect(edgeCount).toBe(segmentCount);
  });
  it('every coastline chain closes into a ring', () => { expect(unclosed).toBe(0); });
});
