import { describe, it, expect } from 'vitest';
import * as d3 from 'd3';
import { greatCircleDistanceKm, sampleGreatCircleArc, computeScaleBar, computeInterruptedScaleBar, niceScaleBarLength } from '../utils/measure';
import { Point } from '../types';

const R = 6371;
const north: Point = { x: 0, y: 1, z: 0 };
const south: Point = { x: 0, y: -1, z: 0 };
const equatorA: Point = { x: 1, y: 0, z: 0 };
const equatorB: Point = { x: 0, y: 0, z: 1 };

describe('greatCircleDistanceKm', () => {
  it('is 0 for identical points', () => {
    expect(greatCircleDistanceKm(equatorA, equatorA, R)).toBeCloseTo(0, 6);
  });

  it('is pi * R for antipodal points', () => {
    expect(greatCircleDistanceKm(north, south, R)).toBeCloseTo(Math.PI * R, 6);
  });

  it('is (pi / 2) * R for points 90 degrees apart', () => {
    expect(greatCircleDistanceKm(equatorA, equatorB, R)).toBeCloseTo((Math.PI / 2) * R, 6);
    expect(greatCircleDistanceKm(equatorA, north, R)).toBeCloseTo((Math.PI / 2) * R, 6);
  });
});

describe('sampleGreatCircleArc', () => {
  it('starts and ends at the input points', () => {
    const arc = sampleGreatCircleArc(equatorA, equatorB, 16);
    expect(arc[0].x).toBeCloseTo(equatorA.x, 6);
    expect(arc[0].y).toBeCloseTo(equatorA.y, 6);
    expect(arc[0].z).toBeCloseTo(equatorA.z, 6);
    const last = arc[arc.length - 1];
    expect(last.x).toBeCloseTo(equatorB.x, 6);
    expect(last.y).toBeCloseTo(equatorB.y, 6);
    expect(last.z).toBeCloseTo(equatorB.z, 6);
  });

  it('produces segments + 1 samples, all unit-length', () => {
    const arc = sampleGreatCircleArc(equatorA, equatorB, 24);
    expect(arc).toHaveLength(25);
    for (const p of arc) {
      const len = Math.hypot(p.x, p.y, p.z);
      expect(len).toBeCloseTo(1, 6);
    }
  });

  it('progresses monotonically toward b (angle from a increases)', () => {
    const arc = sampleGreatCircleArc(equatorA, north, 32);
    let lastAngle = -1;
    for (const p of arc) {
      const dot = Math.max(-1, Math.min(1, equatorA.x * p.x + equatorA.y * p.y + equatorA.z * p.z));
      const angle = Math.acos(dot);
      expect(angle).toBeGreaterThanOrEqual(lastAngle - 1e-9);
      lastAngle = angle;
    }
  });

  it('handles near-identical points without producing NaN', () => {
    const almostA: Point = { x: 0.999999, y: 0.001, z: 0 };
    const arc = sampleGreatCircleArc(equatorA, almostA, 8);
    for (const p of arc) {
      expect(Number.isFinite(p.x)).toBe(true);
      expect(Number.isFinite(p.y)).toBe(true);
      expect(Number.isFinite(p.z)).toBe(true);
    }
  });

  it('handles antipodal points without producing NaN', () => {
    const arc = sampleGreatCircleArc(north, south, 8);
    for (const p of arc) {
      expect(Number.isFinite(p.x)).toBe(true);
      expect(Number.isFinite(p.y)).toBe(true);
      expect(Number.isFinite(p.z)).toBe(true);
      expect(Math.hypot(p.x, p.y, p.z)).toBeCloseTo(1, 3);
    }
  });

  it('defaults to 48 segments', () => {
    const arc = sampleGreatCircleArc(equatorA, equatorB);
    expect(arc).toHaveLength(49);
  });
});

describe('niceScaleBarLength', () => {
  it('picks a value from the 1/2/5 x 10^n family', () => {
    const { km } = niceScaleBarLength(10, 140);
    const normalized = km / Math.pow(10, Math.floor(Math.log10(km)));
    expect([1, 2, 5]).toContain(Math.round(normalized * 1e6) / 1e6);
  });

  it('keeps the bar within maxPixels', () => {
    for (const pixelsPerKm of [0.01, 0.5, 3, 27, 500]) {
      const { km, px } = niceScaleBarLength(pixelsPerKm, 140);
      expect(px).toBeLessThanOrEqual(140);
      expect(km).toBeGreaterThan(0);
      expect(px).toBeCloseTo(km * pixelsPerKm, 6);
    }
  });

  it('picks the largest nice value that fits', () => {
    // maxKm = 140 / 10 = 14 -> largest of {1,2,5,10} <= 14 is 10
    const { km } = niceScaleBarLength(10, 140);
    expect(km).toBe(10);
  });

  it('returns zero for non-finite or non-positive input', () => {
    expect(niceScaleBarLength(0, 140)).toEqual({ km: 0, px: 0 });
    expect(niceScaleBarLength(10, 0)).toEqual({ km: 0, px: 0 });
    expect(niceScaleBarLength(-5, 140)).toEqual({ km: 0, px: 0 });
  });
});

describe('computeScaleBar', () => {
  it('returns a sane, positive pixelsPerKm for a fitted equirectangular projection', () => {
    const width = 2000;
    const height = 1000;
    const projection = d3.geoEquirectangular().fitSize([width, height], { type: 'Sphere' } as d3.GeoPermissibleObjects);
    const result = computeScaleBar(projection, [0, 0], 6371);
    expect(result).not.toBeNull();
    expect(result!.pixelsPerKm).toBeGreaterThan(0);
    expect(Number.isFinite(result!.pixelsPerKm)).toBe(true);

    // Sanity check against the projection's own known scale: at the equator,
    // 1 degree of longitude spans (pi * R / 180) km on a sphere.
    const kmPerDegree = (Math.PI * 6371) / 180;
    const a = projection([-0.5, 0])!;
    const b = projection([0.5, 0])!;
    const pixelsPerDegree = Math.hypot(b[0] - a[0], b[1] - a[1]);
    expect(result!.pixelsPerKm).toBeCloseTo(pixelsPerDegree / kmPerDegree, 6);
  });

  it('returns null when the projection cannot place a sample point', () => {
    const failingProjection = () => null;
    const result = computeScaleBar(failingProjection, [0, 0], 6371);
    expect(result).toBeNull();
  });
});

describe('computeInterruptedScaleBar', () => {
  const width = 2000;
  const height = 1000;
  const smooth = d3.geoEquirectangular().fitSize([width, height], { type: 'Sphere' } as d3.GeoPermissibleObjects);

  it('agrees with computeScaleBar where the plane is not cut', () => {
    const plain = computeScaleBar(smooth, [0, 0], 6371);
    const guarded = computeInterruptedScaleBar(smooth, [0, 0], 6371);
    expect(guarded).not.toBeNull();
    expect(guarded!.pixelsPerKm).toBeCloseTo(plain!.pixelsPerKm, 9);
  });

  it('rejects a centre sitting on a cut', () => {
    // A seam at lon 0: everything east of it is displaced by a face width. This
    // is what a Dymaxion fold does to two points half a degree apart.
    const cut = (lonLat: [number, number]): [number, number] => {
      const [lon, lat] = lonLat;
      const p = smooth(lonLat)!;
      return [p[0] + (lon > 0 ? 400 : 0), p[1] + lat * 0];
    };
    expect(computeInterruptedScaleBar(cut, [0, 0], 6371)).toBeNull();
  });

  it('accepts a centre near, but not on, the same cut', () => {
    const cut = (lonLat: [number, number]): [number, number] => {
      const p = smooth(lonLat)!;
      return [p[0] + (lonLat[0] > 0 ? 400 : 0), p[1]];
    };
    // Both probes land east of the seam, so the step is smooth again.
    expect(computeInterruptedScaleBar(cut, [10, 0], 6371)).not.toBeNull();
  });

  it('catches a fold that leaves the two-point distance looking PLAUSIBLE', () => {
    // The failure the three-point check exists for. Here the displaced half
    // lands back near the other, so `computeScaleBar` alone returns a
    // believable number instead of an absurd one.
    const foldBack = (lonLat: [number, number]): [number, number] => {
      const p = smooth(lonLat)!;
      return [lonLat[0] > 0 ? p[0] - 3 : p[0], p[1]];
    };
    const naive = computeScaleBar(foldBack, [0, 0], 6371);
    expect(naive).not.toBeNull();
    expect(naive!.pixelsPerKm).toBeGreaterThan(0); // plausible, and wrong
    expect(computeInterruptedScaleBar(foldBack, [0, 0], 6371)).toBeNull();
  });

  it('returns null when any of the three samples fails to project', () => {
    expect(computeInterruptedScaleBar(() => null, [0, 0], 6371)).toBeNull();
  });
});
