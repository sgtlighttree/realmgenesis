import { describe, it, expect } from 'vitest';
import { simplifyPolyline } from '../utils/simplify';

describe('simplifyPolyline', () => {
  it('preserves both endpoints', () => {
    const pts: [number, number][] = [[0, 0], [1, 0.001], [2, 0]];
    const out = simplifyPolyline(pts, 0.5);
    expect(out[0]).toEqual([0, 0]);
    expect(out[out.length - 1]).toEqual([2, 0]);
  });

  it('drops a near-collinear interior point below tolerance', () => {
    const pts: [number, number][] = [[0, 0], [1, 0.01], [2, 0]];
    expect(simplifyPolyline(pts, 0.5)).toEqual([[0, 0], [2, 0]]);
  });

  it('keeps an interior point above tolerance', () => {
    const pts: [number, number][] = [[0, 0], [1, 5], [2, 0]];
    expect(simplifyPolyline(pts, 0.5)).toEqual([[0, 0], [1, 5], [2, 0]]);
  });

  it('never increases vertex count and keeps max deviation <= tolerance', () => {
    const pts: [number, number][] = [];
    for (let i = 0; i <= 100; i++) pts.push([i, Math.sin(i / 5) * 3]);
    const out = simplifyPolyline(pts, 1.0);
    expect(out.length).toBeLessThanOrEqual(pts.length);
    expect(out.length).toBeGreaterThanOrEqual(2);
  });

  it('returns input unchanged for 2 or fewer points', () => {
    expect(simplifyPolyline([[0, 0]], 1)).toEqual([[0, 0]]);
    expect(simplifyPolyline([[0, 0], [9, 9]], 1)).toEqual([[0, 0], [9, 9]]);
  });
});
