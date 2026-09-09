import { describe, it, expect } from 'vitest';

import {
  settlementQuota, SETTLEMENTS_PER_STERADIAN, MAX_SETTLEMENT_CELL_FRACTION,
} from '../utils/settlements';

// The old formula, kept here as the calibration reference: quota was
// `landCells * density * 0.02`, which made settlement count a function of MESH
// RESOLUTION rather than of the planet. 200k was the resolution it was tuned at.
const OLD = (landCells: number, density: number) => Math.max(1, Math.round(landCells * density * 0.02));
const REFERENCE_POINTS = 200_000;

describe('settlementQuota', () => {
  it('matches the old cell-count formula at the resolution it was tuned for', () => {
    for (const landFrac of [0.1, 0.155, 0.3, 0.5]) {
      const land = Math.round(REFERENCE_POINTS * landFrac);
      const q = settlementQuota(land, REFERENCE_POINTS, 1.0);
      expect(q, `landFrac ${landFrac}`).toBeCloseTo(OLD(land, 1.0), -1);
    }
  });

  // THE BUG: the same planet at 5k and 200k got wildly different settlement
  // counts, so a low-resolution world was nearly empty and the density slider
  // had almost nothing to scale. Area is a property of the planet; cell count
  // is a property of the tessellation.
  it('is resolution-independent wherever the cell cap does not bind', () => {
    const landFrac = 0.155;
    const counts = [50_000, 100_000, 200_000].map(pts =>
      settlementQuota(Math.round(pts * landFrac), pts, 1.0));
    for (const c of counts) expect(c).toBeCloseTo(counts[0], -1);
  });

  it('gives a low-resolution world far more settlements than the old formula', () => {
    const land = 1553; // measured in-app at 5,000 cells
    expect(OLD(land, 1.0)).toBe(31);
    const now = settlementQuota(land, 5_000, 1.0);
    expect(now).toBeGreaterThan(250);
    // ...but never more than the cap allows.
    expect(now).toBeLessThanOrEqual(land * MAX_SETTLEMENT_CELL_FRACTION * 1.0);
  });

  // The cap scales WITH density, or the slider would go dead at low resolution
  // the moment the cap binds — the same trap as the settlementDensity-0 finding.
  it('keeps density live even where the cap binds', () => {
    const land = 1553;
    const q = (d: number) => settlementQuota(land, 5_000, d);
    expect(q(0.5)).toBeLessThan(q(1.0));
    expect(q(1.0)).toBeLessThan(q(2.0));
    expect(q(2.0) / q(0.5)).toBeCloseTo(4, 0);
  });

  it('keeps density live at high resolution too', () => {
    const land = 31_000;
    const q = (d: number) => settlementQuota(land, 200_000, d);
    expect(q(2.0) / q(0.5)).toBeCloseTo(4, 0);
  });

  it('never exceeds the land available', () => {
    for (const density of [0.5, 1, 2]) {
      for (const [land, pts] of [[100, 1000], [1553, 5000], [31000, 200000]]) {
        expect(settlementQuota(land, pts, density)).toBeLessThanOrEqual(land);
      }
    }
  });

  // Matt's S36 ruling: density 0 is a deliberate capitals-only mode.
  it('floors at one settlement so density 0 stays capitals-only', () => {
    expect(settlementQuota(1553, 5_000, 0)).toBe(1);
  });

  it('exposes its calibration', () => {
    expect(SETTLEMENTS_PER_STERADIAN).toBeCloseTo((REFERENCE_POINTS * 0.02) / (4 * Math.PI), 6);
  });
});
