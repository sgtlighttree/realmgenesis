import { describe, it, expect } from 'vitest';
import {
  DEFAULT_MAX_ELEVATION_M,
  MAX_DEPTH_M,
  elevationMetres,
  formatElevation,
} from '../utils/datum';

// D8a presentation datum. The height field is normalized 0-1 with seaLevel as
// the coastline; these functions turn it into real metres above/below sea level
// against a FIXED maximum (comparable between worlds — see ROADMAP D8a).
describe('elevationMetres', () => {
  const sea = 0.55; // the app default

  it('reads 0 m exactly at the coastline', () => {
    expect(elevationMetres(sea, sea, 9000)).toBe(0);
  });

  it('reads the full datum at the highest peak (height = 1)', () => {
    expect(elevationMetres(1, sea, 9000)).toBeCloseTo(9000, 6);
  });

  it('reads full ocean depth at the deepest floor (height = 0)', () => {
    expect(elevationMetres(0, sea, 9000)).toBeCloseTo(-MAX_DEPTH_M, 6);
  });

  it('applies the quadratic hypsometric curve to land, not raw height', () => {
    // height 0.70 with sea 0.55: frac = (0.70-0.55)/(1-0.55) = 1/3; curved
    // frac^2 = 1/9 of 9000 = 1000 m (linear would have been 3000 m).
    expect(elevationMetres(0.7, sea, 9000)).toBeCloseTo(1000, 6);
  });

  it('scales ocean depth linearly (curve is land-only)', () => {
    // Ocean stays linear: real ocean floor is bottom-heavy abyssal plain, so a
    // concave curve made it implausibly shallow. h 0.50, sea 0.55:
    // (0.50-0.55)/0.55 * 11000 = -1000 m.
    expect(elevationMetres(0.5, sea, 9000)).toBeCloseTo(-1000, 6);
  });

  it('is monotonic increasing across the whole range', () => {
    let prev = -Infinity;
    for (let h = 0; h <= 1.0001; h += 0.05) {
      const m = elevationMetres(h, sea, 9000);
      expect(m).toBeGreaterThan(prev);
      prev = m;
    }
  });

  it('rescales proportionally with the adjustable datum', () => {
    // Doubling the datum doubles the reported land elevation (display-only).
    expect(elevationMetres(0.7, sea, 18000)).toBeCloseTo(2000, 6);
  });

  it('leaves ocean depth unaffected by the elevation datum', () => {
    // MAX_DEPTH_M is fixed; changing the elevation datum must not move depths.
    expect(elevationMetres(0.2, sea, 9000)).toBe(elevationMetres(0.2, sea, 18000));
  });
});

describe('formatElevation', () => {
  const sea = 0.55;

  it('labels the coastline as sea level', () => {
    expect(formatElevation(sea, sea, 9000)).toBe('Sea level');
  });

  it('formats land elevation with thousands separators and a unit', () => {
    // Curved: frac 1/3 -> 1/9 * 9000 = 1000 m.
    expect(formatElevation(0.7, sea, 9000)).toBe('1,000 m');
  });

  it('formats ocean depth as a negative metre value', () => {
    expect(formatElevation(0, sea, 9000)).toBe(`${(-MAX_DEPTH_M).toLocaleString('en-US')} m`);
  });

  it('exports a 9000 m default datum', () => {
    expect(DEFAULT_MAX_ELEVATION_M).toBe(9000);
    expect(elevationMetres(1, sea)).toBeCloseTo(9000, 6);
  });
});
