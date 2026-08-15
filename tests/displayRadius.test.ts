import { describe, it, expect } from 'vitest';
import { displayRadius } from '../utils/displayRadius';

describe('displayRadius', () => {
  it('applies per-cell relief when not smooth (matches the globe mesh)', () => {
    expect(displayRadius(0, false)).toBeCloseTo(1.0);
    expect(displayRadius(1, false)).toBeCloseTo(1.05);
    expect(displayRadius(0.55, false)).toBeCloseTo(1.0275); // sea level datum
  });
  it('collapses every cell to unit radius when smooth', () => {
    expect(displayRadius(0, true)).toBe(1);
    expect(displayRadius(1, true)).toBe(1);
    expect(displayRadius(-0.5, true)).toBe(1);
  });
  it('adds the offset in both modes', () => {
    expect(displayRadius(1, false, 0.005)).toBeCloseTo(1.055);
    expect(displayRadius(1, true, 0.005)).toBeCloseTo(1.005);
  });
});
