import { describe, it, expect } from 'vitest';
import { isVisible } from '../utils/screenProject';

describe('analytic horizon test', () => {
  const cam = [0, 0, 3]; // camera on +Z
  it('front-facing point (near +Z) is visible', () => {
    expect(isVisible(0, 0, 1, cam[0], cam[1], cam[2])).toBe(true);
  });
  it('back-facing point (near −Z) is culled', () => {
    expect(isVisible(0, 0, -1, cam[0], cam[1], cam[2])).toBe(false);
  });
  it('grazing point at the horizon (near +X) is culled within eps', () => {
    // p=(1,0,0): cam−p=(-1,0,3) normalized dot p = -1/sqrt(10) < 0.08 → false
    expect(isVisible(1, 0, 0, cam[0], cam[1], cam[2])).toBe(false);
  });
});
