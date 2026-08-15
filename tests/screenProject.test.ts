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
    // p=(1,0,0): cam−p=(-1,0,3) normalized dot p = -1/sqrt(10) < eps → false
    expect(isVisible(1, 0, 0, cam[0], cam[1], cam[2])).toBe(false);
  });
  it('rendered-radius limb point is visible (F2 radius fix)', () => {
    // A terrain-limb point at r=1.05, ~68° off the camera axis. Under the old
    // unit-radius test with fat eps=0.08 this band was wrongly culled — the grid
    // stopped short of the silhouette. At true radius with a hair eps (0.005) it
    // is visible; the same DIRECTION at unit radius under the old fat eps culls.
    // This is the invariant that guards the queued tenants (roads/contours/etc.).
    const a = (68 * Math.PI) / 180;
    const r = 1.05;
    expect(isVisible(r * Math.sin(a), 0, r * Math.cos(a), 0, 0, 3)).toBe(true);
    expect(isVisible(Math.sin(a), 0, Math.cos(a), 0, 0, 3, 0.08)).toBe(false);
  });
});
