import { describe, it, expect } from 'vitest';
import { readFileSync } from 'node:fs';
import { resolve } from 'node:path';

// Guards the F2 one-frame-lag fix, which depends on MOUNT ORDER and therefore
// cannot be caught by any behavioural test: R3F runs useFrame callbacks in
// subscription order, and React subscribes sibling effects in mount order. If
// <GlobeSpin/> is ever moved after <ScreenOverlay/>, the overlay projects a
// rotation the spin has not yet advanced and every screen-space overlay trails
// the terrain by one frame again.
//
// That regression is invisible when the globe is still, sub-pixel when zoomed
// out, and grows with zoom — it reads as "parallax", which is exactly the
// symptom that consumed sessions 17, 18 and 19 chasing radius math. So the
// invariant is asserted where it can actually break: in the source text.
describe('screen-space overlay frame ordering', () => {
  const src = readFileSync(resolve(__dirname, '../components/WorldViewer.tsx'), 'utf8');

  it('mounts GlobeSpin before ScreenOverlay', () => {
    // Anchor on the props, not the bare tag: both names also appear inside the
    // comments that explain this very invariant.
    const spin = src.indexOf('<GlobeSpin target=');
    const overlay = src.indexOf('<ScreenOverlay world=');
    expect(spin, '<GlobeSpin/> not found in WorldViewer').toBeGreaterThan(-1);
    expect(overlay, '<ScreenOverlay/> not found in WorldViewer').toBeGreaterThan(-1);
    expect(spin).toBeLessThan(overlay);
  });

  it('keeps the globe spin out of a component that renders ScreenOverlay', () => {
    // The spin must not drift back into the parent: a parent's effects subscribe
    // AFTER its children's, so a spin there runs after the overlay again.
    expect(src).toMatch(/const GlobeSpin[^]*?useFrame\(/);
    // Exactly one useFrame may touch rotation.y, and it is GlobeSpin's.
    const spinWrites = src.match(/rotation\.y \+=/g) ?? [];
    expect(spinWrites).toHaveLength(1);
    const globeSpinStart = src.indexOf('const GlobeSpin');
    const globeSpinEnd = src.indexOf('\n};', globeSpinStart);
    const writeAt = src.indexOf('rotation.y +=');
    expect(writeAt).toBeGreaterThan(globeSpinStart);
    expect(writeAt).toBeLessThan(globeSpinEnd);
  });
});

describe('overlay reads the current frame', () => {
  const src = readFileSync(resolve(__dirname, '../components/overlays/ScreenOverlay.tsx'), 'utf8');

  it('forces camera and globe world matrices current before projecting', () => {
    // Three only recomputes world matrices inside render(), which runs after
    // every useFrame callback. Without these the overlay reads last frame's
    // transform — the other half of the same lag.
    expect(src).toContain('camera.updateMatrixWorld()');
    expect(src).toContain('updateWorldMatrix(true, false)');
  });
});
