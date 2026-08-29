import { describe, it, expect } from 'vitest';

import { SvgSubstrate } from '../utils/mapStyle/substrateSvg';
import { PlacedGlyph } from '../utils/mapStyle/types';

const glyph: PlacedGlyph = { x: 10, y: 20, kind: 'mountain', scale: 16, seedRot: 0, cellId: 0 };

describe('SvgSubstrate', () => {
  it('emits a rect for fillRect', () => {
    const s = new SvgSubstrate((() => '') as never, 100, 50);
    s.fillRect(0, 0, 100, 50, '#abcdef');
    expect(s.body()).toContain('<rect');
    expect(s.body()).toContain('#abcdef');
  });

  it('registers a pattern in defs for hatching, and references it', () => {
    const s = new SvgSubstrate((() => '') as never, 100, 50);
    s.hatchRect(0, 0, 100, 50, { color: '#000', spacingPx: 8, widthPx: 1, angleDeg: 45 });
    expect(s.defs()).toContain('<pattern');
    expect(s.body()).toMatch(/fill="url\(#/);
  });

  it('registers a turbulence filter for grain', () => {
    const s = new SvgSubstrate((() => '') as never, 100, 50);
    s.grain({ seed: 'x', opacity: 0.12, scale: 1 });
    expect(s.defs()).toContain('feTurbulence');
  });

  it('emits a stroked, unfilled path for a glyph', () => {
    const s = new SvgSubstrate((() => '') as never, 100, 50);
    s.drawGlyph(glyph, '#3b2f1c', 1);
    expect(s.body()).toContain('fill="none"');
    expect(s.body()).toContain('stroke="#3b2f1c"');
  });

  it('emits fill-opacity rather than an rgba() fill', () => {
    const s = new SvgSubstrate((() => 'M0 0L1 1') as never, 100, 50);
    s.fillFeature({ type: 'Feature', geometry: {} }, '#000000', 0.25);
    expect(s.body()).toContain('fill-opacity="0.250"');
    expect(s.body()).not.toContain('rgba(');
  });

  it('hatches many features as ONE path element', () => {
    const s = new SvgSubstrate((() => 'M0 0L1 1') as never, 100, 50);
    const features = Array.from({ length: 20 }, () => ({ type: 'Feature', geometry: {} }));
    s.hatchFeatures(features, { color: '#000', spacingPx: 8, widthPx: 1, angleDeg: 45 });
    // 13k-17k separate <path> elements would make the SVG unopenable in a
    // vector editor for no visual gain.
    expect(s.body().match(/<path/g)?.length).toBe(1);
    expect(s.body()).toMatch(/fill="url\(#/);
  });

  it('emits nothing for an empty feature set', () => {
    const s = new SvgSubstrate((() => 'M0 0L1 1') as never, 100, 50);
    s.hatchFeatures([], { color: '#000', spacingPx: 8, widthPx: 1, angleDeg: 45 });
    expect(s.body()).toBe('');
  });

  it('counter-transforms glyphs on a mirrored surface', () => {
    // Map2D and exportSVG both draw inside a horizontal flip. A glyph drawn
    // straight into that comes out backwards; a nested flip composes to the
    // identity, and the x-flip puts it back over its own cell.
    const plain = new SvgSubstrate((() => 'M0 0') as never, 100, 50, false);
    plain.drawGlyph(glyph, '#000', 1);
    expect(plain.body()).not.toContain('<g transform=');

    const mirrored = new SvgSubstrate((() => 'M0 0') as never, 100, 50, true);
    mirrored.drawGlyph(glyph, '#000', 1);
    expect(mirrored.body()).toContain('<g transform="translate(100,0) scale(-1,1)">');
    // glyph.x is 10, so the flipped anchor is 100 - 10 = 90.
    expect(mirrored.body()).toContain('90');
  });

  it('reuses one pattern id for identical hatch specs', () => {
    const s = new SvgSubstrate((() => '') as never, 100, 50);
    const spec = { color: '#000', spacingPx: 8, widthPx: 1, angleDeg: 45 };
    s.hatchRect(0, 0, 50, 50, spec);
    s.hatchRect(50, 0, 50, 50, spec);
    expect(s.defs().match(/<pattern/g)?.length).toBe(1);
  });
});
