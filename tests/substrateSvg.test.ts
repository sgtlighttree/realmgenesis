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

  it('reuses one pattern id for identical hatch specs', () => {
    const s = new SvgSubstrate((() => '') as never, 100, 50);
    const spec = { color: '#000', spacingPx: 8, widthPx: 1, angleDeg: 45 };
    s.hatchRect(0, 0, 50, 50, spec);
    s.hatchRect(50, 0, 50, 50, spec);
    expect(s.defs().match(/<pattern/g)?.length).toBe(1);
  });
});
