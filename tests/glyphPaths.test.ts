import { describe, it, expect } from 'vitest';

import { glyphPathData } from '../utils/mapStyle/glyphPaths';
import { GlyphKind, PlacedGlyph } from '../utils/mapStyle/types';

const KINDS: GlyphKind[] = ['mountain', 'hill', 'forest', 'conifer', 'dune', 'marsh'];

const glyph = (kind: GlyphKind, over: Partial<PlacedGlyph> = {}): PlacedGlyph => ({
  x: 100, y: 200, kind, scale: 16, seedRot: 0, cellId: 1, ...over,
});

describe('glyphPathData', () => {
  it('returns a non-empty path for every kind', () => {
    for (const kind of KINDS) {
      const d = glyphPathData(glyph(kind));
      expect(d.length).toBeGreaterThan(0);
      expect(d.startsWith('M')).toBe(true);
      expect(d).not.toMatch(/NaN|Infinity|undefined/);
    }
  });

  it('is deterministic', () => {
    for (const kind of KINDS) {
      expect(glyphPathData(glyph(kind))).toBe(glyphPathData(glyph(kind)));
    }
  });

  it('scales with the glyph scale', () => {
    const small = glyphPathData(glyph('mountain', { scale: 10 }));
    const large = glyphPathData(glyph('mountain', { scale: 20 }));
    expect(small).not.toBe(large);
  });

  it('keeps every coordinate near the glyph origin', () => {
    // A glyph must not stray more than ~2x its own scale from its anchor, or
    // thinning by centre distance would not prevent visible collisions.
    for (const kind of KINDS) {
      const g = glyph(kind);
      const nums = (glyphPathData(g).match(/-?\d+\.?\d*/g) ?? []).map(Number);
      for (let i = 0; i < nums.length; i += 2) {
        expect(Math.abs(nums[i] - g.x)).toBeLessThan(g.scale * 2);
        expect(Math.abs(nums[i + 1] - g.y)).toBeLessThan(g.scale * 2);
      }
    }
  });
});
