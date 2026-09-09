import { describe, it, expect } from 'vitest';

import { LABEL_CONFIG, LabelKind } from '../utils/labels';
import {
  LABEL_GLYPH_KINDS,
  LabelGlyphKind,
  labelGlyphMetrics,
  labelGlyphPathData,
  normalizedLabelGlyph,
} from '../utils/labelGlyphs';

describe('label glyph shapes', () => {
  it('defines a shape for every declared glyph kind', () => {
    for (const kind of LABEL_GLYPH_KINDS) {
      const shape = normalizedLabelGlyph(kind);
      expect(shape.length, `${kind} has no subpaths`).toBeGreaterThan(0);
      for (const sub of shape) {
        expect(sub.length, `${kind} has a degenerate subpath`).toBeGreaterThan(1);
      }
    }
  });

  // The glyph is sized off the label's fontSize, and drawMapLabels grows its
  // declutter rect by the returned metrics. If a shape escaped the unit box the
  // rect would understate it and labels would silently overlap.
  it('normalizes every shape to unit height, reporting its true width', () => {
    for (const kind of LABEL_GLYPH_KINDS) {
      const { width, height } = labelGlyphMetrics(kind);
      // Height is the normalized axis so every glyph carries the same visual
      // weight beside the name; width follows the shape's own aspect.
      expect(height, `${kind} is not normalized to unit height`).toBeCloseTo(1, 6);
      expect(width, `${kind} reports no width`).toBeGreaterThan(0);
      for (const sub of normalizedLabelGlyph(kind)) {
        for (const [x, y] of sub) {
          expect(Math.abs(x), `${kind} x escapes its reported width`).toBeLessThanOrEqual(width / 2 + 1e-9);
          expect(Math.abs(y), `${kind} y escapes its reported height`).toBeLessThanOrEqual(height / 2 + 1e-9);
        }
      }
    }
  });

  // Wide kinds must not come out shorter than square ones — that was the bug
  // long-axis normalization would have shipped: a hamlet ring standing twice
  // as tall as a mountain range's peak.
  it('gives the wide kinds the same height as the square ones', () => {
    for (const kind of LABEL_GLYPH_KINDS) {
      expect(labelGlyphMetrics(kind).height).toBeCloseTo(
        labelGlyphMetrics('hamlet').height, 6);
    }
    expect(labelGlyphMetrics('mountain').width).toBeGreaterThan(1);
  });

  it('every glyph declared in LABEL_CONFIG is a known kind', () => {
    const known = new Set<string>(LABEL_GLYPH_KINDS);
    for (const [kind, config] of Object.entries(LABEL_CONFIG)) {
      if (config.glyph) {
        expect(known.has(config.glyph), `${kind} -> ${config.glyph} is unknown`).toBe(true);
      }
    }
  });

  it('gives every settlement tier and geographic feature a glyph', () => {
    const needsOne: LabelKind[] = [
      'capital', 'city', 'town', 'village', 'hamlet',
      'range', 'desert', 'forest', 'sea', 'ocean', 'island', 'lake',
    ];
    for (const kind of needsOne) {
      expect(LABEL_CONFIG[kind].glyph, `${kind} has no glyph`).toBeTruthy();
    }
    // Territory and user pins stay text-only: a faction name spans a region and
    // has no single symbol, and markers carry their own C4 pin.
    expect(LABEL_CONFIG.faction.glyph).toBeUndefined();
    expect(LABEL_CONFIG.province.glyph).toBeUndefined();
    expect(LABEL_CONFIG.marker.glyph).toBeUndefined();
  });

  it('reads distinct silhouettes across the settlement ladder', () => {
    const tiers: LabelGlyphKind[] = ['capital', 'city', 'town', 'village', 'hamlet'];
    const seen = new Set<string>();
    for (const kind of tiers) {
      const d = labelGlyphPathData(kind, 0, 0, 10);
      expect(seen.has(d), `${kind} duplicates another tier`).toBe(false);
      seen.add(d);
    }
  });
});

describe('labelGlyphPathData', () => {
  it('emits SVG path data anchored at the requested centre', () => {
    const d = labelGlyphPathData('village', 100, 50, 10);
    expect(d.startsWith('M')).toBe(true);
    const { width, height } = labelGlyphMetrics('village');
    const coords = d.slice(1).split(/[ML]/).map(p => p.trim().split(' ').map(Number));
    for (const [x, y] of coords) {
      // Bounds come from the reported metrics, not a hardcoded square: width
      // follows each shape's aspect now that height is the normalized axis.
      expect(Math.abs(x - 100)).toBeLessThanOrEqual(width * 10 / 2 + 0.01);
      expect(Math.abs(y - 50)).toBeLessThanOrEqual(height * 10 / 2 + 0.01);
    }
  });

  it('translates with the anchor and scales with size', () => {
    const a = labelGlyphPathData('town', 0, 0, 10);
    const b = labelGlyphPathData('town', 0, 0, 20);
    expect(a).not.toEqual(b);
    const shifted = labelGlyphPathData('town', 7, 0, 10);
    expect(shifted).not.toEqual(a);
  });
});
