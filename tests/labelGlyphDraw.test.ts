import { describe, it, expect } from 'vitest';

import { drawMapLabels, glyphLayout, MapLabel } from '../utils/labels';
import { labelGlyphPathData } from '../utils/labelGlyphs';
import { DEFAULT_LABEL_VISIBILITY } from '../types';
import { makeFakeCtx } from './helpers/overlayCanvas';

// Screen positions are supplied directly so the geometry under test is the
// label layout, not a projection. `project` maps a label's x to a screen x.
const projectAt = (label: { x: number; y: number; z: number }): [number, number] =>
  [label.x, label.y];

const label = (kind: MapLabel['kind'], name: string, x: number): MapLabel => ({
  kind, name, position: { x, y: 100, z: 0 }, priority: 1,
});

/** Same, but varying the screen Y — the axis the stacked glyph grows along. */
const stacked = (kind: MapLabel['kind'], name: string, y: number): MapLabel => ({
  kind, name, position: { x: 200, y, z: 0 }, priority: 1,
});

// `faction` is the no-glyph control rather than `province`: both lack a glyph,
// but province's ZOOM_THRESHOLDS entry (1.2) filters it at the scale used here,
// so it would test the LOD gate instead of the glyph layout.
const VIS = { ...DEFAULT_LABEL_VISIBILITY, capitals: true, factions: true, geography: true };

const draw = (labels: MapLabel[], opts = {}) => {
  const ctx = makeFakeCtx();
  drawMapLabels(ctx, labels, projectAt, 1, VIS, opts);
  return ctx;
};

describe('label glyphs in drawMapLabels', () => {
  it('draws a glyph for a kind that declares one', () => {
    const ctx = draw([label('capital', 'Aa', 200)]);
    expect(ctx.ops.some(o => o.op === 'moveTo')).toBe(true);
    expect(ctx.ops.filter(o => o.op === 'stroke').length).toBeGreaterThan(0);
  });

  it('draws no glyph for a kind that declares none', () => {
    const ctx = draw([label('faction', 'Aa', 200)]);
    expect(ctx.ops.some(o => o.op === 'moveTo')).toBe(false);
  });

  // The glyph stacks ABOVE the name and the union is centred on the anchor, so
  // the text drops below the anchor and the glyph rises above it. Height-only
  // arithmetic — see the `glyphLayout` doc for why width is unusable here.
  it('stacks the glyph above the name, union centred on the anchor', () => {
    const withGlyph = draw([label('capital', 'Aa', 200)]);
    const without = draw([label('faction', 'Aa', 200)]);
    const textY = (c: ReturnType<typeof draw>) =>
      c.ops.find(o => o.op === 'fillText')!.args[1];
    const textX = (c: ReturnType<typeof draw>) =>
      c.ops.find(o => o.op === 'fillText')!.args[0];

    // The name stays horizontally centred on the anchor in both cases.
    expect(textX(without)).toBeCloseTo(200, 6);
    expect(textX(withGlyph)).toBeCloseTo(200, 6);

    expect(textY(without)).toBeCloseTo(100, 6);
    expect(textY(withGlyph)).toBeGreaterThan(100);

    // Every glyph vertex sits above the name it belongs to.
    const glyphYs = withGlyph.ops
      .filter(o => o.op === 'moveTo' || o.op === 'lineTo')
      .map(o => o.args[1]);
    expect(Math.max(...glyphYs)).toBeLessThan(textY(withGlyph));
  });

  // THE REGRESSION THIS FILE EXISTS FOR. The declutter rectangle was built from
  // text metrics alone; a glyph drawn outside it overlaps neighbouring labels
  // while the algorithm believes there is no collision. The glyph stacks
  // vertically, so the rectangle has to grow in Y: two capitals 26px apart
  // VERTICALLY have non-overlapping text boxes (13.2px + 8 padding = 21.2) but
  // overlapping glyph+text boxes (~10.5 + 2.75 + 13.2 + 8 = 34.4).
  it('reserves the glyph in the declutter rectangle', () => {
    const glyphed = draw([stacked('capital', 'Aa', 100), stacked('capital', 'Bb', 126)]);
    expect(glyphed.ops.filter(o => o.op === 'fillText').map(o => o.text)).toEqual(['Aa']);

    // Same spacing, no glyph: both fit, proving the rejection above came from
    // the glyph's height and not from the spacing being too tight anyway.
    const plain = draw([stacked('faction', 'Aa', 100), stacked('faction', 'Bb', 126)]);
    expect(plain.ops.filter(o => o.op === 'fillText').map(o => o.text)).toEqual(['AA', 'BB']);
  });

  // A3's relief styles already paint mountains, trees and dunes across the
  // land, so a `range` label glyph would land a second mountain beside one of
  // theirs. Settlement glyphs survive — no style paints towns.
  it('suppresses geographic glyphs when the style already paints relief', () => {
    const geo = draw([label('range', 'Aa', 200)], { geoGlyphs: false });
    expect(geo.ops.some(o => o.op === 'moveTo')).toBe(false);
    expect(geo.ops.filter(o => o.op === 'fillText').map(o => o.text)).toEqual(['Aa']);

    const town = draw([label('capital', 'Aa', 200)], { geoGlyphs: false });
    expect(town.ops.some(o => o.op === 'moveTo')).toBe(true);
  });

  it('keeps geographic glyphs by default', () => {
    const ctx = draw([label('range', 'Aa', 200)]);
    expect(ctx.ops.some(o => o.op === 'moveTo')).toBe(true);
  });
});

// The SVG export builds its label glyphs from `glyphLayout` + `labelGlyphPathData`
// rather than from the canvas ops, so nothing structural forces the two to
// agree. This pins the numbers: if drawMapLabels ever stops routing through
// glyphLayout, the raster and vector maps drift and this fails.
describe('raster/vector glyph parity', () => {
  it('draws the canvas glyph at the coordinates the SVG path emits', () => {
    const anchorX = 200;
    const anchorY = 100;
    const ctx = draw([label('capital', 'Aa', anchorX)]);

    // capital: baseFontSize 11, effectiveFontScale min(1,2) = 1.
    const layout = glyphLayout(11, 'capital');
    const d = labelGlyphPathData('capital', anchorX, anchorY + layout.glyphDy, layout.glyphSize);
    const svgPoints = d
      .split(/(?=[ML])/)
      .map(seg => seg.slice(1).trim().split(' ').map(Number));

    const canvasPoints = ctx.ops
      .filter(o => o.op === 'moveTo' || o.op === 'lineTo')
      .map(o => o.args)
      // Halo pass then ink pass draw the same path twice; compare one of them.
      .slice(0, svgPoints.length);

    expect(canvasPoints).toHaveLength(svgPoints.length);
    canvasPoints.forEach(([x, y], i) => {
      expect(x).toBeCloseTo(svgPoints[i][0], 1);
      expect(y).toBeCloseTo(svgPoints[i][1], 1);
    });
  });
});
