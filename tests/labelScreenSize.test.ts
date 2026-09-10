import { describe, it, expect } from 'vitest';

import { drawMapLabels, MapLabel, LABEL_CONFIG } from '../utils/labels';
import { DEFAULT_LABEL_VISIBILITY } from '../types';
import { makeFakeCtx } from './helpers/overlayCanvas';

const label = (kind: MapLabel['kind'], name = 'Aa'): MapLabel => ({
  kind, name, position: { x: 1, y: 0, z: 0 }, priority: 1,
});
const project = (_p: { x: number; y: number; z: number }): [number, number] => [400, 300];
const VIS = { ...DEFAULT_LABEL_VISIBILITY, factions: true, capitals: true };

const fontPx = (opts: object, scale: number, kind: MapLabel['kind'] = 'faction'): number => {
  const ctx = makeFakeCtx();
  const fonts: string[] = [];
  const orig = ctx.fillText;
  (ctx as unknown as { fillText: typeof orig }).fillText = function (t: string, x: number, y: number) {
    fonts.push(ctx.font); return orig.call(ctx, t, x, y);
  };
  drawMapLabels(ctx, [label(kind)], project, scale, VIS, opts);
  const m = /(\d+(?:\.\d+)?)px/.exec(fonts[0] ?? '');
  return m ? parseFloat(m[1]) : NaN;
};

// THE BUG: on the F3 vector path the canvas transform already carries the zoom
// (`setTransform(dpr*scale, ...)` in VectorOverlay), and drawMapLabels was ALSO
// given `scale` as its font scale — so zoom was counted twice and labels grew
// quadratically. Line widths and marker pins already divided by `s`; labels
// did not. `fontScale`/`minFontSize` let a caller under a scaled transform ask
// for constant SCREEN size.
describe('label size under a scaled transform', () => {
  it('keeps a constant on-screen size when the caller divides by its transform scale', () => {
    for (const scale of [1, 2, 4, 6]) {
      const px = fontPx({ fontScale: 1 / scale, minFontSize: 8 / scale }, scale);
      // Multiply back by the transform the caller applies: the on-screen size.
      expect(px * scale, `scale ${scale}`).toBeCloseTo(LABEL_CONFIG.faction.baseFontSize, 6);
    }
  });

  it('applies the floor in SCREEN pixels, not transformed ones', () => {
    // hamlet's base is 7, below the 8px floor: it should clamp to 8 on screen.
    // Scales are above ZOOM_THRESHOLDS.hamlet (3.2) so the LOD does not cull it.
    for (const scale of [4, 6]) {
      const px = fontPx({ fontScale: 1 / scale, minFontSize: 8 / scale }, scale, 'hamlet');
      expect(px * scale, `scale ${scale}`).toBeCloseTo(8, 6);
    }
  });

  it('is unchanged for callers that pass neither option', () => {
    expect(fontPx({}, 1)).toBe(14);
    expect(fontPx({}, 2)).toBe(28);
    // The historical min(scale, 2) saturation still applies.
    expect(fontPx({}, 6)).toBe(28);
  });

  it('still honours an explicit fontScale on its own (export path)', () => {
    expect(fontPx({ fontScale: 2 }, 1)).toBe(28);
  });
});

// The halo and the glyph stroke had the SAME floor-in-the-wrong-unit defect as
// the font size: `Math.max(3, ...)` and `Math.max(1, ...)` are transformed
// units, so under a scaled transform they pinned and the halo rendered at
// `3 * scale` — 16px of black around 2.5px of text at full zoom, which is the
// fat blob behind every name. Floors now derive from `minFontSize`, which the
// caller already expresses in the right unit.
describe('halo and glyph stroke widths', () => {
  const widths = (opts: object, scale: number): number[] => {
    const ctx = makeFakeCtx();
    const seen: number[] = [];
    const origStroke = ctx.stroke;
    const origStrokeText = ctx.strokeText;
    (ctx as unknown as { stroke: typeof origStroke }).stroke = function () {
      seen.push(ctx.lineWidth); return (origStroke as () => void).call(ctx);
    };
    (ctx as unknown as { strokeText: typeof origStrokeText }).strokeText = function (t: string, x: number, y: number) {
      seen.push(ctx.lineWidth); return origStrokeText.call(ctx, t, x, y);
    };
    drawMapLabels(ctx, [label('capital')], project, scale, VIS, opts);
    return seen;
  };

  it('holds every label stroke at a constant on-screen width across zoom', () => {
    const at = (scale: number) =>
      widths({ fontScale: 1 / scale, minFontSize: 8 / scale }, scale).map(w => +(w * scale).toFixed(4));
    const base = at(1);
    expect(base.length).toBeGreaterThan(0);
    for (const scale of [2, 4, 6]) {
      expect(at(scale), `scale ${scale}`).toEqual(base);
    }
  });

  it('leaves the default (unscaled) callers byte-identical', () => {
    // capital at scale 1: fontSize 11. Halo max(3, 11*haloScale); glyph max(1, 0.99).
    const w = widths({}, 1);
    expect(w.some(x => x >= 3)).toBe(true);
    expect(Math.min(...w)).toBeCloseTo(1, 6);
  });
});
