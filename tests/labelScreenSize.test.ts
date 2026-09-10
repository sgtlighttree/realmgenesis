import { describe, it, expect } from 'vitest';

import {
  drawMapLabels, MapLabel, LABEL_CONFIG, LABEL_SIZE_SCALE, smallestVisibleTier,
  ZOOM_THRESHOLDS, MAP2D_MIN_SCALE, MAP2D_MAX_SCALE,
} from '../utils/labels';
import { DEFAULT_LABEL_VISIBILITY } from '../types';
import { makeFakeCtx } from './helpers/overlayCanvas';

const label = (kind: MapLabel['kind'], name = 'Aa'): MapLabel => ({
  kind, name, position: { x: 1, y: 0, z: 0 }, priority: 1,
});
const project = (_p: { x: number; y: number; z: number }): [number, number] => [400, 300];
const VIS = {
  ...DEFAULT_LABEL_VISIBILITY,
  factions: true, capitals: true, cities: true, towns: true, villages: true, hamlets: true,
};

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
      expect(px * scale, `scale ${scale}`)
        .toBeCloseTo(LABEL_CONFIG.faction.baseFontSize * LABEL_SIZE_SCALE, 6);
    }
  });

  it('applies the floor in SCREEN pixels, not transformed ones', () => {
    // The floor only binds when base * LABEL_SIZE_SCALE falls under 8; assert
    // the floor's UNIT (screen px), which holds either way.
    // Scales are above ZOOM_THRESHOLDS.hamlet (3.2) so the LOD does not cull it.
    const expected = Math.max(8, LABEL_CONFIG.hamlet.baseFontSize * LABEL_SIZE_SCALE);
    for (const scale of [4, 6]) {
      const px = fontPx({ fontScale: 1 / scale, minFontSize: 8 / scale }, scale, 'hamlet');
      expect(px * scale, `scale ${scale}`).toBeCloseTo(expected, 6);
    }
  });

  it('scales off baseFontSize for callers that pass neither option', () => {
    const base = LABEL_CONFIG.faction.baseFontSize * LABEL_SIZE_SCALE;
    expect(fontPx({}, 1)).toBeCloseTo(base, 6);
    expect(fontPx({}, 2)).toBeCloseTo(base * 2, 6);
    // The historical min(scale, 2) saturation still applies.
    expect(fontPx({}, 6)).toBeCloseTo(base * 2, 6);
  });

  it('still honours an explicit fontScale on its own (export path)', () => {
    expect(fontPx({ fontScale: 2 }, 1))
      .toBeCloseTo(LABEL_CONFIG.faction.baseFontSize * LABEL_SIZE_SCALE * 2, 6);
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

  // The floors are expressed as fractions of `minFontSize` so they scale with a
  // zoomed caller. At the DEFAULT minFontSize of 8 they must still equal the
  // literals they replaced — 3 for the halo, 1 for the glyph stroke — or every
  // unscaled caller (exports, blit path, globe tenant) silently shifts.
  it('reproduces the old literal floors at the default minFontSize', () => {
    // A tiny fontScale drives fontSize to the floor, so both widths land on
    // their minimums rather than on their size-proportional terms.
    const w = widths({ fontScale: 0.001 }, 1);
    expect(Math.max(...w)).toBeCloseTo(3, 6);   // halo floor
    expect(Math.min(...w)).toBeCloseTo(1, 6);   // glyph-stroke floor
  });
});

// The zoom readout beside the scale bar reports this, so it has to agree with
// the LOD that actually runs inside drawMapLabels.
describe('smallestVisibleTier', () => {
  it('names the finest tier the given zoom admits', () => {
    expect(smallestVisibleTier(MAP2D_MAX_SCALE)).toBe('hamlet');
    expect(smallestVisibleTier(ZOOM_THRESHOLDS.hamlet)).toBe('hamlet');
    expect(smallestVisibleTier(ZOOM_THRESHOLDS.hamlet - 0.01)).toBe('village');
    expect(smallestVisibleTier(ZOOM_THRESHOLDS.village - 0.01)).toBe('town');
    expect(smallestVisibleTier(ZOOM_THRESHOLDS.town - 0.01)).toBe('city');
    expect(smallestVisibleTier(MAP2D_MIN_SCALE)).toBe('capital');
  });

  it('agrees with what drawMapLabels actually draws', () => {
    for (const scale of [0.6, 1, 2, 2.6, 3.2, 6]) {
      const tier = smallestVisibleTier(scale);
      if (!tier) continue;
      const ctx = makeFakeCtx();
      drawMapLabels(ctx, [label(tier)], project, scale, VIS, {});
      expect(ctx.ops.some(o => o.op === 'fillText'), `${tier} at ${scale}`).toBe(true);
    }
  });

  it('returns null below every threshold', () => {
    expect(smallestVisibleTier(MAP2D_MIN_SCALE - 0.01)).toBe(null);
  });
});
