import { describe, it, expect, vi } from 'vitest';

import { Canvas2DSubstrate } from '../utils/mapStyle/substrateCanvas';
import { PlacedGlyph } from '../utils/mapStyle/types';

vi.stubGlobal('Path2D', class { constructor(_d: string) {} });

const makeCtx = () => ({
  save: vi.fn(), restore: vi.fn(), beginPath: vi.fn(), fill: vi.fn(),
  stroke: vi.fn(), fillRect: vi.fn(), clip: vi.fn(), rect: vi.fn(),
  moveTo: vi.fn(), lineTo: vi.fn(), translate: vi.fn(), rotate: vi.fn(),
  createPattern: vi.fn(() => ({})),
  fillStyle: '', strokeStyle: '', lineWidth: 0, globalAlpha: 1,
  lineCap: '', lineJoin: '',
}) as unknown as CanvasRenderingContext2D;

const glyph: PlacedGlyph = { x: 10, y: 20, kind: 'mountain', scale: 16, seedRot: 0, cellId: 0 };

describe('Canvas2DSubstrate', () => {
  it('fills a rect with the given colour', () => {
    const ctx = makeCtx();
    new Canvas2DSubstrate(ctx, (() => {}) as never, 100, 50).fillRect(0, 0, 100, 50, '#abcdef');
    expect(ctx.fillRect).toHaveBeenCalledWith(0, 0, 100, 50);
    expect(ctx.fillStyle).toBe('#abcdef');
  });

  it('strokes a glyph without filling it', () => {
    const ctx = makeCtx();
    new Canvas2DSubstrate(ctx, (() => {}) as never, 100, 50).drawGlyph(glyph, '#3b2f1c', 1);
    expect(ctx.stroke).toHaveBeenCalled();
    expect(ctx.fill).not.toHaveBeenCalled();
  });

  it('hatches at a cost independent of how many features are passed', () => {
    const clips = (n: number) => {
      const ctx = makeCtx();
      const features = Array.from({ length: n }, () => ({ type: 'Feature', geometry: {} }));
      new Canvas2DSubstrate(ctx, (() => {}) as never, 100, 50)
        .hatchFeatures(features, { color: '#000', spacingPx: 8, widthPx: 1, angleDeg: 45 });
      return (ctx.clip as unknown as { mock: { calls: unknown[] } }).mock.calls.length;
    };
    // The features union into ONE composite clip region, so cost does not scale
    // with cell count. An ocean is 13k-17k cells and Map2D re-renders on every
    // viewport change; per-feature clipping is what made that unusable.
    // (Two clips total: the composite region, plus hatchRect's own bounds.)
    expect(clips(5)).toBe(clips(50));
    expect(clips(5)).toBe(2);
  });

  it('draws nothing for an empty feature set', () => {
    const ctx = makeCtx();
    new Canvas2DSubstrate(ctx, (() => {}) as never, 100, 50)
      .hatchFeatures([], { color: '#000', spacingPx: 8, widthPx: 1, angleDeg: 45 });
    expect(ctx.clip).not.toHaveBeenCalled();
  });

  it('paints grain at constant cost, not per output pixel', () => {
    // The first version looped the whole output in 3px steps and called
    // fillRect per speck: 49,500 calls for one 1400x700 Map2D frame (on every
    // pan) and 1,679,374 for an 8192px PNG. A repeating tile is constant cost.
    //
    // A DOM stub is required, not optional: without `document` the tile builder
    // returns null and grain early-returns, so the assertion would pass
    // vacuously on 0 === 0 and prove nothing.
    const tileCtx = {
      createImageData: (w: number, h: number) => ({ data: new Uint8ClampedArray(w * h * 4) }),
      putImageData: vi.fn(),
    };
    vi.stubGlobal('document', {
      createElement: () => ({ width: 0, height: 0, getContext: () => tileCtx }),
    });
    try {
      const calls = (c: CanvasRenderingContext2D) =>
        (c.fillRect as unknown as { mock: { calls: unknown[] } }).mock.calls.length;

      const big = makeCtx();
      new Canvas2DSubstrate(big, (() => {}) as never, 8192, 4096)
        .grain({ seed: 'grain-test', opacity: 0.1, scale: 1 });
      const small = makeCtx();
      new Canvas2DSubstrate(small, (() => {}) as never, 100, 50)
        .grain({ seed: 'grain-test', opacity: 0.1, scale: 1 });

      // The pattern path really ran, and cost is identical at 80x the area.
      expect(big.createPattern).toHaveBeenCalled();
      expect(calls(big)).toBe(1);
      expect(calls(small)).toBe(1);
    } finally {
      vi.unstubAllGlobals();
    }
  });

  it('restores the context for every save', () => {
    const ctx = makeCtx();
    const sub = new Canvas2DSubstrate(ctx, (() => {}) as never, 100, 50);
    sub.hatchRect(0, 0, 100, 50, { color: '#000', spacingPx: 8, widthPx: 1, angleDeg: 45 });
    sub.grain({ seed: 'x', opacity: 0.1, scale: 1 });
    expect((ctx.save as unknown as { mock: { calls: unknown[] } }).mock.calls.length)
      .toBe((ctx.restore as unknown as { mock: { calls: unknown[] } }).mock.calls.length);
  });
});
