import { describe, it, expect, vi } from 'vitest';

import { Canvas2DSubstrate } from '../utils/mapStyle/substrateCanvas';
import { PlacedGlyph } from '../utils/mapStyle/types';

vi.stubGlobal('Path2D', class { constructor(_d: string) {} });

const makeCtx = () => ({
  save: vi.fn(), restore: vi.fn(), beginPath: vi.fn(), fill: vi.fn(),
  stroke: vi.fn(), fillRect: vi.fn(), clip: vi.fn(), rect: vi.fn(),
  moveTo: vi.fn(), lineTo: vi.fn(), translate: vi.fn(), rotate: vi.fn(),
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

  it('restores the context for every save', () => {
    const ctx = makeCtx();
    const sub = new Canvas2DSubstrate(ctx, (() => {}) as never, 100, 50);
    sub.hatchRect(0, 0, 100, 50, { color: '#000', spacingPx: 8, widthPx: 1, angleDeg: 45 });
    sub.grain({ seed: 'x', opacity: 0.1, scale: 1 });
    expect((ctx.save as unknown as { mock: { calls: unknown[] } }).mock.calls.length)
      .toBe((ctx.restore as unknown as { mock: { calls: unknown[] } }).mock.calls.length);
  });
});
