import { describe, it, expect } from 'vitest';
import { Canvas2DSubstrate } from '../utils/mapStyle/substrateCanvas';

if (typeof (globalThis as any).Path2D === 'undefined') {
  (globalThis as any).Path2D = class {} as any;
}

describe('Canvas2DSubstrate.fillCells', () => {
  it('fills each requested cell path with its colour', () => {
    const calls: Array<{ style: string; filled: boolean }> = [];
    const ctx = {
      save() {}, restore() {}, beginPath() {}, closePath() {},
      set fillStyle(v: string) { calls.push({ style: v, filled: false }); },
      get fillStyle() { return ''; },
      set strokeStyle(_v: string) {}, get strokeStyle() { return ''; },
      set lineWidth(_v: number) {}, set globalAlpha(_v: number) {},
      fill(_p?: any) { if (calls.length) calls[calls.length - 1].filled = true; },
      stroke(_p?: any) {}, translate() {}, scale() {},
    } as any;
    const paths = [new Path2D(), new Path2D(), new Path2D()];
    const sub = new Canvas2DSubstrate(ctx, (() => {}) as any, 100, 100, false, paths);
    sub.fillCells([0, 2], ['#112233', '#445566', '#778899']);
    // Two cells filled, with colours indexed by cell id (0 and 2).
    const fills = calls.filter((c) => c.filled);
    expect(fills.map((c) => c.style)).toEqual(['#112233', '#778899']);
  });
});
