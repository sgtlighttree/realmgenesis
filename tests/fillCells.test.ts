import { describe, it, expect } from 'vitest';
import { Canvas2DSubstrate } from '../utils/mapStyle/substrateCanvas';
import { landPass } from '../utils/mapStyle/passes';
import { Substrate } from '../utils/mapStyle/substrate';
import { StyleRenderContext } from '../utils/mapStyle/types';
import { MapGeometryCache } from '../utils/mapCache';
import { ViewMode } from '../types';

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

describe('landPass cache fast path', () => {
  const feature = () => ({ type: 'Feature', geometry: { type: 'Polygon', coordinates: [] } });
  const SEA = 0.55;

  // 3 cells: 2 land (0, 1), 1 ocean (2).
  const cells = [
    { id: 0, height: 0.8, biome: 'Temperate Forest', temperature: 10, moisture: 0.5, center: { x: 1, y: 0, z: 0 } },
    { id: 1, height: 0.6, biome: 'Steppe', temperature: 12, moisture: 0.2, center: { x: 1, y: 0, z: 0 } },
    { id: 2, height: 0.3, biome: 'Ocean', temperature: 8, moisture: 1, center: { x: 1, y: 0, z: 0 } },
  ];

  const colorCache = ['#101010', '#202020', '#303030'];
  const geometryCache = {} as MapGeometryCache;

  const makeCtx = (over: Partial<StyleRenderContext> = {}): StyleRenderContext => ({
    world: {
      cells,
      params: { seaLevel: SEA, season: 0.5, seed: 't' },
      geoJson: { features: cells.map(feature) },
    },
    viewMode: 'height' as ViewMode,
    widthPx: 1024,
    heightPx: 512,
    glyphs: [],
    shadeMap: null,
    colorCtx: { seaLevel: SEA },
    coastlines: [],
    ...over,
  } as unknown as StyleRenderContext);

  const makeSub = () => {
    const fillCellsCalls: Array<{ indices: number[] | Uint32Array; colors: string[] }> = [];
    const fillFeatureCalls: unknown[] = [];
    const sub: Substrate = {
      fillRect: () => {},
      fillFeature: (f) => { fillFeatureCalls.push(f); },
      strokeFeature: () => {},
      strokeSegments: () => {},
      hatchRect: () => {},
      hatchFeatures: () => {},
      grain: () => {},
      withSphereClip: (draw) => draw(),
      drawGlyph: () => {},
      fillCells: (indices, colors) => { fillCellsCalls.push({ indices, colors: colors as string[] }); },
    };
    return { sub, fillCellsCalls, fillFeatureCalls };
  };

  const palette = {
    paper: '#fff', ink: '#000', inkLight: '#333', sea: '#00f', seaHatch: '#00c',
    coast: '#000', shadow: '#000', highlight: '#fff', ice: '#eef', desk: '#888',
    deskShadow: '#444',
  };

  it('RED->GREEN: calls fillCells once with land indices + colour cache when both caches present and mute is 0', () => {
    const { sub, fillCellsCalls, fillFeatureCalls } = makeSub();
    const pass = landPass(() => 'ramp', palette);
    const ctx = makeCtx({ colorCache, geometryCache });
    pass(ctx, sub);
    expect(fillCellsCalls).toHaveLength(1);
    expect(fillCellsCalls[0].indices).toEqual([0, 1]);
    expect(fillCellsCalls[0].colors).toBe(colorCache);
    expect(fillFeatureCalls).toHaveLength(0);
  });

  it('falls back to the per-feature path when caches are absent', () => {
    const { sub, fillCellsCalls, fillFeatureCalls } = makeSub();
    const pass = landPass(() => 'ramp', palette);
    const ctx = makeCtx();
    pass(ctx, sub);
    expect(fillCellsCalls).toHaveLength(0);
    expect(fillFeatureCalls.length).toBeGreaterThan(0);
  });

  it('falls back to the per-feature path for categorical (muted) policy even with caches present', () => {
    const { sub, fillCellsCalls, fillFeatureCalls } = makeSub();
    const pass = landPass(() => 'categorical', palette);
    const ctx = makeCtx({ colorCache, geometryCache });
    pass(ctx, sub);
    expect(fillCellsCalls).toHaveLength(0);
    expect(fillFeatureCalls.length).toBeGreaterThan(0);
  });

  it('falls back to the per-feature path for bare policy even with caches present', () => {
    const { sub, fillCellsCalls } = makeSub();
    const pass = landPass(() => 'bare', palette);
    const ctx = makeCtx({ colorCache, geometryCache });
    pass(ctx, sub);
    expect(fillCellsCalls).toHaveLength(0);
  });
});
