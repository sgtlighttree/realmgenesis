import { describe, it, expect } from 'vitest';

import { ViewMode } from '../types';
import { getMapStyle } from '../utils/mapStyle';
import { runStyle } from '../utils/mapStyle/passes';
import { GrainSpec, HatchSpec, Substrate } from '../utils/mapStyle/substrate';
import { PlacedGlyph, StyleRenderContext } from '../utils/mapStyle/types';

// Recording substrate: every pass writes through this interface, so asserting
// on the call log tests the passes without a canvas or an SVG serializer.
interface Call { op: string; arg?: unknown; opacity?: number; count?: number }

const recorder = () => {
  const calls: Call[] = [];
  const sub: Substrate = {
    fillRect: (_x, _y, _w, _h, fill) => { calls.push({ op: 'fillRect', arg: fill }); },
    fillFeature: (_f, fill, opacity) => { calls.push({ op: 'fillFeature', arg: fill, opacity }); },
    strokeFeature: (_f, stroke) => { calls.push({ op: 'strokeFeature', arg: stroke }); },
    strokeSegments: (segs, stroke) => { calls.push({ op: 'strokeSegments', arg: stroke, count: segs.length }); },
    hatchRect: (_x, _y, _w, _h, spec: HatchSpec) => { calls.push({ op: 'hatchRect', arg: spec.color }); },
    hatchFeatures: (features, spec: HatchSpec) => {
      calls.push({ op: 'hatchFeatures', arg: spec.color, count: features.length });
    },
    grain: (spec: GrainSpec) => { calls.push({ op: 'grain', arg: spec.seed }); },
    drawGlyph: (g: PlacedGlyph, ink) => { calls.push({ op: 'drawGlyph', arg: ink }); void g; },
  };
  return { sub, calls };
};

// Four cells: two land, two ocean. seaLevel 0.55, matching the app default.
const SEA = 0.55;
const feature = () => ({ type: 'Feature', geometry: { type: 'Polygon', coordinates: [] } });

const makeCtx = (viewMode: ViewMode, over: Partial<StyleRenderContext> = {}): StyleRenderContext => {
  const cells = [
    { id: 0, height: 0.8, biome: 'Temperate Forest', temperature: 10, moisture: 0.5 },
    { id: 1, height: 0.6, biome: 'Steppe', temperature: 12, moisture: 0.2 },
    { id: 2, height: 0.3, biome: 'Ocean', temperature: 8, moisture: 1 },
    { id: 3, height: 0.1, biome: 'Deep Ocean', temperature: 6, moisture: 1 },
  ];
  return {
    world: {
      cells,
      params: { seaLevel: SEA, season: 0.5, seed: 't' },
      geoJson: { features: cells.map(feature) },
    },
    viewMode,
    widthPx: 1024,
    heightPx: 512,
    glyphs: [],
    shadeMap: null,
    colorCtx: { seaLevel: SEA },
    coastlines: [],
    ...over,
  } as unknown as StyleRenderContext;
};

const parchment = getMapStyle('parchment');

describe('runStyle', () => {
  it('draws nothing for a style with no passes', () => {
    const { sub, calls } = recorder();
    runStyle(getMapStyle('default'), makeCtx('biome'), sub);
    expect(calls).toEqual([]);
  });

  it('an empty pass list is the signal, not the style id', () => {
    // Render paths branch on passes.length so the invariant lives in one place.
    expect(getMapStyle('default').passes.length).toBe(0);
    expect(parchment.passes.length).toBeGreaterThan(0);
  });
});

describe('parchment passes', () => {
  // REGRESSION: an earlier draft ended the ocean pass with a FULL-BLEED
  // hatchRect. Under the 'bare' policy landPass paints nothing, so that hatch
  // sat directly on the parchment land.
  it('never hatches the whole canvas, only the ocean', () => {
    for (const mode of ['biome', 'satellite', 'height_bw', 'political', 'height'] as ViewMode[]) {
      const { sub, calls } = recorder();
      runStyle(parchment, makeCtx(mode), sub);
      expect(calls.some(c => c.op === 'hatchRect')).toBe(false);
      const hatch = calls.find(c => c.op === 'hatchFeatures');
      expect(hatch?.count).toBe(2); // exactly the two ocean cells
    }
  });

  it('leaves land unpainted in bare modes, and paints it otherwise', () => {
    // 'bare' — only the two ocean cells get a fill.
    const bare = recorder();
    runStyle(parchment, makeCtx('biome'), bare.sub);
    expect(bare.calls.filter(c => c.op === 'fillFeature').length).toBe(2);

    // 'categorical' — two ocean plus two land.
    const cat = recorder();
    runStyle(parchment, makeCtx('political'), cat.sub);
    expect(cat.calls.filter(c => c.op === 'fillFeature').length).toBe(4);

    // 'ramp' — likewise painted; a continuous ramp's information IS its fill,
    // so bare paper would render these modes blank.
    const ramp = recorder();
    runStyle(parchment, makeCtx('temperature'), ramp.sub);
    expect(ramp.calls.filter(c => c.op === 'fillFeature').length).toBe(4);
  });

  it('hatches the ocean AFTER shading it, so relief reads under the hatch', () => {
    const shadeMap = new Float32Array([0.8, 1, 0.75, 1]);
    const { sub, calls } = recorder();
    runStyle(parchment, makeCtx('biome', { shadeMap }), sub);
    const lastShade = calls.map(c => c.op).lastIndexOf('fillFeature');
    const hatchAt = calls.map(c => c.op).indexOf('hatchFeatures');
    expect(hatchAt).toBeGreaterThan(lastShade);
  });

  // REGRESSION: hillshade previously hardcoded #000000/#ffffff, the only
  // colours in the file not drawn from the palette. On bare parchment a neutral
  // grey shadow reads as dirt on the paper rather than as landform.
  it('shades in palette sepia, never pure black or white', () => {
    const shadeMap = new Float32Array([0.8, 1.1, 0.75, 1]);
    const { sub, calls } = recorder();
    runStyle(parchment, makeCtx('biome', { shadeMap }), sub);
    const shades = calls.filter(c => c.op === 'fillFeature' && c.opacity !== undefined);
    expect(shades.length).toBeGreaterThan(0);
    for (const s of shades) {
      expect([parchment.palette.shadow, parchment.palette.highlight]).toContain(s.arg);
      expect(s.arg).not.toBe('#000000');
      expect(s.arg).not.toBe('#ffffff');
    }
  });

  it('passes opacity as a number, never as an rgba() colour string', () => {
    // SVG 1.1 has no rgba() colour syntax; an rgba() fill would render
    // inconsistently in Illustrator and Inkscape.
    const shadeMap = new Float32Array([0.8, 1, 0.75, 1]);
    const { sub, calls } = recorder();
    runStyle(parchment, makeCtx('biome', { shadeMap }), sub);
    for (const c of calls) {
      expect(String(c.arg)).not.toContain('rgba(');
    }
  });

  it('draws every placed glyph in palette ink', () => {
    const glyphs: PlacedGlyph[] = [
      { x: 1, y: 2, kind: 'mountain', scale: 16, seedRot: 0, cellId: 0 },
      { x: 3, y: 4, kind: 'forest', scale: 14, seedRot: 0.1, cellId: 1 },
    ];
    const { sub, calls } = recorder();
    runStyle(parchment, makeCtx('biome', { glyphs }), sub);
    const drawn = calls.filter(c => c.op === 'drawGlyph');
    expect(drawn.length).toBe(2);
    expect(drawn.every(c => c.arg === parchment.palette.ink)).toBe(true);
  });

  it('lays paper down before anything else', () => {
    const { sub, calls } = recorder();
    runStyle(parchment, makeCtx('biome'), sub);
    expect(calls[0].op).toBe('fillRect');
    expect(calls[0].arg).toBe(parchment.palette.paper);
    expect(calls[1].op).toBe('grain');
  });
});
