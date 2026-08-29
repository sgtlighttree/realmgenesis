import { describe, it, expect } from 'vitest';

import { ViewMode } from '../types';
import { getMapStyle } from '../utils/mapStyle';
import { runStyle } from '../utils/mapStyle/passes';
import { GrainSpec, HatchSpec, Substrate } from '../utils/mapStyle/substrate';
import { PlacedGlyph, StyleRenderContext } from '../utils/mapStyle/types';

// Recording substrate: every pass writes through this interface, so asserting
// on the call log tests the passes without a canvas or an SVG serializer.
interface Call { op: string; arg?: unknown; opacity?: number; count?: number; spacing?: number }

const recorder = () => {
  const calls: Call[] = [];
  const sub: Substrate = {
    fillRect: (_x, _y, _w, _h, fill) => { calls.push({ op: 'fillRect', arg: fill }); },
    fillFeature: (_f, fill, opacity) => { calls.push({ op: 'fillFeature', arg: fill, opacity }); },
    strokeFeature: (_f, stroke) => { calls.push({ op: 'strokeFeature', arg: stroke }); },
    strokeSegments: (segs, stroke) => { calls.push({ op: 'strokeSegments', arg: stroke, count: segs.length }); },
    hatchRect: (_x, _y, _w, _h, spec: HatchSpec) => { calls.push({ op: 'hatchRect', arg: spec.color }); },
    hatchFeatures: (features, spec: HatchSpec) => {
      calls.push({
        op: 'hatchFeatures',
        arg: spec.color,
        count: features.length,
        opacity: spec.opacity,
        spacing: spec.spacingPx,
      });
    },
    grain: (spec: GrainSpec) => { calls.push({ op: 'grain', arg: spec.seed }); },
    drawGlyph: (g: PlacedGlyph, ink) => { calls.push({ op: 'drawGlyph', arg: ink }); void g; },
  };
  return { sub, calls };
};

// Four cells: two land, two ocean. seaLevel 0.55, matching the app default.
const SEA = 0.55;
const feature = () => ({ type: 'Feature', geometry: { type: 'Polygon', coordinates: [] } });

// Cell centres are UNIT VECTORS; oceanHatchPass reads latitude off y.
const at = (latDeg: number) => ({
  x: Math.cos((latDeg * Math.PI) / 180),
  y: Math.sin((latDeg * Math.PI) / 180),
  z: 0,
});

const makeCtx = (viewMode: ViewMode, over: Partial<StyleRenderContext> = {}): StyleRenderContext => {
  const cells = [
    { id: 0, height: 0.8, biome: 'Temperate Forest', temperature: 10, moisture: 0.5, center: at(0) },
    { id: 1, height: 0.6, biome: 'Steppe', temperature: 12, moisture: 0.2, center: at(5) },
    { id: 2, height: 0.3, biome: 'Ocean', temperature: 8, moisture: 1, center: at(10) },
    { id: 3, height: 0.1, biome: 'Deep Ocean', temperature: 6, moisture: 1, center: at(-10) },
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

  // REGRESSION (reported by Matt): parchment "eliminates all ice terrain/ice
  // caps from view". Bare paper left ICE_CAP indistinguishable from temperate
  // land — 6.8% of land in a default world — and oceanFillPass painted a flat
  // tone without consulting getCellColor, so D3 seasonal sea ice went too.
  it('paints ice caps even in bare modes', () => {
    const ctx = makeCtx('biome');
    (ctx.world.cells as { biome: string }[])[0].biome = 'Ice Cap';
    const { sub, calls } = recorder();
    runStyle(parchment, ctx, sub);
    const fills = calls.filter(c => c.op === 'fillFeature');
    // Two ocean cells plus the one ice cap; the other land cell stays bare.
    expect(fills.length).toBe(3);
    expect(fills.some(c => c.arg === parchment.palette.ice)).toBe(true);
  });

  it('paints frozen sea as ice, not open water', () => {
    const ctx = makeCtx('biome');
    // Below SEAWATER_FREEZE_C at the neutral season.
    (ctx.world.cells as { temperature: number }[])[2].temperature = -40;
    const { sub, calls } = recorder();
    runStyle(parchment, ctx, sub);
    const fills = calls.filter(c => c.op === 'fillFeature');
    expect(fills.some(c => c.arg === parchment.palette.ice)).toBe(true);
    expect(fills.some(c => c.arg === parchment.palette.sea)).toBe(true);
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

// The globe bakes an equirectangular texture, where every longitude converges
// at the pole. The ocean hatch is a fixed-frequency pattern, so it wound round
// that singularity as a spiral rosette — the defect three rounds of UV fixes
// could only move, never remove. Fading the hatch out leaves plain paper there.
describe('oceanHatchPass polar fade', () => {
  const oceanAt = (...lats: number[]) => {
    const cells = lats.map((lat, i) => ({
      id: i, height: 0.2, biome: 'Ocean', temperature: 8, moisture: 1, center: at(lat),
    }));
    // Cast for the same reason makeCtx does: these passes read four fields, and
    // a full Cell would be noise.
    return { cells, geoJson: { features: cells.map(feature) } } as unknown as
      Pick<StyleRenderContext['world'], 'cells' | 'geoJson'>;
  };

  const hatchCalls = (over: Partial<StyleRenderContext>) => {
    const { sub, calls } = recorder();
    runStyle(parchment, makeCtx('biome', over), sub);
    return calls.filter(c => c.op === 'hatchFeatures');
  };

  it('hatches to the pole when no fade is set — a flat map has no singularity', () => {
    const world = { ...makeCtx('biome').world, ...oceanAt(0, 89) };
    const hatches = hatchCalls({ world });
    expect(hatches.length).toBe(1);
    expect(hatches[0].count).toBe(2);
    expect(hatches[0].opacity).toBeUndefined();
  });

  it('drops ocean cells past the end of the fade band', () => {
    const world = { ...makeCtx('biome').world, ...oceanAt(0, 89) };
    const hatches = hatchCalls({ world, polarHatchFadeDeg: [66, 82] });
    expect(hatches.length).toBe(1);
    expect(hatches[0].count).toBe(1); // the 89-degree cell is gone
    expect(hatches[0].opacity).toBe(1);
  });

  it('ramps opacity through the band instead of cutting at one latitude', () => {
    // A hard cut would just move the artifact onto the cut line.
    const world = { ...makeCtx('biome').world, ...oceanAt(0, 70, 74, 78) };
    const hatches = hatchCalls({ world, polarHatchFadeDeg: [66, 82] });
    const opacities = hatches.map(h => h.opacity);
    expect(opacities.length).toBeGreaterThan(1);
    // Strongest first, so the faded bands read as an edge to it.
    expect([...opacities].sort((a, b) => (b as number) - (a as number))).toEqual(opacities);
    expect(Math.max(...(opacities as number[]))).toBe(1);
    expect(Math.min(...(opacities as number[]))).toBeLessThan(1);
    // Every ocean cell inside the band is still hatched, none dropped.
    expect(hatches.reduce((n, h) => n + (h.count ?? 0), 0)).toBe(4);
  });

  it('treats both hemispheres alike', () => {
    const north = hatchCalls({
      world: { ...makeCtx('biome').world, ...oceanAt(75) },
      polarHatchFadeDeg: [66, 82],
    });
    const south = hatchCalls({
      world: { ...makeCtx('biome').world, ...oceanAt(-75) },
      polarHatchFadeDeg: [66, 82],
    });
    expect(south.map(h => h.opacity)).toEqual(north.map(h => h.opacity));
  });
});

// The globe joins the bake's left and right edges. A hatch is drawn in output
// pixels, so unless its pattern repeats a whole number of times across the
// width, the phase differs either side of the join and the diagonals jog down
// the antimeridian as a thin vertical line. Found by screenshotting at lon 180
// — every earlier screenshot was at lon 0, which is 180 degrees from the seam.
describe('oceanHatchPass seam alignment', () => {
  const spacingFor = (over: Partial<StyleRenderContext>) => {
    const { sub, calls } = recorder();
    runStyle(parchment, makeCtx('biome', over), sub);
    return calls.find(c => c.op === 'hatchFeatures')!.spacing!;
  };

  it('leaves the spacing alone for a flat map, which has no join', () => {
    // 2000px does not divide evenly, so the wrapped case must differ.
    expect(spacingFor({ widthPx: 2000, wrapsHorizontally: true }))
      .not.toBe(spacingFor({ widthPx: 2000 }));
  });

  it('makes the 45-degree hatch repeat a whole number of times across the width', () => {
    for (const widthPx of [1024, 2000, 2048, 4096, 8192]) {
      const spacing = spacingFor({ widthPx, wrapsHorizontally: true, lineScale: 0.5 });
      // A 45-degree hatch repeats every spacing / sin(45) pixels horizontally.
      const repeats = (widthPx * Math.SQRT1_2) / spacing;
      expect(repeats).toBeCloseTo(Math.round(repeats), 6);
    }
  });

  it('barely moves the density — it is the phase that matters, not the spacing', () => {
    const flat = spacingFor({ widthPx: 2048, lineScale: 0.5 });
    const wrapped = spacingFor({ widthPx: 2048, wrapsHorizontally: true, lineScale: 0.5 });
    expect(Math.abs(wrapped - flat) / flat).toBeLessThan(0.02);
  });
});
