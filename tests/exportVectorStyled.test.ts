import { describe, it, expect, beforeAll } from 'vitest';

import { WorldData } from '../types';
import { getMapStyle } from '../utils/mapStyle';
import { exportSVG } from '../utils/exportVector';
import { generateWorld } from '../utils/worldGen';
import { makeParams } from './helpers';

// End-to-end: a real generated world through the whole style chain to an SVG
// string. Every unit up to here was tested in isolation; this is the first
// check that they compose.
describe('exportSVG with the parchment style', () => {
  let world: WorldData;

  beforeAll(async () => {
    world = await generateWorld(makeParams({ points: 1500 }));
  }, 60000);

  it('still produces a valid SVG document', () => {
    const svg = exportSVG(world, 'biome', 'equirectangular', 1024, 'parchment');
    expect(svg.startsWith('<svg')).toBe(true);
    expect(svg.endsWith('</svg>')).toBe(true);
    expect(svg).not.toMatch(/NaN|undefined|Infinity/);
  });

  it('emits the parchment ingredients: paper, grain filter, hatch pattern, glyphs', () => {
    const svg = exportSVG(world, 'biome', 'equirectangular', 1024, 'parchment');
    expect(svg).toContain('feTurbulence');      // paper grain
    expect(svg).toContain('<pattern');          // ocean hatching
    expect(svg).toContain('url(#');             // pattern actually referenced
    expect(svg).toContain('stroke-linecap="round"'); // glyph / coastline strokes
  });

  it('paints paper instead of the black background rect', () => {
    const styled = exportSVG(world, 'biome', 'equirectangular', 1024, 'parchment');
    const plain = exportSVG(world, 'biome', 'equirectangular', 1024, 'default');
    // The default path opens with a full-bleed background rect; the styled path
    // must not, or the paper and its grain sit on top of a black page.
    expect(plain).toContain('fill="#050505"');
    expect(styled).not.toContain('fill="#050505"');
    // Read from the palette, never a literal — a palette tweak must not
    // break a test that is about the background rect.
    expect(styled).toContain(getMapStyle('parchment').palette.paper);
  });

  it('counter-transforms glyphs so they are not mirrored', () => {
    const svg = exportSVG(world, 'biome', 'equirectangular', 1024, 'parchment');
    // The body sits inside one mirror group; each glyph carries a nested flip
    // that composes back to the identity. Without it every glyph renders
    // backwards.
    expect(svg).toContain('<g id="style" transform="translate(1024,0) scale(-1,1)">');
    expect(svg).toContain('<g transform="translate(1024,0) scale(-1,1)"><path');
  });

  it('emits no rgba() FILL — SVG 1.1 has no such colour syntax', () => {
    const svg = exportSVG(world, 'biome', 'equirectangular', 1024, 'parchment');
    // Scoped to fills on purpose. The style system's own risk was baking alpha
    // into a fill colour, which is what `Substrate.fillFeature`'s explicit
    // `opacity` argument exists to prevent; it must emit `fill-opacity`.
    //
    // The document DOES still contain rgba() STROKES, from the pre-A3 river and
    // label groups in exportVector.ts (lines ~85, ~95, ~193). Those predate this
    // work and are out of scope here — recorded in HANDOFF as a follow-up.
    expect(svg).not.toContain('fill="rgba(');
    expect(svg).toContain('fill-opacity=');
  });

  it('suppresses glyphs in continuous-ramp modes', () => {
    // A ramp's information content IS its fill; glyphs would fight it.
    const ramp = exportSVG(world, 'temperature', 'equirectangular', 1024, 'parchment');
    expect(ramp).not.toContain('<g transform="translate(1024,0) scale(-1,1)"><path');
  });

  it('leaves default-style output byte-identical to before the style existed', () => {
    const a = exportSVG(world, 'biome', 'equirectangular', 1024, 'default');
    const b = exportSVG(world, 'biome', 'equirectangular', 1024);
    expect(a).toBe(b);
  });
});
