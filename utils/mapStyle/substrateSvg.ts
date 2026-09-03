import { Point3, toLonLat } from '../geo';
import { glyphPathData } from './glyphPaths';
import { GeoFeatureLike, GrainSpec, HatchSpec, Substrate } from './substrate';
import { PlacedGlyph } from './types';

type PathStringGenerator = (object: unknown) => string | null;

const esc = (s: string): string => s.replace(/[<>&"]/g, c =>
  ({ '<': '&lt;', '>': '&gt;', '&': '&amp;', '"': '&quot;' }[c] as string));

/**
 * Accumulates SVG. `defs()` and `body()` are spliced into the document by the
 * export path.
 *
 * Nothing here is a fallback or an approximation of the raster look: SVG has
 * native equivalents for both parchment ingredients that seem raster-only —
 * hatching is `<pattern>` and paper grain is `<filter><feTurbulence>`.
 */
export class SvgSubstrate implements Substrate {
  private defsParts: string[] = [];
  private bodyParts: string[] = [];
  private patternIds = new Map<string, string>();
  private grainId: string | null = null;
  private sphereClipId: string | null = null;

  /**
   * `mirrored` says the caller wraps `body()` in Map2D's/exportSVG's horizontal
   * flip group. Only glyphs care — see `drawGlyph`.
   */
  constructor(
    private path: PathStringGenerator,
    private width: number,
    private height: number,
    private mirrored = false,
    private cellFeatures: GeoFeatureLike[] = [],
  ) {}

  defs(): string {
    return this.defsParts.length ? `<defs>${this.defsParts.join('')}</defs>` : '';
  }

  body(): string {
    return this.bodyParts.join('');
  }

  fillRect(x: number, y: number, w: number, h: number, fill: string): void {
    this.bodyParts.push(`<rect x="${x}" y="${y}" width="${w}" height="${h}" fill="${esc(fill)}"/>`);
  }

  fillFeature(feature: GeoFeatureLike, fill: string, opacity = 1): void {
    const d = this.path(feature);
    if (!d) return;
    // fill-opacity, NOT an rgba() fill — SVG 1.1 has no rgba() colour syntax.
    //
    // The hairline stroke closes the antialiasing seam between adjacent Voronoi
    // polygons (the pre-A3 renderCellPaths does the same). It MUST carry the
    // same opacity as the fill: without stroke-opacity a translucent pass drew
    // a full-strength outline on every cell, which rendered the whole map as a
    // dark honeycomb mesh — most visible where hillshade fills at 5% opacity.
    const op = opacity < 1
      ? ` fill-opacity="${opacity.toFixed(3)}" stroke-opacity="${opacity.toFixed(3)}"`
      : '';
    this.bodyParts.push(
      `<path d="${d}" fill="${esc(fill)}"${op} stroke="${esc(fill)}" stroke-width="0.5"/>`,
    );
  }

  hatchFeatures(features: GeoFeatureLike[], spec: HatchSpec): void {
    const ds: string[] = [];
    for (const f of features) {
      const d = this.path(f);
      if (d) ds.push(d);
    }
    if (!ds.length) return;
    // One <path> holding every subpath, not one element per feature: an ocean
    // is 13k-17k cells, and that many elements makes the file unopenable in a
    // vector editor for no visual gain.
    const id = this.hatchPatternId(spec);
    this.bodyParts.push(`<path d="${ds.join('')}" fill="url(#${id})"/>`);
  }

  hatchCells(indices: number[] | Uint32Array, spec: HatchSpec): void {
    const feats: GeoFeatureLike[] = [];
    for (let j = 0; j < indices.length; j++) {
      const f = this.cellFeatures?.[indices[j]];
      if (f) feats.push(f);
    }
    this.hatchFeatures(feats, spec);
  }

  strokeFeature(feature: GeoFeatureLike, stroke: string, width: number): void {
    const d = this.path(feature);
    if (!d) return;
    this.bodyParts.push(`<path d="${d}" fill="none" stroke="${esc(stroke)}" stroke-width="${width}"/>`);
  }

  strokeSegments(segments: Array<[Point3, Point3]>, stroke: string, width: number): void {
    // Point3 is a UNIT VECTOR; GeoJSON coordinates are [lon, lat] in DEGREES.
    // Passing the vector straight through made d3 read [x, y] as degrees, so
    // every segment collapsed to a dot near the origin and the coastline simply
    // did not appear. 352 passing tests did not catch it — a screenshot did.
    const ds: string[] = [];
    for (const [a, b] of segments) {
      const d = this.path({ type: 'LineString', coordinates: [toLonLat(a), toLonLat(b)] });
      if (d) ds.push(d);
    }
    if (!ds.length) return;
    this.bodyParts.push(
      `<path d="${ds.join('')}" fill="none" stroke="${esc(stroke)}" stroke-width="${width}" ` +
      `stroke-linecap="round" stroke-linejoin="round"/>`,
    );
  }

  hatchRect(x: number, y: number, w: number, h: number, spec: HatchSpec): void {
    // One <pattern> per distinct spec — repeating a pattern definition per
    // rect would bloat the file for no visual difference.
    const id = this.hatchPatternId(spec);
    this.bodyParts.push(`<rect x="${x}" y="${y}" width="${w}" height="${h}" fill="url(#${id})"/>`);
  }

  /** One <pattern> per distinct spec, shared by hatchRect and hatchFeatures. */
  private hatchPatternId(spec: HatchSpec): string {
    const opacity = spec.opacity ?? 1;
    const key = `${spec.color}|${spec.spacingPx}|${spec.widthPx}|${spec.angleDeg}|${opacity}`;
    let id = this.patternIds.get(key);
    if (!id) {
      id = `hatch${this.patternIds.size}`;
      this.patternIds.set(key, id);
      // stroke-opacity, not an rgba() stroke: SVG 1.1 has no rgba() syntax.
      this.defsParts.push(
        `<pattern id="${id}" width="${spec.spacingPx}" height="${spec.spacingPx}" ` +
        `patternUnits="userSpaceOnUse" patternTransform="rotate(${spec.angleDeg})">` +
        `<line x1="0" y1="0" x2="0" y2="${spec.spacingPx}" ` +
        `stroke="${esc(spec.color)}" stroke-width="${spec.widthPx}"` +
        `${opacity < 1 ? ` stroke-opacity="${opacity}"` : ''}/>` +
        `</pattern>`,
      );
    }
    return id;
  }

  withSphereClip(draw: () => void): void {
    const d = this.path({ type: 'Sphere' });
    // No sphere path (a degenerate projection) means nothing to clip TO, not
    // "clip everything away" — draw unclipped rather than blanking the map.
    if (!d) { draw(); return; }
    if (!this.sphereClipId) {
      this.sphereClipId = 'sphereClip';
      this.defsParts.push(`<clipPath id="${this.sphereClipId}"><path d="${d}"/></clipPath>`);
    }
    this.bodyParts.push(`<g clip-path="url(#${this.sphereClipId})">`);
    draw();
    this.bodyParts.push('</g>');
  }

  grain(spec: GrainSpec): void {
    if (!this.grainId) {
      this.grainId = 'paperGrain';
      this.defsParts.push(
        `<filter id="${this.grainId}" x="0" y="0" width="100%" height="100%">` +
        `<feTurbulence type="fractalNoise" baseFrequency="${0.8 / Math.max(0.1, spec.scale)}" ` +
        `numOctaves="4" seed="${Math.abs(hashSeed(spec.seed)) % 1000}" result="n"/>` +
        `<feColorMatrix type="saturate" values="0" in="n" result="g"/>` +
        `</filter>`,
      );
    }
    this.bodyParts.push(
      `<rect x="0" y="0" width="${this.width}" height="${this.height}" ` +
      `filter="url(#${this.grainId})" opacity="${spec.opacity}"/>`,
    );
  }

  drawGlyph(g: PlacedGlyph, ink: string, width: number): void {
    // Mirror handling matches Canvas2DSubstrate: a nested flip composes with the
    // caller's outer flip to the identity, and the x-flip puts the glyph back
    // over its own cell. Without it every glyph renders backwards.
    const glyph = this.mirrored ? { ...g, x: this.width - g.x } : g;
    const path =
      `<path d="${glyphPathData(glyph)}" fill="none" stroke="${esc(ink)}" stroke-width="${width}" ` +
      `stroke-linecap="round" stroke-linejoin="round"/>`;
    this.bodyParts.push(
      this.mirrored
        ? `<g transform="translate(${this.width},0) scale(-1,1)">${path}</g>`
        : path,
    );
  }

  fillCells(indices: number[] | Uint32Array, colors: string[]): void {
    for (let j = 0; j < indices.length; j++) {
      const i = indices[j];
      const f = this.cellFeatures[i];
      if (!f) continue;
      this.fillFeature(f, colors[i]); // per-cell path, exactly as today
    }
  }
}

const hashSeed = (seed: string): number => {
  let h = 2166136261 >>> 0;
  for (let i = 0; i < seed.length; i++) {
    h ^= seed.charCodeAt(i);
    h = Math.imul(h, 16777619) >>> 0;
  }
  return h | 0;
};
