import { Point3 } from '../geo';
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

  constructor(
    private path: PathStringGenerator,
    private width: number,
    private height: number,
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
    const op = opacity < 1 ? ` fill-opacity="${opacity.toFixed(3)}"` : '';
    this.bodyParts.push(
      `<path d="${d}" fill="${esc(fill)}"${op} stroke="${esc(fill)}" stroke-width="0.5"/>`,
    );
  }

  hatchFeature(feature: GeoFeatureLike, spec: HatchSpec): void {
    const d = this.path(feature);
    if (!d) return;
    const id = this.hatchPatternId(spec);
    this.bodyParts.push(`<path d="${d}" fill="url(#${id})"/>`);
  }

  strokeFeature(feature: GeoFeatureLike, stroke: string, width: number): void {
    const d = this.path(feature);
    if (!d) return;
    this.bodyParts.push(`<path d="${d}" fill="none" stroke="${esc(stroke)}" stroke-width="${width}"/>`);
  }

  strokeSegments(segments: Array<[Point3, Point3]>, stroke: string, width: number): void {
    const ds: string[] = [];
    for (const [a, b] of segments) {
      const d = this.path({ type: 'LineString', coordinates: [a, b] });
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

  /** One <pattern> per distinct spec, shared by hatchRect and hatchFeature. */
  private hatchPatternId(spec: HatchSpec): string {
    const key = `${spec.color}|${spec.spacingPx}|${spec.widthPx}|${spec.angleDeg}`;
    let id = this.patternIds.get(key);
    if (!id) {
      id = `hatch${this.patternIds.size}`;
      this.patternIds.set(key, id);
      this.defsParts.push(
        `<pattern id="${id}" width="${spec.spacingPx}" height="${spec.spacingPx}" ` +
        `patternUnits="userSpaceOnUse" patternTransform="rotate(${spec.angleDeg})">` +
        `<line x1="0" y1="0" x2="0" y2="${spec.spacingPx}" ` +
        `stroke="${esc(spec.color)}" stroke-width="${spec.widthPx}"/>` +
        `</pattern>`,
      );
    }
    return id;
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
    this.bodyParts.push(
      `<path d="${glyphPathData(g)}" fill="none" stroke="${esc(ink)}" stroke-width="${width}" ` +
      `stroke-linecap="round" stroke-linejoin="round"/>`,
    );
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
