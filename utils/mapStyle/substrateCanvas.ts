import { Point3 } from '../geo';
import { glyphPathData } from './glyphPaths';
import { GeoFeatureLike, GrainSpec, HatchSpec, Substrate } from './substrate';
import { PlacedGlyph } from './types';

type PathGeneratorLike = (object: unknown) => unknown;

/** Deterministic 0-1 noise for the grain pass — no RNG import, no shared state. */
const hash01 = (x: number, y: number, seed: string): number => {
  let h = 2166136261 >>> 0;
  const str = `${seed}:${x}:${y}`;
  for (let i = 0; i < str.length; i++) {
    h ^= str.charCodeAt(i);
    h = Math.imul(h, 16777619) >>> 0;
  }
  return (h >>> 8) / 0x1000000;
};

export class Canvas2DSubstrate implements Substrate {
  /**
   * `mirrored` says the context already carries Map2D's horizontal flip
   * (`translate(width,0); scale(-1,1)`). Only glyphs care: polygons and strokes
   * come from the same mirrored path generator and land correctly, but a glyph
   * is drawn from its own coordinates and would come out backwards. See
   * `drawGlyph`.
   */
  constructor(
    private ctx: CanvasRenderingContext2D,
    private path: PathGeneratorLike,
    private width: number,
    private height: number,
    private mirrored = false,
  ) {}

  fillRect(x: number, y: number, w: number, h: number, fill: string): void {
    this.ctx.fillStyle = fill;
    this.ctx.fillRect(x, y, w, h);
  }

  fillFeature(feature: GeoFeatureLike, fill: string, opacity = 1): void {
    this.ctx.save();
    this.ctx.globalAlpha = opacity;
    this.ctx.beginPath();
    this.path(feature);
    this.ctx.fillStyle = fill;
    this.ctx.fill();
    this.ctx.restore();
  }

  hatchFeatures(features: GeoFeatureLike[], spec: HatchSpec): void {
    if (!features.length) return;
    this.ctx.save();
    // ONE composite path: d3's canvas path generator appends to the current
    // path, so calling it repeatedly without an intervening beginPath() unions
    // every feature into a single clip region. One clip, one line sweep.
    this.ctx.beginPath();
    for (const f of features) this.path(f);
    this.ctx.clip();
    this.hatchRect(0, 0, this.width, this.height, spec);
    this.ctx.restore();
  }

  strokeFeature(feature: GeoFeatureLike, stroke: string, width: number): void {
    this.ctx.beginPath();
    this.path(feature);
    this.ctx.strokeStyle = stroke;
    this.ctx.lineWidth = width;
    this.ctx.stroke();
  }

  strokeSegments(segments: Array<[Point3, Point3]>, stroke: string, width: number): void {
    if (!segments.length) return;
    this.ctx.strokeStyle = stroke;
    this.ctx.lineWidth = width;
    this.ctx.lineCap = 'round';
    this.ctx.lineJoin = 'round';
    this.ctx.beginPath();
    for (const [a, b] of segments) {
      this.path({ type: 'LineString', coordinates: [a, b] });
    }
    this.ctx.stroke();
  }

  hatchRect(x: number, y: number, w: number, h: number, spec: HatchSpec): void {
    this.ctx.save();
    this.ctx.beginPath();
    this.ctx.rect(x, y, w, h);
    this.ctx.clip();
    this.ctx.strokeStyle = spec.color;
    this.ctx.lineWidth = spec.widthPx;
    const rad = (spec.angleDeg * Math.PI) / 180;
    const dx = Math.cos(rad);
    const dy = Math.sin(rad);
    const span = Math.hypot(w, h);
    const steps = Math.ceil((span * 2) / Math.max(1, spec.spacingPx));
    this.ctx.beginPath();
    for (let i = -steps; i <= steps; i++) {
      const off = i * spec.spacingPx;
      const cx = x + w / 2 - dy * off;
      const cy = y + h / 2 + dx * off;
      this.ctx.moveTo(cx - dx * span, cy - dy * span);
      this.ctx.lineTo(cx + dx * span, cy + dy * span);
    }
    this.ctx.stroke();
    this.ctx.restore();
  }

  grain(spec: GrainSpec): void {
    this.ctx.save();
    this.ctx.globalAlpha = spec.opacity;
    const step = Math.max(2, Math.round(3 * spec.scale));
    for (let y = 0; y < this.height; y += step) {
      for (let x = 0; x < this.width; x += step) {
        const n = hash01(x, y, spec.seed);
        if (n < 0.55) continue;
        this.ctx.fillStyle = n > 0.85 ? '#ffffff' : '#000000';
        this.ctx.fillRect(x, y, step, step);
      }
    }
    this.ctx.restore();
  }

  drawGlyph(g: PlacedGlyph, ink: string, width: number): void {
    // Path2D accepts SVG path syntax, so glyphs use the SAME path string the
    // SVG substrate embeds. One shape definition, two surfaces.
    //
    // On a mirrored surface a glyph would render backwards — a mountain's
    // shading flick on the wrong side, and any future lettering unreadable.
    // Applying the SAME flip again composes to the identity, so the glyph draws
    // unmirrored; flipping its x puts it back over its own cell. Doing this
    // here keeps every pass free of surface geometry.
    this.ctx.save();
    let glyph = g;
    if (this.mirrored) {
      this.ctx.translate(this.width, 0);
      this.ctx.scale(-1, 1);
      glyph = { ...g, x: this.width - g.x };
    }
    const p = new Path2D(glyphPathData(glyph));
    this.ctx.strokeStyle = ink;
    this.ctx.lineWidth = width;
    this.ctx.lineCap = 'round';
    this.ctx.lineJoin = 'round';
    this.ctx.stroke(p);
    this.ctx.restore();
  }
}
