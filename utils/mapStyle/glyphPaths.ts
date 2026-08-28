import { GlyphKind, PlacedGlyph } from './types';

/**
 * Procedural glyph shapes as SVG path data, in output pixel coordinates.
 *
 * One path string serves both substrates: `SvgSubstrate` embeds it directly and
 * `Canvas2DSubstrate` feeds it to `new Path2D(d)`, which accepts SVG path
 * syntax. That is what guarantees the raster and vector maps draw the SAME
 * shape rather than two drifting approximations.
 *
 * Shapes are defined in a unit box — x and y in [-1, 1], baseline at y = 0,
 * apex toward -y — then rotated by `seedRot`, scaled and translated.
 */

type UnitPoint = [number, number];
type UnitSubpath = UnitPoint[];

const MOUNTAIN: UnitSubpath[] = [
  [[-1, 0], [-0.35, -1], [0, -0.45], [0.35, -1], [1, 0]],
  [[-0.35, -1], [-0.12, -0.62], [0.06, -0.78]], // snow/shading flick
];

const HILL: UnitSubpath[] = [
  [[-1, 0], [-0.55, -0.55], [0, -0.7], [0.55, -0.55], [1, 0]],
];

const FOREST: UnitSubpath[] = [
  [[-0.55, 0], [-0.55, -0.4]],
  [[-0.9, -0.4], [-0.55, -0.95], [-0.2, -0.4]],
  [[0.55, 0], [0.55, -0.35]],
  [[0.2, -0.35], [0.55, -0.85], [0.9, -0.35]],
];

const CONIFER: UnitSubpath[] = [
  [[0, 0], [0, -0.35]],
  [[-0.6, -0.3], [0, -1], [0.6, -0.3]],
  [[-0.42, -0.62], [0, -1], [0.42, -0.62]],
];

const DUNE: UnitSubpath[] = [
  [[-1, 0], [-0.4, -0.42], [0.35, -0.15], [1, -0.32]],
  [[-0.6, 0.28], [0.1, 0.05], [0.85, 0.22]],
];

const MARSH: UnitSubpath[] = [
  [[-0.8, 0], [0.8, 0]],
  [[-0.45, -0.1], [-0.45, -0.75]],
  [[0.05, -0.1], [0.05, -0.9]],
  [[0.55, -0.1], [0.55, -0.7]],
];

const SHAPES: Record<GlyphKind, UnitSubpath[]> = {
  mountain: MOUNTAIN,
  hill: HILL,
  forest: FOREST,
  conifer: CONIFER,
  dune: DUNE,
  marsh: MARSH,
};

const round = (n: number): string => (Math.round(n * 100) / 100).toString();

/**
 * SVG path `d` for a placed glyph, already positioned in output pixels.
 * Open polylines — these are drawn with `strokePath`, never filled.
 */
export const glyphPathData = (g: PlacedGlyph): string => {
  const cos = Math.cos(g.seedRot);
  const sin = Math.sin(g.seedRot);
  const half = g.scale / 2;

  const parts: string[] = [];
  for (const subpath of SHAPES[g.kind]) {
    const points = subpath.map(([ux, uy]) => {
      const sx = ux * half;
      const sy = uy * half;
      return [
        g.x + sx * cos - sy * sin,
        g.y + sx * sin + sy * cos,
      ] as UnitPoint;
    });
    parts.push(
      `M${round(points[0][0])} ${round(points[0][1])}` +
      points.slice(1).map(p => `L${round(p[0])} ${round(p[1])}`).join(''),
    );
  }
  return parts.join('');
};
