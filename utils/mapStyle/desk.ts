import { StylePalette } from './types';

/**
 * The desk — what fills the viewport outside the map.
 *
 * SCREEN SPACE, always. Map2D blits its offscreen bitmap under the pan/zoom
 * transform, so anything painted into that bitmap slides and shrinks with the
 * paper instead of filling the viewport. The desk is furniture, like the scale
 * bar: it belongs to the viewport, not to the map. (It was a style pass for
 * exactly one round; this is why it is not one.)
 *
 * Every style paints a desk, including `default`. A projection rarely fills its
 * canvas — Mercator clipped at +/-85 degrees is square — and raw black in that
 * margin reads as a rendering failure rather than a choice.
 */

/** `#rrggbb` (or `#rgb`) with an alpha, for gradient stops. Non-hex passes through. */
const withAlpha = (hex: string, alpha: number): string => {
  const m = /^#([0-9a-f]{3}|[0-9a-f]{6})$/i.exec(hex.trim());
  if (!m) return hex;
  const h = m[1].length === 3 ? m[1].replace(/./g, (c) => c + c) : m[1];
  const n = parseInt(h, 16);
  return `rgba(${(n >> 16) & 255}, ${(n >> 8) & 255}, ${n & 255}, ${alpha})`;
};

/**
 * Paint the desk over the whole canvas, in device pixels, with the identity
 * transform already set by the caller.
 *
 * Three layers, cheapest first:
 *   1. `desk` — the flat surface colour.
 *   2. `deskTexture` — an optional repeating image (wood, leather, cloth).
 *   3. `deskEdge` — an optional vignette, transparent at the centre and solid
 *      at the corners, so the surface reads as lit rather than as a flat
 *      swatch. Fixed to the viewport: a desk lamp does not pan with the paper.
 *
 * A style that wants a plain fill simply omits the last two.
 */
export const paintDesk = (
  ctx: CanvasRenderingContext2D,
  palette: StylePalette,
  width: number,
  height: number,
  texture?: CanvasImageSource | null,
): void => {
  ctx.fillStyle = palette.desk;
  ctx.fillRect(0, 0, width, height);

  if (texture) {
    const pattern = ctx.createPattern(texture, 'repeat');
    if (pattern) {
      ctx.fillStyle = pattern;
      ctx.fillRect(0, 0, width, height);
    }
  }

  const edge = palette.deskEdge;
  if (!edge) return;
  const cx = width / 2;
  const cy = height / 2;
  const outer = Math.hypot(cx, cy);
  const grad = ctx.createRadialGradient(cx, cy, 0, cx, cy, outer);
  // Deliberately late and shallow. The desk is seen almost entirely in the
  // LETTERBOX MARGIN down each side, and a viewport is wide: the mid-edge sits
  // at ~0.87 of the corner radius, so a gradient that starts early and reaches
  // full strength paints the margin near-black — which is the raw black this
  // whole feature exists to replace. Nothing darkens inside 0.7, and the
  // corners keep just over half the base colour.
  grad.addColorStop(0.7, withAlpha(edge, 0));
  grad.addColorStop(1, withAlpha(edge, 0.55));
  ctx.fillStyle = grad;
  ctx.fillRect(0, 0, width, height);
};
