/**
 * On-screen stroke weights for the 2D vector overlay, in CSS pixels.
 *
 * `VectorOverlay.redraw` sets `setTransform(dpr * scale, ...)`, so everything
 * drawn is already multiplied by the zoom. A width therefore has to be divided
 * by that scale to hold a constant weight on screen — the same reason labels
 * pass `fontScale: 1 / s`.
 *
 * These used to be written inline as `Math.max(FLOOR, W / s)`. That floor is in
 * TRANSFORMED units while it was plainly meant as a screen minimum, so once
 * `W / s` dropped under it the width pinned at the floor and the on-screen
 * weight became `floor * s` — growing linearly with zoom instead of holding.
 * The coastline crossed over at s = 1.5 and reached ~6px by full zoom, as heavy
 * as the 12-14px labels beside it. Naming the widths in screen pixels makes the
 * unit explicit and leaves nowhere for a floor in the wrong unit to hide.
 */

/** Every stroke on the vector overlay, in CSS pixels at any zoom. */
export const SCREEN_STROKE_PX = {
  coast: 1.5,
  lake: 1.2,
  border: 2,
  contour: 1.5,
  graticule: 1,
  river: 1.5,
  road: 1.4,
  seaRoute: 1.2,
  cellOutline: 0.5,
  highlight: 3,
  markerPin: 2.5,
  ruler: 2,
} as const;

/** Sea-route dash pattern, in screen pixels. */
export const SCREEN_DASH_PX: readonly [number, number] = [5, 4];

/**
 * Convert a screen-pixel width into the units the scaled transform expects.
 * Guards a non-positive scale, which would otherwise divide by zero on the
 * first frame before the zoom state settles.
 */
export const screenStroke = (px: number, scale: number): number =>
  px / (scale > 0 ? scale : 1);
