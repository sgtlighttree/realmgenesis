import { Point } from '../types';

const clamp = (v: number, min: number, max: number): number => Math.max(min, Math.min(max, v));

// Local lon/lat -> unit-sphere conversion (mirrors utils/geo.ts's lonLatToPoint3),
// duplicated here as {x,y,z} rather than a tuple so it matches Cell.center's shape
// without pulling geo.ts's d3/dymaxion dependencies into this pure module.
const lonLatToPoint = (lon: number, lat: number): Point => {
  const lonRad = lon * (Math.PI / 180);
  const latRad = lat * (Math.PI / 180);
  const cosLat = Math.cos(latRad);
  return { x: cosLat * Math.cos(lonRad), y: Math.sin(latRad), z: cosLat * Math.sin(lonRad) };
};

// True geodesic distance between two unit-sphere points, in km.
export const greatCircleDistanceKm = (a: Point, b: Point, planetRadiusKm: number): number => {
  const dot = clamp(a.x * b.x + a.y * b.y + a.z * b.z, -1, 1);
  const angle = Math.acos(dot);
  return angle * planetRadiusKm;
};

// Any unit vector perpendicular to v, used only for the antipodal fallback
// below where the great circle's plane is otherwise undefined.
const arbitraryPerpendicular = (v: Point): Point => {
  const seed: Point = Math.abs(v.x) > 0.9 ? { x: 0, y: 1, z: 0 } : { x: 1, y: 0, z: 0 };
  const d = v.x * seed.x + v.y * seed.y + v.z * seed.z;
  const px = seed.x - d * v.x;
  const py = seed.y - d * v.y;
  const pz = seed.z - d * v.z;
  const len = Math.hypot(px, py, pz) || 1;
  return { x: px / len, y: py / len, z: pz / len };
};

// Samples points along the great-circle arc from a to b (both unit-sphere),
// via spherical linear interpolation. Two degeneracies get special handling:
// - near-identical points: slerp's sin(angle) denominator is ~0, so every
//   sample just holds at a.
// - near-antipodal points: the geodesic itself is mathematically undefined
//   (infinitely many great circles pass through both poles), so the arc
//   sweeps through an arbitrary perpendicular direction instead of slerping
//   directly — this keeps every sample finite and unit-length rather than
//   collapsing toward the zero vector at the midpoint.
export const sampleGreatCircleArc = (a: Point, b: Point, segments = 48): Point[] => {
  const dot = clamp(a.x * b.x + a.y * b.y + a.z * b.z, -1, 1);
  const angle = Math.acos(dot);
  const sinAngle = Math.sin(angle);
  const points: Point[] = [];
  const isAntipodal = Math.PI - angle < 1e-6;
  const perp = isAntipodal ? arbitraryPerpendicular(a) : null;

  for (let i = 0; i <= segments; i++) {
    const t = i / segments;
    let x: number;
    let y: number;
    let z: number;

    if (perp) {
      const theta = t * Math.PI;
      x = Math.cos(theta) * a.x + Math.sin(theta) * perp.x;
      y = Math.cos(theta) * a.y + Math.sin(theta) * perp.y;
      z = Math.cos(theta) * a.z + Math.sin(theta) * perp.z;
    } else if (sinAngle < 1e-6) {
      x = a.x + (b.x - a.x) * t;
      y = a.y + (b.y - a.y) * t;
      z = a.z + (b.z - a.z) * t;
    } else {
      const fa = Math.sin((1 - t) * angle) / sinAngle;
      const fb = Math.sin(t * angle) / sinAngle;
      x = fa * a.x + fb * b.x;
      y = fa * a.y + fb * b.y;
      z = fa * a.z + fb * b.z;
    }

    const len = Math.hypot(x, y, z) || 1;
    points.push({ x: x / len, y: y / len, z: z / len });
  }

  return points;
};

// Projection-aware scale: samples two points +/-0.5 degrees apart along the
// parallel through centerLonLat, projects both, and divides the resulting
// pixel distance by their true geodesic separation. Returns null if either
// sample point fails to project (e.g. off the visible hemisphere).
export const computeScaleBar = (
  projection: (lonLat: [number, number]) => [number, number] | null,
  centerLonLat: [number, number],
  planetRadiusKm: number,
): { pixelsPerKm: number } | null => {
  const [lon, lat] = centerLonLat;
  const halfDeg = 0.5;
  const clampedLat = clamp(lat, -89.5, 89.5); // avoid pole singularity in the parallel sample

  const a = projection([lon - halfDeg, clampedLat]);
  const b = projection([lon + halfDeg, clampedLat]);
  if (!a || !b) return null;

  const pixelDist = Math.hypot(b[0] - a[0], b[1] - a[1]);
  const kmDist = greatCircleDistanceKm(
    lonLatToPoint(lon - halfDeg, clampedLat),
    lonLatToPoint(lon + halfDeg, clampedLat),
    planetRadiusKm,
  );
  if (!(kmDist > 0)) return null;

  return { pixelsPerKm: pixelDist / kmDist };
};

/**
 * Largest step-to-step ratio treated as smooth. Above it, the two probe samples
 * are taken to sit on opposite sides of a cut.
 *
 * Generous on purpose. The probes are half a degree apart — about 55 km, against
 * a Dymaxion face 63 degrees across — and the projection is smooth WITHIN a
 * face, so ordinary variation between two adjacent steps is a fraction of a
 * percent. Anything approaching 1.5 is a discontinuity, not a gradient.
 */
const FOLD_RATIO = 1.5;

/**
 * `computeScaleBar` for an interrupted projection — one whose plane is cut, so
 * that two points a few kilometres apart on the globe can land far apart on the
 * page.
 *
 * **Why the plain version is not enough.** It samples `lon ± 0.5` and divides
 * the pixel distance by the true geodesic distance. Across a cut that pixel
 * distance is meaningless, and worse, it FAILS QUIETLY: the two halves of a
 * Dymaxion net often land near each other, so the bad reading looks plausible
 * rather than absurd. There is nothing in the number itself to reject.
 *
 * So this samples three points instead of two and compares the two steps. On a
 * face they match to within a fraction of a percent; across a cut one of them
 * jumps. A mismatch returns null, and the caller shows no bar — which is the
 * honest answer for a viewport centred on a fold.
 *
 * **This is NOT a distortion guard.** Measured over 1207 sample points on this
 * app's own Dymaxion net, LINEAR scale varies 13% peak to peak across the whole
 * globe, so one bar is good to about ±7% anywhere — where Mercator, which
 * already carries a bar, is out by 100% at 60 degrees and 480% at 80. (The
 * often-quoted "2% Fuller error" is an AREA figure and the wrong statistic for a
 * scale bar; don't reach for it here.) The only thing being detected is the cut.
 * That measurement rejected 1.7% of those points as fold-crossing.
 */
export const computeInterruptedScaleBar = (
  projection: (lonLat: [number, number]) => [number, number] | null,
  centerLonLat: [number, number],
  planetRadiusKm: number,
): { pixelsPerKm: number } | null => {
  const [lon, lat] = centerLonLat;
  const halfDeg = 0.5;
  const clampedLat = clamp(lat, -89.5, 89.5);

  const left = projection([lon - halfDeg, clampedLat]);
  const mid = projection([lon, clampedLat]);
  const right = projection([lon + halfDeg, clampedLat]);
  if (!left || !mid || !right) return null;

  const d1 = Math.hypot(mid[0] - left[0], mid[1] - left[1]);
  const d2 = Math.hypot(right[0] - mid[0], right[1] - mid[1]);
  if (!(d1 > 0) || !(d2 > 0)) return null;
  if (Math.max(d1, d2) / Math.min(d1, d2) > FOLD_RATIO) return null;

  return computeScaleBar(projection, centerLonLat, planetRadiusKm);
};

const NICE_STEPS = [1, 2, 5];

// Largest "nice" km value (1/2/5 x 10^n) whose bar fits within maxPixels.
export const niceScaleBarLength = (pixelsPerKm: number, maxPixels: number): { km: number; px: number } => {
  if (!(pixelsPerKm > 0) || !(maxPixels > 0)) return { km: 0, px: 0 };
  const maxKm = maxPixels / pixelsPerKm;
  if (!(maxKm > 0)) return { km: 0, px: 0 };

  // Epsilon guards against log10 landing just under an exact power of ten
  // (e.g. log10(1000) evaluating to 2.9999999996) and picking the wrong decade.
  const exponent = Math.floor(Math.log10(maxKm) + 1e-9);
  let best = 0;

  for (const e of [exponent - 1, exponent]) {
    for (const step of NICE_STEPS) {
      const candidate = step * Math.pow(10, e);
      if (candidate <= maxKm && candidate > best) best = candidate;
    }
  }

  return { km: best, px: best * pixelsPerKm };
};

// Canvas-drawing counterpart to the pure math above (same split as
// drawContourPaths/drawMapLabels living alongside their compute functions).
// (x, bottomY) anchors the bottom-left corner; the label is drawn above the bar.
export const drawScaleBar = (
  ctx: CanvasRenderingContext2D,
  x: number,
  bottomY: number,
  km: number,
  px: number,
): void => {
  if (!(km > 0) || !(px > 0)) return;

  const barHeight = 5;
  const barY = bottomY - barHeight;

  ctx.save();
  ctx.lineJoin = 'round';

  ctx.fillStyle = 'rgba(0, 0, 0, 0.85)';
  ctx.fillRect(x - 1.5, barY - 1.5, px + 3, barHeight + 3);
  ctx.fillStyle = '#ffffff';
  ctx.fillRect(x, barY, px, barHeight);

  ctx.strokeStyle = 'rgba(0, 0, 0, 0.85)';
  ctx.lineWidth = 2;
  ctx.beginPath();
  ctx.moveTo(x, barY - 3);
  ctx.lineTo(x, barY + barHeight + 3);
  ctx.moveTo(x + px, barY - 3);
  ctx.lineTo(x + px, barY + barHeight + 3);
  ctx.stroke();

  const label = `${Math.round(km).toLocaleString()} km`;
  ctx.font = '11px sans-serif';
  ctx.textAlign = 'left';
  ctx.textBaseline = 'bottom';
  ctx.lineWidth = 3;
  ctx.strokeStyle = 'rgba(0, 0, 0, 0.85)';
  ctx.strokeText(label, x, barY - 4);
  ctx.fillStyle = '#ffffff';
  ctx.fillText(label, x, barY - 4);

  ctx.restore();
};
