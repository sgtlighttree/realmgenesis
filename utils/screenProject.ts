// F2 screen-space overlay support: pure projection + occlusion math, no React,
// no THREE component types (kept unit-testable). The globe is a unit sphere, so
// a cell center's surface normal IS the point itself. A point is "visible" when
// the vector from the point to the camera has a positive-enough dot with that
// normal — i.e. the point faces the camera and is not past the horizon. The eps
// margin (0.08) absorbs the ≤1.05× relief grazing edge (spec §2.2).

export interface ProjectedCells {
  x: Float32Array; // screen-pixel x
  y: Float32Array; // screen-pixel y
  visible: Uint8Array; // 1 = on the near hemisphere, 0 = culled
  n: number;
}

export function isVisible(
  px: number, py: number, pz: number,
  camx: number, camy: number, camz: number,
  eps = 0.08,
): boolean {
  let dx = camx - px, dy = camy - py, dz = camz - pz;
  const len = Math.hypot(dx, dy, dz);
  if (len < 1e-9) return true;
  dx /= len; dy /= len; dz /= len;
  // surface normal at a unit-sphere point IS the point itself
  return dx * px + dy * py + dz * pz > eps;
}
