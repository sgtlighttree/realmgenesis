// F2 screen-space overlay support: pure projection + occlusion math, no React,
// no THREE component types (kept unit-testable). A point P is "visible" when the
// (unit) view ray from P to the camera has a positive dot with the outward
// surface normal at P — i.e. P faces the camera and is not past the horizon. The
// normal is P̂ = P/|P|, so this is the EXACT perspective horizon for a sphere of
// radius |P|. Callers MUST pass P at its true rendered radius (the globe mesh
// puts each cell at r = 1 + height·0.05, up to ≈1.05), NOT at unit radius:
// testing a rendered-radius point at r=1 culls a band inside the true limb and
// makes overlays drift off the terrain on zoom. eps is a hair (not the old fat
// 0.08) only to suppress limb flicker.

export interface ProjectedCells {
  x: Float32Array; // screen-pixel x
  y: Float32Array; // screen-pixel y
  visible: Uint8Array; // 1 = on the near hemisphere, 0 = culled
  n: number;
}

export function isVisible(
  px: number, py: number, pz: number,
  camx: number, camy: number, camz: number,
  eps = 0.005,
): boolean {
  let dx = camx - px, dy = camy - py, dz = camz - pz;
  // sqrt, not Math.hypot: hypot's overflow-safe scaling is ~11x slower in V8 and
  // buys nothing here — all inputs are near-unit sphere/camera coords (F4, S29).
  const len = Math.sqrt(dx * dx + dy * dy + dz * dz);
  if (len < 1e-9) return true;
  dx /= len; dy /= len; dz /= len;
  // outward normal at P is P̂ = P/|P| (P at its TRUE radius, not unit)
  const pl = Math.sqrt(px * px + py * py + pz * pz);
  if (pl < 1e-9) return true;
  return (dx * px + dy * py + dz * pz) / pl > eps;
}
