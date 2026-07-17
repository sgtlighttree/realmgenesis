// Shared 2D/3D geometry helpers used by the Dymaxion raster pipelines in
// utils/export.ts and components/Map2D.tsx. Keep these in sync-free: both
// pipelines must rasterize identically or picking drifts from the visible map.

export type Point2 = [number, number];
export type Point3 = [number, number, number];

export const insideTri = (p: Point2, a: Point2, b: Point2, c: Point2): boolean => {
  const v0 = [c[0] - a[0], c[1] - a[1]];
  const v1 = [b[0] - a[0], b[1] - a[1]];
  const v2 = [p[0] - a[0], p[1] - a[1]];
  const dot00 = v0[0] * v0[0] + v0[1] * v0[1];
  const dot01 = v0[0] * v1[0] + v0[1] * v1[1];
  const dot02 = v0[0] * v2[0] + v0[1] * v2[1];
  const dot11 = v1[0] * v1[0] + v1[1] * v1[1];
  const dot12 = v1[0] * v2[0] + v1[1] * v2[1];
  const invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
  const u = (dot11 * dot02 - dot01 * dot12) * invDenom;
  const v = (dot00 * dot12 - dot01 * dot02) * invDenom;
  return u >= -1e-6 && v >= -1e-6 && u + v <= 1 + 1e-6;
};

export const barycentric = (p: Point2, a: Point2, b: Point2, c: Point2): [number, number, number] | null => {
  const v0 = [b[0] - a[0], b[1] - a[1]];
  const v1 = [c[0] - a[0], c[1] - a[1]];
  const v2 = [p[0] - a[0], p[1] - a[1]];
  const d00 = v0[0] * v0[0] + v0[1] * v0[1];
  const d01 = v0[0] * v1[0] + v0[1] * v1[1];
  const d11 = v1[0] * v1[0] + v1[1] * v1[1];
  const d20 = v2[0] * v0[0] + v2[1] * v0[1];
  const d21 = v2[0] * v1[0] + v2[1] * v1[1];
  const denom = d00 * d11 - d01 * d01;
  if (!denom) return null;
  const v = (d11 * d20 - d01 * d21) / denom;
  const w = (d00 * d21 - d01 * d20) / denom;
  const u = 1 - v - w;
  return [u, v, w];
};

export const normalizeVec = (v: Point3): Point3 => {
  const len = Math.hypot(v[0], v[1], v[2]) || 1;
  return [v[0] / len, v[1] / len, v[2] / len];
};

export const toLonLat = (v: Point3): Point2 => {
  const lon = Math.atan2(v[2], v[0]) * (180 / Math.PI);
  const lat = Math.asin(Math.max(-1, Math.min(1, v[1]))) * (180 / Math.PI);
  return [lon, lat];
};

export const lonLatToPoint3 = ([lon, lat]: Point2): Point3 => {
  const lonRad = lon * (Math.PI / 180);
  const latRad = lat * (Math.PI / 180);
  const cosLat = Math.cos(latRad);
  return [cosLat * Math.cos(lonRad), Math.sin(latRad), cosLat * Math.sin(lonRad)];
};
