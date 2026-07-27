export type Point3 = [number, number, number];

// Leaf vector math. Split out of utils/geo.ts so utils/features.ts — and
// therefore the worker bundle — does not pull in all of d3.
export const normalizeVec = (v: Point3): Point3 => {
  const len = Math.hypot(v[0], v[1], v[2]) || 1;
  return [v[0] / len, v[1] / len, v[2] / len];
};
