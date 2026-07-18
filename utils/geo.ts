import * as d3 from 'd3';

import { DymaxionSettings } from '../types';
import { buildDymaxionNet, DymaxionNetFace } from './dymaxion';

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

export const dot3 = (a: Point3, b: Point3): number =>
  a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

export const sub3 = (a: Point3, b: Point3): Point3 =>
  [a[0] - b[0], a[1] - b[1], a[2] - b[2]];

export const cross3 = (a: Point3, b: Point3): Point3 => [
  a[1] * b[2] - a[2] * b[1],
  a[2] * b[0] - a[0] * b[2],
  a[0] * b[1] - a[1] * b[0],
];

export const barycentric3D = (p: Point3, a: Point3, b: Point3, c: Point3): [number, number, number] | null => {
  const v0 = sub3(b, a);
  const v1 = sub3(c, a);
  const v2 = sub3(p, a);
  const d00 = dot3(v0, v0);
  const d01 = dot3(v0, v1);
  const d11 = dot3(v1, v1);
  const d20 = dot3(v2, v0);
  const d21 = dot3(v2, v1);
  const denom = d00 * d11 - d01 * d01;
  if (!denom) return null;
  const v = (d11 * d20 - d01 * d21) / denom;
  const w = (d00 * d21 - d01 * d20) / denom;
  return [1 - v - w, v, w];
};

export const pointInsideSphericalFace = (p: Point3, vertices: Point3[]): boolean => {
  const centroid = normalizeVec([
    vertices[0][0] + vertices[1][0] + vertices[2][0],
    vertices[0][1] + vertices[1][1] + vertices[2][1],
    vertices[0][2] + vertices[1][2] + vertices[2][2],
  ]);
  for (let i = 0; i < 3; i++) {
    const a = vertices[i];
    const b = vertices[(i + 1) % 3];
    let edgeNormal = normalizeVec(cross3(a, b));
    if (dot3(edgeNormal, centroid) < 0) {
      edgeNormal = [-edgeNormal[0], -edgeNormal[1], -edgeNormal[2]];
    }
    if (dot3(edgeNormal, p) < -1e-7) return false;
  }
  return true;
};

export const getDymaxionNetTransform = (layout: DymaxionSettings['layout'], canvasWidth: number, canvasHeight: number) => {
  const net = buildDymaxionNet(layout);
  const faces = net.faces;
  let minX = Infinity;
  let minY = Infinity;
  let maxX = -Infinity;
  let maxY = -Infinity;
  faces.forEach((face) => {
    face.vertices.forEach((v) => {
      minX = Math.min(minX, v[0]);
      minY = Math.min(minY, v[1]);
      maxX = Math.max(maxX, v[0]);
      maxY = Math.max(maxY, v[1]);
    });
  });

  const pad = 8;
  const netWidth = Math.max(1e-6, maxX - minX);
  const netHeight = Math.max(1e-6, maxY - minY);
  const scale = Math.min((canvasWidth - pad * 2) / netWidth, (canvasHeight - pad * 2) / netHeight);
  const offsetX = (canvasWidth - netWidth * scale) / 2 - minX * scale;
  const offsetY = (canvasHeight - netHeight * scale) / 2 - minY * scale;

  return { faces, scale, offsetX, offsetY };
};

// Projects a unit-sphere position into Dymaxion NET space (the un-fitted 2D
// net coordinates). Callers apply their own net→canvas mapping so labels stay
// aligned with whichever raster fit (padding / UV flip) that pipeline uses.
export const projectToDymaxionNet = (
  position: Point3,
  faces: DymaxionNetFace[],
  lon: number,
  lat: number,
  roll: number,
): Point2 | null => {
  const rotate = d3.geoRotation([lon, lat, roll]);
  const [sourceLon, sourceLat] = toLonLat(position);
  const inverted = rotate.invert([-sourceLon, sourceLat]);
  if (!inverted) return null;
  const p3 = lonLatToPoint3(inverted as Point2);

  let selectedFace = faces[0];
  let selectedInside = false;
  let bestScore = -Infinity;

  for (const face of faces) {
    const faceCenter = normalizeVec([
      face.vertices3D[0][0] + face.vertices3D[1][0] + face.vertices3D[2][0],
      face.vertices3D[0][1] + face.vertices3D[1][1] + face.vertices3D[2][1],
      face.vertices3D[0][2] + face.vertices3D[1][2] + face.vertices3D[2][2],
    ]);
    const score = dot3(p3, faceCenter);
    const inside = pointInsideSphericalFace(p3, face.vertices3D);
    if (inside) {
      selectedFace = face;
      selectedInside = true;
      break;
    }
    if (!selectedInside && score > bestScore) {
      bestScore = score;
      selectedFace = face;
    }
  }

  const weights = barycentric3D(
    p3,
    selectedFace.vertices3D[0],
    selectedFace.vertices3D[1],
    selectedFace.vertices3D[2],
  );
  if (!weights) return null;

  const clamped = selectedInside ? weights : weights.map(w => Math.max(0, w)) as [number, number, number];
  const weightSum = clamped[0] + clamped[1] + clamped[2] || 1;
  const u = clamped[0] / weightSum;
  const v = clamped[1] / weightSum;
  const w = clamped[2] / weightSum;
  const x = selectedFace.vertices[0][0] * u + selectedFace.vertices[1][0] * v + selectedFace.vertices[2][0] * w;
  const y = selectedFace.vertices[0][1] * u + selectedFace.vertices[1][1] * v + selectedFace.vertices[2][1] * w;

  return [x, y];
};

export const projectDymaxionPoint = (
  position: Point3,
  layout: DymaxionSettings['layout'],
  canvasWidth: number,
  canvasHeight: number,
  lon: number,
  lat: number,
  roll: number,
): Point2 | null => {
  const { faces, scale, offsetX, offsetY } = getDymaxionNetTransform(layout, canvasWidth, canvasHeight);
  const net = projectToDymaxionNet(position, faces, lon, lat, roll);
  return net ? [net[0] * scale + offsetX, net[1] * scale + offsetY] : null;
};
