import { geoGnomonic } from 'd3';
import { geoPolyhedral } from 'd3-geo-projection';
import type { GeoProjection } from 'd3';
import type { DymaxionLayout } from '../types';

type Vec3 = [number, number, number];
type Face = [number, number, number];

const PHI = (1 + Math.sqrt(5)) / 2;

const RAW_VERTS: Vec3[] = [
  [-1, PHI, 0],
  [1, PHI, 0],
  [-1, -PHI, 0],
  [1, -PHI, 0],
  [0, -1, PHI],
  [0, 1, PHI],
  [0, -1, -PHI],
  [0, 1, -PHI],
  [PHI, 0, -1],
  [PHI, 0, 1],
  [-PHI, 0, -1],
  [-PHI, 0, 1],
];

const FACES: Face[] = [
  [0, 11, 5],
  [0, 5, 1],
  [0, 1, 7],
  [0, 7, 10],
  [0, 10, 11],
  [1, 5, 9],
  [5, 11, 4],
  [11, 10, 2],
  [10, 7, 6],
  [7, 1, 8],
  [3, 9, 4],
  [3, 4, 2],
  [3, 2, 6],
  [3, 6, 8],
  [3, 8, 9],
  [4, 9, 5],
  [2, 4, 11],
  [6, 2, 10],
  [8, 6, 7],
  [9, 8, 1],
];

const normalize = ([x, y, z]: Vec3): Vec3 => {
  const len = Math.hypot(x, y, z) || 1;
  return [x / len, y / len, z / len];
};

const add = (a: Vec3, b: Vec3): Vec3 => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const sub = (a: Vec3, b: Vec3): Vec3 => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const cross = (a: Vec3, b: Vec3): Vec3 => [
  a[1] * b[2] - a[2] * b[1],
  a[2] * b[0] - a[0] * b[2],
  a[0] * b[1] - a[1] * b[0],
];
const dot = (a: Vec3, b: Vec3): number => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

const cartesianToLonLat = ([x, y, z]: Vec3): [number, number] => {
  const lat = Math.asin(Math.max(-1, Math.min(1, y))) * (180 / Math.PI);
  const lon = Math.atan2(z, x) * (180 / Math.PI);
  return [lon, lat];
};

const sphericalToCartesian = (lambda: number, phi: number): Vec3 => {
  const cosphi = Math.cos(phi);
  return [cosphi * Math.cos(lambda), Math.sin(phi), cosphi * Math.sin(lambda)];
};

const buildFaces = () => {
  const verts = RAW_VERTS.map(normalize);
  const vertsLonLat = verts.map(cartesianToLonLat);

  // geoPolyhedral expects clockwise vertex order (as seen from outside).
  // Flip faces whose normals point outward so the winding is clockwise.
  const orientedFacesIdx: Face[] = FACES.map((face) => {
    const [a, b, c] = face;
    const normal = cross(sub(verts[b], verts[a]), sub(verts[c], verts[a]));
    const centroid = normalize(add(add(verts[a], verts[b]), verts[c]));
    return dot(normal, centroid) > 0 ? [a, c, b] : face;
  });

  const facesCart = orientedFacesIdx.map((face) => face.map((idx) => verts[idx]) as Vec3[]);
  const facesLonLat = orientedFacesIdx.map((face) => face.map((idx) => vertsLonLat[idx]) as [number, number][]);
  const normals = facesCart.map((face) => normalize(add(add(face[0], face[1]), face[2])));
  const edgeNormals = facesCart.map((face, fi) => {
    const centroid = normals[fi];
    return face.map((v, i) => {
      const next = face[(i + 1) % face.length];
      let n = cross(v, next);
      if (dot(n, centroid) < 0) n = [-n[0], -n[1], -n[2]];
      return normalize(n);
    });
  });

  return { facesLonLat, normals, edgeNormals, facesIdx: orientedFacesIdx, vertsLonLat, facesCart, verts };
};

export type DymaxionNetFace = {
  index: number;
  vertices: [number, number][];
  localVertices: [number, number][];
  vertices3D: Vec3[];
  transform: [number, number, number, number, number, number];
  inverse: [number, number, number, number, number, number];
};

const inverseMatrix = (m: [number, number, number, number, number, number]) => {
  const det = m[0] * m[4] - m[1] * m[3];
  const invDet = det ? 1 / det : 0;
  return [
    m[4] * invDet,
    -m[1] * invDet,
    (m[1] * m[5] - m[2] * m[4]) * invDet,
    -m[3] * invDet,
    m[0] * invDet,
    (m[2] * m[3] - m[0] * m[5]) * invDet
  ] as [number, number, number, number, number, number];
};

const otherSidePoint = (a: [number, number], b: [number, number], p: [number, number]) => {
  const dx = b[0] - a[0];
  const dy = b[1] - a[1];
  const len = Math.hypot(dx, dy) || 1;
  const h = (Math.sqrt(3) / 2) * len;
  const perp: [number, number] = [-dy / len, dx / len];
  const mid: [number, number] = [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2];
  const c1: [number, number] = [mid[0] + perp[0] * h, mid[1] + perp[1] * h];
  const c2: [number, number] = [mid[0] - perp[0] * h, mid[1] - perp[1] * h];
  const side = (q: [number, number]) => dx * (q[1] - a[1]) - dy * (q[0] - a[0]);
  return side(p) >= 0 ? (side(c1) >= 0 ? c2 : c1) : (side(c1) <= 0 ? c2 : c1);
};

const findSharedEdge = (a: Face, b: Face) => {
  for (let i = 0; i < 3; i++) {
    const a0 = a[i];
    const a1 = a[(i + 1) % 3];
    for (let j = 0; j < 3; j++) {
      const b0 = b[j];
      const b1 = b[(j + 1) % 3];
      if ((a0 === b0 && a1 === b1) || (a0 === b1 && a1 === b0)) {
        return { a0, a1 };
      }
    }
  }
  return null;
};

export const buildDymaxionNet = (layout: DymaxionLayout) => {
  const { facesLonLat, facesIdx, facesCart } = buildFaces();
  const { root, parents } = buildParents(layout, buildAdjacency(FACES));

  const localVertices = facesCart.map((face) => {
    const [a, b, c] = face;
    const xAxis = normalize(sub(b, a));
    const faceNormal = normalize(cross(sub(b, a), sub(c, a)));
    const yAxis = cross(faceNormal, xAxis);
    return face.map((v) => [dot(sub(v, a), xAxis), dot(sub(v, a), yAxis)] as [number, number]);
  });

  const children: number[][] = Array.from({ length: facesLonLat.length }, () => []);
  parents.forEach((p, i) => {
    if (p !== -1) children[p].push(i);
  });

  const faceCoords = new Array(facesLonLat.length).fill(null) as (Record<number, [number, number]> | null)[];
  const rootIds = facesIdx[root];
  const rootLocal = localVertices[root];
  faceCoords[root] = {
    [rootIds[0]]: rootLocal[0],
    [rootIds[1]]: rootLocal[1],
    [rootIds[2]]: rootLocal[2]
  };

  const order: number[] = [root];
  while (order.length) {
    const current = order.shift()!;
    const currentCoords = faceCoords[current];
    if (!currentCoords) continue;
    children[current].forEach((child) => {
      if (faceCoords[child]) return;
      const shared = findSharedEdge(facesIdx[current], facesIdx[child]);
      if (!shared) return;
      const idsParent = facesIdx[current];
      const idsChild = facesIdx[child];
      const parentThird = idsParent.find((id) => id !== shared.a0 && id !== shared.a1) as number;
      const childThird = idsChild.find((id) => id !== shared.a0 && id !== shared.a1) as number;
      const A = currentCoords[shared.a0];
      const B = currentCoords[shared.a1];
      const P = currentCoords[parentThird];
      const C = otherSidePoint(A, B, P);
      faceCoords[child] = {
        [shared.a0]: A,
        [shared.a1]: B,
        [childThird]: C
      };
      order.push(child);
    });
  }

  const netFaces: DymaxionNetFace[] = facesLonLat.map((face, i) => {
    const coords = faceCoords[i];
    const ids = facesIdx[i];
    const vertices = coords
      ? (ids.map((id) => coords[id]) as [number, number][])
      : localVertices[i];
    const transform = [1, 0, 0, 0, 1, 0] as [number, number, number, number, number, number];
    return {
      index: i,
      vertices,
      localVertices: localVertices[i],
      vertices3D: facesCart[i],
      transform,
      inverse: inverseMatrix(transform)
    };
  });

  return { faces: netFaces };
};

const buildAdjacency = (faces: Face[]) => {
  const edgeMap = new Map<string, number[]>();
  const adjacency: number[][] = Array.from({ length: faces.length }, () => []);
  const edgeKey = (a: number, b: number) => (a < b ? `${a}-${b}` : `${b}-${a}`);

  faces.forEach((face, fi) => {
    for (let i = 0; i < 3; i++) {
      const a = face[i];
      const b = face[(i + 1) % 3];
      const key = edgeKey(a, b);
      const list = edgeMap.get(key) || [];
      list.push(fi);
      edgeMap.set(key, list);
    }
  });

  edgeMap.forEach((facesIdx) => {
    if (facesIdx.length === 2) {
      const [a, b] = facesIdx;
      adjacency[a].push(b);
      adjacency[b].push(a);
    }
  });

  return adjacency;
};

const buildParents = (layout: DymaxionLayout, adjacency: number[][]) => {
  // Explicit spanning trees for known Dymaxion-style nets.
  const parents = new Array(adjacency.length).fill(-1);

  if (layout === 'classic') {
    const root = 5;
    parents[15] = 5;
    parents[6] = 15;
    parents[16] = 6;
    parents[7] = 16;
    parents[17] = 7;
    parents[8] = 17;
    parents[18] = 8;
    parents[9] = 18;
    parents[19] = 9;

    parents[1] = 5;
    parents[0] = 6;
    parents[4] = 7;
    parents[3] = 8;
    parents[2] = 9;

    parents[10] = 15;
    parents[11] = 16;
    parents[12] = 17;
    parents[13] = 18;
    parents[14] = 19;

    return { root, parents };
  }

  const root = 6;
  parents[16] = 6;
  parents[7] = 16;
  parents[17] = 7;
  parents[8] = 17;
  parents[18] = 8;
  parents[9] = 18;
  parents[19] = 9;
  parents[5] = 19;
  parents[15] = 5;

  parents[0] = 6;
  parents[4] = 7;
  parents[3] = 8;
  parents[2] = 9;
  parents[1] = 5;

  parents[11] = 16;
  parents[12] = 17;
  parents[13] = 18;
  parents[14] = 19;
  parents[10] = 15;

  return { root, parents };
};

export const createDymaxionProjection = (layout: DymaxionLayout): GeoProjection => {
  const { facesLonLat, normals, edgeNormals } = buildFaces();
  const adjacency = buildAdjacency(FACES);
  const { root, parents } = buildParents(layout, adjacency);

  const nodes = facesLonLat.map((face, i) => {
    const c = cartesianToLonLat(normals[i]);
    return {
      face,
      project: geoGnomonic().scale(1).translate([0, 0]).rotate([-c[0], -c[1]]),
    };
  });

  parents.forEach((parent, i) => {
    if (parent === -1) return;
    const parentNode = nodes[parent];
    parentNode.children ||= [];
    parentNode.children.push(nodes[i]);
  });

  const faceSelector = (lambda: number, phi: number) => {
    const p = sphericalToCartesian(lambda, phi);
    const eps = 1e-9;
    for (let i = 0; i < edgeNormals.length; i++) {
      const edges = edgeNormals[i];
      if (dot(edges[0], p) >= -eps && dot(edges[1], p) >= -eps && dot(edges[2], p) >= -eps) {
        return nodes[i];
      }
    }
    let bestIndex = 0;
    let bestDot = -Infinity;
    for (let i = 0; i < normals.length; i++) {
      const d = dot(normals[i], p);
      if (d > bestDot) {
        bestDot = d;
        bestIndex = i;
      }
    }
    return nodes[bestIndex];
  };

  return geoPolyhedral(nodes[root], faceSelector);
};
