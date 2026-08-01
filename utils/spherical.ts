import { Point } from '../types';
import { RNG } from './rng';

// Quaternion rotation on the unit sphere
export function quatFromAxisAngle(axis: Point, angle: number): [number, number, number, number] {
  const half = angle / 2;
  const s = Math.sin(half);
  const len = Math.sqrt(axis.x * axis.x + axis.y * axis.y + axis.z * axis.z);
  if (len === 0) return [1, 0, 0, 0];
  const nx = axis.x / len, ny = axis.y / len, nz = axis.z / len;
  return [Math.cos(half), nx * s, ny * s, nz * s];
}

export function quatRotate(q: [number, number, number, number], v: Point): Point {
  const [w, x, y, z] = q;
  const vx = v.x, vy = v.y, vz = v.z;
  const tx = 2 * (y * vz - z * vy);
  const ty = 2 * (z * vx - x * vz);
  const tz = 2 * (x * vy - y * vx);
  return {
    x: vx + w * tx + (y * tz - z * ty),
    y: vy + w * ty + (z * tx - x * tz),
    z: vz + w * tz + (x * ty - y * tx),
  };
}

// Build an Euler pole: a unit axis + angular rate (radians per timestep)
export function randomEulerPole(rng: RNG): { axis: Point; rate: number } {
  const theta = 2 * Math.PI * rng.next();
  const phi = Math.acos(2 * rng.next() - 1);
  return {
    axis: {
      x: Math.sin(phi) * Math.cos(theta),
      y: Math.sin(phi) * Math.sin(theta),
      z: Math.cos(phi),
    },
    rate: 0.001 + rng.next() * 0.009,
  };
}

// Chord distance between two points on the unit sphere
export function chordDistance(a: Point, b: Point): number {
  const dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
  return Math.sqrt(dx * dx + dy * dy + dz * dz);
}

// 3D vector operations
export function cross3(a: Point, b: Point): Point {
  return {
    x: a.y * b.z - a.z * b.y,
    y: a.z * b.x - a.x * b.z,
    z: a.x * b.y - a.y * b.x,
  };
}

export function dot3(a: Point, b: Point): number {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

export function sub3(a: Point, b: Point): Point {
  return { x: a.x - b.x, y: a.y - b.y, z: a.z - b.z };
}

export function scale3(v: Point, s: number): Point {
  return { x: v.x * s, y: v.y * s, z: v.z * s };
}

export function magnitude(v: Point): number {
  return Math.sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

export function normalizeVec(v: Point): Point {
  const m = magnitude(v);
  if (m === 0) return { x: 0, y: 0, z: 0 };
  return { x: v.x / m, y: v.y / m, z: v.z / m };
}

// Generate evenly-distributed points on the sphere (Fibonacci sphere)
export function generateFibonacciSphere(samples: number, rng: RNG, jitter: number): Point[] {
  const points: Point[] = [];
  const phi = Math.PI * (3 - Math.sqrt(5));
  const spacing = Math.sqrt(4 * Math.PI / samples);

  for (let i = 0; i < samples; i++) {
    const y = 1 - (i / (samples - 1)) * 2;
    const radius = Math.sqrt(1 - y * y);
    const theta = phi * i;

    let x = Math.cos(theta) * radius;
    let z = Math.sin(theta) * radius;
    let py = y;

    if (jitter > 0) {
      x += (rng.next() - 0.5) * jitter * spacing * 1.5;
      py += (rng.next() - 0.5) * jitter * spacing * 1.5;
      z += (rng.next() - 0.5) * jitter * spacing * 1.5;
      const len = Math.sqrt(x*x + py*py + z*z);
      x /= len; py /= len; z /= len;
    }
    points.push({ x, y: py, z });
  }
  return points;
}