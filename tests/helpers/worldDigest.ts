import { createHash } from 'node:crypto';
import { WorldData } from '../../types';

export type WorldDigest = Record<string, string>;

// Canonical Cell field list. Deliberately explicit rather than Object.keys:
// a field silently disappearing across the worker boundary must show up as a
// digest change, and Object.keys would simply stop emitting it on both sides.
export const CELL_FIELDS = [
  'id', 'height', 'plateId', 'temperature', 'moisture', 'biome', 'flux',
  'regionId', 'provinceId', 'isCapital', 'isTown', 'population',
  'cultureId', 'religionId',
] as const;

// Exact IEEE-754 bits, not toFixed. A Float64 -> Float32 downcast of 0.5123456
// survives toFixed(6) but changes these bits, which is precisely the class of
// bug this stage exists to catch.
const f64 = new DataView(new ArrayBuffer(8));
const bits = (n: number): string => {
  f64.setFloat64(0, n);
  return f64.getBigUint64(0).toString(16).padStart(16, '0');
};

// undefined, null, -0 and NaN each get a distinct, stable encoding.
const enc = (v: unknown): string => {
  if (v === undefined) return 'u';
  if (v === null) return 'n';
  if (typeof v === 'number') return Object.is(v, -0) ? '-0' : bits(v);
  if (typeof v === 'boolean') return v ? 'T' : 'F';
  if (typeof v === 'string') return 's' + v;
  if (Array.isArray(v)) return '[' + v.map(enc).join(',') + ']';
  if (typeof v === 'object') {
    const o = v as Record<string, unknown>;
    return '{' + Object.keys(o).sort().map(k => k + ':' + enc(o[k])).join(',') + '}';
  }
  return 'x' + String(v);
};

const hash = (s: string): string => createHash('sha256').update(s).digest('hex').slice(0, 32);

export const digestWorld = (world: WorldData): WorldDigest => {
  const d: WorldDigest = {};
  d['cellCount'] = hash(String(world.cells.length));

  for (const field of CELL_FIELDS) {
    d[`cell.${field}`] = hash(
      world.cells.map(c => enc((c as unknown as Record<string, unknown>)[field])).join(';'),
    );
  }
  d['cell.center'] = hash(world.cells.map(c => enc(c.center)).join(';'));
  d['cell.vertices'] = hash(world.cells.map(c => enc(c.vertices)).join(';'));
  d['cell.neighbors'] = hash(world.cells.map(c => enc(c.neighbors)).join(';'));

  d['geoJson'] = hash(enc(world.geoJson as unknown));
  d['params'] = hash(enc(world.params as unknown));
  d['rivers'] = hash(enc(world.rivers as unknown));
  d['lakes'] = hash(enc(world.lakes as unknown));
  d['features'] = hash(enc(world.features as unknown));
  d['cultures'] = hash(enc(world.cultures as unknown));
  d['religions'] = hash(enc(world.religions as unknown));
  d['civData'] = hash(enc(world.civData as unknown));
  d['markers'] = hash(enc(world.markers as unknown));
  d['routes'] = hash(enc(world.routes as unknown));
  return d;
};

export const diffDigests = (a: WorldDigest, b: WorldDigest): string[] => {
  const keys = new Set([...Object.keys(a), ...Object.keys(b)]);
  return [...keys].filter(k => a[k] !== b[k]).sort();
};
