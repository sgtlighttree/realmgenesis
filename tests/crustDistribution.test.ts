import { describe, it, expect } from 'vitest';
import { seedCrustField } from '../utils/crust';
import { RNG, SimplexNoise } from '../utils/rng';
import { DEFAULT_PARAMS } from '../utils/defaultParams';
import { Point, LandStyle } from '../types';

// Regression guard for D9 (pangea bias, 2026-08-22). The crust field once
// sampled a single low-frequency noise octave, collapsing all land into one
// hemisphere-sized component on EVERY seed. These tests assert the Continents
// field spreads land across several separate masses, and that the Pangea style
// still intentionally clumps. See utils/crust.ts for the tuning rationale.

function fibSphere(n: number): Point[] {
  const pts: Point[] = [];
  const phi = Math.PI * (3 - Math.sqrt(5));
  for (let i = 0; i < n; i++) {
    const y = 1 - (i / (n - 1)) * 2;
    const r = Math.sqrt(1 - y * y);
    const t = phi * i;
    pts.push({ x: Math.cos(t) * r, y, z: Math.sin(t) * r });
  }
  return pts;
}

// Undirected kNN=6 graph (mirrors buildMacroNeighborGraph's connectivity).
function knn6(pts: Point[]): number[][] {
  const n = pts.length;
  const k = 6;
  const bd = Array.from({ length: n }, () => new Float64Array(k).fill(Infinity));
  const bi = Array.from({ length: n }, () => new Int32Array(k).fill(-1));
  const upd = (d: Float64Array, id: Int32Array, cand: number, cd: number) => {
    let wi = 0, wd = d[0];
    for (let t = 1; t < k; t++) if (d[t] > wd) { wd = d[t]; wi = t; }
    if (wd > cd) { d[wi] = cd; id[wi] = cand; }
  };
  for (let i = 0; i < n; i++) {
    const pi = pts[i];
    for (let j = i + 1; j < n; j++) {
      const dx = pi.x - pts[j].x, dy = pi.y - pts[j].y, dz = pi.z - pts[j].z;
      const d = dx * dx + dy * dy + dz * dz;
      upd(bd[i], bi[i], j, d);
      upd(bd[j], bi[j], i, d);
    }
  }
  const sets: Set<number>[] = Array.from({ length: n }, () => new Set<number>());
  for (let i = 0; i < n; i++) for (const j of bi[i]) if (j >= 0) { sets[i].add(j); sets[j].add(i); }
  return sets.map((s) => Array.from(s));
}

// Connected-component sizes of continental cells, descending.
function componentSizes(nb: number[][], cont: Uint8Array): number[] {
  const n = cont.length;
  const comp = new Int32Array(n).fill(-1);
  const sizes: number[] = [];
  for (let i = 0; i < n; i++) {
    if (cont[i] !== 1 || comp[i] >= 0) continue;
    let sz = 0;
    const st = [i];
    comp[i] = sizes.length;
    while (st.length) {
      const c = st.pop()!;
      sz++;
      for (const m of nb[c]) if (cont[m] === 1 && comp[m] < 0) { comp[m] = sizes.length; st.push(m); }
    }
    sizes.push(sz);
  }
  return sizes.sort((a, b) => b - a);
}

const N = 6000;
const pts = fibSphere(N);
const nb = knn6(pts);
const SEEDS = ['realmgenesis', 'abc', 'xyz123', 'seed42', 'hello'];

function measure(style: LandStyle) {
  let landFrac = 0, majorMasses = 0, largest = 0;
  for (const seed of SEEDS) {
    const simplex = new SimplexNoise(new RNG(seed + '_crust'));
    const field = seedCrustField(pts, { ...DEFAULT_PARAMS, landStyle: style }, simplex, new RNG(seed + '_crust_jit'));
    const cont = field.crustTypes;
    let land = 0;
    for (let i = 0; i < N; i++) if (cont[i] === 1) land++;
    const sizes = componentSizes(nb, cont);
    const total = land || 1;
    landFrac += land / N;
    majorMasses += sizes.filter((s) => s >= 0.03 * total).length;
    largest += (sizes[0] || 0) / total;
  }
  const d = SEEDS.length;
  return { landFrac: landFrac / d, majorMasses: majorMasses / d, largest: largest / d };
}

describe('crust distribution (D9 pangea guard)', () => {
  it('Continents spreads land across several separate masses, not one blob', () => {
    const m = measure('Continents');
    // A pangea = one component holding ~all land (largest ~1.0, majorMasses ~1).
    expect(m.largest).toBeLessThan(0.75);
    expect(m.majorMasses).toBeGreaterThanOrEqual(3);
    // Sanity: keep a plausible ocean-dominated land fraction.
    expect(m.landFrac).toBeGreaterThan(0.15);
    expect(m.landFrac).toBeLessThan(0.55);
  });

  it('Pangea style still produces one dominant supercontinent', () => {
    const m = measure('Pangea');
    expect(m.largest).toBeGreaterThan(0.7);
  });

  it('Islands style fragments land into many small masses', () => {
    const m = measure('Islands');
    expect(m.largest).toBeLessThan(0.3);
    expect(m.majorMasses).toBeGreaterThan(8);
  });
});
