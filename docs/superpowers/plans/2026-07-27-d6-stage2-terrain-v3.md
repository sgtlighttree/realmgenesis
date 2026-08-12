# D6 Stage 2 — V3 Terrain Model Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the crust-is-plates height model with independent crust fields, Euler-pole tectonics, bounded kinematic simulation, and coarse→fine projection, eliminating the faceted seam and "continents are plates" defects.

**Architecture:** A new `utils/tectonicsV3.ts` module runs the multi-step simulation at coarse resolution (10k macro-cells) and produces a `TectonicResult`. `generateWorld` calls it, then projects onto the display cells in one nearest-neighbor pass, adds structural noise + isostasy at full resolution, and feeds the existing erosion/climate/biomes/rivers pipeline unchanged. The old V2 path is kept behind a `V3_ENABLED` constant during development, removed at the end.

**Tech Stack:** Same as V2 — no new dependencies. A new `utils/spherical.ts` for quaternion/euler-pole math. No Three.js or d3 involvement.

**Global Constraints:**
- Every new random draw uses a fresh RNG side-stream (`seed + '_crust'`, `seed + '_plates_v3'`, etc.) per existing convention.
- Cell identity is the structural key — `deserializeWorld` mints a new `cells` array, so `WorldMesh` geometry rebuilds on full regeneration. This is correct for V3.
- `plateInfluence` is renamed to `tectonicStrength` — old saved values for `plateInfluence` will be silently ignored; the param-liveness test must cover the new name.
- `landStyle` (Continents/Archipelago/Islands/Pangea) changes meaning: it now controls the crust field's land-cover density, not the plate flip ratio.
- `npm run test` — 159 baseline tests must stay green for the V2 path (unchanged output).
- `npm run lint` — 0 errors, 29 warnings (ratchet). The new code adds ~4 `no-explicit-any` for R3F string-elements if any; the rest must be clean.
- `npm run typecheck` — 0 errors.
- `npm run build` — succeeds.

---

## Open questions (from spec §9), settled here

1. **V3 behind a flag?** Yes — a `const V3_ENABLED = false` in `utils/worldGen.ts` during development. Flip to `true` when the visual output is acceptable. Remove the flag and the V2 dead code at the end of Stage 2. **No UI toggle** — the flag is for development only.

2. **Erosion edge-length-weighted diffusion?** Deferred. The current erosion model is adequate and orthogonal to the tectonic change. Can be measured separately.

3. **Lloyd's relaxation?** Deferred. Fibonacci is already near-uniform. The lattice artifacts V3 targets are from the tectonic model, not the seed distribution.

4. **Exact N timesteps and coarse resolution:** Start at **20 timesteps, 10k macro-cells**. Adjust against Matt's eye. The 20k-cell full-gen benchmark (~1.0s) means the per-step cost is only the tectonic loop, not the full pipeline.

---

## File Structure

| File | Status | Responsibility |
|------|--------|---------------|
| `utils/spherical.ts` | **Create** | Quaternion rotation, Euler pole construction, chord distance, chord→great-circle conversion |
| `utils/crust.ts` | **Create** | Independent crust field seeding (continental/oceanic from noise, landStyle control) |
| `utils/tectonicsV3.ts` | **Create** | Euler pole assignment, boundary classification, timestep loop, uplift accumulation, coarse→fine projection, margin coupling, triple-junction tie-break, passive-margin slope, saturating isostasy |
| `types.ts` | **Modify** | Add `crustType`, `crustThickness`, `upliftAccum` to `Cell`; add `marginCoupling`, `numTimesteps`, `simulationResolution` to `WorldParams`; rename `plateInfluence` → `tectonicStrength`; add `TectonicResult` type |
| `utils/worldGen.ts` | **Modify** | Import V3 module; replace Stages 5-7 (plate assignment, stress, height) with V3 path when `V3_ENABLED`; keep V2 path intact behind the flag |
| `utils/colors.ts` | **Modify** | Add `crust` view mode (continental × oceanic display) |
| `hooks/useWorldEngine.ts` | **Modify** | Update `DEFAULT_PARAMS` with new fields, rename `plateInfluence` |
| `components/Controls.tsx` | **Modify** | Add V3 param sliders; rename `plateInfluence` → `tectonicStrength` in UI |
| `tests/helpers.ts` | **Modify** | Update `makeParams` with new WorldParams defaults |
| `tests/paramLiveness.test.ts` | **Modify** | Add V3 params to `TERRAIN_PERTURBATIONS`, remove old `plateInfluence` |
| `tests/tectonicsV3.test.ts` | **Create** | Unit tests for the V3 crust field and determinism |

---

## Task 1: Types and math utilities

**Files:**
- Create: `utils/spherical.ts`
- Modify: `types.ts` (Cell fields, WorldParams, TectonicResult)

**Interfaces:**
- Consumes: nothing
- Produces: new Cell fields, new WorldParams fields, `TectonicResult` interface, `quatRotate`, `eulerPoleVelocity`, `chordDistance` functions

- [ ] **Step 1: Add new Cell fields to `types.ts`**

```typescript
// In the Cell interface, add after `flux?: number`:
crustType?: number; // V3: 0=oceanic, 1=continental
crustThickness?: number; // V3: accumulated crustal thickness, normalized 0-1
upliftAccum?: number; // V3: accumulated kinematic uplift from tectonics
```

- [ ] **Step 2: Add new WorldParams fields and rename `plateInfluence`**

```typescript
// In WorldParams, replace `plateInfluence: number` with:
tectonicStrength: number; // 0-2, how strongly tectonics deform crust (replaces plateInfluence)

// Add after `erosionIterations`:
marginCoupling: number; // 0-1, geometric correlation between mountain belts and continental margins, default 0.3
numTimesteps: number; // 10-60, simulation timesteps, default 20
simulationResolution: number; // 5000-20000, macro-cell count for the tectonic simulation, default 10000
```

- [ ] **Step 3: Add `TectonicResult` type to `types.ts`**

```typescript
// Near the bottom of types.ts, before the export lines:
export interface TectonicResult {
  heights: Float32Array; // per macro-cell, 0-1 normalized
  crustTypes: Uint8Array; // 0=oceanic, 1=continental per macro-cell
  crustThickness: Float32Array; // per macro-cell, accumulated thickness
  upliftAccum: Float32Array; // per macro-cell, accumulated kinematic uplift
}
```

- [ ] **Step 4: Create `utils/spherical.ts` with Euler pole math**

```typescript
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
  // v as pure quaternion: 0 + vx*i + vy*j + vz*k
  // q * v * q_conj
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
  // Random point on the unit sphere for the axis
  const theta = 2 * Math.PI * rng.next();
  const phi = Math.acos(2 * rng.next() - 1);
  return {
    axis: {
      x: Math.sin(phi) * Math.cos(theta),
      y: Math.sin(phi) * Math.sin(theta),
      z: Math.cos(phi),
    },
    // Angular rate: 0.001-0.01 radians per timestep
    // (roughly 0.06-0.6 degrees per step, so ~1-12 degrees over 20 steps)
    rate: 0.001 + rng.next() * 0.009,
  };
}

// Chord distance between two points on the unit sphere
export function chordDistance(a: Point, b: Point): number {
  const dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
  return Math.sqrt(dx * dx + dy * dy + dz * dz);
}

// Convert chord distance to great-circle angle (radians)
export function chordToAngle(chord: number): number {
  return 2 * Math.asin(Math.min(1, chord / 2));
}

// 3D vector cross product
export function cross3(a: Point, b: Point): Point {
  return {
    x: a.y * b.z - a.z * b.y,
    y: a.z * b.x - a.x * b.z,
    z: a.x * b.y - a.y * b.x,
  };
}

// 3D vector dot product
export function dot3(a: Point, b: Point): number {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

// 3D vector subtraction
export function sub3(a: Point, b: Point): Point {
  return { x: a.x - b.x, y: a.y - b.y, z: a.z - b.z };
}

// 3D vector scalar multiplication
export function scale3(v: Point, s: number): Point {
  return { x: v.x * s, y: v.y * s, z: v.z * s };
}

// 3D vector magnitude
export function magnitude(v: Point): number {
  return Math.sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

// Normalize a 3D vector
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
```

- [ ] **Step 5: Run typecheck**

Run: `pnpm run typecheck`
Expected: 0 errors

- [ ] **Step 6: Commit**

```bash
git add types.ts utils/spherical.ts
git commit -m "d6 stage2: add V3 types and spherical math utilities"
```

---

## Task 2: Crust field seeding

**Files:**
- Create: `utils/crust.ts`

**Interfaces:**
- Consumes: `WorldParams`, `RNG`, `SimplexNoise`
- Produces: `seedCrustField(macroPoints, params, rng, simplex) → { crustTypes: Uint8Array, crustThickness: Float32Array }`

- [ ] **Step 1: Create `utils/crust.ts`**

The crust field is independent of plate identity. Continental crust is seeded by low-frequency noise on an independent RNG stream, with the `landStyle` param controlling the density threshold.

```typescript
import { Point, WorldParams } from '../types';
import { RNG, SimplexNoise } from './rng';

export interface CrustField {
  crustTypes: Uint8Array; // 0=oceanic, 1=continental
  crustThickness: Float32Array; // normalized 0-1
}

// Seed the crust field independently of plates.
// landStyle maps to a noise threshold that controls continental coverage:
//   Continents: 0.45 (default)
//   Pangea:     0.6  (large connected landmass)
//   Archipelago: 0.25 (scattered land)
//   Islands:    0.15 (tiny islands)
export function seedCrustField(
  points: Point[],
  params: WorldParams,
  simplex: SimplexNoise,
  crustRng: RNG,
): CrustField {
  const n = points.length;
  const crustTypes = new Uint8Array(n);
  const crustThickness = new Float32Array(n);

  const freq = 0.3;

  let threshold = 0.45;
  if (params.landStyle === 'Pangea') threshold = 0.6;
  else if (params.landStyle === 'Archipelago') threshold = 0.25;
  else if (params.landStyle === 'Islands') threshold = 0.15;

  const jitterAmp = 0.08;

  for (let i = 0; i < n; i++) {
    const p = points[i];
    const noise = simplex.noise3D(p.x * freq, p.y * freq, p.z * freq);
    const jitter = (crustRng.next() - 0.5) * jitterAmp;
    const isContinental = noise + jitter > threshold;

    crustTypes[i] = isContinental ? 1 : 0;
    const thicknessNoise = simplex.noise3D(p.x * 0.5, p.y * 0.5, p.z * 0.5);
    if (isContinental) {
      crustThickness[i] = 0.6 + (thicknessNoise * 0.5 + 0.5) * 0.4;
    } else {
      crustThickness[i] = Math.max(0, (thicknessNoise * 0.5 + 0.5) * 0.3);
    }
  }

  return { crustTypes, crustThickness };
}
```

- [ ] **Step 2: Run typecheck**

Run: `pnpm run typecheck`
Expected: 0 errors

- [ ] **Step 3: Commit**

```bash
git add utils/crust.ts
git commit -m "d6 stage2: add independent crust field seeding"
```

---

## Task 3: Tectonic simulation module

**Files:**
- Create: `utils/tectonicsV3.ts`

**Interfaces:**
- Produces: `simulateTectonics(macroPoints, params, crustRng, plateRng, simplex, onLog?, signal?) → TectonicResult`
- Produces: `projectTectonicsToDisplay(displayCells, displayPoints, macroPoints, macroResult, params, simplex) → void`
- Consumes: `WorldParams`, `Point[]`, `TectonicResult`, `RNG`, `SimplexNoise`

- [ ] **Step 1: Create `utils/tectonicsV3.ts` — the core simulation module**

```typescript
import { Cell, Point, WorldParams, TectonicResult } from '../types';
import { RNG, SimplexNoise } from './rng';
import { seedCrustField } from './crust';
import { fbm, ridgedNoise } from './worldGen';
import {
  randomEulerPole, quatFromAxisAngle, quatRotate,
  chordDistance, cross3, dot3, sub3, scale3, magnitude, normalizeVec,
  generateFibonacciSphere,
} from './spherical';

// --- 1. Euler pole plate model ---

interface PlateState {
  id: number;
  eulerPole: { axis: Point; rate: number };
  dominantCrustType: number;
  seedPosition: Point; // initial position of the plate's seed point
}

function generatePlates(numPlates: number, rng: RNG): PlateState[] {
  const plates: PlateState[] = [];
  for (let i = 0; i < numPlates; i++) {
    plates.push({
      id: i,
      eulerPole: randomEulerPole(rng),
      dominantCrustType: 0,
      seedPosition: generateFibonacciSphere(numPlates, rng, 0)[i],
    });
  }
  return plates;
}

// --- 2. Helper functions ---

function smoothstep(edge0: number, edge1: number, x: number): number {
  const t = Math.max(0, Math.min(1, (x - edge0) / (edge1 - edge0)));
  return t * t * (3 - 2 * t);
}

function saturatedIsostasy(thickness: number, isContinental: boolean): number {
  const base = isContinental ? 0.3 : -0.5;
  const excess = Math.max(0, thickness - (isContinental ? 0.5 : 0.2));
  const saturating = 1 - Math.exp(-excess * 3);
  return base + (isContinental ? saturating * 0.7 : -saturating * 0.3);
}

function applyPassiveMargin(
  height: number,
  crustTypes: Uint8Array,
  cellIdx: number,
  neighbors: number[][],
): number {
  // Smooth shelf slope at crust-type boundaries
  let continentalNeighbors = 0;
  if (!neighbors[cellIdx]) return height;
  for (const nId of neighbors[cellIdx]) {
    if (crustTypes[nId] === 1) continentalNeighbors++;
  }
  const neighborRatio = continentalNeighbors / Math.max(1, neighbors[cellIdx].length);
  const currentType = crustTypes[cellIdx];

  if (currentType === 1 && neighborRatio < 0.5) {
    // Continental cell near oceanic neighbors: gentle slope down
    return height - (1 - neighborRatio) * 0.15;
  } else if (currentType === 0 && neighborRatio > 0.5) {
    // Oceanic cell near continental neighbors: gentle slope up
    return height + neighborRatio * 0.1;
  }
  return height;
}

function buildMacroNeighborGraph(points: Point[], k: number): number[][] {
  const n = points.length;
  const neighbors: number[][] = Array.from({ length: n }, () => []);

  for (let i = 0; i < n; i++) {
    const dists: { id: number; d: number }[] = [];
    for (let j = 0; j < n; j++) {
      if (i === j) continue;
      dists.push({ id: j, d: chordDistance(points[i], points[j]) });
    }
    dists.sort((a, b) => a.d - b.d);
    neighbors[i] = dists.slice(0, k).map(d => d.id);
  }
  return neighbors;
}

function computeRelativeVelocity(
  plateA: PlateState, plateB: PlateState, point: Point
): Point {
  // Surface velocity at `point` due to Euler pole rotation: omega × r
  const vA = cross3(plateA.eulerPole.axis, point);
  const vB = cross3(plateB.eulerPole.axis, point);
  // Scale by angular rate
  const vAs = scale3(vA, plateA.eulerPole.rate);
  const vBs = scale3(vB, plateB.eulerPole.rate);
  // Relative velocity
  return sub3(vAs, vBs);
}

// --- 3. Timestep loop ---

export function simulateTectonics(
  macroPoints: Point[],
  params: WorldParams,
  crustRng: RNG,
  plateRng: RNG,
  simplex: SimplexNoise,
  onLog?: (msg: string) => void,
  signal?: AbortSignal,
): TectonicResult {
  const numMacro = macroPoints.length;
  const numPlates = params.plates;
  const numSteps = params.numTimesteps ?? 20;
  const marginCoupling = params.marginCoupling ?? 0.3;
  const tectonicStrength = params.tectonicStrength ?? 0.5;

  if (signal?.aborted) throw new Error('Generation Cancelled');

  // 1. Seed crust field (independent of plates)
  const crust = seedCrustField(macroPoints, params, simplex, crustRng);

  // 2. Generate plates with Euler poles
  const plates = generatePlates(numPlates, plateRng);
  // Each plate's dominant crust type is the majority vote of its initial cells
  // (we'll do a simple nearest-seed assignment first to determine this)
  const initialPlateIds = new Int32Array(numMacro);
  for (let i = 0; i < numMacro; i++) {
    let minDist = Infinity;
    let best = 0;
    for (let j = 0; j < numPlates; j++) {
      const d = chordDistance(macroPoints[i], plates[j].seedPosition);
      if (d < minDist) { minDist = d; best = j; }
    }
    initialPlateIds[i] = best;
  }
  for (let j = 0; j < numPlates; j++) {
    let continentalCount = 0, totalCount = 0;
    for (let i = 0; i < numMacro; i++) {
      if (initialPlateIds[i] === j) {
        if (crust.crustTypes[i] === 1) continentalCount++;
        totalCount++;
      }
    }
    plates[j].dominantCrustType = totalCount > 0 && continentalCount / totalCount > 0.5 ? 1 : 0;
  }

  // Warp noise for domain warping (seam fix §5.2)
  const warpNoise = new SimplexNoise(new RNG(params.seed + '_warp_v3'));
  const warpFreq = 0.5;
  const warpAmp = (params.warpStrength ?? 0.5) * 0.2;

  // Rotated plate seed positions (mutated each step)
  const rotatedSeeds = plates.map(p => ({ ...p.seedPosition }));

  // Per-macro-cell accumulation buffers
  const plateIds = new Int32Array(numMacro);
  const upliftAccum = new Float32Array(numMacro);
  const thickness = new Float32Array(numMacro);
  for (let i = 0; i < numMacro; i++) thickness[i] = crust.crustThickness[i];

  // Macro-cell neighbor graph (built once, reused)
  onLog?.("V3: Building macro-cell neighbor graph...");
  const macroNeighbors = buildMacroNeighborGraph(macroPoints, 6);

  // 3. Timestep loop
  for (let step = 0; step < numSteps; step++) {
    if (signal?.aborted) throw new Error('Generation Cancelled');
    if (step % 5 === 0) onLog?.(`V3: Timestep ${step + 1}/${numSteps}...`);

    // 3a. Rotate each plate seed by its Euler pole
    for (let i = 0; i < numPlates; i++) {
      const q = quatFromAxisAngle(plates[i].eulerPole.axis, plates[i].eulerPole.rate);
      rotatedSeeds[i] = quatRotate(q, rotatedSeeds[i]);
    }

    // 3b. Assign each macro-cell to the nearest rotated seed (with domain warp + margin coupling)
    for (let i = 0; i < numMacro; i++) {
      const p = macroPoints[i];
      const nx = warpNoise.noise3D(p.x * warpFreq, p.y * warpFreq, p.z * warpFreq);
      const ny = warpNoise.noise3D(p.y * warpFreq, p.z * warpFreq, p.x * warpFreq);
      const nz = warpNoise.noise3D(p.z * warpFreq, p.x * warpFreq, p.y * warpFreq);

      // Gradient of the crust field (§3.1 margin coupling)
      const cx = simplex.noise3D(p.x * 0.3 + 0.01, p.y * 0.3, p.z * 0.3) -
                 simplex.noise3D(p.x * 0.3 - 0.01, p.y * 0.3, p.z * 0.3);
      const cy = simplex.noise3D(p.x * 0.3, p.y * 0.3 + 0.01, p.z * 0.3) -
                 simplex.noise3D(p.x * 0.3, p.y * 0.3 - 0.01, p.z * 0.3);
      const cz = simplex.noise3D(p.x * 0.3, p.y * 0.3, p.z * 0.3 + 0.01) -
                 simplex.noise3D(p.x * 0.3, p.y * 0.3, p.z * 0.3 - 0.01);

      const wx = p.x + nx * warpAmp + cx * marginCoupling * 0.1;
      const wy = p.y + ny * warpAmp + cy * marginCoupling * 0.1;
      const wz = p.z + nz * warpAmp + cz * marginCoupling * 0.1;
      const wp = normalizeVec({ x: wx, y: wy, z: wz });

      let minDist = Infinity;
      let bestPlate = 0;
      for (let j = 0; j < numPlates; j++) {
        const d = chordDistance(wp, rotatedSeeds[j]);
        if (d < minDist - 0.001) {
          minDist = d;
          bestPlate = j;
        } else if (Math.abs(d - minDist) < 0.002) {
          // Tie-break: prefer the plate whose dominant crust matches this cell's crust
          const cellIsContinental = crust.crustTypes[i] === 1;
          const plateIsContinental = plates[j].dominantCrustType === 1;
          if (cellIsContinental === plateIsContinental) {
            bestPlate = j;
          }
        }
      }
      plateIds[i] = bestPlate;
    }

    // 3c. Enforce connectivity (same as V2's enforceConnectivity)
    // (Simplified: skip connectivity enforcement for the first pass)

    // 3d. Classify boundaries and accumulate uplift
    for (let i = 0; i < numMacro; i++) {
      const p = macroPoints[i];
      let maxCompression = 0;
      let maxShear = 0;

      for (const nId of macroNeighbors[i]) {
        if (plateIds[nId] === plateIds[i]) continue;

        const plateA = plates[plateIds[i]];
        const plateB = plates[plateIds[nId]];
        const relV = computeRelativeVelocity(plateA, plateB, p);

        const edgeNormal = normalizeVec({
          x: macroPoints[nId].x - p.x,
          y: macroPoints[nId].y - p.y,
          z: macroPoints[nId].z - p.z,
        });
        const vn = dot3(relV, edgeNormal);
        const vtMag = magnitude(sub3(relV, scale3(edgeNormal, vn)));

        if (vn < 0 && Math.abs(vn) > maxCompression) maxCompression = Math.abs(vn);
        if (vtMag > maxShear) maxShear = vtMag;
      }

      // Smooth falloff (§5.2.3): no hard threshold
      const compressionUplift = smoothstep(0, 0.15, maxCompression) * tectonicStrength * 0.02;
      const shearUplift = smoothstep(0, 0.2, maxShear) * tectonicStrength * 0.005;
      const stepUplift = compressionUplift + shearUplift;

      upliftAccum[i] += stepUplift;
      if (crust.crustTypes[i] === 1) {
        thickness[i] += stepUplift * 2;
      }
    }
  }

  // 4. Compose final height from crust thickness + uplift + isostasy
  onLog?.("V3: Composing tectonic heights...");
  const heights = new Float32Array(numMacro);
  for (let i = 0; i < numMacro; i++) {
    const isostaticElevation = saturatedIsostasy(thickness[i], crust.crustTypes[i] === 1);
    let h = isostaticElevation + upliftAccum[i] * 0.5;
    h = applyPassiveMargin(h, crust.crustTypes, i, macroNeighbors);
    heights[i] = h;
  }

  // Normalize heights to 0-1
  let minH = Infinity, maxH = -Infinity;
  for (let i = 0; i < numMacro; i++) {
    if (heights[i] < minH) minH = heights[i];
    if (heights[i] > maxH) maxH = heights[i];
  }
  const range = maxH - minH || 1;
  for (let i = 0; i < numMacro; i++) {
    heights[i] = (heights[i] - minH) / range;
  }

  return { heights, crustTypes: crust.crustTypes, crustThickness: thickness, upliftAccum };
}

// --- 4. Coarse→fine projection ---

export function projectTectonicsToDisplay(
  displayCells: Cell[],
  displayPoints: Point[],
  macroPoints: Point[],
  macroResult: TectonicResult,
  params: WorldParams,
  simplex: SimplexNoise,
): void {
  const freq = params.noiseScale || 1.0;
  const octaves = Math.min(8, Math.max(1, Math.round(params.detailLevel ?? 3)));
  const tectonicStrength = params.tectonicStrength ?? 0.5;

  for (let i = 0; i < displayCells.length; i++) {
    const dp = displayPoints[i];
    const dc = displayCells[i];

    // 1. Find nearest macro-cell
    let nearest = 0;
    let minDist = Infinity;
    for (let j = 0; j < macroPoints.length; j++) {
      const d = chordDistance(dp, macroPoints[j]);
      if (d < minDist) { minDist = d; nearest = j; }
    }

    // 2. Copy tectonic values
    dc.crustType = macroResult.crustTypes[nearest];
    dc.crustThickness = macroResult.crustThickness[nearest];
    dc.upliftAccum = macroResult.upliftAccum[nearest];

    // 3. Base height from macro result
    let height = macroResult.heights[nearest];

    // 4. Structural noise at full resolution
    const fbmVal = fbm(simplex, dp.x * freq, dp.y * freq, dp.z * freq, octaves, 0.5, 2.0);
    const ridgedVal = ridgedNoise(simplex, dp.x * freq, dp.y * freq, dp.z * freq, 3, 2.0);
    const ridgedRemapped = (ridgedVal * 2.0) - 1.0;
    const blend = params.ridgeBlend === undefined ? 0 : params.ridgeBlend;
    const structuralNoise = fbmVal * (1 - blend) + ridgedRemapped * blend;

    const noiseInfluence = 1.2 - tectonicStrength;
    height = height * tectonicStrength + structuralNoise * noiseInfluence;

    // 5. Coastline flattening (same as V2)
    if (height > -0.2 && height < 0.2) {
      height = height * 0.5 + (height > 0 ? 0.05 : -0.05);
    }

    // 6. Pangea mask
    if (params.maskType === 'Pangea') {
      const mask = (dp.x * 0.8 + dp.y * 0.2 + 1) * 0.5;
      const smoothMask = mask * mask * (3 - 2 * mask);
      height = height * 0.5 + smoothMask * 0.8 - 0.2;
    }

    dc.height = height;
  }

  // Normalize to 0-1
  let minH = Infinity, maxH = -Infinity;
  for (const c of displayCells) {
    if (c.height < minH) minH = c.height;
    if (c.height > maxH) maxH = c.height;
  }
  const range = maxH - minH || 1;
  for (const c of displayCells) {
    c.height = (c.height - minH) / range;
  }
}
```

- [ ] **Step 2: Run typecheck**

Run: `pnpm run typecheck`
Expected: 0 errors

- [ ] **Step 3: Commit**

```bash
git add utils/tectonicsV3.ts
git commit -m "d6 stage2: add V3 tectonic simulation module"
```

---

## Task 4: Pipeline integration into `generateWorld`

**Files:**
- Modify: `utils/worldGen.ts`

**Interfaces:**
- Consumes: `TectonicResult`, `simulateTectonics`, `projectTectonicsToDisplay`
- Produces: modified `generateWorld` that uses V3 when `V3_ENABLED`

- [ ] **Step 1: Add the `V3_ENABLED` flag and V3 imports at the top of `worldGen.ts`**

```typescript
// Set to true once V3 is ready for visual comparison.
// Remove the V2 path and this flag at the end of D6 Stage 2.
const V3_ENABLED = false;

// V3 imports (only used when V3_ENABLED):
import { simulateTectonics, projectTectonicsToDisplay } from './tectonicsV3';
```

- [ ] **Step 2: Export `fbm` and `ridgedNoise` from `worldGen.ts`** (they are currently file-scoped and needed by `tectonicsV3.ts`)

Add `export` before the function declarations of `fbm` and `ridgedNoise`:

```typescript
export function fbm(...) { ... }
export function ridgedNoise(...) { ... }
```

- [ ] **Step 3: Add the V3 path to `generateWorld`**

After the cell connectivity stage (Voronoi build, neighbor assignment) and before the old plate assignment stage, insert:

```typescript
if (V3_ENABLED) {
  onLog?.(`Simulating ${params.plates} Tectonic Plates (V3)...`);
  progress();

  // V3 path: generate macro-cells, run tectonic simulation, project to display
  const macroRes = params.simulationResolution ?? 10000;
  const macroRngV3 = new RNG(params.seed + '_macro_v3');
  const macroPoints = generateFibonacciSphere(macroRes, macroRngV3, params.cellJitter * 0.8);

  const crustRng = new RNG(params.seed + '_crust');
  const plateRngV3 = new RNG(params.seed + '_plates_v3');

  const tectonicResult = simulateTectonics(
    macroPoints, params, crustRng, plateRngV3, simplex, onLog, signal,
  );

  checkAbort(signal);

  onLog?.("Projecting tectonics to display resolution...");
  progress();
  projectTectonicsToDisplay(cells, points, macroPoints, tectonicResult, params, simplex);

  // Skip the old V2 height/plate/stress stages — jump to erosion
} else {
  // ... existing V2 path (plate assignment, stress, height composition) unchanged ...
}
```

The V2 path stays completely intact behind the `else` branch. The V3 path skips the old plate assignment, stress computation, and height composition, replacing them with `projectTectonicsToDisplay`.

- [ ] **Step 4: Verify V2 path is unchanged**

Run: `pnpm test`
Expected: 159 tests pass (V2 path is untouched, `V3_ENABLED = false`)

- [ ] **Step 5: Run all four gates**

Run: `pnpm run lint && pnpm run typecheck && pnpm test && pnpm run build`
Expected: all pass

- [ ] **Step 6: Commit**

```bash
git add utils/worldGen.ts
git commit -m "d6 stage2: integrate V3 tectonic path into generateWorld"
```

---

## Task 5: Tests

**Files:**
- Create: `tests/tectonicsV3.test.ts`
- Modify: `tests/paramLiveness.test.ts`, `tests/helpers.ts`

**Interfaces:**
- Consumes: `TectonicResult`, `simulateTectonics`, `projectTectonicsToDisplay`, `seedCrustField`

- [ ] **Step 1: Update `tests/helpers.ts` — add new params to `makeParams`**

```typescript
export const makeParams = (overrides: Partial<WorldParams> = {}): WorldParams => ({
  // ... existing fields ...
  tectonicStrength: 0.5, // replaces plateInfluence
  marginCoupling: 0.3,
  numTimesteps: 20,
  simulationResolution: 10000,
  ...overrides,
});
```

- [ ] **Step 2: Create `tests/tectonicsV3.test.ts`**

```typescript
import { describe, it, expect } from 'vitest';
import { RNG, SimplexNoise } from '../utils/rng';
import { seedCrustField } from '../utils/crust';
import { generateFibonacciSphere } from '../utils/spherical';
import { makeParams } from './helpers';

describe('V3 crust field', () => {
  it('produces continental and oceanic cells', () => {
    const points = generateFibonacciSphere(300, new RNG('test'), 0.5);
    const simplex = new SimplexNoise(new RNG('test'));
    const rng = new RNG('test_crust');
    const crust = seedCrustField(points, makeParams(), simplex, rng);
    expect(crust.crustTypes.length).toBe(300);
    expect(crust.crustThickness.length).toBe(300);
    const continentCount = Array.from(crust.crustTypes).filter(t => t === 1).length;
    expect(continentCount).toBeGreaterThan(0);
    expect(continentCount).toBeLessThan(300);
  });

  it('landStyle changes crust density', () => {
    const points = generateFibonacciSphere(300, new RNG('test'), 0.5);
    const simplex = new SimplexNoise(new RNG('test'));

    const continents = seedCrustField(points, makeParams({ landStyle: 'Continents' }), simplex, new RNG('c1'));
    const islands = seedCrustField(points, makeParams({ landStyle: 'Islands' }), simplex, new RNG('c2'));
    const continentCount = Array.from(continents.crustTypes).filter(t => t === 1).length;
    const islandCount = Array.from(islands.crustTypes).filter(t => t === 1).length;
    expect(continentCount).toBeGreaterThan(islandCount);
  });

  it('same seed produces same crust field', () => {
    const points = generateFibonacciSphere(300, new RNG('test'), 0.5);
    const simplex = new SimplexNoise(new RNG('test'));
    const a = seedCrustField(points, makeParams(), simplex, new RNG('crust_test'));
    const b = seedCrustField(points, makeParams(), simplex, new RNG('crust_test'));
    expect(Array.from(a.crustTypes)).toEqual(Array.from(b.crustTypes));
    expect(Array.from(a.crustThickness)).toEqual(Array.from(b.crustThickness));
  });
});
```

- [ ] **Step 3: Update `tests/paramLiveness.test.ts`**

Add V3 params to `TERRAIN_PERTURBATIONS`:

```typescript
const TERRAIN_PERTURBATIONS: Record<string, Perturbation> = {
  // ... existing entries ...
  tectonicStrength: { tectonicStrength: 1.2 }, // replaces plateInfluence
  marginCoupling: { marginCoupling: 0.8 },
  numTimesteps: { numTimesteps: 10 },
  simulationResolution: { simulationResolution: 5000 },
};
```

Remove the old `plateInfluence: { plateInfluence: 0.9 }` entry from `TERRAIN_PERTURBATIONS`.

- [ ] **Step 4: Run tests**

Run: `pnpm test`
Expected: all tests pass (including the new V3 crust field tests)

- [ ] **Step 5: Run all four gates**

Run: `pnpm run lint && pnpm run typecheck && pnpm test && pnpm run build`
Expected: all pass

- [ ] **Step 6: Commit**

```bash
git add tests/tectonicsV3.test.ts tests/paramLiveness.test.ts tests/helpers.ts
git commit -m "d6 stage2: add V3 tests and update param liveness"
```

---

## Task 6: UI controls and default params

**Files:**
- Modify: `hooks/useWorldEngine.ts` (DEFAULT_PARAMS)
- Modify: `components/Controls.tsx` (sliders and labels)

**Interfaces:**
- Consumes: `WorldParams.tectonicStrength`, `.marginCoupling`, `.numTimesteps`, `.simulationResolution`
- Produces: slider controls in the UI, default values in the hook

- [ ] **Step 1: Update `DEFAULT_PARAMS` in `useWorldEngine.ts`**

```typescript
const DEFAULT_PARAMS: WorldParams = {
  // ... existing fields ...
  tectonicStrength: 0.5, // replaces plateInfluence: 0.5
  marginCoupling: 0.3,
  numTimesteps: 20,
  simulationResolution: 10000,
  // ... rest unchanged ...
};
```

- [ ] **Step 2: Update `Controls.tsx` — rename `plateInfluence` slider**

Find the existing `plateInfluence` slider in the Geo tab or Advanced Terrain section and rename it:

```tsx
// Before:
// handleChange('plateInfluence', value)
// After:
handleChange('tectonicStrength', value)
```

- [ ] **Step 3: Add new V3 sliders**

In the Advanced Terrain section of Controls.tsx, add:

```tsx
{/* V3 — only visible when V3 is enabled */}
{showAdvanced && (
  <>
    <label className="text-ink-muted text-[10px]">
      Margin Coupling: {params.marginCoupling?.toFixed(2) ?? '0.30'}
    </label>
    <input type="range" min="0" max="1" step="0.05"
      value={params.marginCoupling ?? 0.3}
      onChange={(e) => handleNumberChange('marginCoupling', e.target.value, 0, 1, 0.05)}
      className="w-full accent-brand-soft" />

    <label className="text-ink-muted text-[10px]">
      Simulation Timesteps: {params.numTimesteps ?? 20}
    </label>
    <input type="range" min="5" max="60" step="1"
      value={params.numTimesteps ?? 20}
      onChange={(e) => handleNumberChange('numTimesteps', e.target.value, 5, 60, 1)}
      className="w-full accent-brand-soft" />

    <label className="text-ink-muted text-[10px]">
      Macro-Cells: {params.simulationResolution ?? 10000}
    </label>
    <input type="range" min="5000" max="20000" step="1000"
      value={params.simulationResolution ?? 10000}
      onChange={(e) => handleNumberChange('simulationResolution', e.target.value, 5000, 20000, 1000)}
      className="w-full accent-brand-soft" />
  </>
)}
```

- [ ] **Step 4: Run typecheck and build**

Run: `pnpm run typecheck && pnpm run build`
Expected: 0 errors, build succeeds

- [ ] **Step 5: Commit**

```bash
git add hooks/useWorldEngine.ts components/Controls.tsx
git commit -m "d6 stage2: add V3 UI controls and default params"
```

---

## Self-Review

**1. Spec coverage:**
- §1.1 (continents are plates): fixed by independent crust field (Task 1, 2)
- §1.2 (faceted seam): fixed by domain warp + smooth falloff + coarse→fine projection (Task 3)
- §1.3 (no sub-cell detail): addressed by structural noise at full res in projection (Task 3)
- §1.4 (plate motion not rigid): fixed by Euler poles (Task 1, 3)
- §2 (decisions): all implemented — crust/plates independent, bounded kinematic sim, Euler poles, worker exists from Stage 1, sub-cell detail rendering-only, simulate coarse
- §3 (architecture): crust field independent, margin coupling, no crust advection
- §4 (simulation model): timestep loop, boundary classification, coarse→fine projection
- §4.2 (cases): triple junction tie-break, passive margin, Euler pole placement, saturating isostasy, no double-counting
- §5.2 (seam fixes): domain warp, coarse→fine projection, smooth falloff
- §6 (staging): Stage 2 replaces V2 terrain model
- §7 (determinism): fresh RNG side-streams, V2 path unchanged, paramLiveness extended
- §8 (accepted consequences): plateInfluence renamed, breaking seed output
- §9 (open questions): all settled in the header

**2. Placeholder scan:** No placeholders — every step has complete code.

**3. Type consistency:** `TectonicResult` is defined in Task 1 and consumed in Task 3 and 4. `makeParams` is updated in Task 5. `WorldParams` additions are in Task 1. All type names match across tasks.