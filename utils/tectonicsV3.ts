import { Cell, Point, WorldParams, TectonicResult } from '../types';
import { RNG, SimplexNoise } from './rng';
import { seedCrustField } from './crust';
import {
  randomEulerPole, quatFromAxisAngle, quatRotate,
  chordDistance, cross3, dot3, sub3, scale3, magnitude, normalizeVec,
  generateFibonacciSphere,
} from './spherical';

// Noise helpers (duplicated from worldGen.ts to avoid circular imports)
// The V3 module will import from worldGen once V3_ENABLED becomes true
// and the old V2 path is removed — for now these are standalone copies.
function fbm(simplex: SimplexNoise, x: number, y: number, z: number, octaves: number, persistence: number, lacunarity: number): number {
    let total = 0;
    let frequency = 1;
    let amplitude = 1;
    let maxValue = 0;
    for(let i=0;i<octaves;i++) {
        total += simplex.noise3D(x * frequency, y * frequency, z * frequency) * amplitude;
        maxValue += amplitude;
        amplitude *= persistence;
        frequency *= lacunarity;
    }
    return total / maxValue;
}

function ridgedNoise(simplex: SimplexNoise, x: number, y: number, z: number, octaves: number, lacunarity: number): number {
    let total = 0;
    let frequency = 1;
    let amplitude = 1;
    let weight = 1;
    let max = 0;
    for (let i = 0; i < octaves; i++) {
        let signal = simplex.noise3D(x * frequency, y * frequency, z * frequency);
        signal = 1.0 - Math.abs(signal);
        signal *= signal;
        signal *= weight;
        weight = signal * 2;
        if (weight > 1) weight = 1;
        if (weight < 0) weight = 0;
        total += signal * amplitude;
        max += amplitude;
        amplitude *= 0.5;
        frequency *= lacunarity;
    }
    return total / max;
}

// --- 1. Euler pole plate model ---

interface PlateState {
  id: number;
  eulerPole: { axis: Point; rate: number };
  dominantCrustType: number;
  seedPosition: Point;
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
  let continentalNeighbors = 0;
  if (!neighbors[cellIdx]) return height;
  for (const nId of neighbors[cellIdx]) {
    if (crustTypes[nId] === 1) continentalNeighbors++;
  }
  const neighborRatio = continentalNeighbors / Math.max(1, neighbors[cellIdx].length);
  const currentType = crustTypes[cellIdx];

  if (currentType === 1 && neighborRatio < 0.5) {
    return height - (1 - neighborRatio) * 0.15;
  } else if (currentType === 0 && neighborRatio > 0.5) {
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
  const vA = cross3(plateA.eulerPole.axis, point);
  const vB = cross3(plateB.eulerPole.axis, point);
  const vAs = scale3(vA, plateA.eulerPole.rate);
  const vBs = scale3(vB, plateB.eulerPole.rate);
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

  // Determine each plate's dominant crust type by initial nearest-seed assignment
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
  const rotatedSeeds = plates.map(p => ({ x: p.seedPosition.x, y: p.seedPosition.y, z: p.seedPosition.z }));

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

    // 3c. Classify boundaries and accumulate uplift
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