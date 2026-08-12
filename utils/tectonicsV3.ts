import { Cell, Point, WorldParams, TectonicResult } from '../types';
import { RNG, SimplexNoise } from './rng';
import { MinHeap } from './pathfinding';
import { seedCrustField } from './crust';
import {
  randomEulerPole, quatFromAxisAngle, quatRotate,
  chordDistance, cross3, dot3, sub3, scale3, magnitude, normalizeVec,
  generateFibonacciSphere,
} from './spherical';

// Noise helpers (duplicated from worldGen.ts to avoid circular imports)
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

// Generate plates with irregular sizes via hierarchical merging.
// Seed 2.5x desired count of proto-plates, apply plateJitter perturbation,
// assign macro-cells, then merge the smallest plates into neighbors.
function generatePlates(numPlates: number, rng: RNG, plateJitter: number): PlateState[] {
  // Remove Euler-pole draws from the RNG that generatePlates used to consume
  // — we generate poles for the final set, not the proto-set.
  const protoCount = Math.ceil(numPlates * 2.5);
  // Proto-plate seed positions from Fibonacci, then jittered
  let protoSeeds = generateFibonacciSphere(protoCount, rng, 0);
  if (plateJitter > 0) {
    protoSeeds = protoSeeds.map(s => {
      const theta = rng.next() * 2 * Math.PI * plateJitter;
      const phi = (rng.next() - 0.5) * Math.PI * plateJitter;
      const q = quatFromAxisAngle(
        { x: Math.cos(theta), y: Math.sin(theta), z: 0 },
        phi
      );
      return quatRotate(q, s);
    });
  }

  // We'll assign macro-cells to proto-plates below (in simulateTectonics),
  // but for merge we just need the seed positions. The actual merge happens
  // after the initial assignment in simulateTectonics.
  const plates: PlateState[] = [];
  for (let i = 0; i < numPlates; i++) {
    plates.push({
      id: i,
      eulerPole: randomEulerPole(rng),
      dominantCrustType: 0,
      seedPosition: protoSeeds[i],
    });
  }
  return plates;
}

function mergeSmallPlates(
  macroNeighbors: number[][],
  plateIds: Int32Array,
  plates: PlateState[],
  minCellThreshold: number,
): void {
  const numPlates = plates.length;
  const numCells = plateIds.length;
  // Count cells per plate
  const counts = new Int32Array(numPlates);
  for (let i = 0; i < numCells; i++) {
    counts[plateIds[i]]++;
  }

  // Repeatedly dissolve the smallest sub-threshold plate. Because Dijkstra
  // assignment makes each plate a single connected region, merging the WHOLE
  // region into its most-common adjacent plate keeps the result connected —
  // the old per-cell nearest-seed reassignment could scatter the dying plate
  // across several non-adjacent neighbors, re-introducing fragments.
  let changed = true;
  while (changed) {
    changed = false;
    let smallestId = -1;
    let smallestCount = Infinity;
    for (let j = 0; j < numPlates; j++) {
      if (counts[j] > 0 && counts[j] < minCellThreshold && counts[j] < smallestCount) {
        smallestCount = counts[j];
        smallestId = j;
      }
    }
    if (smallestId < 0) break;

    // Vote for the most-common adjacent plate across the whole dying region.
    const votes = new Int32Array(numPlates);
    for (let i = 0; i < numCells; i++) {
      if (plateIds[i] !== smallestId) continue;
      for (const nId of macroNeighbors[i]) {
        const np = plateIds[nId];
        if (np !== smallestId && counts[np] > 0) votes[np]++;
      }
    }
    let target = -1;
    let bestVotes = 0;
    for (let j = 0; j < numPlates; j++) {
      if (votes[j] > bestVotes) { bestVotes = votes[j]; target = j; }
    }

    if (target >= 0) {
      for (let i = 0; i < numCells; i++) {
        if (plateIds[i] === smallestId) plateIds[i] = target;
      }
      counts[target] += counts[smallestId];
    }
    counts[smallestId] = 0;
    changed = true;
  }
}

// --- 2. Helper functions ---

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

  // Track the k nearest for each cell without sorting all N.
  // For each pair (i, j), update both cells' top-k arrays.
  const bestDists: Float64Array[] = Array.from({ length: n }, () => new Float64Array(k).fill(Infinity));
  const bestIds: Int32Array[] = Array.from({ length: n }, () => new Int32Array(k).fill(-1));

  for (let i = 0; i < n; i++) {
    const pi = points[i];
    for (let j = i + 1; j < n; j++) {
      const d = chordDistance(pi, points[j]);
      // Update i's top-k
      updateTopK(bestDists[i], bestIds[i], j, d, k);
      // Update j's top-k (symmetric)
      updateTopK(bestDists[j], bestIds[j], i, d, k);
    }
    // Convert to regular array after all candidates processed
    neighbors[i] = Array.from(bestIds[i]).filter(id => id >= 0);
  }

  // Symmetrize: the top-k relation is directional (i may keep j while j's own
  // top-k is full of closer cells), but Dijkstra region-growth needs an
  // undirected graph or plate fronts can leak through one-way edges. Add the
  // reciprocal of every edge, deduped, in ascending index order for
  // determinism.
  const sets: Set<number>[] = neighbors.map(list => new Set(list));
  for (let i = 0; i < n; i++) {
    for (const j of neighbors[i]) sets[j].add(i);
  }
  for (let i = 0; i < n; i++) {
    neighbors[i] = Array.from(sets[i]).sort((a, b) => a - b);
  }
  return neighbors;
}

function updateTopK(dists: Float64Array, ids: Int32Array, candId: number, candDist: number, k: number): void {
  // Find the slot with the largest distance (the "worst" in a top-k)
  let worstIdx = 0;
  let worstDist = dists[0];
  for (let t = 1; t < k; t++) {
    if (dists[t] > worstDist) {
      worstDist = dists[t];
      worstIdx = t;
    }
  }
  // Replace the worst only if the candidate beats it
  if (worstDist > candDist) {
    dists[worstIdx] = candDist;
    ids[worstIdx] = candId;
  }
}

// Static per-edge growth costs for plate region-growing. A plate front pays
// `chordDistance * noiseMul * marginMul` to cross each macro edge:
//   - noiseMul (from boundaryRoughness, sampled at a warp-displaced edge
//     midpoint) makes fronts advance unevenly, so boundaries stop reading as
//     clean Voronoi bisectors. Keeps `boundaryRoughness` AND `warpStrength` live.
//   - marginMul makes crust-type boundaries cheaper to reach, so plate margins
//     are attracted to continental/oceanic transitions. Keeps `marginCoupling` live.
// Returned parallel to `macroNeighbors`: costs[i][t] is the cost of edge
// (i -> macroNeighbors[i][t]). The edge is undirected — both directions share
// the midpoint sample and the symmetric crust-mismatch test.
function computeEdgeCosts(
  macroPoints: Point[],
  macroNeighbors: number[][],
  crustTypes: Uint8Array,
  boundaryRoughness: number,
  warpStrength: number,
  marginCoupling: number,
  warpNoise: SimplexNoise,
  edgeNoise: SimplexNoise,
): number[][] {
  const n = macroPoints.length;
  const warpAmp = warpStrength * 0.2;
  const costs: number[][] = new Array(n);
  for (let i = 0; i < n; i++) {
    const pi = macroPoints[i];
    const list = macroNeighbors[i];
    const row = new Array<number>(list.length);
    for (let t = 0; t < list.length; t++) {
      const nId = list[t];
      const pn = macroPoints[nId];
      // Edge midpoint, warp-displaced so the roughness pattern isn't aligned
      // to the lattice.
      let mx = (pi.x + pn.x) * 0.5;
      let my = (pi.y + pn.y) * 0.5;
      let mz = (pi.z + pn.z) * 0.5;
      if (warpAmp > 0) {
        mx += warpNoise.noise3D(mx * 0.5, my * 0.5, mz * 0.5) * warpAmp;
        my += warpNoise.noise3D(my * 0.5, mz * 0.5, mx * 0.5) * warpAmp;
        mz += warpNoise.noise3D(mz * 0.5, mx * 0.5, my * 0.5) * warpAmp;
      }
      const noiseMul = boundaryRoughness > 0
        ? Math.max(0.15, 1 + boundaryRoughness * 0.8 * edgeNoise.noise3D(mx * 4, my * 4, mz * 4))
        : 1;
      const marginMul = crustTypes[i] !== crustTypes[nId]
        ? Math.max(0.2, 1 - marginCoupling * 0.5)
        : 1;
      row[t] = chordDistance(pi, pn) * noiseMul * marginMul;
    }
    costs[i] = row;
  }
  return costs;
}

// Multi-source Dijkstra region-growth. Every plate grows outward from the macro
// cell nearest its (rotated) seed, following graph edges — so each plate's
// territory is a single connected region BY CONSTRUCTION. This is what kills
// the enclave/exclave artifacts the old per-cell argmin-with-noise produced:
// there is no way for a cell to end up owned by a plate it has no path to.
// `plateSpeeds[j]` scales plate j's growth (faster plates claim more territory),
// giving the size irregularity the old proto-plate-merge scheme was meant to.
function assignPlatesDijkstra(
  macroPoints: Point[],
  macroNeighbors: number[][],
  edgeCosts: number[][],
  rotatedSeeds: Point[],
  plateSpeeds: number[],
  activeMask: (j: number) => boolean,
  plateIds: Int32Array,
): void {
  const n = macroPoints.length;
  const dist = new Float64Array(n).fill(Infinity);
  plateIds.fill(-1);

  interface Node { cell: number; plate: number; dist: number; }
  // Score baked with tiny index/plate terms so pop order is fully deterministic
  // on distance ties (identical seeds must give identical assignments).
  const heap = new MinHeap<Node>(node => node.dist + node.cell * 1e-9 + node.plate * 1e-12);

  // Seed each active plate at its nearest macro cell.
  for (let j = 0; j < rotatedSeeds.length; j++) {
    if (!activeMask(j)) continue;
    let best = -1;
    let bestD = Infinity;
    const s = rotatedSeeds[j];
    for (let i = 0; i < n; i++) {
      const d = chordDistance(macroPoints[i], s);
      if (d < bestD) { bestD = d; best = i; }
    }
    if (best >= 0) heap.push({ cell: best, plate: j, dist: 0 });
  }

  while (heap.size() > 0) {
    const { cell, plate, dist: d } = heap.pop()!;
    if (plateIds[cell] !== -1) continue; // already settled by a cheaper front
    plateIds[cell] = plate;
    dist[cell] = d;
    const speed = plateSpeeds[plate];
    const list = macroNeighbors[cell];
    const costRow = edgeCosts[cell];
    for (let t = 0; t < list.length; t++) {
      const nId = list[t];
      if (plateIds[nId] !== -1) continue;
      const nd = d + costRow[t] / speed;
      if (nd < dist[nId]) {
        dist[nId] = nd;
        heap.push({ cell: nId, plate, dist: nd });
      }
    }
  }

  // Fallback for any cell unreachable through the graph (isolated component):
  // hand it to the nearest active seed so nothing stays -1.
  for (let i = 0; i < n; i++) {
    if (plateIds[i] !== -1) continue;
    let best = 0;
    let bestD = Infinity;
    for (let j = 0; j < rotatedSeeds.length; j++) {
      if (!activeMask(j)) continue;
      const d = chordDistance(macroPoints[i], rotatedSeeds[j]);
      if (d < bestD) { bestD = d; best = j; }
    }
    plateIds[i] = best;
  }
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
  const plateJitter = params.plateJitter ?? 0.3;

  if (signal?.aborted) throw new Error('Generation Cancelled');

  // 1. Seed crust field (independent of plates)
  const crust = seedCrustField(macroPoints, params, simplex, crustRng);

  // 2. Generate plates with Euler poles (proto-plates + merge for irregular sizes)
  const plates = generatePlates(numPlates, plateRng, plateJitter);

  // Macro-cell neighbor graph (built once; drives BOTH region growth and the
  // boundary classification below). Hoisted ahead of assignment because
  // Dijkstra needs it.
  onLog?.("V3: Building macro-cell neighbor graph...");
  const macroNeighbors = buildMacroNeighborGraph(macroPoints, 6);

  // Noise streams: warpNoise displaces the edge-midpoint samples so roughness
  // isn't lattice-aligned; edgeNoise is an independent side-stream for the
  // per-edge roughness cost.
  const warpNoise = new SimplexNoise(new RNG(params.seed + '_warp_v3'));
  const edgeNoise = new SimplexNoise(new RNG(params.seed + '_edge_v3'));
  const boundaryRoughness = params.boundaryRoughness ?? 0.3;
  const warpStrength = params.warpStrength ?? 0.5;

  // Static per-edge growth costs — independent of plate motion, so computed
  // once and reused every timestep.
  const edgeCosts = computeEdgeCosts(
    macroPoints, macroNeighbors, crust.crustTypes,
    boundaryRoughness, warpStrength, marginCoupling, warpNoise, edgeNoise,
  );

  // Per-plate growth speed ∈ [0.75, 1.3] — faster plates claim more territory,
  // giving the power-law size spread the old proto-plate merge approximated.
  const plateSpeeds = plates.map(() => plateRng.range(0.75, 1.3));

  // Rotated plate seed positions (mutated each step; start at the raw seeds).
  const rotatedSeeds = plates.map(p => ({ x: p.seedPosition.x, y: p.seedPosition.y, z: p.seedPosition.z }));

  // 3. Initial plate assignment (Dijkstra region-growth) + merge of small plates
  const plateIds = new Int32Array(numMacro);
  assignPlatesDijkstra(
    macroPoints, macroNeighbors, edgeCosts, rotatedSeeds, plateSpeeds,
    () => true, plateIds,
  );

  // Merge plates below 0.5% cell threshold
  const minCellThreshold = Math.max(1, Math.floor(numMacro * 0.005));
  mergeSmallPlates(macroNeighbors, plateIds, plates, minCellThreshold);

  // Re-tally plate counts after merge, and set dominant crust types
  const postMergeCounts = new Int32Array(plates.length);
  for (let i = 0; i < numMacro; i++) {
    postMergeCounts[plateIds[i]]++;
  }
  for (let j = 0; j < plates.length; j++) {
    let continentalCount = 0;
    for (let i = 0; i < numMacro; i++) {
      if (plateIds[i] === j) {
        if (crust.crustTypes[i] === 1) continentalCount++;
      }
    }
    plates[j].dominantCrustType = postMergeCounts[j] > 0 && continentalCount / Math.max(1, postMergeCounts[j]) > 0.5 ? 1 : 0;
  }

  const activePlateCount = postMergeCounts.filter(c => c > 0).length;
  onLog?.(`V3: ${activePlateCount} active plates after merging`);

  // Per-macro-cell accumulation buffers
  const upliftAccum = new Float32Array(numMacro);
  const thickness = new Float32Array(numMacro);
  for (let i = 0; i < numMacro; i++) thickness[i] = crust.crustThickness[i];

  // 4. Timestep loop
  for (let step = 0; step < numSteps; step++) {
    if (signal?.aborted) throw new Error('Generation Cancelled');
    if (step % 5 === 0) onLog?.(`V3: Timestep ${step + 1}/${numSteps}...`);

    // 4a. Rotate each active plate seed by its Euler pole
    for (let i = 0; i < plates.length; i++) {
      if (postMergeCounts[i] === 0) continue;
      const q = quatFromAxisAngle(plates[i].eulerPole.axis, plates[i].eulerPole.rate);
      rotatedSeeds[i] = quatRotate(q, rotatedSeeds[i]);
    }

    // 4b. Re-grow plate regions from the rotated seeds. Region-growth over the
    // macro graph keeps every plate connected each step — no argmin lottery, so
    // no enclaves/exclaves. Roughness and margin attraction live in the static
    // edge costs; plate motion enters only through the moving seeds.
    assignPlatesDijkstra(
      macroPoints, macroNeighbors, edgeCosts, rotatedSeeds, plateSpeeds,
      j => postMergeCounts[j] > 0, plateIds,
    );

    // 4c. Classify boundaries and accumulate uplift — velocity-scaled, asymmetric
    for (let i = 0; i < numMacro; i++) {
      const p = macroPoints[i];
      let maxCollisionUplift = 0;
      let maxSubductionUplift = 0;
      let maxSubductionTrench = 0;
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
        const vnMag = Math.abs(dot3(relV, edgeNormal));
        const vn = dot3(relV, edgeNormal);
        const vtMag = magnitude(sub3(relV, scale3(edgeNormal, vn)));

        const crustA = crust.crustTypes[i];
        const crustB = crust.crustTypes[nId];

        if (vn < -0.0005) {
          // Convergent
          if (crustA === 1 && crustB === 1) {
            // Continental collision — broad, massive
            if (vnMag > maxCollisionUplift) maxCollisionUplift = vnMag;
          } else if (crustA === 0 && crustB === 1) {
            // This cell is oceanic, neighbor is continental — this side subducts
            if (vnMag > maxSubductionTrench) maxSubductionTrench = vnMag;
          } else if (crustA === 1 && crustB === 0) {
            // This cell is continental, neighbor is oceanic — arc on this side
            if (vnMag > maxSubductionUplift) maxSubductionUplift = vnMag;
          } else {
            // Oceanic-oceanic — symmetric trench + arc
            if (vnMag > maxSubductionTrench) maxSubductionTrench = vnMag;
            if (vnMag > maxSubductionUplift) maxSubductionUplift = vnMag;
          }
        } else if (vn > 0.0005) {
          // Divergent — rift (no positive uplift contribution)
          const riftAmount = vn * tectonicStrength * 10;
          upliftAccum[i] -= riftAmount;
          thickness[i] = Math.max(0, thickness[i] - riftAmount);
        }

        if (vtMag > maxShear) maxShear = vtMag;
      }

      // Per-boundary noise modulation (segmented mountain belts)
      const noiseVal = simplex.noise3D(p.x * 2, p.y * 2, p.z * 2);
      const noiseFactor = 0.3 + (noiseVal * 0.5 + 0.5) * 0.7;

      // Apply uplift — velocity-scaled, no smoothstep saturation
      const collisionUplift = maxCollisionUplift * tectonicStrength * 60 * noiseFactor;
      upliftAccum[i] += collisionUplift;
      if (collisionUplift > 0) thickness[i] += collisionUplift * 1.5;

      const subductionArcUplift = maxSubductionUplift * tectonicStrength * 30 * noiseFactor;
      upliftAccum[i] += subductionArcUplift;
      if (subductionArcUplift > 0) thickness[i] += subductionArcUplift * 2;

      const subductionTrench = maxSubductionTrench * tectonicStrength * 20 * noiseFactor;
      upliftAccum[i] -= subductionTrench;

      const shearUplift = maxShear * tectonicStrength * 5 * noiseFactor;
      upliftAccum[i] += shearUplift;
    }
  }

  // 5. Compose final height from crust thickness + uplift + isostasy
  onLog?.("V3: Composing tectonic heights...");
  const heights = new Float32Array(numMacro);
  for (let i = 0; i < numMacro; i++) {
    const isostaticElevation = saturatedIsostasy(thickness[i], crust.crustTypes[i] === 1);
    let h = isostaticElevation + upliftAccum[i] * 0.5;
    h = applyPassiveMargin(h, crust.crustTypes, i, macroNeighbors);
    heights[i] = h;
  }

  // Count how many plates actually have cells assigned for logging
  const plateCounts = new Map<number, number>();
  for (let i = 0; i < numMacro; i++) {
    plateCounts.set(plateIds[i], (plateCounts.get(plateIds[i]) ?? 0) + 1);
  }
  onLog?.(`V3: ${plateCounts.size} plates with macro-cells assigned after simulation`);

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

  return { heights, crustTypes: crust.crustTypes, crustThickness: thickness, upliftAccum, plateIds };
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
    dc.plateId = macroResult.plateIds[nearest];

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