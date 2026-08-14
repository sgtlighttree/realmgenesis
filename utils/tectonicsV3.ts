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

// Multi-source Dijkstra region-growth. Every plate grows outward from its
// seed chain (a connected run of macro cells — see the `elongation` seeding
// below), following graph edges — so each plate's MACRO-cell territory is a
// single connected region BY CONSTRUCTION. This is what kills the enclave/
// exclave artifacts the old per-cell argmin-with-noise produced at the macro
// level: there is no way for a macro cell to end up owned by a plate it has
// no path to. NOTE: this guarantee is at the macro-graph level only — the
// nearest-macro-cell downsample onto the fine mesh (`dc.plateId =
// macroResult.plateIds[nearest]`) can still pinch a connected macro region
// into disconnected fine-mesh cells; measured pre-existing at elongation 0.0
// too (see tests/plateConnectivity.test.ts). `plateSpeeds[j]` scales plate
// j's growth (faster plates claim more territory), giving the size
// irregularity the old proto-plate-merge scheme was meant to.
function assignPlatesDijkstra(
  macroPoints: Point[],
  macroNeighbors: number[][],
  edgeCosts: number[][],
  rotatedSeeds: Point[],
  plateSpeeds: number[],
  activeMask: (j: number) => boolean,
  plateIds: Int32Array,
  plates: PlateState[],
  elongation: number,
): void {
  const n = macroPoints.length;
  const dist = new Float64Array(n).fill(Infinity);
  plateIds.fill(-1);

  interface Node { cell: number; plate: number; dist: number; }
  // Score baked with tiny index/plate terms so pop order is fully deterministic
  // on distance ties (identical seeds must give identical assignments).
  const heap = new MinHeap<Node>(node => node.dist + node.cell * 1e-9 + node.plate * 1e-12);

  // Seed each active plate at a chain of macro cells grown from the nearest
  // unclaimed cell along the plate's velocity direction — elongated seeds
  // produce elongated plates instead of round Voronoi blobs.
  const chainLen = 1 + Math.round(elongation * 4);
  // Cells already used as a dist-0 source by an EARLIER plate's chain. Disjoint
  // source sets across plates are what keep connectivity by construction: two
  // plates must never both seed the same cell, or one plate's chain gets severed
  // into an exclave. A dist-0 pop always beats any positive-dist front, so every
  // unclaimed source settles for its own plate.
  const claimed = new Uint8Array(n);
  for (let j = 0; j < rotatedSeeds.length; j++) {
    if (!activeMask(j)) continue;
    // nearest UNCLAIMED macro cell to the (rotated) seed
    let best = -1, bestD = Infinity;
    const s = rotatedSeeds[j];
    for (let i = 0; i < n; i++) {
      if (claimed[i]) continue;
      const d = chordDistance(macroPoints[i], s);
      if (d < bestD) { bestD = d; best = i; }
    }
    if (best < 0) continue;

    // chain axis = plate velocity direction at the seed (tangent to sphere)
    const vel = scale3(cross3(plates[j].eulerPole.axis, s), plates[j].eulerPole.rate);
    const axis = magnitude(vel) > 1e-9 ? normalizeVec(vel) : null;

    // walk a connected chain of UNCLAIMED cells forward along the axis; all become dist-0 sources
    claimed[best] = 1;
    heap.push({ cell: best, plate: j, dist: 0 });
    let tip = best;
    let chainSize = 1;
    while (axis && chainSize < chainLen) {
      let nextCell = -1, nextDot = 0; // require forward progress (dot > 0)
      for (const nId of macroNeighbors[tip]) {
        if (claimed[nId]) continue;
        const dir = normalizeVec({
          x: macroPoints[nId].x - macroPoints[tip].x,
          y: macroPoints[nId].y - macroPoints[tip].y,
          z: macroPoints[nId].z - macroPoints[tip].z,
        });
        const d = dot3(dir, axis);
        if (d > nextDot || (d === nextDot && nextCell >= 0 && nId < nextCell)) { nextDot = d; nextCell = nId; }
      }
      if (nextCell < 0) break;
      claimed[nextCell] = 1;
      heap.push({ cell: nextCell, plate: j, dist: 0 });
      tip = nextCell;
      chainSize++;
    }
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

// GDH1 seafloor depth (Stein & Stein 1992), meters, as a function of crust age
// in Ma. Ridge crest ~2500 m, flattening toward ~5651 m for old floor.
function gdh1Depth(ageMa: number): number {
  const t = Math.max(0, ageMa);
  return t < 20 ? 2500 + 350 * Math.sqrt(t) : 5651 - 2473 * Math.exp(-0.0278 * t);
}

// Map a GDH1 depth (m) into the engine's oceanic height band: ridge ≈ −0.5,
// oldest floor ≈ −0.85. Kept inside the existing ocean band so the global
// normalization and the Stage-9b oceanDepth remap keep behaving.
function depthToBandHeight(depthM: number): number {
  const f = Math.min(1, Math.max(0, (depthM - 2500) / (5651 - 2500)));
  return -0.5 - f * 0.35;
}

// Seafloor age via multi-source Dijkstra distance-to-ridge over OCEANIC macro
// cells only (continents block propagation). Ridges = oceanic cells on a
// divergent boundary in the final state. age = distance / spreadRate, capped.
// Returns Ma per cell; −1 for continental/non-ocean. Empty-ridge worlds return
// all −1 so the caller falls back to isostasy.
const MAX_SEAFLOOR_AGE = 180; // Ma; oldest surviving oceanic crust on Earth
function computeSeafloorAge(
  macroPoints: Point[],
  macroNeighbors: number[][],
  plateIds: Int32Array,
  crustTypes: Uint8Array,
  plates: PlateState[],
  rotatedSeeds: Point[],
  spreadRate: number,
  seed: string,
): Float32Array {
  const ageNoise = new SimplexNoise(new RNG(seed + '_agenoise_v3'));
  const n = macroPoints.length;
  const age = new Float32Array(n).fill(-1);
  const dist = new Float64Array(n).fill(Infinity);
  const heap = new MinHeap<{ cell: number; dist: number }>(x => x.dist + x.cell * 1e-9);

  for (let i = 0; i < n; i++) {
    if (crustTypes[i] !== 0) continue; // oceanic only
    const p = macroPoints[i];
    for (const nId of macroNeighbors[i]) {
      if (plateIds[nId] === plateIds[i]) continue;
      const relV = computeRelativeVelocity(plates[plateIds[i]], plates[plateIds[nId]], p);
      const edgeNormal = normalizeVec({
        x: macroPoints[nId].x - p.x, y: macroPoints[nId].y - p.y, z: macroPoints[nId].z - p.z,
      });
      if (dot3(relV, edgeNormal) > 0.0005) { // divergent → ridge
        dist[i] = 0;
        heap.push({ cell: i, dist: 0 });
        break;
      }
    }
  }
  if (heap.size() === 0) return age; // no ridges (e.g. Pangea) → all −1, isostasy fallback

  while (heap.size() > 0) {
    const { cell, dist: d } = heap.pop()!;
    if (d > dist[cell]) continue;
    for (const nId of macroNeighbors[cell]) {
      if (crustTypes[nId] !== 0) continue; // stay in the ocean
      const nd = d + chordDistance(macroPoints[cell], macroPoints[nId]);
      if (nd < dist[nId]) { dist[nId] = nd; heap.push({ cell: nId, dist: nd }); }
    }
  }

  const rate = Math.max(1e-4, spreadRate);
  for (let i = 0; i < n; i++) {
    if (crustTypes[i] !== 0) continue;
    // Unreached oceanic basin (ringed by land): saturate to max age.
    const d = dist[i] === Infinity ? MAX_SEAFLOOR_AGE * rate : dist[i];
    const raw = d / rate;
    const p = macroPoints[i];
    const perturbed = raw * (1 + 0.1 * ageNoise.noise3D(p.x * 2, p.y * 2, p.z * 2));
    age[i] = Math.max(0, Math.min(MAX_SEAFLOOR_AGE, perturbed));
  }
  return age;
}

// Microplate injection (D7 part 2, Goal 1). Instead of peeling cells post-hoc
// (which desyncs plate lookups), seed new small plates at high-shear boundary
// segments and let the existing per-timestep Dijkstra re-assignment grow them
// as connected regions — so the 0-exclave invariant holds by construction.
// Draws only from `microRng` so plateRng's draw order (determinism) is untouched.
function injectMicroplates(
  macroPoints: Point[],
  macroNeighbors: number[][],
  plateIds: Int32Array,
  plates: PlateState[],
  plateSpeeds: number[],
  rotatedSeeds: Point[],
  crustTypes: Uint8Array,
  intensity: number,
  microRng: RNG,
): void {
  if (intensity <= 0) return;
  const n = macroPoints.length;
  // Shear magnitude at each boundary cell (tangential relative velocity).
  const shear = new Float64Array(n).fill(-1);
  for (let i = 0; i < n; i++) {
    const p = macroPoints[i];
    let maxShear = 0;
    for (const nId of macroNeighbors[i]) {
      if (plateIds[nId] === plateIds[i]) continue;
      const relV = computeRelativeVelocity(plates[plateIds[i]], plates[plateIds[nId]], p);
      const edgeNormal = normalizeVec({
        x: macroPoints[nId].x - p.x, y: macroPoints[nId].y - p.y, z: macroPoints[nId].z - p.z,
      });
      const vn = dot3(relV, edgeNormal);
      const vt = magnitude(sub3(relV, scale3(edgeNormal, vn)));
      if (vt > maxShear) maxShear = vt;
    }
    if (maxShear > 0) shear[i] = maxShear;
  }
  // Rank boundary cells by shear, take the strongest, spaced apart so microplates
  // don't clump. Count scales with intensity and plate count.
  const ranked = [...Array(n).keys()].filter(i => shear[i] > 0).sort((a, b) =>
    shear[b] - shear[a] || a - b);
  const target = Math.round(intensity * Math.max(2, plates.length * 0.6));
  const chosen: number[] = [];
  const minSepSq = 0.25 * 0.25; // ~chord separation between microplate seeds
  for (const c of ranked) {
    if (chosen.length >= target) break;
    const pc = macroPoints[c];
    if (chosen.some(o => {
      const q = macroPoints[o];
      return (pc.x - q.x) ** 2 + (pc.y - q.y) ** 2 + (pc.z - q.z) ** 2 < minSepSq;
    })) continue;
    chosen.push(c);
  }
  for (const c of chosen) {
    const micro: PlateState = {
      id: plates.length,
      eulerPole: randomEulerPole(microRng),
      dominantCrustType: crustTypes[c],
      seedPosition: { x: macroPoints[c].x, y: macroPoints[c].y, z: macroPoints[c].z },
    };
    plates.push(micro);
    plateSpeeds.push(microRng.range(0.4, 0.7)); // slow → small
    rotatedSeeds.push({ x: macroPoints[c].x, y: macroPoints[c].y, z: macroPoints[c].z });
    plateIds[c] = micro.id; // plant one cell so it's active; re-assignment grows it
  }
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
  const spreadRate = params.spreadRate ?? 0.008;
  const microplateIntensity = params.microplateIntensity ?? 0;

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
  const plateElongation = params.plateElongation ?? 0.4;

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
    () => true, plateIds, plates, plateElongation,
  );

  // Merge plates below 0.5% cell threshold
  const minCellThreshold = Math.max(1, Math.floor(numMacro * 0.005));
  mergeSmallPlates(macroNeighbors, plateIds, plates, minCellThreshold);

  // D7 part 2: inject shear-driven microplates AFTER merge (so the 0.5% cutoff
  // doesn't eat them). They grow to connected regions in the timestep loop.
  const microRng = new RNG(params.seed + '_micro_v3');
  injectMicroplates(
    macroPoints, macroNeighbors, plateIds, plates, plateSpeeds, rotatedSeeds,
    crust.crustTypes, microplateIntensity, microRng,
  );

  // Re-tally plate counts after merge + microplate injection, and set dominant
  // crust types (sizes to the possibly-grown plates array).
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
      j => postMergeCounts[j] > 0, plateIds, plates, plateElongation,
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
        } else if (vn > 0.0005 && crustA === 1 && crustB === 1) {
          // Divergent CONTINENTAL rift only — a rift valley. Oceanic divergence
          // is a spreading ridge; its elevation comes from GDH1 age→depth below,
          // so lowering it here would fight the bathymetry model.
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

  // 4d. Seafloor age → bathymetry (D7 part 2, Goal 2). Final-state ridges only.
  onLog?.("V3: Computing seafloor age (bathymetry)...");
  const seafloorAge = computeSeafloorAge(
    macroPoints, macroNeighbors, plateIds, crust.crustTypes, plates, rotatedSeeds, spreadRate, params.seed,
  );
  // Abyssal-hill texture amplitude, baked at the former seafloorDetail default
  // (0.5). The Seafloor Detail slider was repurposed into the seafloorDepth datum
  // (Stage 9b, worldGen.ts); this internal texture is no longer user-exposed.
  const seafloorDetail = 0.5;
  const abyssalNoise = new SimplexNoise(new RNG(params.seed + '_abyssal_v3'));

  // 5. Compose final height. Oceanic floor with a valid age follows GDH1
  // subsidence (deeper away from ridges) instead of flat isostasy; continents
  // and ridgeless basins keep isostasy. Uplift (subduction trenches/arcs) still
  // adds on top so convergent features survive.
  onLog?.("V3: Composing tectonic heights...");
  const heights = new Float32Array(numMacro);
  for (let i = 0; i < numMacro; i++) {
    let base: number;
    if (crust.crustTypes[i] === 0 && seafloorAge[i] >= 0) {
      const p = macroPoints[i];
      const hill = abyssalNoise.noise3D(p.x * 6, p.y * 6, p.z * 6) * seafloorDetail * 0.06;
      base = depthToBandHeight(gdh1Depth(seafloorAge[i])) + hill;
    } else {
      base = saturatedIsostasy(thickness[i], crust.crustTypes[i] === 1);
    }
    let h = base + upliftAccum[i] * 0.5;
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
    // Roughness scales sub-cell structural relief. Centered so the default
    // (0.5) is ×1.0 — seeds at default roughness stay byte-identical while the
    // Terrain Roughness slider stays meaningful under V3 (it drove V2's fBm
    // amplitude; here it drives the display-resolution relief amplitude).
    const roughness = params.roughness ?? 0.5;
    const structuralNoise = (fbmVal * (1 - blend) + ridgedRemapped * blend) * (0.5 + roughness);

    // Deep-ocean cells carry GDH1 bathymetry; damp the structural noise there
    // so the age→depth gradient isn't washed out. Damping baked at the former
    // seafloorDetail default (0.5 → factor 0.675); the slider is now seafloorDepth.
    let noiseInfluence = 1.2 - tectonicStrength;
    if (dc.crustType === 0 && macroResult.heights[nearest] < 0.5) {
      noiseInfluence *= 1 - 0.65 * 0.5;
    }
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