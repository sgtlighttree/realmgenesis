import { geoVoronoi } from 'd3-geo-voronoi';
import { Cell, Point, WorldData, WorldParams, BiomeType, CivData, FactionData, ProvinceData, TownData } from '../types';
import { RNG, SimplexNoise } from './rng';
import { BIOME_COLORS } from './colors';

// --- DATA STRUCTURES ---

class MinHeap<T> {
    private heap: T[];
    private scoreFunction: (t: T) => number;

    constructor(scoreFunction: (t: T) => number) {
        this.heap = [];
        this.scoreFunction = scoreFunction;
    }

    push(node: T) {
        this.heap.push(node);
        this.bubbleUp(this.heap.length - 1);
    }

    pop(): T | undefined {
        if (this.heap.length === 0) return undefined;
        const top = this.heap[0];
        const bottom = this.heap.pop();
        if (this.heap.length > 0 && bottom !== undefined) {
            this.heap[0] = bottom;
            this.sinkDown(0);
        }
        return top;
    }

    size(): number { return this.heap.length; }

    private bubbleUp(index: number) {
        while (index > 0) {
            const parentIndex = Math.floor((index - 1) / 2);
            if (this.scoreFunction(this.heap[index]) >= this.scoreFunction(this.heap[parentIndex])) break;
            [this.heap[index], this.heap[parentIndex]] = [this.heap[parentIndex], this.heap[index]];
            index = parentIndex;
        }
    }

    private sinkDown(index: number) {
        const length = this.heap.length;
        const element = this.heap[index];
        const elemScore = this.scoreFunction(element);

        while (true) {
            let leftChildIdx = 2 * index + 1;
            let rightChildIdx = 2 * index + 2;
            let leftScore, rightScore;
            let swap = null;

            if (leftChildIdx < length) {
                leftScore = this.scoreFunction(this.heap[leftChildIdx]);
                if (leftScore < elemScore) swap = leftChildIdx;
            }
            if (rightChildIdx < length) {
                rightScore = this.scoreFunction(this.heap[rightChildIdx]);
                if (swap === null) {
                    if (rightScore < elemScore) swap = rightChildIdx;
                } else {
                    if (rightScore < leftScore!) swap = rightChildIdx;
                }
            }

            if (swap === null) break;
            [this.heap[index], this.heap[swap]] = [this.heap[swap], this.heap[index]];
            index = swap;
        }
    }
}

// --- MATH HELPERS ---

function toSpherical(x: number, y: number, z: number): [number, number] {
  const r = Math.sqrt(x * x + y * y + z * z);
  if (r === 0) return [0, 0];
  let lat = Math.asin(Math.max(-1, Math.min(1, y / r))) * (180 / Math.PI);
  let lon = Math.atan2(z, x) * (180 / Math.PI);
  return [lat, lon];
}

function generateFibonacciSphere(samples: number, rng: RNG, jitter: number): Point[] {
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

function randomVector(rng: RNG): Point {
    const u = rng.next();
    const v = rng.next();
    const theta = 2 * Math.PI * u;
    const phi = Math.acos(2 * v - 1);
    return {
        x: Math.sin(phi) * Math.cos(theta),
        y: Math.sin(phi) * Math.sin(theta),
        z: Math.cos(phi)
    };
}

// --- NOISE ALGORITHMS ---

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

// --- EROSION ---

// Helper to check for abort signal
const checkAbort = (signal?: AbortSignal) => {
    if (signal?.aborted) {
        throw new Error("Generation Cancelled");
    }
};

async function applyHydraulicErosion(cells: Cell[], iterations: number, seaLevel: number, signal?: AbortSignal): Promise<void> {
    cells.forEach(c => c.flux = 0);
    const sorted = [...cells].sort((a, b) => b.height - a.height);
    const erosionRate = 0.02;
    const depositionRate = 0.01;
    const rainAmount = 0.1;

    // Yield every few iterations to keep UI responsive
    const chunkSize = 5;

    for (let iter = 0; iter < iterations; iter++) {
        if (iter % chunkSize === 0) {
            await new Promise(r => setTimeout(r, 0));
            checkAbort(signal);
        }

        // Only rain on land
        sorted.forEach(c => c.flux = c.height >= seaLevel ? rainAmount : 0);
        
        sorted.forEach(c => {
            let lowestH = c.height;
            let targetId = -1;
            for (const nId of c.neighbors) {
                const n = cells[nId];
                if (n.height < lowestH) {
                    lowestH = n.height;
                    targetId = nId;
                }
            }
            if (targetId !== -1) {
                const target = cells[targetId];
                target.flux! += c.flux!;
                const slope = c.height - lowestH;
                const streamPower = c.flux! * slope * 10; 
                const erosion = streamPower * erosionRate;
                const safeErosion = Math.min(erosion, slope * 0.9);
                c.height -= safeErosion;
                target.height += safeErosion * depositionRate; 
            }
        });
    }
}

async function applyThermalErosion(cells: Cell[], iterations: number, signal?: AbortSignal) {
    const talus = 0.008; // Min slope diff
    const rate = 0.2; 
    const chunkSize = 5;

    for(let iter=0; iter<iterations; iter++) {
        if (iter % chunkSize === 0) {
            await new Promise(r => setTimeout(r, 0));
            checkAbort(signal);
        }
        cells.forEach(c => {
            let maxDiff = 0;
            let lowestNIndex = -1;
            for(const nId of c.neighbors) {
                const diff = c.height - cells[nId].height;
                if (diff > maxDiff) {
                    maxDiff = diff;
                    lowestNIndex = nId;
                }
            }
            if (maxDiff > talus && lowestNIndex !== -1) {
                const transfer = (maxDiff - talus) * rate;
                c.height -= transfer;
                cells[lowestNIndex].height += transfer;
            }
        });
    }
}

// --- RIVER GENERATION ---

async function generateRivers(cells: Cell[], seaLevel: number, params: WorldParams, onProgress?: (msg: string) => void, signal?: AbortSignal): Promise<Point[][]> {
    const numCells = cells.length;
    onProgress?.("Rivers: Initializing drainage map...");
    await new Promise(r => setTimeout(r, 0));
    checkAbort(signal);
    
    // 1. Depression Filling (Drainage Enforcement)
    // CRITICAL FIX: Use Float64Array to prevent infinite loops caused by precision mismatch 
    // between 32-bit storage and 64-bit JS calculations.
    const waterLevel = new Float64Array(numCells).fill(Infinity);
    const downstream = new Int32Array(numCells).fill(-1);
    
    const heap = new MinHeap<{id: number, lvl: number}>(x => x.lvl);
    
    let oceanCount = 0;
    cells.forEach(c => {
        if (c.height < seaLevel) {
            waterLevel[c.id] = c.height;
            heap.push({id: c.id, lvl: c.height});
            oceanCount++;
        }
    });

    if (oceanCount === 0) {
        onProgress?.("Warning: No ocean found. River generation skipped.");
        return [];
    }

    let processed = 0;
    
    onProgress?.(`Rivers: Propagating water levels...`);

    while(heap.size() > 0) {
        // Safety break and log update
        if (++processed % 2000 === 0) {
            await new Promise(r => setTimeout(r, 0));
            checkAbort(signal);
            onProgress?.(`Rivers: Drainage processed ${processed} cells...`);
        }

        const {id, lvl} = heap.pop()!;
        
        if (lvl > waterLevel[id]) continue; 

        const c = cells[id];
        for(const nId of c.neighbors) {
            const n = cells[nId];
            const targetLvl = Math.max(n.height, lvl);
            
            if (targetLvl < waterLevel[nId]) {
                waterLevel[nId] = targetLvl;
                downstream[nId] = id; 
                heap.push({id: nId, lvl: targetLvl});
            }
        }
    }

    onProgress?.("Rivers: Accumulating flux...");
    await new Promise(r => setTimeout(r, 0));
    checkAbort(signal);

    // 2. Accumulate Flux
    const sortedIndices = Array.from({length: numCells}, (_, i) => i)
                               .sort((a,b) => waterLevel[b] - waterLevel[a]);
    
    const flux = new Float32Array(numCells).fill(0);
    
    for(const idx of sortedIndices) {
        const c = cells[idx];
        if (c.height < seaLevel) continue;
        const precip = c.moisture * (params.rainfallMultiplier || 1.0);
        flux[idx] += precip;
        const target = downstream[idx];
        if (target !== -1) flux[target] += flux[idx];
    }
    cells.forEach((c, i) => c.flux = flux[i]);

    onProgress?.("Rivers: Tracing paths...");
    // 3. Trace Rivers
    const threshold = 1.0; 
    const visited = new Set<number>();
    const riverPaths: Point[][] = [];
    
    const candidates = sortedIndices.filter(i => flux[i] > threshold && cells[i].height >= seaLevel);

    const getRenderPoint = (c: Cell) => {
        const r = 1 + c.height * 0.05 + 0.005;
        return { x: c.center.x * r, y: c.center.y * r, z: c.center.z * r };
    };

    processed = 0;
    for (const startId of candidates) {
        if (visited.has(startId)) continue;
        if (++processed % 500 === 0) {
             await new Promise(r => setTimeout(r, 0));
             checkAbort(signal);
        }
        
        const path: Point[] = [];
        let curr = startId;
        let safety = 0;
        
        while(curr !== -1 && safety++ < 2000) {
            path.push(getRenderPoint(cells[curr]));
            visited.add(curr);
            
            const next = downstream[curr];
            if (next === -1) break; 
            
            if (cells[next].height < seaLevel) {
                path.push(getRenderPoint(cells[next]));
                break;
            }
            if (visited.has(next)) {
                path.push(getRenderPoint(cells[next]));
                break;
            }
            curr = next;
        }
        
        if (path.length >= 2) riverPaths.push(path);
    }
    
    onProgress?.(`Rivers: Generated ${riverPaths.length} segments.`);
    return riverPaths;
}

// --- BIOME ---

function determineBiome(height: number, temp: number, moisture: number, seaLevel: number): BiomeType {
  if (height < seaLevel) {
    if (height < seaLevel * 0.6) return BiomeType.DEEP_OCEAN;
    return BiomeType.OCEAN;
  }
  const landH = (height - seaLevel) / (1 - seaLevel);
  if (landH < 0.02 && temp > 15) return BiomeType.BEACH;
  if (landH > 0.85 && temp > -5) return BiomeType.VOLCANIC;
  if (temp < -10) return BiomeType.ICE_CAP;
  if (temp < 0) return BiomeType.TUNDRA;
  
  if (moisture < 0.15) {
      if (temp > 25) return BiomeType.HOT_DESERT;
      if (temp > 10) return BiomeType.STEPPE;
      return BiomeType.COLD_DESERT;
  }
  if (moisture < 0.4) {
      if (temp > 25) return BiomeType.TROPICAL_SAVANNA;
      if (temp > 10) return BiomeType.MEDITERRANEAN;
      return BiomeType.STEPPE;
  }
  
  if (temp > 25) return BiomeType.TROPICAL_RAINFOREST;
  if (temp > 15) return BiomeType.TEMPERATE_RAINFOREST;
  if (temp > 5) return BiomeType.TEMPERATE_FOREST;
  return BiomeType.BOREAL_FOREST;
}

// --- TECTONIC HELPERS ---

function enforceConnectivity(cells: Cell[], numPlates: number) {
    const componentId = new Int32Array(cells.length).fill(-1);
    const compSize: number[] = [];
    const compPlate: number[] = [];
    
    let compCount = 0;
    
    for(let i=0; i<cells.length; i++) {
        if(componentId[i] !== -1) continue;
        
        const pid = cells[i].plateId;
        const q = [i];
        componentId[i] = compCount;
        let size = 0;
        
        let head = 0;
        while(head < q.length) {
            const curr = q[head++];
            size++;
            for(const nId of cells[curr].neighbors) {
                if(componentId[nId] === -1 && cells[nId].plateId === pid) {
                    componentId[nId] = compCount;
                    q.push(nId);
                }
            }
        }
        compSize.push(size);
        compPlate.push(pid);
        compCount++;
    }

    const largestCompForPlate = new Int32Array(numPlates).fill(-1);
    const maxS = new Int32Array(numPlates).fill(-1);
    
    for(let c=0; c<compCount; c++) {
        const pid = compPlate[c];
        if(compSize[c] > maxS[pid]) {
            maxS[pid] = compSize[c];
            largestCompForPlate[pid] = c;
        }
    }

    const isOrphan = (cIdx: number) => {
        const pid = compPlate[cIdx];
        return largestCompForPlate[pid] !== cIdx;
    }

    const compCells: number[][] = Array.from({length: compCount}, () => []);
    for(let i=0; i<cells.length; i++) {
        compCells[componentId[i]].push(i);
    }

    const orphanIndices: number[] = [];
    for(let c=0; c<compCount; c++) {
        if(isOrphan(c)) orphanIndices.push(c);
    }
    orphanIndices.sort((a,b) => compSize[a] - compSize[b]);

    orphanIndices.forEach(cIdx => {
        const myCells = compCells[cIdx];
        const neighborCounts = new Map<number, number>();
        
        for(const cellId of myCells) {
            for(const nId of cells[cellId].neighbors) {
                const nComp = componentId[nId];
                if(nComp !== cIdx) {
                    const nPlate = cells[nId].plateId;
                    neighborCounts.set(nPlate, (neighborCounts.get(nPlate) || 0) + 1);
                }
            }
        }

        let bestP = -1;
        let maxCount = -1;
        neighborCounts.forEach((count, pid) => {
            if(count > maxCount) { maxCount = count; bestP = pid; }
        });

        if(bestP !== -1) {
            for(const cellId of myCells) {
                cells[cellId].plateId = bestP;
            }
        }
    });
}

// --- GEOGRAPHY GENERATION ---

export async function generateWorld(params: WorldParams, onLog?: (msg: string) => void, signal?: AbortSignal): Promise<WorldData> {
  onLog?.(`Initializing Grid (${params.points} cells)...`);
  const macroRng = new RNG(params.seed + '_macro');
  const simplex = new SimplexNoise(new RNG(params.seed));
  
  const points = generateFibonacciSphere(params.points, macroRng, params.cellJitter * 0.8);
  const geoPoints: [number, number][] = points.map(p => {
     const [lat, lon] = toSpherical(p.x, p.y, p.z);
     return [lon, lat]; 
  });
  
  onLog?.("Computing Connectivity...");
  await new Promise(r => setTimeout(r, 10)); 
  checkAbort(signal);

  const voronoi = geoVoronoi(geoPoints);
  const polygons = voronoi.polygons();
  const links = voronoi.links().features;

  const cells: Cell[] = points.map((p, i) => {
     const feature = polygons.features[i];
     let vertices: Point[] = [];
     if (feature && feature.geometry) {
        vertices = feature.geometry.coordinates[0].map((coord: any) => {
            const lon = (coord[0] * Math.PI) / 180;
            const lat = (coord[1] * Math.PI) / 180;
            return { x: Math.cos(lat) * Math.cos(lon), y: Math.sin(lat), z: Math.cos(lat) * Math.sin(lon) };
        });
        if (vertices.length > 0) vertices.pop();
     }
     return { id: i, center: p, vertices, neighbors: [], height: 0, plateId: 0, temperature: 0, moisture: 0, biome: BiomeType.OCEAN };
  });

  const coordIdMap = new Map<string, number>();
  const getKey = (coord: number[]) => `${coord[0].toFixed(4)},${coord[1].toFixed(4)}`;
  geoPoints.forEach((p, i) => coordIdMap.set(getKey(p), i));

  links.forEach((link: any) => {
     const p0 = link.geometry.coordinates[0];
     const p1 = link.geometry.coordinates[1];
     const i0 = coordIdMap.get(getKey(p0));
     const i1 = coordIdMap.get(getKey(p1));
     if (i0 !== undefined && i1 !== undefined && i0 !== i1) {
         cells[i0].neighbors.push(i1);
         cells[i1].neighbors.push(i0);
     }
  });
  cells.forEach(c => c.neighbors = [...new Set(c.neighbors)]);

  onLog?.(`Simulating ${params.plates} Tectonic Plates...`);
  await new Promise(r => setTimeout(r, 0));
  checkAbort(signal);

  const numPlates = params.plates;
  const plateRng = new RNG(params.seed + '_plates_loc'); 
  
  const plateVectors: Point[] = [];
  for(let i=0; i<numPlates; i++) {
      plateVectors.push(randomVector(plateRng));
  }

  const warpNoise = new SimplexNoise(new RNG(params.seed + '_warp'));
  const warpFreq = 0.5; 
  const warpAmp = (params.warpStrength ?? 0.5) * 0.2; 

  cells.forEach(cell => {
      const nx = warpNoise.noise3D(cell.center.x * warpFreq, cell.center.y * warpFreq, cell.center.z * warpFreq);
      const ny = warpNoise.noise3D(cell.center.y * warpFreq, cell.center.z * warpFreq, cell.center.x * warpFreq);
      const nz = warpNoise.noise3D(cell.center.z * warpFreq, cell.center.x * warpFreq, cell.center.y * warpFreq);
      
      const wx = cell.center.x + nx * warpAmp;
      const wy = cell.center.y + ny * warpAmp;
      const wz = cell.center.z + nz * warpAmp;
      
      let minDist = Infinity;
      let bestPlate = 0;
      
      for(let i=0; i<numPlates; i++) {
          const p = plateVectors[i];
          const d = (wx - p.x)**2 + (wy - p.y)**2 + (wz - p.z)**2;
          if (d < minDist) {
              minDist = d;
              bestPlate = i;
          }
      }
      cell.plateId = bestPlate;
  });

  enforceConnectivity(cells, numPlates);

  const moveRng = new RNG(params.seed + '_plates_move');
  const plateDrift = plateVectors.map(() => ({ 
      x: moveRng.next() - 0.5, 
      y: moveRng.next() - 0.5, 
      z: moveRng.next() - 0.5 
  }));

  const cellStress = new Float32Array(cells.length).fill(0); 
  const distToEdge = new Float32Array(cells.length).fill(0);

  cells.forEach(c => {
      let isBoundary = false;
      let maxStress = 0;
      
      for (const nId of c.neighbors) {
          const n = cells[nId];
          if (n.plateId !== c.plateId) {
              isBoundary = true;
              const driftA = plateDrift[c.plateId % plateDrift.length];
              const driftB = plateDrift[n.plateId % plateDrift.length];
              const dx = n.center.x - c.center.x;
              const dy = n.center.y - c.center.y;
              const dz = n.center.z - c.center.z;
              const rvx = driftA.x - driftB.x;
              const rvy = driftA.y - driftB.y;
              const rvz = driftA.z - driftB.z;
              const dot = (rvx*dx + rvy*dy + rvz*dz) * 10; 
              if (Math.abs(dot) > Math.abs(maxStress)) maxStress = dot;
          }
      }
      
      if (isBoundary) {
          cellStress[c.id] = maxStress;
          distToEdge[c.id] = 0; 
      } else {
          distToEdge[c.id] = 1.0; 
      }
  });

  const spreadIterations = Math.max(2, Math.floor(4 * Math.sqrt(params.points / 4000)));
  
  const nextStress = new Float32Array(cells.length);
  const nextDist = new Float32Array(cells.length);

  for(let i=0; i<spreadIterations; i++) {
      cells.forEach(c => {
          let stressSum = cellStress[c.id];
          let distSum = distToEdge[c.id];
          let count = 1;
          c.neighbors.forEach(nId => {
              stressSum += cellStress[nId];
              distSum += distToEdge[nId];
              count++;
          });
          nextStress[c.id] = stressSum / count;
          nextDist[c.id] = distSum / count + 0.1; 
      });
      nextStress.forEach((v,k) => cellStress[k] = v);
      nextDist.forEach((v,k) => distToEdge[k] = v);
  }

  onLog?.("Applying Height & Noise...");
  const freq = params.noiseScale || 1.0;
  const plateInf = (params.plateInfluence === undefined ? 0.5 : params.plateInfluence); 

  const plateHeights = new Float32Array(numPlates);
  const pRng = new RNG(params.seed + '_plates_h');
  
  let landChance = 0.45;
  let landLevel = 0.2;
  let oceanLevel = -0.5;

  if (params.landStyle === 'Archipelago') { landChance = 0.25; landLevel = 0.1; oceanLevel = -0.3; }
  if (params.landStyle === 'Islands') { landChance = 0.15; landLevel = 0.2; oceanLevel = -0.6; }

  for (let i = 0; i < numPlates; i++) {
      const isLand = pRng.next() < landChance;
      plateHeights[i] = isLand ? (landLevel + pRng.next() * 0.3) : (oceanLevel + pRng.next() * 0.3);
  }

  cells.forEach(c => {
      const fbmVal = fbm(simplex, c.center.x * freq, c.center.y * freq, c.center.z * freq, 3, 0.5, 2.0);
      const ridgedVal = ridgedNoise(simplex, c.center.x * freq, c.center.y * freq, c.center.z * freq, 3, 2.0);
      const ridgedRemapped = (ridgedVal * 2.0) - 1.0;
      const blend = params.ridgeBlend === undefined ? 0 : params.ridgeBlend;
      const structuralNoise = fbmVal * (1 - blend) + ridgedRemapped * blend;
      let baseSum = 0; 
      let bCount = 0;
      c.neighbors.forEach(n => { baseSum += plateHeights[cells[n].plateId]; bCount++; });
      baseSum += plateHeights[c.plateId]; bCount++;
      const avgBase = baseSum / bCount;
      const influence = Math.min(1, Math.max(0.1, plateInf));
      let height = avgBase * influence + structuralNoise * (1.2 - influence);
      const stress = cellStress[c.id]; 
      const edgeProx = Math.max(0, 1.0 - distToEdge[c.id] * 0.5); 
      if (edgeProx > 0) {
          if (stress > 0.05) {
              const mtnHeight = stress * edgeProx * 1.5;
              const ridge = ridgedNoise(simplex, c.center.x * freq, c.center.y * freq, c.center.z * freq, 4, 2.5);
              height += mtnHeight + (ridge * 0.3 * mtnHeight);
          } else if (stress < -0.05) {
              height -= Math.abs(stress) * edgeProx * 1.0;
          }
      }
      const detail = fbm(simplex, c.center.x * 6, c.center.y * 6, c.center.z * 6, 2, 0.5, 2.5);
      height += detail * params.roughness * 0.15;
      if (height > -0.2 && height < 0.2) {
          height = height * 0.5 + (height > 0 ? 0.05 : -0.05);
      }
      if (params.maskType === 'Pangea') {
          const mask = (c.center.x * 0.8 + c.center.y * 0.2 + 1) * 0.5;
          const smoothMask = mask * mask * (3 - 2 * mask);
          height = height * 0.5 + smoothMask * 0.8 - 0.2;
      }
      c.height = height;
  });
  
  let minH = Infinity, maxH = -Infinity;
  cells.forEach(c => { if (c.height < minH) minH = c.height; if (c.height > maxH) maxH = c.height; });
  let range = maxH - minH || 1;
  cells.forEach(c => c.height = (c.height - minH) / range);

  // EROSION
  if (params.erosionIterations > 0) {
      onLog?.(`Eroding Terrain (${params.erosionIterations} iter)...`);
      await new Promise(r => setTimeout(r, 0));
      checkAbort(signal);

      const resFactor = Math.sqrt(params.points / 5000);
      const hydraulicSteps = Math.ceil(params.erosionIterations * 2 * resFactor);
      const thermalSteps = Math.ceil(params.erosionIterations * 0.5 * resFactor);
      
      await applyHydraulicErosion(cells, hydraulicSteps, params.seaLevel, signal); 
      await applyThermalErosion(cells, thermalSteps, signal);
      
      minH = Infinity; maxH = -Infinity;
      cells.forEach(c => { if (c.height < minH) minH = c.height; if (c.height > maxH) maxH = c.height; });
      range = maxH - minH || 1;
      cells.forEach(c => c.height = (c.height - minH) / range);
  }

  onLog?.("Calculating Climate (Wind & Rain)...");
  
  const windVectors = cells.map(c => {
      const tiltRad = (params.axialTilt || 0) * (Math.PI / 180);
      const cosT = Math.cos(tiltRad);
      const sinT = Math.sin(tiltRad);
      const rotY = c.center.y * cosT - c.center.x * sinT; 
      const lat = Math.asin(Math.max(-1, Math.min(1, rotY))); 
      const latDeg = lat * (180 / Math.PI);
      let dir = 1; 
      if (Math.abs(latDeg) < 30) dir = -1; 
      else if (Math.abs(latDeg) < 60) dir = 1; 
      else dir = -1; 

      const len = Math.sqrt(c.center.x*c.center.x + c.center.z*c.center.z);
      if (len === 0) return {x:0, y:0, z:0};
      
      return { x: (-c.center.z / len) * dir, y: 0, z: (c.center.x / len) * dir };
  });

  cells.forEach(c => { 
      if (c.height < params.seaLevel) c.moisture = 1.0; 
      else c.moisture = 0.1 * params.rainfallMultiplier; 
  });
  
  const moistureMix = params.moistureTransport === undefined ? 0.5 : params.moistureTransport;
  
  for(let pass=0; pass<8; pass++) {
      const newMoisture = new Float32Array(cells.length);
      cells.forEach((c, i) => {
          if (c.height < params.seaLevel) { 
              newMoisture[i] = 1.0; 
              return; 
          }
          let incomingMoisture = 0; 
          let count = 0;
          c.neighbors.forEach(nId => {
             const n = cells[nId];
             const wind = windVectors[nId];
             const dx = c.center.x - n.center.x; 
             const dy = c.center.y - n.center.y;
             const dz = c.center.z - n.center.z;
             const dot = dx*wind.x + dy*wind.y + dz*wind.z; 
             
             if (dot > 0) { 
                 let carry = n.moisture;
                 const heightDiff = c.height - n.height;
                 if (heightDiff > 0.02) carry *= 1.5;
                 else if (heightDiff < -0.02) carry *= 0.2; 
                 incomingMoisture += carry; 
                 count++; 
             }
          });
          if (count === 0) { newMoisture[i] = c.moisture * 0.95; return; }
          incomingMoisture /= count; 
          newMoisture[i] = c.moisture * (1 - moistureMix) + incomingMoisture * moistureMix;
          if (c.height > params.seaLevel) newMoisture[i] *= 0.98; 
      });
      cells.forEach((c, i) => c.moisture = newMoisture[i]);
  }
  
  const tempVariance = params.temperatureVariance === undefined ? 5 : params.temperatureVariance;
  cells.forEach(c => {
      const tiltRad = (params.axialTilt || 0) * (Math.PI / 180);
      const rotY = c.center.y * Math.cos(tiltRad) - c.center.x * Math.sin(tiltRad);
      const lat = Math.asin(Math.max(-1, Math.min(1, rotY)));
      const latRatio = Math.abs(lat) / (Math.PI / 2); 
      let temp = params.baseTemperature * (1 - latRatio * latRatio) + params.poleTemperature * (latRatio * latRatio);
      const elevation = Math.max(0, c.height - params.seaLevel);
      temp -= elevation * 60;
      if (tempVariance > 0) temp += simplex.noise3D(c.center.x * 5, c.center.y * 5, c.center.z * 5) * tempVariance;
      c.temperature = temp;
      c.moisture = Math.max(0, Math.min(1, c.moisture * params.rainfallMultiplier));
      c.biome = determineBiome(c.height, c.temperature, c.moisture, params.seaLevel);
  });

  const rivers = await generateRivers(cells, params.seaLevel, params, onLog, signal);
  const world: WorldData = { cells, params, geoJson: polygons, rivers };
  
  return recalculateCivs(world, params, onLog);
}

// ... (recalculateCivs and recalculateProvinces remain unchanged as they are synchronous)
export function recalculateCivs(world: WorldData, params: WorldParams, onLog?: (msg: string) => void): WorldData {
    onLog?.(`Forging ${params.numFactions} Civilizations...`);
    
    // Reset
    world.cells.forEach(c => {
        c.regionId = undefined;
        c.provinceId = undefined;
        c.isCapital = false;
        c.isTown = false;
        c.population = 0;
    });

    const civRng = new RNG(params.civSeed);
    const numFactions = params.numFactions;
    const factions: FactionData[] = [];
    
    const suitable = world.cells.filter(c => 
        c.height >= params.seaLevel && 
        c.biome !== BiomeType.ICE_CAP &&
        c.biome !== BiomeType.VOLCANIC
    );
    
    const candidates = suitable.length > 0 ? suitable : world.cells.filter(c => c.height >= params.seaLevel);
    
    if (candidates.length === 0) {
        world.civData = { factions: [] };
        return world;
    }

    const capitals: number[] = [];
    const minDist = (world.cells.length / numFactions) * params.capitalSpacing * 0.5; 
    
    let attempts = 0;
    while(capitals.length < numFactions && attempts < 1000) {
        attempts++;
        const candidate = candidates[Math.floor(civRng.next() * candidates.length)];
        
        let tooClose = false;
        for(const capId of capitals) {
            const cap = world.cells[capId];
            const d = (candidate.center.x - cap.center.x)**2 + (candidate.center.y - cap.center.y)**2 + (candidate.center.z - cap.center.z)**2;
            if (d < minDist * 0.0001) { 
                tooClose = true;
                break;
            }
        }
        
        if (!tooClose) {
            capitals.push(candidate.id);
            candidate.isCapital = true;
            candidate.regionId = capitals.length - 1;
            factions.push({
                id: capitals.length - 1,
                name: `Faction ${capitals.length}`,
                color: '#ffffff', 
                capitalId: candidate.id,
                provinces: [],
                totalPopulation: 0
            });
        }
    }

    const pq = new MinHeap<{id: number, cost: number, region: number}>(x => x.cost);
    const costs = new Map<number, number>();
    
    capitals.forEach((capId, idx) => {
        pq.push({ id: capId, cost: 0, region: idx });
        costs.set(capId, 0);
    });

    const waterCost = (params.waterCrossingCost || 0.5) * 50; 
    const landCost = 1;
    const territorialRange = (params.territorialWaters || 0.2) * 50; 

    while(pq.size() > 0) {
        const { id, cost, region } = pq.pop()!;
        if (world.cells[id].regionId !== undefined && world.cells[id].regionId !== region) continue;
        world.cells[id].regionId = region;
        if (cost > 200) continue; 
        const currCell = world.cells[id];
        for(const nId of currCell.neighbors) {
            const nCell = world.cells[nId];
            let moveCost = landCost;
            if (nCell.biome === BiomeType.ICE_CAP) moveCost *= 4;
            if (nCell.biome === BiomeType.HOT_DESERT) moveCost *= 2;
            if (nCell.biome === BiomeType.VOLCANIC) moveCost *= 5;
            const slope = Math.abs(nCell.height - currCell.height);
            moveCost += slope * 20;
            const isWater = nCell.height < params.seaLevel;
            if (isWater) moveCost = waterCost;
            moveCost *= (1 + (civRng.next() * params.borderRoughness)); 
            const newCost = cost + moveCost;
            if (isWater && newCost > territorialRange) continue;
            if (!costs.has(nId) || newCost < costs.get(nId)!) {
                costs.set(nId, newCost);
                if (world.cells[nId].regionId === undefined) {
                     pq.push({ id: nId, cost: newCost, region });
                }
            }
        }
    }
    
    const palette = [
        '#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'
    ];
    for (let i = palette.length - 1; i > 0; i--) {
        const j = Math.floor(civRng.next() * (i + 1));
        [palette[i], palette[j]] = [palette[j], palette[i]];
    }
    
    factions.forEach((f, i) => f.color = palette[i % palette.length]);
    
    world.civData = { factions };
    
    return recalculateProvinces(world, params);
}

export function recalculateProvinces(world: WorldData, params: WorldParams): WorldData {
    if (!world.civData) return world;
    const provRng = new RNG(params.civSeed + '_prov');

    world.cells.forEach(c => {
        let suitability = 0;
        if (c.height < params.seaLevel) {
            c.population = 0;
            return;
        }
        switch(c.biome) {
            case BiomeType.TROPICAL_RAINFOREST: suitability = 0.4; break;
            case BiomeType.TROPICAL_SAVANNA: suitability = 0.7; break;
            case BiomeType.HOT_DESERT: suitability = 0.1; break;
            case BiomeType.COLD_DESERT: suitability = 0.1; break;
            case BiomeType.TEMPERATE_FOREST: suitability = 0.9; break;
            case BiomeType.TEMPERATE_RAINFOREST: suitability = 0.8; break;
            case BiomeType.MEDITERRANEAN: suitability = 1.0; break; 
            case BiomeType.STEPPE: suitability = 0.5; break;
            case BiomeType.BOREAL_FOREST: suitability = 0.4; break;
            case BiomeType.TUNDRA: suitability = 0.2; break;
            case BiomeType.ICE_CAP: suitability = 0.0; break;
            case BiomeType.VOLCANIC: suitability = 0.1; break;
            case BiomeType.BEACH: suitability = 0.6; break;
            default: suitability = 0.5;
        }
        if ((c.flux || 0) > 0.5) suitability += 0.3;
        if ((c.flux || 0) > 2.0) suitability += 0.2; 
        let coast = false;
        for(const nId of c.neighbors) {
            if (world.cells[nId].height < params.seaLevel) { coast = true; break; }
        }
        if (coast) suitability += 0.3;
        if (c.height > 0.6) suitability -= (c.height - 0.6) * 2;
        c.population = Math.floor(suitability * 10000 * (0.8 + provRng.next() * 0.4));
    });

    world.civData.factions.forEach(faction => {
        faction.provinces = [];
        faction.totalPopulation = 0;
        const landCells = world.cells.filter(c => c.regionId === faction.id && c.height >= params.seaLevel);
        if (landCells.length === 0) return;
        const density = params.provinceSize || 0.5; 
        const targetSize = 20 + density * 100; 
        let numProvinces = Math.max(1, Math.ceil(landCells.length / targetSize));
        
        const townIds: number[] = [faction.capitalId];
        const capitalCell = world.cells[faction.capitalId];
        capitalCell.provinceId = 0; 
        
        let attempts = 0;
        while(townIds.length < numProvinces && attempts < 500) {
            attempts++;
            const candidate = landCells[Math.floor(provRng.next() * landCells.length)];
            if (townIds.includes(candidate.id)) continue;
            let tooClose = false;
            for(const tId of townIds) {
                const t = world.cells[tId];
                const d = (candidate.center.x - t.center.x)**2 + (candidate.center.y - t.center.y)**2 + (candidate.center.z - t.center.z)**2;
                if (d < 0.005 * density) { 
                    tooClose = true;
                    break;
                }
            }
            if (!tooClose) townIds.push(candidate.id);
        }
        
        townIds.forEach((tId, idx) => {
            const tCell = world.cells[tId];
            tCell.isTown = true;
            tCell.population = (tCell.population || 0) * 5; 
            faction.provinces.push({
                id: idx,
                name: idx === 0 ? "Capital Region" : `Province ${idx}`,
                towns: [{ name: idx === 0 ? "Capital City" : "Town", cellId: tId, population: tCell.population || 0, isCapital: tId === faction.capitalId }],
                totalPopulation: 0
            });
        });

        const pq = new MinHeap<{id: number, cost: number, provIdx: number}>(x => x.cost);
        const costs = new Map<number, number>();
        townIds.forEach((tId, idx) => {
            pq.push({ id: tId, cost: 0, provIdx: idx });
            costs.set(tId, 0);
        });

        const claimed = new Set<number>();
        while(pq.size() > 0) {
            const { id, cost, provIdx } = pq.pop()!;
            if (world.cells[id].regionId !== faction.id) continue;
            if (world.cells[id].provinceId !== undefined && world.cells[id].provinceId !== provIdx) continue;
            world.cells[id].provinceId = provIdx;
            claimed.add(id);
            faction.provinces[provIdx].totalPopulation += (world.cells[id].population || 0);
            for(const nId of world.cells[id].neighbors) {
                if (world.cells[nId].regionId !== faction.id) continue;
                if (world.cells[nId].provinceId !== undefined) continue;
                let moveCost = 1;
                moveCost += Math.abs(world.cells[nId].height - world.cells[id].height) * 10;
                const newCost = cost + moveCost;
                if (!costs.has(nId) || newCost < costs.get(nId)!) {
                    costs.set(nId, newCost);
                    pq.push({ id: nId, cost: newCost, provIdx });
                }
            }
        }
        faction.totalPopulation = faction.provinces.reduce((sum, p) => sum + p.totalPopulation, 0);
    });

    return world;
}