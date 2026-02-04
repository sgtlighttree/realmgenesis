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

function applyHydraulicErosion(cells: Cell[], iterations: number) {
    cells.forEach(c => c.flux = 0);
    const sorted = [...cells].sort((a, b) => b.height - a.height);
    // Tuned constants for decent looking rivers
    const erosionRate = 0.02;
    const depositionRate = 0.01;
    const rainAmount = 0.1;

    for (let iter = 0; iter < iterations; iter++) {
        sorted.forEach(c => c.flux = rainAmount);
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
                const safeErosion = Math.min(erosion, slope * 0.9); // Don't dig pits
                c.height -= safeErosion;
                target.height += safeErosion * depositionRate; 
            }
        });
    }
}

function applyThermalErosion(cells: Cell[], iterations: number) {
    const talus = 0.008; // Min slope diff
    const rate = 0.2; 
    for(let iter=0; iter<iterations; iter++) {
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

// --- BIOME ---

function determineBiome(height: number, temp: number, moisture: number, seaLevel: number): BiomeType {
  if (height < seaLevel) {
    if (height < seaLevel * 0.6) return BiomeType.DEEP_OCEAN;
    return BiomeType.OCEAN;
  }
  const landH = (height - seaLevel) / (1 - seaLevel);
  if (landH < 0.02 && temp > 5) return BiomeType.BEACH;
  if (landH > 0.85 && temp > -5) return BiomeType.VOLCANIC;
  if (temp < -5) return BiomeType.ICE_CAP;
  if (temp < 5) return BiomeType.TUNDRA;
  
  const aridityThreshold = (temp + 10) / 100; 
  if (moisture < aridityThreshold) {
      if (moisture < aridityThreshold * 0.5) return temp > 18 ? BiomeType.HOT_DESERT : BiomeType.COLD_DESERT;
      else return BiomeType.STEPPE;
  }
  if (temp > 18) {
      if (moisture > 0.6) return BiomeType.TROPICAL_RAINFOREST;
      return BiomeType.TROPICAL_SAVANNA;
  }
  if (temp < 12) return BiomeType.BOREAL_FOREST;
  if (moisture < 0.5) return BiomeType.MEDITERRANEAN;
  if (moisture > 0.75) return BiomeType.TEMPERATE_RAINFOREST;
  return BiomeType.TEMPERATE_FOREST;
}

// --- TECTONIC HELPERS ---

function enforceConnectivity(cells: Cell[], numPlates: number) {
    // Robust Region Merging to remove enclaves
    // 1. Label connected components via BFS
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

    // 2. Identify "Main Body" (largest component) for each Plate ID
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
        // If a plate has no components (rare), this check is safe
        return largestCompForPlate[pid] !== cIdx;
    }

    // 3. Collect orphans and adopt them to dominant neighbors
    // We group cells by component for easy reassignment
    const compCells: number[][] = Array.from({length: compCount}, () => []);
    for(let i=0; i<cells.length; i++) {
        compCells[componentId[i]].push(i);
    }

    // Sort orphans by size (smallest first) to dissolve specks into larger regions
    const orphanIndices: number[] = [];
    for(let c=0; c<compCount; c++) {
        if(isOrphan(c)) orphanIndices.push(c);
    }
    orphanIndices.sort((a,b) => compSize[a] - compSize[b]);

    orphanIndices.forEach(cIdx => {
        const myCells = compCells[cIdx];
        const neighborCounts = new Map<number, number>();
        
        // Scan border to find dominant neighbor plate
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
            // Reassign
            for(const cellId of myCells) {
                cells[cellId].plateId = bestP;
                // Note: We don't update componentId/compSize on the fly, 
                // assuming adoption by a massive plate makes intermediate state irrelevant.
            }
        }
    });
}

// --- GEOGRAPHY GENERATION ---

export async function generateWorld(params: WorldParams, onProgress?: (msg: string, pct: number) => void): Promise<WorldData> {
  onProgress?.("Initializing Grid...", 10);
  const macroRng = new RNG(params.seed + '_macro');
  const simplex = new SimplexNoise(new RNG(params.seed));
  
  const points = generateFibonacciSphere(params.points, macroRng, params.cellJitter * 0.8);
  const geoPoints: [number, number][] = points.map(p => {
     const [lat, lon] = toSpherical(p.x, p.y, p.z);
     return [lon, lat]; 
  });
  
  onProgress?.("Computing Connectivity...", 20);
  await new Promise(r => setTimeout(r, 0));
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

  onProgress?.("Simulating Tectonics...", 40);
  const numPlates = params.plates;
  const plateRng = new RNG(params.seed + '_plates_loc'); // Dedicated seed for plate locations
  
  // 1. Plate Centers (Resolution Independent)
  // We place vectors on the sphere. These are the "soul" of the plates.
  const plateVectors: Point[] = [];
  for(let i=0; i<numPlates; i++) {
      plateVectors.push(randomVector(plateRng));
  }

  // 2. Assign Plates (Warped Voronoi) - O(N * Plates)
  // Much faster than Dijkstra for high N, and Resolution Independent
  const warpNoise = new SimplexNoise(new RNG(params.seed + '_warp'));
  
  // Lower frequency = smoother larger curves
  const warpFreq = 0.5; 
  // Lower amplitude scaling = prevents extreme folds/enclaves
  const warpAmp = (params.warpStrength ?? 0.5) * 0.2; 

  cells.forEach(cell => {
      // Warped coordinates
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

  // 3. Clean up artifacts (Merge disjoint enclaves)
  enforceConnectivity(cells, numPlates);

  // 4. Plate Movement & Stress
  // These determine mountains/rifts
  const moveRng = new RNG(params.seed + '_plates_move');
  const plateDrift = plateVectors.map(() => ({ 
      x: moveRng.next() - 0.5, 
      y: moveRng.next() - 0.5, 
      z: moveRng.next() - 0.5 
  }));

  const cellStress = new Float32Array(cells.length).fill(0); 
  const distToEdge = new Float32Array(cells.length).fill(0);

  // Calculate stress and "Edge proximity"
  // We do one pass to find boundary cells
  cells.forEach(c => {
      let isBoundary = false;
      let maxStress = 0;
      
      for (const nId of c.neighbors) {
          const n = cells[nId];
          if (n.plateId !== c.plateId) {
              isBoundary = true;
              
              // Stress calculation
              const driftA = plateDrift[c.plateId % plateDrift.length];
              const driftB = plateDrift[n.plateId % plateDrift.length];
              const dx = n.center.x - c.center.x;
              const dy = n.center.y - c.center.y;
              const dz = n.center.z - c.center.z;
              const rvx = driftA.x - driftB.x;
              const rvy = driftA.y - driftB.y;
              const rvz = driftA.z - driftB.z;
              // Converging = Positive, Diverging = Negative
              const dot = (rvx*dx + rvy*dy + rvz*dz) * 10; 
              if (Math.abs(dot) > Math.abs(maxStress)) maxStress = dot;
          }
      }
      
      if (isBoundary) {
          cellStress[c.id] = maxStress;
          distToEdge[c.id] = 0; // It is on the edge
      } else {
          distToEdge[c.id] = 1.0; // Interior
      }
  });

  // Diffuse Edge Distance and Stress inwards
  // Scaling iterations by resolution so mountains have consistent width
  // Base 4000 points = ~2 iters. 40000 points = ~6 iters
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
          nextDist[c.id] = distSum / count + 0.1; // Add distance penalty
      });
      nextStress.forEach((v,k) => cellStress[k] = v);
      nextDist.forEach((v,k) => distToEdge[k] = v);
  }

  onProgress?.("Generating Terrain...", 60);
  const featureFreq = params.noiseScale || 1.0;
  const plateInf = (params.plateInfluence === undefined ? 0.5 : params.plateInfluence); 

  // --- Assign Base Height based on Plate ID ---
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
      // 1. Structural Noise (Fractal)
      const structuralNoise = fbm(simplex, c.center.x, c.center.y, c.center.z, 3, 0.5, 1.5 * featureFreq);
      
      // 2. Base Plate Height 
      // We perform a local average to smooth the sharp transition between plates
      let baseSum = 0; 
      let bCount = 0;
      c.neighbors.forEach(n => { baseSum += plateHeights[cells[n].plateId]; bCount++; });
      baseSum += plateHeights[c.plateId]; bCount++;
      const avgBase = baseSum / bCount;

      const influence = Math.min(1, Math.max(0.1, plateInf));
      
      // Mix Plate Base and Noise
      let height = avgBase * influence + structuralNoise * (1.2 - influence);
      
      // 3. Tectonic Interaction (Mountains/Rifts)
      const stress = cellStress[c.id]; // -1 to 1 typically
      const edgeProx = Math.max(0, 1.0 - distToEdge[c.id] * 0.5); // 1.0 at edge, 0.0 deep inside
      
      if (edgeProx > 0) {
          if (stress > 0.05) {
              // Converging: Mountains
              // We multiply by edgeProx so mountains only appear near borders
              const mtnHeight = stress * edgeProx * 1.5;
              const ridge = ridgedNoise(simplex, c.center.x, c.center.y, c.center.z, 4, 2.5);
              height += mtnHeight + (ridge * 0.3 * mtnHeight);
          } else if (stress < -0.05) {
              // Diverging: Rifts
              height -= Math.abs(stress) * edgeProx * 1.0;
          }
      }

      // 4. Detail & Roughness
      const detail = fbm(simplex, c.center.x * 6, c.center.y * 6, c.center.z * 6, 2, 0.5, 2.5);
      height += detail * params.roughness * 0.15;

      // 5. Continental Shelf Curve
      // Flatten values near sea level to create shelves
      // Soft clamp
      if (height > -0.2 && height < 0.2) {
          height = height * 0.5 + (height > 0 ? 0.05 : -0.05);
      }

      // Pangea Mask Logic
      if (params.maskType === 'Pangea') {
          const mask = (c.center.x * 0.8 + c.center.y * 0.2 + 1) * 0.5;
          const smoothMask = mask * mask * (3 - 2 * mask);
          height = height * 0.5 + smoothMask * 0.8 - 0.2;
      }
      
      c.height = height;
  });
  
  // Normalize Height
  let minH = Infinity, maxH = -Infinity;
  cells.forEach(c => { if (c.height < minH) minH = c.height; if (c.height > maxH) maxH = c.height; });
  let range = maxH - minH || 1;
  cells.forEach(c => c.height = (c.height - minH) / range);

  if (params.erosionIterations > 0) {
      onProgress?.("Eroding Terrain...", 70);
      await new Promise(r => setTimeout(r, 0));
      // Scale erosion steps by resolution 
      // At higher res, water travels less distance per step, so we need more steps
      const resFactor = Math.sqrt(params.points / 5000);
      const hydraulicSteps = Math.ceil(params.erosionIterations * 2 * resFactor);
      const thermalSteps = Math.ceil(params.erosionIterations * 0.5 * resFactor);
      
      applyHydraulicErosion(cells, hydraulicSteps);
      applyThermalErosion(cells, thermalSteps);

      // Renormalize after erosion
      minH = Infinity; maxH = -Infinity;
      cells.forEach(c => { if (c.height < minH) minH = c.height; if (c.height > maxH) maxH = c.height; });
      range = maxH - minH || 1;
      cells.forEach(c => c.height = (c.height - minH) / range);
  }

  onProgress?.("Calculating Climate...", 80);
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

  cells.forEach(c => { if (c.height < params.seaLevel) c.moisture = 1.0; else c.moisture = 0.1 * params.rainfallMultiplier; });
  const moistureMix = params.moistureTransport === undefined ? 0.5 : params.moistureTransport;
  for(let pass=0; pass<6; pass++) {
      const newMoisture = new Float32Array(cells.length);
      cells.forEach((c, i) => {
          if (c.height < params.seaLevel) { newMoisture[i] = 1.0; return; }
          let incoming = 0; let count = 0;
          c.neighbors.forEach(nId => {
             const n = cells[nId]; const dx = c.center.x - n.center.x; const dz = c.center.z - n.center.z;
             const wind = windVectors[nId]; const dot = dx*wind.x + 0 + dz*wind.z; 
             if (dot > 0) { let carry = n.moisture; if (c.height > n.height + 0.05) carry *= 0.5; incoming += carry; count++; }
          });
          if (count === 0) { newMoisture[i] = c.moisture * 0.95; return; }
          incoming /= count; newMoisture[i] = c.moisture * (1 - moistureMix) + incoming * moistureMix;
          if (c.height > params.seaLevel + 0.2) newMoisture[i] *= 0.8; 
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
      temp -= Math.max(0, c.height - params.seaLevel) * 50; 
      if (tempVariance > 0) temp += simplex.noise3D(c.center.x * 5, c.center.y * 5, c.center.z * 5) * tempVariance;
      c.temperature = temp;
      c.moisture = Math.max(0, Math.min(1, c.moisture * params.rainfallMultiplier));
      c.biome = determineBiome(c.height, c.temperature, c.moisture, params.seaLevel);
  });

  const world: WorldData = { cells, params, geoJson: polygons };
  return recalculateCivs(world, params, onProgress);
}

export function recalculateCivs(world: WorldData, params: WorldParams, onProgress?: (msg: string, pct: number) => void): WorldData {
    onProgress?.("Placing Civilizations...", 90);
    
    world.cells.forEach(c => {
        c.regionId = undefined;
        c.provinceId = undefined;
        c.isCapital = false;
        c.isTown = false;
        c.population = 0;
    });

    const rng = new RNG(params.civSeed);
    const numFactions = params.numFactions;
    const capitals: number[] = [];
    
    const landCells = world.cells.filter(c => c.height >= params.seaLevel && c.biome !== BiomeType.ICE_CAP);
    if (landCells.length === 0) return world;

    // Poisson Disk-ish for Capitals
    for(let i=0; i<numFactions; i++) {
        let bestCell = -1;
        let maxDist = -1;
        for(let attempt=0; attempt<20; attempt++) {
            const candidate = landCells[Math.floor(rng.next() * landCells.length)];
            let minDist = Infinity;
            for(const capId of capitals) {
                const cap = world.cells[capId];
                // Euclidean dist on sphere is fine for relative comparisons
                const d = (candidate.center.x - cap.center.x)**2 + (candidate.center.y - cap.center.y)**2 + (candidate.center.z - cap.center.z)**2;
                if (d < minDist) minDist = d;
            }
            if (capitals.length === 0) minDist = Infinity;
            if (minDist > maxDist) { maxDist = minDist; bestCell = candidate.id; }
        }
        if (bestCell !== -1) {
            capitals.push(bestCell);
            world.cells[bestCell].isCapital = true;
            world.cells[bestCell].regionId = i;
        }
    }

    const factionColors = ['#ef4444', '#f97316', '#f59e0b', '#84cc16', '#10b981', '#06b6d4', '#3b82f6', '#6366f1', '#8b5cf6', '#d946ef', '#f43f5e', '#e11d48'];
    const factions: FactionData[] = capitals.map((capId, i) => ({
        id: i,
        name: `Faction ${i+1}`,
        color: factionColors[i % factionColors.length],
        capitalId: capId,
        provinces: [],
        totalPopulation: 0
    }));

    // Dijkstra Expansion using MinHeap
    const costs = new Float32Array(world.cells.length).fill(Infinity);
    const owner = new Int32Array(world.cells.length).fill(-1);
    const frontier = new MinHeap<{id: number, cost: number, owner: number, waterDist: number}>(x => x.cost);

    capitals.forEach(id => {
        costs[id] = 0;
        owner[id] = world.cells[id].regionId!;
        frontier.push({ id, cost: 0, owner: world.cells[id].regionId!, waterDist: 0 });
    });
    
    while(frontier.size() > 0) {
        const current = frontier.pop()!;
        if (current.cost > costs[current.id]) continue;
        
        const cell = world.cells[current.id];
        
        for(const nId of cell.neighbors) {
            const neighbor = world.cells[nId];
            let moveCost = 1;
            
            const isNeighborWater = neighbor.height < params.seaLevel;
            let newWaterDist = 0;

            const dist = Math.sqrt((cell.center.x - neighbor.center.x)**2 + (cell.center.y - neighbor.center.y)**2 + (cell.center.z - neighbor.center.z)**2);

            if (isNeighborWater) {
                newWaterDist = current.waterDist + dist;
                const limit = params.territorialWaters ?? 0.2;
                if (newWaterDist > limit) continue;
                if (params.waterCrossingCost > 0.8) moveCost = 1000; 
                else moveCost = dist * (5 + (params.waterCrossingCost * 20)); // Scale cost by distance (resolution independent)
            } else {
                 newWaterDist = 0;
                 const hDiff = Math.abs(neighbor.height - cell.height);
                 moveCost = dist * (1 + hDiff * 50); // Scale cost by distance
            }
            
            const newCost = current.cost + moveCost;
            if (newCost < costs[nId]) {
                costs[nId] = newCost;
                owner[nId] = current.owner;
                // Optimization: Don't push if cost is excessively high
                if (newCost < 500) { 
                    frontier.push({ id: nId, cost: newCost, owner: current.owner, waterDist: newWaterDist });
                }
            }
        }
    }
    
    world.cells.forEach((c, i) => { if (owner[i] !== -1) c.regionId = owner[i]; });
    world.civData = { factions };
    return recalculateProvinces(world, params);
}

export function recalculateProvinces(world: WorldData, params: WorldParams): WorldData {
    if (!world.civData) return world;
    const rng = new RNG(params.civSeed + '_provs');
    world.cells.forEach(c => { c.provinceId = undefined; c.isTown = false; });

    world.civData.factions.forEach(faction => {
        const factionCells = world.cells.filter(c => c.regionId === faction.id && c.height >= params.seaLevel);
        if (factionCells.length === 0) return;

        faction.provinces = [];
        const towns: number[] = [faction.capitalId]; 
        const numProvinces = Math.max(1, Math.floor(factionCells.length * (params.provinceSize * 0.05))); 
        
        for(let i=1; i<numProvinces; i++) {
             const candidate = factionCells[Math.floor(rng.next() * factionCells.length)];
             towns.push(candidate.id);
             world.cells[candidate.id].isTown = true;
        }

        factionCells.forEach(c => {
            let nearestTown = -1;
            let minDist = Infinity;
            towns.forEach((tId, idx) => {
                const t = world.cells[tId];
                const d = (c.center.x - t.center.x)**2 + (c.center.y - t.center.y)**2 + (c.center.z - t.center.z)**2;
                if (d < minDist) { minDist = d; nearestTown = idx; }
            });
            c.provinceId = nearestTown;
        });
        
        towns.forEach((tId, idx) => {
            const isCap = tId === faction.capitalId;
            const pop = Math.floor(rng.next() * 10000 + (isCap ? 50000 : 5000));
            faction.provinces.push({
                id: idx,
                name: isCap ? `Capital Province` : `Province ${idx}`,
                towns: [{ name: isCap ? `Capital City` : `Town ${idx}`, cellId: tId, population: pop, isCapital: isCap }],
                totalPopulation: pop
            });
        });
        faction.totalPopulation = faction.provinces.reduce((s, p) => s + p.totalPopulation, 0);
    });
    return world;
}