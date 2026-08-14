import { geoVoronoi } from 'd3-geo-voronoi';
import { Cell, LakeData, Point, WorldData, WorldParams, BiomeType, FactionData, CultureData, ReligionData } from '../types';
import { RNG, SimplexNoise } from './rng';
import { FACTION_COLORS, CULTURE_COLORS, RELIGION_COLORS, FOLK_COLORS } from './palette';
import { createNameGenerator, NameGenerator, NameStyle, NAME_STYLES } from './namegen';
import { detectFeatures } from './features';
import { MinHeap, landTerrainStepCost } from './pathfinding';
import { simulateTectonics, projectTectonicsToDisplay } from './tectonicsV3';

// --- DATA STRUCTURES ---

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

// --- NOISE ALGORITHMS ---

export function fbm(simplex: SimplexNoise, x: number, y: number, z: number, octaves: number, persistence: number, lacunarity: number): number {
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

export function ridgedNoise(simplex: SimplexNoise, x: number, y: number, z: number, octaves: number, lacunarity: number): number {
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

async function applyHydraulicErosion(cells: Cell[], iterations: number, seaLevel: number): Promise<void> {
    cells.forEach(c => c.flux = 0);
    const sorted = [...cells].sort((a, b) => b.height - a.height);
    const erosionRate = 0.02;
    const depositionRate = 0.01;
    const rainAmount = 0.1;

    // Each iteration rains on land, then routes flux and erodes/deposits
    // downhill across the whole cell set. Runs inside a worker, so no
    // per-iteration yielding is needed to keep the UI responsive.
    for (let iter = 0; iter < iterations; iter++) {
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

async function applyThermalErosion(cells: Cell[], iterations: number) {
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

// Lakes are land-height depressions, not ocean. Civ suitability, drainage, and
// coast checks must treat them as water even though height >= seaLevel.
export function isLakeCell(cell: Cell): boolean {
  return cell.biome === BiomeType.LAKE || cell.biome === BiomeType.SALT_LAKE;
}

// --- RIVER GENERATION ---

async function generateRivers(cells: Cell[], seaLevel: number, params: WorldParams, onProgress?: (msg: string) => void): Promise<{ rivers: Point[][]; lakes: LakeData[] }> {
    const numCells = cells.length;
    onProgress?.("Rivers: Initializing drainage map...");
    
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
        return { rivers: [], lakes: [] };
    }

    let processed = 0;
    
    onProgress?.(`Rivers: Propagating water levels...`);

    while(heap.size() > 0) {
        // Safety break and log update
        if (++processed % 2000 === 0) {
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

    // 1b. Surface depressions as lakes.
    // A land cell whose enforced water level sits above its own terrain is
    // flooded — the fill pooled water there. Group contiguous flooded land
    // cells into discrete lakes, then classify and stamp their biome so it
    // wins over the earlier terrestrial assignment. EPS guards float equality
    // (a freely-draining cell settles at exactly its own height).
    const EPS = 1e-9;
    const lakes: LakeData[] = [];
    const isFlooded = (id: number) => cells[id].height >= seaLevel && waterLevel[id] > cells[id].height + EPS;
    const lakeVisited = new Uint8Array(numCells);

    for (let s = 0; s < numCells; s++) {
        if (lakeVisited[s] || !isFlooded(s)) continue;
        // BFS the flooded component on the cell adjacency graph.
        const group: number[] = [];
        const stack = [s];
        lakeVisited[s] = 1;
        while (stack.length > 0) {
            const id = stack.pop()!;
            group.push(id);
            for (const nId of cells[id].neighbors) {
                if (!lakeVisited[nId] && isFlooded(nId)) {
                    lakeVisited[nId] = 1;
                    stack.push(nId);
                }
            }
        }

        const groupSet = new Set(group);
        let surfaceLevel = 0;
        let sumTemp = 0;
        let sumMoist = 0;
        // Outflow: the basin spills through one sill cell whose downstream
        // pointer leaves the lake toward the ocean. That neighbour is the pour
        // point just outside the lake (null only if the basin never drains).
        let outflowCellId: number | null = null;
        for (const id of group) {
            surfaceLevel = Math.max(surfaceLevel, waterLevel[id]);
            sumTemp += cells[id].temperature;
            sumMoist += cells[id].moisture;
            const d = downstream[id];
            if (outflowCellId === null && d !== -1 && !groupSet.has(d)) outflowCellId = d;
        }

        const meanTemp = sumTemp / group.length;
        const meanMoist = sumMoist / group.length;
        // Salt lakes form where a closed/arid basin loses water to evaporation
        // faster than inflow, concentrating dissolved salts. Model that as a
        // hot + dry basin: mean temp above ~18C and mean moisture below ~0.3.
        const isSalt = meanTemp > 18 && meanMoist < 0.3;
        // Priority-Flood gives every filled basin a spill point, so a truly
        // outflow-less basin is rare; a salt (evaporative) basin is a terminal
        // sink regardless, so it also reads as endorheic.
        const isEndorheic = outflowCellId === null || isSalt;

        const biome = isSalt ? BiomeType.SALT_LAKE : BiomeType.LAKE;
        for (const id of group) cells[id].biome = biome;

        lakes.push({ id: lakes.length, cellIds: group, surfaceLevel, outflowCellId, isEndorheic, isSalt });
    }

    onProgress?.("Rivers: Accumulating flux...");

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
    
    // Rivers must not start inside a lake — flux routes through the basin to
    // the outflow, where a fresh segment restarts on its own downstream path.
    const candidates = sortedIndices.filter(i => flux[i] > threshold && cells[i].height >= seaLevel && !isLakeCell(cells[i]));

    const getRenderPoint = (c: Cell) => {
        const r = 1 + c.height * 0.05 + 0.005;
        return { x: c.center.x * r, y: c.center.y * r, z: c.center.z * r };
    };

    for (const startId of candidates) {
        if (visited.has(startId)) continue;

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
            // Cut the polyline at the lake shore; the interior is open water.
            if (isLakeCell(cells[next])) {
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
    
    onProgress?.(`Rivers: Generated ${riverPaths.length} segments, ${lakes.length} lakes.`);
    return { rivers: riverPaths, lakes };
}

// --- BIOME ---

export function determineBiome(height: number, temp: number, moisture: number, seaLevel: number): BiomeType {
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

// --- GEOGRAPHY GENERATION ---

export async function generateWorld(params: WorldParams, onLog?: (msg: string) => void, signal?: AbortSignal, onProgress?: (stage: number, total: number) => void): Promise<WorldData> {
  // The only abort checkpoint left. The nine setTimeout(0) yields this used to
  // ride on existed to keep the MAIN thread barely responsive between stages;
  // generation now runs in a worker (utils/worldGenClient.ts), where they cost
  // time and buy nothing. Mid-run cancellation is worker.terminate() — a
  // synchronous generation loop cannot drain the message queue, so a
  // message-based abort could only ever be seen at a yield, which is precisely
  // what was removed. This entry check stays because generateWorld is still
  // directly callable (the test suite, dev/goldenCompare.html) and the
  // determinism suite asserts an already-aborted signal throws before any work.
  checkAbort(signal);

  // Must equal the number of progress() calls below (the erosion tick fires
  // even when erosion is skipped) so the bar reaches 100%
  const TOTAL_STAGES = 7;
  let stage = 0;
  const progress = () => onProgress?.(++stage, TOTAL_STAGES);

  onLog?.(`Initializing Grid (${params.points} cells)...`);
  progress();
  const macroRng = new RNG(params.seed + '_macro');
  const simplex = new SimplexNoise(new RNG(params.seed));
  
  const points = generateFibonacciSphere(params.points, macroRng, params.cellJitter * 0.8);
  const geoPoints: [number, number][] = points.map(p => {
     const [lat, lon] = toSpherical(p.x, p.y, p.z);
     return [lon, lat]; 
  });
  
  onLog?.("Computing Connectivity...");

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

  let minH = Infinity, maxH = -Infinity, range = 1;

  onLog?.(`Simulating ${params.plates} Tectonic Plates...`);
  progress();

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

  // EROSION
  progress();
  if (params.erosionIterations > 0) {
      onLog?.(`Eroding Terrain (${params.erosionIterations} iter)...`);

      const resFactor = Math.sqrt(params.points / 5000);
      const hydraulicSteps = Math.ceil(params.erosionIterations * 2 * resFactor);
      const thermalSteps = Math.ceil(params.erosionIterations * 0.5 * resFactor);
      
      await applyHydraulicErosion(cells, hydraulicSteps, params.seaLevel); 
      await applyThermalErosion(cells, thermalSteps);
      
      minH = Infinity; maxH = -Infinity;
      cells.forEach(c => { if (c.height < minH) minH = c.height; if (c.height > maxH) maxH = c.height; });
      range = maxH - minH || 1;
      cells.forEach(c => c.height = (c.height - minH) / range);
  }

  // Post-normalization height remap: independently scale peaks and trenches
  // Uses a power curve so heights stay in [0,1] and seaLevel position is preserved
  {
      const mh = params.mountainHeight ?? 1.0;
      const od = params.oceanDepth ?? 1.0;
      // seafloorDepth: linear ocean-floor datum (mean water depth). Multiplies
      // the reshaped depth, holding the coastline fixed. 1.0 is a no-op (shaped
      // ∈ [0,1], so min(1, shaped·1) === shaped) → byte-identical to before.
      const sd = params.seafloorDepth ?? 1.0;
      if (mh !== 1.0 || od !== 1.0 || sd !== 1.0) {
          const sl = params.seaLevel;
          cells.forEach(c => {
              if (c.height >= sl) {
                  const t = (c.height - sl) / Math.max(1e-6, 1 - sl);
                  c.height = sl + (1 - sl) * Math.pow(Math.max(0, t), 1 / Math.max(0.1, mh));
              } else {
                  const t = (sl - c.height) / Math.max(1e-6, sl);
                  const shaped = Math.pow(Math.max(0, t), 1 / Math.max(0.1, od));
                  c.height = sl - sl * Math.min(1, shaped * sd);
              }
          });
      }
  }

  onLog?.("Calculating Climate (Wind & Rain)...");
  progress();
  
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

  progress();
  const { rivers, lakes } = await generateRivers(cells, params.seaLevel, params, onLog);
  const world: WorldData = { cells, params, geoJson: polygons, rivers, lakes };

  // Named geographic features are terrain-derived (B3) — detect them before
  // civs so they never depend on, and are never mutated by, the civ passes.
  world.features = detectFeatures(world);

  progress();
  return recalculateCivs(world, params, onLog);
}

// Biome groups used by the culture terrain-affinity profiles below.
const FOREST_BIOMES = new Set<BiomeType>([
  BiomeType.TROPICAL_RAINFOREST,
  BiomeType.TEMPERATE_FOREST,
  BiomeType.TEMPERATE_RAINFOREST,
  BiomeType.BOREAL_FOREST,
]);
const DESERT_BIOMES = new Set<BiomeType>([
  BiomeType.HOT_DESERT,
  BiomeType.COLD_DESERT,
  BiomeType.STEPPE,
]);

// Looks up the naming style a cell's culture speaks, falling back to the
// world's base nameStyle for cells with no culture assignment yet (e.g. an
// oceanic world where recalculateCultures found zero land candidates).
function getCultureNameStyle(world: WorldData, params: WorldParams, cellId: number): NameStyle {
  const cell = world.cells[cellId];
  const culture = cell?.cultureId !== undefined ? world.cultures?.find(c => c.id === cell.cultureId) : undefined;
  return culture?.nameStyle ?? (params.nameStyle ?? 'fantasy');
}

// One NameGenerator per NameStyle, lazily created on first use and cached for
// the life of the caller. Each style gets its own RNG side-stream (seedPrefix
// + style) so faction/province naming stays fully deterministic regardless of
// which cultures happen to touch which capitals/towns.
function makeStyleNameGenCache(seedPrefix: string): (style: NameStyle) => NameGenerator {
  const cache = new Map<NameStyle, NameGenerator>();
  return (style: NameStyle) => {
    let gen = cache.get(style);
    if (!gen) {
      const rng = new RNG(seedPrefix + style);
      gen = createNameGenerator(style, () => rng.next());
      cache.set(style, gen);
    }
    return gen;
  };
}

// C1: the culture layer. Distinct from (and computed before) factions —
// culture regions spread by terrain affinity and give faction/province/town
// naming a consistent regional flavor. Deterministic from civSeed alone.
//
// DETERMINISM CONSTRAINT: this function must never read from or advance
// civRng (created afterwards in recalculateCivs) or any other existing RNG
// stream. It owns a single new side-stream, civSeed + '_cultures', so every
// existing seed keeps byte-identical faction/province geometry — only the
// generated names change (see recalculateCivs/recalculateProvinces below).
export function recalculateCultures(world: WorldData, params: WorldParams, onLog?: (msg: string) => void): WorldData {
    const numCultures = Math.max(1, params.numCultures || 4);
    onLog?.(`Weaving ${numCultures} Cultures...`);

    world.cells.forEach(c => { c.cultureId = undefined; });

    const cultureRng = new RNG(params.civSeed + '_cultures');

    // Land = same definition used everywhere else in the civ passes.
    const candidates = world.cells.filter(c => c.height >= params.seaLevel && !isLakeCell(c));
    if (candidates.length === 0) {
        world.cultures = [];
        return world;
    }

    // --- Seed home cells, spaced like faction capitals (squared-chord minimum) ---
    // No dedicated spacing slider exists yet for cultures, so a fixed mid
    // spacing factor (matching the capitalSpacing default) is used.
    const CULTURE_SPACING = 0.5;
    const minChordSq = CULTURE_SPACING ** 2 * (4 / numCultures);
    const homes: number[] = [];
    let attempts = 0;
    while (homes.length < numCultures && attempts < 1000) {
        attempts++;
        const candidate = candidates[Math.floor(cultureRng.next() * candidates.length)];
        let tooClose = false;
        for (const homeId of homes) {
            const h = world.cells[homeId];
            const d = (candidate.center.x - h.center.x) ** 2 + (candidate.center.y - h.center.y) ** 2 + (candidate.center.z - h.center.z) ** 2;
            if (d < minChordSq) { tooClose = true; break; }
        }
        if (!tooClose) homes.push(candidate.id);
    }

    // Coastal land cells, precomputed once — reused both by the coastal
    // affinity profile and by classifying a culture's own home cell.
    const coastalIds = new Set<number>();
    for (const c of candidates) {
        for (const nId of c.neighbors) {
            const n = world.cells[nId];
            if (n.height < params.seaLevel || isLakeCell(n)) { coastalIds.add(c.id); break; }
        }
    }

    // Culture 0 always gets params.nameStyle (keeps the nameStyle param live
    // even at numCultures=1); the rest rotate through NAME_STYLES from there.
    const startIdx = NAME_STYLES.indexOf(params.nameStyle ?? 'fantasy');
    const cultures: CultureData[] = homes.map((homeId, i) => {
        const style = NAME_STYLES[(startIdx + i) % NAME_STYLES.length];
        const gen = createNameGenerator(style, () => cultureRng.next());
        return {
            id: i,
            name: gen.faction(),
            color: CULTURE_COLORS[i % CULTURE_COLORS.length],
            nameStyle: style,
            homeCellId: homeId,
        };
    });

    // --- Per-culture terrain-affinity cost profile ---
    // Each profile returns an additive cost term for a candidate cell during
    // expansion (summed with slope + roughness below): a low value pulls the
    // frontier toward that terrain, a higher value mildly resists it so the
    // bias actually shapes the region instead of being lost in noise. Priority
    // order mirrors the design: coastal > highland > forest/desert > generic.
    const affinityFns: Array<(cell: Cell) => number> = homes.map((homeId) => {
        const home = world.cells[homeId];
        if (coastalIds.has(homeId)) {
            // Coastal culture: cheap to hug the shoreline, pricier inland.
            return (cell: Cell) => (coastalIds.has(cell.id) ? 0.4 : 1.1);
        }
        if (home.height > 0.7) {
            // Highland culture: cheap at altitude, pricier in the lowlands.
            return (cell: Cell) => (cell.height > 0.7 ? 0.4 : 1.1);
        }
        if (FOREST_BIOMES.has(home.biome)) {
            // Forest culture: cheap in matching forest biomes.
            return (cell: Cell) => (FOREST_BIOMES.has(cell.biome) ? 0.4 : 1.1);
        }
        if (DESERT_BIOMES.has(home.biome)) {
            // Desert/steppe culture: cheap in matching arid biomes.
            return (cell: Cell) => (DESERT_BIOMES.has(cell.biome) ? 0.4 : 1.1);
        }
        // Generic culture: no terrain preference — flat cost roughly midway
        // between the "favored" (0.4) and "disfavored" (1.1) cases above so
        // it neither dominates nor loses every frontier race.
        return () => 0.75;
    });

    // --- Multi-source Dijkstra over the full cell graph ---
    // Water is impassable except a bounded coastal hop so islands still pick
    // up a culture; land cells beyond the hop's reach are caught by the
    // straggler snap pass below instead of chaining indefinitely across water.
    const pq = new MinHeap<{ id: number; cost: number; culture: number }>(x => x.cost);
    const costs = new Map<number, number>();
    homes.forEach((homeId, idx) => {
        pq.push({ id: homeId, cost: 0, culture: idx });
        costs.set(homeId, 0);
    });

    const waterHopCost = (params.waterCrossingCost || 0.5) * 2;
    const maxWaterCost = waterHopCost * 3;

    while (pq.size() > 0) {
        const { id, cost, culture } = pq.pop()!;
        const cell = world.cells[id];
        const isWater = cell.height < params.seaLevel || isLakeCell(cell);
        if (!isWater) {
            if (cell.cultureId !== undefined && cell.cultureId !== culture) continue;
            cell.cultureId = culture;
        } else if (cost > maxWaterCost) {
            continue;
        }
        for (const nId of cell.neighbors) {
            const nCell = world.cells[nId];
            const nIsWater = nCell.height < params.seaLevel || isLakeCell(nCell);
            let moveCost: number;
            if (nIsWater) {
                moveCost = waterHopCost;
            } else {
                const slope = Math.abs(nCell.height - cell.height);
                moveCost = slope + affinityFns[culture](nCell) + cultureRng.next() * 0.3;
            }
            const newCost = cost + moveCost;
            if (nIsWater && newCost > maxWaterCost) continue;
            if (!costs.has(nId) || newCost < costs.get(nId)!) {
                costs.set(nId, newCost);
                if (nCell.cultureId === undefined) {
                    pq.push({ id: nId, cost: newCost, culture });
                }
            }
        }
    }

    // Stragglers: land cells the bounded hop never reached (isolated islands)
    // snap to whichever assigned cell is nearest by chord distance, so every
    // land cell ends this pass with a cultureId.
    const unassigned = candidates.filter(c => c.cultureId === undefined);
    if (unassigned.length > 0) {
        const assignedCells = candidates.filter(c => c.cultureId !== undefined);
        for (const c of unassigned) {
            let bestCultureId: number | undefined;
            let bestDistSq = Infinity;
            for (const a of assignedCells) {
                const d = (c.center.x - a.center.x) ** 2 + (c.center.y - a.center.y) ** 2 + (c.center.z - a.center.z) ** 2;
                if (d < bestDistSq) { bestDistSq = d; bestCultureId = a.cultureId; }
            }
            if (bestCultureId !== undefined) c.cultureId = bestCultureId;
        }
    }

    world.cultures = cultures;
    return world;
}

// Naming templates for the religion layer (C2). Each takes the fresh name
// drawn from the relevant culture's style generator (X) and, where useful,
// the culture's own name; which template is picked is itself an rng draw so
// two runs of the same seed still agree on phrasing.
const FOLK_TEMPLATES: Array<(x: string, cultureName: string) => string> = [
  (x) => `${x} Rites`,
  (_x, cultureName) => `Old Faith of ${cultureName}`,
  (x) => `${x} Creed`,
];
const ORGANIZED_TEMPLATES: Array<(x: string) => string> = [
  (x) => `Church of ${x}`,
  (x) => `${x} Path`,
  (x) => `Temple of ${x}`,
];

// C2: the religion layer. Two kinds share the world: exactly one 'folk'
// religion per culture (covering every land cell by default) and a handful
// of 'organized' religions that spread outward from holy cities (the
// highest-population towns) and convert folk cells within a limited budget.
// Called at the end of recalculateProvinces (see below) because holy-city
// selection needs towns/populations, which only exist once provinces have
// run — recalculateCivs's own last statement is `return
// recalculateProvinces(...)`, so hooking there still covers every civ
// recalculation, plus the standalone "Update Provinces" path in App.tsx.
//
// DETERMINISM CONSTRAINT: mirrors recalculateCultures — this function owns
// a single new side-stream, civSeed + '_religions', and never reads from or
// advances civRng/provRng/cultureRng or any cell field other than
// religionId, so numCultures/numFactions/provinceSize signatures
// (terrainSignature/civSignature/cultureSignature in tests/helpers.ts) stay
// byte-identical — only world.religions and cell.religionId are new.
export function recalculateReligions(world: WorldData, params: WorldParams, onLog?: (msg: string) => void): WorldData {
    onLog?.('Kindling Religions...');

    world.cells.forEach(c => { c.religionId = undefined; });

    if (!world.cultures || world.cultures.length === 0) {
        world.religions = [];
        return world;
    }

    const religionRng = new RNG(params.civSeed + '_religions');

    // Land = same definition used everywhere else in the civ passes.
    const candidates = world.cells.filter(c => c.height >= params.seaLevel && !isLakeCell(c));
    if (candidates.length === 0) {
        world.religions = [];
        return world;
    }

    // Names are drawn directly from this single religionRng stream rather
    // than spinning a fresh per-style sub-stream (contrast
    // makeStyleNameGenCache, used for faction/province naming) — the design
    // calls for "a fresh name from the culture's style generator on the
    // religions stream," so every religion name, of any style, shares one
    // deterministic draw order.
    const relNameGenCache = new Map<NameStyle, NameGenerator>();
    const getRelNameGen = (style: NameStyle): NameGenerator => {
        let gen = relNameGenCache.get(style);
        if (!gen) {
            gen = createNameGenerator(style, () => religionRng.next());
            relNameGenCache.set(style, gen);
        }
        return gen;
    };

    // --- Folk religions: exactly one per culture ---
    // Folk id === culture id (cultures are indexed 0..n-1 in creation
    // order), which lets the final stamping pass below fall back to
    // `cell.cultureId` directly as the matching folk religion id.
    const folkReligions: ReligionData[] = world.cultures.map((culture) => {
        const gen = getRelNameGen(culture.nameStyle);
        const template = FOLK_TEMPLATES[Math.floor(religionRng.next() * FOLK_TEMPLATES.length)];
        return {
            id: culture.id,
            name: template(gen.faction(), culture.name),
            kind: 'folk',
            color: FOLK_COLORS[culture.id % FOLK_COLORS.length],
            cultureId: culture.id,
            holyCellId: null,
        };
    });

    // --- Organized religions: seeded at holy cities ---
    // Count derives from numCultures rather than a new param, per design.
    const organizedCount = Math.max(1, Math.floor((params.numCultures || 4) / 2));

    const townCells = world.cells.filter(c => c.isTown);
    let organizedReligions: ReligionData[] = [];
    let holyCityIds: number[] = [];

    if (townCells.length > 0) {
        // rng-tiebreak for equal-population towns, pre-drawn per town (in
        // world.cells id order) so the number of religionRng draws never
        // depends on how the sort algorithm orders its comparisons.
        const tiebreak = new Map<number, number>();
        for (const t of townCells) tiebreak.set(t.id, religionRng.next());
        const sortedTowns = [...townCells].sort((a, b) => {
            const diff = (b.population || 0) - (a.population || 0);
            if (diff !== 0) return diff;
            return tiebreak.get(a.id)! - tiebreak.get(b.id)!;
        });

        // Spaced like faction capitals / culture homes (squared-chord minimum).
        const HOLY_CITY_SPACING = 0.5;
        const minChordSq = HOLY_CITY_SPACING ** 2 * (4 / organizedCount);
        for (const cell of sortedTowns) {
            if (holyCityIds.length >= organizedCount) break;
            let tooClose = false;
            for (const hId of holyCityIds) {
                const h = world.cells[hId];
                const d = (cell.center.x - h.center.x) ** 2 + (cell.center.y - h.center.y) ** 2 + (cell.center.z - h.center.z) ** 2;
                if (d < minChordSq) { tooClose = true; break; }
            }
            if (!tooClose) holyCityIds.push(cell.id);
        }
        // Fallback for maps with too few well-spaced towns: fill remaining
        // slots from the same population-ranked list, ignoring spacing, so
        // organizedCount is honored whenever at least that many towns exist.
        if (holyCityIds.length < organizedCount) {
            for (const cell of sortedTowns) {
                if (holyCityIds.length >= organizedCount) break;
                if (!holyCityIds.includes(cell.id)) holyCityIds.push(cell.id);
            }
        }

        organizedReligions = holyCityIds.map((holyCellId, i) => {
            const style = getCultureNameStyle(world, params, holyCellId);
            const gen = getRelNameGen(style);
            const template = ORGANIZED_TEMPLATES[Math.floor(religionRng.next() * ORGANIZED_TEMPLATES.length)];
            return {
                id: folkReligions.length + i,
                name: template(gen.faction()),
                kind: 'organized',
                color: RELIGION_COLORS[i % RELIGION_COLORS.length],
                cultureId: null,
                holyCellId,
            };
        });
    }

    // --- Multi-source Dijkstra spread from holy cities over land ---
    // (bounded water hops copied from recalculateCultures so organized
    // faiths can still reach nearby islands, same as culture regions do).
    const converted = new Map<number, number>(); // land cellId -> organizedReligions index

    if (organizedReligions.length > 0) {
        const holyCultureIds = holyCityIds.map(id => world.cells[id].cultureId);

        // Per-religion expansion budget is a CELL-COUNT quota, not a cost
        // ceiling: Dijkstra cost grows with map size (more hops to cross a
        // bigger map), so a cost-based cap would either convert nothing or
        // everything depending on map scale. A count quota scales with the
        // map instead, guaranteeing roughly half of all land stays folk at
        // any resolution. Base ~ half of an equal per-religion share of
        // land cells, varied 0.7-1.3x via religionRng so organized faiths
        // don't all claim identically-sized territories.
        const baseBudget = (candidates.length / organizedReligions.length) * 0.5;
        const budgets = organizedReligions.map(() =>
            Math.max(1, Math.round(baseBudget * (0.7 + religionRng.next() * 0.6)))
        );
        const convertedCount = organizedReligions.map(() => 0);

        const pq = new MinHeap<{ id: number; cost: number; religion: number }>(x => x.cost);
        const costs = new Map<number, number>();
        holyCityIds.forEach((cellId, idx) => {
            pq.push({ id: cellId, cost: 0, religion: idx });
            costs.set(cellId, 0);
        });

        const waterHopCost = (params.waterCrossingCost || 0.5) * 2;
        const maxWaterCost = waterHopCost * 3;

        while (pq.size() > 0) {
            const { id, cost, religion } = pq.pop()!;
            const cell = world.cells[id];
            const isWater = cell.height < params.seaLevel || isLakeCell(cell);

            if (isWater) {
                if (cost > maxWaterCost) continue;
            } else {
                if (converted.has(id)) continue;
                // Budget exhausted: this religion's frontier stops
                // expanding from here, so the cell (and everything beyond
                // it) stays folk. Deliberate contrast with C1 cultures,
                // which snap every straggler — organized faiths plateau.
                if (convertedCount[religion] >= budgets[religion]) continue;
                converted.set(id, religion);
                convertedCount[religion]++;
            }

            for (const nId of cell.neighbors) {
                const nCell = world.cells[nId];
                const nIsWater = nCell.height < params.seaLevel || isLakeCell(nCell);
                if (!nIsWater && converted.has(nId)) continue;
                let moveCost: number;
                if (nIsWater) {
                    moveCost = waterHopCost;
                } else {
                    const slope = Math.abs(nCell.height - cell.height);
                    // Cheaper to spread within the holy city's own culture
                    // region — an organized faith travels easiest among
                    // its founding people.
                    const withinHomeCulture = nCell.cultureId === holyCultureIds[religion];
                    moveCost = slope + (withinHomeCulture ? 0.35 : 0.6) + religionRng.next() * 0.3;
                }
                const newCost = cost + moveCost;
                if (nIsWater && newCost > maxWaterCost) continue;
                if (!costs.has(nId) || newCost < costs.get(nId)!) {
                    costs.set(nId, newCost);
                    pq.push({ id: nId, cost: newCost, religion });
                }
            }
        }
    }

    // --- Stamp every land cell: converted -> its organized religion, else folk ---
    // (folk id === culture id, so cell.cultureId is the direct fallback).
    for (const cell of candidates) {
        const relIdx = converted.get(cell.id);
        cell.religionId = relIdx !== undefined ? organizedReligions[relIdx].id : cell.cultureId;
    }

    world.religions = [...folkReligions, ...organizedReligions];
    return world;
}

// ... (recalculateCivs and recalculateProvinces remain unchanged as they are synchronous)
export function recalculateCivs(world: WorldData, params: WorldParams, onLog?: (msg: string) => void): WorldData {
    // Cultures are computed first and are purely additive: recalculateCultures
    // owns its own RNG stream and never touches civRng below, so existing
    // seeds keep byte-identical faction/province geometry.
    recalculateCultures(world, params, onLog);

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
    // Faction names draw from a per-NameStyle cache, each on its own RNG
    // side-stream (civSeed + '_facnames_' + style) — so naming a faction
    // after its capital's culture never touches civRng and stays fully
    // deterministic regardless of culture layout.
    const getFacNameGen = makeStyleNameGenCache(params.civSeed + '_facnames_');
    const numFactions = params.numFactions;
    const factions: FactionData[] = [];
    
    const suitable = world.cells.filter(c =>
        c.height >= params.seaLevel &&
        !isLakeCell(c) &&
        c.biome !== BiomeType.ICE_CAP &&
        c.biome !== BiomeType.VOLCANIC
    );

    const candidates = suitable.length > 0 ? suitable : world.cells.filter(c => c.height >= params.seaLevel && !isLakeCell(c));
    
    if (candidates.length === 0) {
        world.civData = { factions: [] };
        return world;
    }

    const capitals: number[] = [];
    // Scale-independent minimum separation: at spacing 1.0 the squared-chord
    // threshold is 4/numFactions (approaching an even spread over the sphere);
    // at 0 capitals may land anywhere. The old formula scaled with cell count
    // and effectively never rejected candidates on smaller maps.
    const minChordSq = (params.capitalSpacing ?? 0.5) ** 2 * (4 / numFactions);

    let attempts = 0;
    while(capitals.length < numFactions && attempts < 1000) {
        attempts++;
        const candidate = candidates[Math.floor(civRng.next() * candidates.length)];

        let tooClose = false;
        for(const capId of capitals) {
            const cap = world.cells[capId];
            const d = (candidate.center.x - cap.center.x)**2 + (candidate.center.y - cap.center.y)**2 + (candidate.center.z - cap.center.z)**2;
            if (d < minChordSq) {
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
                name: getFacNameGen(getCultureNameStyle(world, params, candidate.id)).faction(),
                color: '#ffffff',
                capitalId: candidate.id,
                provinces: [],
                totalPopulation: 0
            });
        }
    }

    // civSizeVariance: per-faction competitive movement-cost scaling.
    // A faction with a larger size factor pays proportionally less per cell,
    // so it wins Dijkstra frontier races against neighbors and grows bigger.
    // Cost scaling (rather than an absolute budget) works at any map
    // resolution because expansion is competitive, not distance-capped.
    const sizeVariance = params.civSizeVariance ?? 0;
    const costMult = capitals.map(() =>
        1 / Math.min(2, Math.max(0.25, 1 + (civRng.next() * 2 - 1) * sizeVariance))
    );

    const pq = new MinHeap<{id: number, cost: number, region: number}>(x => x.cost);
    const costs = new Map<number, number>();

    capitals.forEach((capId, idx) => {
        pq.push({ id: capId, cost: 0, region: idx });
        costs.set(capId, 0);
    });

    const waterCost = (params.waterCrossingCost || 0.5) * 50;
    const territorialRange = (params.territorialWaters || 0.2) * 50;

    while(pq.size() > 0) {
        const { id, cost, region } = pq.pop()!;
        if (world.cells[id].regionId !== undefined && world.cells[id].regionId !== region) continue;
        world.cells[id].regionId = region;
        if (cost > 200) continue;
        const currCell = world.cells[id];
        for(const nId of currCell.neighbors) {
            const nCell = world.cells[nId];
            let moveCost = landTerrainStepCost(currCell, nCell);
            const isWater = nCell.height < params.seaLevel || isLakeCell(nCell);
            if (isWater) moveCost = waterCost;
            moveCost *= (1 + (civRng.next() * params.borderRoughness));
            moveCost *= costMult[region];
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
    
    // Shuffle FACTION_COLORS with civRng so colors are varied but deterministic per seed
    const shuffled = [...FACTION_COLORS];
    for (let i = shuffled.length - 1; i > 0; i--) {
        const j = Math.floor(civRng.next() * (i + 1));
        [shuffled[i], shuffled[j]] = [shuffled[j], shuffled[i]];
    }
    factions.forEach((f, i) => f.color = shuffled[i % shuffled.length]);
    
    world.civData = { factions };
    
    return recalculateProvinces(world, params);
}

export function recalculateProvinces(world: WorldData, params: WorldParams): WorldData {
    if (!world.civData) return world;
    const provRng = new RNG(params.civSeed + '_prov');
    // Province/town names draw from a per-NameStyle cache, keyed by each
    // town cell's own culture, on dedicated RNG side-streams (civSeed +
    // '_provnames_' + style) — separate from recalculateCivs's faction
    // stream and from provRng, so naming never perturbs geometry.
    const getProvNameGen = makeStyleNameGenCache(params.civSeed + '_provnames_');

    world.cells.forEach(c => {
        let suitability = 0;
        if (c.height < params.seaLevel || isLakeCell(c)) {
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
            const n = world.cells[nId];
            if (n.height < params.seaLevel || isLakeCell(n)) { coast = true; break; }
        }
        if (coast) suitability += 0.3;
        if (c.height > 0.6) suitability -= (c.height - 0.6) * 2;
        // High peaks could drive suitability negative, producing negative
        // populations that silently deflated faction totals
        suitability = Math.max(0, suitability);
        c.population = Math.floor(suitability * 10000 * (0.8 + provRng.next() * 0.4));
    });

    world.civData.factions.forEach(faction => {
        faction.provinces = [];
        faction.totalPopulation = 0;
        const landCells = world.cells.filter(c => c.regionId === faction.id && c.height >= params.seaLevel && !isLakeCell(c));
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
            const nameGen = getProvNameGen(getCultureNameStyle(world, params, tId));
            faction.provinces.push({
                id: idx,
                name: nameGen.province(),
                towns: [{ name: nameGen.town(), cellId: tId, population: tCell.population || 0, isCapital: tId === faction.capitalId }],
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

    // C2: religions need towns/populations for holy-city selection, which
    // this function just finished computing above — see the determinism
    // and hook-point comment on recalculateReligions itself.
    recalculateReligions(world, params);

    // C3: routes are derived but computed LAZILY at the App level, gated on the
    // Roads & Routes toggle — computeRoutes is O(towns · A*) and runs several
    // seconds near the 200k-cell cap, so a routes-off generation must pay zero.
    // Clearing here invalidates any stale routes after a civ/province recalc so
    // the lazy pass recomputes against the new town graph.
    world.routes = undefined;

    return world;
}
