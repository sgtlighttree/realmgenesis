import { geoVoronoi } from 'd3-geo-voronoi';
import { Cell, Point, WorldData, WorldParams, BiomeType, CivData, FactionData, ProvinceData, TownData } from '../types';
import { RNG, SimplexNoise } from './rng';
import { BIOME_COLORS } from './colors';

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

function findClosestCellIndex(cells: Cell[], p: Point): number {
    let minDist = Infinity;
    let index = -1;
    for (let i = 0; i < cells.length; i++) {
        const c = cells[i].center;
        const d = (c.x - p.x)**2 + (c.y - p.y)**2 + (c.z - p.z)**2;
        if (d < minDist) {
            minDist = d;
            index = i;
        }
    }
    return index;
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

function domainWarp(simplex: SimplexNoise, x: number, y: number, z: number, strength: number): number {
    if (strength <= 0.01) return fbm(simplex, x, y, z, 4, 0.5, 2.0);

    const qx = fbm(simplex, x, y, z, 2, 0.5, 2.0);
    const qy = fbm(simplex, x + 5.2, y + 1.3, z + 2.8, 2, 0.5, 2.0);
    const qz = fbm(simplex, x + 9.2, y + 2.8, z + 11.0, 2, 0.5, 2.0);

    return fbm(simplex, x + strength * qx, y + strength * qy, z + strength * qz, 6, 0.5, 2.0);
}

// --- EROSION ---

function applyHydraulicErosion(cells: Cell[], iterations: number) {
    cells.forEach(c => c.flux = 0);
    const sorted = [...cells].sort((a, b) => b.height - a.height);
    const erosionRate = 0.025;
    const depositionRate = 0.002;
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
                const streamPower = c.flux! * slope;
                const erosion = streamPower * erosionRate;
                const safeErosion = Math.min(erosion, slope * 0.9);
                c.height -= safeErosion;
                target.height += safeErosion * depositionRate; 
            }
        });
    }
}

function applyThermalErosion(cells: Cell[], iterations: number) {
    const talus = 0.01; 
    const rate = 0.1; 
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

// --- GEOGRAPHY GENERATION ---

export async function generateGeography(params: WorldParams, onProgress?: (msg: string, pct: number) => void): Promise<WorldData> {
  onProgress?.("Initializing Grid...", 10);
  const macroRng = new RNG(params.seed + '_macro');
  const simplex = new SimplexNoise(new RNG(params.seed));
  const points = generateFibonacciSphere(params.points, macroRng, params.cellJitter);
  const geoPoints: [number, number][] = points.map(p => {
     const [lat, lon] = toSpherical(p.x, p.y, p.z);
     return [lon, lat]; 
  });
  onProgress?.("Computing Voronoi...", 20);
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
  const plateVectors: Point[] = [];
  for(let i=0; i<numPlates; i++) plateVectors.push(randomVector(macroRng));
  const plateCenters: number[] = plateVectors.map(v => findClosestCellIndex(cells, v));
  const plateDrift = plateVectors.map(() => ({ x: macroRng.next() - 0.5, y: macroRng.next() - 0.5, z: macroRng.next() - 0.5 }));
  
  const queue: {id: number, plateIdx: number}[] = [];
  cells.forEach(c => c.plateId = -1);
  plateCenters.forEach((id, idx) => {
    if (id >= 0 && id < cells.length) {
        cells[id].plateId = idx;
        queue.push({ id, plateIdx: idx });
    }
  });
  let head = 0;
  while(head < queue.length) {
     const { id, plateIdx } = queue[head++];
     const cell = cells[id];
     for (const nId of cell.neighbors) {
        if (cells[nId].plateId === -1) {
             cells[nId].plateId = plateIdx;
             queue.push({ id: nId, plateIdx });
        }
     }
  }
  cells.forEach(c => {
    if (c.plateId === -1 && c.neighbors.length > 0) {
        c.plateId = cells[c.neighbors[0]].plateId;
    }
  });

  const cellStress = new Float32Array(cells.length).fill(0); 
  cells.forEach(c => {
      let maxStress = 0;
      for (const nId of c.neighbors) {
          const n = cells[nId];
          if (n.plateId !== c.plateId) {
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
      cellStress[c.id] = maxStress;
  });

  const nextStress = new Float32Array(cells.length);
  for(let i=0; i<3; i++) {
      cells.forEach(c => {
          let sum = cellStress[c.id];
          c.neighbors.forEach(nId => sum += cellStress[nId]);
          nextStress[c.id] = sum / (c.neighbors.length + 1);
      });
      nextStress.forEach((v,k) => cellStress[k] = v);
  }

  onProgress?.("Generating Terrain...", 60);
  const featureFreq = params.noiseScale || 1.0;
  const warpStr = (params.warpStrength === undefined ? 0.5 : params.warpStrength) * 0.5;
  const plateInf = (params.plateInfluence === undefined ? 0.5 : params.plateInfluence) * 2.0; 
  const userRidgeBias = params.ridgeBlend; 

  cells.forEach(c => {
      let continental = fbm(simplex, c.center.x * featureFreq, c.center.y * featureFreq, c.center.z * featureFreq, 2, 0.5, 2.0);
      const stress = cellStress[c.id];
      const stressFactor = Math.min(1, Math.max(0, stress * plateInf)); 
      const mix = stressFactor + (1 - stressFactor) * userRidgeBias;
      const localRoughness = params.roughness + stressFactor * 0.5;
      const warpedX = c.center.x * 2.0;
      const warpedY = c.center.y * 2.0;
      const warpedZ = c.center.z * 2.0;
      const fbmNoise = domainWarp(simplex, warpedX, warpedY, warpedZ, warpStr);
      const ridgeNoise = ridgedNoise(simplex, warpedX, warpedY, warpedZ, 5, 2.0) * 2 - 0.5;
      let detail = fbmNoise * (1 - mix) + ridgeNoise * mix;
      if (params.points > 5000) {
          const micro = simplex.noise3D(c.center.x * 12, c.center.y * 12, c.center.z * 12);
          detail += micro * 0.1 * params.roughness;
      }
      let height = continental; 
      height += detail * localRoughness * 0.5;
      if (stress > 0) height += stress * 0.8; 
      else if (stress < 0) height += stress * 1.5;
      if (params.maskType === 'Pangea') {
          const mask = (c.center.x * 0.8 + c.center.y * 0.2 + 1) * 0.5;
          const smoothMask = mask * mask * (3 - 2 * mask);
          height = height * 0.5 + smoothMask * 0.8 - 0.3;
      }
      c.height = height;
  });
  
  let minH = Infinity, maxH = -Infinity;
  cells.forEach(c => { if (c.height < minH) minH = c.height; if (c.height > maxH) maxH = c.height; });
  let range = maxH - minH || 1;
  cells.forEach(c => c.height = (c.height - minH) / range);

  if (params.erosionIterations > 0) {
      onProgress?.("Eroding Terrain...", 70);
      await new Promise(r => setTimeout(r, 0));
      applyHydraulicErosion(cells, Math.ceil(params.erosionIterations * 2));
      applyThermalErosion(cells, Math.ceil(params.erosionIterations / 3));
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

// --------------------------------------------------------
// --- POLITICAL & PROVINCE GENERATION ---
// --------------------------------------------------------

const PLATE_COLORS = ['#ef5350', '#ab47bc', '#7e57c2', '#5c6bc0', '#42a5f5', '#29b6f6', '#26c6da', '#26a69a', '#66bb6a', '#9ccc65', '#d4e157', '#ffee58', '#ffca28', '#ffa726', '#ff7043', '#8d6e63', '#bdbdbd', '#78909c'];

export function recalculateProvinces(world: WorldData, params: WorldParams): WorldData {
    const civRng = new RNG(Math.random().toString(36));
    const numFactions = params.numFactions;
    const factions: FactionData[] = [];
    const factionCells = Array.from({length: numFactions}, () => [] as number[]);
    const capitalIds = new Array(numFactions).fill(-1);

    world.cells.forEach(c => {
        c.provinceId = undefined;
        c.isTown = false;
        if (c.regionId !== undefined) {
            factionCells[c.regionId].push(c.id);
            if (c.isCapital) capitalIds[c.regionId] = c.id;
        }
    });

    const provinceSizeParam = params.provinceSize || 0.5;

    factionCells.forEach((cellsInFaction, fId) => {
        if (cellsInFaction.length === 0) return;
        const existingFaction = world.civData?.factions.find(f => f.id === fId);
        const factionName = existingFaction ? existingFaction.name : `Faction ${fId + 1}`;
        const factionColor = existingFaction ? existingFaction.color : PLATE_COLORS[fId % PLATE_COLORS.length];
        const targetCellsPerProvince = 50 + (provinceSizeParam * 450); 
        const numProvinces = Math.max(1, Math.floor(cellsInFaction.length / targetCellsPerProvince));
        const provCenters: number[] = [];
        const capId = capitalIds[fId] !== -1 ? capitalIds[fId] : cellsInFaction[0];
        provCenters.push(capId);
        
        for(let k=1; k<numProvinces; k++) {
            let bestC = -1; let maxMinDist = -1;
            for(let j=0; j<20; j++) {
                const candId = cellsInFaction[Math.floor(civRng.next() * cellsInFaction.length)];
                let minDist = Infinity;
                for(const pc of provCenters) {
                    const dist = (world.cells[candId].center.x - world.cells[pc].center.x)**2 + (world.cells[candId].center.y - world.cells[pc].center.y)**2 + (world.cells[candId].center.z - world.cells[pc].center.z)**2;
                    if (dist < minDist) minDist = dist;
                }
                if (minDist > maxMinDist) { maxMinDist = minDist; bestC = candId; }
            }
            if (bestC !== -1) provCenters.push(bestC);
        }
        
        const provincesData: ProvinceData[] = provCenters.map((pid, i) => ({
            id: i,
            name: (existingFaction && existingFaction.provinces[i]) ? existingFaction.provinces[i].name : `Province ${i+1}`,
            towns: [],
            totalPopulation: 0,
        }));

        cellsInFaction.forEach(cid => {
             const c = world.cells[cid];
             let minDist = Infinity; let pIdx = 0;
             provCenters.forEach((pid, idx) => {
                 const pCell = world.cells[pid];
                 const d = (c.center.x - pCell.center.x)**2 + (c.center.y - pCell.center.y)**2 + (c.center.z - pCell.center.z)**2;
                 if (d < minDist) { minDist = d; pIdx = idx; }
             });
             c.provinceId = pIdx;
             provincesData[pIdx].totalPopulation += (c.population || 0);
        });
        
        provincesData.forEach(prov => {
            let maxPop = -1; let townId = -1;
            cellsInFaction.forEach(cid => {
                const c = world.cells[cid];
                // CRITICAL: Force towns to be on land (height >= seaLevel)
                if (c.provinceId === prov.id && c.height >= params.seaLevel) {
                    if ((c.population || 0) > maxPop) { maxPop = c.population || 0; townId = cid; }
                }
            });
            // Fallback to highest pop cell if whole province is underwater (unlikely but possible)
            if (townId === -1 && cellsInFaction.length > 0) {
               cellsInFaction.forEach(cid => {
                  const c = world.cells[cid];
                  if (c.provinceId === prov.id && (c.population || 0) > maxPop) { maxPop = c.population || 0; townId = cid; }
               });
            }
            if (townId !== -1) {
                world.cells[townId].isTown = true;
                const isCap = world.cells[townId].isCapital || false;
                let tName = isCap ? `Capital` : `Town`;
                if (existingFaction && existingFaction.provinces[prov.id] && existingFaction.provinces[prov.id].towns.length > 0) {
                     const oldT = existingFaction.provinces[prov.id].towns.find(t => t.isCapital === isCap);
                     if (oldT) tName = oldT.name;
                }
                prov.towns.push({ name: tName, cellId: townId, population: maxPop, isCapital: isCap });
            }
        });
        factions.push({ id: fId, name: factionName, color: factionColor, capitalId: capId, provinces: provincesData, totalPopulation: provincesData.reduce((sum, p) => sum + p.totalPopulation, 0), description: existingFaction?.description });
    });
    world.civData = { factions };
    return world;
}

export function recalculateCivs(world: WorldData, params: WorldParams, onProgress?: (msg: string, pct: number) => void): WorldData {
    onProgress?.("Simulating History...", 90);
    world.cells.forEach(c => { c.regionId = undefined; c.provinceId = undefined; c.isCapital = false; c.isTown = false; c.population = 0; });
    world.cells.forEach(c => {
        if (c.height < params.seaLevel) return;
        let score = 10;
        if (c.flux && c.flux > 0.5) score += 50; 
        if (c.neighbors.some(n => world.cells[n].height < params.seaLevel)) score += 30;
        if (c.temperature > 0 && c.temperature < 30) score += 20;
        else if (c.temperature < -10 || c.temperature > 40) score -= 20;
        let variance = 0; c.neighbors.forEach(n => variance += Math.abs(c.height - world.cells[n].height));
        if (variance < 0.05) score += 20; else score -= 10;
        c.population = Math.max(0, Math.floor(score * 100));
    });
    const civRng = new RNG(params.civSeed || params.seed + '_civs');
    const numFactions = params.numFactions;
    const capitalIds: number[] = [];
    const factionPowers: number[] = [];
    const sizeVar = params.civSizeVariance === undefined ? 0.5 : params.civSizeVariance;
    const capitalSpacing = params.capitalSpacing === undefined ? 0.5 : params.capitalSpacing;
    const minCapitalDistSq = Math.pow(0.1 + capitalSpacing * 1.0, 2);
    for(let i=0; i<numFactions; i++) factionPowers.push(Math.max(0.2, 1.0 + (civRng.next() - 0.5) * 2 * sizeVar));
    for (let i=0; i<numFactions; i++) {
        let bestId = -1; let maxScore = -Infinity;
        for (let k=0; k<50; k++) {
            const vec = randomVector(civRng); const id = findClosestCellIndex(world.cells, vec);
            if (id === -1) continue;
            const c = world.cells[id];
            if (c.height < params.seaLevel || capitalIds.includes(id)) continue;
            let score = (c.population || 0);
            for (const cap of capitalIds) {
                const dSq = (c.center.x - world.cells[cap].center.x)**2 + (c.center.y - world.cells[cap].center.y)**2 + (c.center.z - world.cells[cap].center.z)**2;
                if (dSq < minCapitalDistSq) score -= 50000;
            }
            if (score > maxScore) { maxScore = score; bestId = id; }
        }
        if (bestId !== -1) { capitalIds.push(bestId); world.cells[bestId].isCapital = true; world.cells[bestId].population! += 50000; }
    }
    const queue: {id: number, regionId: number, cost: number}[] = [];
    const claimed = new Map<number, number>(); 
    capitalIds.forEach((id, idx) => { world.cells[id].regionId = idx; queue.push({ id, regionId: idx, cost: 0 }); claimed.set(id, 0); });
    const noise = new SimplexNoise(civRng);
    const waterCostBase = params.waterCrossingCost === undefined ? 0.8 : params.waterCrossingCost;
    const waterPenalty = 1 + (1.0 - waterCostBase) * 200; 
    while(queue.length > 0) {
        const { id, regionId, cost } = queue.shift()!;
        if (claimed.has(id) && claimed.get(id)! < cost) continue;
        const current = world.cells[id];
        const powerMultiplier = 1.0 / factionPowers[regionId];
        for (const nId of current.neighbors) {
            const neighbor = world.cells[nId];
            let stepCost = 1.0;
            if (neighbor.height < params.seaLevel) stepCost = waterPenalty;
            else { if (neighbor.height > current.height + 0.05) stepCost += 10; if (neighbor.biome === BiomeType.ICE_CAP) stepCost += 15; if (neighbor.biome === BiomeType.HOT_DESERT) stepCost += 8; if (neighbor.biome === BiomeType.TROPICAL_RAINFOREST) stepCost += 5; }
            const nVal = noise.noise3D(neighbor.center.x * 10, neighbor.center.y * 10, neighbor.center.z * 10);
            stepCost += nVal * params.borderRoughness * 10;
            if (stepCost < 1) stepCost = 1;
            const totalCost = cost + (stepCost * powerMultiplier);
            if (!claimed.has(nId) || totalCost < claimed.get(nId)!) { claimed.set(nId, totalCost); neighbor.regionId = regionId; queue.push({ id: nId, regionId, cost: totalCost }); }
        }
        if (queue.length > 0 && queue.length % 500 === 0) queue.sort((a,b) => a.cost - b.cost);
    }
    const result = recalculateProvinces(world, params);
    onProgress?.("Finalizing...", 100);
    return result;
}

export async function generateWorld(params: WorldParams, onProgress?: (msg: string, pct: number) => void): Promise<WorldData> {
    return generateGeography(params, onProgress);
}