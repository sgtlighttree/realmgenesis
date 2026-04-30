import * as d3 from 'd3';
import { geoWinkel3, geoRobinson, geoMollweide } from 'd3-geo-projection';
import { WorldData, ViewMode, WorldParams, CivData, DymaxionSettings } from '../types';
import { getCellColor } from './colors';
import { createDymaxionProjection, buildDymaxionNet } from './dymaxion';

export type ExportResolution = 4096 | 8192 | 16384 | 32768;
export type ProjectionType = 'equirectangular' | 'mercator' | 'winkeltripel' | 'orthographic' | 'robinson' | 'mollweide' | 'dymaxion';
export type DymaxionExportSettings = Pick<DymaxionSettings, 'layout' | 'lon' | 'lat' | 'roll'>;

const renderEquirectangular = (
  world: WorldData,
  viewMode: ViewMode,
  width: number,
  height: number
) => {
  const canvas = document.createElement('canvas');
  canvas.width = width;
  canvas.height = height;
  const ctx = canvas.getContext('2d');
  if (!ctx) return null;

  ctx.fillStyle = viewMode === 'satellite' || viewMode === 'biome' ? '#050505' : '#000000';
  ctx.fillRect(0, 0, width, height);

  ctx.save();
  ctx.translate(width, 0);
  ctx.scale(-1, 1);

  const projection = d3.geoEquirectangular().fitSize([width, height], { type: 'Sphere' } as any);
  const pathGenerator = d3.geoPath(projection, ctx);

  world.cells.forEach((cell, i) => {
    const feature = world.geoJson.features[i];
    if (!feature) return;
    const threeColor = getCellColor(cell, viewMode, world.params.seaLevel);
    ctx.beginPath();
    pathGenerator(feature);
    ctx.fillStyle = '#' + threeColor.getHexString();
    ctx.fill();
  });

  ctx.restore();
  return canvas;
};

const DEBUG_DYMAXION = false;

const exportDymaxionRaster = (
  world: WorldData,
  viewMode: ViewMode,
  width: number,
  height: number,
  dymaxionSettings?: DymaxionExportSettings
) => {
  const srcWidth = width;
  const srcHeight = Math.round(width / 2);
  const source = renderEquirectangular(world, viewMode, srcWidth, srcHeight);
  if (!source) return;
  const srcCtx = source.getContext('2d');
  if (!srcCtx) return;
  const srcImage = srcCtx.getImageData(0, 0, srcWidth, srcHeight);
  const srcData = srcImage.data;

  const canvas = document.createElement('canvas');
  canvas.width = width;
  canvas.height = height;
  const ctx = canvas.getContext('2d');
  if (!ctx) return;

  ctx.fillStyle = viewMode === 'satellite' || viewMode === 'biome' ? '#050505' : '#000000';
  ctx.fillRect(0, 0, width, height);

  const layout = dymaxionSettings?.layout || 'classic';
  const net = buildDymaxionNet(layout);
  const faces = net.faces;
  const isBlender = layout === 'blender';

  // For the Blender UV net, UV coords map directly to pixels:
  //   px = u * width,  py = (1 - v) * height
  // (V is flipped because image y=0 is top, UV v=0 is bottom.)
  // For the classic net, auto-fit the net to fill the canvas with padding.
  let scale = 1;
  let offsetX = 0;
  let offsetY = 0;

  if (!isBlender) {
    let minX = Infinity, minY = Infinity, maxX = -Infinity, maxY = -Infinity;
    faces.forEach((face) => {
      face.vertices.forEach((v) => {
        minX = Math.min(minX, v[0]);
        minY = Math.min(minY, v[1]);
        maxX = Math.max(maxX, v[0]);
        maxY = Math.max(maxY, v[1]);
      });
    });
    const pad = 12;
    const netWidth = Math.max(1e-6, maxX - minX);
    const netHeight = Math.max(1e-6, maxY - minY);
    scale = Math.min((width - pad * 2) / netWidth, (height - pad * 2) / netHeight);
    offsetX = (width - netWidth * scale) / 2 - minX * scale;
    offsetY = (height - netHeight * scale) / 2 - minY * scale;
  }

  const rotate = dymaxionSettings ? d3.geoRotation([dymaxionSettings.lon, dymaxionSettings.lat, dymaxionSettings.roll]) : null;

  const output = ctx.getImageData(0, 0, width, height);
  const outData = output.data;

  const insideTri = (p: [number, number], a: [number, number], b: [number, number], c: [number, number]) => {
    const v0 = [c[0] - a[0], c[1] - a[1]];
    const v1 = [b[0] - a[0], b[1] - a[1]];
    const v2 = [p[0] - a[0], p[1] - a[1]];
    const dot00 = v0[0] * v0[0] + v0[1] * v0[1];
    const dot01 = v0[0] * v1[0] + v0[1] * v1[1];
    const dot02 = v0[0] * v2[0] + v0[1] * v2[1];
    const dot11 = v1[0] * v1[0] + v1[1] * v1[1];
    const dot12 = v1[0] * v2[0] + v1[1] * v2[1];
    const invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
    const u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    const v = (dot00 * dot12 - dot01 * dot02) * invDenom;
    return u >= -1e-6 && v >= -1e-6 && u + v <= 1 + 1e-6;
  };

  const applyMatrix = (m: [number, number, number, number, number, number], p: [number, number]) => ([
    m[0] * p[0] + m[1] * p[1] + m[2],
    m[3] * p[0] + m[4] * p[1] + m[5]
  ] as [number, number]);

  const barycentric = (p: [number, number], a: [number, number], b: [number, number], c: [number, number]) => {
    const v0 = [b[0] - a[0], b[1] - a[1]];
    const v1 = [c[0] - a[0], c[1] - a[1]];
    const v2 = [p[0] - a[0], p[1] - a[1]];
    const d00 = v0[0] * v0[0] + v0[1] * v0[1];
    const d01 = v0[0] * v1[0] + v0[1] * v1[1];
    const d11 = v1[0] * v1[0] + v1[1] * v1[1];
    const d20 = v2[0] * v0[0] + v2[1] * v0[1];
    const d21 = v2[0] * v1[0] + v2[1] * v1[1];
    const denom = d00 * d11 - d01 * d01;
    if (!denom) return null;
    const v = (d11 * d20 - d01 * d21) / denom;
    const w = (d00 * d21 - d01 * d20) / denom;
    const u = 1 - v - w;
    return [u, v, w] as [number, number, number];
  };

  const normalizeVec = (v: [number, number, number]) => {
    const len = Math.hypot(v[0], v[1], v[2]) || 1;
    return [v[0] / len, v[1] / len, v[2] / len] as [number, number, number];
  };

  const toLonLat = (v: [number, number, number]) => {
    const lon = Math.atan2(v[2], v[0]) * (180 / Math.PI);
    const lat = Math.asin(Math.max(-1, Math.min(1, v[1]))) * (180 / Math.PI);
    return [lon, lat] as [number, number];
  };

  faces.forEach((face) => {
    const verts = isBlender
      ? face.vertices.map((v) => [v[0] * width, (1 - v[1]) * height]) as [number, number][]
      : face.vertices.map((v) => [v[0] * scale + offsetX, v[1] * scale + offsetY]) as [number, number][];
    const [a, b, c] = verts;
    const minBX = Math.max(0, Math.floor(Math.min(a[0], b[0], c[0])));
    const maxBX = Math.min(width - 1, Math.ceil(Math.max(a[0], b[0], c[0])));
    const minBY = Math.max(0, Math.floor(Math.min(a[1], b[1], c[1])));
    const maxBY = Math.min(height - 1, Math.ceil(Math.max(a[1], b[1], c[1])));

    for (let y = minBY; y <= maxBY; y++) {
      for (let x = minBX; x <= maxBX; x++) {
        const p: [number, number] = [x + 0.5, y + 0.5];
        if (!insideTri(p, a, b, c)) continue;
        const netPoint: [number, number] = isBlender
          ? [p[0] / width, 1 - p[1] / height]
          : [(p[0] - offsetX) / scale, (p[1] - offsetY) / scale];
        const weights = barycentric(netPoint, face.vertices[0], face.vertices[1], face.vertices[2]);
        if (!weights) continue;
        const [u, v, w] = weights;
        const v0 = face.vertices3D[0];
        const v1 = face.vertices3D[1];
        const v2 = face.vertices3D[2];
        const p3 = normalizeVec([
          u * v0[0] + v * v1[0] + w * v2[0],
          u * v0[1] + v * v1[1] + w * v2[1],
          u * v0[2] + v * v1[2] + w * v2[2]
        ]);
        const lonLat = toLonLat(p3);
        const rotated = rotate ? rotate(lonLat) : lonLat;
        const lon = rotated[0];
        const lat = rotated[1];
        const srcX = Math.min(srcWidth - 1, Math.max(0, Math.floor((lon + 180) / 360 * srcWidth)));
        const srcY = Math.min(srcHeight - 1, Math.max(0, Math.floor((90 - lat) / 180 * srcHeight)));
        const srcIdx = (srcY * srcWidth + srcX) * 4;
        const outIdx = (y * width + x) * 4;
        outData[outIdx] = srcData[srcIdx];
        outData[outIdx + 1] = srcData[srcIdx + 1];
        outData[outIdx + 2] = srcData[srcIdx + 2];
        outData[outIdx + 3] = 255;
      }
    }
  });

  ctx.putImageData(output, 0, 0);

  if (DEBUG_DYMAXION) {
    ctx.save();
    ctx.lineWidth = Math.max(1, Math.round(width / 1024));
    ctx.strokeStyle = 'rgba(255,255,255,0.85)';
    ctx.fillStyle = 'rgba(255,255,255,0.9)';
    ctx.font = `${Math.max(12, Math.round(width / 200))}px sans-serif`;
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    faces.forEach((face) => {
      const verts = face.vertices.map((v) => [v[0] * scale + offsetX, v[1] * scale + offsetY]) as [number, number][];
      ctx.beginPath();
      ctx.moveTo(verts[0][0], verts[0][1]);
      ctx.lineTo(verts[1][0], verts[1][1]);
      ctx.lineTo(verts[2][0], verts[2][1]);
      ctx.closePath();
      ctx.stroke();
      const cx = (verts[0][0] + verts[1][0] + verts[2][0]) / 3;
      const cy = (verts[0][1] + verts[1][1] + verts[2][1]) / 3;
      ctx.fillText(String(face.index), cx, cy);
    });
    ctx.restore();
  }

  const link = document.createElement('a');
  const mapName = world.params.mapName || 'map';
  const seed = world.params.seed;
  const layoutSuffix = isBlender ? 'blender' : 'dymaxion';
  link.download = `realmgenesis_${mapName}_${seed}_${viewMode}_${layoutSuffix}_${width}x${height}.png`;
  link.href = canvas.toDataURL('image/png', 0.8);
  link.click();
};

export const exportMap = async (
  world: WorldData, 
  viewMode: ViewMode, 
  resolution: ExportResolution = 4096,
  projectionType: ProjectionType = 'equirectangular',
  dymaxionSettings?: DymaxionExportSettings
) => {
  const width = resolution;
  let height = resolution / 2;
  if (projectionType === 'mercator') height = resolution; 
  if (projectionType === 'orthographic') height = resolution; 
  if (projectionType === 'dymaxion') {
    height = dymaxionSettings?.layout === 'blender' ? resolution : Math.round(resolution * 0.6);
  }

  if (projectionType === 'dymaxion') {
    exportDymaxionRaster(world, viewMode, width, height, dymaxionSettings);
    return;
  }

  const canvas = document.createElement('canvas');
  canvas.width = width;
  canvas.height = height;
  const ctx = canvas.getContext('2d');
  if (!ctx) return;

  // Background color depending on mode
  if (viewMode === 'satellite' || viewMode === 'biome') {
     ctx.fillStyle = '#050505'; // Space/Dark
  } else {
     ctx.fillStyle = '#000000';
  }
  ctx.fillRect(0, 0, width, height);

  ctx.save();
  ctx.translate(width, 0);
  ctx.scale(-1, 1);

  let projection: d3.GeoProjection;
  switch (projectionType) {
      case 'mercator': projection = d3.geoMercator(); break;
      case 'winkeltripel': projection = geoWinkel3(); break;
      case 'robinson': projection = geoRobinson(); break;
      case 'mollweide': projection = geoMollweide(); break;
      case 'orthographic': projection = d3.geoOrthographic(); break;
      case 'dymaxion': projection = createDymaxionProjection(dymaxionSettings?.layout || 'classic'); break;
      case 'equirectangular': default: projection = d3.geoEquirectangular(); break;
  }
  if (projectionType === 'dymaxion' && dymaxionSettings) {
    projection.rotate([dymaxionSettings.lon, dymaxionSettings.lat, dymaxionSettings.roll]);
  }
  projection.fitSize([width, height], { type: "Sphere" } as any);
  const pathGenerator = d3.geoPath(projection, ctx);

  world.cells.forEach((cell, i) => {
    const feature = world.geoJson.features[i];
    if (!feature) return;
    const threeColor = getCellColor(cell, viewMode, world.params.seaLevel);
    ctx.beginPath();
    pathGenerator(feature);
    ctx.fillStyle = '#' + threeColor.getHexString();
    ctx.fill();
  });
  
  ctx.restore();

  const link = document.createElement('a');
  // realmgenesis_mapName_seedValue_viewLayer_projection_resolution.png
  const mapName = world.params.mapName || 'map';
  const seed = world.params.seed;
  link.download = `realmgenesis_${mapName}_${seed}_${viewMode}_${projectionType}_${width}x${height}.png`;
  link.href = canvas.toDataURL('image/png', 0.8); 
  link.click();
};

// --- CONFIG SAVE/LOAD ---

export interface LoadedMap {
    params: WorldParams;
    civData?: CivData;
}

export const saveMapConfig = (params: WorldParams, world?: WorldData) => {
  const date = new Date();
  const dateStr = date.toISOString().split('T')[0].replace(/-/g, ''); // YYYYMMDD
  
  const content = {
      version: "1.4",
      date: date.toISOString(),
      params,
      // We only save the metadata (lore/names). 
      // The geometry (borders/provinces) will be regenerated deterministically from the seed.
      civData: world?.civData || null,
  };

  const dataStr = JSON.stringify(content, null, 2);
  const dataUri = 'data:application/json;charset=utf-8,'+ encodeURIComponent(dataStr);
  
  const filename = `realmgenesis_${params.mapName || 'map'}_${dateStr}_${params.seed}.json`;
  
  const linkElement = document.createElement('a');
  linkElement.setAttribute('href', dataUri);
  linkElement.setAttribute('download', filename);
  linkElement.click();
};

const validateWorldParams = (params: unknown): params is Record<string, unknown> => {
    if (typeof params !== 'object' || params === null || Array.isArray(params)) {
        return false;
    }
    const p = params as Record<string, unknown>;
    const numericBounds: Record<string, [number, number]> = {
        points: [2000, 1000000],
        plates: [2, 50],
        seaLevel: [0.1, 0.9],
        roughness: [0, 1],
        noiseScale: [0.1, 5.0],
        ridgeBlend: [0, 1],
        warpStrength: [0, 2.0],
        plateInfluence: [0, 2.0],
        erosionIterations: [0, 50],
        baseTemperature: [-10, 50],
        poleTemperature: [-50, 20],
        rainfallMultiplier: [0, 3],
        moistureTransport: [0, 1],
        temperatureVariance: [0, 20],
        numFactions: [2, 20],
        capitalSpacing: [0, 1],
        provinceSize: [0.1, 1.0],
        civSizeVariance: [0, 1],
        waterCrossingCost: [0.1, 1.0],
        territorialWaters: [0.01, 1.0],
        axialTilt: [-90, 90],
        cellJitter: [0, 1],
        borderRoughness: [0, 1],
        detailLevel: [0, 10],
        planetRadius: [1000, 20000],
    };
    for (const [key, [min, max]] of Object.entries(numericBounds)) {
        if (key in p) {
            const val = p[key];
            if (typeof val !== 'number' || isNaN(val) || !isFinite(val) || val < min || val > max) {
                return false;
            }
        }
    }
    if ('mapName' in p && typeof p.mapName !== 'string') return false;
    if ('seed' in p && typeof p.seed !== 'string') return false;
    if ('civSeed' in p && typeof p.civSeed !== 'string') return false;
    if ('landStyle' in p && typeof p.landStyle !== 'string') return false;
    if ('maskType' in p && typeof p.maskType !== 'string') return false;
    if ('loreLevel' in p) {
        const ll = p.loreLevel;
        if (typeof ll !== 'number' || ![1, 2, 3].includes(ll)) return false;
    }
    return true;
};

export const loadMapConfig = async (file: File): Promise<LoadedMap | null> => {
    return new Promise((resolve) => {
        const reader = new FileReader();
        reader.onload = (event) => {
            try {
                const json = JSON.parse(event.target?.result as string);
                
                if (json.params) {
                    if (!validateWorldParams(json.params)) {
                        console.error("Invalid or out-of-bounds params in config file");
                        resolve(null);
                        return;
                    }
                    resolve({
                        params: json.params as unknown as WorldParams,
                        civData: json.civData
                    });
                } else if (json.points) {
                    if (!validateWorldParams(json)) {
                        console.error("Invalid or out-of-bounds params in legacy config file");
                        resolve(null);
                        return;
                    }
                    resolve({ params: json as unknown as WorldParams });
                } else {
                    throw new Error("Invalid structure");
                }
            } catch (e) {
                console.error("Failed to parse config", e);
                resolve(null);
            }
        };
        reader.readAsText(file);
    });
};

// --- LOCAL STORAGE MANAGER ---

const LS_KEY = 'realmgenesis_saves';

export interface SavedMapEntry {
    name: string;
    date: number; // timestamp
    params: WorldParams;
    civData?: CivData;
}

export const getSavedMaps = (): SavedMapEntry[] => {
    try {
        const raw = localStorage.getItem(LS_KEY);
        return raw ? JSON.parse(raw) : [];
    } catch(e) {
        return [];
    }
};

export const saveMapToBrowser = (name: string, params: WorldParams, civData?: CivData) => {
    try {
        const current = getSavedMaps();
        const existingIdx = current.findIndex(m => m.name === name);
        const entry: SavedMapEntry = { name, date: Date.now(), params, civData };
        
        if (existingIdx >= 0) {
            current[existingIdx] = entry;
        } else {
            current.push(entry);
        }
        localStorage.setItem(LS_KEY, JSON.stringify(current));
        return true;
    } catch (e) {
        return false;
    }
};

export const deleteSavedMap = (name: string) => {
    const current = getSavedMaps();
    const filtered = current.filter(m => m.name !== name);
    localStorage.setItem(LS_KEY, JSON.stringify(filtered));
};
