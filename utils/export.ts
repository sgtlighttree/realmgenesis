import * as d3 from 'd3';
import { geoWinkel3, geoRobinson, geoMollweide } from 'd3-geo-projection';
import { WorldData, ViewMode, WorldParams, CivData } from '../types';
import { getCellColor } from './colors';

export type ExportResolution = 4096 | 8192 | 16384 | 32768;
export type ProjectionType = 'equirectangular' | 'mercator' | 'winkeltripel' | 'orthographic' | 'robinson' | 'mollweide';

export const exportMap = async (
  world: WorldData, 
  viewMode: ViewMode, 
  resolution: ExportResolution = 4096,
  projectionType: ProjectionType = 'equirectangular'
) => {
  const width = resolution;
  let height = resolution / 2;
  if (projectionType === 'mercator') height = resolution; 
  if (projectionType === 'orthographic') height = resolution; 

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
      case 'equirectangular': default: projection = d3.geoEquirectangular(); break;
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
  // terraforge_mapName_seedValue_viewLayer_projection_resolution.png
  const mapName = world.params.mapName || 'map';
  const seed = world.params.seed;
  link.download = `terraforge_${mapName}_${seed}_${viewMode}_${projectionType}_${width}x${height}.png`;
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
  
  const filename = `terraforge_${params.mapName || 'map'}_${dateStr}_${params.seed}.json`;
  
  const linkElement = document.createElement('a');
  linkElement.setAttribute('href', dataUri);
  linkElement.setAttribute('download', filename);
  linkElement.click();
};

export const loadMapConfig = async (file: File): Promise<LoadedMap | null> => {
    return new Promise((resolve) => {
        const reader = new FileReader();
        reader.onload = (event) => {
            try {
                const json = JSON.parse(event.target?.result as string);
                
                // Validate minimal structure
                if (json.params) {
                    resolve({
                        params: json.params,
                        civData: json.civData
                    });
                } else if (json.points) {
                    // Legacy format support
                    resolve({ params: json });
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

const LS_KEY = 'terraforge_saves';

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