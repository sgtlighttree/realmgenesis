export interface Point {
  x: number;
  y: number;
  z: number;
}

export interface Cell {
  id: number;
  center: Point;
  vertices: Point[]; 
  neighbors: number[]; 
  
  // Physical properties
  height: number; // 0-1 normalized
  plateId: number;
  temperature: number; // Celsius
  moisture: number; // 0-1
  biome: BiomeType;
  flux?: number; // Water flux for erosion
  
  // Political/Cultural
  regionId?: number; // Faction ID
  provinceId?: number; // Local ID within faction
  isCapital?: boolean;
  isTown?: boolean;
  population?: number;
}

export enum BiomeType {
  OCEAN = 'Ocean',
  DEEP_OCEAN = 'Deep Ocean',
  
  // E - Polar
  ICE_CAP = 'Ice Cap',
  TUNDRA = 'Tundra',
  
  // B - Dry
  HOT_DESERT = 'Hot Desert',
  COLD_DESERT = 'Cold Desert',
  STEPPE = 'Steppe',
  
  // A - Tropical
  TROPICAL_RAINFOREST = 'Tropical Rainforest',
  TROPICAL_SAVANNA = 'Tropical Savanna',
  
  // C - Temperate
  MEDITERRANEAN = 'Mediterranean',
  TEMPERATE_FOREST = 'Temperate Forest',
  TEMPERATE_RAINFOREST = 'Temperate Rainforest',
  
  // D - Continental
  BOREAL_FOREST = 'Boreal Forest', // Taiga

  // Special
  BEACH = 'Beach',
  VOLCANIC = 'Volcanic',
}

export type LandStyle = 'Continents' | 'Archipelago' | 'Islands' | 'Pangea' | 'Custom';
export type MaskType = 'None' | 'Pangea';

export interface WorldParams {
  // System
  mapName: string;
  points: number;
  seed: string;
  planetRadius: number; // km
  axialTilt: number; // degrees (visual/climate)
  
  // Geography
  landStyle: LandStyle;
  cellJitter: number; // 0-1 randomization of grid
  
  // Advanced Terrain Controls
  noiseScale: number; // Feature Frequency
  ridgeBlend: number; // 0 = Rounded (FBM), 1 = Linear/Spikey (Ridged)
  maskType: MaskType;
  warpStrength: number; // 0-2
  plateInfluence: number; // 0-2
  erosionIterations: number; // 0-50

  plates: number;
  seaLevel: number;
  roughness: number; // 0-1
  detailLevel: number; 
  
  // Climate
  baseTemperature: number; // Equator
  poleTemperature: number; // Pole
  rainfallMultiplier: number;
  moistureTransport: number;
  temperatureVariance: number;
  
  // Political
  numFactions: number;
  civSeed: string; 
  borderRoughness: number; 
  civSizeVariance: number; 
  waterCrossingCost: number; 
  capitalSpacing: number; 
  provinceSize: number; // 0.1 (Small) to 1.0 (Huge)
  
  // Meta
  loreLevel: 1 | 2 | 3;
}

export interface TownData {
    name: string;
    cellId: number;
    population: number;
    isCapital: boolean;
}

export interface ProvinceData {
    id: number; // local id
    name: string;
    towns: TownData[];
    totalPopulation: number;
    color?: string;
}

export interface FactionData {
    id: number;
    name: string;
    color: string;
    capitalId: number;
    provinces: ProvinceData[];
    totalPopulation: number;
    description?: string; // Level 3
}

export interface CivData {
    factions: FactionData[];
}

export interface WorldData {
  cells: Cell[];
  params: WorldParams;
  geoJson: any; // Cached for export
  civData?: CivData;
}

export type ViewMode = 'biome' | 'height' | 'height_bw' | 'temperature' | 'moisture' | 'plates' | 'political' | 'population' | 'province' | 'satellite';

export interface LoreData {
  name: string;
  description: string;
}