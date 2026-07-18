
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
export type DymaxionLayout = 'classic' | 'blender';
export type DymaxionControlMode = 'planet' | 'overlay';

export interface DymaxionSettings {
  layout: DymaxionLayout;
  lon: number;
  lat: number;
  roll: number;
  showOverlay: boolean;
  mode: DymaxionControlMode;
}

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
  mountainHeight: number; // 0.5-2.0 power-curve scale for heights above seaLevel
  oceanDepth: number; // 0.5-2.0 power-curve scale for depths below seaLevel
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
  territorialWaters: number; // Max distance from land to claim water
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

// Minimal GeoJSON shape for the cached d3-geo-voronoi polygons.
// Feature index i corresponds to WorldData.cells[i]. Structurally
// compatible with d3's ExtendedFeature so features can be passed to
// d3.geoPath generators without casts.
export interface GeoJsonFeature {
  type: 'Feature';
  geometry: { type: 'Polygon'; coordinates: number[][][] } | null;
  properties: Record<string, unknown> | null;
}

export interface GeoJsonCollection {
  type: 'FeatureCollection';
  features: GeoJsonFeature[];
}

export interface WorldData {
  cells: Cell[];
  params: WorldParams;
  geoJson: GeoJsonCollection; // Cached for export
  civData?: CivData;
  rivers?: Point[][]; // Array of paths for smooth river rendering
}

export type DisplayMode = 'globe' | 'mercator' | 'dymaxion';
export type InspectMode = 'click' | 'hover' | 'off';
export type ViewMode = 'biome' | 'height' | 'height_bw' | 'temperature' | 'moisture' | 'plates' | 'political' | 'population' | 'province' | 'satellite';
export type EditMode = 'off' | 'terrain-raise' | 'terrain-lower' | 'terrain-flatten' | 'terrain-smooth' | 'biome' | 'political' | 'world-edit';
export type PaintStyle = 'adaptive' | 'freeform';
export const POLITICAL_ERASER_ID = -1;

export interface UndoSnapshot {
    cells: Map<number, { height: number; biome: BiomeType; regionId?: number; provinceId?: number }>;
}

export interface LoreData {
  name: string;
  description: string;
}

export interface LabelVisibility {
  factions: boolean;
  capitals: boolean;
  towns: boolean;
  provinces: boolean;
  borders: boolean;
}

export const DEFAULT_LABEL_VISIBILITY: LabelVisibility = {
  factions: true,
  capitals: true,
  towns: false,
  provinces: false,
  borders: true,
};
