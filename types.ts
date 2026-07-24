import { NameStyle } from './utils/namegen';

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
  cultureId?: number; // Culture ID (C1) — undefined on water cells
  religionId?: number; // Religion ID (C2) — undefined on water cells
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

  // Hydrology (derived — assigned after biome classification, never by classifyBiome)
  LAKE = 'Lake',
  SALT_LAKE = 'Salt Lake',
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

  numCultures: number; // 2-8, home count for the culture layer (C1)

  nameStyle: NameStyle;

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

// A culture region (C1): a terrain-affinity-driven population cluster,
// distinct from (and predating) faction borders. Factions inherit their
// naming style from the culture their capital sits in.
export interface CultureData {
  id: number;
  name: string;
  color: string;
  nameStyle: NameStyle;
  homeCellId: number;
}

// The religion layer (C2). 'folk' faiths are one-per-culture and cover
// every land cell by default; 'organized' faiths spread outward from a
// holy city (a town cell) and convert folk cells within a limited budget,
// leaving unreached land as folk. See recalculateReligions in worldGen.ts.
export type ReligionKind = 'folk' | 'organized';

export interface ReligionData {
  id: number;
  name: string;
  kind: ReligionKind;
  color: string;
  cultureId: number | null; // set for folk religions, null for organized
  holyCellId: number | null; // set for organized religions, null for folk
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

// A filled depression surfaced as a first-class water body. Derived from the
// Priority-Flood drainage solve in generateRivers — deterministic from terrain.
export interface LakeData {
  id: number;
  cellIds: number[]; // contiguous flooded land cells forming the lake
  surfaceLevel: number; // filled water level (spill elevation), normalized 0-1
  outflowCellId: number | null; // land cell just outside the sill the lake spills to; null if none
  isEndorheic: boolean; // no surface drainage to the ocean (arid evaporative sink or no spill)
  isSalt: boolean; // hot + arid basin — concentrates salt via evaporation
}

// Auto-detected, named geographic feature (B3). Purely terrain-derived —
// clustered from the cell graph after biomes/lakes, never touched by civ passes.
export type GeoFeatureKind = 'range' | 'desert' | 'forest' | 'sea' | 'ocean' | 'island' | 'lake';

export interface GeoFeature {
  id: number;
  kind: GeoFeatureKind;
  name: string; // filled by the naming pass; seeded from the TERRAIN seed
  cellIds: number[]; // member cells (contiguous over the neighbor graph)
  anchor: Point; // unit-sphere centroid of the member cells (label position)
  size: number; // cellIds.length
}

// User-placed point of interest (C4). Kind is an open-ended flavor tag, not
// tied to generation data.
export type MarkerKind = 'dungeon' | 'ruin' | 'battlefield' | 'portal' | 'poi';

// Anchored to the unit sphere (not a cellId) so it survives regeneration —
// terrain/civ data is rebuilt from the seed, but a sphere position is stable.
// A marker may end up over water after terrain changes; that's acceptable.
export interface MarkerData {
  id: number;
  kind: MarkerKind;
  name: string;
  note: string;
  position: Point;
}

export interface WorldData {
  cells: Cell[];
  params: WorldParams;
  geoJson: GeoJsonCollection; // Cached for export
  civData?: CivData;
  rivers?: Point[][]; // Array of paths for smooth river rendering
  lakes?: LakeData[];
  features?: GeoFeature[]; // named geographic features (B3)
  markers?: MarkerData[]; // user-placed POIs (C4)
  cultures?: CultureData[]; // culture layer (C1)
  religions?: ReligionData[]; // religion layer (C2)
  routes?: RouteData[]; // roads & sea trade routes (C3)
}

// C3: a land road or sea trade route between two towns. Derived from civData +
// terrain, never serialized — regenerates deterministically like rivers.
export interface RouteData {
  path: Point[];
  kind: 'road' | 'searoute';
  fromCellId: number;
  toCellId: number;
}

export type DisplayMode = 'globe' | 'mercator' | 'dymaxion';
export type InspectMode = 'click' | 'hover' | 'off';
export type ViewMode = 'biome' | 'height' | 'height_bw' | 'temperature' | 'moisture' | 'plates' | 'political' | 'population' | 'province' | 'satellite' | 'culture' | 'religion';
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
  geography: boolean; // single toggle for all auto-detected geographic labels (B3)
  markers: boolean; // user-placed POI labels/pins (C4)
}

export const DEFAULT_LABEL_VISIBILITY: LabelVisibility = {
  factions: true,
  capitals: true,
  towns: false,
  provinces: false,
  borders: true,
  geography: true,
  markers: true,
};
