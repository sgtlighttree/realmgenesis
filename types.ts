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

  // V3 terrain model fields (independent crust + Euler-pole tectonics)
  crustType?: number; // V3: 0=oceanic, 1=continental
  crustThickness?: number; // V3: accumulated crustal thickness, normalized 0-1
  upliftAccum?: number; // V3: accumulated kinematic uplift from tectonics
  
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
  tectonicStrength: number; // 0-2, how strongly tectonics deform crust (replaces plateInfluence)
  erosionIterations: number; // 0-50

  // V3 terrain model params
  marginCoupling: number; // 0-1, geometric correlation between mountain belts and continental margins
  numTimesteps: number; // 10-60, simulation timesteps
  simulationResolution: number; // 5000-20000, macro-cell count for tectonic simulation
  plateJitter: number; // 0-3, how irregularly plate seeds are distributed (0=uniform, 1=chaotic, >1=strongly varied sizes)
  boundaryRoughness: number; // 0-3, how jagged/fractal plate boundaries are (0=straight arcs, 1=highly irregular, >1=extreme fracture)

  // D7 part 2 — geophysical realism (all inert at their zero-ish defaults on old seeds)
  spreadRate: number; // 0.004-0.02, chord-per-Ma seafloor spreading rate; smaller = older/deeper floor
  seafloorDepth: number; // 0.3-2.0, linear ocean-floor depth datum (mean water depth); 1.0 = unchanged, <1 shallower seas, >1 deeper abyss. Coastline (seaLevel) held fixed.
  microplateIntensity: number; // 0-1, how many shear-driven microplates to inject (0 = none, byte-identical)
  plateElongation: number; // 0–1: seed-chain length → plate elongation (0 = round blobs)

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
  season: number; // D1: orbital position 0-1, 0.5 = neutral (annual mean). Render-only, not a generation param.
  starClass: StarClass; // D5: host star spectral class — scales global insolation/temperature (G = Sun-like default).
  currentStrength: number; // D2: ocean-current intensity (0-2). 0 = disabled/byte-identical; scales Coriolis, coastal temperature moderation, and warm-current evaporation.

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

export interface OceanCurrentData {
  vx: Float32Array;   // per-cell tangential velocity (unit-sphere frame)
  vy: Float32Array;
  vz: Float32Array;
  sst: Float32Array;  // per-cell SST anomaly, °C (ocean cells; 0 elsewhere)
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
  currents?: OceanCurrentData; // D2 field, gen-time only; never serialized (regenerates like rivers/routes)
}

// C3: a land road or sea trade route between two towns. Derived from civData +
// terrain, never serialized — regenerates deterministically like rivers.
export interface RouteData {
  path: Point[];
  // Cell id per path point, parallel to `path`. Routes are built by walking the
  // cell graph, so every point IS a cell center — keeping the ids lets the F2
  // overlay drape each point at its own terrain radius with no nearest-cell
  // search. Never serialized (routes regenerate deterministically).
  cellIds: number[];
  kind: 'road' | 'searoute';
  fromCellId: number;
  toCellId: number;
}

// V3 tectonic simulation result: per-macro-cell fields projected onto the
// display cells in one nearest-neighbor pass at the end of the simulation.
export interface TectonicResult {
  heights: Float32Array; // per macro-cell, 0-1 normalized
  crustTypes: Uint8Array; // 0=oceanic, 1=continental per macro-cell
  crustThickness: Float32Array; // per macro-cell
  upliftAccum: Float32Array; // per macro-cell, accumulated kinematic uplift
  plateIds: Int32Array; // per macro-cell, assigned plate ID at end of simulation
}

export type DisplayMode = 'globe' | 'mercator' | 'dymaxion';
export type InspectMode = 'click' | 'hover' | 'off';
export type ViewMode = 'biome' | 'height' | 'height_bw' | 'temperature' | 'moisture' | 'plates' | 'political' | 'population' | 'province' | 'satellite' | 'culture' | 'religion';
// D5: host star spectral class (hottest → coolest). Scales global insolation.
export type StarClass = 'O' | 'B' | 'A' | 'F' | 'G' | 'K' | 'M';
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
