import { BiomeType, Cell, CivData, CultureData, GeoJsonCollection, LakeData, GeoFeature, MarkerData, Point, ReligionData, RouteData, WorldData, WorldParams } from '../types';

// Stable biome <-> byte mapping. Order is the declaration order of the enum and
// must never be reordered — it is a wire format, not a display list.
const BIOME_LIST = Object.values(BiomeType) as BiomeType[];
const BIOME_INDEX = new Map<BiomeType, number>(BIOME_LIST.map((b, i) => [b, i]));

// Presence bits. Every optional Cell field needs one: live code tests
// `c.regionId === undefined`, so undefined and a sentinel must round-trip
// distinctly. A shared bitfield beats one array per optional.
const P_FLUX = 1, P_REGION = 2, P_PROVINCE = 4, P_POPULATION = 8,
      P_CULTURE = 16, P_RELIGION = 32, P_CAPITAL = 64, P_TOWN = 128;
// Value bits for the two booleans (presence says "is set", this says "is true").
const V_CAPITAL = 1, V_TOWN = 2;

export interface Ragged { offsets: Uint32Array; data: Float64Array }
export interface RaggedI32 { offsets: Uint32Array; data: Int32Array }

export interface WorldPayload {
  cellCount: number;
  // per-cell scalars
  height: Float64Array;
  temperature: Float64Array;
  moisture: Float64Array;
  flux: Float64Array;
  population: Float64Array;
  plateId: Int32Array;
  regionId: Int32Array;
  provinceId: Int32Array;
  cultureId: Int32Array;
  religionId: Int32Array;
  biome: Uint8Array;
  presence: Uint8Array;
  bools: Uint8Array;
  center: Float64Array;        // 3 per cell
  // ragged per-cell
  vertOffsets: Uint32Array; vertData: Float64Array;   // 3 per vertex
  nbrOffsets: Uint32Array;  nbrData: Int32Array;
  ringOffsets: Uint32Array; ringData: Float64Array;   // 2 per coord (lon, lat)
  geomPresent: Uint8Array;                            // geoJson geometry !== null
  // geoJson feature properties. VERIFIED against d3-geo-voronoi v2.1: every
  // feature carries { site: [lon,lat], sitecoordinates: [lon,lat],
  // neighbours: number[] }. It is per-CELL-scale, not roster-scale, so it gets
  // the SoA treatment like everything else at 200k. `neighbours` is d3's own
  // adjacency and is NOT the same array as Cell.neighbors (which is built from
  // links(), deduped, and differently ordered) — transfer both, never alias.
  propsPresent: Uint8Array;
  geoSite: Float64Array;                              // 2 per cell
  geoSiteCoords: Float64Array;                        // 2 per cell
  geoNbrOffsets: Uint32Array; geoNbrData: Int32Array;
  // rivers: array of paths of Points
  riverOffsets: Uint32Array; riverData: Float64Array; // 3 per point
  hasRivers: boolean;
  // ragged cellIds for lakes and features
  lakeIdOffsets: Uint32Array; lakeIdData: Int32Array;
  featIdOffsets: Uint32Array; featIdData: Int32Array;
  // roster-scale: plain structured clone
  params: WorldParams;
  civData?: CivData;
  cultures?: CultureData[];
  religions?: ReligionData[];
  markers?: MarkerData[];
  routes?: RouteData[];
  lakesMeta?: Omit<LakeData, 'cellIds'>[];
  featuresMeta?: Omit<GeoFeature, 'cellIds'>[];
}

const packRagged = <T>(rows: T[][], stride: number, write: (t: T, out: Float64Array, at: number) => void): Ragged => {
  const offsets = new Uint32Array(rows.length + 1);
  let total = 0;
  for (let i = 0; i < rows.length; i++) { offsets[i] = total; total += rows[i].length; }
  offsets[rows.length] = total;
  const data = new Float64Array(total * stride);
  let at = 0;
  for (const row of rows) for (const item of row) { write(item, data, at); at += stride; }
  return { offsets, data };
};

const packRaggedI32 = (rows: number[][]): RaggedI32 => {
  const offsets = new Uint32Array(rows.length + 1);
  let total = 0;
  for (let i = 0; i < rows.length; i++) { offsets[i] = total; total += rows[i].length; }
  offsets[rows.length] = total;
  const data = new Int32Array(total);
  let at = 0;
  for (const row of rows) for (const n of row) data[at++] = n;
  return { offsets, data };
};

const writePoint = (p: Point, out: Float64Array, at: number) => { out[at] = p.x; out[at + 1] = p.y; out[at + 2] = p.z; };
const readPoint = (d: Float64Array, at: number): Point => ({ x: d[at], y: d[at + 1], z: d[at + 2] });

export const serializeWorld = (world: WorldData): { payload: WorldPayload; transfer: ArrayBuffer[] } => {
  const cells = world.cells;
  const n = cells.length;

  const height = new Float64Array(n), temperature = new Float64Array(n), moisture = new Float64Array(n);
  const flux = new Float64Array(n), population = new Float64Array(n);
  const plateId = new Int32Array(n), regionId = new Int32Array(n), provinceId = new Int32Array(n);
  const cultureId = new Int32Array(n), religionId = new Int32Array(n);
  const biome = new Uint8Array(n), presence = new Uint8Array(n), bools = new Uint8Array(n);
  const center = new Float64Array(n * 3);

  for (let i = 0; i < n; i++) {
    const c = cells[i];
    height[i] = c.height; temperature[i] = c.temperature; moisture[i] = c.moisture;
    plateId[i] = c.plateId;
    biome[i] = BIOME_INDEX.get(c.biome)!;
    writePoint(c.center, center, i * 3);
    let p = 0, b = 0;
    if (c.flux !== undefined) { flux[i] = c.flux; p |= P_FLUX; }
    if (c.regionId !== undefined) { regionId[i] = c.regionId; p |= P_REGION; }
    if (c.provinceId !== undefined) { provinceId[i] = c.provinceId; p |= P_PROVINCE; }
    if (c.population !== undefined) { population[i] = c.population; p |= P_POPULATION; }
    if (c.cultureId !== undefined) { cultureId[i] = c.cultureId; p |= P_CULTURE; }
    if (c.religionId !== undefined) { religionId[i] = c.religionId; p |= P_RELIGION; }
    if (c.isCapital !== undefined) { p |= P_CAPITAL; if (c.isCapital) b |= V_CAPITAL; }
    if (c.isTown !== undefined) { p |= P_TOWN; if (c.isTown) b |= V_TOWN; }
    presence[i] = p; bools[i] = b;
  }

  const verts = packRagged(cells.map(c => c.vertices), 3, writePoint);
  const nbrs = packRaggedI32(cells.map(c => c.neighbors));

  const geomPresent = new Uint8Array(n);
  const rings = world.geoJson.features.map((f, i) => {
    if (!f.geometry) return [];
    geomPresent[i] = 1;
    return f.geometry.coordinates[0];
  });
  const ringPacked = packRagged(rings as number[][][], 2, (coord, out, at) => { out[at] = coord[0]; out[at + 1] = coord[1]; });

  const propsPresent = new Uint8Array(n);
  const geoSite = new Float64Array(n * 2);
  const geoSiteCoords = new Float64Array(n * 2);
  const geoNbrRows: number[][] = world.geoJson.features.map((f, i) => {
    const props = f.properties as { site?: number[]; sitecoordinates?: number[]; neighbours?: number[] } | null;
    if (!props) return [];
    propsPresent[i] = 1;
    if (props.site) { geoSite[i * 2] = props.site[0]; geoSite[i * 2 + 1] = props.site[1]; }
    if (props.sitecoordinates) { geoSiteCoords[i * 2] = props.sitecoordinates[0]; geoSiteCoords[i * 2 + 1] = props.sitecoordinates[1]; }
    return props.neighbours ?? [];
  });
  const geoNbrs = packRaggedI32(geoNbrRows);

  const riverPacked = packRagged(world.rivers ?? [], 3, writePoint);
  const lakeIds = packRaggedI32((world.lakes ?? []).map(l => l.cellIds));
  const featIds = packRaggedI32((world.features ?? []).map(f => f.cellIds));

  const payload: WorldPayload = {
    cellCount: n,
    height, temperature, moisture, flux, population,
    plateId, regionId, provinceId, cultureId, religionId,
    biome, presence, bools, center,
    vertOffsets: verts.offsets, vertData: verts.data,
    nbrOffsets: nbrs.offsets, nbrData: nbrs.data,
    ringOffsets: ringPacked.offsets, ringData: ringPacked.data, geomPresent,
    propsPresent, geoSite, geoSiteCoords,
    geoNbrOffsets: geoNbrs.offsets, geoNbrData: geoNbrs.data,
    riverOffsets: riverPacked.offsets, riverData: riverPacked.data,
    hasRivers: world.rivers !== undefined,
    lakeIdOffsets: lakeIds.offsets, lakeIdData: lakeIds.data,
    featIdOffsets: featIds.offsets, featIdData: featIds.data,
    params: world.params,
    civData: world.civData,
    cultures: world.cultures,
    religions: world.religions,
    markers: world.markers,
    routes: world.routes,
    lakesMeta: world.lakes?.map(({ cellIds: _cellIds, ...rest }) => rest),
    featuresMeta: world.features?.map(({ cellIds: _cellIds, ...rest }) => rest),
  };

  const transfer: ArrayBuffer[] = [];
  for (const v of Object.values(payload as unknown as Record<string, unknown>)) {
    if (ArrayBuffer.isView(v)) transfer.push(v.buffer as ArrayBuffer);
  }
  return { payload, transfer };
};

export const deserializeWorld = (p: WorldPayload): WorldData => {
  const n = p.cellCount;
  const cells: Cell[] = new Array(n);

  for (let i = 0; i < n; i++) {
    const vs = p.vertOffsets[i], ve = p.vertOffsets[i + 1];
    const vertices: Point[] = new Array(ve - vs);
    for (let k = vs; k < ve; k++) vertices[k - vs] = readPoint(p.vertData, k * 3);

    const ns = p.nbrOffsets[i], ne = p.nbrOffsets[i + 1];
    const neighbors: number[] = new Array(ne - ns);
    for (let k = ns; k < ne; k++) neighbors[k - ns] = p.nbrData[k];

    const pr = p.presence[i], bo = p.bools[i];
    const c: Cell = {
      id: i,
      center: readPoint(p.center, i * 3),
      vertices,
      neighbors,
      height: p.height[i],
      plateId: p.plateId[i],
      temperature: p.temperature[i],
      moisture: p.moisture[i],
      biome: BIOME_LIST[p.biome[i]],
    };
    if (pr & P_FLUX) c.flux = p.flux[i];
    if (pr & P_REGION) c.regionId = p.regionId[i];
    if (pr & P_PROVINCE) c.provinceId = p.provinceId[i];
    if (pr & P_POPULATION) c.population = p.population[i];
    if (pr & P_CULTURE) c.cultureId = p.cultureId[i];
    if (pr & P_RELIGION) c.religionId = p.religionId[i];
    if (pr & P_CAPITAL) c.isCapital = (bo & V_CAPITAL) !== 0;
    if (pr & P_TOWN) c.isTown = (bo & V_TOWN) !== 0;
    cells[i] = c;
  }

  const unpackIds = (offsets: Uint32Array, data: Int32Array, i: number): number[] => {
    const s = offsets[i], e = offsets[i + 1];
    const out: number[] = new Array(e - s);
    for (let k = s; k < e; k++) out[k - s] = data[k];
    return out;
  };

  const geoJson: GeoJsonCollection = {
    type: 'FeatureCollection',
    features: new Array(n).fill(null).map((_, i) => {
      let geometry: { type: 'Polygon'; coordinates: number[][][] } | null = null;
      if (p.geomPresent[i]) {
        const s = p.ringOffsets[i], e = p.ringOffsets[i + 1];
        const ring: number[][] = new Array(e - s);
        for (let k = s; k < e; k++) ring[k - s] = [p.ringData[k * 2], p.ringData[k * 2 + 1]];
        geometry = { type: 'Polygon', coordinates: [ring] };
      }
      const properties = p.propsPresent[i]
        ? {
            site: [p.geoSite[i * 2], p.geoSite[i * 2 + 1]],
            sitecoordinates: [p.geoSiteCoords[i * 2], p.geoSiteCoords[i * 2 + 1]],
            neighbours: unpackIds(p.geoNbrOffsets, p.geoNbrData, i),
          }
        : null;
      return { type: 'Feature' as const, geometry, properties };
    }),
  };

  const rivers: Point[][] | undefined = p.hasRivers
    ? new Array(p.riverOffsets.length - 1).fill(null).map((_, i) => {
        const s = p.riverOffsets[i], e = p.riverOffsets[i + 1];
        const path: Point[] = new Array(e - s);
        for (let k = s; k < e; k++) path[k - s] = readPoint(p.riverData, k * 3);
        return path;
      })
    : undefined;

  const lakes: LakeData[] | undefined = p.lakesMeta?.map((m, i) => ({ ...m, cellIds: unpackIds(p.lakeIdOffsets, p.lakeIdData, i) }));
  const features: GeoFeature[] | undefined = p.featuresMeta?.map((m, i) => ({ ...m, cellIds: unpackIds(p.featIdOffsets, p.featIdData, i) }));

  const world: WorldData = { cells, params: p.params, geoJson };
  if (p.civData !== undefined) world.civData = p.civData;
  if (rivers !== undefined) world.rivers = rivers;
  if (lakes !== undefined) world.lakes = lakes;
  if (features !== undefined) world.features = features;
  if (p.cultures !== undefined) world.cultures = p.cultures;
  if (p.religions !== undefined) world.religions = p.religions;
  if (p.markers !== undefined) world.markers = p.markers;
  if (p.routes !== undefined) world.routes = p.routes;
  return world;
};
