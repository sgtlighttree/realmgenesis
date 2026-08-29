import * as d3 from 'd3';
import { geoWinkel3, geoRobinson, geoMollweide } from 'd3-geo-projection';
import { WorldData, ViewMode, Point } from '../types';
import { buildFactionColorMap, buildCultureColorMap, buildReligionColorMap, getCellColor } from './colors';
import { seasonalTemperatureDelta } from './seasons';
import { toLonLat, Point3 } from './geo';
import { computeCoastlineSegments, computeFactionBorderSegments } from './boundaries';
import { computeShadeMap } from './shading';
import { getMapStyle, MapStyleId } from './mapStyle';
import { OverlayInk } from './mapStyle/overlayInk';
import { runStyle } from './mapStyle/passes';
import { placeGlyphs } from './mapStyle/placeGlyphs';
import { SvgSubstrate } from './mapStyle/substrateSvg';
import { collectLabels, LABEL_CONFIG } from './labels';
import { ProjectionType } from './export';
import { elevationMetres } from './datum';

// 'dymaxion' is a raster-only path (triangular net rasterization, not a
// d3 projection) — SVG export reuses the same projection select as PNG
// export but excludes it at the type level.
export type VectorProjectionType = Exclude<ProjectionType, 'dymaxion'>;

const toPoint3 = (p: Point): Point3 => [p.x, p.y, p.z];

const escapeXml = (value: string): string =>
  value
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;')
    .replace(/'/g, '&apos;');

// --- SVG EXPORT ---

const buildProjection = (projectionType: VectorProjectionType, width: number, height: number): d3.GeoProjection => {
  let projection: d3.GeoProjection;
  switch (projectionType) {
    case 'mercator': projection = d3.geoMercator(); break;
    case 'winkeltripel': projection = geoWinkel3(); break;
    case 'robinson': projection = geoRobinson(); break;
    case 'mollweide': projection = geoMollweide(); break;
    case 'orthographic': projection = d3.geoOrthographic(); break;
    case 'equirectangular': default: projection = d3.geoEquirectangular(); break;
  }
  projection.fitSize([width, height], { type: 'Sphere' as const });
  return projection;
};

const renderCellPaths = (world: WorldData, viewMode: ViewMode, pathGenerator: d3.GeoPath): string => {
  const factionColors = buildFactionColorMap(world.civData);
  const cultureColors = buildCultureColorMap(world.cultures);
  const religionColors = buildReligionColorMap(world.religions);
  const seaLevel = world.params.seaLevel;
  const parts: string[] = [];
  world.cells.forEach((cell, i) => {
    const feature = world.geoJson.features[i];
    if (!feature?.geometry) return;
    const d = pathGenerator(feature);
    if (!d) return;
    const hex = '#' + getCellColor(cell, viewMode, {
      seaLevel,
      factionColors,
      cultureColors,
      religionColors,
      seasonalDelta: seasonalTemperatureDelta(cell, world.params),
    }).getHexString();
    parts.push(`<path d="${d}" fill="${hex}" stroke="${hex}" stroke-width="1"/>`);
  });
  return parts.join('');
};

const renderSegmentPathData = (segments: Array<[Point3, Point3]>, pathGenerator: d3.GeoPath): string[] => {
  const paths: string[] = [];
  segments.forEach(([a, b]) => {
    const d = pathGenerator({ type: 'LineString' as const, coordinates: [toLonLat(a), toLonLat(b)] });
    if (d) paths.push(d);
  });
  return paths;
};

const buildCoastlinesGroup = (world: WorldData, pathGenerator: d3.GeoPath, width: number): string => {
  const dList = renderSegmentPathData(computeCoastlineSegments(world), pathGenerator);
  const strokeWidth = Math.max(1, width / 2048);
  const paths = dList.map(d => `<path d="${d}"/>`).join('');
  return `<g id="coastlines" transform="translate(${width},0) scale(-1,1)" fill="none" ` +
    `stroke="rgba(15,23,42,0.9)" stroke-width="${strokeWidth}">${paths}</g>`;
};

const buildBordersGroup = (world: WorldData, pathGenerator: d3.GeoPath, width: number): string => {
  const dList = renderSegmentPathData(computeFactionBorderSegments(world), pathGenerator);
  const strokeScale = width / 1024;
  const underWidth = 2.5 * strokeScale;
  const overWidth = 1 * strokeScale;
  const paths = dList.map(d => `<path d="${d}"/>`).join('');
  return `<g id="borders" transform="translate(${width},0) scale(-1,1)">` +
    `<g fill="none" stroke="rgba(2,6,23,0.9)" stroke-width="${underWidth}" ` +
    `stroke-linejoin="round" stroke-linecap="round">${paths}</g>` +
    `<g fill="none" stroke="rgba(248,250,252,0.95)" stroke-width="${overWidth}" ` +
    `stroke-linejoin="round" stroke-linecap="round">${paths}</g>` +
    `</g>`;
};

// Builds a single "M x,y L x,y M x,y ..." path string, breaking into a new
// subpath whenever consecutive points jump the antimeridian — same detection
// as Map2D's river pass, since a straight line between the projected points
// either side of a jump would otherwise streak across the whole map.
const buildRiverPathData = (path: Point[], projection: d3.GeoProjection): string => {
  const segments: string[] = [];
  let lastLon: number | null = null;
  let penDown = false;
  for (const p of path) {
    const [lon, lat] = toLonLat(toPoint3(p));
    const isJump = lastLon !== null && Math.abs(lon - lastLon) > 180;
    const projected = projection([lon, lat]);
    lastLon = lon;
    if (!projected || !Number.isFinite(projected[0]) || !Number.isFinite(projected[1])) {
      penDown = false;
      continue;
    }
    if (!penDown || isJump) {
      segments.push(`M${projected[0]},${projected[1]}`);
      penDown = true;
    } else {
      segments.push(`L${projected[0]},${projected[1]}`);
    }
  }
  return segments.join(' ');
};

const buildRiversGroup = (world: WorldData, projection: d3.GeoProjection, width: number, ink: OverlayInk): string => {
  const strokeWidth = Math.max(0.5, 1.5 * (width / 2048) * ink.riverWidthScale);
  const paths: string[] = [];
  (world.rivers ?? []).forEach(path => {
    if (path.length < 2) return;
    const d = buildRiverPathData(path, projection);
    if (d) paths.push(`<path d="${d}"/>`);
  });
  return `<g id="rivers" transform="translate(${width},0) scale(-1,1)" fill="none" ` +
    // stroke-opacity, never an rgba() stroke: SVG 1.1 has no rgba() syntax.
    `stroke="${ink.river}" stroke-width="${strokeWidth}" stroke-opacity="0.8">${paths.join('')}</g>`;
};

// C3: roads (solid tan) and sea routes (dashed teal) as two styled sub-groups,
// reusing the same antimeridian-aware polyline builder as rivers.
const buildRoutesGroup = (world: WorldData, projection: d3.GeoProjection, width: number, ink: OverlayInk): string => {
  if (!world.routes || world.routes.length === 0) return '';
  const scale = width / 2048;
  const roadWidth = Math.max(0.5, 1.4 * scale);
  const seaWidth = Math.max(0.5, 1.2 * scale);
  const roadPaths: string[] = [];
  const seaPaths: string[] = [];
  world.routes.forEach(route => {
    if (route.path.length < 2) return;
    const d = buildRiverPathData(route.path, projection);
    if (!d) return;
    (route.kind === 'road' ? roadPaths : seaPaths).push(`<path d="${d}"/>`);
  });
  return `<g id="routes" transform="translate(${width},0) scale(-1,1)" fill="none" stroke-opacity="0.9">` +
    `<g stroke="${ink.road}" stroke-width="${roadWidth}">${roadPaths.join('')}</g>` +
    `<g stroke="${ink.seaRoute}" stroke-width="${seaWidth}" stroke-dasharray="${5 * scale} ${4 * scale}">${seaPaths.join('')}</g>` +
    `</g>`;
};

// Labels are drawn OUTSIDE the mirrored geographic groups: text inside a
// horizontally-flipped <g> would render backwards, so we counter-mirror by
// projecting normally and placing at x = width - projected[0] (matching the
// post-restore canvas label pass in export.ts/Map2D.tsx).
const buildLabelsGroup = (
  world: WorldData,
  projection: d3.GeoProjection,
  projectionType: VectorProjectionType,
  width: number,
): string => {
  const fontScale = width / 1024;
  const parts: string[] = [];
  collectLabels(world).forEach(label => {
    // Orthographic clips paths to the front hemisphere but projection() still
    // maps back-hemisphere points through the disc — skip them (see export.ts).
    if (projectionType === 'orthographic' && label.position.x <= 0.02) return;
    const projected = projection(toLonLat(toPoint3(label.position)));
    if (!projected) return;
    const x = width - projected[0];
    const y = projected[1];
    if (!Number.isFinite(x) || !Number.isFinite(y)) return;

    const config = LABEL_CONFIG[label.kind];
    const fontSize = Math.max(8, config.baseFontSize * fontScale);
    const text = config.uppercase ? label.name.toUpperCase() : label.name;
    const fill = config.fill ?? '#f8fafc';
    const styleAttr = config.italic ? ' font-style="italic"' : '';

    parts.push(
      `<text x="${x.toFixed(2)}" y="${y.toFixed(2)}" font-family="Inter, ui-sans-serif, sans-serif" ` +
      `font-size="${fontSize.toFixed(2)}" font-weight="${config.fontWeight}"${styleAttr} ` +
      `fill="${fill}" fill-opacity="${config.alpha}" stroke="rgba(2,6,23,0.95)" ` +
      `stroke-width="${Math.max(2, fontSize * 0.22).toFixed(2)}" paint-order="stroke" ` +
      `text-anchor="middle" dominant-baseline="middle">${escapeXml(text)}</text>`,
    );
  });
  return `<g id="labels">${parts.join('')}</g>`;
};

// Builds the SVG markup string. Pure — no `document` usage — so it can be
// exercised directly in jsdom-free tests; downloadSVG below handles the DOM.
export const exportSVG = (
  world: WorldData,
  viewMode: ViewMode,
  projectionType: VectorProjectionType = 'equirectangular',
  width = 2048,
  mapStyleId: MapStyleId = 'default',
): string => {
  let height = width / 2;
  if (projectionType === 'mercator' || projectionType === 'orthographic') height = width;

  const projection = buildProjection(projectionType, width, height);
  const pathGenerator = d3.geoPath(projection);

  const background = viewMode === 'satellite' || viewMode === 'biome' ? '#050505' : '#000000';
  const mirror = `translate(${width},0) scale(-1,1)`;

  const labels = buildLabelsGroup(world, projection, projectionType, width);

  const style = getMapStyle(mapStyleId);
  if (style.passes.length > 0) {
    // `mirrored: true` — the body is emitted inside the mirror group below, and
    // the substrate counter-transforms glyphs so they are not drawn backwards.
    // Do NOT hoist glyphs into a separate unmirrored group: that would flip
    // them twice.
    const sub = new SvgSubstrate(
      pathGenerator as unknown as (o: unknown) => string | null, width, height, true,
    );
    const glyphs = style.fillPolicy(viewMode) === 'ramp'
      ? []
      // The OUTPUT width, so a large export is the same map at higher
      // resolution rather than a denser one.
      : placeGlyphs(world.cells, projection, width, {
        seaLevel: world.params.seaLevel,
        seed: world.params.seed,
      });
    runStyle(style, {
      world, viewMode, widthPx: width, heightPx: height,
      glyphs,
      shadeMap: computeShadeMap(world.cells, world.params.seaLevel),
      coastlines: computeCoastlineSegments(world),
      colorCtx: {
        seaLevel: world.params.seaLevel,
        factionColors: buildFactionColorMap(world.civData),
        cultureColors: buildCultureColorMap(world.cultures),
        religionColors: buildReligionColorMap(world.religions),
      },
    }, sub);
    // No background rect: paperPass paints the ground, and a black rect beneath
    // it would show through the grain.
    return (
      `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">` +
      sub.defs() +
      `<g id="style" transform="${mirror}">${sub.body()}</g>` +
      labels +
      `</svg>`
    );
  }

  const cellPaths = renderCellPaths(world, viewMode, pathGenerator);
  const coastlines = buildCoastlinesGroup(world, pathGenerator, width);
  const rivers = buildRiversGroup(world, projection, width, style.overlayInk);
  const routes = buildRoutesGroup(world, projection, width, style.overlayInk);
  const borders = buildBordersGroup(world, pathGenerator, width);

  return (
    `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">` +
    `<rect width="${width}" height="${height}" fill="${background}"/>` +
    `<g id="cells" transform="${mirror}">${cellPaths}</g>` +
    coastlines +
    rivers +
    routes +
    borders +
    labels +
    `</svg>`
  );
};

export const downloadSVG = (
  world: WorldData,
  viewMode: ViewMode,
  projectionType: VectorProjectionType = 'equirectangular',
  width = 2048,
  mapStyleId: MapStyleId = 'default',
): void => {
  const svg = exportSVG(world, viewMode, projectionType, width, mapStyleId);
  const blob = new Blob([svg], { type: 'image/svg+xml' });
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  const mapName = world.params.mapName || 'map';
  const seed = world.params.seed;
  link.download = `realmgenesis_${mapName}_${seed}_${viewMode}_${projectionType}.svg`;
  link.href = url;
  link.click();
  URL.revokeObjectURL(url);
};

// --- GEOJSON EXPORT ---

interface PolygonGeometry { type: 'Polygon'; coordinates: number[][][]; }
interface LineStringGeometry { type: 'LineString'; coordinates: number[][]; }
interface PointGeometry { type: 'Point'; coordinates: number[]; }
type AnyGeometry = PolygonGeometry | LineStringGeometry | PointGeometry;

interface GeoJsonFeatureOut<G extends AnyGeometry> {
  type: 'Feature';
  geometry: G;
  properties: Record<string, unknown>;
}

interface GeoJsonFeatureCollectionOut {
  type: 'FeatureCollection';
  features: GeoJsonFeatureOut<AnyGeometry>[];
}

const buildFactionNameMap = (world: WorldData): Map<number, string> =>
  new Map((world.civData?.factions ?? []).map(f => [f.id, f.name]));

const buildProvinceNameMap = (world: WorldData): Map<string, string> => {
  const map = new Map<string, string>();
  (world.civData?.factions ?? []).forEach(faction => {
    faction.provinces.forEach(province => {
      map.set(`${faction.id}-${province.id}`, province.name);
    });
  });
  return map;
};

const buildCellFeatures = (world: WorldData): GeoJsonFeatureOut<PolygonGeometry>[] => {
  const factionNames = buildFactionNameMap(world);
  const provinceNames = buildProvinceNameMap(world);
  const features: GeoJsonFeatureOut<PolygonGeometry>[] = [];
  // D8a: export real metres relative to sea level, not the raw 0-1 field.
  // Geometry is already geodesic lon/lat; this makes the vertical genuine too,
  // so QGIS/Blender open the file in consistent units.
  const { seaLevel, maxElevationM } = world.params;

  world.cells.forEach((cell, i) => {
    const feature = world.geoJson.features[i];
    if (!feature?.geometry) return;

    const properties: Record<string, unknown> = {
      id: cell.id,
      height: Math.round(elevationMetres(cell.height, seaLevel, maxElevationM)),
      biome: cell.biome,
      temperature: cell.temperature,
      moisture: cell.moisture,
    };
    if (cell.regionId !== undefined) {
      properties.regionId = cell.regionId;
      const factionName = factionNames.get(cell.regionId);
      if (factionName !== undefined) properties.factionName = factionName;
      if (cell.provinceId !== undefined) {
        properties.provinceId = cell.provinceId;
        const provinceName = provinceNames.get(`${cell.regionId}-${cell.provinceId}`);
        if (provinceName !== undefined) properties.provinceName = provinceName;
      }
    }
    if (cell.population !== undefined) properties.population = cell.population;
    if (cell.isCapital !== undefined) properties.isCapital = cell.isCapital;
    if (cell.isTown !== undefined) properties.isTown = cell.isTown;

    features.push({
      type: 'Feature',
      geometry: { type: 'Polygon', coordinates: feature.geometry.coordinates },
      properties,
    });
  });

  return features;
};

const buildRiverFeatures = (world: WorldData): GeoJsonFeatureOut<LineStringGeometry>[] => {
  const features: GeoJsonFeatureOut<LineStringGeometry>[] = [];
  (world.rivers ?? []).forEach((path, index) => {
    if (path.length < 2) return;
    features.push({
      type: 'Feature',
      geometry: { type: 'LineString', coordinates: path.map(p => toLonLat(toPoint3(p))) },
      properties: { kind: 'river', index },
    });
  });
  return features;
};

const buildRouteFeatures = (world: WorldData): GeoJsonFeatureOut<LineStringGeometry>[] => {
  const features: GeoJsonFeatureOut<LineStringGeometry>[] = [];
  (world.routes ?? []).forEach((route, index) => {
    if (route.path.length < 2) return;
    features.push({
      type: 'Feature',
      geometry: { type: 'LineString', coordinates: route.path.map(p => toLonLat(toPoint3(p))) },
      properties: { kind: route.kind, index },
    });
  });
  return features;
};

const buildSegmentFeatures = (
  segments: Array<[Point3, Point3]>,
  kind: 'border' | 'coastline',
): GeoJsonFeatureOut<LineStringGeometry>[] => {
  const features: GeoJsonFeatureOut<LineStringGeometry>[] = [];
  segments.forEach(([a, b]) => {
    features.push({
      type: 'Feature',
      geometry: { type: 'LineString', coordinates: [toLonLat(a), toLonLat(b)] },
      properties: { kind },
    });
  });
  return features;
};

const buildLabelFeatures = (world: WorldData): GeoJsonFeatureOut<PointGeometry>[] => {
  const features: GeoJsonFeatureOut<PointGeometry>[] = [];
  collectLabels(world).forEach(label => {
    features.push({
      type: 'Feature',
      geometry: { type: 'Point', coordinates: toLonLat(toPoint3(label.position)) },
      properties: { kind: label.kind, name: label.name },
    });
  });
  return features;
};

// Builds the GeoJSON as a JSON string (not an object) so this stays a pure
// string-builder alongside exportSVG, callable from tests without a DOM.
// RFC 7946-clean: no top-level "crs" or extra members beyond type/features.
export const exportGeoJSON = (world: WorldData): string => {
  const features: GeoJsonFeatureOut<AnyGeometry>[] = [
    ...buildCellFeatures(world),
    ...buildRiverFeatures(world),
    ...buildRouteFeatures(world),
    ...buildSegmentFeatures(computeFactionBorderSegments(world), 'border'),
    ...buildSegmentFeatures(computeCoastlineSegments(world), 'coastline'),
    ...buildLabelFeatures(world),
  ];
  const collection: GeoJsonFeatureCollectionOut = { type: 'FeatureCollection', features };
  return JSON.stringify(collection);
};

export const downloadGeoJSON = (world: WorldData): void => {
  const json = exportGeoJSON(world);
  const blob = new Blob([json], { type: 'application/geo+json' });
  const url = URL.createObjectURL(blob);
  const link = document.createElement('a');
  const mapName = world.params.mapName || 'map';
  const seed = world.params.seed;
  link.download = `realmgenesis_${mapName}_${seed}.geojson`;
  link.href = url;
  link.click();
  URL.revokeObjectURL(url);
};
