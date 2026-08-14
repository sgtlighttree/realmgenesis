import * as d3 from 'd3';
import { geoWinkel3, geoRobinson, geoMollweide } from 'd3-geo-projection';
import { WorldData, ViewMode, WorldParams, CivData, DymaxionSettings, LabelVisibility, DEFAULT_LABEL_VISIBILITY, MarkerData, MarkerKind, RouteData } from '../types';
import { buildFactionColorMap, buildCultureColorMap, buildReligionColorMap, getCellColor } from './colors';
import { buildDymaxionNet } from './dymaxion';
import { insideTri, barycentric, normalizeVec, toLonLat, projectToDymaxionNet, Point2 } from './geo';
import { collectLabels, drawMapLabels } from './labels';
import { computeShadeMap, computeContourSegments, drawContourPaths } from './shading';
import { computeScaleBar, niceScaleBarLength, drawScaleBar } from './measure';
import { NAME_STYLES, NameStyle } from './namegen';

// Older saved configs predate params added after the config format was
// established (nameStyle, then numCultures for C1); default them here so the
// generator and every UI control bound to these params always receive a
// valid value instead of `undefined` (which would e.g. turn the Cultures
// slider into an uncontrolled input).
const withParamDefaults = (params: WorldParams): WorldParams => ({
  ...params,
  nameStyle: NAME_STYLES.includes(params.nameStyle as NameStyle) ? params.nameStyle : 'fantasy',
  numCultures: typeof params.numCultures === 'number' && isFinite(params.numCultures) ? params.numCultures : 4,
});

// Matches the options offered in the Export tab. 16K+ exceeded browser
// canvas limits on most devices and was removed from the UI.
export type ExportResolution = 2048 | 4096 | 8192;
export type ProjectionType = 'equirectangular' | 'mercator' | 'winkeltripel' | 'orthographic' | 'robinson' | 'mollweide' | 'dymaxion';
export type DymaxionExportSettings = Pick<DymaxionSettings, 'layout' | 'lon' | 'lat' | 'roll'>;

// C3: draw roads (solid) + sea routes (dashed) onto an already-mirrored ctx,
// using the same antimeridian-jump handling as the on-screen Map2D route pass.
const drawRoutesOnCtx = (
  ctx: CanvasRenderingContext2D,
  projection: d3.GeoProjection,
  routes: RouteData[] | undefined,
  lineScale: number,
) => {
  if (!routes) return;
  ctx.save();
  ctx.globalAlpha = 0.9;
  routes.forEach(route => {
    if (route.path.length < 2) return;
    ctx.strokeStyle = route.kind === 'road' ? '#c8a25a' : '#5eb8c8';
    ctx.lineWidth = (route.kind === 'road' ? 1.4 : 1.2) * lineScale;
    ctx.setLineDash(route.kind === 'searoute' ? [5 * lineScale, 4 * lineScale] : []);
    ctx.beginPath();
    let lastLon: number | null = null;
    route.path.forEach((p, i) => {
      const lon = Math.atan2(p.z, p.x) * (180 / Math.PI);
      const lat = Math.asin(Math.max(-1, Math.min(1, p.y))) * (180 / Math.PI);
      const isJump = lastLon !== null && Math.abs(lon - lastLon) > 180;
      const pt = projection([lon, lat]);
      if (pt) {
        if (i === 0 || isJump) ctx.moveTo(pt[0], pt[1]);
        else ctx.lineTo(pt[0], pt[1]);
      }
      lastLon = lon;
    });
    ctx.stroke();
  });
  ctx.setLineDash([]);
  ctx.restore();
};

const renderEquirectangular = (
  world: WorldData,
  viewMode: ViewMode,
  width: number,
  height: number,
  showHillshade = false,
  showContours = false,
  showRoutes = false
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
  const factionColors = buildFactionColorMap(world.civData);
  const cultureColors = buildCultureColorMap(world.cultures);
  const religionColors = buildReligionColorMap(world.religions);
  const shadeMap = showHillshade ? computeShadeMap(world.cells, world.params.seaLevel) : null;

  world.cells.forEach((cell, i) => {
    const feature = world.geoJson.features[i];
    if (!feature) return;
    const threeColor = getCellColor(cell, viewMode, world.params.seaLevel, factionColors, cultureColors, religionColors);
    if (shadeMap) threeColor.multiplyScalar(shadeMap[cell.id]);
    const hexColor = '#' + threeColor.getHexString();
    ctx.beginPath();
    pathGenerator(feature);
    ctx.fillStyle = hexColor;
    ctx.strokeStyle = hexColor;
    ctx.lineWidth = 1;
    ctx.fill();
    ctx.stroke();
  });

  if (showContours) {
    drawContourPaths(ctx, pathGenerator, computeContourSegments(world.cells, world.params.seaLevel, 0.1), Math.max(1, width / 2048));
  }

  if (showRoutes) {
    drawRoutesOnCtx(ctx, projection, world.routes, Math.max(1, width / 2048));
  }

  ctx.restore();
  return canvas;
};

const DEBUG_DYMAXION = false;

const exportDymaxionRaster = (
  world: WorldData,
  viewMode: ViewMode,
  width: number,
  height: number,
  dymaxionSettings?: DymaxionExportSettings,
  labelVisibility: LabelVisibility = DEFAULT_LABEL_VISIBILITY,
  showHillshade = false,
  showContours = false,
  showRoutes = false
) => {
  const srcWidth = width;
  const srcHeight = Math.round(width / 2);
  const source = renderEquirectangular(world, viewMode, srcWidth, srcHeight, showHillshade, showContours, showRoutes);
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

  const labels = collectLabels(world);
  if (labels.length > 0) {
    // Net→canvas mapping must mirror this raster's own fit (blender UV flip /
    // pad-12 classic fit), not getDymaxionNetTransform's pad-8 screen fit.
    const netToCanvas: (p: Point2) => Point2 = isBlender
      ? (p) => [p[0] * width, (1 - p[1]) * height]
      : (p) => [p[0] * scale + offsetX, p[1] * scale + offsetY];
    drawMapLabels(
      ctx,
      labels,
      (position) => {
        const net = projectToDymaxionNet(
          [position.x, position.y, position.z],
          faces,
          dymaxionSettings?.lon ?? 0,
          dymaxionSettings?.lat ?? 0,
          dymaxionSettings?.roll ?? 0,
        );
        return net ? netToCanvas(net) : null;
      },
      2.5,
      labelVisibility,
      width / 1024,
    );
  }

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
  dymaxionSettings?: DymaxionExportSettings,
  labelVisibility: LabelVisibility = DEFAULT_LABEL_VISIBILITY,
  showHillshade = false,
  showContours = false,
  showRoutes = false
) => {
  const width = resolution;
  let height = resolution / 2;
  if (projectionType === 'mercator') height = resolution; 
  if (projectionType === 'orthographic') height = resolution; 
  if (projectionType === 'dymaxion') {
    height = dymaxionSettings?.layout === 'blender' ? resolution : Math.round(resolution * 0.6);
  }

  if (projectionType === 'dymaxion') {
    exportDymaxionRaster(world, viewMode, width, height, dymaxionSettings, labelVisibility, showHillshade, showContours, showRoutes);
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

  // 'dymaxion' is handled by the raster path above and never reaches here
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
  const factionColors = buildFactionColorMap(world.civData);
  const cultureColors = buildCultureColorMap(world.cultures);
  const religionColors = buildReligionColorMap(world.religions);
  const shadeMap = showHillshade ? computeShadeMap(world.cells, world.params.seaLevel) : null;

  world.cells.forEach((cell, i) => {
    const feature = world.geoJson.features[i];
    if (!feature) return;
    const threeColor = getCellColor(cell, viewMode, world.params.seaLevel, factionColors, cultureColors, religionColors);
    if (shadeMap) threeColor.multiplyScalar(shadeMap[cell.id]);
    const hexColor = '#' + threeColor.getHexString();
    ctx.beginPath();
    pathGenerator(feature);
    ctx.fillStyle = hexColor;
    ctx.strokeStyle = hexColor;
    ctx.lineWidth = 1;
    ctx.fill();
    ctx.stroke();
  });

  if (showContours) {
    drawContourPaths(ctx, pathGenerator, computeContourSegments(world.cells, world.params.seaLevel, 0.1), Math.max(1, width / 2048));
  }

  if (showRoutes) {
    drawRoutesOnCtx(ctx, projection, world.routes, Math.max(1, width / 2048));
  }

  ctx.restore();

  const labels = collectLabels(world);
  if (labels.length > 0) {
    drawMapLabels(
      ctx,
      labels,
      (position) => {
        // Orthographic clips paths to the front hemisphere but projection()
        // still maps back-hemisphere points through the disc — skip them.
        if (projectionType === 'orthographic' && position.x <= 0.02) return null;
        const projected = projection(toLonLat([position.x, position.y, position.z]));
        // Cells were drawn under a horizontal mirror; labels draw post-restore.
        return projected ? [width - projected[0], projected[1]] : null;
      },
      2.5,
      labelVisibility,
      width / 1024,
    );
  }

  // Scale bar: bottom-left, post-restore/screen space, for every d3
  // projection reaching this function (dymaxion is the separate raster path
  // above and never gets here). maxPixels scales with resolution so the bar
  // stays proportionally sized across export sizes.
  const scaleBarCenter = projection.invert?.([width / 2, height / 2]);
  if (scaleBarCenter && Number.isFinite(scaleBarCenter[0]) && Number.isFinite(scaleBarCenter[1])) {
    const scaleBarInfo = computeScaleBar(projection, scaleBarCenter, world.params.planetRadius);
    if (scaleBarInfo) {
      const { km, px } = niceScaleBarLength(scaleBarInfo.pixelsPerKm, width * 0.12);
      if (km > 0) {
        const margin = width * 0.02;
        drawScaleBar(ctx, margin, height - margin, km, px);
      }
    }
  }

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
    markers?: MarkerData[];
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
      // Sphere-anchored, so they carry over unchanged across regeneration.
      markers: world?.markers || null,
  };

  const dataStr = JSON.stringify(content, null, 2);
  const dataUri = 'data:application/json;charset=utf-8,'+ encodeURIComponent(dataStr);
  
  const filename = `realmgenesis_${params.mapName || 'map'}_${dateStr}_${params.seed}.json`;
  
  const linkElement = document.createElement('a');
  linkElement.setAttribute('href', dataUri);
  linkElement.setAttribute('download', filename);
  linkElement.click();
};

export const validateWorldParams = (params: unknown): params is Record<string, unknown> => {
    if (typeof params !== 'object' || params === null || Array.isArray(params)) {
        return false;
    }
    const p = params as Record<string, unknown>;
    const numericBounds: Record<string, [number, number]> = {
        // Upper bound matches the UI maximum: generation is single-threaded
        // and larger worlds freeze the tab
        points: [2000, 200000],
        plates: [2, 50],
        seaLevel: [0.1, 0.9],
        roughness: [0, 1],
        noiseScale: [0.1, 5.0],
        ridgeBlend: [0, 1],
        warpStrength: [0, 2.0],
        tectonicStrength: [0, 2.0],
        erosionIterations: [0, 50],
        baseTemperature: [-10, 50],
        poleTemperature: [-50, 20],
        rainfallMultiplier: [0, 3],
        moistureTransport: [0, 1],
        temperatureVariance: [0, 20],
        numFactions: [2, 20],
        numCultures: [2, 8],
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
    if ('nameStyle' in p && typeof p.nameStyle !== 'string') return false;
    if ('loreLevel' in p) {
        const ll = p.loreLevel;
        if (typeof ll !== 'number' || ![1, 2, 3].includes(ll)) return false;
    }
    return true;
};

// Shape-check imported civData so a hand-edited or corrupt file can't crash
// the restore loop or downstream rendering. Params are bounds-checked above;
// civData is metadata (names/colors) so a failed check degrades gracefully
// to "load terrain without metadata" rather than rejecting the file.
export const validateCivData = (civData: unknown): civData is CivData => {
    if (typeof civData !== 'object' || civData === null) return false;
    const c = civData as { factions?: unknown };
    if (!Array.isArray(c.factions)) return false;
    return c.factions.every((f: unknown) => {
        if (typeof f !== 'object' || f === null) return false;
        const fac = f as Record<string, unknown>;
        if (typeof fac.id !== 'number' || typeof fac.name !== 'string') return false;
        if (typeof fac.color !== 'string' || typeof fac.capitalId !== 'number') return false;
        if (!Array.isArray(fac.provinces)) return false;
        return fac.provinces.every((p: unknown) => {
            if (typeof p !== 'object' || p === null) return false;
            const prov = p as Record<string, unknown>;
            return typeof prov.id === 'number' && typeof prov.name === 'string' && Array.isArray(prov.towns);
        });
    });
};

const MARKER_KINDS: MarkerKind[] = ['dungeon', 'ruin', 'battlefield', 'portal', 'poi'];

// Sanitizer, not a type guard: markers are cosmetic overlays, so a malformed
// individual entry is dropped rather than failing the whole load (same
// "degrade gracefully" posture as validateCivData). Returns undefined only
// when the top-level shape isn't an array at all.
export const validateMarkers = (markers: unknown): MarkerData[] | undefined => {
    if (!Array.isArray(markers)) return undefined;
    const result: MarkerData[] = [];
    for (const m of markers) {
        if (typeof m !== 'object' || m === null) continue;
        const rec = m as Record<string, unknown>;
        if (typeof rec.id !== 'number' || !isFinite(rec.id)) continue;
        if (typeof rec.kind !== 'string' || !MARKER_KINDS.includes(rec.kind as MarkerKind)) continue;
        if (typeof rec.name !== 'string') continue;
        if (typeof rec.note !== 'string') continue;
        const pos = rec.position as Record<string, unknown> | undefined;
        if (typeof pos !== 'object' || pos === null) continue;
        const { x, y, z } = pos;
        if (typeof x !== 'number' || !isFinite(x)) continue;
        if (typeof y !== 'number' || !isFinite(y)) continue;
        if (typeof z !== 'number' || !isFinite(z)) continue;
        result.push({ id: rec.id, kind: rec.kind as MarkerKind, name: rec.name, note: rec.note, position: { x, y, z } });
    }
    return result;
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
                    let civData: CivData | undefined;
                    if (json.civData != null) {
                        if (validateCivData(json.civData)) {
                            civData = json.civData;
                        } else {
                            console.error("Ignoring malformed civData in config file; loading terrain only");
                        }
                    }
                    const markers = json.markers != null ? validateMarkers(json.markers) : undefined;
                    resolve({
                        params: withParamDefaults(json.params as unknown as WorldParams),
                        civData,
                        markers,
                    });
                } else if (json.points) {
                    if (!validateWorldParams(json)) {
                        console.error("Invalid or out-of-bounds params in legacy config file");
                        resolve(null);
                        return;
                    }
                    resolve({ params: withParamDefaults(json as unknown as WorldParams) });
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
    markers?: MarkerData[];
}

export const getSavedMaps = (): SavedMapEntry[] => {
    try {
        const raw = localStorage.getItem(LS_KEY);
        return raw ? JSON.parse(raw) : [];
    } catch {
        return [];
    }
};

export const saveMapToBrowser = (name: string, params: WorldParams, civData?: CivData, markers?: MarkerData[]) => {
    try {
        const current = getSavedMaps();
        const existingIdx = current.findIndex(m => m.name === name);
        const entry: SavedMapEntry = { name, date: Date.now(), params, civData, markers };
        
        if (existingIdx >= 0) {
            current[existingIdx] = entry;
        } else {
            current.push(entry);
        }
        localStorage.setItem(LS_KEY, JSON.stringify(current));
        return true;
    } catch {
        return false;
    }
};

export const deleteSavedMap = (name: string) => {
    const current = getSavedMaps();
    const filtered = current.filter(m => m.name !== name);
    localStorage.setItem(LS_KEY, JSON.stringify(filtered));
};
