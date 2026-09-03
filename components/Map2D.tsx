import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import * as d3 from 'd3';
import { WorldData, ViewMode, InspectMode, DymaxionSettings, EditMode, LabelVisibility, DEFAULT_LABEL_VISIBILITY, Point, MarkerData } from '../types';
import { getCellColor } from '../utils/colors';
import { computeCoastlineSegments, computeFactionBorderSegments } from '../utils/boundaries';
import { buildMapGeometryCache } from '../utils/mapCache';
import { buildCellColorCache } from '../utils/mapColorCache';
import { buildCellQuadtree, findCellIdAtPoint } from '../utils/mapPick';
import { getMapStyle, MapStyleId } from '../utils/mapStyle';
import { OverlayInk } from '../utils/mapStyle/overlayInk';
import { buildProjection, FlatProjectionType } from '../utils/projections';
import { DISPLAY_MODES } from './ViewControls';
import { useLabelFonts } from '../utils/mapStyle/useLabelFonts';
import { useDeskTexture } from '../utils/mapStyle/useDeskTexture';
import { paintDesk } from '../utils/mapStyle/desk';
import { runStyle } from '../utils/mapStyle/passes';
import { placeGlyphs } from '../utils/mapStyle/placeGlyphs';
import { Canvas2DSubstrate } from '../utils/mapStyle/substrateCanvas';
import { seasonalTemperatureDelta } from '../utils/seasons';
import { insideTri, barycentric, normalizeVec, toLonLat, lonLatToPoint3, getDymaxionNetTransform, projectDymaxionPoint, Point2, Point3 } from '../utils/geo';
import { collectLabels, drawMapLabels } from '../utils/labels';
import { computeShadeMap, computeContourSegments, drawContourPaths, contourInterval } from '../utils/shading';
import { drawCurrents2D } from './overlays/currents2D';
import { computeScaleBar, computeInterruptedScaleBar, niceScaleBarLength } from '../utils/measure';

type Size = { width: number; height: number };

const clamp = (v: number, min: number, max: number) => Math.max(min, Math.min(max, v));
const INTERACTION_DPR = 1;
const MAX_SHARP_DPR = 3;
const MAX_SHARP_SCALE = 2.5;
const SETTLE_MS = 200;

/**
 * Rivers, in projected space.
 *
 * ONE routine, called by both the Dymaxion source raster and the direct
 * projection path. It used to be two byte-identical blocks differing only in
 * their context variable and DPR divisor, each with its own hardcoded
 * `#38bdf8` — two of the five sites that made the river colour unreachable by
 * a map style. Collapsed here so the next style cannot diverge the same way.
 */
const drawRiverPaths = (
  ctx: CanvasRenderingContext2D,
  rivers: Point[][],
  project: (p: [number, number]) => [number, number] | null,
  lineWidth: number,
  ink: OverlayInk,
) => {
  ctx.strokeStyle = ink.river;
  ctx.lineWidth = lineWidth;
  ctx.globalAlpha = 0.8;
  for (const path of rivers) {
    if (path.length < 2) continue;
    ctx.beginPath();
    let lastLon: number | null = null;
    path.forEach((p, i) => {
      const lon = Math.atan2(p.z, p.x) * (180 / Math.PI);
      const lat = Math.asin(Math.max(-1, Math.min(1, p.y))) * (180 / Math.PI);
      // Antimeridian crossing: break the subpath rather than draw across the map.
      const isJump = lastLon !== null && Math.abs(lon - lastLon) > 180;
      const pt = project([lon, lat]);
      if (pt) {
        if (i === 0 || isJump) ctx.moveTo(pt[0], pt[1]);
        else ctx.lineTo(pt[0], pt[1]);
      }
      lastLon = lon;
    });
    ctx.stroke();
  }
  ctx.globalAlpha = 1;
};

const drawFactionBorders = (ctx: CanvasRenderingContext2D, pathGenerator: d3.GeoPath, borders: Array<[Point3, Point3]>, lineWidth: number, ink: OverlayInk) => {
  if (borders.length === 0) return;

  const drawPass = () => {
    borders.forEach(([a, b]) => {
      ctx.beginPath();
      pathGenerator({
        type: 'LineString',
        coordinates: [toLonLat(a), toLonLat(b)],
      } as any);
      ctx.stroke();
    });
  };

  ctx.save();
  ctx.lineJoin = 'round';
  ctx.lineCap = 'round';
  ctx.globalAlpha = 0.95;
  ctx.strokeStyle = ink.borderCasing;
  ctx.lineWidth = lineWidth * ink.borderWidthScale * 2.5;
  ctx.setLineDash([]);
  drawPass();
  ctx.strokeStyle = ink.border;
  ctx.lineWidth = lineWidth * ink.borderWidthScale;
  ctx.setLineDash(ink.borderDash);
  drawPass();
  ctx.setLineDash([]);
  ctx.restore();
};

// Small amber diamond per marker, drawn in the same projected space as the
// adjacent drawMapLabels call so both share one projection lambda. halfSize
// is in the caller's coordinate space (CSS px for the transformed mercator
// ctx, device px for the raw dymaxion raster) — see call sites.
const drawMarkerPins = (
  ctx: CanvasRenderingContext2D,
  markers: MarkerData[],
  project: (position: { x: number; y: number; z: number }) => [number, number] | null,
  halfSize: number,
): void => {
  if (markers.length === 0) return;
  ctx.save();
  ctx.fillStyle = '#f59e0b';
  ctx.strokeStyle = 'rgba(28, 18, 7, 0.9)';
  ctx.lineWidth = Math.max(0.75, halfSize * 0.3);
  for (const marker of markers) {
    const projected = project(marker.position);
    if (!projected) continue;
    const [x, y] = projected;
    ctx.beginPath();
    ctx.moveTo(x, y - halfSize);
    ctx.lineTo(x + halfSize, y);
    ctx.lineTo(x, y + halfSize);
    ctx.lineTo(x - halfSize, y);
    ctx.closePath();
    ctx.fill();
    ctx.stroke();
  }
  ctx.restore();
};

export const rasterizeDymaxionSource = ({
  sourceData,
  sourceWidth,
  sourceHeight,
  outputData,
  outputWidth,
  outputHeight,
  canvasWidth,
  canvasHeight,
  renderDpr,
  layout,
  lon,
  lat,
  roll,
}: {
  sourceData: Uint8ClampedArray;
  sourceWidth: number;
  sourceHeight: number;
  outputData: Uint8ClampedArray;
  outputWidth: number;
  outputHeight: number;
  canvasWidth: number;
  canvasHeight: number;
  renderDpr: number;
  layout: DymaxionSettings['layout'];
  lon: number;
  lat: number;
  roll: number;
}) => {
  const { faces, scale, offsetX, offsetY } = getDymaxionNetTransform(layout, canvasWidth, canvasHeight);
  const rotate = d3.geoRotation([lon, lat, roll]);

  faces.forEach((face) => {
    const verts = face.vertices.map((v) => [v[0] * scale + offsetX, v[1] * scale + offsetY]) as Point2[];
    const [a, b, c] = verts;
    const minBX = Math.max(0, Math.floor(Math.min(a[0], b[0], c[0])));
    const maxBX = Math.min(canvasWidth - 1, Math.ceil(Math.max(a[0], b[0], c[0])));
    const minBY = Math.max(0, Math.floor(Math.min(a[1], b[1], c[1])));
    const maxBY = Math.min(canvasHeight - 1, Math.ceil(Math.max(a[1], b[1], c[1])));

    const startOY = Math.floor(minBY * renderDpr);
    const endOY = Math.min(outputHeight - 1, Math.ceil(maxBY * renderDpr));
    const startOX = Math.floor(minBX * renderDpr);
    const endOX = Math.min(outputWidth - 1, Math.ceil(maxBX * renderDpr));

    for (let oy = startOY; oy <= endOY; oy++) {
      for (let ox = startOX; ox <= endOX; ox++) {
        const x = ox / renderDpr;
        const y = oy / renderDpr;
        if (!insideTri([x, y], a, b, c)) continue;

        const netPoint: Point2 = [(x - offsetX) / scale, (y - offsetY) / scale];
        const weights = barycentric(netPoint, face.vertices[0], face.vertices[1], face.vertices[2]);
        if (!weights) continue;
        const [u, v, w] = weights;
        const v0 = face.vertices3D[0];
        const v1 = face.vertices3D[1];
        const v2 = face.vertices3D[2];
        const p3 = normalizeVec([
          u * v0[0] + v * v1[0] + w * v2[0],
          u * v0[1] + v * v1[1] + w * v2[1],
          u * v0[2] + v * v1[2] + w * v2[2],
        ]);
        const [rotatedLon, rotatedLat] = rotate(toLonLat(p3));
        const srcX = Math.min(sourceWidth - 1, Math.max(0, Math.floor((rotatedLon + 180) / 360 * sourceWidth)));
        const srcY = Math.min(sourceHeight - 1, Math.max(0, Math.floor((90 - rotatedLat) / 180 * sourceHeight)));
        const srcIdx = (srcY * sourceWidth + srcX) * 4;
        const outIdx = (oy * outputWidth + ox) * 4;
        if (outIdx >= 0 && outIdx < outputData.length - 3) {
          outputData[outIdx] = sourceData[srcIdx];
          outputData[outIdx + 1] = sourceData[srcIdx + 1];
          outputData[outIdx + 2] = sourceData[srcIdx + 2];
          outputData[outIdx + 3] = 255;
        }
      }
    }
  });
};

const Map2D: React.FC<{
  world: WorldData | null;
  viewMode: ViewMode;
  inspectMode: InspectMode;
  onInspect: (cellId: number | null) => void;
  highlightCellId?: number | null;
  /**
   * Which flat view to draw. Everything except `dymaxion` is a d3 projection
   * from `utils/projections.ts`; Dymaxion is an unfolded net with its own path.
   */
  projectionType?: FlatProjectionType | 'dymaxion';
  dymaxionSettings?: DymaxionSettings;
  showGrid?: boolean;
  showRivers?: boolean;
  showRoutes?: boolean;
  showHillshade?: boolean;
  mapStyleId?: MapStyleId;
  showContours?: boolean;
  showCurrents?: boolean;
  labelVisibility?: LabelVisibility;
  editMode?: EditMode;
  onPaint?: (cellId: number, phase: 'start' | 'stroke' | 'end', isRightClick?: boolean) => void;
  factionColors?: Map<number, string>;
  cultureColors?: Map<number, string>;
  religionColors?: Map<number, string>;
  brushSize?: number;
  rulerArc?: Point[] | null;
}> = ({ world, viewMode, inspectMode, onInspect, highlightCellId = null, projectionType = 'mercator', dymaxionSettings, showGrid = false, showRivers = true, showRoutes = false, showHillshade = false, mapStyleId = 'default', showContours = false, showCurrents = false, labelVisibility = DEFAULT_LABEL_VISIBILITY, editMode = 'off', onPaint, factionColors, cultureColors, religionColors, brushSize = 1, rulerArc = null }) => {
  const containerRef = useRef<HTMLDivElement>(null);
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const offscreenRef = useRef<HTMLCanvasElement | null>(null);
  const pickCanvasRef = useRef<HTMLCanvasElement | null>(null);
  const pickCtxRef = useRef<CanvasRenderingContext2D | null>(null);
  /**
   * Bumped whenever the Dymaxion pick buffer is rebuilt.
   *
   * The buffer lives in a REF, filled by an effect — invisible to anything that
   * reads it during render. The scale bar does exactly that, and without this it
   * computed once, before the effect had run, saw a null ref, and returned null
   * for good: nothing in its dependency array ever changed again. The bar simply
   * never appeared, with no error anywhere.
   *
   * Same lesson as `useLabelFonts`: a value produced after the first render has
   * to change some IDENTITY, or its readers never recompute.
   */
  const [pickVersion, setPickVersion] = useState(0);
  const [size, setSize] = useState<Size>({ width: 0, height: 0 });
  const [scale, setScale] = useState(1);
  const [offset, setOffset] = useState({ x: 0, y: 0 });
  const dragging = useRef(false);
  const dragDistance = useRef(0);
  const lastPos = useRef({ x: 0, y: 0 });
  const isPaintingRef = useRef(false);
  const lastPaintCell = useRef<number | null>(null);
  const [isSpaceHeld, setIsSpaceHeld] = useState(false);
  const [mousePos, setMousePos] = useState<{ x: number; y: number } | null>(null);
  const dymaxionLayout = dymaxionSettings?.layout || 'blender';
  const dymaxionLon = dymaxionSettings?.lon ?? 0;
  const dymaxionLat = dymaxionSettings?.lat ?? 0;
  const dymaxionRoll = dymaxionSettings?.roll ?? 0;
  // Sourced from the shared boundary scan (utils/boundaries.ts), matching how
  // coastlines are sourced — was a byte-for-byte private copy. The
  // labelVisibility.borders gate is preserved here (computeFactionBorderSegments
  // has no visibility argument): dropping it would draw borders while toggled off.
  const factionBorders = useMemo(
    () => (world && labelVisibility.borders ? computeFactionBorderSegments(world) : []),
    [world, labelVisibility.borders],
  );
  const mapLabels = useMemo(
    () => (world ? collectLabels(world) : []),
    [world],
  );
  // Keyed on world identity: heights mutate in place on paint strokes and
  // WorldData is shallow-copied, matching the other paint-sensitive memos.
  const shadeMap = useMemo(
    () => (world && showHillshade ? computeShadeMap(world.cells, world.params.seaLevel) : null),
    [world, showHillshade],
  );
  const contourSegments = useMemo(
    () => (world && showContours ? computeContourSegments(world.cells, world.params.seaLevel, contourInterval(world.cells, world.params.seaLevel)) : []),
    [world, showContours],
  );

  useEffect(() => {
    const onDown = (e: KeyboardEvent) => { if (e.code === 'Space' && !e.repeat) { e.preventDefault(); setIsSpaceHeld(true); } };
    const onUp = (e: KeyboardEvent) => { if (e.code === 'Space') setIsSpaceHeld(false); };
    window.addEventListener('keydown', onDown);
    window.addEventListener('keyup', onUp);
    return () => { window.removeEventListener('keydown', onDown); window.removeEventListener('keyup', onUp); };
  }, []);
  const [isInteracting, setIsInteracting] = useState(false);
  const [qualityDpr, setQualityDpr] = useState(1);
  const rafId = useRef<number | null>(null);
  const pendingOffset = useRef<{ x: number; y: number } | null>(null);
  const settleTimer = useRef<number | null>(null);
  const wheelTimer = useRef<number | null>(null);
  const [renderCount, setRenderCount] = useState(0);
  const inspectEnabled = inspectMode === 'click';

  useEffect(() => {
    if (!containerRef.current) return;
    const ro = new ResizeObserver(entries => {
      const entry = entries[0];
      setSize({
        width: Math.floor(entry.contentRect.width),
        height: Math.floor(entry.contentRect.height)
      });
    });
    ro.observe(containerRef.current);
    return () => { ro.disconnect(); };
  }, []);

  useEffect(() => {
    if (!size.width || !size.height) return;
    setScale(1);
    setOffset({ x: 0, y: 0 });
  }, [size.width, size.height, world?.params.seed]);

  useEffect(() => {
    if (settleTimer.current) {
      window.clearTimeout(settleTimer.current);
      settleTimer.current = null;
    }
    if (isInteracting) {
      setQualityDpr(INTERACTION_DPR);
      return;
    }
    settleTimer.current = window.setTimeout(() => {
      const baseDpr = Math.min(2, window.devicePixelRatio || 1);
      const scaled = baseDpr * Math.min(scale, MAX_SHARP_SCALE);
      const target = Math.min(MAX_SHARP_DPR, scaled);
      setQualityDpr(target);
      settleTimer.current = null;
    }, SETTLE_MS);
  }, [isInteracting, scale]);

  const scaleRef = useRef(scale);
  const offsetRef = useRef(offset);

  useEffect(() => {
    scaleRef.current = scale;
  }, [scale]);

  useEffect(() => {
    offsetRef.current = offset;
  }, [offset]);

  const projection = useMemo(() => {
    if (!size.width || !size.height) return null;
    if (projectionType === 'dymaxion') return null;
    // Same factory the exports use, so screen and file cannot drift apart.
    return buildProjection(projectionType, size.width, size.height);
  }, [size.width, size.height, projectionType]);

  // F3 Phase 1 caches. Built ONCE per world/projection — NEVER keyed on
  // scale/offset/qualityDpr — so pan/zoom settle redraws reuse them instead of
  // reprojecting + recolouring every cell (the ~2.2s stall). Non-dymaxion only:
  // the Dymaxion path renders a source raster with its own per-feature loop and
  // pick buffer, and must not receive these (see invariant 1). Built from the
  // SAME memoized `projection` + `size` the rest of the component uses, so cached
  // geometry aligns with everything drawn through that projection.
  const geometryCache = useMemo(
    () => (world && projection && projectionType !== 'dymaxion'
      ? buildMapGeometryCache(world, projection, size.width, size.height)
      : null),
    [world, projection, projectionType, size.width, size.height],
  );
  const colorCache = useMemo(
    // shadeMap is intentionally NOT baked in (passed null). The ONLY consumer is
    // landPass's fast path, whose per-cell fallback does not apply shadeMap
    // either — hillshadePass overlays relief shading separately. Baking it here
    // would double-shade ramp modes with hillshade on. shadeMap stays in the dep
    // key (per invariant 4) so a future shaded variant would rebuild correctly.
    () => (world ? buildCellColorCache(world, viewMode, {
      seaLevel: world.params.seaLevel, factionColors, cultureColors, religionColors,
    }, null) : null),
    // shadeMap/season/seaLevel are kept in the key deliberately (invariant 4):
    // shadeMap is not read here today (passed null) but must trigger a rebuild if
    // a future shaded variant bakes it in; season/seaLevel change the seasonal
    // delta and sea/land split the cache depends on, even though `world` usually
    // co-changes. The lint rule sees these as redundant; they are intentional.
    // eslint-disable-next-line react-hooks/exhaustive-deps
    [world, viewMode, factionColors, cultureColors, religionColors, shadeMap,
     world?.params.season, world?.params.seaLevel],
  );
  const pickQuadtree = useMemo(
    () => (world && projection && projectionType !== 'dymaxion'
      ? buildCellQuadtree(world, projection, size.width, size.height) : null),
    [world, projection, projectionType, size.width, size.height],
  );

  const style = useMemo(() => getMapStyle(mapStyleId), [mapStyleId]);
  // Identity changes once the style's webfonts resolve, which repaints the
  // canvas through the ordinary dependency rules. See useLabelFonts.
  const labelTheme = useLabelFonts(style.labelTheme);
  const overlayInk = style.overlayInk;
  // null until the image decodes, so the desk falls back to its colour rather
  // than painting nothing. See useDeskTexture.
  const deskTexture = useDeskTexture(style.palette.deskTexture);

  // A style with no passes draws nothing, so the legacy per-cell loop runs
  // instead. This is the ONE test for that — never a comparison against the
  // style id, which would put the same invariant in two places.
  const styled = style.passes.length > 0;

  // Coastline geometry is only needed by a style, and the scan is O(cells x
  // neighbours), so it stays behind that flag.
  const coastlines = useMemo(
    () => (world && styled ? computeCoastlineSegments(world) : []),
    [world, styled],
  );

  // Placement depends on the projection and the output width, so it re-runs
  // when either changes — but not per frame. Ramp modes get an empty list:
  // glyphs would fight a continuous fill.
  const glyphs = useMemo(() => {
    if (!world || !styled || !projection) return [];
    if (style.fillPolicy(viewMode) === 'ramp') return [];
    return placeGlyphs(world.cells, projection, size.width, {
      seaLevel: world.params.seaLevel,
      seed: world.params.seed,
    });
  }, [world, styled, style, viewMode, projection, size.width]);

  // The Dymaxion path rasterizes an equirectangular SOURCE buffer at its own
  // width, then re-projects it onto triangles. Glyphs must be placed against
  // THAT projection and THAT width, not the screen's Mercator — otherwise they
  // land in the wrong place and at the wrong size.
  const dymaxionGlyphs = useMemo(() => {
    if (!world || !styled || projectionType !== 'dymaxion') return [];
    if (style.fillPolicy(viewMode) === 'ramp') return [];
    const srcWidth = Math.max(1, Math.floor(size.width * qualityDpr));
    const srcHeight = Math.max(1, Math.round(srcWidth / 2));
    const srcProjection = d3.geoEquirectangular()
      .fitSize([srcWidth, srcHeight], { type: 'Sphere' } as d3.GeoPermissibleObjects);
    return placeGlyphs(world.cells, srcProjection, srcWidth, {
      seaLevel: world.params.seaLevel,
      seed: world.params.seed,
    });
  }, [world, styled, style, viewMode, projectionType, size.width, qualityDpr]);

  useEffect(() => {
    if (!world || !size.width || !size.height) return;
    const offscreen = offscreenRef.current ?? document.createElement('canvas');
    offscreenRef.current = offscreen;

    const renderDpr = qualityDpr;
    const ctx = offscreen.getContext('2d');
    if (!ctx) return;

    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, offscreen.width, offscreen.height);
    ctx.setTransform(renderDpr, 0, 0, renderDpr, 0, 0);
    ctx.translate(size.width, 0);
    ctx.scale(-1, 1);

    if (projectionType === 'dymaxion') {
      const srcWidth = Math.max(1, Math.floor(size.width * renderDpr));
      const srcHeight = Math.max(1, Math.round(srcWidth / 2));
      const source = document.createElement('canvas');
      source.width = srcWidth;
      source.height = srcHeight;
      const srcCtx = source.getContext('2d');
      if (!srcCtx) return;
      const projection = d3.geoEquirectangular().fitSize([srcWidth, srcHeight], { type: 'Sphere' } as d3.GeoPermissibleObjects);
      const pathGenerator = d3.geoPath(projection, srcCtx);
      srcCtx.fillStyle = viewMode === 'satellite' || viewMode === 'biome' ? '#050505' : '#000000';
      srcCtx.fillRect(0, 0, srcWidth, srcHeight);
      srcCtx.save();
      srcCtx.translate(srcWidth, 0);
      srcCtx.scale(-1, 1);
      if (styled) {
        const srcSub = new Canvas2DSubstrate(
          srcCtx, pathGenerator as unknown as (o: unknown) => unknown, srcWidth, srcHeight, true,
        );
        runStyle(style, {
          world, viewMode, widthPx: srcWidth, heightPx: srcHeight,
          glyphs: dymaxionGlyphs, shadeMap, coastlines,
          colorCtx: {
            seaLevel: world.params.seaLevel, factionColors, cultureColors, religionColors,
          },
        }, srcSub);
      } else {
      for (let i = 0; i < world.cells.length; i++) {
        const feature = world.geoJson?.features?.[i];
        if (!feature || !feature.geometry) continue;
        const color = getCellColor(world.cells[i], viewMode, {
          seaLevel: world.params.seaLevel,
          factionColors,
          cultureColors,
          religionColors,
          seasonalDelta: seasonalTemperatureDelta(world.cells[i], world.params),
        });
        if (shadeMap) color.multiplyScalar(shadeMap[i]);
        const hexColor = '#' + color.getHexString();
        srcCtx.beginPath();
        pathGenerator(feature);
        srcCtx.fillStyle = hexColor;
        srcCtx.strokeStyle = hexColor;
        srcCtx.lineWidth = 1;
        srcCtx.fill();
        srcCtx.stroke();
      }
      }

      drawContourPaths(srcCtx, pathGenerator, contourSegments, Math.max(1, renderDpr));

      if (highlightCellId !== null) {
        const feature = world.geoJson?.features?.[highlightCellId];
        if (feature && feature.geometry) {
          srcCtx.save();
          srcCtx.strokeStyle = '#f9a8a8';
          srcCtx.lineWidth = Math.max(2, 3 * renderDpr);
          srcCtx.lineJoin = 'round';
          srcCtx.lineCap = 'round';
          srcCtx.beginPath();
          pathGenerator(feature);
          srcCtx.stroke();
          srcCtx.restore();
        }
      }

      // Draw Grid on source equirectangular canvas
      if (showGrid) {
        srcCtx.strokeStyle = 'rgba(255,255,255,0.15)';
        srcCtx.lineWidth = 1;
        srcCtx.beginPath();
        pathGenerator(d3.geoGraticule().step([10, 10])());
        srcCtx.stroke();
      }

      // Draw Rivers on source equirectangular canvas
      if (showRivers && world.rivers) {
        drawRiverPaths(srcCtx, world.rivers, projection, Math.max(0.5, 1.5 / renderDpr), overlayInk);
      }

      // Draw Routes on source equirectangular canvas (C3)
      if (showRoutes && world.routes) {
        srcCtx.globalAlpha = 0.9;
        const rw = Math.max(0.5, 1.4 / renderDpr);
        world.routes.forEach(route => {
          if (route.path.length < 2) return;
          srcCtx.strokeStyle = route.kind === 'road' ? '#c8a25a' : '#5eb8c8';
          srcCtx.lineWidth = rw;
          srcCtx.setLineDash(route.kind === 'searoute' ? [rw * 4, rw * 3] : []);
          srcCtx.beginPath();
          let lastLon: number | null = null;
          route.path.forEach((p, i) => {
            const lon = Math.atan2(p.z, p.x) * (180 / Math.PI);
            const lat = Math.asin(Math.max(-1, Math.min(1, p.y))) * (180 / Math.PI);
            const isJump = lastLon !== null && Math.abs(lon - lastLon) > 180;
            const pt = projection([lon, lat]);
            if (pt) {
              if (i === 0 || isJump) srcCtx.moveTo(pt[0], pt[1]);
              else srcCtx.lineTo(pt[0], pt[1]);
            }
            lastLon = lon;
          });
          srcCtx.stroke();
        });
        srcCtx.setLineDash([]);
        srcCtx.globalAlpha = 1.0;
      }

      drawFactionBorders(
        srcCtx,
        pathGenerator,
        factionBorders,
        Math.max(1.5, 2 * renderDpr),
        overlayInk,
      );

      srcCtx.restore();

      const srcImage = srcCtx.getImageData(0, 0, srcWidth, srcHeight);
      const srcData = srcImage.data;

      const canvasWidth = size.width;
      const canvasHeight = size.height;
      const outWidth = Math.floor(canvasWidth * renderDpr);
      const outHeight = Math.floor(canvasHeight * renderDpr);

      // Create a temporary canvas for heavy rendering to avoid clearing the screen too early
      const tempCanvas = document.createElement('canvas');
      tempCanvas.width = outWidth;
      tempCanvas.height = outHeight;
      const tCtx = tempCanvas.getContext('2d');
      if (!tCtx) return;

      tCtx.fillStyle = viewMode === 'satellite' || viewMode === 'biome' ? '#050505' : '#000000';
      tCtx.fillRect(0, 0, outWidth, outHeight);

      const output = tCtx.getImageData(0, 0, outWidth, outHeight);
      const outData = output.data;

      rasterizeDymaxionSource({
        sourceData: srcData,
        sourceWidth: srcWidth,
        sourceHeight: srcHeight,
        outputData: outData,
        outputWidth: outWidth,
        outputHeight: outHeight,
        canvasWidth,
        canvasHeight,
        renderDpr,
        layout: dymaxionLayout,
        lon: dymaxionLon,
        lat: dymaxionLat,
        roll: dymaxionRoll,
      });

      tCtx.putImageData(output, 0, 0);

      const dymaxionMarkerProject = (position: { x: number; y: number; z: number }): [number, number] | null => {
        const projected = projectDymaxionPoint(
          [position.x, position.y, position.z],
          dymaxionLayout,
          canvasWidth,
          canvasHeight,
          dymaxionLon,
          dymaxionLat,
          dymaxionRoll,
        );
        return projected ? [projected[0] * renderDpr, projected[1] * renderDpr] : null;
      };

      if (labelVisibility.markers && world?.markers) {
        drawMarkerPins(tCtx, world.markers, dymaxionMarkerProject, Math.max(2, 2.5 * renderDpr));
      }

      drawMapLabels(
        tCtx,
        mapLabels,
        dymaxionMarkerProject,
        renderDpr,
        labelVisibility,
        { theme: labelTheme },
      );

      // Ruler arc: projected per-sample (not baked into the source equirect
      // raster like rivers) so a jump across net faces can break the polyline
      // instead of drawing a spurious line straight across the canvas.
      if (rulerArc && rulerArc.length > 1) {
        tCtx.save();
        tCtx.strokeStyle = '#fbbf24';
        tCtx.lineWidth = Math.max(1.5, 2 * renderDpr);
        tCtx.globalAlpha = 0.9;
        tCtx.beginPath();
        let hasPoint = false;
        let lastPt: [number, number] | null = null;
        const jumpThreshold = canvasWidth * renderDpr * 0.15;
        rulerArc.forEach((p) => {
          const projected = projectDymaxionPoint(
            [p.x, p.y, p.z],
            dymaxionLayout,
            canvasWidth,
            canvasHeight,
            dymaxionLon,
            dymaxionLat,
            dymaxionRoll,
          );
          if (!projected) {
            lastPt = null;
            return;
          }
          const pt: [number, number] = [projected[0] * renderDpr, projected[1] * renderDpr];
          const isJump = lastPt !== null && Math.hypot(pt[0] - lastPt[0], pt[1] - lastPt[1]) > jumpThreshold;
          if (!hasPoint || isJump) {
            tCtx.moveTo(pt[0], pt[1]);
            hasPoint = true;
          } else {
            tCtx.lineTo(pt[0], pt[1]);
          }
          lastPt = pt;
        });
        tCtx.stroke();
        tCtx.restore();
      }

      // Update the main offscreen canvas only once the heavy rendering is complete
      offscreen.width = outWidth;
      offscreen.height = outHeight;
      const finalCtx = offscreen.getContext('2d');
      if (finalCtx) {
        finalCtx.drawImage(tempCanvas, 0, 0);
        setRenderCount(c => c + 1);
      }
      return;
    }

    if (!projection) return;
    offscreen.width = Math.max(1, Math.floor(size.width * renderDpr));
    offscreen.height = Math.max(1, Math.floor(size.height * renderDpr));
    // Reset ctx state after resize
    ctx.setTransform(renderDpr, 0, 0, renderDpr, 0, 0);
    ctx.translate(size.width, 0);
    ctx.scale(-1, 1);
    
    const pathGenerator = d3.geoPath(projection, ctx);

    if (styled) {
      // The context carries the horizontal flip applied just above, so the
      // substrate is told it is mirrored — only glyphs need to compensate. The
      // cached cellPaths are un-flipped CSS-px and drawn THROUGH this same
      // flipped/DPR context, so the mirror and DPR apply automatically (see
      // invariant 2 — no flip is baked into the cache).
      const sub = new Canvas2DSubstrate(
        ctx, pathGenerator as unknown as (o: unknown) => unknown, size.width, size.height, true,
        geometryCache?.cellPaths,
      );
      runStyle(style, {
        world, viewMode, widthPx: size.width, heightPx: size.height,
        glyphs, shadeMap, coastlines,
        colorCtx: {
          seaLevel: world.params.seaLevel, factionColors, cultureColors, religionColors,
        },
        geometryCache: geometryCache ?? undefined,
        colorCache: colorCache ?? undefined,
      }, sub);
    } else {
    for (let i = 0; i < world.cells.length; i++) {
        const feature = world.geoJson?.features?.[i];
      if (!feature || !feature.geometry) continue;
        const color = getCellColor(world.cells[i], viewMode, {
          seaLevel: world.params.seaLevel,
          factionColors,
          cultureColors,
          religionColors,
          seasonalDelta: seasonalTemperatureDelta(world.cells[i], world.params),
        });
        if (shadeMap) color.multiplyScalar(shadeMap[i]);
        const hexColor = '#' + color.getHexString();
      ctx.beginPath();
      pathGenerator(feature);
      ctx.fillStyle = hexColor;
      ctx.strokeStyle = hexColor;
      ctx.lineWidth = 1;
      ctx.fill();
      ctx.stroke();
    }
    }

    drawContourPaths(ctx, pathGenerator, contourSegments, Math.max(0.75, 1.5 / qualityDpr));

    // Draw Grid
    if (showGrid) {
      ctx.strokeStyle = 'rgba(255,255,255,0.15)';
      ctx.lineWidth = 1;
      ctx.beginPath();
      pathGenerator(d3.geoGraticule().step([10, 10])());
      ctx.stroke();
    }

    // Draw Rivers
    if (showRivers && world.rivers) {
      drawRiverPaths(ctx, world.rivers, projection, 1.5 / qualityDpr, overlayInk);
    }

    // Draw Routes (C3) — same antimeridian-jump pattern as rivers above.
    if (showRoutes && world.routes) {
      ctx.globalAlpha = 0.9;
      world.routes.forEach(route => {
        if (route.path.length < 2) return;
        ctx.strokeStyle = route.kind === 'road' ? '#c8a25a' : '#5eb8c8';
        ctx.lineWidth = (route.kind === 'road' ? 1.4 : 1.2) / qualityDpr;
        ctx.setLineDash(route.kind === 'searoute' ? [5 / qualityDpr, 4 / qualityDpr] : []);
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
      ctx.globalAlpha = 1.0;
    }

    // Draw ocean currents (F2) — arrows over ocean cells, warm/cold SST tint.
    if (showCurrents && world.currents && projection) {
      drawCurrents2D(ctx, world, projection, qualityDpr);
    }

    // Ruler arc — same antimeridian-jump projection pattern as rivers above.
    if (rulerArc && rulerArc.length > 1) {
      ctx.strokeStyle = '#fbbf24';
      ctx.lineWidth = 2 / qualityDpr;
      ctx.globalAlpha = 0.9;
      ctx.beginPath();
      let lastLon: number | null = null;
      rulerArc.forEach((p, i) => {
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
      ctx.globalAlpha = 1.0;
    }

    drawFactionBorders(
      ctx,
      pathGenerator,
      factionBorders,
      Math.max(1, 2 / qualityDpr),
      overlayInk,
    );

    if (viewMode === 'political') {
      ctx.save();
      ctx.globalAlpha = 0.25;
      ctx.strokeStyle = '#ffffff';
      ctx.lineWidth = 0.5;
      for (let i = 0; i < world.cells.length; i++) {
        const feature = world.geoJson?.features?.[i];
        if (!feature || !feature.geometry) continue;
        ctx.beginPath();
        pathGenerator(feature);
        ctx.stroke();
      }
      ctx.restore();
    }

    if (highlightCellId !== null) {
      const feature = world.geoJson?.features?.[highlightCellId];
      if (feature && feature.geometry) {
        ctx.save();
        ctx.strokeStyle = '#f9a8a8';
        ctx.lineWidth = 3 / qualityDpr;
        ctx.lineJoin = 'round';
        ctx.lineCap = 'round';
        ctx.beginPath();
        pathGenerator(feature);
        ctx.stroke();
        ctx.restore();
      }
    }

    const mercatorMarkerProject = (position: { x: number; y: number; z: number }): [number, number] | null => {
      const projected = projection(toLonLat([position.x, position.y, position.z]));
      return projected ? [size.width - projected[0], projected[1]] : null;
    };

    // mercatorMarkerProject applies its own manual mirror (size.width - x), so
    // it must run under the same *absolute* (un-mirrored) transform as the
    // label pass below — the ambient ctx transform above still carries the
    // translate/scale(-1,1) cell-drawing mirror, which would double it up.
    if (labelVisibility.markers && world?.markers) {
      ctx.save();
      ctx.setTransform(renderDpr, 0, 0, renderDpr, 0, 0);
      drawMarkerPins(ctx, world.markers, mercatorMarkerProject, Math.max(1.5, 2.5 / qualityDpr));
      ctx.restore();
    }

    if (mapLabels.length > 0) {
      ctx.save();
      ctx.setTransform(renderDpr, 0, 0, renderDpr, 0, 0);
      drawMapLabels(
        ctx,
        mapLabels,
        mercatorMarkerProject,
        // Read via ref: the effect re-runs on the qualityDpr settle cycle, so
        // LOD tracks settled zoom without redrawing every cell per wheel tick.
        scaleRef.current,
        labelVisibility,
        { theme: labelTheme },
      );
      ctx.restore();
    }

    // Signal the blit effect (deps include renderCount) to copy the freshly
    // drawn offscreen to the visible canvas. The dymaxion path bumps this at
    // line ~585; the mercator path must too, or overlay toggles (currents, grid,
    // rivers, …) redraw the offscreen but never reach the screen until a pan/zoom.
    setRenderCount(c => c + 1);
  }, [
    // A3: without these the canvas keeps the previous style's pixels — the
    // offscreen buffer is redrawn only when this effect re-runs.
    style,
    labelTheme,
    overlayInk,
    styled,
    glyphs,
    dymaxionGlyphs,
    coastlines,
    projection,
    size.width,
    size.height,
    world,
    viewMode,
    qualityDpr,
    highlightCellId,
    projectionType,
    dymaxionLayout,
    dymaxionLon,
    dymaxionLat,
    dymaxionRoll,
    showGrid,
    showRivers,
    showRoutes,
    showCurrents,
    factionBorders,
    mapLabels,
    labelVisibility,
    shadeMap,
    contourSegments,
    factionColors,
    cultureColors,
    religionColors,
    rulerArc,
    geometryCache,
    colorCache,
  ]);

  useEffect(() => {
    if (!world || !size.width || !size.height || projectionType !== 'dymaxion') {
      pickCtxRef.current = null;
      return;
    }
    const pickCanvas = pickCanvasRef.current ?? document.createElement('canvas');
    pickCanvasRef.current = pickCanvas;
    pickCanvas.width = size.width;
    pickCanvas.height = size.height;
    const ctx = pickCanvas.getContext('2d', { willReadFrequently: true });
    if (!ctx) return;
    pickCtxRef.current = ctx;

    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, size.width, size.height);

    const srcWidth = Math.max(1, Math.floor(size.width));
    const srcHeight = Math.max(1, Math.round(srcWidth / 2));
    const source = document.createElement('canvas');
    source.width = srcWidth;
    source.height = srcHeight;
    const srcCtx = source.getContext('2d', { willReadFrequently: true });
    if (!srcCtx) return;

    const sourceProjection = d3.geoEquirectangular().fitSize(
      [srcWidth, srcHeight],
      { type: 'Sphere' } as d3.GeoPermissibleObjects,
    );
    const pathGenerator = d3.geoPath(sourceProjection, srcCtx);

    srcCtx.fillStyle = 'rgb(0,0,0)';
    srcCtx.fillRect(0, 0, srcWidth, srcHeight);
    srcCtx.save();
    srcCtx.translate(srcWidth, 0);
    srcCtx.scale(-1, 1);

    for (let i = 0; i < world.cells.length; i++) {
      const feature = world.geoJson?.features?.[i];
      if (!feature || !feature.geometry) continue;
      const id = i + 1;
      const r = id & 255;
      const g = (id >> 8) & 255;
      const b = (id >> 16) & 255;
      const pickColor = `rgb(${r},${g},${b})`;
      srcCtx.beginPath();
      pathGenerator(feature);
      srcCtx.fillStyle = pickColor;
      srcCtx.strokeStyle = pickColor;
      srcCtx.lineWidth = 1;
      srcCtx.fill();
      srcCtx.stroke();
    }

    srcCtx.restore();

    ctx.fillStyle = 'rgb(0,0,0)';
    ctx.fillRect(0, 0, size.width, size.height);
    const srcImage = srcCtx.getImageData(0, 0, srcWidth, srcHeight);
    const output = ctx.getImageData(0, 0, size.width, size.height);
    rasterizeDymaxionSource({
      sourceData: srcImage.data,
      sourceWidth: srcWidth,
      sourceHeight: srcHeight,
      outputData: output.data,
      outputWidth: size.width,
      outputHeight: size.height,
      canvasWidth: size.width,
      canvasHeight: size.height,
      renderDpr: 1,
      layout: dymaxionLayout,
      lon: dymaxionLon,
      lat: dymaxionLat,
      roll: dymaxionRoll,
    });
    ctx.putImageData(output, 0, 0);
    setPickVersion(v => v + 1); // see the ref's own note: readers cannot see a ref
    // The pick buffer encodes cell IDs from the world's structure (cells +
    // geoJson are stable references across paint strokes), so keying on
    // world.cells instead of world identity skips a full-canvas per-pixel
    // reprojection on every stroke event while painting
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [
    size.width,
    size.height,
    world?.cells,
    projectionType,
    dymaxionLayout,
    dymaxionLon,
    dymaxionLat,
    dymaxionRoll,
  ]);

  useEffect(() => {
    const canvas = canvasRef.current;
    const offscreen = offscreenRef.current;
    if (!canvas || !offscreen || !size.width || !size.height) return;
    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    const displayDpr = qualityDpr;
    canvas.width = Math.max(1, Math.floor(size.width * displayDpr));
    canvas.height = Math.max(1, Math.floor(size.height * displayDpr));

    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, canvas.width, canvas.height);

    // The desk the map lies on. SCREEN SPACE, like the scale bar below — the
    // offscreen bitmap is blitted under the pan/zoom transform just after this,
    // so a backdrop painted into it would shrink and slide with the paper
    // instead of filling the viewport. (It was a style pass first, for exactly
    // one round; this is why it is not.)
    //
    // Unconditional: `styled` gates the STYLE PASSES, and the desk is not one.
    // Every style, `default` included, letterboxes under a projection that does
    // not fill its canvas, and every style should fill that margin.
    paintDesk(ctx, style.palette, canvas.width, canvas.height, deskTexture);

    ctx.setTransform(displayDpr * scale, 0, 0, displayDpr * scale, displayDpr * offset.x, displayDpr * offset.y);
    // Drop shadow under the sheet, so it rests on the desk rather than floating
    // in it. Set on the blit itself: canvas shadows apply to drawImage. Also
    // unconditional — it is the only thing separating a near-black default map
    // from the desk it sits on.
    ctx.shadowColor = style.palette.deskShadow;
    ctx.shadowBlur = 28;
    ctx.shadowOffsetY = 6;
    ctx.drawImage(offscreen, 0, 0, size.width, size.height);
    ctx.shadowColor = 'transparent';
    ctx.shadowBlur = 0;
    ctx.shadowOffsetY = 0;

  }, [size.width, size.height, scale, offset.x, offset.y, qualityDpr, viewMode, world?.params.seed, world, projectionType, renderCount, projection, rulerArc, styled, style.palette, deskTexture]);

  /**
   * Screen point -> cell id, by reading the Dymaxion pick buffer.
   *
   * Declared here, above the scale bar, rather than down with the other pointer
   * helpers: the scale bar needs it too. It is the ONLY inverse the Dymaxion net
   * has — see the scale bar's own note.
   */
  const getPickBufferCellId = useCallback((mapX: number, mapY: number) => {
    if (!pickCtxRef.current || !world) return null;
    const data = pickCtxRef.current.getImageData(Math.floor(mapX), Math.floor(mapY), 1, 1).data;
    const id = data[0] + (data[1] << 8) + (data[2] << 16);
    if (id <= 0 || id > world.cells.length) return null;
    return id - 1;
    // pickVersion is not read here — it is the whole point. Listing it gives
    // this callback a new identity when the buffer is rebuilt, which is what
    // makes the scale bar's memo recompute against a buffer that now exists.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [world, pickVersion]);

  /**
   * Ground scale at the centre of the viewport.
   *
   * Was drawn onto the canvas at a hardcoded x of 224, chosen to clear a
   * BiomeLegend panel that has since moved: the shell's panels float ABOVE a
   * full-bleed canvas, so the bar sat underneath the left sidebar and could not
   * be seen. As DOM in the same corner as the pan/zoom hint it cannot be
   * occluded by a panel laid out in the same coordinate space, and it stops
   * being a magic number that goes stale whenever the chrome moves.
   *
   * The value DOES change as you pan north or south, and that is correct rather
   * than a glitch — Mercator and equirectangular both stretch with latitude, so
   * a bar valid at the equator is wrong at 60 degrees. The latitude it applies
   * to is shown beside it, so the change reads as information.
   *
   * **Dymaxion has one too, and this used to say it could not.** The old reason
   * — "net distortion varies per face, so one figure would be a lie" — had it
   * backwards. Measured on this app's own net, linear scale varies 13% peak to
   * peak over the whole globe; Mercator, which already carried a bar, is out by
   * 100% at 60 degrees. A single figure is MORE defensible on the Dymaxion than
   * on the view that had one. What actually blocked it was mechanical.
   * `computeScaleBar` needs a centre lon/lat, every other view gets one from
   * `projection.invert`, and Dymaxion has no analytic inverse.
   *
   * It has a per-pixel one: the pick buffer built for cell selection. Read it at
   * the viewport centre, take that cell's centre, and the missing lon/lat is
   * there. The forward direction is `projectDymaxionPoint`, in the same CSS
   * pixel space the pick buffer and `centerMapX/Y` already use.
   *
   * The genuine hazard is the CUT, not the distortion, and it is handled by
   * `computeInterruptedScaleBar` — see there for why a fold fails quietly rather
   * than absurdly. Centre the viewport on a fold, or on the background between
   * faces, and no bar is drawn, which is the honest answer.
   */
  const scaleBar = useMemo(() => {
    if (!world) return null;
    const centerMapX = (size.width / 2 - offset.x) / scale;
    const centerMapY = (size.height / 2 - offset.y) / scale;

    if (projectionType === 'dymaxion') {
      if (centerMapX < 0 || centerMapY < 0 || centerMapX >= size.width || centerMapY >= size.height) return null;
      const cellId = getPickBufferCellId(centerMapX, centerMapY);
      if (cellId === null) return null; // background between faces
      const c = world.cells[cellId].center;
      const centerLonLat = toLonLat([c.x, c.y, c.z]);
      const forward = (lonLat: [number, number]): [number, number] | null =>
        projectDymaxionPoint(
          lonLatToPoint3(lonLat), dymaxionLayout,
          size.width, size.height,
          dymaxionLon, dymaxionLat, dymaxionRoll,
        );
      const dymInfo = computeInterruptedScaleBar(forward, centerLonLat, world.params.planetRadius);
      if (!dymInfo) return null;
      const dym = niceScaleBarLength(dymInfo.pixelsPerKm * scale, 140);
      if (dym.km <= 0) return null;
      const dymLat = Math.round(centerLonLat[1]);
      return { km: dym.km, px: dym.px, lat: `${Math.abs(dymLat)}\u00b0${dymLat < 0 ? 'S' : 'N'}` };
    }

    if (!projection) return null;
    const centerLonLat = projection.invert?.([size.width - centerMapX, centerMapY]);
    if (!centerLonLat || !Number.isFinite(centerLonLat[0]) || !Number.isFinite(centerLonLat[1])) return null;
    const info = computeScaleBar(projection, centerLonLat, world.params.planetRadius);
    if (!info) return null;
    const { km, px } = niceScaleBarLength(info.pixelsPerKm * scale, 140);
    if (km <= 0) return null;
    const lat = Math.round(centerLonLat[1]);
    return { km, px, lat: `${Math.abs(lat)}\u00b0${lat < 0 ? 'S' : 'N'}` };
  }, [
    world, projection, projectionType, size.width, size.height, offset.x, offset.y, scale,
    getPickBufferCellId, dymaxionLayout, dymaxionLon, dymaxionLat, dymaxionRoll,
  ]);

  const handleWheel = useCallback((event: WheelEvent) => {
    event.preventDefault();
    setIsInteracting(true);
    if (wheelTimer.current) window.clearTimeout(wheelTimer.current);
    wheelTimer.current = window.setTimeout(() => {
      setIsInteracting(false);
      wheelTimer.current = null;
    }, 180);
    const rect = canvasRef.current?.getBoundingClientRect();
    if (!rect) return;
    const mx = event.clientX - rect.left;
    const my = event.clientY - rect.top;

    const prevScale = scaleRef.current;
    const zoomFactor = event.deltaY < 0 ? 1.1 : 0.9;
    const nextScale = clamp(prevScale * zoomFactor, 0.6, 6);

    const prevOffset = offsetRef.current;
    const worldX = (mx - prevOffset.x) / prevScale;
    const worldY = (my - prevOffset.y) / prevScale;

    const nextOffsetX = mx - worldX * nextScale;
    const nextOffsetY = my - worldY * nextScale;

    setScale(nextScale);
    setOffset({ x: nextOffsetX, y: nextOffsetY });
  }, []);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const listener = (event: WheelEvent) => { handleWheel(event); };
    canvas.addEventListener('wheel', listener, { passive: false });
    return () => { canvas.removeEventListener('wheel', listener); };
  }, [handleWheel]);

  const isPaintMode = editMode !== 'off' && editMode !== 'world-edit';

  const getCellIdAtMapPoint = useCallback((mapX: number, mapY: number) => {
    if (!world) return null;

    if (projectionType === 'dymaxion') {
      return getPickBufferCellId(mapX, mapY);
    }

    if (!projection || !pickQuadtree) return null;
    // Un-flip x the same way the old invert path did (size.width - mapX), then
    // query the quadtree of projected cell centres. Nearest-projected-centre
    // replaces the O(cells) geodesic scan (getNearestCellId).
    return findCellIdAtPoint(pickQuadtree, size.width - mapX, mapY);
  }, [getPickBufferCellId, projection, projectionType, size.width, world, pickQuadtree]);

  const paintPickAt = useCallback((clientX: number, clientY: number, phase: 'start' | 'stroke' | 'end') => {
    if (!canvasRef.current || !onPaint) return;
    const rect = canvasRef.current.getBoundingClientRect();
    const mapX = (clientX - rect.left - offset.x) / scale;
    const mapY = (clientY - rect.top - offset.y) / scale;
    if (mapX < 0 || mapY < 0 || mapX >= size.width || mapY >= size.height) return;
    const cellId = getCellIdAtMapPoint(mapX, mapY);
    if (cellId === null) return;
    lastPaintCell.current = cellId;
    onPaint(cellId, phase);
  }, [getCellIdAtMapPoint, offset.x, offset.y, scale, size.width, size.height, onPaint]);

  const handleMouseDown = (e: React.MouseEvent<HTMLCanvasElement>) => {
    if (e.button === 1) e.preventDefault();
    dragging.current = true;
    dragDistance.current = 0;
    lastPos.current = { x: e.clientX, y: e.clientY };
    setIsInteracting(true);

    if (isPaintMode && !isSpaceHeld && e.button !== 1) {
      const isRight = e.button === 2;
      if (isRight) {
        // Right-click: send start with isRightClick=true, don't begin a stroke
        if (!canvasRef.current) return;
        const rect = canvasRef.current.getBoundingClientRect();
        const mapX = (e.clientX - rect.left - offset.x) / scale;
        const mapY = (e.clientY - rect.top - offset.y) / scale;
        if (mapX >= 0 && mapY >= 0 && mapX < size.width && mapY < size.height) {
          const cellId = getCellIdAtMapPoint(mapX, mapY);
          if (cellId !== null && onPaint) onPaint(cellId, 'start', true);
        }
        return;
      }
      isPaintingRef.current = true;
      paintPickAt(e.clientX, e.clientY, 'start');
      paintPickAt(e.clientX, e.clientY, 'stroke');
    }
  };

  const handleMouseMove = (e: React.MouseEvent<HTMLCanvasElement>) => {
    // Track mouse for brush circle overlay
    if (isPaintMode) {
      const rect = canvasRef.current?.getBoundingClientRect();
      if (rect) setMousePos({ x: e.clientX - rect.left, y: e.clientY - rect.top });
    }

    if (!dragging.current) return;
    const dx = e.clientX - lastPos.current.x;
    const dy = e.clientY - lastPos.current.y;
    lastPos.current = { x: e.clientX, y: e.clientY };
    dragDistance.current += Math.abs(dx) + Math.abs(dy);

    if (isPaintMode && isPaintingRef.current && !isSpaceHeld) {
      paintPickAt(e.clientX, e.clientY, 'stroke');
      return; // skip pan during paint stroke
    }

    const next = { x: offset.x + dx, y: offset.y + dy };
    pendingOffset.current = next;
    if (rafId.current === null) {
      rafId.current = requestAnimationFrame(() => {
        rafId.current = null;
        if (pendingOffset.current) {
          setOffset(pendingOffset.current);
          pendingOffset.current = null;
        }
      });
    }
  };

  const pickAt = (clientX: number, clientY: number, clearOnMiss = true) => {
    if (!canvasRef.current) return;
    const rect = canvasRef.current.getBoundingClientRect();
    const mapX = (clientX - rect.left - offset.x) / scale;
    const mapY = (clientY - rect.top - offset.y) / scale;
    if (mapX < 0 || mapY < 0 || mapX >= size.width || mapY >= size.height) {
      if (clearOnMiss) onInspect(null);
      return;
    }
    const cellId = getCellIdAtMapPoint(mapX, mapY);
    if (cellId === null) {
      if (clearOnMiss) onInspect(null);
      return;
    }
    onInspect(cellId);
  };

  const handleHover = (e: React.MouseEvent<HTMLCanvasElement>) => {
    if (inspectMode !== 'hover' || dragging.current) return;
    pickAt(e.clientX, e.clientY);
  };

  const endDrag = () => {
    if (isPaintMode && isPaintingRef.current) {
      if (lastPaintCell.current !== null && onPaint) {
        onPaint(lastPaintCell.current, 'end');
      }
      isPaintingRef.current = false;
      lastPaintCell.current = null;
    }
    dragging.current = false;
    setIsInteracting(false);
  };

  const handleMouseUp = (e: React.MouseEvent<HTMLCanvasElement>) => {
    endDrag();
    if (inspectEnabled && dragDistance.current < 6) {
      pickAt(e.clientX, e.clientY, false);
    }
  };

  return (
    <div ref={containerRef} className="w-full h-full bg-black relative">
      <canvas
        ref={canvasRef}
        className="w-full h-full cursor-grab active:cursor-grabbing"
        onMouseDown={handleMouseDown}
        onMouseMove={(e) => { handleMouseMove(e); handleHover(e); }}
        onMouseUp={handleMouseUp}
        onMouseLeave={() => { endDrag(); setMousePos(null); if (inspectMode === 'hover') onInspect(null); }}
        onAuxClick={e => { if (e.button === 1) e.preventDefault(); }}
        onContextMenu={e => e.preventDefault()}
      />
      {/* Brush size circle overlay */}
      {isPaintMode && mousePos && world && (() => {
        const avgCellPx = Math.sqrt((size.width * size.height) / world.cells.length / Math.PI);
        const brushPx = Math.max(8, avgCellPx * (brushSize + 0.5) * scale);
        return (
          <div
            className="pointer-events-none absolute border-2 border-white/70 rounded-full"
            style={{
              left: mousePos.x - brushPx,
              top: mousePos.y - brushPx,
              width: brushPx * 2,
              height: brushPx * 2,
            }}
          />
        );
      })()}
      {!world && (
        <div className="absolute inset-0 flex items-center justify-center text-ink-strong/50">
          Forging World...
        </div>
      )}
      {scaleBar && (
        <div className="absolute bottom-11 right-4 flex items-center gap-2 text-[10px] bg-black/60 border border-white/10 px-2 py-1 text-ink-soft">
          <span
            className="inline-block border-x border-b border-current h-1.5"
            style={{ width: `${scaleBar.px}px` }}
          />
          <span>{scaleBar.km.toLocaleString()} km at {scaleBar.lat}</span>
        </div>
      )}
      {world && (
        <div className="absolute bottom-4 right-4 text-[10px] bg-black/60 border border-white/10 px-2 py-1 text-ink-soft">
          {DISPLAY_MODES.find(m => m.mode === projectionType)?.label ?? '2D Map'} • Scroll to zoom, drag to pan
        </div>
      )}
    </div>
  );
};

export default Map2D;
