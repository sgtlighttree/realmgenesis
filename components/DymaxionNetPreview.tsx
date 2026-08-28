import React, { useEffect, useRef } from 'react';
import * as d3 from 'd3';

import { WorldData, ViewMode, DymaxionSettings } from '../types';
import { buildFactionColorMap, buildCultureColorMap, buildReligionColorMap, getCellColor } from '../utils/colors';
import { seasonalTemperatureDelta } from '../utils/seasons';
import { rasterizeDymaxionSource } from './Map2D';

interface Props {
  world: WorldData | null;
  viewMode: ViewMode;
  settings: DymaxionSettings;
}

// Source equirectangular size (plate carrée, 2:1). rasterizeDymaxionSource
// samples it as lon∈[-180,180]→x, lat∈[90,-90]→y, so the aspect must be 2:1
// and the projection must fit the full sphere with no padding.
const SRC_W = 360;
const SRC_H = 180;
const OUT_W = 240; // CSS px
const OUT_H = 130;

/**
 * A compact, live preview of the ACTUAL Dymaxion net (the same projection the
 * full 2D Dymaxion view uses), not a Mercator with an overlay. It renders the
 * world to an offscreen equirectangular source, then warps that source onto the
 * icosahedron net via the shared rasterizer, keyed off the cage lon/lat/roll so
 * it updates as the cage moves. Debounced like the minimap — hosts must UNMOUNT
 * it to stop the redraw.
 */
const DymaxionNetPreview: React.FC<Props> = ({ world, viewMode, settings }) => {
  const canvasRef = useRef<HTMLCanvasElement>(null);

  useEffect(() => {
    if (!world || !canvasRef.current) return;
    const timer = window.setTimeout(() => {
      const canvas = canvasRef.current;
      if (!canvas) return;
      const dpr = Math.min(2, window.devicePixelRatio || 1);
      const outW = Math.round(OUT_W * dpr);
      const outH = Math.round(OUT_H * dpr);
      canvas.width = outW;
      canvas.height = outH;
      const ctx = canvas.getContext('2d');
      if (!ctx) return;

      // 1) Offscreen equirectangular source (no flip — the rasterizer assumes
      //    a plain plate carrée). fitSize on a 2:1 canvas fits the sphere exactly.
      const src = document.createElement('canvas');
      src.width = SRC_W;
      src.height = SRC_H;
      const sctx = src.getContext('2d');
      if (!sctx) return;
      sctx.fillStyle = '#111';
      sctx.fillRect(0, 0, SRC_W, SRC_H);
      const projection = d3.geoEquirectangular().fitSize([SRC_W, SRC_H], { type: 'Sphere' } as d3.GeoPermissibleObjects);
      const path = d3.geoPath(projection, sctx);
      const factionColors = buildFactionColorMap(world.civData);
      const cultureColors = buildCultureColorMap(world.cultures);
      const religionColors = buildReligionColorMap(world.religions);
      world.cells.forEach((cell, i) => {
        const feature = world.geoJson?.features?.[i];
        if (!feature || !feature.geometry || feature.geometry.coordinates.length === 0) return;
        const color = getCellColor(cell, viewMode, {
          seaLevel: world.params.seaLevel,
          factionColors,
          cultureColors,
          religionColors,
          seasonalDelta: seasonalTemperatureDelta(cell, world.params),
        });
        sctx.beginPath();
        path(feature);
        sctx.fillStyle = '#' + color.getHexString();
        sctx.fill();
      });
      const sourceImage = sctx.getImageData(0, 0, SRC_W, SRC_H);

      // 2) Warp the source onto the icosahedron net.
      const out = ctx.createImageData(outW, outH);
      for (let i = 0; i < out.data.length; i += 4) {
        out.data[i] = 17; out.data[i + 1] = 17; out.data[i + 2] = 17; out.data[i + 3] = 255;
      }
      rasterizeDymaxionSource({
        sourceData: sourceImage.data,
        sourceWidth: SRC_W,
        sourceHeight: SRC_H,
        outputData: out.data,
        outputWidth: outW,
        outputHeight: outH,
        canvasWidth: OUT_W,
        canvasHeight: OUT_H,
        renderDpr: dpr,
        layout: settings.layout,
        lon: settings.lon,
        lat: settings.lat,
        roll: settings.roll,
      });
      ctx.putImageData(out, 0, 0);
    }, 150);
    return () => { window.clearTimeout(timer); };
  }, [world, viewMode, settings.layout, settings.lon, settings.lat, settings.roll]);

  if (!world) return null;

  return (
    <div>
      <div className="text-[10px] text-ink-muted mb-1">Dymaxion net (live)</div>
      <canvas
        ref={canvasRef}
        className="block w-full opacity-90 hover:opacity-100 transition-opacity max-w-full"
      />
    </div>
  );
};

export default DymaxionNetPreview;
