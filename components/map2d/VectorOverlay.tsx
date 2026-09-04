import React, { useCallback, useEffect, useRef } from 'react';
import * as d3 from 'd3';

import { WorldData, LabelVisibility, Point } from '../../types';
import { MapLabel } from '../../utils/labels';
import { OverlayInk } from '../../utils/mapStyle/overlayInk';
import { LabelTheme } from '../../utils/mapStyle/labelTheme';
import { ContourSegment } from '../../utils/shading';
import {
  paintVectorOverlay, PaintVectorOverlayParams, VectorOverlayBoundaryPaths, VectorOverlayToggles,
} from './paintVectorOverlay';

export interface VectorOverlayHandle {
  redraw: () => void;
}

type Size = { width: number; height: number };

export interface VectorOverlayProps {
  size: Size;
  dpr: number;
  /** Live pan/zoom, applied on top of dpr — same values the GL fill's `setTransform` receives (Task 6). */
  scale: number;
  offsetX: number;
  offsetY: number;
  world: WorldData;
  projection: d3.GeoProjection;
  /**
   * Cached, LOD-swappable boundary Path2Ds from `MapGeometryCache`
   * (`coast`/`borders`/`lakes`). The host swaps these on zoom-settle for
   * versions rebuilt at a lower simplification tolerance — see
   * `VectorOverlayBoundaryPaths` in `paintVectorOverlay.ts` for why this is
   * geometry LOD, not pixel LOD.
   */
  boundaryPaths: VectorOverlayBoundaryPaths;
  contourSegments: ContourSegment[];
  mapLabels: MapLabel[];
  labelVisibility: LabelVisibility;
  labelTheme: LabelTheme;
  overlayInk: OverlayInk;
  /** `style.palette.coast` — see the field doc on `PaintVectorOverlayParams.coastColor`. */
  coastColor: string;
  toggles: VectorOverlayToggles;
  /** A5 ruler arc; null when inactive. See `PaintVectorOverlayParams.rulerArc`. */
  rulerArc?: Point[] | null;
  /** Selected/hover cell's cached ring Path2D. See `PaintVectorOverlayParams.highlightPath`. */
  highlightPath?: Path2D | null;
  /** Combined cell-outline path, political mode only. See `PaintVectorOverlayParams.cellOutlinePath`. */
  cellOutlinePath?: Path2D | null;
  /** Exposes `redraw()` to the host — same pattern as `GLFillSurface`'s `handleRef`. */
  handleRef?: React.MutableRefObject<VectorOverlayHandle | null>;
}

/**
 * Transparent per-frame Canvas2D layer over the WebGL cell-fill surface
 * (`GLFillSurface`). Owns its own `<canvas>`, sized `size.width*dpr ×
 * size.height*dpr` with `size.width × size.height` CSS dimensions, and
 * exposes `redraw()` via `handleRef` for the host to call on pan/zoom
 * (Task 6 wires the event handling; this component only draws).
 *
 * `redraw()` sets the SAME composed transform order Map2D's offscreen build
 * uses (mirror ∘ dpr — see docs/map-styles.md §"The mirror trap"), extended
 * with live pan/zoom, then delegates the actual draw sequence to
 * `paintVectorOverlay` so the sequence itself lives in exactly one place.
 */
export const VectorOverlay: React.FC<VectorOverlayProps> = ({
  size, dpr, scale, offsetX, offsetY, world, projection, boundaryPaths, contourSegments,
  mapLabels, labelVisibility, labelTheme, overlayInk, coastColor, toggles,
  rulerArc = null, highlightPath = null, cellOutlinePath = null, handleRef,
}) => {
  const canvasRef = useRef<HTMLCanvasElement>(null);

  const redraw = useCallback(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    // pan/zoom (display) ∘ dpr, matching the current blit's displayDpr*scale
    // form:
    ctx.setTransform(dpr * scale, 0, 0, dpr * scale, dpr * offsetX, dpr * offsetY);
    // mirror (same as the offscreen build): flip X about width.
    ctx.translate(size.width, 0);
    ctx.scale(-1, 1);

    const params: PaintVectorOverlayParams = {
      world, projection, size, toggles, boundaryPaths, contourSegments,
      mapLabels, labelVisibility, labelTheme, overlayInk, coastColor,
      // Line widths divide by `scale` (not `qualityDpr`) so strokes keep
      // constant screen weight across zoom — the vector transform above
      // already carries `dpr`.
      strokeScale: scale,
      rulerArc, highlightPath, cellOutlinePath,
    };
    paintVectorOverlay(ctx, params);
  }, [
    dpr, scale, offsetX, offsetY, size, world, projection, toggles, boundaryPaths,
    contourSegments, mapLabels, labelVisibility, labelTheme, overlayInk, coastColor,
    rulerArc, highlightPath, cellOutlinePath,
  ]);

  // Backing store sized size.width*dpr × size.height*dpr; CSS size stays
  // size.width × size.height via the inline style below.
  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas) return;
    const w = Math.max(1, Math.floor(size.width * dpr));
    const h = Math.max(1, Math.floor(size.height * dpr));
    if (canvas.width !== w) canvas.width = w;
    if (canvas.height !== h) canvas.height = h;
  }, [size.width, size.height, dpr]);

  // Redraw whenever anything the draw sequence reads changes. Task 6 also
  // calls `handleRef.current.redraw()` directly during live pan/zoom, ahead
  // of this effect's own dependency-driven redraw — both paths land on the
  // same function.
  useEffect(() => {
    redraw();
  }, [redraw]);

  useEffect(() => {
    if (!handleRef) return undefined;
    handleRef.current = { redraw };
    return () => {
      if (handleRef.current?.redraw === redraw) {
        handleRef.current = null;
      }
    };
  }, [handleRef, redraw]);

  return (
    <canvas
      ref={canvasRef}
      style={{
        position: 'absolute',
        top: 0,
        left: 0,
        width: size.width,
        height: size.height,
        pointerEvents: 'none',
      }}
    />
  );
};
