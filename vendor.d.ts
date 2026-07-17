// Local declaration shims for untyped d3 plugins.
// @types/d3-geo-projection is not published for this registry and
// d3-geo-voronoi ships no types at all; keep these shims minimal and
// typed at the call sites instead.

declare module 'd3-geo-projection' {
  import type { GeoProjection, GeoRawProjection } from 'd3';

  export function geoWinkel3(): GeoProjection;
  export function geoRobinson(): GeoProjection;
  export function geoMollweide(): GeoProjection;

  export interface PolyhedralNode {
    face: [number, number][];
    project: GeoProjection;
    children?: PolyhedralNode[];
  }

  export function geoPolyhedral(
    root: PolyhedralNode,
    faceSelector: (lambda: number, phi: number) => PolyhedralNode,
  ): GeoProjection;

  export function geoGnomonicRaw(): GeoRawProjection;
}

declare module 'd3-geo-voronoi' {
  import type { GeoJsonCollection } from './types';

  export interface GeoVoronoi {
    polygons(): GeoJsonCollection;
    links(): { features: Array<{ geometry: { coordinates: [number, number][] } }> };
  }

  export function geoVoronoi(points: [number, number][]): GeoVoronoi;
}
