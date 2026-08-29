import * as d3 from 'd3';
import { geoWinkel3, geoRobinson, geoMollweide } from 'd3-geo-projection';

import { ProjectionType } from './export';

/** Every projection except Dymaxion, which is an unfolded net, not a d3 one. */
export type FlatProjectionType = Exclude<ProjectionType, 'dymaxion'>;

/**
 * The one place a d3 projection is constructed.
 *
 * The export path had this list and the on-screen map did not, so `Map2D` could
 * only show Mercator while `exportSVG` already offered six. Two registries would
 * drift — and had, silently, for a while.
 *
 * `fitSize` against the sphere is part of the contract: it is what makes the
 * shader's `(lon+180)/360` assumption and the substrate's sphere clip agree with
 * whatever comes out of here.
 */
export const buildProjection = (
  type: FlatProjectionType,
  width: number,
  height: number,
): d3.GeoProjection => {
  let projection: d3.GeoProjection;
  switch (type) {
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

/**
 * Projections offered as an on-screen 2D view, with their labels.
 *
 * A subset of what the exporter can produce, on purpose. Orthographic shows one
 * hemisphere, which the 3D globe already does better and interactively;
 * Robinson and Mollweide are close enough to Winkel Tripel that three of them in
 * a menu is noise rather than choice.
 *
 * Every one here MUST have a working `invert`, because cell picking projects the
 * pointer back to a lon/lat. Verified for Winkel Tripel: geoWinkel3().invert()
 * round-trips (37.5, -12.25) exactly.
 */
export const SCREEN_PROJECTIONS: ReadonlyArray<{ id: FlatProjectionType; label: string }> = [
  { id: 'equirectangular', label: 'Equirectangular' },
  { id: 'mercator', label: 'Mercator' },
  { id: 'winkeltripel', label: 'Winkel Tripel' },
];
