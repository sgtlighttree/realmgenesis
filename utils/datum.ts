/**
 * D8a — presentation datum (units).
 *
 * The height field is a normalized 0-1 value with `seaLevel` as the coastline.
 * This module turns it into real metres above/below sea level so the app reads
 * in kilometres horizontally AND metres vertically, instead of kilometres
 * horizontally and percent vertically (the inconsistency ROADMAP D8 names).
 *
 * The datum is a FIXED maximum that `mountainHeight` distributes terrain within,
 * NOT a value derived from the world — that is deliberate. A datum that varied
 * with the world would report a different altitude for the same cell when a
 * slider moved, making altitudes incomparable between worlds. `maxElevationM` is
 * user-adjustable but defaults to a fixed 9000, and adjusting it is DISPLAY-ONLY:
 * it rescales every readout and changes no terrain (that coupling is D8b).
 *
 * `MAX_DEPTH_M` is a plain constant, NOT adjustable, on purpose: `seafloorDepth`
 * already exists as a generation-side ocean-depth datum, so a second adjustable
 * depth here would duplicate it.
 */

/** Default presentation datum in metres. Single source of truth. */
export const DEFAULT_MAX_ELEVATION_M = 9000;

/** Fixed ocean-depth datum in metres (~Mariana Trench). Not user-adjustable. */
export const MAX_DEPTH_M = 11000;

/**
 * Convert a normalized 0-1 cell height to metres relative to sea level.
 * Land (height > seaLevel) scales toward +maxElevationM at height 1.
 * Ocean (height < seaLevel) scales toward -MAX_DEPTH_M at height 0.
 * The coastline (height === seaLevel) is exactly 0.
 */
export const elevationMetres = (
  height: number,
  seaLevel: number,
  maxElevationM: number = DEFAULT_MAX_ELEVATION_M,
): number => {
  if (height >= seaLevel) {
    return ((height - seaLevel) / (1 - seaLevel)) * maxElevationM;
  }
  return ((height - seaLevel) / seaLevel) * MAX_DEPTH_M;
};

/**
 * Human-readable elevation for a cell: "3,000 m", "-2,400 m", or "Sea level".
 */
export const formatElevation = (
  height: number,
  seaLevel: number,
  maxElevationM: number = DEFAULT_MAX_ELEVATION_M,
): string => {
  const m = Math.round(elevationMetres(height, seaLevel, maxElevationM));
  if (m === 0) return 'Sea level';
  return `${m.toLocaleString('en-US')} m`;
};
