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
 * Hypsometric power curve applied to the above-sea LAND fraction before scaling
 * to metres. The normalized height field sits mostly at mid-range (median land
 * frac ~0.23), so a LINEAR datum reported a median land elevation of ~2 km —
 * Earth is ~840 m with ~70% of land under 1 km. A quadratic curve (frac^2)
 * matches Earth almost exactly (measured: mean 824 m, 72% under 1 km) without
 * touching the height field, generation, or the 3D relief.
 *
 * Ocean is deliberately LEFT LINEAR: real ocean hypsometry is bottom-heavy
 * (most of the floor is 3-6 km abyssal plain), the opposite of land, so the same
 * concave curve made the ocean implausibly shallow (measured median -1,059 m vs
 * Earth's abyssal ~-3,700 m). The land curve alone already softens coastal cliffs
 * from the land side. Verified with scripts/queryWorld.mjs hypsometry.
 */
export const HYPSOMETRIC_EXPONENT = 2.0;

/**
 * Convert a normalized 0-1 cell height to metres relative to sea level.
 * Land (height > seaLevel) rises toward +maxElevationM at height 1, through the
 * HYPSOMETRIC_EXPONENT curve so the metre distribution matches Earth's.
 * Ocean (height < seaLevel) descends LINEARLY toward -MAX_DEPTH_M at height 0.
 * The coastline (height === seaLevel) is exactly 0.
 */
export const elevationMetres = (
  height: number,
  seaLevel: number,
  maxElevationM: number = DEFAULT_MAX_ELEVATION_M,
): number => {
  if (height >= seaLevel) {
    const frac = (height - seaLevel) / (1 - seaLevel);
    return Math.pow(frac, HYPSOMETRIC_EXPONENT) * maxElevationM;
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
