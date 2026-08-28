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
 * D8b — physical climate coupling constants (gated by `physicalClimate`).
 *
 * `LAPSE_RATE_C_PER_KM` is the standard environmental lapse rate: air cools
 * 6.5 °C for every kilometre of altitude. It replaces the invented `* 60`
 * multiplier on normalized height in worldGen's temperature finalize.
 *
 * The OROG_* constants scale the orographic (rain-shadow) moisture term by the
 * real barrier metres crossed. Windward ascent boosts rain (capped so coasts do
 * not saturate); leeward descent dries it (floored so a single edge cannot lose
 * everything). Tuned empirically against TWO measured targets, 3 seeds at 20k:
 *   (1) land moisture<0.15 share in the ~30-36% accept band (spec §6);
 *   (2) a functioning rain shadow — cells behind an upwind barrier (>500 m) are
 *       measurably drier than exposed cells (a windward/leeward contrast metric).
 * Final: windward 0.85, leeward floor 0.5 / per-km 0.3 → 3-seed avg 35.4% dry
 * share AND a shadow contrast of ~0.135 (shadowed ~0.26-0.34 vs exposed
 * ~0.42-0.44 mean moisture).
 *
 * NOTE — an earlier Task 5 tune hit 32.4% by nearly DISABLING leeward drying
 * (floor 0.95 / per-km 0.02), which killed rain shadows (contrast ~0.08 — no
 * better than the incidental inland-dryness gradient). That aggregate-only target
 * was smeared, not located. This tune trades ~3 pts of dry share (32→35%) for
 * real, correctly-placed rain-shadow aridity — steppe rises only to ~18-25%
 * (still below pre-D8b ~26.5-29%). The root cause of inland dryness is the 8-pass
 * transport under-delivering moisture (worldGen.ts:611 land decay + fixed pass
 * count), which the orographic knob can locate but not add water to — addressing
 * that is a follow-up, not part of this pair.
 */
export const LAPSE_RATE_C_PER_KM = 6.5;
export const OROG_WINDWARD_PER_KM = 0.85;
export const OROG_WINDWARD_CAP = 2.5;
export const OROG_LEEWARD_PER_KM = 0.3;
export const OROG_LEEWARD_FLOOR = 0.5;

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
