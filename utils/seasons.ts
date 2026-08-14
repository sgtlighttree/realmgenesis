import { Cell, WorldParams } from '../types';

// D1 — Seasonal cycle. Pure, THREE-free, worldGen-free helpers shared by the
// generator (annual-mean temperature) and the render/color layer (per-season
// display excursion). See docs/superpowers/specs/2026-08-14-d1-seasonal-cycle-design.md.
//
// The season slider is orbital position s ∈ [0,1). s = 0.5 is the neutral point
// (declination 0 → the shown world equals the annual mean). axialTilt is the
// amplitude of the seasonal excursion, no longer a static climate offset.

// D3: seawater freezing point (°C). A physical constant, not a creative knob —
// the iciness levers are poleTemperature/baseTemperature/axialTilt (all live +
// regenerating). Sea-ice is a render overlay: water cells whose *seasonal*
// temperature falls below this render as ice (see utils/colors.ts).
export const SEAWATER_FREEZE_C = -2;

// Samples used to orbit-average the latitude temperature curve. A smooth
// sinusoid-composed function converges fast; 96 keeps the annual mean stable to
// ~1e-3 °C so the seasonal excursion is symmetric to well under a paramLiveness
// tick, at a per-cell cost of ~96 flops (trivial even at 200k cells).
const ORBIT_SAMPLES = 96;

// Latitude temperature curve, matching worldGen's climate pass exactly.
// ψ is geometric latitude in radians. Quadratic falloff equator → pole.
const tLat = (psiRad: number, base: number, pole: number): number => {
  const r = Math.abs(psiRad) / (Math.PI / 2);
  return base * (1 - r * r) + pole * (r * r);
};

// Subsolar declination (radians) at orbital position `season`, for a given
// axial tilt in degrees. δ(s) = tilt·sin(2π·s): 0 at the equinox neutral point
// (s = 0.5 and s = 0), ±tilt at the solstices (s = 0.25 / 0.75).
export const seasonalDeclination = (season: number, axialTiltDeg: number): number => {
  const tilt = (axialTiltDeg || 0) * (Math.PI / 180);
  return tilt * Math.sin(2 * Math.PI * season);
};

// Orbit-averaged annual-mean latitude temperature at geometric latitude φ.
// Because tLat is quadratic, this shifts with axial tilt (Jensen's inequality) —
// which is what keeps axialTilt a live generation parameter after the seasonal
// excursion is moved into the render layer.
export const annualMeanLatTemp = (phiRad: number, params: WorldParams): number => {
  const tilt = (params.axialTilt || 0) * (Math.PI / 180);
  const base = params.baseTemperature;
  const pole = params.poleTemperature;
  if (tilt === 0) return tLat(phiRad, base, pole);
  let sum = 0;
  for (let k = 0; k < ORBIT_SAMPLES; k++) {
    const decl = tilt * Math.sin((2 * Math.PI * k) / ORBIT_SAMPLES);
    sum += tLat(phiRad - decl, base, pole);
  }
  return sum / ORBIT_SAMPLES;
};

// Seasonal temperature excursion (°C) added to the stored cell.temperature by
// the render layer to show the climate at orbital position params.season.
//
// Anchored to the EQUINOX (δ=0), not the annual mean: delta = Tlat(φ−δ(s)) −
// Tlat(φ). This is exactly 0 — and continuously so — at the equinox seasons
// (s = 0.5 default, and s = 0), so the neutral view reproduces the canonical
// annual-mean world with no discontinuous "pop" when the slider leaves neutral.
// (Anchoring to the annual mean instead would make every cell jump ~C·tilt²/2 °C
// the instant you nudge off neutral — see the D1 spec.) The stored temperature
// remains the tilt-dependent orbit mean, so axialTilt stays a live generation param.
export const seasonalTemperatureDelta = (cell: Cell, params: WorldParams): number => {
  const season = params.season ?? 0.5;
  if (!params.axialTilt || season === 0.5) return 0; // neutral: canonical world
  const decl = seasonalDeclination(season, params.axialTilt);
  if (decl === 0) return 0;
  const phi = Math.asin(Math.max(-1, Math.min(1, cell.center.y)));
  const base = params.baseTemperature;
  const pole = params.poleTemperature;
  return tLat(phi - decl, base, pole) - tLat(phi, base, pole);
};
