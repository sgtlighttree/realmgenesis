import { StarClass } from '../types';

// D5 — Planetary parameters. Only `starClass` ships with a mechanical hook this
// iteration (day length / gravity / moons have no principled hook in the current
// model and are deferred — see docs/ENGINEERING-NOTES.md and ROADMAP D5).
//
// Star class shifts the world's temperature: a hotter class warms the whole
// planet, a cooler one cools it. Applied as a Stefan-Boltzmann-style scale in
// KELVIN (T ∝ L^¼), NOT a multiply on °C — the Kelvin form drives the (negative)
// pole temperature further down for a dimmer star instead of wrongly warming it.
//
// The factors are a STYLIZED insolation range, not literal stellar luminosity
// (a real O star is ~10^5× the Sun and would vaporize any planet outside a
// razor-thin habitable zone). The range is kept narrow so biomes stay varied at
// both extremes (verified by biome-variety census); G = 1.0 is an exact no-op so
// default worlds are byte-identical.
const STAR_TEMP_FACTOR: Record<StarClass, number> = {
  O: 1.07, // blue giant — hottest
  B: 1.05,
  A: 1.035,
  F: 1.018,
  G: 1.0, // Sun-like — default, no-op
  K: 0.965,
  M: 0.93, // red dwarf — coolest
};

export const STAR_CLASSES: StarClass[] = ['O', 'B', 'A', 'F', 'G', 'K', 'M'];

// Human-readable labels for UI / lore.
export const STAR_CLASS_LABELS: Record<StarClass, string> = {
  O: 'O — blue giant (scorching)',
  B: 'B — blue-white',
  A: 'A — white',
  F: 'F — yellow-white',
  G: 'G — yellow, Sun-like',
  K: 'K — orange dwarf',
  M: 'M — red dwarf (frigid)',
};

// Scale a Celsius temperature by the star class's insolation factor, in Kelvin.
// G-class (or unknown) is an exact identity — no ULP perturbation.
export const applyStarClass = (tempC: number, starClass: StarClass | undefined): number => {
  const f = STAR_TEMP_FACTOR[starClass ?? 'G'] ?? 1.0;
  if (f === 1) return tempC;
  return (tempC + 273.15) * f - 273.15;
};
