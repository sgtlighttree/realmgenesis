import { describe, it, expect } from 'vitest';
import { BiomeType, Cell } from '../types';
import {
  seasonalDeclination,
  annualMeanLatTemp,
  seasonalTemperatureDelta,
  SEAWATER_FREEZE_C,
} from '../utils/seasons';
import { makeParams } from './helpers';

// A minimal cell at a given geometric latitude (radians). Only `center` matters
// to the seasonal helpers; the rest is filler to satisfy the Cell shape.
const cellAtLat = (latRad: number): Cell => ({
  id: 0,
  center: { x: 0, y: Math.sin(latRad), z: Math.cos(latRad) },
  vertices: [],
  neighbors: [],
  height: 0.6,
  plateId: 0,
  temperature: 0,
  moisture: 0.5,
  biome: BiomeType.OCEAN,
});

const SEASONS = [0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875];

describe('seasonal declination', () => {
  it('is zero at the equinox neutral point (season 0.5) and at season 0', () => {
    const p = makeParams({ axialTilt: 23.5 });
    expect(seasonalDeclination(0.5, p.axialTilt)).toBeCloseTo(0, 10);
    expect(seasonalDeclination(0, p.axialTilt)).toBeCloseTo(0, 10);
  });

  it('peaks at ±tilt at the solstices (season 0.25 / 0.75)', () => {
    const tiltDeg = 23.5;
    const tiltRad = (tiltDeg * Math.PI) / 180;
    expect(seasonalDeclination(0.25, tiltDeg)).toBeCloseTo(tiltRad, 10);
    expect(seasonalDeclination(0.75, tiltDeg)).toBeCloseTo(-tiltRad, 10);
  });
});

describe('seasonalTemperatureDelta', () => {
  it('is exactly 0 at the neutral season (0.5) everywhere', () => {
    const p = makeParams({ season: 0.5, axialTilt: 23.5 });
    for (const lat of [-1.2, -0.5, 0, 0.5, 1.2]) {
      expect(seasonalTemperatureDelta(cellAtLat(lat), p)).toBe(0);
    }
  });

  it('is identically 0 for all seasons when axialTilt is 0', () => {
    for (const s of SEASONS) {
      const p = makeParams({ season: s, axialTilt: 0 });
      for (const lat of [-1.0, 0, 0.8]) {
        expect(seasonalTemperatureDelta(cellAtLat(lat), p)).toBe(0);
      }
    }
  });

  it('is 0 at both equinox seasons (0 and 0.5) — anchored to the canonical world', () => {
    const p0 = makeParams({ axialTilt: 23.5 });
    for (const lat of [-1.1, -0.6, 0.3, 0.9]) {
      expect(seasonalTemperatureDelta(cellAtLat(lat), { ...p0, season: 0.5 })).toBe(0);
      // s=0 also has δ=0; the excursion is ~0 there (continuous, no detent).
      expect(seasonalTemperatureDelta(cellAtLat(lat), { ...p0, season: 0 })).toBeCloseTo(0, 6);
    }
  });

  it('is continuous through the neutral point (no "pop" leaving 0.5)', () => {
    const p0 = makeParams({ axialTilt: 23.5 });
    const lat = 0.8;
    const justOff = seasonalTemperatureDelta(cellAtLat(lat), { ...p0, season: 0.5 + 1e-4 });
    // A tiny step off neutral must produce a tiny excursion, not a finite jump.
    expect(Math.abs(justOff)).toBeLessThan(0.05);
  });

  it('flips sign between hemispheres at a solstice', () => {
    const p = makeParams({ season: 0.25, axialTilt: 23.5 });
    const north = seasonalTemperatureDelta(cellAtLat(0.7), p);
    const south = seasonalTemperatureDelta(cellAtLat(-0.7), p);
    // Northern-summer solstice: the subsolar point is north, so the northern
    // mid-latitude warms above its annual mean and the southern one cools.
    expect(north).toBeGreaterThan(0);
    expect(south).toBeLessThan(0);
    expect(Math.sign(north)).not.toBe(Math.sign(south));
  });
});

describe('D3 sea-ice edge moves with season', () => {
  // An ocean cell's stored temperature ≈ annual-mean latitude temp (no elevation
  // lapse underwater). Shown temp = stored + seasonal excursion. Somewhere in the
  // mid-high latitudes there must be a cell that is open water in its hemisphere's
  // summer but frozen (below the seawater freeze constant) in winter — i.e. the
  // sea-ice edge migrates. This is the physics D3's render overlay reads.
  it('has a latitude that freezes in winter but thaws in summer', () => {
    const p = makeParams({ axialTilt: 23.5, poleTemperature: -30, baseTemperature: 30 });
    const shown = (lat: number, s: number) =>
      annualMeanLatTemp(lat, p) + seasonalTemperatureDelta(cellAtLat(lat), { ...p, season: s });
    let flips = 0;
    for (let lat = 0.6; lat < 1.5; lat += 0.02) {
      const nSummer = shown(lat, 0.25); // subsolar north
      const nWinter = shown(lat, 0.75); // subsolar south
      if (nWinter < SEAWATER_FREEZE_C && nSummer >= SEAWATER_FREEZE_C) flips++;
    }
    expect(flips).toBeGreaterThan(0);
  });

  it('winter is colder than summer in the northern hemisphere', () => {
    const p = makeParams({ axialTilt: 23.5 });
    const lat = 1.1;
    const shown = (s: number) =>
      annualMeanLatTemp(lat, p) + seasonalTemperatureDelta(cellAtLat(lat), { ...p, season: s });
    expect(shown(0.75)).toBeLessThan(shown(0.25));
  });
});

describe('annualMeanLatTemp', () => {
  it('equals the plain latitude curve when tilt is 0', () => {
    const p = makeParams({ axialTilt: 0, baseTemperature: 30, poleTemperature: -30 });
    // Equator (lat 0): r=0 → base. Pole (lat π/2): r=1 → pole.
    expect(annualMeanLatTemp(0, p)).toBeCloseTo(30, 6);
    expect(annualMeanLatTemp(Math.PI / 2, p)).toBeCloseTo(-30, 6);
  });

  it('shifts with axialTilt (Jensen on the quadratic) so tilt stays a live param', () => {
    const equator = 0;
    const flat = makeParams({ axialTilt: 0, baseTemperature: 30, poleTemperature: -30 });
    const tilted = makeParams({ axialTilt: 60, baseTemperature: 30, poleTemperature: -30 });
    const dFlat = annualMeanLatTemp(equator, flat);
    const dTilt = annualMeanLatTemp(equator, tilted);
    // The equatorial annual mean must move by well more than the paramLiveness
    // 3-decimal temperature threshold.
    expect(Math.abs(dTilt - dFlat)).toBeGreaterThan(0.5);
  });
});
