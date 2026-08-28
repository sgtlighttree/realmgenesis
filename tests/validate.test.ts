import { describe, it, expect } from 'vitest';
import { validateWorldParams, validateCivData, withParamDefaults } from '../utils/export';
import { makeParams } from './helpers';

// makeParams uses 300 points for engine-test speed, below the validator's
// UI-aligned floor of 2000 — bump to a realistic count for validation tests.
const validParams = (overrides = {}) => makeParams({ points: 5000, ...overrides });

describe('validateWorldParams', () => {
  it('accepts the default parameter set', () => {
    expect(validateWorldParams(validParams())).toBe(true);
  });

  it('accepts every land style preset', () => {
    for (const landStyle of ['Continents', 'Archipelago', 'Islands', 'Pangea', 'Custom'] as const) {
      expect(validateWorldParams(validParams({ landStyle }))).toBe(true);
    }
  });

  it('rejects out-of-bounds numeric values', () => {
    expect(validateWorldParams(validParams({ points: 999999999 }))).toBe(false);
    expect(validateWorldParams(validParams({ seaLevel: 2 }))).toBe(false);
    expect(validateWorldParams(validParams({ numFactions: 100 }))).toBe(false);
    expect(validateWorldParams(validParams({ erosionIterations: -1 }))).toBe(false);
  });

  it('rejects NaN/Infinity and wrong types', () => {
    expect(validateWorldParams(validParams({ seaLevel: NaN }))).toBe(false);
    expect(validateWorldParams(validParams({ roughness: Infinity }))).toBe(false);
    expect(validateWorldParams({ ...validParams(), seed: 5 })).toBe(false);
    expect(validateWorldParams({ ...validParams(), loreLevel: 7 })).toBe(false);
  });

  it('rejects non-objects', () => {
    expect(validateWorldParams(null)).toBe(false);
    expect(validateWorldParams([])).toBe(false);
    expect(validateWorldParams('params')).toBe(false);
  });
});

describe('physicalClimate param', () => {
  it('defaults a missing physicalClimate to true (old saves get grounded climate)', () => {
    const p = makeParams();
    delete (p as unknown as Record<string, unknown>).physicalClimate;
    expect(withParamDefaults(p).physicalClimate).toBe(true);
  });
  it('preserves an explicit false', () => {
    const p = makeParams({ physicalClimate: false });
    expect(withParamDefaults(p).physicalClimate).toBe(false);
  });
  it('rejects a non-boolean physicalClimate', () => {
    // Use validParams() (points: 5000, clears the [2000, 200000] floor) so
    // physicalClimate is the ONLY reason validateWorldParams can reject —
    // makeParams()'s default 300 points would reject on its own.
    const bad = { ...validParams(), physicalClimate: 'yes' } as unknown;
    expect(validateWorldParams(bad)).toBe(false);
  });
  it('accepts a boolean physicalClimate', () => {
    // makeParams() defaults to 300 points (engine-test speed), below the
    // validator's UI-aligned 2000 floor — use validParams() so this
    // assertion isolates physicalClimate, not an unrelated points failure.
    expect(validateWorldParams(validParams({ physicalClimate: true }))).toBe(true);
  });
});

describe('validateCivData', () => {
  const validCivData = {
    factions: [
      {
        id: 0,
        name: 'Faction 1',
        color: '#e53935',
        capitalId: 12,
        totalPopulation: 5000,
        provinces: [
          { id: 0, name: 'Capital Region', totalPopulation: 5000, towns: [{ name: 'Capital City', cellId: 12, population: 5000, isCapital: true }] },
        ],
      },
    ],
  };

  it('accepts well-formed civData', () => {
    expect(validateCivData(validCivData)).toBe(true);
  });

  it('accepts an empty faction list', () => {
    expect(validateCivData({ factions: [] })).toBe(true);
  });

  it('rejects structural corruption', () => {
    expect(validateCivData(null)).toBe(false);
    expect(validateCivData({})).toBe(false);
    expect(validateCivData({ factions: 'nope' })).toBe(false);
    expect(validateCivData({ factions: [{ id: 'zero', name: 'x', color: '#fff', capitalId: 1, provinces: [] }] })).toBe(false);
    expect(validateCivData({ factions: [{ id: 0, name: 'x', color: '#fff', capitalId: 1, provinces: [{ id: 0, name: 'p' }] }] })).toBe(false);
  });
});
