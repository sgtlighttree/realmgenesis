import { describe, it, expect } from 'vitest';
import { RNG, SimplexNoise } from '../utils/rng';
import { seedCrustField } from '../utils/crust';
import { generateFibonacciSphere, chordDistance } from '../utils/spherical';
import { makeParams } from './helpers';

describe('V3 crust field', () => {
  it('produces continental and oceanic cells', () => {
    const points = generateFibonacciSphere(300, new RNG('test'), 0.5);
    const simplex = new SimplexNoise(new RNG('test'));
    const rng = new RNG('test_crust');
    const crust = seedCrustField(points, makeParams(), simplex, rng);
    expect(crust.crustTypes.length).toBe(300);
    expect(crust.crustThickness.length).toBe(300);
    const continentCount = Array.from(crust.crustTypes).filter(t => t === 1).length;
    expect(continentCount).toBeGreaterThan(0);
    expect(continentCount).toBeLessThan(300);
  });

  it('landStyle changes crust density', () => {
    const points = generateFibonacciSphere(300, new RNG('test'), 0.5);
    const simplex = new SimplexNoise(new RNG('test'));

    const continents = seedCrustField(
      points, makeParams({ landStyle: 'Continents' }), simplex, new RNG('c1')
    );
    const islands = seedCrustField(
      points, makeParams({ landStyle: 'Islands' }), simplex, new RNG('c2')
    );
    const continentCount = Array.from(continents.crustTypes).filter(t => t === 1).length;
    const islandCount = Array.from(islands.crustTypes).filter(t => t === 1).length;
    expect(continentCount).toBeGreaterThan(islandCount);
  });

  it('same seed produces same crust field', () => {
    const points = generateFibonacciSphere(300, new RNG('test'), 0.5);
    const simplex = new SimplexNoise(new RNG('test'));
    const a = seedCrustField(points, makeParams(), simplex, new RNG('crust_test'));
    const b = seedCrustField(points, makeParams(), simplex, new RNG('crust_test'));
    expect(Array.from(a.crustTypes)).toEqual(Array.from(b.crustTypes));
    expect(Array.from(a.crustThickness)).toEqual(Array.from(b.crustThickness));
  });

  it('oceanic crust is thinner than continental', () => {
    const points = generateFibonacciSphere(300, new RNG('test'), 0.5);
    const simplex = new SimplexNoise(new RNG('test'));
    const crust = seedCrustField(points, makeParams(), simplex, new RNG('crust_test'));

    let avgContinental = 0, avgOceanic = 0;
    let contCount = 0, oceanCount = 0;
    for (let i = 0; i < crust.crustTypes.length; i++) {
      if (crust.crustTypes[i] === 1) {
        avgContinental += crust.crustThickness[i];
        contCount++;
      } else {
        avgOceanic += crust.crustThickness[i];
        oceanCount++;
      }
    }
    avgContinental /= contCount;
    avgOceanic /= oceanCount;
    expect(avgContinental).toBeGreaterThan(avgOceanic);
  });
});

describe('V3 spherical math', () => {
  it('chord distance is zero for same point', () => {
    const p = { x: 1, y: 0, z: 0 };
    expect(chordDistance(p, p)).toBe(0);
  });

  it('chord distance between opposite poles is 2 on unit sphere', () => {
    const north = { x: 0, y: 1, z: 0 };
    const south = { x: 0, y: -1, z: 0 };
    expect(chordDistance(north, south)).toBeCloseTo(2, 5);
  });
});