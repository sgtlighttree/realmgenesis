import { describe, it, expect } from 'vitest';
import { RNG, SimplexNoise } from '../utils/rng';

describe('RNG (mulberry32)', () => {
  it('is deterministic for the same seed', () => {
    const a = new RNG('hello');
    const b = new RNG('hello');
    for (let i = 0; i < 100; i++) {
      expect(a.next()).toBe(b.next());
    }
  });

  it('produces different sequences for different seeds', () => {
    const a = new RNG('seed-a');
    const b = new RNG('seed-b');
    const seqA = Array.from({ length: 10 }, () => a.next());
    const seqB = Array.from({ length: 10 }, () => b.next());
    expect(seqA).not.toEqual(seqB);
  });

  it('stays within [0, 1)', () => {
    const rng = new RNG('bounds');
    for (let i = 0; i < 10000; i++) {
      const v = rng.next();
      expect(v).toBeGreaterThanOrEqual(0);
      expect(v).toBeLessThan(1);
    }
  });

  it('range() respects min/max', () => {
    const rng = new RNG('range');
    for (let i = 0; i < 1000; i++) {
      const v = rng.range(-5, 5);
      expect(v).toBeGreaterThanOrEqual(-5);
      expect(v).toBeLessThan(5);
    }
  });
});

describe('SimplexNoise', () => {
  it('is deterministic per seed', () => {
    const a = new SimplexNoise(new RNG('noise'));
    const b = new SimplexNoise(new RNG('noise'));
    for (let i = 0; i < 50; i++) {
      const x = i * 0.13, y = i * 0.29, z = i * 0.41;
      expect(a.noise3D(x, y, z)).toBe(b.noise3D(x, y, z));
    }
  });

  it('outputs roughly within [-1, 1]', () => {
    const noise = new SimplexNoise(new RNG('bounds'));
    for (let i = 0; i < 5000; i++) {
      const v = noise.noise3D(i * 0.017, i * 0.031, i * 0.047);
      expect(Math.abs(v)).toBeLessThanOrEqual(1.05);
    }
  });
});
