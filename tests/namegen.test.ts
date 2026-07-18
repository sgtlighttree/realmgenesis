import { describe, it, expect } from 'vitest';
import { RNG } from '../utils/rng';
import { createNameGenerator, NAME_STYLES, NameStyle } from '../utils/namegen';

// Draw a fixed-length sequence of names from a fresh generator seeded from a
// deterministic rng stream.
const sequence = (style: NameStyle, seed: string, count: number): string[] => {
  const rng = new RNG(seed);
  const gen = createNameGenerator(style, () => rng.next());
  const out: string[] = [];
  for (let i = 0; i < count; i++) {
    // Rotate the three semantic methods so all are exercised.
    if (i % 3 === 0) out.push(gen.faction());
    else if (i % 3 === 1) out.push(gen.province());
    else out.push(gen.town());
  }
  return out;
};

describe('namegen', () => {
  it('is deterministic for the same style + seed', () => {
    for (const style of NAME_STYLES) {
      const a = sequence(style, 'seed_alpha', 30);
      const b = sequence(style, 'seed_alpha', 30);
      expect(a).toEqual(b);
    }
  });

  it('produces different sequences for different seeds', () => {
    const a = sequence('fantasy', 'seed_one', 30);
    const b = sequence('fantasy', 'seed_two', 30);
    expect(a).not.toEqual(b);
  });

  it('produces mostly distinct output across styles', () => {
    // Same rng seed, different style — the wordlist should dominate the result.
    const styles = NAME_STYLES;
    for (let i = 0; i < styles.length; i++) {
      for (let j = i + 1; j < styles.length; j++) {
        const a = sequence(styles[i], 'shared_seed', 40);
        const b = sequence(styles[j], 'shared_seed', 40);
        const shared = a.filter(n => b.includes(n)).length;
        // Short 4-char names collide occasionally; require a clear majority
        // to be distinct rather than perfect separation.
        expect(shared).toBeLessThan(a.length * 0.4);
      }
    }
  });

  it('emits sane names: non-empty, capitalized, within length bounds', () => {
    for (const style of NAME_STYLES) {
      for (const name of sequence(style, 'sanity_seed', 60)) {
        expect(name.length).toBeGreaterThanOrEqual(3);
        expect(name.length).toBeLessThanOrEqual(14);
        expect(name).toMatch(/^[A-Z]/);
        expect(name.trim()).not.toBe('');
      }
    }
  });

  it('deduplicates names within a single generator instance', () => {
    const rng = new RNG('dedup_seed');
    const gen = createNameGenerator('fantasy', () => rng.next());
    const names = new Set<string>();
    for (let i = 0; i < 50; i++) names.add(gen.faction());
    expect(names.size).toBe(50);
  });
});
