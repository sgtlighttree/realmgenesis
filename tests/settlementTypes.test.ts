import { describe, it, expect } from 'vitest';
import { TIER_ORDER, TIER_MULTIPLIER } from '../utils/settlements';
import { DEFAULT_PARAMS } from '../utils/defaultParams';

describe('settlement tiers & params', () => {
  it('orders capital→hamlet with descending multipliers', () => {
    expect(TIER_ORDER).toEqual(['capital', 'city', 'town', 'village', 'hamlet']);
    const m = TIER_ORDER.map(t => TIER_MULTIPLIER[t]);
    expect(m).toEqual([...m].sort((a, b) => b - a));
    expect(TIER_MULTIPLIER.hamlet).toBe(1);
  });
  it('defaults settlementDensity to 1.0', () => {
    expect(DEFAULT_PARAMS.settlementDensity).toBe(1.0);
  });
});
