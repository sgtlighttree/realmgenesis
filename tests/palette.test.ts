import { describe, it, expect } from 'vitest';
import { CULTURE_COLORS, FOLK_COLORS } from '../utils/palette';
import { darkenForFolk } from '../utils/colors';

describe('palette', () => {
  it('FOLK_COLORS matches the THREE-based darkenForFolk exactly', () => {
    expect(FOLK_COLORS).toEqual(CULTURE_COLORS.map(darkenForFolk));
  });

  it('FOLK_COLORS is index-aligned to CULTURE_COLORS', () => {
    expect(FOLK_COLORS.length).toBe(CULTURE_COLORS.length);
  });
});
