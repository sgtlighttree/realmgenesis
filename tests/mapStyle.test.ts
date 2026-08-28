import { describe, it, expect } from 'vitest';

import { getMapStyle, MAP_STYLES } from '../utils/mapStyle';

describe('map style registry', () => {
  it('exposes default and parchment', () => {
    expect(Object.keys(MAP_STYLES).sort()).toEqual(['default', 'parchment']);
  });

  it('returns the default style for an unknown id', () => {
    // @ts-expect-error deliberately probing the runtime fallback
    expect(getMapStyle('nonsense').id).toBe('default');
  });

  it('gives every ViewMode a fill policy in every style', () => {
    const modes = [
      'biome', 'height', 'height_bw', 'temperature', 'moisture', 'plates',
      'political', 'population', 'province', 'satellite', 'culture', 'religion',
    ] as const;
    for (const style of Object.values(MAP_STYLES)) {
      for (const mode of modes) {
        expect(['bare', 'categorical', 'ramp']).toContain(style.fillPolicy(mode));
      }
    }
  });

  it('keeps every mode on the ramp policy in the default style', () => {
    // The default style is the pre-A3 look: every mode paints its own fill.
    expect(getMapStyle('default').fillPolicy('political')).toBe('ramp');
    expect(getMapStyle('default').fillPolicy('satellite')).toBe('ramp');
  });
});
