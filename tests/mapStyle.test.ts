import { describe, it, expect } from 'vitest';

import { getMapStyle, MAP_STYLES } from '../utils/mapStyle';

describe('map style registry', () => {
  it('exposes every registered style', () => {
    // An exact inventory on purpose: registering a style is a deliberate act,
    // and this is the line that makes an accidental one visible. Adding a style
    // means updating this list — that is the test working, not the test being
    // in the way.
    expect(Object.keys(MAP_STYLES).sort()).toEqual(['blueprint', 'boardgame', 'default', 'inkwash', 'parchment']);
  });

  it('inks every line overlay in an OPAQUE colour', () => {
    // The SVG export consumes these and SVG 1.1 has no rgba() syntax. The
    // default style is exempt: it is the grandfathered pre-A3 look and its
    // rgba() values reproduce it exactly. Every style added since must be
    // opaque, and `contour`/`currents` may be null to keep their own ramps.
    for (const style of Object.values(MAP_STYLES)) {
      if (style.id === 'default') continue;
      for (const [key, value] of Object.entries(style.overlayInk)) {
        if (typeof value !== 'string') continue;
        expect(`${style.id}.${key}=${value}`).not.toContain('rgba(');
      }
    }
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
