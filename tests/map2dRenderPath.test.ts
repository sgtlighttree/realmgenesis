import { describe, it, expect } from 'vitest';
import { shouldUseGLFill } from '../utils/map2dRenderPath';

describe('shouldUseGLFill', () => {
  const base = { styleId: 'default', projectionType: 'mercator', webglAvailable: true, contextLost: false };
  it('uses GL for default style on a non-dymaxion projection with webgl', () => {
    expect(shouldUseGLFill(base)).toBe(true);
  });
  it('falls back for parchment', () => {
    expect(shouldUseGLFill({ ...base, styleId: 'parchment' })).toBe(false);
  });
  it('falls back for dymaxion', () => {
    expect(shouldUseGLFill({ ...base, projectionType: 'dymaxion' })).toBe(false);
  });
  it('falls back when webgl unavailable', () => {
    expect(shouldUseGLFill({ ...base, webglAvailable: false })).toBe(false);
  });
  it('falls back after context loss', () => {
    expect(shouldUseGLFill({ ...base, contextLost: true })).toBe(false);
  });
});
