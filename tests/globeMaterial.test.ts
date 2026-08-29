import { describe, it, expect } from 'vitest';
import * as THREE from 'three';

import { createStyleTexture, createStyledGlobeMaterial } from '../utils/mapStyle/globeMaterial';

// CanvasTexture only stores the image; it touches no DOM, so a stub is enough.
const fakeCanvas = () => ({ width: 2048, height: 1024 }) as unknown as HTMLCanvasElement;

describe('createStyleTexture', () => {
  // THE bug this file exists for. three.js sets UNPACK_FLIP_Y_WEBGL from
  // texture.flipY, which DEFAULTS TO TRUE, so the canvas's top row uploads to
  // v = 1. Row 0 of an equirectangular canvas is the NORTH pole, and both
  // d3.geoEquirectangular and toLonLat put north at v = 0 — so the default
  // renders the world upside down, southern hemisphere over northern mesh.
  //
  // It survived three sessions of UV fixes across two sessions, because every
  // fix changed the UV formula (which was right) and never the upload. Proved
  // with `renderGlobePreview --marker --views=90,-90`: the canvas bottom came
  // out on the north pole.
  it('does NOT flip the texture on upload', () => {
    expect(createStyleTexture(fakeCanvas()).flipY).toBe(false);
  });

  it('repeats in u so the antimeridian blends, and clamps in v because the poles are the edge', () => {
    const tex = createStyleTexture(fakeCanvas());
    expect(tex.wrapS).toBe(THREE.RepeatWrapping);
    expect(tex.wrapT).toBe(THREE.ClampToEdgeWrapping);
  });

  // atan jumps a full turn across the antimeridian, so the screen-space
  // derivative explodes there and mip selection falls to the coarsest level —
  // a blurred band down the seam. Linear filtering needs no derivative.
  it('carries no mipmap chain', () => {
    const tex = createStyleTexture(fakeCanvas());
    expect(tex.generateMipmaps).toBe(false);
    expect(tex.minFilter).toBe(THREE.LinearFilter);
  });
});

describe('createStyledGlobeMaterial', () => {
  // Run the injection the way three.js does, against the include markers it
  // would supply, so a rename of any include is caught here rather than as a
  // silently untextured globe.
  const compile = () => {
    const mat = createStyledGlobeMaterial(createStyleTexture(fakeCanvas()));
    const shader = {
      vertexShader: '#include <common>\nvoid main(){\n#include <begin_vertex>\n}',
      fragmentShader: '#include <common>\nvoid main(){\n#include <map_fragment>\n}',
    };
    mat.onBeforeCompile(shader as never, null as never);
    return { mat, ...shader };
  };

  it('passes the undisplaced sphere direction to the fragment shader', () => {
    const { vertexShader, fragmentShader } = compile();
    expect(vertexShader).toContain('attribute vec3 sphereDir;');
    expect(vertexShader).toContain('vSphereDir = sphereDir;');
    expect(fragmentShader).toContain('varying vec3 vSphereDir;');
    // Interpolating unit vectors shortens them inside the triangle.
    expect(fragmentShader).toContain('normalize( vSphereDir )');
  });

  // Equirectangular, north at v = 0 — the same convention as
  // d3.geoEquirectangular and toLonLat. It pairs with flipY === false above:
  // change one without the other and the world goes upside down again.
  it('derives an equirectangular coordinate with north at v = 0', () => {
    const { fragmentShader } = compile();
    expect(fragmentShader).toContain('atan( sphereN.z, sphereN.x )');
    expect(fragmentShader).toContain('0.5 - asin( clamp( sphereN.y, -1.0, 1.0 ) )');
    expect(fragmentShader).toContain('texture2D( map, sphereUv )');
  });

  it('carries a program cache key, or three.js may reuse an uninjected program', () => {
    const { mat } = compile();
    expect(typeof mat.customProgramCacheKey).toBe('function');
    expect(mat.customProgramCacheKey()).toBeTruthy();
  });
});
