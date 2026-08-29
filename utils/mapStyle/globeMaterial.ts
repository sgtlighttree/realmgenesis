import * as THREE from 'three';

/**
 * The styled globe's material: unlit, and it derives its own equirectangular
 * texture coordinate PER FRAGMENT from an interpolated sphere direction.
 *
 * **Why not a `uv` attribute.** An equirectangular u IS a longitude, and
 * longitude is not linear across a triangle on a sphere. Per-vertex UVs are
 * interpolated linearly, so every triangle sampled a slightly wrong strip of
 * texture, and cells near a pole — where a small triangle legitimately spans a
 * wide longitude range — sampled an enormous one. Wrapping each vertex to the u
 * nearest its cell centre took the worst measured triangle u-span from 0.633 to
 * 0.000, but only by collapsing polar cells onto a single meridian, which traded
 * dark polar streaks for a spiral rosette at the pole. Trading one artifact for
 * another is all per-vertex UVs can do here.
 *
 * Computing the coordinate in the fragment shader removes the interpolation
 * error outright — every fragment gets its own exact spherical coordinate — so
 * there is no seam wrap and no polar collapse to get wrong.
 *
 * The direction comes from a `sphereDir` attribute holding the UNDISPLACED unit
 * sphere position, not from `position`: position carries per-cell height, so
 * normalizing it would make neighbouring cells at different heights sample
 * different content either side of their shared edge — a UV jog on every cell
 * boundary.
 *
 * The pole still pinches, because equirectangular content genuinely converges
 * there. That is a texture-content problem, handled in the bake (the ocean hatch
 * fades out at high latitude), not a geometry one.
 *
 * UNLIT, like political mode, and for the same reason: the baked texture already
 * contains the hillshade pass, so a lit material shades it a second time and
 * warm paper renders as dark grey-brown. Unlit also makes the globe match the 2D
 * map and the exports exactly, which is the point of putting the style here.
 * Relief still reads — the mesh keeps its displacement, and the baked shading
 * travels with the texture.
 *
 * The caller owns disposal.
 */
export const createStyledGlobeMaterial = (map: THREE.Texture): THREE.MeshBasicMaterial => {
  const mat = new THREE.MeshBasicMaterial({ map, toneMapped: false, side: THREE.FrontSide });

  mat.onBeforeCompile = (shader) => {
    shader.vertexShader = shader.vertexShader
      .replace(
        '#include <common>',
        '#include <common>\nattribute vec3 sphereDir;\nvarying vec3 vSphereDir;',
      )
      .replace(
        '#include <begin_vertex>',
        '#include <begin_vertex>\n\tvSphereDir = sphereDir;',
      );

    shader.fragmentShader = shader.fragmentShader
      .replace('#include <common>', '#include <common>\nvarying vec3 vSphereDir;')
      .replace('#include <map_fragment>', [
        // Interpolating three unit vectors gives a shorter vector inside the
        // triangle, so it has to be renormalized before the angles are read.
        'vec3 sphereN = normalize( vSphereDir );',
        // 0.15915494 = 1/(2*PI); 0.31830989 = 1/PI. Same mapping as
        // d3.geoEquirectangular uses for the bake, so texel and cell agree.
        'vec2 sphereUv = vec2(',
        '  atan( sphereN.z, sphereN.x ) * 0.15915494 + 0.5,',
        '  0.5 - asin( clamp( sphereN.y, -1.0, 1.0 ) ) * 0.31830989',
        ');',
        'diffuseColor *= texture2D( map, sphereUv );',
      ].join('\n'));
  };

  // Without this, three.js can hand back a cached program compiled from the
  // same material feature set but WITHOUT the injection above.
  mat.customProgramCacheKey = () => 'globe-sphere-uv';
  return mat;
};

/**
 * Prepare a baked style canvas for use as the globe texture.
 *
 * NO MIPMAPS, deliberately. The fragment shader computes u with `atan`, which
 * jumps by a full turn across the antimeridian, so the screen-space derivative
 * explodes on that one line and mip selection falls to the coarsest level — a
 * blurred band down the seam. Linear filtering needs no derivative at all. The
 * texture is 2048 wide and the globe shows roughly half a turn across the
 * viewport, so this sits near 1:1, not deep in minification. (Anisotropy goes
 * with the mipmap chain; it has nothing left to filter.)
 *
 * RepeatWrapping on u so the antimeridian blends across the seam instead of
 * clamping: the texel at u = 0 is the neighbour of the texel at u = 1.
 * ClampToEdge on v — the default — is correct, because the poles ARE the edge.
 *
 * **`flipY = false` is the whole ball game.** three.js sets
 * `UNPACK_FLIP_Y_WEBGL` from `texture.flipY`, which defaults to TRUE, so the
 * canvas's top row is uploaded to v = 1. Row 0 of an equirectangular canvas is
 * the NORTH pole, and both `d3.geoEquirectangular` and `toLonLat` put north at
 * v = 0 — so the default flip renders the world UPSIDE DOWN, with the southern
 * hemisphere painted over the northern mesh.
 *
 * It had been that way since `buildGlobeUVs`, whose `vOf` used the same
 * (correct) formula, which is why no UV fix ever touched it: coastlines sat on
 * the wrong hemisphere, mountain glyphs landed on ocean, and the "polar swirl"
 * chased across three sessions was the SOUTH pole's content drawn on the north.
 * Proved with `renderGlobePreview --marker --views=90,-90`, which paints the
 * canvas top red and bottom blue: blue came out on the north pole.
 *
 * Corrected here, at the upload, rather than by negating v in the shader. The
 * shader's `v = 0.5 - asin(y)/PI` agrees with the projection it mirrors; a
 * negated copy would disagree with it, and the next person re-derives the lot.
 */
export const createStyleTexture = (canvas: HTMLCanvasElement): THREE.CanvasTexture => {
  const tex = new THREE.CanvasTexture(canvas);
  tex.flipY = false;
  tex.wrapS = THREE.RepeatWrapping;
  tex.colorSpace = THREE.SRGBColorSpace;
  tex.generateMipmaps = false;
  tex.minFilter = THREE.LinearFilter;
  return tex;
};
