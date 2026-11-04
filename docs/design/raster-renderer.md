# R2: raster-displacement globe (design, s03, 2026-11-04)

Status: **proposed**, fresh-context review pending. Parent: `docs/design/v2-direction-s03.md` (R2).

## Goal

Heights come from a **baked raster on the GPU**, not from per-vertex CPU sampling:
- tiles become geometry only (sphere positions + face uv);
- the vertex shader displaces by the raster;
- the fragment shader lights from a baked gradient raster, so normals are per pixel at any zoom.

That removes the CPU height path (detail network, gullies, per-tile sampling), its latency and its UI stalls. It also gives later work one place to write heights: R4's regional bakes override the raster inside their footprint, and R5's exports read the same raster.

**Expectation (for Matt):** at the same data, R2 looks about like today's macro-only Ultra. The macro is the data. R2 buys smoothness at every zoom, no stalls, and the plumbing for regions and exports. *Detail* comes from R1 (routing) and R4 (bakes).

## What stays

Built and working, kept as is:
- the quadtree selection, cache/LRU, culling and horizon test (`rg-render/src/terrain.rs`);
- camera-relative instances (f64 centre − eye on the CPU);
- skirts and the exaggeration uniform, with its orbit → near-ground fade and "keep relief";
- the equirect layer texture for colour (plus its `height_tex`, used by the 2D map);
- the camera and tilt rules.

## Data: `DisplayRaster` (CPU, then GPU)

- **Layout:** 6 faces × (S+1)² samples, **corner-aligned like `CubeField`** (face-edge samples duplicated). Bilinear and bicubic then agree at face edges without cross-face lookups, except for the 1-texel bicubic footprint (see seams).
- **S = 2048** for every preset: ~4.9 km per texel at Earth radius. One code path and predictable memory. For presets coarser than 2048 it's a bicubic upsample.
- **Height:** baked from `TerrainSampler::macro_height` in macro-only mode: Catmull-Rom, with the coast sign taken from the drained surface (s02 fix). This is a *display* cache derived from committed fields; it isn't committed, and it isn't stored in the world file. Region bakes (R4) will be committed and stored.
- **Gradient:** the tangent-plane gradient of h (m/m) in a local (east, north) frame. It's baked by 3D central differences on the sphere (`macro_height(dir ± ε·east)`, etc.), so it crosses face edges correctly; that's the cross-face apron the reviewer asked for. The frame is singular exactly at the poles: there, fall back to a fixed tangent (`normalize_or`); it only affects a pixel or so.
- **GPU formats:**
  - height `R32Float` in a `texture_2d_array` (6 layers), read with `textureLoad` (needs no filterable-float feature);
  - gradient `Rg16Float` array, hardware-filtered (filterable everywhere).
  - Memory: ~100 MB + ~100 MB at S = 2048. The CPU copy is kept (it serves CPU height queries: camera clearance, hover).
  - Later options if memory bites: `R16Uint` heights (0.3 m steps over −11…+9 km), or S = 1024 on Draft/Normal.
- **Bake cost:** 25 M samples. Run it on the generation thread with rayon, as a stage after climate ("display raster" in the console). Estimate 2–5 s on the M1; to be measured. It reruns on load and on relief/sea-level changes that alter the macro.

## Tiles (CPU, trivial)

- The `TileVertex` becomes:
  - `pos` (sea-level sphere point minus the tile centre, f32, as today);
  - `uv` (face uv in raster texel space);
  - `skirt`.
- Per instance: today's `offset`/`morph`/`skirt_m`/`cell_m`, plus the **face layer**.
- `TileGenerator::generate` shrinks to positions + uv (4.4k `uv_to_dir` calls, well under 1 ms), so tiles can generate synchronously in `select` order. Worker threads become optional; keep the queue for now.
- **Error and split:**
  - `child_error_m` comes from the raster: the max |bicubic − bilinear-between-vertices| at cell centres, computed once per tile on the CPU from the `DisplayRaster` copy.
  - **No splitting below `cell ≤ texel / 2`**: past that the mesh only re-samples a smooth spline.
  - Today's τ-px rule stays.

## Shaders (`planet.wgsl`, globe path)

- **Height:** Catmull-Rom over 4×4 `textureLoad`s at the vertex uv (bilinear would crease the geometry at texel scale, the s02 "triangular" look). Indices are clamped within the face. The duplicated edge samples keep same-face edges exact; across a face edge, the one-sided stencil leaves a sub-texel kink, which is hidden by the gradient normals.
- **Geomorph:** as today, morph only the height toward the parent mesh. The parent height at an odd vertex is the mean of the heights at its two parent-grid neighbours' uvs (parity from `vertex_index` on the 65-wide grid; neighbour uv = uv ± one tile cell). That's up to two extra bicubic samples per vertex: about 48 loads, ~2.5 M vertices, fine on the M1 (to be measured).
- **Normal (fragment):** sample the gradient texture (hardware bilinear) at the interpolated uv. Normal = normalize(up − exag·(gₑ·east + gₙ·north)), with east/north from `up` and the same pole fallback. The vertex normal attribute goes away.
- **Rivers:** the existing layer-alpha river overlay. Its zoom fade (built so the detail network's carved valleys took over) is **removed**, so rivers stay visible at all zooms. A vector river layer is a separate later step (it also feeds the SVG export).

## Zoom floor (coverage-based)

- min altitude = max(terrain clearance, texel_m × (viewport_h_px / (2·tan(fov/2))) / M), where M = the allowed magnification (screen px per data texel).
- **texel_m** = the *data* texel: max(raster texel, macro cell). So Normal (31 km cells) floors higher than Ultra.
- **M defaults to 8** (reviewer's rule; ~800 km over a 5 km texel at 1080 px, 45° FOV). It's a slider in the view settings so Matt can push it; R4 regions lower texel_m locally.

## Removed from the app (code stays in git history)

- The "Streamed detail (experimental)" toggle, the `DetailNet` build in `start_tiles`, network root tiles, and the per-vertex sampler call in tiles.
- `rg-terrain`'s detail code (`detail_net.rs`, `erosion_filter.rs`, gullies) stays for the probes until R1 lands. Deleting it is a separate commit after Matt signs off.

## Steps and gates

1. **`DisplayRaster` + bake** (`rg-terrain`). Unit tests:
   - equals `macro_height` at nodes;
   - identical across face-edge twins;
   - gradient vs an analytic field (`d.z · 3000`) within 2%.
2. **Render path:** texture arrays, new vertex/instance formats, WGSL height/geomorph/normal.
   - Gate: same-level neighbours share edge heights (port `shared_edges_match_between_neighbours` to uv sampling), plus `cargo test`/clippy.
3. **App:** drop the detail toggle/network, sync tile generation, split cap, zoom floor + slider, river fade removed, bake stage in the console.
4. **Measure and look** (seed 1 Ultra, Matt's s02 viewpoints):
   - fps orbit/close, idle CPU, RSS (today: 0% idle, 625 MB at Normal);
   - bake time;
   - screenshots at 3 altitudes, compared with s02.
   - Matt judges.

## Risks

- The vertex-shader cost of 16–48 loads: if fps drops, cache the parent heights in a per-tile buffer (CPU-sampled, 17 KB per tile) instead.
- A gradient/height mismatch at face edges (lighting seam): the 3D-difference bake should prevent it; check the 12 edges visually.
- Pole frame: check the poles visually.
- Display S = 2048 vs Ultra's 1024 macro: double bicubic smoothing, harmless.
