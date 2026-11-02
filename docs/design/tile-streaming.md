# Tile streaming, first increment (M1a)

Status: draft (s01, 2026-11-02). Goal: zoom from orbit down to about 10 m/px on the globe, with terrain detail appearing as you zoom. The ARCHITECTURE §4–5 constraints apply.

## Scope of M1a

- **In:** quadtree LOD over the cube-sphere; CPU tile generation within a per-frame time budget (web-safe, single-threaded); camera-relative rendering; procedural detail constrained by the committed macro field; crack hiding; LRU cache; a camera that can get close and tilt.
- **Out (M1b and later):** GPU compute tile generation (`rg-gpu`); integer-lattice noise in WGSL; texture detail (albedo stays the global equirect layer); 2D map tiles; sculpt deltas; `Job`-ified world generation.

## Render tiles vs data tiles

ARCHITECTURE pins **data** tiles at 256 cells / 257² vertices. Render tiles are separate: **64 cells / 65² vertices** plus a skirt ring. They're cheap to generate on the CPU (about 4.2k samples) and fine-grained for LOD.
- Render tile level L: cell ≈ (π/2·R) / (64·2^L). At R = 6371 km: L0 ≈ 156 km, L8 ≈ 610 m, **L14 ≈ 9.5 m (the floor)**.
- When data tiles arrive (sculpt deltas, M6), a render tile at level L samples the data tile at level L − 2 (since 256 = 4 × 64), so the two pyramids line up.

## Height of a sample

`h(dir) = macro(dir) + detail(dir, macro)`

- `macro`: bilinear sample of the committed `CubeField` (`tectonics.height_m`). Smooth, continuous across faces.
- `detail`: fBm + ridged noise from `rg-noise`, evaluated on the CPU in **f64** at the sample's unit direction scaled by frequency.
  - f64 on the CPU sidesteps the f32 lattice problem for M1a. The GPU port (M1b) must switch to integer lattice coordinates from the tile address.
  - Octaves start just above the macro field's Nyquist (wavelength ≈ 2 macro cells) and continue to about 2 × the floor cell.
  - Octave count is fixed per sample, not per tile level, so the same point gets the same height at every LOD. **This rule is what prevents popping beyond the usual geomorphing.** A tile at level L only *evaluates* octaves whose wavelength is at least 2 of its cells. Higher octaves are below its sampling rate anyway.
  - That means h does differ slightly between levels (missing octaves), so it's a band-limited version of the same function. Popping is bounded by the amplitude of the dropped octaves, which is tiny when the gain is 0.5.
- **Coast constraint:** the detail amplitude is scaled by `smoothstep(0, 150 m, |macro|)` and clamped so `sign(h) == sign(macro)`. The coastline is exactly the committed macro contour at every zoom; the vector layer stays honest.
- **Context:** amplitude grows with macro elevation (mountains rough, plains gentle) and is small on the ocean floor. Ridged octaves are weighted by elevation.

## Selection (per frame)

- Start from the 6 roots and recurse.
- **Cull:** a frustum test against the tile's bounding sphere, plus a horizon test against the planet sphere (radius = R + max height).
- **Split** if `projected_error_px > τ`, where `error = cell_size_m × k`, projected by distance to the camera (eye to nearest point on the bounding sphere). τ ≈ 2 px for 65² tiles.
- A tile is drawn only if it is resident. If a child isn't ready, the parent is drawn instead (and its generation is requested), so no holes appear while streaming.

## Generation budget

- Requests are queued and sorted by (level, distance).
- Each frame generates tiles until **6 ms (native) or 4 ms (web)** have elapsed. 65² samples × about 12 octaves in f64 takes about 1–2 ms per tile native; to be measured.
- Generated: positions relative to the **tile centre** (f32), normals from central differences over a 1-sample apron (so normals match across seams), the macro-height attribute for colouring, and skirt vertices dropped by `cell × 4` radially.
- **Cache:** `HashMap<TileId, GpuTile>`, LRU by last-used frame. Caps: 1200 tiles native, 600 web. Shared index buffer for the 65² grid + skirts.

## Camera-relative rendering

- Each tile has a centre `c` (f64, metres, ECEF scaled by R).
- Per draw, the vertex shader gets `offset = (c − eye)` as f32 (computed in f64 on the CPU) and a rotation-only view-projection.
- Per-tile offsets go in one storage buffer indexed by `instance_index`: draw with `first_instance = slot`, so there are no dynamic uniforms.
- The near plane is about altitude × 0.1 (clamped ≥ 0.5 m), with an infinite reversed-Z far plane.

## Camera

- Orbit camera with **distance in metres above the surface**. Min altitude is about 30 m above the terrain under the camera (sampled from the macro field + detail).
- **Tilt:** right-drag or two-finger vertical drag pitches the view toward the horizon, Google Earth-style. Tilt is capped by altitude: none at orbit, up to 75° near the ground.
- Zoom speed scales with altitude.

## Colour

Tiles sample the existing equirect layer texture by lat/lon (as now). It's blurry close up (about 20 km/px). Lighting from the detailed normals carries the near-field read for M1a. A detail-albedo pass (slope/altitude rock and snow, noise-varied biome tint) is M1c.

## Open questions for review

1. Is f64-CPU detail plus band-limited octaves acceptable as M1a, given the integer-lattice requirement is deferred to the GPU port?
2. Skirts vs stitching for cracks?
3. Is a storage buffer indexed by `instance_index` fine on WebGPU and Metal, or would push constants or dynamic uniforms be better?
4. Is the CPU per-frame budget OK on the web (single thread) until `rg-gpu`?

## Review outcome (advisor stand-in, s01): approved with changes, all adopted unless noted

1. **Exaggeration is a uniform**, not baked into vertices. Vertices store the unrelieved sphere position, true `h`, the parent-level `h` and the exaggeration-1 normal. The shader displaces by `max(h,0)·exag` and tilts normals by exag. Effective exaggeration = the user value at orbit, fading to about 1.5× below about 50 km. It's driven by camera distance from the target, not height above terrain (avoids feedback).
2. **Horizon cull:** the occluder is the sphere of radius R (oceans render at sea level). Tile bounds sit at R + h_max·exag. Occluded = fully inside the eye's tangent cone AND fully behind the horizon plane.
3. **Popping:** every vertex carries `h_parent` (the same point at the parent's octave set), and tiles morph `h_parent → h` over 250 ms after upload (by time, not distance, since tiles arrive late). Normals may seam between levels in M1a.
4. **Cost:** detail is skipped where macro ≤ 0 (rendered flat anyway). The native build uses worker threads (a pure `generate` fn); the web build generates in-frame within a budget (at least 1 tile per frame). Use `web-time` (`std::time::Instant` panics on wasm32). If a web dive takes > 5 s to settle, bring `rg-gpu` forward.
5. **Fallbacks:** draw the parent or all four children (culled children count as satisfied), never a mix. L0–L2 are pinned. Selection touches every visited ancestor, so the LRU never evicts one.
6. **Split metric:** a per-tile measured error = the max |h(level) − h(child level)| × exag + the chord sag, projected to physical px; τ ≈ 2 px. Budget about 600 drawn tiles.
7. **Camera:** Earth-style, with target lat/lon, distance (m), heading and tilt; tilt is capped by altitude. The minimum clearance is checked against displayed terrain. Picking hits the sphere for M1a (tile-height picking is later).
8. **Coast:** a smoothstep fade on |macro| plus a sign clamp is kept for M1a. *Known:* coasts are bilinear segments one macro cell long, and flat terraces can form at sea level. A coast SDF band + soft knee comes later.
9. **Bicubic macro sampling** (lighting creases at macro-cell scale): deferred.
10. **Memory:** about 36 B/vertex × 4.4k ≈ 160 KB/tile. Caps: 800 native, 500 web. Packing (octahedral normals) is deferred.
- **Open question 1 resolution:** f64 CPU noise is OK for M1a, tagged `DETAIL_GENERATOR_VERSION = 0` (not persistent). The GPU port must use a **3D** integer lattice (`floor(dir·f)` as i32 + an f32 offset); a 2D face-local lattice would seam along the 12 cube edges.
- **Question 3:** an instance-rate vertex buffer with `first_instance = slot` (works on Metal, WebGPU and DX12-emulated). No push constants on WebGPU.
