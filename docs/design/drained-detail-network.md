# Drained detail network (s01 design, v1)

**Why:** stateless noise detail can be steep or drained, not both (s01 measurements: hills put 13–16% of mountain land in closed depressions; gullies alone reach 14–17° mean slopes only at 768 macro res and still dam 2–5%; at Draft res the band from the macro cell down to 4 km is empty and mountains average 6–8°, against Earth's 15–30°). Matt also finds the gully filter "a bit swirly". The fix is to make the sub-grid detail a **channel network**, so drainage holds by construction and steepness comes from hillslope physics, not from grid resolution.

**Lineage:** Dendry (Gaillard, Guérin, Génevaux, Galin 2019, "Dendry: a procedural model for dendritic patterns", I3D) builds a point-evaluated hierarchy of jittered control points connected into branches; Génevaux et al. 2013 ("Terrain generation using procedural models based on hydrology", SIGGRAPH) build terrain from a river network plus hillslope profiles. This design combines the two and ties both to the macro stream-power solve. Original code; no code from either paper is used.

## The surface

For a land point p:

```
h(p) = min( U(p),  min over channel segments s near p of  W_s(p) )
W_s(p) = z_s(π_s(p)) + H(d_s(p))
```

- `π_s(p)`: the nearest point on segment s; `d_s(p)` the distance to it; `z_s` the channel elevation, linear along s, **strictly decreasing downstream**.
- `H(d) = S_h · w · tanh(d / w)`: the hillslope profile. It is strictly increasing, so a wedge descends towards its channel everywhere, with slope `S_h` at the channel and convex hilltops (diffusion-like) at about `w`.
- `U(p) = M(p) + R(p)`: a cap, the drained macro interpolant plus a ridge allowance, so the network never builds relief the macro doesn't have. When it binds, it gives flat-topped divides of macro height.

### Drainage guarantee

1. A local minimum of `min(f_i)` is a local minimum of the active `f_i`. (Where `f_a` is active and has a descent direction, the min decreases along it too.)
2. A wedge `W_s` has only one local minimum: the downstream end of s. That's because `H` strictly increases with d and `z_s` strictly decreases downstream.
3. The downstream end of a level-k segment lies on a level-(k−1) segment, with the same z. So the parent wedge is ≤ there and keeps descending.
4. Level 0 is the macro network (node → receiver), and it reaches the sea.

So the network part has no closed depressions. `U` adds at most the macro interpolant's residual (0.05–3% today, and only where U binds). Flats can only occur where `H` saturates numerically; `tanh` keeps it strict.

## The hierarchy

- **Level 0:** macro nodes and their receiver edges (from `Hydrology::flow_dir`), z = macro node height, k_sn from the macro edge: `k_sn = S · A^θ` (stream-power steady state, θ = m = 0.45).
- **Level k ≥ 1:** one jittered point per cell of a lattice with spacing `ℓ_k = ℓ_0 / 3^k` (ℓ_0 = macro cell), down to about 100 m (4–5 levels).
  - Cells are addressed per cube face, in (u, v); a neighbour cell past a face edge is resolved through its direction to the owning face. So every point has a single identity (face, level, i, j), and nothing seams.
  - Each point q **attaches** to the nearest point on the level-(k−1) segments within `2 ℓ_{k−1}`. Its channel runs straight from q to there (v1; curvature later).
- **Channel elevation:** `z(q) = z_attach + max(ε·L, M(q) − D_k(q) − z_attach)`, where L is the segment length.
  - So channel heads follow the macro surface minus a valley depth `D_k`, yet always rise upstream (ε = 1e-3, strict).
  - `D_k = min(S_h · ℓ_k / 2, channel relief from k_sn)`: hillslope-limited in steep terrain, channel-limited on plains.
- **Hillslope gradient `S_h`** from the erosion rate, saturating at the threshold angle (Roering et al. 2007; Montgomery & Brandon 2002): `S_h = S_c · (1 − exp(−E / E_0))`, with `E = K · k_sn` (n = 1), `S_c = tan 33°`. On plains (low k_sn) hillslopes stay gentle. In fast-uplift belts they reach 25–35°.

## Evaluation and cost

- **Point-evaluated and deterministic:** level points come from a hash of (seed, face, level, i, j). A per-thread memo (cell → point, attach, z) avoids recomputing parents, and it's only a cache. Per sample: about 16 macro segments plus about 25 points per level × 5 levels. Target: no slower than today's gully filter (about 40 ms per deep tile on the CPU).
- **LOD:** a level is used only if `ℓ_k ≥ 2 ×` the tile cell (band-limited like the old octaves). Parent-level heights for the CDLOD morph come from dropping the finest level.
- `river_corridor` is replaced: the level-0 wedges *are* the macro river valleys.
- The gully filter, the hills and their `RG_*` knobs are retired once this passes the scorecard. The MPL port stays in its file, unused, until removed deliberately.

## Validation (`probe_realism`, final surface)

- Depressions < 1% everywhere (mountain and lowland boxes), at Draft, Normal and 768.
- Mountain 90 m slopes: mean 15–30°, p99 ≤ ~50°. Plains < 4°. `terrain_realism` passes unchanged.
- Hypsometry of the final surface within the land curve tolerances. Mean land elevation may shift by `mean(H)`; if it does, calibrate the remap target, not the detail.
- `river_map` with `RG_MAP_FULL=1`: dendritic valleys, no swirls, no cell-edge lines, no seams at face edges.

## Open questions (for the reviewer)

1. Does `U` binding break drainage in practice? It's the macro interpolant, with a 0.05–3% residual.
2. Is `D_k = S_h ℓ_k / 2` the right valley-depth scaling, or should channel relief follow `k_sn ∫ A^−θ`?
3. Straight child segments will look like a polyline at the finest level; Dendry uses midpoint displacement. Is that needed in v1?

## Review outcome (advisor stand-in, s01): go-with-changes → v1.1

The lemmas hold. The surface as first specified didn't, so these are the changes:

1. **Reach wall:** `H(d) = (S_h − S_t)·w·tanh(g/w) + S_t·g + P·max(0, d − ρ_k)²` with `g = √(d² + f²) − f` (rounded floor of width f).
   - The wall P makes a wedge exceed any land height beyond `r_k`. So a candidate set of "all segments within `r_k`" is exact: it depends only on p, so tiles agree and there are no seams.
   - A parent passes through its child's junction, so it's always in reach there, and there's no pit at the junction.
   - `S_h, S_t, w, f` are fixed per segment, not sampled at p (lemma 2 needs that).
2. **Coast:**
   - Level-0 edges end at z = 0 at sea.
   - The cap `U = T(p) + R(T(p))`, where T is the **triangle-split** macro (proven drained) and R(0) = 0. U is then a strictly increasing function of a drained surface, so it's drained.
   - Land h ∈ (0, U], so the coastline is exact without `coast_rule`.
3. **Elevation, junction up:** a child channel rises from its junction with gradient `S_k = min(S_cap, k_sn·A_k^−θ)`, where `A_k ≈ ℓ_k²` is the area of a level-k cell and `k_sn` is inherited from the parent (macro edges: `k_sn = S·A^θ`).
   - The head is capped at `U(q) − hollow_k`, always ≥ `z_a + εL` (ε = 1e-3).
   - So valleys are shallow at the head and deepen downstream, as real ones do, and relief comes from the stream-power law at every level.
   - Hillslopes: `S_h = S_c·(1 − e^(−E/E_0))` with `E = K·k_sn`, `S_c = tan 33°`, `E_0 = 0.3 mm/yr`. `S_t = 0.3·S_h`. `w = ℓ_k/4`, `f = ℓ_k/8` (level 0: hydraulic width).
4. **Attach:**
   - To the nearest point on the union of all coarser levels, shifted downstream (about 60° junctions).
   - Skip targets the head would have to climb to (U at the midpoint above U at the head).
   - The attach depends only on q's own cell window (memo-safe).
5. **Evaluation per patch:** gather the segments within reach of a patch (a tile, or a probe box) coarse to fine into flat arrays, sorted by `z_min`. Per vertex, stop once `z_min ≥` the current best.
6. **LOD:** a level fades in by blending its profile toward a steep ramp from the junction. That's strictly monotone at every blend, so it stays drained, with no pops.
   - Probes and readouts use unmorphed heights.
   - Distances in 3D (chords), never in (u, v).
7. **Later:** midpoint-displaced (curved) segments, stored receivers (schema bump), and a strictly increasing remap of the final surface's hypsometry with g(0) = 0.
