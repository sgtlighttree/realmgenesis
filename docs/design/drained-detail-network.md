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
