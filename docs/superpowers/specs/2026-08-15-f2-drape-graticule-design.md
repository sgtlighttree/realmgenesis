# F2 — Draped graticule (grid follows terrain relief)

**Status:** design · **Date:** 2026-08-15 · **Branch:** `f2-drape-graticule`

## Problem

Session 17 shipped smooth-globe mode: turning on the Lat/Long grid flattens the
globe so the graticule can't parallax against relief. That works, but it forces
the user off the raised (3D-relief) look to read a grid. This is the deferred
"drape" tier: make the graticule **follow the terrain surface** on the *raised*
globe so it has zero parallax AND keeps relief — the Google-Earth read.

The residual parallax (S17 finding) was because the graticule sits at ONE radius
while terrain spans r=1.0..1.05; a grid line and a mountain at the same lat/lon
project to different pixels under perspective. The fix: project each graticule
sample at the **terrain radius at that lat/lon**, not a fixed radius.

## Design

Scope is the **3D graticule tenant only** (`components/overlays/tenants.ts`
`drawGraticuleTenant`). Currents are already cell-bound (parallax-free); the 2D
map has no parallax; smooth-globe mode is unchanged.

### Draped radius

For each graticule sample (a unit direction on the sphere), find the nearest
display cell and project the point at that cell's rendered radius, clamped to sea
level:

```
r = 1 + max(cell.height, seaLevel) * 0.05
```

The sea-level clamp means over ocean the grid sits at the constant sea-level
radius (identical to the S17 shipped behavior — the ocean-limb screenshot stays a
valid regression anchor), and over land it rides the terrain, meeting it exactly
at coastlines (same formula the mesh uses).

Applies **only in the raised path**. When `smooth` is true the graticule stays at
`r = 1` (S17 behavior) — there is no relief to follow.

### Nearest cell = greedy neighbor walk (no spatial index)

The display cell graph is a **symmetric, complete Voronoi adjacency**
(`utils/worldGen.ts:445-446` pushes both directions), so a hill-climb over
`dot(dir, cell.center)` cannot dead-end and terminates at the true nearest cell
(Voronoi ⇒ the nearest-center function is unimodal on the graph).

`nearestCellWalk(cells, dx, dy, dz, startId) → cellId`: from `startId`, repeatedly
hop to the neighbor with the highest `dot(dir, center)`; stop when no neighbor
improves. Cap iterations (e.g. `cells.length`) as a cycle guard (dot is strictly
increasing, so the cap is only defensive).

**Seeding:** the graticule samples ~5100 points (≈53 lines × 97 samples). Carry a
single `startId` across the ENTIRE draw: brute-force scan for the first sample
only (one O(n) per redraw), then every subsequent sample (including each line's
first) walks from the running `startId`. Consecutive samples on a line are ~3.75°
apart → 1–3 hops; line-to-line jumps cost a few more, bounded by graph diameter
(~√n). One O(n) per redraw + ~5100 short walks — cheap even at the 200k cap, and
the overlay redraw is gated (not per-frame). No index build, no per-world cache.

### What stays

- **Grid→smooth coupling stays** (advisor: it's the safety net for the queued
  ScreenOverlay tenants — roads/routes, contours, borders, labels — none of which
  drape; and the fallback if drape ever looks wrong). Drape makes the
  *raised+grid* path parallax-free for users who turn smooth off manually; it does
  not replace the coupling. Dropping the coupling is a one-line change Matt can
  make later with eyes on screen.
- `smoothGlobe` toggle, `displayRadius`, currents, all S17 behavior — untouched.

### Redraw gate

Draped radii depend on `world` (heights mutate in place on paint). The gate key is
camera+globe matrix + tenant ids + smooth flag — a paint stroke doesn't bump it,
BUT `ScreenOverlay`'s `[nCells, world]` effect already resets `lastKey.current=''`
on any `world` change (S16), so paint-driven height changes force a redraw.
**Verify** this covers it rather than adding a world term to the key.

## Testing

- **Unit** (`tests/nearestCell.test.ts`): on a small hand-built cell graph, the
  walk returns the true nearest center from several start points; monotonic dot;
  terminates.
- **In-browser** (the discriminating views, raised globe, grid on + smooth OFF):
  1. grid crossing a **coastline** at high zoom — meets the terrain at the
     shoreline, no step-off;
  2. grid over a **mountain range** at high zoom — the line ripples with relief,
     doesn't cut through peaks or float above them;
  3. ocean limb unchanged vs S17.
  Plus: 0 console errors, perf fine at 5k, gates green.

## Files

- `utils/nearestCell.ts` (new) — `nearestCellWalk`.
- `tests/nearestCell.test.ts` (new).
- `components/overlays/tenants.ts` — `drawGraticuleTenant` draped raised path.
