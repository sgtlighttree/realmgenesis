# RealmGenesis v2: roadmap

Current state and evidence: `HANDOFF.md`. Detail and history: `docs/session-logs/`. Terrain design: `docs/design/terrain-model.md` and `docs/design/drained-detail-network.md`.

## Where terrain stands (end of s01)

The terrain is physically grounded and checked against Earth by the same code that scores ours (`examples/earth_reference.rs`):
- **Drainage:** closed depressions 0.00% on the final surface, by construction.
- **Slopes:** mountains 10–18° mean (Appalachian-class), lowlands 2–5°.
- **Hypsometry (emergent, not forced):** > 2 km 9.1%, > 4 km 2.0% (Earth 6–9% / ~2%).
- **Oceans:** mean depth 3626–3709 m across seeds (Earth 3688).
- **Climate:** land runoff 333 mm/yr (Earth ~300).
- **Lakes and erosion:** climate-driven erosion; lakes and closed basins from a water balance; exact mass balance.

## Next: terrain (choose one per session)

1. **Lake rims + peaks** *(recommended next: visible artifacts, moderate effort)*
   - Closed-lake rims leak: valleys of a neighbouring basin cut below the lake level across the macro divide, so water stands as a wall on part of the shore. Fix: inside closed catchments, only wedges draining to that lake count, with a continuous rim at catchment boundaries (beware: min/max compositions can create pits).
   - Peaks top out near 5 km (Earth 8.8). Real summits stand 2–3 km above smoothed plateau levels through deep valley and glacial incision.
2. **Sediment & deltas** *(big scientific step, multi-session)*: transport-limited deposition (Yuan et al. 2019, implicit ξ–q). Alluvial plains, fans and deltas fed by the erosion flux, filled basins, and margin sediment tied to river mouths (today GlobSed's wedge is uniform).
3. **Crustal isostasy** *(fixes "heights fitted" at its root, multi-session)*: crustal thickness from shortening with a 65–70 km cap, Airy isostasy, then flexure (Te 20–40 km). Rock uplift = thickening rate (replaces the exp(−d/L) belts and the 7 km elevation cap). Hypsometry stays a check only.
4. **Endorheic share:** we get 2.3–4.5% of land vs Earth's ~20%. 90–95% of ours is arid (the mechanism is right). Needs a fair Earth measurement at macro resolution (HydroSHEDS-conditioned DEM; the 0.25° test is dominated by resampling artifacts on both sides), a better arid PET than Thornthwaite, and tectonically maintained basins.
5. **Slope–area θ:** ours is −0.04–0.27 vs Earth's 0.19–0.45. Channel heads are relief-capped.
6. **Glacial regime:** ELA "buzzsaw", U-shaped valleys, cirques, fjords. **Seafloor fabric:** ridge–transform segmentation, axial morphology vs spreading rate, ridge-parallel abyssal hills.
7. **Validation debt:** the bathymetry age remap and dynamic-topography std were tuned on the bands they're checked against; Horton/Tokunaga ratios and per-basin HI on the final surface; shelf width and continental-slope gradient.

## Next: the "alive world" layer

From `docs/research/2026-11-02-alive-world-features.md` (licences verified there). In dependency order:
1. Population / habitability.
2. Names (our own phoneme grammars; don't copy FMG's name lists).
3. Burgs, capitals, ports.
4. Cultures (cost-weighted flood fill).
5. States / borders.

Then labels, markers (waterfalls from knickpoints, volcanoes on arcs), roads and trade, provinces, religions, and later military, diplomacy and heraldry.

Licence rules: FMG (MIT), mewo2 (MIT) and mapgen4 (Apache-2.0) may be ported with attribution. watabou TownGenerator and planet_heightmap_generation are GPL-3.0: ideas only. Armoria has no LICENSE file, so hold off. WappenWiki charges are CC BY-NC: blocked.

## Engineering

- **Tile cost:** network tiles cost 85–240 ms on the CPU (old noise detail: 9–33 ms). Native copes on worker threads; the web build streams slowly. Plan: the GPU port (`rg-gpu`), since the per-patch channel buffers port directly; per-level window reuse across tiles in the meantime.
- **Generation time:** ~6 s at Normal (two climate passes + lakes per step).
- **Cleanup:** retire the old hills/gully detail once the app no longer needs the fallback (the MPL port stays in its own file until removed deliberately). Store macro receivers instead of snapping `flow_dir`.
- **App:** M1c items (close-up albedo, LOD normal seams, picking), tier B climate, vector layer, sculpting.
