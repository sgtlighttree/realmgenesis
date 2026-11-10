# s08 (2026-11-10, local M1): equirect patch projection, continuous fidelity

Picks up s07's queue: item 6 (equirect patch projection, "counter distort") and 7b (continuous fidelity + tolerance sliders).

## Setup
- On battery, Blender open (TCHQ session `claude-playground` active), OBS running. Builds capped (`taskpolicy -b -j 4`). TCHQ session messaged: it tells me before a render/bake; I check Blender's CPU before each build. Disk 24 GB free at start.
- Registry had lost crates (`accesskit` etc.); on home Wi-Fi, so `cargo fetch --target aarch64-apple-darwin` (allowed by the SSID rule). `cargo metadata --offline` still fails on Linux-only crates; that's expected with a target-filtered fetch and doesn't affect builds.
- Matt (this session): on a Fable session the advisor is Fable too, so consult it less often. One consult taken at the start (below).

## Advisor verdict (commitment boundary: Placement field across ~8 files + solver change)
**Go on both, one trap to fix:**
1. `solver.rs` uses the one `guided` weight for two jobs: the uplift blend and the strength of the final clamp (`h += (clamped − h)·g`). A fractional freedom fed through it would enforce the tolerance only partly. → Blend with freedom × pin; clamp at pin strength wherever freedom > 0; keep the `g > 0` guard so Exact stays bit-identical.
2. Tolerance needs its own per-node field beside `guided`; fill uncovered nodes with 0.10 (the field is interpolated).
3. Gate: `world_hash` of an Exact and a Guided import before/after. Exact identical; Guided at freedom 1 / tolerance 0.10 identical or explained.
4. Equirect: `pixel`/`dir_of_pixel` exact inverses; wrap Δlon; clamp latitude in `dir_of_pixel` (overlay and handles call it past the poles).
5. Don't default the new field to azimuthal at pass-through sites (`import.rs` rescale, `svg.rs`, `coasts.rs`, `app.rs`): it would silently reset the projection.
6. New JSON keys omitted at defaults; old `"fidelity"` strings still read; no schema bump (the `peak_boost` precedent).

## Decisions
- **Equirect patch = normal-aspect plate carrée about the patch centre**, equal degrees per pixel on both axes, `km_per_px` measured along a meridian (= at the equator for longitude). Rotation turns the image in its own plane before the mapping. **Rejected: oblique equirect** (the patch centre as its own equator): it moves rigidly, but it isn't what an Azgaar map is; a map placed at 50° N must stretch as plate carrée does there. Assumption (not looked up): Azgaar's globe treats the map as plate carrée with equal degrees per pixel.
- **Fidelity → `freedom` 0..1 + `tolerance`** (fraction of height, 30 m floor kept). `freedom` 0 = Exact, 1 = today's Guided. The `Fidelity` enum goes.

## Log
- **Baseline hashes (HEAD `e333cf4`, Klerga V1 patch 20,10 at 0.5 km/px, sea grey 0.196, smooth 3, blank ocean, seed 1, Normal):** Exact `height_m 2b41e227df32d751`, drainage `7179450036100856`, lakes `6c44c4418540a1f8`, 122 rivers. Guided `2201ac4cd8bda070`, `5c7938b06081be18`, `ee2aed5686daf040`, 118 rivers. Command: `RG_MASK=<png> RG_MASK_PATCH=20,10,0.5 RG_MASK_THRESHOLD=0.196 RG_MASK_OUTSIDE=0 RG_NO_CACHE=1 RG_MASK_HEIGHT=4000,3,<0|1> world_hash 1`.
- **Built (untested until the build below reports):**
  - `mask::Projection { Azimuthal, Equirect }` as a field of `Placement::Patch`. `pixel`/`dir_of_pixel` branch on it; Δlon wraps about the centre; `dir_of_pixel` stops at the poles. The overlay, outline and handles go through `dir_of_pixel`, so they follow. The scale handle measures the real corner distance for equirect.
  - Pass-through sites carry the projection: `import::user_placement`/`raster_placement` (SVG and step rescale), `svg::coverage_from_svg`, the panel (`CoastUi.projection`), `RG_MASK_PROJ=equirect` (app and probes).
  - `Elevation.fidelity` → `freedom` (0..1) + `tolerance`. `Imported.tolerance` → `Snapshot.guided_tolerance` → solver. The clamp's strength is now the pin weight; it runs only where freedom > 0.
  - World file: `placement.projection: "equirect"` only when set; `elevation.fidelity` still "exact"/"guided" at the ends, `elevation.freedom` only in between, `elevation.tolerance` only when not 10%. No schema bump.
  - UI: projection radio ("regional" / "flat lat/lon"), span shown in degrees plus a **"fit 360°"** button for equirect (my addition: it's the "scaled to world size" step); "erosion freedom" and "tolerance %" sliders replace the exact/guided radio. New heightmap layers still start at freedom 0.
  - Test edits: existing tests got the new constructor field (`projection: Azimuthal`, `freedom: 0.0`/`1.0`) only; no assertion changed. New tests: `equirect_patch_round_trips_and_is_lat_lon_linear`, `freedom_and_tolerance_are_honoured`, `projection_freedom_tolerance_round_trip`.
- **Known limits (equirect):** bilinear sampling clamps at the image edge, so a map fitted to exactly 360° has a one-pixel seam opposite its centre (whole-world placement wraps; a patch doesn't). The mip footprint uses the meridian pixel size, as whole-world does, so high latitudes are sampled a little sharp. With a rotation, an image wider than 360° maps each point once (the wrapped one).
- **TCHQ coordination:** its session reported swap at 6.7 of 7.2 GB (RAM is the tight resource), so builds run at `-j 2`. It asked for a ~10 min hold on new builds during a 306 MB OBJ import; held.
- **Finding (test, n=1 map: the dome at res 64):** with the straight uplift blend `u = u_inv + (smooth − u_inv)·g`, freedom 0.3 moved pinned land a mean 183 m against 180 m at freedom 1: no gradation, both sit on the 10% clamp. Read: the inverted uplift varies by orders of magnitude between neighbouring nodes, so any linear share of the smoothed rate swamps it. **Fix: blend in log space** for 0 < g < 1 (`u = exact·(smooth/exact)^g`, `exact` floored at 0.1% of the base rate); g ≥ 1 keeps the old formula so a full Guided layer's interior is unchanged. **Side effect:** a Guided layer's feathered edge (g = pin < 1) now uses the log blend too, so existing Guided worlds change in that band. Mask tests: 10/10 pass.
- **Refuted: the log-space uplift blend.** Same test: freedom 0.3 moved 185 m against 180 m. Any share of the smoothed uplift takes the land past the clamp, so the uplift is not a usable continuous knob (n=1 map, two blend formulas). Removed.
- **Decision: freedom = the share of the guided change that is kept.** Any freedom > 0 runs the full guided solve (uplift blend weight = pin, as before s08); then `h = h0 + freedom · clamp(h_guided − h0, ±tolerance)`. Predictable (same valleys, shallower), monotone, never past the tolerance, and freedom 1 is today's Guided exactly. Rejected: scaling the tolerance by the freedom (two sliders doing one job). Cost: a blend of two drained surfaces can leave small pits; the final routing fills them as it does for any surface. This departs from the advisor's item 1 (clamp at pin strength): that would make every freedom > 0 identical.
- **Disk: 18 GB free** mid-session (24 at start). `target/` is still 1.6 GB, so it isn't this project; swap went 7.2 → 9.2 GB while Blender worked. Flagged to Matt.
- **Tests:** `rg-pipeline import::` 10/10 with the new freedom meaning (`freedom_and_tolerance_are_honoured` passes: freedom 0.3 moves less than freedom 1, all three cases inside the tolerance).
- **Hash gate (verified, same command as the baseline):** Exact (freedom 0) is **bit-identical** on all four lines. Guided (freedom 1, tolerance 10%) first came out with the same drainage and river count but different height bits; cause confirmed by the fix: the tolerance arrived as f32 0.1. The solver now uses the exact constant when the tolerance is the default, and Guided is **bit-identical** too (`2201ac4cd8bda070` / `5c7938b06081be18` / `ee2aed5686daf040` / 118). Freedom 0.5 gives a third world (`74309c16d0f5c3e2`, 121 rivers).
- TCHQ asked for two more ~10 min holds (EEVEE renders); held each time.
