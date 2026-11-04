# R3: coastline conditioning and the plausibility pass (design, s04)

Status: **draft, awaiting fresh review and Matt's approval.** No code yet.
Log: `docs/session-logs/2026-11-04-s04-r3-design.md`. Roadmap row: `ROADMAP.md` → R3.

## Goal

Coastlines become an **input**: procedural (as now) or imported (raster or vector), placed on the sphere, and edited with a brush. The tectonics are then **fitted to the coasts**, not the other way round:
- Plate boundaries lie where the coasts say they should:
  - active coasts get a trench offshore and a cordillera or arc;
  - passive coasts sit inside a plate, with a wide shelf and a sediment wedge;
  - mid-ocean ridges run between rifted coasts, mirroring them.
- Coastal bathymetry follows from the same classification.
- Erosion and climate run as today. Editing means changing the mask, then regenerating with the same seed.

Decided earlier (s03, don't re-litigate): continents first, plates fitted by a plausibility pass. Active coasts are automatic by default and can be repainted ("both"). No history simulation.

## What changes in the pipeline

Today (`rg-tectonics`): `crust` (noise continents) → `plates` (random seeds, random Euler poles, kind read off the motion) → `relief` (profiles by distance to boundaries) → sea-level bisection (the coast emerges from the heights) → `finish` (shelves from the emergent coast).

R3:

```
0. mask      land mask: procedural, imported or edited → per-node land coverage, coast SDF
1. coast     landmasses → coast segments → seaways (medial axis) → class per seaway/segment
2. crust     continental crust = land + a margin whose width depends on the coast class
3. plates    partition fitted to the classes; big landmasses split along sutures
4. motion    Euler vectors solved by least squares to the intended boundary kinds
5. relief    today's profiles on the fitted boundaries (+ passive-margin escarpments)
6. conform   heights forced to the mask's sign (strict) or allowed to add islands offshore
7. finish    shelves/slopes from the coast class (replaces the `activity` noise)
```

New modules: `rg-tectonics/src/mask.rs`, `coast.rs`, `fit.rs` (partition + motion solve). `crust.rs` keeps the procedural continent noise. `plates.rs` keeps boundary measurement and `BoundaryNode`, but no longer rolls the poles or decides the kinds. File decoding (PNG, GeoJSON) lives in `rg-pipeline/src/import.rs`; the core crate only sees a decoded `LandMask`.

## 0. The land mask

```rust
pub enum LandMask {
    /// Today's continent noise, thresholded to `land_fraction`.
    Procedural,
    /// A coverage raster (0..1, 0.5 = coast) with a placement on the sphere.
    Raster { data: Vec<f32>, w: u32, h: u32, place: Placement },
    /// Polygons in lon/lat (degrees), holes allowed, with an optional rigid rotation.
    Polygons { rings: Vec<Ring>, rotate: Rotation },
}
pub enum Placement {
    /// 2:1 image covering the whole sphere, plate carrée.
    WorldEquirect,
    /// A regional map: centre lat/lon, km per pixel, rotation (degrees clockwise from north).
    Patch { center_deg: [f64; 2], km_per_px: f64, rotation_deg: f64 },
}
```

- **Per-node coverage** `c(x) ∈ [0, 1]` (bilinear in the raster; 0/1 inside or outside a polygon, with a sub-cell distance blend at edges). Land = `c > 0.5`. Sampling a continuous coverage at node positions, rather than rasterising at one fixed resolution, keeps coasts **resolution-coherent** (the s03 lesson).
- **Patch projection: azimuthal equidistant** about the centre. Distances from the centre are exact. The transverse scale error is θ/sin θ: +1.7% at 2,000 km from the centre, +4% at 3,000 km, +11% at 5,000 km. The UI should warn above ~3,000 km radius; at that size a world-equirect image fits better.
- **Outside a patch:** either procedural continents with a guard band (no procedural land within ~800 km of the patch's land), or ocean. Recommended default: procedural, so the planet isn't one continent in an empty ocean.
- **Polygons:** GeoJSON `Polygon`/`MultiPolygon`. Edges are densified to ≤ 50 km great-circle steps. Inside/outside uses a winding test in a local azimuthal projection per polygon; polygons bigger than a hemisphere are split. `Rotation` moves a chosen point to a chosen point and spins about it, so a real-world outline can be relocated.
- **SVG (later):** the same `Patch` placement, with viewBox units scaled to km.
- **Brush strokes**, applied in order on top of the base mask: `Stroke { path: Vec<DVec3>, radius_km, op: Land | Sea | Active | Passive | Unpaint }`. They are resolution-independent (a distance-to-polyline test) and small enough to store as JSON. The brush UI is R3.3.
- Islands below a size floor (default 200,000 km²; Britain 209,000, Iceland 103,000) are **islands**, not landmasses. Ones on a landmass's shelf join its crust. Others get oceanic crust (volcanic: an apron, no shelf).

## 1. Coast classification: seaways first, coasts inherit

Classifying each coast segment independently with noise gives short active/passive/active runs, which Earth doesn't show. Instead:

1. **Landmasses** = connected land components above the size floor.
2. **Seaways**: run a multi-source distance over the ocean from every coast node, labelled by landmass. Where two labels' fronts meet is the **medial axis** of the ocean between them. The pair of coasts that face each other across it is a **seaway** `(A, B)` with half-width `D` (coast to axis). Coast that faces no other landmass within ~6,000 km faces "open ocean", a seaway with itself.
3. **Class per seaway**, rolled per seed:
   - **Opening** (Atlantic): both coasts passive; a ridge on the medial axis.
   - **Closing** (Pacific): both rims active, with oceanic plates in between.
   - **One-sided** (Mediterranean, Caribbean; the default when `D` is under ~300 km): one rim active and the other passive, so the passive side's plate dives under the active side's.

   The odds of "closing" rise with `D` (wide oceans are Pacific-like). One per-world knob, **active share** (default 0.3 of continental coast length; see the Earth numbers below), steers the dice. The rolls are made per seaway, so they stay put across resolutions.
4. **Segments:** a seaway's coast is cut into segments of at least ~1,000 km. Within a closing seaway, a segment may roll passive (default 15%; e.g. a rifted microcontinent), so long rims aren't uniform.
5. **Painted overrides** (`Active`/`Passive` strokes) win over the rolls. They re-split segments at the stroke ends; the minimum length doesn't apply to painted runs.

Earth guide (for the default share; **verify before quoting**):
- modern passive margins total ~105,000 km (Bradley 2008, *Earth-Science Reviews*);
- subduction zones total ~55,000 km (Lallemand et al. 2005), a large part of it intra-oceanic.

Active continental coast is therefore roughly 25–35% of continental-margin length.

Output: a per-node `margin` field (0 none, 1 passive, 2 active). It is stored (schema 5) and drawn as a new **Margins** layer with a legend.

## 2. Crust from the mask

- Continental crust = land, plus a margin offshore:
  - **passive:** shelf + slope, ~150–250 km (varied by the existing `margin_w` noise);
  - **active:** the forearc, up to the trench offset minus ~20 km.
- `edge_m` (the signed distance to the crust edge) keeps its meaning for every consumer (sediment wedge, plate costs, relief weights).
- A new `coast_m` is the signed distance to the **mask coast**. The continental base height ramps up from the coast rather than the crust edge, so the coast sits at 0 by construction. Today's `-130 + interior·ease(edge_m)` would put the coast near +600 m once the crust edge is 200 km offshore.

## 3. Plate partition fitted to the classes

Two phases, both multi-source growth on the existing `Grid::grow_regions` (noisy costs, so borders aren't geometric):

1. **Continental plates.** Each landmass is a plate. Large ones are split; see the sutures below. A plate claims:
   - its land;
   - its ocean up to the medial axis of every **opening** seaway (the ridge becomes the boundary between the two continental plates, the South America/Africa pattern);
   - its ocean up to the **trench offset** (default 120 km) off every **active** segment. Fronts cannot pass this line.
2. **Oceanic plates.** The ocean left unclaimed (the interiors of closing seaways and open ocean) is split among oceanic seeds:
   - count = `plate_count` minus the continental plates, at least one per closing seaway;
   - placed at random places (`nearest_node(unit_vec)`, the s03 rule).

   Boundaries between oceanic plates take their kind from the solved motion (ridge, transform or intra-oceanic subduction with an arc), as today.

**Sutures inside large landmasses** (the advisor's point: imported continents would otherwise have no collision belts, and land relief is what Matt is happy with):
- Landmasses over a size threshold (default 8 M km²; Australia 7.7, Eurasia/Africa/the Americas above) are split into 2–3 plates, grown from seeds inside them.
- Each internal boundary rolls a kind:
  - **collision** (the default, 60%): Himalaya/Alps;
  - **rift** (25%): East Africa;
  - **dormant** (15%): an old suture; it becomes an old belt, with no plate boundary motion target.
- `plate_count` is a target, not exact. The UI shows the achieved count.

Microcontinents and islands join the plate around them.

## 4. Motion: solved, then used for rates only

Unknowns: one rotation vector `w_k = omega_k · pole_k` per plate (3 numbers each). At a boundary sample `p` between plates `a` and `b`, with unit normal `n` from `a` to `b`:

```
convergence rate  c(p) = ((w_a − w_b) × p) · n      // linear in the w's
```

- **Targets**, weighted by boundary length:
  - active/subduction: `c = +55 mm/yr`;
  - collision: `+30`;
  - ridge: `−40` (full rate);
  - rift: `−10`;
  - free boundaries between oceanic plates: no target.
- **Regularisation:** `λ Σ |w_k − w_k⁰|²`, with `w_k⁰` rolled as today. This gives each plate its own direction and obliquity, instead of the minimum-norm answer.
- **Frame:** no-net-rotation (area-weighted Σ w = 0). Only hotspot track directions depend on it.
- Size: 3K unknowns (K ≤ 60), so the normal equations are at most 180 × 180. Plain Gaussian elimination, deterministic, no new dependency.

**The kind comes from the classification, not from the sign of the solved `c`** (the advisor's point 1). Otherwise we're back to "active if the dice agree", just through the solver's residuals. The solve supplies:
- **rate**, feeding `strength`, `speed_mm_yr` and uplift. Clamped to a minimum (subduction ≥ 15 mm/yr, ridge ≥ 10) where the residual has the wrong sign or is too small;
- **obliquity**. On ridges, high obliquity splits the boundary into transform pieces (fed to R3b). Active margins stay subduction when oblique (strain partitioning, as at Sumatra).

Wrong-sign samples are counted and printed by the probe; target < 10% of prescribed boundary length.

`Plate::omega > 0 && < 0.05` (a test) still holds: 55 mm/yr is 0.0086 rad/Myr at Earth radius.

## 5. Relief on the fitted boundaries

`relief.rs` mostly runs unchanged: the profiles key on distance to boundary nodes of each `Kind`. Changes:
- **Subduction geometry follows the coast.** With the boundary at the trench offset, the trench lies ~120–150 km offshore, and the cordillera (170 km from the boundary) lands ~50 km inland. Andes check: crest ~150–250 km from the trench axis.
- **Passive-margin escarpments** (new, on a share of passive segments, default 30%: southern Africa, the Western Ghats, eastern Australia, southeast Brazil): a seaward-facing escarpment 50–200 km inland, 600–1,500 m of plateau behind it, and a low coastal plain. These are high-elevation passive margins; most passive margins are low.
- **Seafloor age** keeps the GDH1 law and the Parsons area–age calibration. Opening seaways put the ridge on the medial axis, so the oldest floor lies at the passive coasts (Atlantic pattern). The sediment wedge stays tied to `edge_m`.
- **Interior old belts** (noise) stay; dormant sutures add to them.

## 6. Conforming to the mask

- **Strict** (the default for imports): no sea-level bisection.
  - Land nodes: `h = max(h, ramp(coast_m))`, a small positive ramp (~2 m at the coast, rising ~1 m/km, out to ~20 km).
  - Sea nodes: heights above −200 m are softly compressed into [−200, −1] m, so an arc or hotspot crest becomes a shallow bank, not a flat sheet at −1 m.
  - Land fraction = the mask's; the slider is disabled and shows the value.
- **Islands allowed** (the procedural default, an option for imports): tectonic features (arcs, hotspots) may rise above sea level, but only beyond a guard distance from the mask coast (default 100 km). Drawn coasts never move.
- **Procedural mode:** the mask is the continent noise thresholded straight to `land_fraction` (not 1.4×), plus one coast-detail fBm octave, because today the bays and islands come from height noise crossing zero. **Every procedural world changes again**, coastlines included.

## 7. Coastal bathymetry

`finish` swaps the `activity` noise for the `margin` class:
- **passive:** shelf 80–200 km (mean ~80 km on Earth; `margin_w` varies it), slope, then the rise with the existing sediment wedge;
- **active:** shelf 10–30 km, steep slope into the trench;
- **volcanic islands:** apron only.

Mid-ocean ridges come from step 3. Open-ocean fabric (transforms, fracture zones, abyssal hills) is R3b.

## Gates (invariants to test, not assumptions)

Unverified until a test says so (log: "n=0"):
1. **Coast survival (strict):** sign(final h) = mask at the committed resolution after the full pipeline (erosion, lakes, hypsometric calibration). Target ≥ 99.9% of nodes, and every miss within one cell of the coast. Risks:
   - resampling to the erosion grid flips coastal nodes;
   - hypsometric calibration pushes low land ≤ 0;
   - `ocean = h ≤ 0` (`solver.rs:593`) turns masked land below sea level into ocean.

   Decision needed: **ocean comes from the mask** (carried in the snapshot), not from the sign of h. Until then, land is clamped ≥ 1 m, so there's no Dead Sea or Caspian (a documented limitation).
2. **Plausibility:**
   - ≥ 90% of active coast length has a trench within 300 km;
   - no trench within 400 km of a passive coast;
   - ridge nodes of opening seaways have `|d_A − d_B| / (d_A + d_B) < 0.1`;
   - seafloor age rises toward passive coasts;
   - wrong-sign boundary length < 10%.
3. **No regressions:**
   - `probe_bathymetry`: mean depth within 3,688 ± 60 m;
   - `probe_realism`: hypsometry within today's bands;
   - `probe_res_coherence`: height r ≥ 0.9 at 320 vs 640, same plate ≥ 95%.
4. **Cost:** tectonics time at Normal and Ultra, measured before and after. Expect +20–40% (3–5 extra distance passes; the solve itself is trivial).

## Tests that will change (need Matt's OK before editing)

- `rg-tectonics`:
  - `land_fraction_hits_target`: procedural still targets `land_fraction`, but islands add a little. Tolerance may need widening.
  - `heights_are_plausible`: probably passes; re-check.
- Goldens: `rg-pipeline`/`rg-world` determinism goldens for regenerated worlds will change (every world changes).
- New tests: the gates above, `mask` sampling (equirect, patch, polygon, resolution coherence of the coast), the seaway classification (deterministic per seed, minimum segment length), and the motion solve (prescribed signs on a synthetic two-plate case).

## World file: schema 5

- Header `conditioning`:
  - `mode` (procedural | imported);
  - the `placement` or `rotation`;
  - `strict`, `islands_guard_km`, `active_share`, `trench_offset_km`, `escarpment_share`, `suture_area_km2`, `landmass_min_km2`;
  - `strokes` (JSON).
- Chunks:
  - `conditioning.source` (the original PNG or GeoJSON bytes, so the world regenerates anywhere);
  - `tectonics.margin` (u8 per sample, committed).
- Schema-4 files load as `mode: procedural`, with `margin` absent. The Margins layer is then hidden.

## Build phases (each gated)

- **R3.1, the plausibility pass on procedural continents.** Steps 0–7 without import or UI, plus an "active share" slider and the Margins layer. Gate: the plausibility and no-regression gates, plus Matt's look in the app. This is where the science gets judged.
- **R3.2, import + placement + strict conformance + schema 5.**
  - PNG first (world equirect and patch), then GeoJSON.
  - The "ocean from the mask" change and the coast-survival gate.
  - UI: source picker, placement fields (centre, km/px, rotation), the strict / islands option.
- **R3.3, editing.** The brush on the globe (needs ray picking against the terrain, an M1c item), land/sea and active/passive strokes, regenerate.
- Then **R3b**, open-ocean fabric, on R3's ridges and obliquity.

## Rejected alternatives

- **Mask imposed after tectonics** (override the crust, adjust uplift and sea level): plates ignore the coast, so trenches and ranges land in the wrong places. That's the problem R3 exists to fix.
- **Per-segment noise classification:** short active/passive alternations; nothing ties conjugate coasts together.
- **Random poles, rerolled until the classes match:** no guarantee, and slow. The least-squares solve is direct.
- **Kind read off the solved velocity:** the classification would no longer be in charge (see §4).
- **Full plate-history simulation** (Cortial et al. 2019 style): rejected in s03; a long-term stretch goal.
- **Keeping today's emergent coastlines in procedural mode** (run today's tectonics once to get a coast, then fit): it doubles tectonics time (Ultra ≈ 19 s). The thresholded noise plus a coast-detail octave looks similar.

## Open questions for Matt (recommended answer first)

1. First import format: **PNG mask (world equirect + placed patch)**, then GeoJSON; SVG later?
2. Outside an imported patch: **procedural continents with a guard band**, or empty ocean?
3. Imports: **islands allowed beyond 100 km of drawn coasts**, or strictly none?
4. Procedural coastlines change again (and every world with them). OK?
5. Land below sea level in imports (Dead Sea, Caspian): **not in R3.2** (land ≥ 1 m), later with "ocean from the mask"?
6. OK to update the tests and goldens listed above once R3.1's gates pass?
7. Default active share **0.3**, and closing seaways get oceanic plates between active rims. Does the Pacific/Atlantic split read right for fantasy worlds, or should "all passive" (a young, rifting world) be a one-click preset?
