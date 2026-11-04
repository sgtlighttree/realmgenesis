# R3: coastline conditioning and the plausibility pass (design, s04)

Status: **built in s05** (R3.1a, R3.1b, R3.2 PNG + SVG, R3.3 MVP; deviations and gate results in `docs/session-logs/2026-11-05-s05-r3-build.md`). Approved by Matt (s04), design v2 (fresh review: *go with changes*, all adopted, marked "review").

Matt's answers to the open questions (s04, 17:00):
1. **PNG first, SVG next, GeoJSON later.** This changes the R3.2 order.
2. Outside an imported patch: **procedural continents and other landmass types, kept a distance away** (default).
3. Islands beyond 100 km of drawn coasts: **yes** ("so it doesn't feel too lonely").
4. Procedural worlds changing: OK.
5. All land stays **≥ 1 m above sea level**.
6. Updating the affected tests and goldens: OK ("project is in flux").
7. Active share 0.25, plus an **"all passive" preset**: yes.
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
- **Brush strokes** (R3.3), applied in order on top of the base mask: `Stroke { path: Vec<DVec3>, radius_km, op: Land | Sea | Active | Passive | Unpaint }`. They are resolution-independent (a distance-to-polyline test) and small enough to store as JSON. The brush UI is R3.3.
- Islands below a size floor (default 200,000 km²; Britain 209,000, Iceland 103,000) are **islands**, not landmasses. Ones on a landmass's shelf join its crust. Others get oceanic crust (volcanic: an apron, no shelf).

## 1. Coast classification: seaways first, coasts inherit

Classifying each coast segment independently with noise gives short active/passive/active runs, which Earth doesn't show. Instead:

1. **Landmasses** = connected land components above the size floor.
   - Components and areas come from the coverage at a **fixed canonical resolution** (res 256, coverage-weighted area), not the generation grid. Thin isthmuses then split or merge the same way at 320 and 640 (review).
   - Landmasses are numbered by sorted centroid position, not scan order.
2. **Seaways**:
   - `nearest_source` over the ocean from every coast node; each `src` maps back to its landmass.
   - Where two labels meet is the **medial axis**.
   - A seaway is one **connected component of the medial axis** (flood fill over the axis nodes), not a landmass pair. The Africa/South America pair meets twice, across the Atlantic and the long way round across the Pacific, and those must roll separately (review).
   - The seaway's key (for its dice) is its quantised axis centroid, a position, never an index.
   - Its half-width `D` is the median coast-to-axis distance.
   - Coast facing no other landmass within ~6,000 km faces "open ocean", treated as a closing seaway.
3. **Class per seaway**, rolled per seed:
   - **Opening** (Atlantic): both coasts passive; a ridge on the medial axis.
   - **Closing** (Pacific): both rims active, with oceanic plates in between, **and a ridge on the medial axis**. A basin consumed from both sides needs a spreading centre, and seafloor age needs a ridge to count from (review).
   - **One-sided** (Mediterranean, Caribbean; the default when `D` is under ~300 km): one rim active and the other passive, so the passive side's plate dives under the active side's.

   The odds of "closing" rise with `D` (wide oceans are Pacific-like). One per-world knob, **active share** (default 0.25 of continental coast length; see the Earth numbers below), steers the dice. The rolls are made per seaway, so they stay put across resolutions.
4. **Segments:** a seaway's coast is cut where a position-keyed 1D noise along the coast peaks, with a minimum of ~1,000 km. The cuts then fall in the same places at every resolution; cutting "every 1,000 km from a start node" would make the start node arbitrary (review). Within a closing seaway, a segment may roll passive (default 15%; e.g. a rifted microcontinent), so long rims aren't uniform.
5. **Painted overrides** (`Active`/`Passive` strokes) win over the rolls. They re-split segments at the stroke ends; the minimum length doesn't apply to painted runs.

Earth guide (s04 fact-check):
- **Verified** (Bradley 2008, *Earth-Science Reviews* 91, p. 1–26): "The present-day passive margins have an aggregate length of 105,000 km, even longer than the spreading ridges (65,000 km) or the convergent plate boundaries (53,000 km)."
- **Unverified (search snippet only):** intra-oceanic systems are "some 40% of the subduction margins" (Leat & Larter 2003).
- Together: continental active margins come to ~32,000 km, which is **~0.23** of continental-margin length. Default active share: **0.25**.

Output: a per-node `margin` field (0 none, 1 passive, 2 active). It is stored (schema 5) and drawn as a new **Margins** layer with a legend.

## 2. Crust from the mask

- Continental crust = land, plus a margin offshore:
  - **passive:** shelf + slope, ~150–250 km (varied by the existing `margin_w` noise), **capped at 0.4·D**. Otherwise the two margins of a narrow seaway overlap and its ridge lands in continental crust, which `plates.rs` would call a `Rift` (review);
  - **active:** the forearc, up to the trench offset minus ~20 km.
- All widths are continuous distances (`nearest_source` arc lengths), never cell counts, so they mean the same at every resolution.
- `edge_m` (the signed distance to the crust edge) keeps its meaning for every consumer (sediment wedge, plate costs, relief weights).
- A new `coast_m` is the signed distance to the **mask coast**. The continental base height ramps up from the coast rather than the crust edge, so the coast sits at 0 by construction. Today's `-130 + interior·ease(edge_m)` would put the coast near +600 m once the crust edge is 200 km offshore.

## 3. Plate partition fitted to the classes

`Grid::grow_regions` can't do this: it has one global cost array, with no pre-labelled or blocked nodes, and a high-cost "wall" only delays fronts (review, `grid.rs:281`). This is **new code**, not reuse.

1. **Continental plates** (no growth needed):
   - land nodes are pre-labelled with their landmass, or its sub-plate (see sutures);
   - ocean ownership comes from the same labelled `nearest_source` as the seaways, with distances warped by low-amplitude noise so ridges aren't geometric;
   - a plate keeps its ocean up to the medial axis of every **opening** seaway (the ridge is the boundary between the two continental plates: the South America/Africa pattern), and up to the **trench offset** (default 120 km) off every **active** segment.
2. **Oceanic plates:** the ocean left unclaimed (inside closing seaways, open ocean) is split by a new `grow_labelled(owner, blocked, seeds, cost)`, a plain multi-source flood restricted to unclaimed nodes.
   - At least one seed goes in **every connected unclaimed component**, including small pockets behind barrier lines.
   - Every closing seaway gets at least two, so its ridge is a plate boundary.
   - Seeds sit at random places (`nearest_node(unit_vec)`, the s03 rule). No "extra seeds per region" (that trick breaks connectivity here).
   - Connectivity is asserted (plus `every_plate_is_connected`).
   - Boundaries between oceanic plates take their kind from the solved motion, except a closing seaway's medial axis, which is prescribed as a ridge.

**Sutures inside large landmasses** (the advisor's point: imported continents would otherwise have no collision belts, and land relief is what Matt is happy with):
- Landmasses over a size threshold (default 8 M km²; Australia 7.7, Eurasia/Africa/the Americas above) are split into 2–3 plates, grown from seeds inside them.
- Each internal boundary rolls a kind, keyed by the landmass centroid:
  - **collision** (the default, ~70%): Himalaya/Alps;
  - **rift** (~30%): East Africa.
- An "old suture" isn't a plate boundary. It merges into one plate and paints an old belt instead (review: a dormant boundary costs a plate id and a special case).
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
- **The fit is soft.** A rigid rotation makes `c` vary sinusoidally along a boundary, so a subduction run spanning ≳60° of arc can't hold +55 mm/yr from end to end: its ends go oblique or the wrong sign (review). Earth's rates vary along strike too, so that's fine as long as the classification stays in charge (below).
- **`omega ≥ ε`:** no-net-rotation can leave a plate with `w ≈ 0`, which has no defined pole. Clamp `omega` to a small minimum and pick a deterministic pole (the rolled `w_k⁰` direction). This keeps `plates_are_valid_rotations` true.

**The kind comes from the classification, not from the sign of the solved `c`** (the advisor's point 1). Otherwise we're back to "active if the dice agree", just through the solver's residuals. The solve supplies:
- **rate**, feeding `strength`, `speed_mm_yr` and uplift. Clamped to a minimum (subduction ≥ 15 mm/yr, ridge ≥ 10) where the residual has the wrong sign or is too small.
  - **Uplift must read the same clamped rate.** Today `uplift.rs:~120` recomputes convergence from the plates' poles and omega itself, so a wrong-sign sample would get a mountain belt in the relief but no erosion uplift (review, confirmed in the code).
  - Fix: the snapshot gains a per-sample `convergence_mm_yr` field (clamped, committed, schema 5), and `uplift.rs` uses it when present. Old files keep the pole path.
- **polarity**: which plate dives (`sub_plate`) and whether the overrider is continental (`over_cont`) come from the classification. The active coast's plate overrides. Today they come from a crust-majority vote in a 3-hop ball (~100 km at res 320), and with a forearc narrower than a cell that vote is a coin toss that differs between resolutions (review).
- **obliquity**. On ridges, high obliquity splits the boundary into transform pieces (fed to R3b). Active margins stay subduction when oblique (strain partitioning, as at Sumatra).

Gate on the **length-weighted residual** (RMS of `c − target` over prescribed boundaries). Report the wrong-sign share as a diagnostic, not a pass/fail. Wide closing seaways will have some.

`Plate::omega > 0 && < 0.05` (a test) still holds: 55 mm/yr is 0.0086 rad/Myr at Earth radius.

## 5. Relief on the fitted boundaries

`relief.rs` mostly runs unchanged: the profiles key on distance to boundary nodes of each `Kind`. Changes:
- **Subduction geometry follows the coast.** With the boundary at the trench offset, the trench lies ~120–150 km offshore, and the cordillera (170 km from the boundary) lands ~50 km inland. The trench axis sits 30 km seaward of the boundary (`relief.rs`), so trench to crest is ~200 km. Andes check: the nearest volcanoes are ~240 km from the trench in southern Peru (search snippet, unverified).
- **Passive-margin escarpments: deferred.** Great escarpments (southern Africa, the Western Ghats, eastern Australia) are rift-flank features. Rolling them on a random share of coasts isn't physical (review). Revisit with the rift-shoulder profile once opening seaways carry a rift age.
- **Seafloor age** keeps the GDH1 law and the Parsons area–age calibration. Opening seaways put the ridge on the medial axis, so the oldest floor lies at the passive coasts (Atlantic pattern). The sediment wedge stays tied to `edge_m`.
- **Interior old belts** (noise) stay; merged "old sutures" add to them.
- **Seafloor age in closing seaways** counts from their own prescribed ridge (§1), not from a far-away one.

## 6. Conforming to the mask

- **Strict** (the default for imports): no sea-level bisection.
  - Land nodes: `h = max(h, 1 m)`. A finer coastal ramp would be moot: hypsometric calibration rewrites land heights by rank, so the lowest ~10% of land becomes 1–67 m whatever the input (review; `hypsometry.rs:86` floors at 1 m, the linear path at 1 mm, `solver.rs:714`).
  - Sea nodes: heights above −200 m are softly compressed into [−200, −1] m, so an arc or hotspot crest becomes a shallow bank, not a flat sheet at −1 m.
  - Land fraction = the mask's; the slider is disabled and shows the value.
  - **Ocean comes from the mask, not from the sign of h.** Today erosion resamples heights to its own grid (bilinear, possibly another resolution, jittered when `node_jitter` > 0) and sets `ocean = h0 ≤ 0` (`solver.rs:592`). A +1 m land node next to an −8 m shelf node puts the interpolated zero crossing ~10–20% of a cell off the mask coast, so coastal nodes flip. The snapshot therefore carries the mask (as a coverage function), and erosion evaluates it at its own node positions.
- **Islands allowed** (an option for imports): tectonic features (arcs, hotspots) may rise above sea level, but only beyond a guard distance from the mask coast (default 100 km). Drawn coasts never move. The land fraction is then the mask's plus the islands', and the slider stays disabled.
- **Procedural mode** (R3.1a): today's crust noise gives the continents. The coast still **emerges** from the heights through the sea-level bisection, so the land-fraction slider works as now. The plausibility pass classifies the coasts of the crust edge. A procedural strict mask (the noise thresholded straight to `land_fraction`, plus a coast-detail octave) is an R3.1b option, decided once the strict path exists. **Every procedural world still changes** (new plates).

## 7. Coastal bathymetry

`finish` swaps the `activity` noise for the `margin` class:
- **passive:** shelf mean ~88 km (`margin_w` varies it, ~40–250 km), slope, then the rise with the existing sediment wedge;
- **active:** shelf mean ~31 km, steep slope into the trench;
- **volcanic islands:** apron only.

Both shelf means are from Harris et al. 2014, *Marine Geology* 352, verified via the abstract: "continental shelf width of passive margins (88 km) is nearly three times that of active margins (31 km)".

Mid-ocean ridges come from step 3. Open-ocean fabric (transforms, fracture zones, abyssal hills) is R3b.

## Gates (invariants to test, not assumptions)

Unverified until a test says so (log: "n=0"):
1. **Coast survival (strict), two checks:**
   - erosion's `ocean` equals the mask at the erosion nodes (exact, by construction);
   - the committed height's sign matches the mask at the committed resolution, with every miss within one cell of the coast.

   Hypsometric calibration is safe (review: it floors land at 1 m / 1 mm). **Unchecked:** whether `lakes::compute` can leave a lake bed ≤ 0 on mask land. Below-sea-level land (Dead Sea, Caspian) is out of scope until ocean-from-mask has landed and is tested.
2. **Plausibility:**
   - ≥ 90% of active coast length has a trench within 300 km;
   - no trench within 400 km of a passive coast;
   - ridge nodes of opening seaways have `|d_A − d_B| / (d_A + d_B) < 0.15` (relaxed for the warp noise);
   - seafloor age rises toward passive coasts, and toward both rims of closing seaways;
   - motion fit: length-weighted RMS residual under a threshold set from the first runs; wrong-sign share reported;
   - same classification (seaway classes, segment cuts) at res 320 and 640.
3. **No regressions:**
   - `probe_bathymetry`: mean depth within 3,688 ± 60 m;
   - `probe_realism`: hypsometry within today's bands;
   - `probe_res_coherence`: height r ≥ 0.9 at 320 vs 640, same plate ≥ 95%.
4. **Cost:** tectonics time at Normal and Ultra, measured before and after. Expect +20–40% (3–5 extra distance passes; the solve itself is trivial).

## Tests that will change (need Matt's OK before editing)

- `rg-tectonics`:
  - `land_fraction_hits_target`: should still pass in R3.1a (the bisection stays); strict mode needs its own test.
  - `heights_are_plausible`: probably passes; re-check.
  - `plates_are_valid_rotations`: kept true by the `omega ≥ ε` clamp.
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
  - `tectonics.margin` (u8 per sample, committed);
  - `tectonics.convergence_mm_yr` (f32 per boundary sample, clamped; read by `uplift.rs`).
- Schema-4 files load as `mode: procedural`, with `margin` absent. The Margins layer is then hidden.

## Build phases (each gated)

R3.1 is split so that "did plausibility improve" is judged apart from "do coasts survive" (review).

- **R3.1a, the plausibility pass on procedural continents.**
  - Includes: coast classification, the fitted partition, prescribed polarity, the motion solve, the committed `convergence_mm_yr` read by uplift, ridges in closing seaways, class-driven shelves, the "active share" slider and the Margins layer.
  - The sea-level bisection stays.
  - Gate: the plausibility and no-regression gates, plus Matt's look in the app.
- **R3.1b, strict conformance.** The mask carried in the snapshot, ocean from the mask in erosion, the coast-survival gate, and the lake-bed check.
- **R3.2, import, placement, schema 5.**
  - PNG first (world equirect and patch), then SVG (Matt, s04), GeoJSON later.
  - UI: source picker, placement fields (centre, km/px, rotation), the strict / islands option.
- **R3.3, editing.** The brush on the globe (needs ray picking against the terrain, an M1c item), land/sea and active/passive strokes, regenerate.
- **Deferred:** GeoJSON, polygon rotation, passive-margin escarpments.
- Then **R3b**, open-ocean fabric, on R3's ridges and obliquity.

## Rejected alternatives

- **Mask imposed after tectonics** (override the crust, adjust uplift and sea level): plates ignore the coast, so trenches and ranges land in the wrong places. That's the problem R3 exists to fix.
- **Per-segment noise classification:** short active/passive alternations; nothing ties conjugate coasts together.
- **Random poles, rerolled until the classes match:** no guarantee, and slow. The least-squares solve is direct.
- **Plate partition by `grow_regions` with cost walls at trenches** (the draft): walls only delay fronts, so noise leaks plates across them, and there's no way to pre-label land (review).
- **Seaway = landmass pair:** merges oceans that the same two continents face on opposite sides (review).
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
7. Default active share **0.25**, and closing seaways get oceanic plates between active rims. Does the Pacific/Atlantic split read right for fantasy worlds, or should "all passive" (a young, rifting world) be a one-click preset?
