# RealmGenesis 3D — Feature Roadmap

This document lays out the feature direction for turning RealmGenesis from a
planetary-simulation viewer into a full fantasy map / planet creation tool.
Azgaar's Fantasy Map Generator is the reference point for scope; the goal is
not feature parity but complementarity — RealmGenesis is sphere-first and
simulation-first, and the roadmap leans into that.

**Where we stand vs. Azgaar.** RealmGenesis already leads on the planetary
side: true spherical Voronoi cells, plate tectonics with stress, hydraulic and
thermal erosion, latitude/tilt-aware climate, Dymaxion and multi-projection
export, GLB export. The gaps are in cartographic *presentation* (nothing on
the map is labeled or styled like a map), worldbuilding depth beyond
factions/provinces/towns, and a handful of simulation features that only make
sense on a sphere — which flat-map tools cannot do correctly at all.

Themes are ordered roughly by leverage. Within each theme, items are ordered
by suggested implementation order. File references point at the current code
that each feature would build on.

**Status key:** ✅ DONE · 🟡 PARTIAL (see note under the item) · ⬜ TODO.

As of 2026-08-21 (Sessions 18–20, merged to `main`): the
F2 ScreenOverlay migration now covers the **graticule, ocean currents,
roads/routes, contours, faction borders, and rivers**; only labels remain. Overlay
parallax is fixed and confirmed — the residual cause turned out to be a
one-frame lag in the render loop, not geometry. **D8 World Datum** was added and
scoped after auditing where heights are actually consumed.

As of 2026-08-12: the "Map identity" milestone (A1/A2/B1/B3), A4/A5, C1–C4,
E1/E2, and the D6 terrain overhaul are shipped. D7 part 1 (connected plates)
landed Session 9; the F1 shell is now the default route.

Note from Matt: While the ROADMAP isn't finished apart from UI stuff, don't be afraid to break "backwards compatibility" with existing seeds, no one else is using this project, and far from production ready anyway."

---

## A. Cartographic presentation

The biggest gap. A generated world currently *renders* but does not *read* as
a map: no names appear anywhere, and every view mode is a data visualization
rather than a cartographic style.

### A1. Map labels & typography  —  ✅ DONE
> _core shipped; per-glyph curved text still deferred_

Render names — factions, towns, rivers, lakes, geographic features — on the
2D map and the 3D globe. Names already exist in `CivData`
(`types.ts`) but are only visible in the Inspector and the 3D curved faction
labels; nothing else is ever labeled.

- 2D (`components/Map2D.tsx`): canvas text drawn in the composite pass, sized
  by importance (faction > ocean > capital > range/sea > town > lake > river),
  with greedy priority decluttering and zoom-based level of detail. Elongated
  features (rivers, ranges) get angled labels along their principal axis;
  per-glyph curved text is a stretch goal.
- 3D (`components/WorldViewer.tsx`): billboarded sprite labels with
  back-of-globe and distance culling. Factions keep the existing curved
  textured-mesh labels.
- Toggleable per label kind; included in PNG export (`utils/export.ts`).
- The Dymaxion point projector currently local to `Map2D.tsx` gets promoted
  into shared geo code so export projects label anchors identically
  (single-source-of-truth invariant, `utils/geo.ts`).

This is the single most transformative feature for "proper fantasy map" feel.

### A2. Offline name generation (namebases)  —  ✅ DONE

Seeded Markov-chain / syllable namebases so every world gets names with no
API key. Azgaar's namebase model is the reference: several built-in styles
(generic fantasy, Norse, Latin-imperial, desert, elven, …), each a word list
that trains a small char-level Markov model.

- Deterministic via fresh RNG side-streams from `utils/rng.ts`
  (`seed + '_names'`), so existing seeds keep their terrain and civ geometry.
- Hooks where placeholder names are currently assigned in `recalculateCivs` /
  `recalculateProvinces` (`utils/worldGen.ts`); geographic names key off the
  terrain seed so re-rolling civs doesn't rename mountains.
- Gemini (`services/gemini.ts`) demotes to an optional lore *enhancer* that
  overwrites the offline names when invoked — no longer a requirement for a
  named world.
- New `nameStyle` param + Civ-tab selector.

### A3. Map style system (fantasy rendering)  —  ✅ DONE (parchment)
> _shipped 2026-08-29 (Session 27b). Spec + plan in `docs/superpowers/`; reference doc `docs/map-styles.md`._

A style layer **orthogonal to `ViewMode`** — `viewMode` decides what the map
shows, a style decides how it is drawn, so parchment + political is a valid and
intended combination.

**Shipped:** one **Parchment** preset — paper tone and grain, hatched ocean,
strong ink coastlines, procedural relief glyphs (mountain / hill / forest /
conifer / dune / marsh) over a soft land-only hillshade, bare-paper land.
Selected from a **Map Style** control; Map2D, PNG export and SVG export all
drive the same pass list.

**Scope taken deliberately:** Map2D + exports only. The **3D globe is
excluded** — it paints per-cell vertex colours, so hatching and paper grain
would need a baked texture or custom shaders: roughly double the work, with M1
performance risk, for a view used to inspect rather than to publish.

**Architecture:** a style is an ordered list of passes written once against a
`Substrate` adapter with Canvas2D and SVG implementations, so the raster and
vector maps cannot drift. Glyph placement is a pure function of
`(cells, projection, widthPx)` — screen-space decisions kept out of drawing.
Full detail, including the mirror trap that affects **every** 2D render path,
is in `docs/map-styles.md`.

**Fill-policy rule (the non-obvious part):** continuous-ramp modes (`height`,
`temperature`, `moisture`, `population`) keep their own fill and suppress
glyphs. A ramp's entire information content *is* its fill, so bare paper would
render those modes blank.

**Preset backlog — logged, not built:**
- **Blueprint** — cyan paper, white line work, technical annotation, grid always on.
- **Ink & wash** — monochrome linework with a single wash colour for elevation.
- **Antique nautical** — rhumb lines, compass roses, sea monsters, heavy graticule.
- **Satellite-clean** — the current satellite look with better typography and no chrome.
- **Risk/boardgame** — flat saturated territory fills, thick black borders, province numbers.
- **Topographic** — contour-led (A4 contours exist), minimal fill, survey typography.

### A4. Hillshading & contour lines  —  ✅ DONE

Cheap wins from existing per-cell height: a hillshade pass (light-vector dot
cell normal) usable in any view, and elevation isolines. Both are derived
data — no pipeline changes.

### A5. Scale bar & great-circle measurement  —  ✅ DONE

`planetRadius` is already a parameter, so true geodesic distance is one
formula away. A ruler tool measuring great-circle distance between clicked
points (with an optional travel-time estimate reusing the Dijkstra terrain
costs from civ expansion in `utils/worldGen.ts`), plus a projection-aware
scale bar on exports. Flat-map tools fundamentally cannot do this correctly —
it's a headline geodesic feature.

### A6. Grid overlays  —  🟡 PARTIAL
> _3D + 2D graticule shipped; configurable density + geodesic hex grid TODO_

The 10° graticule exists on the globe; add configurable density, 2D-view
graticules, and a geodesic hex grid (icosahedral subdivision) overlay for
hex-crawl RPG use.

---

## B. Physical geography completeness

### B1. Lakes as first-class features  —  ✅ DONE

The Priority-Flood depression fill in `generateRivers()`
(`utils/worldGen.ts`) already computes a per-cell `waterLevel`; land cells
where `waterLevel > height` are exactly the filled depressions — currently
discarded. Surfacing them:

- `Lake` entity: member cells, surface level, outflow river or endorheic
  flag, fresh vs. salt (arid + hot basins become endorheic salt lakes).
- New `LAKE` / `SALT_LAKE` biomes; lake cells render as water in all views
  (`utils/colors.ts` is the single coloring function, so one change
  propagates to globe, Map2D, minimap, PNG, and GLB).
- Rivers terminate at lake entry and restart at the outflow.
- Civ integration: no capitals/towns/population on lakes; lake adjacency
  counts as coast.

### B2. River & lake editing / naming  —  ⬜ TODO

Extend the existing paint/edit system (`utils/paintUtils.ts`) to hydrology:
rename rivers and lakes, delete or redirect river segments, adjust lake
levels. Painting terrain over a lake clears it (snapshot semantics), with a
"resurvey" action to recompute hydrology-derived features.

### B3. Auto-detected & named geographic features  —  ✅ DONE

Cluster cells to detect mountain ranges (contiguous high elevation), peaks,
deserts, forests, oceans/seas (water bodies by size), and islands (small land
components); name them via A2 and store anchors (plus a principal axis for
elongated features) for labeling. All O(n) BFS over the existing cell graph —
fine at the 200k-cell cap. Surfaced in the Inspector and consumed by A1.

---

## C. Worldbuilding depth (civilization layer)

### C1. Cultures layer  —  ✅ DONE

Distinct from factions: culture regions spread by terrain affinity, each tied
to a namebase from A2. Factions then inherit (and can span) cultures —
enabling believable multi-ethnic states and consistent regional naming.

### C2. Religions layer  —  ✅ DONE

Folk religions derived per culture, plus organized religions spreading from
holy cities along the faction graph and (later) trade routes.

### C3. Roads & trade routes  —  ✅ DONE

A*/Dijkstra paths between towns over the existing terrain-cost model (reuse
the cost logic in `recalculateCivs`), plus sea routes between ports. Route
connectivity should feed back into town importance/population.

### C4. Markers & notes (POIs)  —  ✅ DONE

User-placed pins (dungeon, ruin, battlefield, portal, …) with free-text
notes, persisted in the JSON save; optional Gemini lore per marker.

### C5. Editor completeness  —  🟡 PARTIAL
> _merge/rename/capital-relocation done; faction split deferred_

Already tracked in `HANDOFF.md` next tasks: merge/split factions, province
management, bulk rename, town/capital relocation.

### C5b. Territorial waters are dead at default params  —  ✅ DONE
> _fixed 2026-08-22 (`c9b70b7`); peaks half verified not-reproduced 2026-08-22 (S25)_

**The bug:** `recalculateCivs` charged **40** for one water step against a total
territorial-water budget of **7.5**, so no water cell was ever claimed at default
params and `territorialWaters` did nothing. Interior lakes were folded into the
same test, making them impassable walls rather than claimable territory — high
mountains starved separately against the `cost > 200` cap.

**The fix (`c9b70b7`):** territorial waters are now counted in **ocean steps from
the coast** (`maxWaterSteps = round(territorialWaters * 12)`), not a cost budget a
single step already blows; reaching land resets the allowance, so a realm
island-hops one strait at a time. **Lakes count as land** (`isOcean = height <
seaLevel && !isLakeCell`), so an enclosed lake belongs to its realm instead of
walling it off. The `cost > 200` frontier cap is **gone**, so land expansion runs
to completion. Guarded by `tests/civClaim.test.ts` (4 cases).

**The peaks half — verified not-reproduced (S25).** The cap removal was a
*reasoned* fix for stranded high peaks that never reproduced at 3k/30k; Matt saw
them at 200k. Probed 5 seeds up to 200k cells: **`reachableGaps = 0` in every run**
(no unclaimed land ever borders claimed land — expansion can no longer die partway
across a landmass). Unclaimed peaks appear only on isolated landmasses beyond the
water-step allowance, which is correct. Structural guarantee: no cap → no
partway-death → any unclaimed peak sits on an unreachable mass.

**Seed compatibility broken deliberately** — Matt authorised it in ROADMAP:34.

### C6. Diplomacy & relations (later)  —  ⬜ TODO

Inter-faction relations, vassals, alliances, war status — mostly flavor and
map overlays, cheap once C1–C3 exist.

---

## D. Planet-scale simulation (the geodesic differentiators)

Features that exploit the sphere-native architecture — the category where
RealmGenesis can do things no flat-map tool can.

### D1. Seasonal cycle  —  ✅ DONE
> _shipped Session 14; wind/moisture-monsoon variant deferred to ENGINEERING-NOTES_

`axialTilt` was a static climate offset; it's now the amplitude of a seasonal
excursion. A `season` slider (orbital position, neutral at 0.5) shifts
temperature, snow line, and biome edges through the year — render-only, no
regeneration. Stored `cell.temperature` is the tilt-dependent orbit-averaged
annual mean (keeps `axialTilt` a live generation param); the per-season excursion
is a closed-form per-cell formula (`utils/seasons.ts`) applied in the color path,
anchored to the equinox so neutral reproduces the canonical annual-mean world.
Biomes shift for display; civs/export/labels stay pinned to the annual world.
Wind + moisture stay annual (the monsoon variant is deferred). Spec:
`docs/superpowers/specs/2026-08-14-d1-seasonal-cycle-design.md`.

### D2. Ocean currents  —  ✅ DONE
> _shipped Session 15; annual/steady-state, semi-simulated gyres; seasonal current reversal deferred_

Currents deflected by continents and Coriolis, feeding the existing 8-pass
moisture transport and temperature model (warm currents moderating coastal
climates — think Gulf Stream). A major realism upgrade to climate.

**Shipped: `currentStrength`** (0–2, default-on 1.0). A fixed-pass, fixed-order
relaxation (`utils/currents.ts`, deterministic, no RNG, no Poisson solve — the
divergence-free tier would break determinism + the browser budget) seeds ocean
velocity from wind stress, then relaxes with Coriolis deflection (∝ sin lat),
advective smoothing, and net-land-normal boundary tangency → emergent gyres +
western-boundary currents. A heat-advection pass yields an SST anomaly that
**moderates coastal temperature** (ocean cells + a 1-ring land coastal blend) and
**boosts warm-current evaporation** into the 8-pass moisture transport. Cascades
through D1 biomes + D3 sea-ice for free. **0 = byte-identical escape hatch**
(stage short-circuits). Coupling into **both temperature and moisture** pulls part
of the deferred D1 seasonal-moisture work forward — see the composability seam in
ENGINEERING-NOTES. Spec:
`docs/superpowers/specs/2026-08-15-d2-ocean-currents-design.md`.

### D3. Ice caps & glaciers  —  ✅ DONE
> _sea-ice shipped Session 14; sea-level coupling deferred to ENGINEERING-NOTES_

Seasonal sea-ice: open-ocean cells below the seawater freeze point (−2°C, a
physical constant) render as ice in the satellite + biome views — a render
overlay reading D1's seasonal temperature, so polar caps form at neutral and
migrate hemisphere-to-hemisphere with the season. `cell.biome` stays `OCEAN`
(no civ/nav impact). Land ice stays with the existing `ICE_CAP` biome + snow.
Sea-level coupling (ice ↔ coastline) deferred (it's a generation-stage change,
not a render overlay). Spec: `docs/superpowers/specs/2026-08-15-d3-sea-ice-design.md`.

### D4. Regional zoom / submap generation  —  ⬜ TODO

Re-run generation at higher point density inside a selected spherical cap for
local/regional maps — Azgaar's "submap" concept, but projection-correct, with
parent-world values as boundary constraints.
Idea: A "base resolution" inside one world, but certain spots can be selected with a bounding box
to upscale/regenerate cells at much higher densities

### D5. Planetary parameters  —  🟡 PARTIAL
> _host star class shipped Session 14; day length / gravity / moons deferred_

Day length, gravity, star class, moons — mostly lore/export metadata at
first, but D1 (seasons) and tides give them mechanical hooks over time.

**Shipped: host star class** (`starClass`, O–M, default G). Scales global
insolation → temperature (Kelvin-space Stefan-Boltzmann scale), which cascades
through D1's biomes and D3's sea-ice — a K-class world ices over, an F-class one
bakes. G is a byte-identical no-op. Surfaced in Controls, lore, and save.
**Deferred** (no principled hook yet — see ENGINEERING-NOTES): day length (needs
D2/Coriolis), gravity (would be a relief fudge duplicating `mountainHeight`),
moons (needs a tide model). Spec:
`docs/superpowers/specs/2026-08-15-d5-star-class-design.md`.

### D6.  (added by Matt/Maintainer)  —  ✅ DONE
> _V3 terrain model shipped & live (Session 8/9)_

Another overhaul of terrain generaton algorithm to make it more realistic
and get rid of seam lines on plate boundaries, and more detailed heightmap
rendering and calculation without increasing cell count

### D10. Flat sea bed — no bathymetric relief  —  ✅ DONE
> _root-caused + fixed 2026-08-28 (Session 27). Cause was NOT the missing slider._

Matt reported the sea bed showed little height variation next to land at any
resolution, and that Sea Depth / Trench Depth did not change it.

The sea bed was **not** flat in the data — ocean depth already spanned ~8,000 m
against land's ~2,700 m. The deficit was purely high-frequency. Measured across
3 seeds at 20k, per-neighbour relief in normalized height: **ocean 0.0153 vs land
0.0250, a 1.64x texture gap.** Big smooth swells, no bumps. An aggregate depth
histogram cannot see this — the diagnosis needed a two-scale relief metric
(neighbour delta vs neighbour-of-neighbour).

**Root causes (two, neither of them the missing control):**
1. **`applyThermalErosion` had no sea-level check** (`utils/worldGen.ts`). It is a
   smoothing operator — any slope above `talus` sheds to its lowest neighbour —
   and it ran on every cell, planing the sea bed flat regardless of what
   generation produced. Ablation: with erosion off, a relief sweep moved ocean
   texture 0.0117→0.0140; with erosion on, 0.0162→0.0159, fully erased. Talus
   erosion models freeze-thaw and gravity scree, a **subaerial** process;
   submarine slopes hold far steeper angles. `thermalSteps` scales with
   `sqrt(points)`, which is the "at any resolution" part of the report.
2. **`computeShadeMap` short-circuited every water cell to `shade = 1.0`**
   (`utils/shading.ts`), so the sea floor had zero relief shading in Map2D and
   every export whatever the bathymetry said.

The retired `seafloorDetail` knob would not have helped: it was **inverted** at
the display site (`noiseInfluence *= 1 - 0.65 * seafloorDetail`), so turning
"detail" up damped fine noise harder. It grew macro swells while flattening
texture.

**Fix:**
- Split the talus constant — subaerial `0.008`, submarine `0.12`. Ocean cells stay
  donors, so nothing piles at the coast (verified: no deposition rim).
- Water hillshades off its **water** neighbours only (a land neighbour would put
  the full land/ocean step into the gradient and rim every coastline). Land
  untouched, so existing maps do not shift.
- New **`seafloorRelief`** param (0–2.0, default 1.0) applied in a new **Stage 9c**
  — after normalization, after erosion, at display resolution. 4-octave fBm,
  zero-mean, shelf-tapered, clamped below `seaLevel` so **land fraction is
  invariant by construction**. Distinct from `seafloorDepth`, which is a
  mean-depth datum and cannot change roughness.

| seafloorRelief | ocean texture | land texture | ratio | land cells |
|----------------|---------------|--------------|-------|------------|
| 0   | 0.0167 | 0.0250 | 1.50x | 4,558 |
| 1.0 (default) | 0.0253 | 0.0251 | **0.99x** | 4,565 |
| 2.0 | 0.0373 | 0.0251 | 0.67x | 4,593 |

Default 1.0 means the sea bed carries the same texture as land. Holds at 80k
(1.53x → 0.96x). Seed movement was authorized by Matt; there is no byte-identical
hatch.

**Rejected:** driving the param from the tectonicsV3 display-noise damping factor.
Implemented, measured, reverted — un-damping that land-tuned noise pushed the
global pre-normalization minimum down, and renormalization against a fixed
`seaLevel` inflated land-cell count ~52% (4,410 → 6,692). It also saturated at
relief 1.0, so half the slider did nothing. See HANDOFF S27 for the full record.

### D9. Pangea bias — terrain concentrates into supercontinents  —  ✅ DONE
> _root-caused + fixed 2026-08-22 (Session 23). Cause was NOT plate seeding._

Since the tectonic-plate improvements, generated habitable terrain gathered into
one part of the sphere — pangea and superoceans — and the existing sliders did
not break it up.

**Root cause (the plate-seeding suspects were all wrong):** `seedCrustField`
(`utils/crust.ts`) decided continental vs oceanic crust by thresholding a
**single** simplex octave at base frequency **0.3** on the unit sphere. At that
scale the whole planet sits inside one noise lobe, so the threshold split the
sphere into one land cap and one ocean cap. Measured on the real field (10k macro
points, 5 seeds): **one connected component holding 100% of land**, mean-position
clump metric **0.74**. The crust field is seeded independently of plates
(`tectonicsV3.ts:582`), so plate count / jitter / elongation never touched it.

**Fix:** the field is now **fractal (fBm) at base frequency ~1.0**, with
per-`landStyle` parameters, so continental crust breaks into several separate
masses at continental scale. Measured (mean of 5 seeds):

| landStyle | masses | largest mass | land frac |
|-----------|--------|--------------|-----------|
| Continents (default, + Custom) | ~7 | ~47% | ~35% |
| Pangea (intentional) | 1 dominant + satellite | ~86% | ~31% |
| Archipelago | ~21 | ~21% | ~35% |
| Islands | ~27 | ~10% | ~34% |

Continents now matches Earth (~6 masses, largest ~40–57%). `Pangea` is preserved
as a deliberate user choice. Guarded by `tests/crustDistribution.test.ts`.

**No escape hatch, deliberately.** Unlike D2/D5 (where prior output was fine and
the feature was additive), here the prior output *is* the defect — a "reproduce
pangea" hatch would only preserve the bug. ROADMAP:34 authorizes seed breakage.
Rejected alternative: a `crustFractalOctaves=1` flag reproducing the old field.

### D8. World datum (units)  —  🟡 PARTIAL
> _surfaced Session 19d; scoped Session 19e; **D8a shipped Session 25** (`21aaded`); D8b still TODO_

Heights are a normalized 0–1 field with no real-world meaning, while horizontal
distance is already genuine (A5 measures great-circle **kilometres** from
`planetRadius`). **The app therefore measures horizontally in real units and
vertically in percent** — that inconsistency is the core argument, not DEM export.

It splits into two tiers with very different costs.

#### D8a — presentation datum (cheap, no seed changes)  —  ✅ DONE (S25, `21aaded`)

`maxElevationM` (default 9000) is a **display-only `WorldParams` key** — like
`planetRadius`/`season`, user-adjustable but read by nothing in generation, so
seeds stay byte-identical. `utils/datum.ts` is the single source:
`elevationMetres`/`formatElevation` scale **above sea level** against the FIXED
maximum (comparable between worlds); `MAX_DEPTH_M` is a fixed constant because
`seafloorDepth` already owns generation-side depth. A `season`-style sync effect
pushes the live value into `world.params` so the slider takes effect with no
regenerate. Guarded by `tests/datum.test.ts`.

- Inspector `Elev: 70%` → `1,000 m` (S25b: a quadratic hypsometric curve, `frac^2`,
  replaced the linear map — linear reported a ~2 km median land elevation vs Earth's
  ~840 m. Curved matches Earth: mean 824 m, 72% under 1 km. `HYPSOMETRIC_EXPONENT`
  in `utils/datum.ts`; verified with `scripts/queryWorld.mjs hypsometry`). **DONE.**
- Contour labels → metres via `contourLabel()`. Still DORMANT (labels pulled S19e);
  updated as the documented single change point so it never re-earns the drift. **DONE.**
- GeoJSON export: `height` **replaced** with rounded metres (Matt's call — the
  vertical is now genuine, matching the already-geodesic lon/lat). **DONE.**
- Lore prompts: **deferred** — `services/gemini.ts` has no existing terrain-height
  text to convert; adding it is net-new prompt work, out of D8a scope.

**Caveat that decides the design:** if the datum is *derived* from
`mountainHeight`, the same cell reports a different altitude when that slider
moves, so values are not comparable between worlds. It must be a fixed maximum
that `mountainHeight` distributes terrain within.

#### D8b — simulation coupling (`physicalClimate`)  —  ✅ DONE (merged to `main` `16ee4ce`, 2026-08-28; rain-shadow retune `da2dccf`; pending only Matt visual sign-off)

**Shipped (S26, branch `d8b-climate-coupling`).** `physicalClimate` (boolean,
default **ON**) gates two sites in `utils/worldGen.ts`; off = **byte-identical**
to pre-D8b `main` (verified per-cell: all generation output identical, only the
`params` echo differs). New physical constants live in `utils/datum.ts`.

- **Lapse rate** — grounded at **6.5 °C/km** on real datum metres (`LAPSE_RATE_C_PER_KM`).
  ICE_CAP stayed ~7% (the `frac^2` curve keeps land low, so only real peaks cool).
- **Orographic rain shadow** — scaled by real barrier metres (land-side above-sea,
  ocean → 0), tuned constants in `utils/datum.ts`.
- **Snow line** — emergent for free; `determineBiome` unchanged.
- **Moisture retune + rain-shadow fix** (folded in; retuned `da2dccf`): land
  `moisture<0.15` share **~35.4%** across 3 seeds (in the 30-36% band); steppe
  ~18-25% (below pre-D8b 26.5-29%), grassland healthy, no collapse. Final
  constants: windward 0.85, leeward floor 0.5 / per-km 0.3. An earlier tune hit
  32.4% by nearly disabling leeward drying, which KILLED rain shadows; the retune
  restores them — verified with a windward/leeward moisture-contrast metric (0.08
  → 0.135). Deeper root cause (base moisture under-delivers inland: spec §6, the
  8-pass transport's land decay + fixed pass count) is a follow-up the orographic
  knob can locate but not fix.
- **`maxElevationM` is now a GENERATION param** (Matt's call) — it drives lapse +
  orographic, so it regenerates to apply; its D8a live display-only sync was removed.
- One existing test (`lakes.test.ts` salt-lake seed) pinned to `physicalClimate:false`:
  the retune wets s17's basin so it freshens instead of forming an arid salt lake —
  the authorized behavior change, testing SALT_LAKE hydrology needs the classic climate.

Spec: `docs/superpowers/specs/2026-08-23-d8b-climate-coupling-design.md`; plan:
`docs/superpowers/plans/2026-08-28-d8b-climate-coupling.md`. Deferred (out of scope,
recorded): air-temperature orographic factor, volcanic-vs-temperature decoupling,
D5 gravity, D3 sea-level coupling.

**Original design notes (kept for rationale):**

**Decisions locked S25b (2026-08-23):** flag `physicalClimate`, default **ON**
(off = byte-identical old formulas); **all 3 couplings under it**; datum pick
**physical lapse 6.5 °C/km + datum 9000 + the `frac^2` curve** — measured stable
(ICE_CAP 7.0→6.7%) because the curve keeps land low. Matt accepts civ-layout
movement on existing seeds. **Moisture dryness folds in here:** 42% of land is
moisture <0.15 (Earth ~33% arid+semiarid) — the real cause of steppe dominance —
so retune the inland rain-shadow when rewriting the orographic pass, once. The
GRASSLAND biome added S25b will populate as interiors get wetter. Full context:
HANDOFF "D8b — IN DESIGN".

The larger payoff, because several tuned magic constants are really physical
quantities wearing normalized clothes:

- **Lapse rate.** `utils/worldGen.ts` does `temp -= elevation * 60`. That is a
  lapse rate with an invented 60. In metres it is ~6.5 °C/km — grounded, and it
  stops fighting `mountainHeight` (raise peaks today and you also silently
  refrigerate them by an unrelated amount).
- **Orographic precipitation.** The 8-pass moisture transport uses
  `heightDiff > 0.02 → carry *= 1.5`, `< -0.02 → *= 0.2`. Rain shadow strength
  physically depends on barrier height and air temperature, both of which a
  datum supplies.
- **Snow line.** `determineBiome` decides ice and tundra from temperature alone,
  and volcanic from `landH > 0.85`. A real snow line is an *altitude* that
  varies with latitude.
- **Unblocks D5 gravity**, which was deferred precisely because it "would be a
  relief fudge duplicating `mountainHeight`". With a datum it gets a principled
  hook: isostasy caps how high mountains can stand (roughly ∝ 1/g), so gravity
  shapes maximum elevation instead of duplicating a slider.
- **Unblocks D3 sea-level coupling**, deferred as a generation-stage change —
  "the caps melting raises sea level N metres" is only expressible with a datum.

**Cost:** every item above is a generation input, so it changes output for
existing seeds. The project has treated that as a hard line before (D2 ships
`currentStrength = 0` as a byte-identical escape hatch; D5's G-class star is an
exact no-op). D8b should follow the same discipline or it invalidates saved
worlds.

**Verdict:** D8a is close to free and fixes a real internal inconsistency. D8b is
a genuine simulation feature whose benefit is grounding three tuned constants and
unblocking two deferred items — worth doing, but scoped and escape-hatched, not
folded into D8a.

### D7. More realistic tectonic plates  —  ✅ DONE (heuristic tier); Cortial rebuild = deliberate NO-GO
> _part 1 (connected plates) Session 9; part 2 (seafloor age→bathymetry + microplates) Session 10; part 3 (plate-shape polish) Session 12. The Cortial boundary-curve rebuild was evaluated and declined — see below._

D6 set the foundation; parts 1–3 completed the heuristic-tier plate model:
- **Part 1 (S9):** multi-source Dijkstra region-growth → 0 enclaves/exclaves by construction.
- **Part 2 (S10):** seafloor age→GDH1 bathymetry, shear-driven microplates.
- **Part 3 (S12):** the real de-blob levers turned out to be **plateJitter** (size/position
  variety) and **boundaryRoughness** (jagged boundaries), both range-extended to 0–3 and
  defaulted to 1.5; plus band/chain seeding (`plateElongation`, mild) and seafloor age noise.

**Cortial boundary-curve rebuild — NO-GO** (both advisors + render agreed): it would destroy
the 0-exclave-by-construction invariant for weeks of sphere geometric-robustness work, and
the research report's own recommendation for a browser budget is the heuristic tier already
shipped. **Shelved (not abandoned)** next levers, recorded in HANDOFF Session 12 for a future
`docs/ENGINEERING-NOTES.md`: anisotropic Dijkstra growth (real shape lever), transform-edge
fracture, and the Cortial path itself. Also open: fine-mesh (display) plate connectivity is
not guaranteed by the macro→fine downsample (pre-existing).

---

## E. Interoperability & export

### E1. SVG export  —  ✅ DONE

Vector export of coastlines, borders, rivers, and labels for post-processing
in Inkscape/Illustrator — a key Azgaar strength and a frequent request for
any mapping tool. The cached GeoJSON cell geometry (`WorldData.geoJson`) plus
d3-geo path generation makes this tractable.

### E2. GeoJSON export  —  ✅ DONE

Cells, coastlines, borders, and rivers as GeoJSON. The data is genuinely
geodesic (lon/lat on a sphere), so exports open cleanly in QGIS and web-GIS
tooling with zero fudging — worth advertising.

### E3. Azgaar `.map` import (stretch)  —  ⬜ TODO
> _stretch goal_

Project an Azgaar flat map onto the sphere (equirectangular assumption,
re-tessellate onto the cell graph). Lossy by nature; stretch goal.

### E4. Blender-accurate UV mapping support for Dymaxion  —  ✅ DONE
> _verified against live Blender 2026-08-22; net was already exact, default fixed_

Original note: "I recall importing Blender Icosphere model data to tune the
Dymaxion export UVs, but the current 2D dymaxion probably isn't exact to what
Blender shows."

**Outcome:** it WAS exact. Dumped the default icosphere's UV + vertex data live
from Blender 5.1 and compared: `buildBlenderNet` matched to 5 decimals (the Feb
extraction held). The export rasterizes at `px = u·W, py = (1−v)·H` on a square
canvas — Blender's own UV sampling formula — so a Blender-layout PNG drops onto
the icosphere with zero tweaking (confirmed by Matt: "drops perfectly"). The
black upper band is the icosphere's unused vertical UV (v only reaches 0.472),
not a defect.

The only real bug was the **default layout**: everything defaulted to the
`classic` staircase, so users exported a net that can never map onto the
icosphere. Default is now the `blender` sawtooth (view, preview, export). Ground
truth is frozen in `tests/fixtures/dymaxionBlenderUV.json` +
`tests/dymaxionBlenderNet.test.ts`; full writeup in `docs/dymaxion.md`.

---

## F. Frontend/Rendering/UI Overhaul

# F1. Redesigning the UI/Frontend  —  🟡 PARTIAL
> _desktop foundational done + F1 shell now the default route (Session 9); F1b brand pass + mobile polish pending_
The UI has become a mess in an attempt to add more features, especially on mobile.
A full redesign and rearchitecture is warranted here.
Can come *before or alongside* D6 with the roadmap in mind.

# F2. 3D Mode Presentation  —  ✅ DONE (overlay migration); permanent-3D items noted below
> _ScreenOverlay foundation + ocean-current viz + graticule Session 16; smooth-globe + graticule drape Sessions 17-18; roads/routes Session 19; contours Session 19b; faction borders Session 20; rivers Session 21; point labels Session 22; ruler Session 23; **Dymaxion cage Session 24 (last named tenant)**. Selection/highlight rings declined in writing (S23). Faction labels stay 3D permanently by design. **The migration is complete.**_
Part of redesign is figuring out how the planet is presented; overlays like borders, rivers and roads and routes and the lat/lon grid are also 3D objects, not 2D overlays simply composited over the 3D globe, which affects visibility and accuracy. Or maybe make the globe entirely smooth by default, instead of applying height per cell. Perhaps 3D mode should more like Google Earth Pro in this respect.

- **Ocean-current visualization** (from D2): a currents overlay drawing the
  `computeOceanCurrents` velocity field on the globe and 2D map, with warm/cold tint
  from the SST anomaly. Static arrows for v1 (animated particle advection deferred —
  thermals on the M1 Air). The velocity field is computed at generation but currently
  discarded; F2 persists it optionally on `WorldData`.

- **`ScreenOverlay` layer (the foundation).** Today every globe overlay (borders,
  rivers, graticule, Dymaxion, rulers, selection) is a *physical R3F 3D object* —
  which is exactly the "affects visibility and accuracy" problem above. F2 introduces
  a generic screen-space overlay layer: a canvas sibling to the WebGL canvas that
  projects visible cells (analytic horizon test, not depth-buffer readback) and draws
  overlays in pure 2D, occluding the far hemisphere without physical geometry on the
  globe. v1 lands **two tenants** to prove the abstraction generalizes: the currents
  field **and the graticule** migrated off its 3D `lineSegments`.
  **Migrated so far:** currents + graticule (S16), roads/routes (S19), contours
  (S19b), faction borders (S20), rivers (S21), point labels (S22, on branch
  `f2-labels-tenant`, verified but unmerged).

  **Status of the last three (updated 2026-08-22, S23):**
  - **`RulerArc`** — ✅ MIGRATED. Now the `ruler` ScreenOverlay tenant
    (`components/overlays/tenants.ts` `drawRulerTenant`), limb-broken polyline +
    endpoint dots at the fixed ruler radius. Guarded by `tests/rulerTenant.test.ts`.
  - **Cell highlight + selection ring** — ⛔ DECLINED IN WRITING (Matt's call,
    S23). Unlike Dymaxion's fixed 1.12 float, these already drape via
    `displayRadius` (a tiny 0.012 lift), so they barely have the parallax the seam
    removes — and a selection ring should read as *attached* to the selected
    cell's geometry, not composited flat over it. Flattening them is a downgrade,
    the way `CurvedFactionLabel` is. Applies to `CellSelectionOverlay`,
    `CellHighlightOutline`, and `BrushRing` together. Leave them 3D.
  - **`DymaxionOverlay`** — ✅ MIGRATED (S24). Now the `dymaxion` ScreenOverlay
    tenant (`drawDymaxionTenant`), **edges only, amber faces dropped** as Matt
    chose. Geometry (12 unit verts + 30 edges, exact `THREE.IcosahedronGeometry`
    orientation) lives in `utils/dymaxionCage.ts`; `cageEdges(settings)` rotates
    it by the same YXZ euler the old 3D cage used. Fixed r = 1.12 in both globe
    modes (a reference frame, not draped). Guarded by
    `tests/dymaxionTenant.test.ts`. The 3D `DymaxionOverlay` component + JSX are
    deleted.

    **Occlusion — the S23 endpoint-only plan was REVERSED (refuted).** The S23
    brief said to cull each edge on its two r = 1.12 endpoints, reasoning that a
    per-sample test would wrongly cull the chord's middle. That reasoning was
    inverted: `isVisible` tests each point against a sphere of the point's OWN
    radius, so the mid-chord (r ≈ 0.95) is MORE permissive than its endpoints,
    not less. Endpoint-only culling drops most of the cage — measured 3.6 / 30
    edges at camDist 2.5, 0.6 / 30 zoomed to 1.5 (the cage nearly vanishes). The
    tenant instead samples each chord (16 pts) and breaks the polyline at the
    horizon — the routes/ruler idiom — which clips edges at the limb (~13 / 30 at
    2.5). Perspective preserves straight lines, so collinear samples still draw
    the exact straight edge. "Back edges faint" stays out of scope (it needs a
    non-culling projector variant on the seam).

  **F2's overlay migration is COMPLETE.** Every named 3D overlay has migrated or
  been declined in writing: currents, graticule, routes, contours, borders,
  rivers, labels, ruler, and the Dymaxion cage are tenants; the selection/
  highlight rings are declined (S23). Faction labels (`CurvedFactionLabel`) stay
  3D permanently — Canvas2D cannot reproduce curved textured meshes without the
  per-glyph text A1 deferred; flattening them is a downgrade, not a migration.
  `CityMarkers`, `MarkerPins`, and `TiltAxisLine` were never in the F2 named
  scope and stay 3D (migrating them would fix the accepted overpaint nit as a
  side effect — a separate deliberate decision, not silent inclusion).

  **Even then the queue does NOT fully empty.** Faction labels stay 3D permanently:
  `CurvedFactionLabel` renders curved textured meshes following the sphere's
  curvature, which Canvas2D cannot reproduce without the per-glyph curved text A1
  deferred. Flattening them is a downgrade, not a migration — do not "finish the
  job" by removing them.

  Every tenant must render correctly in BOTH globe modes, keyed off `smoothGlobe`:
  flat on the unit sphere when smooth, draped over relief when raised. A tenant's
  radius must equal the terrain mesh's — `displayRadius(cell.height, smooth)`, with
  the height RAW — and must assert that equality in a test. "Parallax-free by
  construction" without such a test is exactly how the S18 graticule bug shipped.

  Matching radii is necessary but NOT sufficient: the overlay must also read the
  same FRAME as the renderer. See `ScreenOverlay`'s forced matrix update and
  `<GlobeSpin/>`'s mount position — three sessions of radius fixes chased a
  symptom whose real cause was a one-frame lag in the render loop.

# F3. True 2D vector map  —  ⬜ TODO
Make it a true vector map like most web mapping apps, but keep it as optimized as possible

# F4. Performance Optimizations  —  ⬜ TODO
Self-explanatory, optimize wherever while keeping visual fidelity. Can come last, but better to make efficient renderers/frontends and code wherever possible.


## Recommended first milestone — "Map identity"  —  ✅ DONE
> _A1 + A2 + B1 + B3 all shipped (pre-D6 tier, Sessions 3–4)._

**A1 + A2 + B1 + B3: labels, offline namebases, lakes, and named geographic
features.** Together these turn the simulation viewer into something that
reads as a *map*, and none of them require new simulation machinery — lakes
fall out of the existing depression fill, and labels/names consume data the
pipeline already computes or can derive in one pass.

A detailed implementation design exists for this milestone. Summary:

- **Phasing** (~2–3 weeks total): (1) namebase engine + civ naming hooks;
  (2) hydrology refactor + `Lake` entity + lake biomes/colors; (3) feature
  detection stage + geographic naming + Inspector integration; (4) label
  model + 2D/3D renderers + declutter + export labels. Phases 1 and 2 are
  independent and can land as separate PRs.
- **Data model**: `Lake`, `GeoFeature`, `RiverInfo` types; `Cell.lakeId`;
  `WorldParams.nameStyle`; `LAKE`/`SALT_LAKE` biomes; JSON save schema bumps
  to v1.5 with a names-only `geoNames` container (geometry is never
  serialized — it regenerates deterministically, matching the existing
  save philosophy in `utils/export.ts`).
- **Determinism rules**: every new random draw uses a fresh RNG side-stream
  (`new RNG(seed + '_<purpose>')`) so existing seeds keep their terrain and
  civ geometry byte-identical; geographic names key off the terrain seed,
  civ names off `civSeed`.
- **Top risks**: label decluttering complexity (scoped to greedy
  priority + bbox rejection; curved text deferred), label paths crossing
  Dymaxion net seams (straight anchored labels only there), and 200k-cell
  performance (O(n) typed-array clustering, labels drawn in the cheap
  composite pass, allocation-free 3D culling).

## Suggested order after that

1. **A3 + A4** (map styles, hillshading/contours) — presentation payoff
   compounds with labels.
2. **A5 + E1/E2** (geodesic ruler, SVG/GeoJSON export) — small, high-value,
   and differentiating.
3. **C1 → C3** (cultures → religions → routes) — worldbuilding depth, in
   dependency order.
4. **D1 → D2** (seasons, currents) — the sphere-native simulation showcase.
5. **D4** (submaps) and the rest as demand dictates.
