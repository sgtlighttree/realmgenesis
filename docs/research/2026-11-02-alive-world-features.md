# Alive-world layer: research report

- **Date:** 2026-11-02
- **Scope:** the features that make a generated planet feel inhabited (population, cultures, states, settlements, roads, religions, provinces, names, military, markers, diplomacy, history, heraldry, labels), how Azgaar's Fantasy Map Generator (FMG) and other open projects do them, and how each would map onto our cube-sphere world.
- **Status:** research only. No code changed, no other file touched, nothing committed or pushed.

---

## 0. Summary

**Top 5 to build first** (in dependency order):

1. **Population and habitability.** One pass over the cells. It drives almost everything else.
2. **Burgs (settlements, capitals, ports).** Placed by suitability with a minimum spacing.
3. **Names.** Our own generator, with invented syllable inventories per culture. No real-world name lists.
4. **Cultures.** A cost-weighted flood fill from culture centres. This is the expansion engine that states and religions reuse.
5. **States and borders.** The same flood fill, seeded at capitals.

**Next tier:** labels (they arrive with the vector layer, M7), markers (cheap and physically motivated), a minimal history (founding years and a short event list), roads and trade routes, provinces, religions.

**Later:** military, diplomacy, full chronicle, heraldry.

**Licence red flags:**

- **Do not copy** `watabou/TownGeneratorOS` or `raguilar011095/planet_heightmap_generation`. Both are GPL-3.0.
- **Armoria** declares MIT in its README and package.json, but has no LICENSE file in any of its 1,154 commits. Its complex charge art is CC BY-NC 3.0. FMG bundles that same art. **Never copy any charge art.**
- **Name lists with unknown provenance** (FMG's "real-world" name bases, langgen's sample lists) must not be copied.
- Our climate research doc (`docs/research/2026-11-02-climate-tier-b.md`) calls `planet_heightmap_generation` "licence unverified". It is now verified as GPL-3.0. That doc needs a correction (see §8).

---

## 1. Method and limits

- **Repos:** I shallow-cloned each public repo with `git` (the GitHub API was not enabled for these repos in this session). I read each LICENSE file from the clone. Commit hashes are in §2.
- **FMG source:** I read `Azgaar/Fantasy-Map-Generator` at commit `d279257` (package version 1.156.0). Its `docs/wiki/` folder explains the rules in plain words.
- **Not reachable from this session.** WebFetch could not resolve `www.naturalearthdata.com`, `www.geonames.org` or `www.dwarffortresswiki.org`. The session reports that these hosts must be added to its allowed network domains. Claims that depend on them come from search snippets and are marked **secondary** or **unverified**.
- **Not checked online.** Dwarf Fortress's own terms, Wonderdraft, Inkarnate, HYDE, the Beck Köppen dataset licence, OpenMoji, game-icons.net and JSZip. Each one is listed in §8 with the URL to check.

---

## 2. Licences of the projects checked

| Project | Commit | Licence (where stated) | May we port code? | Source |
|---|---|---|---|---|
| Azgaar/Fantasy-Map-Generator (FMG) | `d279257` | **MIT.** LICENSE: "Copyright 2017-2024 Max Haniyeu (Azgaar)" | Yes, keeping the MIT notice. **Not** its charge art (see §7.4) | https://github.com/Azgaar/Fantasy-Map-Generator/blob/master/LICENSE |
| Azgaar/Armoria (heraldry) | `b0bb9b8` | **MIT declared** in `package.json` and the README "License" section. **No LICENSE file** in any of 1,154 commits (full-history check) | Not yet. Ask the author for a LICENSE file first. **Never** its charge art: the README says complex charges are from WappenWiki under CC BY-NC 3.0 | https://github.com/Azgaar/Armoria (README "License"; `package.json`) |
| mewo2/terrain | `a2511d6` | **MIT** (LICENSE.md). Maintainer says the project is finished | Yes, with notice | https://github.com/mewo2/terrain/blob/master/LICENSE.md |
| mewo2/naming-language | `4abc8d3` | **MIT** (LICENSE.md) | Yes, with notice | https://github.com/mewo2/naming-language |
| redblobgames/mapgen4 | `c1d8cb0` | **Apache-2.0** (LICENSE; README "License") | Yes, with the Apache notice and change notes | https://github.com/redblobgames/mapgen4 |
| Mindwerks/worldengine | `80dd33e` | **MIT** (LICENSE.txt) | Yes, with notice | https://github.com/Mindwerks/worldengine |
| ftomassetti/civs (civilisation sim, Clojure) | `dc76bff` | **Apache-2.0** (LICENSE) | Yes, with notice. In practice an idea source, since we would rewrite it in Rust | https://github.com/ftomassetti/civs |
| ftomassetti/langgen (names, Python) | `794c0a1` | **Apache-2.0** (LICENSE) | Code yes. Its sample lists: **no provenance**, do not ship (see §7.4) | https://github.com/ftomassetti/langgen |
| watabou/TownGeneratorOS (Medieval Fantasy City Generator source) | `7fbc87a` | **GPL-3.0** (LICENSE) | **No.** Ideas only | https://github.com/watabou/TownGeneratorOS |
| raguilar011095/planet_heightmap_generation (climate, JS) | `5e434b3` | **GPL-3.0** (LICENSE) | **No.** Ideas only | https://github.com/raguilar011095/planet_heightmap_generation |
| Watabou's other repos | n/a | Not verified. The profile lists 8 public repos (6 shown on the page: pixel-dungeon, TownGeneratorOS, PD-classes, RuneGeneratorOS, CompassOS, pixel-dungeon-gdx). RuneGeneratorOS and CompassOS were not checked | Check each before use | https://github.com/watabou |
| Dwarf Fortress | n/a | **Not open source.** Freely redistributable binaries, per two secondary sources (Debian non-free archive; Fedora packaging). Bay 12's own page was not fetched | Ideas only | https://metadata.ftp-master.debian.org/changelogs/non-free/d/dwarf-fortress/stable_copyright ; https://forums.fedora-fr.org/viewtopic.php?id=70502 |
| Wonderdraft, Inkarnate | n/a | Closed commercial tools. Not checked online | Ideas only | n/a |

**Notes.**

- FMG bundles third-party libraries listed in `docs/wiki/Dependencies.md` (d3, Delaunator, polylabel, lineclip, Three.js, Alea, JSZip, Dropbox SDK, RgbQuant and others). I did not check their licences. If we copy a file that imports one of them, check that library first. https://github.com/Azgaar/Fantasy-Map-Generator/blob/master/docs/wiki/Dependencies.md
- FMG's own policy page says the source is MIT, maps are free to use, and "some bundled assets (fonts, icon sets, coat-of-arms artwork) carry their own licences". https://github.com/Azgaar/Fantasy-Map-Generator/blob/master/docs/wiki/Policy.md

---

## 3. Our inputs and how the simulations should run

**Inputs we already plan to commit** (ARCHITECTURE §6 to §8): Köppen class and monthly temperature and precipitation per cell, elevation, the river network with flux, lakes, the coastline and coast distance, plate boundaries, and erosion outputs.

**Scale.** FMG runs its generators on a Voronoi graph. Its default density is 10,000 cells, with a range of 1,000 to 100,000 (`src/data/graph-density.ts`). Our committed macro is about 25 million samples per field at L3, which is far too big for graph searches.

**Recommended society grid (v1):**

- **Use the L0 cube-sphere cells.** That is about 393,000 cells (6 faces × 256²) at about 39 km. That spacing suits towns; villages would need L1 (about 20 km) later.
- The grid's 4 or 8 neighbours produce blocky borders. Fix this with deterministic cost noise and by smoothing the border contours in the vector layer.
- **Alternative:** a spherical Fibonacci site set with a stored spherical Delaunay neighbour graph. It avoids cross-face seams but needs a deterministic sphere Delaunay. Decide this at the architecture review.

**Determinism rules** (ARCHITECTURE §5):

- **Integer costs** (fixed-point) in every Dijkstra or flood fill, with ties broken by cell id. Float ordering differs across platforms; integer costs don't.
- **Randomness** from an integer hash of (world seed, stream id, cell id). No `Math.random`.
- **Stored, not re-derived.** Everything the layer commits (cultures, states, burgs, routes, religions and so on) goes into the world file. Regeneration is an explicit command.
- **Locks.** FMG keeps "locked" cultures, states and burgs when the rest regenerates (`cultures-generator.ts`, `states-generator.ts`). That is a good pattern to copy as an idea.

**Cost estimate.** A flood fill over about 400,000 cells should take well under a second natively. This is an estimate, not a measurement. On wasm, run it as stepped `Job`s in chunks under 100 ms (ARCHITECTURE §5b). Single-threaded wasm is fine at this size; `rayon` is only a native speed-up.

**Köppen to habitability (first guess, to tune).** This table is our proposal, not FMG's. FMG uses its own biome table (`biomes-generator.ts`, line 51: values from 0 for water and ice-type biomes up to 100 for a temperate forest biome).

| Köppen group | Example classes | Habitability (0 to 100) |
|---|---|---|
| Tropical wet | Af, Am | 60 |
| Tropical, dry winter | Aw | 50 |
| Hot desert | BWh | 0 |
| Cold desert | BWk | 5 |
| Steppe | BSh, BSk | 25 |
| Mediterranean | Csa, Csb | 85 |
| Humid subtropical, oceanic | Cfa, Cfb | 90 |
| Monsoon subtropical | Cwa | 80 |
| Humid continental, warm summer | Dfa, Dfb | 65 |
| Subarctic | Dfc, Dfd | 20 |
| Tundra | ET | 5 |
| Ice cap | EF | 0 |

---

## 4. Feature by feature

Each entry gives FMG's method, the inputs it uses, how we would map it, and whether the idea is ours to reimplement or a code port.

### 4.1 Population and habitability

- **FMG does:** `src/model/society/population-generator.ts`, `rankCells()`. Each land cell gets a suitability score:
  - the biome's habitability;
  - up to +250 for big rivers and confluences (normalised flux);
  - minus (height − 50) / 5 for elevation;
  - coast bonuses: estuary +15, ocean coast +5, safe harbour +20, freshwater lake +30, salt lake +10, dry −5, lava −30.

  Rural population is suitability × cell area ÷ mean cell area. Burg population is about suitability ÷ 8, multiplied by 1.3 for capitals and ports, with random noise.
- **Inputs:** Köppen or biome, flux, elevation, coast flag, lake type, harbour flag.
- **Our mapping:** the Köppen table in §3 for the base score; our river flux (log-scaled); elevation penalty in metres; a slope penalty (FMG has none); coast and lake bonuses from our coastline and lake outputs. Store the score per cell and the rural population per cell.
- **Idea or port:** idea. Reimplement from this description. The numbers are tuning parameters, not expression.

### 4.2 Names

- **FMG does:** `src/model/society/names-generator.ts`. `calculateChain()` splits each sample name into pseudo-syllables and records which syllable can follow which letter. `getBase()` walks that chain to build a word within a min and max length, then cleans up doubled letters and capitals. Each culture points to a "name base": a comma-separated list of sample names (`src/data/name-bases.ts`, 43 default bases, including fantasy ones such as Elven, Dwarven, Orc and Goblin).
- **mewo2/naming-language (MIT)** takes a different approach. `language.js` builds a phoneme grammar: syllable templates with optional parts, phoneme classes, a random skew exponent, rejection regexes that ban certain syllables, and a per-letter spelling map. https://github.com/mewo2/naming-language ; method notes: http://mewo2.com/notes/naming-language/
- **langgen (Apache-2.0)** builds a generator from sample lists, using the same family of sequence model.
- **Our mapping:** we want a grammar generator (mewo2's shape, written by us in Rust) with per-culture phoneme inventories, plus a Markov fallback. Each culture gets its own inventory. Use the same generator for burgs, states, rivers, mountains and religions, with the culture as the input.
- **Idea or port:** idea for the Markov chain (a standard technique). For the phoneme grammar, reimplement from the idea, or port mewo2's MIT code with a notice. **Do not** copy FMG's or langgen's name lists (see §7.4).

### 4.3 Burgs: settlements, capitals and ports

- **FMG does:** `src/model/political/burgs-generator.ts`.
  - **Capitals:** rank populated land cells by suitability with random jitter. Take the best ones that are at least a minimum spacing apart, where spacing = (width + height) / 2 ÷ capital count. If there are not enough, shrink the spacing by ÷ 1.2 and retry.
  - **Towns:** the same method, with smaller randomised spacing that halves on each pass.
  - **Ports:** coastal candidates are chosen so ports spread out. The burg moves to the shore beside its "haven" (the water cell it opens onto).
  - **Burg type** from geography: Naval (port), Lake, Highland (height above 60), River (flux at least 100), Nomadic or Hunting (low population).
- **Inputs:** suitability, coast and haven, rivers with flux, elevation, Köppen.
- **Our mapping:** the same spacing rules on the L0 grid, measured by great-circle distance. Ports are coastal cells that touch the sea; the haven is the sea neighbour. River ports use the river's mouth or a navigable reach.
- **Idea or port:** idea.

### 4.4 Cultures

- **FMG does:** `src/model/society/cultures-generator.ts`, `expand()`.
  - Culture centres are picked from a named set (`src/data/culture-sets.ts`, for example "All-world" with up to 32 cultures).
  - Each culture gets a type from its geography. The rules are in FMG's own wiki (`docs/wiki/Culture-types.md`): Nomadic (hot or cold desert, or grassland, below height 70), Highland (height above 50), Lake (around a lake over 5 cells), Naval (coastal), River (flux at least 100), Hunting (savanna, rainforest, taiga, tundra or wetland, more than two cells from the coast), and Generic otherwise.
  - Expansion is a multi-source Dijkstra from the centres. Each step costs: a biome-change penalty of 20, a biome cost by culture type, a height cost (crossing sea, hills or mountains), a river cost, and a coastline cost, all divided by the culture's expansionism. A cap stops runaway growth. Only populated cells take the culture, but the flood keeps running through empty land.
- **Inputs:** biome or Köppen, elevation, rivers, coast, population.
- **Our mapping:** the same flood fill. Culture types come from our Köppen groups, elevation and flux. Costs come from Köppen change, slope and river crossings.
- **Idea or port:** idea. The cost table is a tuning choice.

### 4.5 States and borders

- **FMG does:** `src/model/political/states-generator.ts`.
  - `createStates()`: each capital becomes a state. Expansionism is random × size variety + 1.
  - `expandStates()`: the same flood fill, seeded at capitals, with a growth cap (cells ÷ 2 × growth rate). A step inside the state's own culture costs −9 against 100 for another culture. Sea crossing costs 1,000. Mountains cost 2,200, hills 300. Empty land costs 5,000, populated land costs little.
  - Then `normalize()`, `findNeighbors()`, `assignColors()`, `generateCampaigns()` and `generateDiplomacy()`.
- **Inputs:** capitals, cultures, population, elevation, rivers, biome.
- **Our mapping:** the same flood fill as cultures, seeded at capitals. Borders are the cell boundaries, smoothed into contours in the vector layer (ARCHITECTURE §8). Borders must stay put when the user zooms, so they belong in the committed vector features.
- **Idea or port:** idea.

### 4.6 Provinces

- **FMG does:** `src/model/political/provinces-generator.ts`. Each state gets a number of provinces equal to its burg count × ratio (minimum 2). Province seeds are its burgs: capital first, then by population with random weight. Each province grows inside its state with the same cost-based flood.
- **Inputs:** states, burgs, population.
- **Our mapping:** after states exist, run the same flood inside each state. Province labels come from the vector layer.
- **Idea or port:** idea. It's cheap once states exist.

### 4.7 Roads and trade routes

- **FMG does:** `src/model/economy/routes-generator.ts`.
  - **Main roads:** capitals on each landmass are joined by an Urquhart graph (a Delaunay triangulation with the longest edge of each triangle removed). Then least-cost paths are found.
  - **Land cost** per step: squared distance × habitability factor (1 to 1.1) × height factor (about 2 to 4 for heights 0 to 100; the code comment says 1 to 3, but the formula gives 2 to 4) × 0.5 on an existing road × 3 when there is no burg in the cell. Uninhabitable cells are impassable.
  - **Trails** join the other burgs. **Sea routes** and **river routes** use a separate water cost; rivers are navigable above a flux threshold.
  - Routes are drawn along the cell path, with meanders on river routes.
- **Inputs:** burgs, habitability, height, rivers (navigable reaches), coast, the sea.
- **Our mapping:** the same Urquhart and least-cost approach on the L0 grid, with great-circle distance. Sphere Delaunay for a few hundred capitals can use a simple 3D convex hull. River routes use our flux network.
- **Idea or port:** idea.

### 4.8 Religions

- **FMG does:** `src/model/society/religions-generator.ts`.
  - **Folk religions** are one per culture (or per ancestry). They simply cover their culture's cells.
  - **Organised religions** are seeded at centres and spread with a passage cost. On land, roads cost 1, other routes cost a third of the biome cost, and cells without a route cost the biome cost. Water costs 50 with a route and 500 without. Expansion is capped.
  - **Heresies** split off from organised religions with expansionism at least 3, on their boundary cells.
  - **Names** come from forms such as "Culture + ism", "Faith of + deity" or "Place + ian + type", with the culture's name generator.
- **Inputs:** cultures, states, roads, biome.
- **Our mapping:** the same flood as cultures, with road costs from our route layer. Folk religion follows the culture; organised religion follows roads and states.
- **Idea or port:** idea.

### 4.9 Markers and points of interest

- **FMG does:** `src/model/story/markers-generator.ts`. A config table lists about 25 marker types (for example volcanoes, waterfalls, mines, bridges, inns, lighthouses, caves, ruins, battlefields and portals). Each type has a candidate function (for example `listVolcanoes`, `listWaterfalls`, `listInns`) and parameters (minimum, per-candidate rate, multiplier).
- **Inputs:** rivers, coast, roads, volcanic and tectonic features, mountains.
- **Our mapping: physics, not random.** This is where our terrain model helps most.
  - **Waterfalls:** knickpoints in the river long profile, which the stream-power solver already produces.
  - **Volcanoes:** convergent plate boundaries with an arc, from tectonics.
  - **Hot springs:** near volcanoes and rifts.
  - **Lighthouses:** convex headlands on the coastline contour.
  - **Bridges:** where a road crosses a river above a flux threshold.
  - **Inns:** along roads at about a day's travel (25 to 35 km).
  - **Mines:** mountain belts. There is no ore geology yet, so these are flavour only.
- **Idea or port:** idea.

### 4.10 Military

- **FMG does:** `src/model/political/military-generator.ts`. Each state's army is rural population ÷ 100 × a unit share × a terrain modifier, plus urban population × urbanisation × a share. Everything is multiplied by the state's "alert" (its expansion share × diplomacy rate × neighbour pressure, clamped 0.1 to 5). Units are grouped into regiments of about three platoons, each with a name and an emblem.
- **Inputs:** states, population, diplomacy, terrain.
- **Our mapping:** a per-state force count from our population layer, with garrison markers at border burgs. Simplest form first.
- **Idea or port:** idea. Later version.

### 4.11 Diplomacy

- **FMG does:** `states-generator.ts`, `generateDiplomacy()`. Each pair of states gets a relation from a weighted table: Ally, Friendly, Neutral, Suspicion, Rival, Enemy, Vassal, Suzerain. The weights depend on being a neighbour, a neighbour of a neighbour, or far away. Vassals copy their suzerain's relations. Campaigns (wars, conquests, invasions, rebellions and so on) come from neighbour pairs, with random names and start years.
- **Inputs:** the neighbour graph, state expansionism, culture, religion.
- **Our mapping:** the same weights on our state adjacency. Campaigns feed the history (§4.12).
- **Idea or port:** idea. Later version.

### 4.12 History and notes

- **FMG does:** free-text notes on any entity (`src/model/story/notes.ts`). Diplomacy changes and campaigns feed a chronicle (`states-generator.ts`, `getChronicle()`). Journeys are generated along routes (`src/model/story/journeys/`). Event zones such as invasion, rebels, plague, eruption, flood and tsunami cover cells (`zones-generator.ts`).
- **Our mapping, minimal first:** a founding year per burg and state, plus a short event list. Full chronicle later. Event zones can come from our tectonics (eruptions, faults) and rivers (floods).
- **Idea or port:** idea. Dwarf Fortress's legends mode is the best-known example of generated history (sites, figures, wars and artefacts). Its details were not checked online this session (the DF wiki was blocked), so treat that as idea-only general knowledge.

### 4.13 Emblems and heraldry

- **FMG does:** `src/model/political/emblems-generator.ts`. Each state, culture, burg and religion gets a generated coat of arms. Tinctures follow the "rule of tincture" (metal on colour, and so on). Shields, ordinaries and charges come from fixed lists. The charge art is in `src/data/emblems/`.
- **Armoria** is the heraldry generator and editor, used by FMG. Its code is MIT per its README and package.json (see §2). Its charge art is not (see §7.4).
- **Our mapping:** write our own blazon generator (tincture rules, a few ordinaries, simple charges). Draw our own CC0 charges. Heraldic ideas and the rule of tincture are common practice and free to use. The written blazon is not copyrightable, but a drawn render may be (see §7.1).
- **Idea or port:** idea for the generator. Art: draw our own. Later version.

### 4.14 Labels

- **FMG does:** label groups per type (`src/model/annotations/labels-generator.ts`). Each group has a style. `src/controllers/annotations/labels/label-spread.ts` moves overlapping labels, and path labels run along rivers and routes.
- **Our mapping:** labels live in the vector layer (ARCHITECTURE §8, M7). Use the committed features: river polylines, ranges from elevation maxima, burg sizes from population.
- **Idea or port:** idea. The label-spread code is MIT, so porting is also an option, with a notice.

---

## 5. Our build order and how features depend on each other

```
population (4.1) ─┬─> burgs (4.3) ─┬─> roads (4.7) ─> religions (4.8), bridges and inns (4.9)
                  │                └─> states (4.5) ─┬─> provinces (4.6), diplomacy (4.11), military (4.10)
names (4.2) ──────┴─> cultures (4.4) ─────────────────┘
                                                      history (4.12), emblems (4.13)
labels (4.14) and markers (4.9) come with the vector layer
```

---

## 6. Ranking by aliveness per effort

Scores are my judgement for a first version. **A** = aliveness (1 to 5, higher is better). **E** = effort (1 to 5, higher is harder). The ratio is A ÷ E.

| Rank | Feature | A | E | A ÷ E | Depends on | Notes |
|---|---|---|---|---|---|---|
| 1 | Population and habitability | 4 | 1 | 4.0 | Köppen, rivers, coast | Drives all others; visible as density shading |
| 2 | Names | 4 | 1.5 | 2.7 | none (per culture later) | Every label gets a name |
| 3 | Burgs, capitals, ports | 5 | 2 | 2.5 | population | The first thing people read on a map |
| 4 | Labels | 4 | 2 | 2.0 | features | Needed anyway for M7 |
| 5 | States and borders | 5 | 3 | 1.7 | burgs, cultures | Biggest visual change; reuses the flood fill |
| 6 | Markers (physics-based) | 3 | 2 | 1.5 | rivers, tectonics, roads | Cheap if derived from terrain |
| 7 | Minimal history | 3 | 2 | 1.5 | burgs, states | Founding years and a short event list |
| 8 | Cultures | 4 | 3 | 1.3 | population, names | Prerequisite for states and religions |
| 9 | Roads and trade | 4 | 3 | 1.3 | burgs, habitability | Shows connection; needs a path search |
| 10 | Provinces | 3 | 2.5 | 1.2 | states | Cheap once states exist |
| 11 | Religions | 3 | 3 | 1.0 | cultures, roads | Good depth; lower priority |
| 12 | Heraldry | 3 | 4.5 | 0.7 | states | Needs our own art |
| 13 | Military | 2 | 3 | 0.7 | population, states | Mostly numbers until the map shows it |
| 14 | Diplomacy | 2 | 3.5 | 0.6 | states | Mostly text until history uses it |

**Recommendation.** Build 1 to 5 first. Cultures (rank 8) must come before states, and names (rank 2) before cultures are shown, so the order is population, names, burgs, cultures, states. Labels arrive with M7.

---

## 7. Legal and safety

### 7.1 Safe to reimplement from the idea

Ideas, methods, algorithms and mechanics are not protected by copyright. In US law, 17 U.S.C. §102(b) says so directly: no protection for "any idea, procedure, process, system, method of operation, concept, principle, or discovery". https://www.govinfo.gov/link/uscode/17/102 (found by search, text not fetched in full). Other jurisdictions have similar rules; check before publishing in one of them.

Safe to reimplement from the description in this report:

- suitability from habitability, rivers and coast;
- spacing rules for burgs and capitals;
- cost-weighted flood fills for cultures, states, provinces and religions;
- the Urquhart graph and least-cost roads;
- Markov chains and phoneme grammars for names (the technique, not any specific list);
- folk and organised religion rules, heresies;
- diplomacy weights and campaign generation;
- the rule of tincture and the blazon idea;
- label overlap avoidance.

**Rule of thumb.** Our own cost numbers and weights are fine. Copying a file's text, its structure, its comments or its data is not.

**Renders.** A written blazon is not copyrightable, but a drawn coat of arms may be. A Wikimedia Commons-style licence notice (search result, secondary) says the same: the written description is free of copyright, but a drawn render may be protected. We draw our own.

### 7.2 Code we could port, with obligations

| Licence | Projects | Obligations | Practical note |
|---|---|---|---|
| MIT | FMG, mewo2/terrain, mewo2/naming-language, worldengine | Keep the copyright and licence notice in copies or substantial portions | Our rules allow MIT with attribution |
| Apache-2.0 | mapgen4, civs, langgen | Include the licence text and NOTICE, mark changed files, do not use the names as trademarks | Gives a patent grant too |

- Porting between languages (JS, Clojure or Python to Rust) is a rewrite. Reimplementing from the description is usually simpler and gives cleaner provenance.
- If we port any file, its header should give the source URL, commit hash, licence and what we changed.
- Armoria code: **hold** until it has a LICENSE file. A package.json field and a README line are a strong signal, but not a licence text. Ask the author to add the LICENSE file.

### 7.3 Do not copy or mix

- **GPL-3.0:** `watabou/TownGeneratorOS` and `raguilar011095/planet_heightmap_generation`. Copying code or a near-verbatim translation would pull GPL obligations into the project. Ideas only.
- **Watabou's other repos:** not checked. Treat as GPL until a LICENSE file says otherwise.
- **Armoria:** no licence file (see §7.2).
- **Any charge art or heraldry artwork from FMG or Armoria** (see §7.4).
- **FMG's public images, heightmap PNGs and icon files:** licences not checked. Do not copy.

### 7.4 Art, data and name lists: red flags

| Item | Problem | Action |
|---|---|---|
| FMG's name bases (`src/data/name-bases.ts`) | Labelled "real-world bases by Azgaar". No sources listed. Many entries are real place names | Do not copy. Write invented syllable inventories |
| langgen sample lists (`lang_samples/`) | Welsh, Norse, Greek, Roman and other lists. The Apache licence covers the code, not these lists. Provenance not given | Do not ship |
| FMG and Armoria complex charge art | Armoria's README: complex renders from WappenWiki, **CC BY-NC 3.0** (non-commercial). FMG's PRD says its charge art is under mixed licences, including CC BY-NC 3.0 for WappenWiki renders. Armoria README: https://github.com/Azgaar/Armoria ; FMG PRD: https://github.com/Azgaar/Fantasy-Map-Generator/blob/master/docs/prd/emblems-icon-library.md | **Never copy.** Our rules block CC BY-NC-SA; CC BY-NC is also incompatible with an open, commercially usable project. Confirm with Matt (§7.7) |
| FMG's simple charges | Self-made, **CC0** (Armoria README; FMG PRD says the same for simple charges) | Usable, but we should still draw our own |
| Natural Earth | Public domain (secondary sources; the official page was not reachable). Some contributed layers carry separate notices. https://en.wikipedia.org/wiki/Natural_Earth | OK for real-world reference. Check each layer's notice. We don't need it for a fantasy planet |
| GeoNames | **CC BY 4.0** per secondary sources (some older copies say 3.0). Credit required. https://geonames.org/about.html (search result, not fetched) | If used, credit GeoNames. Prefer not to ship real place names in a fantasy tool |
| Wikipedia text and lists | **CC BY-SA 4.0** and GFDL (Wikimedia terms). Share-alike and attribution. https://meta.wikimedia.org/wiki/Terms_of_use/Creative_Commons_4.0 (search result) | Do not copy text or lists. Share-alike would spread into our data |
| HYDE 3.x (population history) | Licence not verified. Utrecht record: https://public.yoda.uu.nl/geo/UU01/8K9D7F.html (search result) | Check the record before use |
| Beck et al. 2018 Köppen maps | Licence not shown in the Figshare search result. https://figshare.com/articles/dataset/Present_and_future_K_ppen-Geiger_climate_classification_maps_at_1-km_resolution/6396959 | Check the Figshare licence field. Fine for validation if the licence allows |
| OpenMoji (emoji art) | Believed to be CC BY-SA 4.0. **Not verified.** https://openmoji.org | Check before shipping any of it |
| game-icons.net (FMG's free-icon source, `src/controllers/icon-picker/index.ts`) | Believed CC BY 3.0. **Not verified.** https://game-icons.net | Check before use |
| Wikimedia Commons images | Each file has its own licence, many CC BY-SA | Check each file |
| Unicode emoji (FMG marker glyphs, e.g. the volcano) | Glyph rendering varies. egui's default fonts lack some emoji (HANDOFF.md) | Draw our own icons |

### 7.5 Method notes for the legal read

- I verified LICENSE files from clones at the commits in §2. The licence claims for those projects are direct.
- Natural Earth, GeoNames, Wikipedia and the statute are from search results or secondary sources. Re-check them against the primary sources before we rely on them.

### 7.6 Trademarks and names to keep out of product and UI

- Product or marketing names: **Wonderdraft**, **Inkarnate**, **Dwarf Fortress**, **Watabou**, **Azgaar**, **Armoria**. Mentioning them in research notes is fine.
- **Tolkien** and **Middle-earth** names, and **Dungeons & Dragons** terms (Wizards of the Coast). Check before using any of them in a product name.
- **Google Earth**. Our architecture doc uses the phrase "Google Earth, but fantasy" internally. Keep it out of public copy.

I did not check trademark registrations. This is a caution list, not a clearance.

### 7.7 Decisions for Matt

1. **Our project licence** (still TBD). MIT or Apache-2.0 accepts every MIT and Apache input in §2.
2. **CC BY-NC.** I recommend treating it as blocked, like CC BY-NC-SA. This affects FMG's charge art and Armoria's renders.
3. **Armoria.** Ask Azgaar for a LICENSE file before any code is ported. Until then, reimplement from the idea.
4. **Name lists.** Confirm that we write our own invented inventories, with no real-world lists.

---

## 8. Open items and corrections needed in other files

I did not change these files, per the task instructions.

- **`docs/research/2026-11-02-climate-tier-b.md`**, lines 3 and 11. It says `raguilar011095/planet_heightmap_generation` has an unverified licence and recommends reimplementing from docs. That is now **verified GPL-3.0** (commit `5e434b3`). Keep the no-copy rule and update the note.
- **`docs/research/2026-11-02-terrain-erosion.md`**, line 12. It names `bevy_erosion_filter` as MIT. The s01 session log corrected this to **MPL-2.0**. The research doc still says MIT.
- **Unverified items to check before relying on them:** Dwarf Fortress's official terms, Wonderdraft and Inkarnate, HYDE, the Beck Köppen dataset licence, OpenMoji, game-icons.net, JSZip (FMG dependency), Watabou's RuneGeneratorOS and CompassOS, the official Natural Earth and GeoNames terms.
- **Sandbox access.** `www.naturalearthdata.com`, `www.geonames.org` and `www.dwarffortresswiki.org` need to be added to the session's allowed network domains to verify the dataset claims directly.
- **Session logs and handoff.** The project rules (CLAUDE.md) ask for a session-log entry and a HANDOFF update, followed by a commit and push. This task said not to modify other files and not to commit or push, so I did neither. The launching agent should decide whether to record this report in the session log.
- **Pre-existing change.** `crates/rg-pipeline/examples/probe_realism.rs` has an uncommitted modification that predates this task. I left it alone.

---

## 9. Sources

Repos (cloned; LICENSE and README read at the commits in §2):

- https://github.com/Azgaar/Fantasy-Map-Generator (LICENSE; `src/model/**`; `docs/wiki/Culture-types.md`, `docs/wiki/Policy.md`, `docs/wiki/Dependencies.md`, `docs/prd/emblems-icon-library.md`)
- https://github.com/Azgaar/Armoria (README "License"; `package.json`; full-history LICENSE check)
- https://github.com/mewo2/terrain (LICENSE.md; `terrain.js`)
- https://github.com/mewo2/naming-language (LICENSE.md; `language.js`)
- https://github.com/redblobgames/mapgen4 (LICENSE; `README.org`)
- https://github.com/Mindwerks/worldengine (LICENSE.txt; README)
- https://github.com/ftomassetti/civs (LICENSE; README; `src/civs/`)
- https://github.com/ftomassetti/langgen (LICENSE; README; `lang_samples/`)
- https://github.com/watabou/TownGeneratorOS (LICENSE; README)
- https://github.com/raguilar011095/planet_heightmap_generation (LICENSE)
- https://github.com/watabou (repo list)

Search results (secondary; not fetched in full): Natural Earth, GeoNames, Wikimedia terms, 17 U.S.C. §102, HYDE and the Beck Köppen dataset, Dwarf Fortress licence status (Debian and Fedora packaging pages).

Blocked from fetch this session: www.naturalearthdata.com, www.geonames.org, www.dwarffortresswiki.org.
