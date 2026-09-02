# RealmGenesis 3D — Documentation

This folder is the **settled** knowledge about RealmGenesis: what the system
is, how it works, and *why* it was built that way. It is rebuilt from the code,
not from older docs — every factual claim should cite a `file:line` you can
check.

## The three-document rule

Knowledge about this project lives in exactly one of three places, by its
time-tense. When you're unsure where something belongs, this rule decides:

| Where | Tense | Holds |
|-------|-------|-------|
| **`docs/`** (here) | settled | What is, and why it was decided that way. Architecture, data model, invariants, test philosophy, engineering decisions. |
| **`HANDOFF.md`** (root) | live | Recent-session state, unconfirmed findings, traps, "not pushed yet". Perishable. |
| **`ROADMAP.md`** (root) | future | Features not yet built, milestone ordering. |

A fact graduates from HANDOFF to `docs/` once it stops changing — a decision
that's been made and won't be re-litigated, an invariant that holds. Until
then it stays in HANDOFF at the confidence level it actually has.

`CLAUDE.md` and `AGENTS.md` stay at the repo root: they are harness/agent
convention entrypoints, not reference docs.

## Contents

| Doc | Covers | Status |
|-----|--------|--------|
| [architecture.md](architecture.md) | System overview, state ownership (`useWorldEngine`), shell routing, the Web-Worker generation boundary, module map | ✅ current |
| [generation-pipeline.md](generation-pipeline.md) | The worker → `generateWorld` stage pipeline as it runs today | ✅ current |
| [tectonics-v3.md](tectonics-v3.md) | The V3 terrain model: crust fields, Euler-pole kinematics, Dijkstra plate growth, GDH1 bathymetry, microplates | ✅ current |
| [data-model.md](data-model.md) | `Cell`, `WorldData`, `WorldParams`, `BiomeType` (17), `ViewMode` (12), and the six worldbuilding entity types | ✅ current |
| [params-reference.md](params-reference.md) | Every `WorldParams` key with range/default/consumer, and the param-liveness contract | ✅ current |
| [rendering.md](rendering.md) | 3D globe (`WorldViewer`), 2D map (`Map2D`), Dymaxion, coloring, labels, hillshading | ✅ current |
| [map-styles.md](map-styles.md) | The A3 style layer: substrate adapter, parchment passes, glyph placement, the fill-policy rule, the mirror trap, how to add a style | ✅ current |
| [civilization.md](civilization.md) | Cultures → religions → factions → provinces → towns → routes → markers → AI lore | ✅ current |
| [export.md](export.md) | PNG raster, SVG, GeoJSON, GLB, and save/load persistence | ✅ current |
| [dymaxion.md](dymaxion.md) | Dymaxion projection, the classic vs Blender net layouts, and the Blender-UV interop (extraction script + why exports drop on cleanly) | ✅ current |
| [performance-findings.md](performance-findings.md) | The F4 **measured baseline** and ranked implementation queue — globe tenant costs, the 2.2 s Map2D gesture, and three optimisations whose correctness numbers already exist | ✅ current |
| [performance-workflow.md](performance-workflow.md) | The F4 method: the four measurement harnesses, the before/after/**correctness** rule, and how to delegate perf work without paying orientation cost N times | ✅ current |
| [invariants.md](invariants.md) | Non-obvious facts that break things if violated — each re-verified against current code | ✅ current |
| [testing.md](testing.md) | Test philosophy, the three determinism instruments, why there is no golden fixture, the param-liveness contract | ✅ current |
| [ENGINEERING-NOTES.md](ENGINEERING-NOTES.md) | Shelved-not-abandoned levers, refuted hypotheses, decisions with their rationale | ✅ current |
| [archive/](archive/) | Point-in-time documents kept for the record, not current truth (dated) | — |

Status key: ✅ current · 🟡 pending/stub · 📦 archived.

## Archive

- [`archive/ARCHITECTURE-legacy.md`](archive/ARCHITECTURE-legacy.md) — the
  original monolithic architecture doc (last accurate ~Session 6). **Superseded
  by this folder and known to drift** (it predates the Web Worker, the V3
  terrain model, the F1 shell, and `tectonicStrength`). Kept only so old links
  resolve and to mine for topics not yet re-covered.
- [`archive/audit-2026-07-17.md`](archive/audit-2026-07-17.md) — a technical
  audit from before CI, tests, the worker, and V3 existed. Most of its findings
  are now resolved; read the header before trusting any item.

`research/` and `superpowers/` (specs + plans) are unchanged.
