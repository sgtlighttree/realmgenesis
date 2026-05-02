<!--
  HANDOFF DOCUMENT
  ─────────────────────────────────────────────────────────────────────────────
  When the user asks for a new handoff for a new session, the contents of this
  file can be fully replaced with the updated state. The instruction to clear it
  should come from the user explicitly ("new handoff" / "update the handoff").
  ─────────────────────────────────────────────────────────────────────────────
-->

# RealmGenesis — Agent Handoff

This document gives an incoming agent the minimum context needed to continue
work on RealmGenesis without reading the full git history. Read
[`ARCHITECTURE.md`](./ARCHITECTURE.md) for deep technical detail on any
specific system; this file covers *what was done recently*, *why*, and *what
is immediately actionable*.

---

## Project snapshot

**RealmGenesis 3D** — browser-based procedural fantasy world generator. Runs
entirely client-side (React 19, Three.js, d3-geo-voronoi). No backend.

- Dev: `npm run dev` (port 3000, or 3001 if 3000 is taken)
- Build check: `npm run build` — must pass with zero TS errors before any commit
- Lint: `npm run lint` — zero errors enforced
- Deployment: Netlify (static SPA, `_redirects` in `public/`)
- AI lore: Google Gemini BYOK (ephemeral key, never persisted)

---

## Work completed in this session

### Bug fixes

| Commit | What was fixed |
|--------|---------------|
| `8afb3be` | **Globe blanks on Borders view** — `<CountryLabels>` uses `<Text>` from drei which calls `suspend()` for async font loading. R3F's Canvas wraps all children in a top-level `<Suspense>`; when CountryLabels suspended, the entire canvas went blank. Fixed by wrapping `<CountryLabels>` in its own `<React.Suspense fallback={null}>` inside `WorldMesh` so only the labels hide during font loading. |
| `8afb3be` | **Seed shuffle button had no visible effect** — `handleRandomizeSeed` called `handleChange` twice in sequence; both calls closed over the same stale `params`, so the second call overwrote the first's seed change. Fixed by merging both updates into a single `setParams(prev => ...)` functional updater. |
| `a26ad2f` | **Horizontal scan-lines in 2D map and all exports** — Canvas 2D anti-aliases path fill edges, leaving the dark background visible as ~1px gaps at Voronoi cell boundaries. In equirectangular projection these gaps align horizontally and appear as continuous lines. Fixed by stroking each cell polygon with its own fill color (`lineWidth = 1`) immediately after filling. Applied to all four cell-render loops (Map2D Mercator, Map2D Dymaxion source canvas, `export.ts renderEquirectangular`, `export.ts exportMap`). |

### Features / improvements (minor improvements batch, committed earlier)

- World stats panel (Sys tab, collapsible)
- Keyboard shortcuts: G = generate, R = randomize seed, Escape = close inspector, 1–5 = switch tab
- Generation progress bar (thin blue bar above Generate button)
- Copy-seed button with 1.5 s checkmark feedback
- Heightmap PNG export button in Exp tab
- Dymaxion orientation presets (N.Pole / Pacific / Atlantic / Asia)
- Slider tooltips on ambiguous parameters

### Terrain generation improvements (latest commit `b6f2723`)

1. **Tectonic Strength slider corrected** — was 0–2.0 but worldGen.ts clamps `plateInfluence` to [0.1, 1.0] internally; slider now reflects the true range (0.1–1.0, step 0.05).

2. **Mountain Heights + Sea/Trench Depth sliders** — two new `WorldParams` fields (`mountainHeight`, `oceanDepth`, both default 1.0). After all height normalization (Stage 9b in the pipeline) a power-curve remap is applied independently to land cells and ocean cells. `mountainHeight > 1` pushes peaks taller; `oceanDepth > 1` pushes trenches deeper. Heights stay in [0,1]. Sliders live in the Geo tab immediately below Ridge Intensity (yellow-300 / amber-600 accent colors).

3. **Archipelago and Islands presets reworked** — the root cause of both looking like large continents was `plateInfluence` being at or above 1.0 (max tectonic structure) and `seaLevel` defaulting to 0.55. Both presets now set `plateInfluence` low (0.25 / 0.4), `seaLevel` high (0.65 / 0.60), include `plates` and `roughness`, and use the new `mountainHeight` / `oceanDepth` fields. Continents and Pangea also reset the new params to 1.0 when selected.

---

## Files most recently touched

| File | What changed |
|------|-------------|
| `types.ts` | Added `mountainHeight: number` and `oceanDepth: number` to `WorldParams` |
| `App.tsx` | Added `mountainHeight: 1.0` and `oceanDepth: 1.0` to `DEFAULT_PARAMS` |
| `utils/worldGen.ts` | Inserted Stage 9b height remap block after post-erosion normalization |
| `components/Controls.tsx` | Fixed Tectonic Strength slider range; added Mountain Heights + Sea/Trench Depth sliders; reworked Archipelago and Islands preset values; `handleRandomizeSeed` stale-closure fix |
| `components/WorldViewer.tsx` | Wrapped `<CountryLabels>` in `<React.Suspense fallback={null}>` |
| `components/Map2D.tsx` | Added same-color stroke after fill in all cell-render loops |
| `utils/export.ts` | Added same-color stroke after fill in both cell-render loops |
| `ARCHITECTURE.md` | Updated WorldParams table, pipeline stages, Controls tab description, added invariants 16–19 |

---

## Key invariants to know before touching anything

These have all bitten us or been explicitly noted — read the full list in
`ARCHITECTURE.md § Key Invariants & Gotchas`:

- **`plateInfluence`** is hard-clamped to [0.1, 1.0] in worldGen.ts. Slider max is now 1.0 to match.
- **`mountainHeight` / `oceanDepth` remap** must stay *after* all normalization and *before* climate/biome assignment (Stage 10). If you add a new normalization pass, move Stage 9b after it.
- **`CountryLabels` Suspense boundary** — do not move CountryLabels outside its `<React.Suspense>` wrapper or the 3D globe will blank on first switch to Borders view.
- **Multi-field `setParams`** — always use a single functional updater `setParams(prev => ({ ...prev, a, b }))` when updating more than one field. Two sequential non-functional calls lose the first update (stale closure).
- **Cell boundary seam prevention** — every cell polygon must be filled AND stroked with the same color. Do not add new cell-rendering loops without this pattern.
- **No test framework** — quality gate is `npm run build` + `npm run lint` zero errors + manual browser testing.
- **All state lives in App.tsx** — no Context/Redux/Zustand. Trace any value up to `App.tsx`.
- **Ephemeral Gemini key** — `apiKey` state is never persisted; do not add persistence.

---

## Potential next tasks (not committed to)

The user has not explicitly queued these, but they came up in conversation or
are natural follow-ons:

- **Preset fine-tuning after visual testing** — the Archipelago/Islands preset values are best-effort without live testing. The user may want to adjust `seaLevel`, `noiseScale`, or `plateInfluence` after seeing results.
- **Map painting — terrain and political** — user-requested major feature. Two sub-modes:
  - *Terrain paint*: click or drag over cells in the 3D globe or 2D Mercator view to raise/lower height or force-assign a biome. After a stroke, re-run `determineBiome` on affected cells and optionally re-propagate moisture/temperature for realism. The cell graph (`WorldData.cells`) is mutable in place; a shallow `setWorld({ ...world })` triggers a re-render.
  - *Political paint*: reassign `cell.regionId` / `cell.provinceId` by dragging a brush over cells, effectively redrawing faction borders. No terrain recalculation needed — just color/border re-render.
  - Both modes should be gated behind an "Edit Mode" toggle (separate from `inspectMode`) so normal orbit/pan interaction isn't disrupted. The 2D pick buffer (`pickCanvasRef` in Map2D.tsx) is already set up for per-cell hit detection and is the natural entry point for 2D painting. In 3D, `faceMap` in WorldViewer.tsx maps ray-cast triangle index → cell ID.

- **World editor — names, demographics, and faction data** — user wants granular per-cell and per-faction editing beyond what AI lore provides:
  - *Cell / city editing*: rename towns (`TownData.name`), adjust `population`, reassign `isCapital` / `isTown` flags.
  - *Faction editing*: rename factions (`FactionData.name`, `color`, `description`), adjust `totalPopulation`, merge or split factions.
  - *Province editing*: rename provinces, reassign towns between provinces.
  - All of this data already lives in `WorldData.civData` (see Political Hierarchy in ARCHITECTURE.md). An editor panel (likely extending the existing `Inspector` HUD or adding a dedicated Edit sidebar tab) that exposes these fields as text inputs / color pickers / number inputs would surface them. Changes are in-place mutations followed by `setWorld({ ...world })`. No generation re-run is required.

- **Additional `landStyle` masks** — currently only `'None'` and `'Pangea'`; a `'Ring'` or `'Twin Continents'` mask is a natural extension of the existing `maskType` mechanism.
- **Update ARCHITECTURE.md** — the user typically asks for doc updates after each batch of changes; this was done at the end of this session.

---

## Workflow notes

- The user commits manually and asks the agent to stage/commit/push when satisfied. **Do not push without explicit request.**
- Commit messages follow the 50/72 rule, single-line imperative.
- All commits include `Co-Authored-By: <model-name> <noreply@anthropic.com>` (substitute appropriate model identity).
- The user runs the dev server themselves via `! npm run dev` when testing. The agent should not start the dev server unless explicitly asked.
- `HANDOFF.md` (this file) is updated on request at the end of a session. When the user says "new handoff," replace this file's content entirely with the new session's state.
