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
work without reading the full git history. Read [`ARCHITECTURE.md`](./ARCHITECTURE.md)
for deep technical detail; this file covers *what was done recently*, *why*, and
*what is immediately actionable*.

---

## Project snapshot

**RealmGenesis 3D** — browser-based procedural fantasy world generator. Runs
entirely client-side (React 19, Three.js, d3-geo-voronoi). No backend.

- Dev: `npm run dev` (port 3000, or 3001 if taken)
- Build check: `npm run build` — zero TS errors required before any commit
- Lint: `npm run lint` — zero errors enforced
- Deployment: Netlify (static SPA)
- AI lore: Google Gemini BYOK (ephemeral key, never persisted)

---

## Work completed in this session

### Major feature: Edit Mode

A full interactive map editing system was added across multiple rounds of
implementation and polish. Two top-level capabilities:

**1. Map Painting** — click/drag in the 3D globe or 2D Mercator view to
modify the world:

| Mode | Behaviour |
|------|-----------|
| `terrain-raise` | Increase cell heights; Adaptive or Freeform style |
| `terrain-lower` | Decrease cell heights |
| `terrain-flatten` | Right-click to sample elevation, drag to flatten to that level |
| `terrain-smooth` | Average cell heights with neighbours |
| `biome` | Force-assign a biome type to cells |
| `political` | Reassign faction borders (`cell.regionId`) |

Controls: brush size (BFS ring radius 0–5), strength 0.1–1.0, adaptive biomes
toggle (gates whether terrain strokes re-run `determineBiome`). Space+drag
restores orbit (3D) or pan (2D) while painting.

Undo: Ctrl+Z, up to 20 strokes. Snapshot pushed at stroke **start** with a
Map reference; the Map grows in-place throughout the drag. This ensures undo
works even if the 'end' event is missed.

**2. World Editor** — `'world-edit'` mode; click a cell → Inspector HUD shows
editable inputs for faction name, color, description, province name, town name,
and town population. Changes are in-place mutations followed by
`setWorld({ ...world })`.

Faction editing also exposed in the **CIV tab** as a new "Factions" section
with inline color picker + name field per faction.

### Bug fixes this session

| Bug | Root cause | Fix |
|-----|-----------|-----|
| 3D globe blanks on Borders view | `<CountryLabels>` calls `suspend()` for font; R3F top-level `<Suspense>` blanked entire canvas | Wrapped `<CountryLabels>` in its own `<React.Suspense fallback={null}>` inside WorldMesh |
| Seed shuffle had no effect | Two sequential `setParams` calls with stale closure — second overwrote first | Merged into `setParams(prev => ...)` single functional updater |
| 2D horizontal scan-lines | Canvas 2D fill anti-aliasing left background visible at cell edges | Same-color `ctx.stroke()` after every `ctx.fill()` in all four cell-render loops |
| Faction colors mismatched map | Chips showed `f.color` (random palette); map rendered `PLATE_COLORS[regionId]` — different arrays | New `FACTION_COLORS` array used for both; `getCellColor` accepts optional `factionColors` map |
| First faction chip unresponsive | `paintFaction` defaulted to `0`; faction IDs are capital cell indices (large numbers, not 0) | EditToolbar initialises `paintFaction` to `factions[0].id` on civData load |
| Undo not working | Snapshot pushed at 'end' phase; 'end' unreliable if pointer released outside element | Push at 'start' with Map reference; Map grows during stroke; 'end' just nulls ref |
| Political brush painting distant cells (partial) | Pick canvas used fill-only → anti-aliased boundary pixels decoded to wrong cell IDs | Same-color stroke added to pick canvas loop; chord-distance cap in BFS |

### Other improvements

- **Terrain sliders** — Tectonic Strength range corrected (0.1–1.0); Mountain Heights + Sea/Trench Depth sliders added (power-curve height remap, Stage 9b)
- **Presets reworked** — Archipelago and Islands now set `plateInfluence` low, `seaLevel` high, and include `roughness`/`plates`
- **Faction colors** — `FACTION_COLORS` palette (18 perceptually distinct hues); shuffled with `civRng` for variety while remaining deterministic per seed
- **Brush visualization** — white ring on globe surface; CSS circle overlay in 2D
- **Cell highlight** — yellow outline on hovered cell in all edit modes
- **EditToolbar layout** — sectioned: Off | Terrain group | Biome | Political | Edit
- **CIV tab** — "Generation Parameters" collapsible; "Factions" editor section

---

## Known issue (not yet fixed)

**Political brush paints distant cells** — even with the pick-canvas stroke fix
and chord-distance BFS cap, at low zoom levels in 2D Mercator, painting near
cell boundaries occasionally registers the wrong cell ID. The pick buffer's
RGB color blending at boundaries can still produce IDs that decode to unrelated
cells. The issue is most visible in Map2D at default zoom; the 3D globe is
unaffected (uses ray-cast triangle index, not pick buffer). A fuller fix would
require integer-safe pick rendering (e.g., no anti-aliasing on boundary paths)
or a nearest-cell lookup from cursor coordinates.

---

## Files most recently touched

| File | What changed |
|------|-------------|
| `types.ts` | Added `EditMode`, `PaintStyle`, `UndoSnapshot`; `mountainHeight`/`oceanDepth` to WorldParams |
| `App.tsx` | 10 new state vars + 3 new handlers (`handlePaint`, `handleUndo`, `handleEditFaction`, `handleEditWorldData`); `factionColors` useMemo; Ctrl+Z useEffect |
| `utils/colors.ts` | Exported `PLATE_COLORS`; added `FACTION_COLORS`; `getCellColor` accepts optional `factionColors` 4th arg |
| `utils/worldGen.ts` | Exports `determineBiome`; uses `FACTION_COLORS` (shuffled per civRng) for faction assignment; Stage 9b height remap |
| `utils/paintUtils.ts` | New file — `getCellsInRadius` (BFS + chord-distance cap), `applyTerrainStroke`, `applySmoothStroke`, `applyFlattenStroke`, `applyPoliticalStroke`, `applyBiomeStroke`, `refreshBiomes`, `snapshotCells` |
| `components/EditToolbar.tsx` | New file — floating paint/edit HUD with mode buttons, sub-controls, undo |
| `components/Inspector.tsx` | Extended for world-edit mode: editable faction/province/town fields |
| `components/WorldViewer.tsx` | Orbit gated by editMode; Space key re-enables orbit; `BrushRing` component; cell highlight; right-click for flatten sample; `factionColors` prop |
| `components/Map2D.tsx` | Same-color stroke on pick canvas; Space+pan; brush circle overlay; right-click flatten; `factionColors` prop |
| `components/Controls.tsx` | CIV tab: "Generation Parameters" collapsible; "Factions" inline editor; `onEditFaction` prop |
| `ARCHITECTURE.md` | Updated LLM nav guide, project structure, app state table, new invariants 20–24 |

---

## Key invariants (see ARCHITECTURE.md § Key Invariants for full list)

- **Pick canvas must stroke**: same-color `ctx.stroke()` after `ctx.fill()` in the pick loop or boundary pixels decode to wrong cell IDs (invariant 22)
- **Undo snapshot is a shared Map reference**: pushed at 'start', mutated during 'stroke', never pushed again (invariant 20)
- **`factionColors` useMemo deps on `world`** (not `world?.civData`) — the shallow spread `{ ...world }` creates a new world reference but reuses civData; only the world ref triggers the memo (invariant 21)
- **`FACTION_COLORS` not `PLATE_COLORS` for political/faction** — PLATE_COLORS hues cluster (invariant 23)
- **`CountryLabels` must stay in its own `<Suspense>`** — see invariant 16
- **Multi-field `setParams` must use functional updater** — see invariant 17
- **All state in App.tsx, prop-drilled** — no Context/Redux

---

## Potential next tasks

- **Fix political brush distant-cell painting** — fuller fix requires integer-safe pick buffer (no anti-aliasing at boundaries). Options: render pick canvas with a lower-quality but unblended path (e.g., integer-aligned cells only), or replace RGB-pick with a nearest-cell distance lookup from cursor position.
- **Painting — terrain and political** — terrain painting works but Adaptive vs Freeform visual difference is subtle (Freeform has 2.4× stronger delta + no edge blending). May need further delta tuning.
- **Additional `landStyle` masks** — `'None'` and `'Pangea'` only; `'Ring'` or `'Twin Continents'` is a natural extension.
- **World editor enhancements** — merge/split factions, reorder provinces, bulk-rename.

---

## Workflow notes

- Do not push without explicit user request.
- Commit messages: 50/72 rule, single-line imperative.
- All commits: `Co-Authored-By: <model-name> <noreply@anthropic.com>`.
- User runs dev server via `! npm run dev`. Do not start it unless asked.
- Replace this file entirely when user says "new handoff."
