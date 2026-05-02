# RealmGenesis 3D - Agent Handoff

This file is the quick state transfer for the next session. Read
`ARCHITECTURE.md` for deeper technical context and `AGENTS.md` for repo
workflow/style rules.

---

## Project Snapshot

RealmGenesis 3D is a browser-only procedural fantasy world generator built with
React 19, Three.js/R3F, Canvas2D, d3 geo tooling, and Gemini BYOK lore support.
There is no backend.

- Dev server: `npm run dev` on port 3000
- Quality gates: `npm run build` and `npm run lint`
- No test framework or formatter is configured
- Deployment target: Netlify static SPA

---

## Work Completed In This Session

### Inspector And Selection Performance

- 3D globe no longer raycasts or fetches cell info on pointer move when the
  inspector is disabled and edit mode is off.
- Click inspection now selects cells on click/up only; dragging the globe no
  longer accidentally selects a new cell.
- Selected cells remain highlighted with a brighter fill and thicker outline
  until another cell is selected.
- Hover interaction stays limited to explicit hover-inspect mode.

### 2D Interaction And Dymaxion Picking

- 2D Mercator and Dymaxion click selection now respects the drag threshold:
  dragging pans, clicking selects.
- Empty clicks do not clear the current inspector selection.
- Dymaxion cell picking now uses a hidden color-ID raster generated through the
  same Dymaxion raster pipeline as the visible map. This fixed the drift where
  selecting a Dymaxion cell could select the corresponding Mercator-area cell
  instead.
- Mercator painting/selection avoids blended RGB boundary IDs by inverting the
  projection and using nearest-cell lookup.

### Province-Aware Political Painting

- Political colors again include deterministic province shade variations inside
  each faction.
- Political painting now writes both `regionId` and `provinceId`.
- Political brush strokes lock the target province at stroke start:
  - If the stroke starts inside the selected faction and a valid province, that
    province is used for the whole stroke.
  - Otherwise, the nearest existing province in the selected faction is used.
- Political eraser brush was added; it clears `regionId` and `provinceId`.
- Undo snapshots include `provinceId`, and undo restores faction/province
  ownership together.
- Faction/province population totals are recalculated after political strokes
  and undo.

### Brush And Navigation UX

- 3D edit-mode hover lag was reduced: paint raycasting happens during active
  strokes, not idle hover.
- Middle mouse now rotates the 3D globe.
- Middle mouse now pans 2D Mercator/Dymaxion maps.
- Space+left-drag still remains available as an alternate navigation path while
  painting.

### Faction Overlay

- Added a `Faction Overlay` checkbox in the Sys tab, independent of the active
  view layer.
- Faction borders and labels can now render over biome, satellite, height,
  climate, plates, political, etc.
- 2D Mercator and Dymaxion draw faction borders and labels over the raster map.
- 3D globe draws faction borders as line segments and faction names as curved
  textured mesh patches projected just above the globe surface, so they rotate
  with the planet and are naturally occluded by the mesh.

---

## Files Most Recently Touched

| File | Notes |
|------|-------|
| `App.tsx` | Added `showFactionOverlay`; wired selected-cell highlight, province-aware political strokes, eraser behavior, totals recalculation, undo restore. |
| `components/Controls.tsx` | Added `Faction Overlay` checkbox. |
| `components/EditToolbar.tsx` | Added numbered faction chips and political eraser chip. |
| `components/Map2D.tsx` | Added Dymaxion pick raster, drag-safe selection, middle-mouse pan, 2D faction overlay borders/labels. |
| `components/WorldViewer.tsx` | Added click-only inspection, native paint raycasting, middle-mouse rotate, selected-cell overlay, curved 3D faction labels, toggleable faction borders. |
| `types.ts` | Added `POLITICAL_ERASER_ID`; undo snapshots now include `provinceId`. |
| `utils/colors.ts` | Restored province-derived political color variants. |
| `utils/paintUtils.ts` | Political stroke writes `regionId` and `provinceId`; snapshots include province ownership. |
| `README.md` | Updated feature summary. |
| `ARCHITECTURE.md` | Updated rendering/state/invariants for faction overlay, picking, and province-aware painting. |
| `AGENTS.md` | Added a git workflow note recommending 50/72 commit messages. |

---

## Validation

Last successful checks:

- `npm.cmd run build`
- `npm.cmd run lint`

Lint currently reports warnings only; there are zero lint errors.

---

## Important Invariants

- All app-level state remains in `App.tsx` and is prop-drilled.
- 3D inspection should stay click-only unless `inspectMode === 'hover'`.
- 3D paint raycasting should only run during active strokes to avoid idle hover
  lag.
- Dymaxion picking must use the same raster pipeline as the visible Dymaxion map.
- Political brush undo must preserve `height`, `biome`, `regionId`, and
  `provinceId`.
- Political brush should not create provinces; it assigns cells to an existing
  selected-faction province or clears ownership when using the eraser.
- `getCellColor(..., factionColors)` must receive the live `factionColors` map
  on render paths that display political ownership.
- 3D faction labels are curved textured meshes, not HTML overlays or billboard
  sprites.

---

## Potential Next Tasks

- Manual browser regression pass across 3D, 2D Mercator, and 2D Dymaxion after
  this commit lands.
- Tune 3D curved faction label sizing/placement if very long faction names are
  introduced.
- Consider reducing existing lint warnings when there is a quiet maintenance
  window.
- Broaden editor features later: merge/split factions, province management,
  bulk rename, town/capital relocation.

---

## Workflow Notes

- Do not push without explicit user request.
- Git commit messages are strongly recommended to follow the 50/72 rule:
  subject line at or under 50 characters, body wrapped near 72 characters.
- Keep commits focused and imperative.
- The user normally runs the dev server; only start it when useful for manual
  verification or when asked.
