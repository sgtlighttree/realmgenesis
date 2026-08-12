# F1 — Docked bucket model (ShellApp → WideShell/NarrowShell)

Date: 2026-07-25. Branch: `redesign`. Follows Phase 1 (`useWorldEngine`
extraction) and Phase 2 (ShellApp v1 reframe).

## The fork, decided

Matt's call, 2026-07-25: **migrate ShellApp onto WideShell/NarrowShell** (not
grow the bespoke frame), and **strip each floater's own chrome, wrap in the
shell `Panel`**. Rationale: both branches need the same two hard pieces
(positioning surgery + a real `view` slot), so the cost gap is small, and the
migration keeps layout separate from wiring instead of tangling both in one
~450-line file. The alternative — growing ShellApp in place — was rejected
because it leaves WideShell/NarrowShell to rot as unused prototypes.

## Constraints verified before planning (not assumed)

- **`inspectorCollapsed` is pure presentation.** Only consumer is `Inspector`'s
  `collapsed` prop (App.tsx:131, ShellApp.tsx:172); nothing in `useWorldEngine`
  branches on it. The docked path may pass `collapsed={false}` and ignore it.
- **`sidebarOpen` must stay in the hook.** Classic `App` still drives its drawer
  from it (App.tsx:22,58) even after `NarrowShell` takes over that metaphor in
  the shell path. Deleting ShellApp's hamburger does NOT mean deleting the state.
- **MiniMap's `isCollapsed` is a performance gate, not just presentation** —
  it early-returns the d3 redraw and sits in the effect dep array
  (MiniMap.tsx:17,45). Docking must **conditionally render (unmount)** the card,
  never CSS-hide it, or the 5k-path redraw runs on every paint stroke.
- **MiniMap returns `null` when `!world`** (MiniMap.tsx:47) — build the Read card
  array conditionally or the rail shows an empty titled box.
- **`EditToolbar` is a vertical stack of two already-horizontal rows** — the
  sub-controls row (`flex-wrap max-w-[560px]`) and the mode-button row
  (EditToolbar.tsx:107,111,220). So the Do bar is a chrome-strip + width
  problem, not the layout rewrite it looked like from the class strings alone.

## Explicit assumption (verify in browser, revisit if wrong)

**No per-card collapse in the docked shells.** Collapse exists on Legend/MiniMap
because they float over and occlude the globe; in `WideShell` the Read slot is
already a scrolling column and in `NarrowShell` it is a dismissable sheet, so the
pressure is gone. `Panel` gets NO new collapse mechanism this pass — that would
be adding a feature under cover of removing one. If the rail reads too tall with
Inspector + Biomes + MiniMap stacked, collapse goes into `Panel`'s existing
`right` slot later, implemented once. Classic App keeps collapse unchanged via
the `className` default.

## Steps

### A. `ViewControls` extraction (fills the real `view` slot)

Do this FIRST. Wiring the migration with the `ViewPanel` stub would put a
placeholder in a live app path and force an immediate redo.

New `components/ViewControls.tsx` exporting **shared primitives**, not an
orientation-prop component (an orientation flag is the two-personality smell):

- `DisplayButton` / `ViewButton` / `LayerToggle` — extracted verbatim from
  Controls.tsx:387–412 and the toggle rows at 514–605.
- `VIEW_LAYERS` (the 12 `ViewButton` entries, Controls.tsx:624–635) and
  `LAYER_TOGGLES` (Grid/Rivers/Routes/Hillshade/Contours) as data arrays.
- `ViewStrip` — the horizontal composition for `WideShell`'s top strip.

Each host composes its own layout from those primitives: `Controls` keeps its
current vertical markup in the Sys tab (Render Mode at 435, toggles at 514,
View Layer at 621 — non-contiguous, which is exactly why layout stays per-host),
the strip composes horizontally. Gen params, Map Name, Seed, AI settings, and
World Stats all stay in Make.

Acceptance: classic App renders byte-identically (visual diff, not just gates).

### B. Positioning + chrome surgery on the four floaters

Pattern: **`className` prop defaulting to the current positioning string** — NOT
a `bare` boolean. Classic App and ShellApp v1 keep working with zero changes.

- `Legend.tsx:9` and `MiniMap.tsx:50` — position AND chrome sit on the same root
  div, so one `className` prop covers both axes. Easy.
- `Inspector.tsx:141` — position is on the root, chrome + the `w-28 / w-64 /
  min-w-[220px]` width switching is on the INNER card (142). Two knobs, not one.
  The root is `pointer-events-none` with the inner `pointer-events-auto`; that
  pair is self-consistent, but verify marker/rename/relocate controls still
  respond after docking — a `pointer-events-none` wrapper in a rail is exactly
  the shape that silently kills interaction.
- `EditToolbar.tsx:107` — strip the `absolute bottom-20 left-1/2` root and the
  `bg-black/85 backdrop-blur border shadow-xl` on the two inner rows (111, 220).
  Watch total width against `max-w-[calc(100%_-_18rem)]`.

### C. Migrate `ShellApp` to the shells

`ShellApp` becomes wiring only: `useWorldEngine()` → compose finished elements →
hand to `WideShell` or `NarrowShell` as named slots. The shell never sees props.

- Keep `editing`, the Esc listener, and `closeEdit`'s `setEditMode('off')` in
  `ShellApp`; shells only call `onSetEditing`. `NarrowShell` already derives
  editing from `openTab === 'do'`, so that composes without special-casing.
- Delete ShellApp's hamburger + `sidebarOpen` drawer (NarrowShell owns it).
  Do NOT remove `sidebarOpen` from the hook — classic still uses it.
- Build the `read` array conditionally (MiniMap only when `world`).

### D. Verification — the bar is higher than "gates green"

This pass restructures `EditToolbar`, which drives `handlePaint` / `handleUndo` —
the paths with **zero test and zero browser coverage**. The interactive smoke is
therefore the verification for THIS pass, not a later gate before deleting
classic:

1. Four gates: typecheck 0, lint 0 errors / ≤30 warnings, 138 tests, build OK.
2. Fresh build on `:4180` (`vite preview`) — the standing `:3000` dev server's
   Tailwind JIT does not scan brand-new files' class strings. **Matt's `:3000`
   is his; do not kill it.**
3. Interactive smoke in BOTH folds: paint a stroke, undo it, from the docked Do
   bar. Plus inspect a cell, toggle a layer, switch render mode.
4. Classic App (`/` with no `?shell`) still renders and generates unchanged.
