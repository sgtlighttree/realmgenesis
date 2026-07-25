# RealmGenesis 3D - Agent Handoff

This file is the quick state transfer for the next session. Read
`ARCHITECTURE.md` for deeper technical context and `AGENTS.md` for repo
workflow/style rules.

## Matt's (maintainer's) notes

> Matt's scratchpad and notes for things observed outside an active coding session. If an item is addressed, click the checkbox, and/or add a ~~strikethrough~~ for emphasis.

- [ ] Make a true vector 2D mode instead of raster, but keep it optimized
- [ ] V3 of terrain generation algorithm. Goal is to make plate boundaries far more realistic, make part of Milestone D.
- [ ] Major UI/frontend/rendering overhaul (Milestone F), use skill `/impeccable` for visual UI review
- [ ] Major feature, for much, MUCH later: World Formats: Planet, Flat Earth (Disc, Rectangle, etc.)
- [ ] Add a favicon just to clear the constant 404'ing

---

## ⚡ NEW-THREAD PICKUP (updated 2026-07-24, end of Session 4)

**C3 (roads & trade routes) SHIPPED this session** — the last pre-D6
additive feature. The whole C-tier and pre-D6 batch are now done. A fresh
thread picks up at the **big-rock planning phase**:

1. **D6 / vector-2D / A3 as ONE rendering-contract decision** (see the D6/F1
   sequencing analysis below — that framing still holds). This is a
   COMMITMENT BOUNDARY: brainstorm + advisor-consult before writing code.
2. **F1 (UI overhaul)** — may come before/alongside D6; needs Matt's design
   input, use `/impeccable`. C-tier UI was kept deliberately minimal (C3
   added exactly one "Roads & Routes" toggle) precisely to limit F1 rework.

The spec + plan for C3 live at `docs/superpowers/specs/2026-07-24-c3-*.md`
and `docs/superpowers/plans/2026-07-24-c3-*.md` (brainstorming → writing-plans
→ executing-plans workflow; useful template for the next feature).

**Execution-mode note for Session 4:** Matt directed inline/self execution
(no subagent delegation) because C3 was a serial one-file-at-a-time chain,
and codified that as a new CLAUDE.md clause. Delegation stays the default
for parallelizable work; skip it when serial.

### Session 3 delegation protocol (working policies, also in memory)

- **Sonnet 5 subagents by default** — Matt's directive. Opus only if
  unavoidable (and then he wants 4.6; the Agent tool can't pin versions, so
  flag to him instead of silently using another Opus). Subagent spend limit
  was hit once mid-session (killed the A4 agent mid-task) — if agents fail
  with a spend-limit error, finish the work inline from the brief.
- **Briefs carry ALL design decisions** (exact files, integration points,
  acceptance criteria incl. "lint ratchet exactly 30 warnings, add none",
  "do not commit", "do not touch HANDOFF/CLAUDE/ROADMAP"). Sonnet's
  literalness is an asset with a complete brief.
- **One agent at a time** — every feature funnels through App.tsx/
  Controls.tsx (prop-drilled architecture); parallel agents collide.
- **Fallback heartbeat**: alongside each agent launch, arm a ~40-min Monitor
  timer (`sleep 2400; echo ...`). Agent finishes first → TaskStop the timer.
  Timer fires first → SendMessage the agent for status. Never heartbeat with
  no agent running. (Cache economics: at this context size one miss ≈ ten
  warm turns.)
- **Orchestrator verifies everything**: re-run all four gates yourself,
  read the key diffs, browser-verify via Playwright (dev server on :3000;
  synthetic clicks need MouseEvent not PointerEvent for Map2D picking),
  commit in logical chunks with 50/72 messages. Do NOT push (standing rule).

### Post-milestone tier — SHIPPED this session (commits 47ef94f..f0459a0)

| Feature | Commits | Notes |
|---|---|---|
| A4 hillshading + contours | 47ef94f | Relief-only Lambert shade map (no terminator), cell-edge isolines, toggles + export. Agent died at spend limit; finished inline. |
| A5 geodesic ruler + scale bar | a78c60c | measure.ts pure math; ruler intercepts onInspect (children untouched); projection-aware scale bar (project-2-points method); agent caught a Map2D blit-deps bug itself. |
| E1/E2 SVG + GeoJSON export | 3461990 | Layered SVG (mirror on geo groups, counter-mirrored text); RFC 7946 FeatureCollection; validated with xmllint + python beyond the suite. |
| C4 markers/POIs | b429b83 | Sphere-position-anchored (survive regen), 'marker' LabelKind through shared pipeline, save/load with sanitizer; agent caught a pin double-mirroring bug. |
| C5 civ editor ops | 117a0d5 | mergeFactions (full province-id map built BEFORE cell rewrite), renames, capital relocation (dual isCapital flag pair). Split deferred. |
| C1 cultures | 09f4bdf, f0459a0 | Terrain-affinity Dijkstra cultures on '_cultures' stream (civRng untouched — liveness-proven); per-culture namebase styles drive faction/town naming by capital's culture. Browser-verified: NAJRA/ZAGHATI (desert) beside VESTAD/Isgard (norse). |

Suite: 52 → 119 tests across the tier. Every feature: typecheck 0, lint
0 errors / exactly 30 warnings (ratchet — do not exceed), build OK.

### D6 / F1 sequencing analysis (agreed with Matt)

- D6 (terrain V3: realistic plate boundaries, sub-cell heightmap detail)
  breaks VALUES not INTERFACES for most features — derived layers (civs,
  cultures, routes) regenerate by design. True wait-list: D4 submaps
  (reuses the generator), A3 raster-heavy styling, B2 resurvey semantics,
  D1–D3 tuning. The D6 planning phase should absorb THREE things as one
  rendering-contract decision: terrain V3 + Matt's vector-2D note + A3.
- F1 (UI overhaul, Matt's addition): may come before or alongside D6.
  Deliberately NOT started in full-auto — needs Matt's design input. C-tier
  UI additions were kept minimal (buttons/selects) to limit rework.
- Pre-D6 batch order was: A5 → E1/E2 → C4 → C5 → C1 → C2 → C3 (all shipped
  except C2 in-flight, C3 next).

---

## Session 6b (2026-07-25) — impeccable audit + critique, then fixes

Same session, after the docked bucket model landed. Commits `3d5a500`..`c4baa75`.
Gates throughout: typecheck 0, lint 0 errors / 29 warnings, 138 tests, build OK.
Critique snapshot: `.impeccable/critique/2026-07-25T03-53-02Z__components-shell.md`.

**Scores: Design Health 22/40, Audit Health 10/20.** Slop verdict: not slop at a
glance, borderline under 30s of clicking — the tells are interaction-level
(three control vocabularies for peer decisions), not visual.

**The detector's 9 findings were ALL false positives.** Every one is a ternary
`active ? 'bg-blue-600 text-white' : 'bg-gray-800 text-gray-400'`; the
`gray-on-color` rule scans the whole template literal and pairs the inactive
branch's gray text with the active branch's saturated bg — a pairing that never
renders. Don't "fix" these; the rule is wrong, not the code.

**Decisions + rationale:**

- **Themed form controls are NOT a violation of the product register's "don't
  reinvent standard affordances (custom scrollbars, weird form controls)" ban.**
  That ban protects standard *behaviour*. A native `<select>` on a near-black app
  renders a light OS menu — that is a theming defect, not an affordance being
  honoured. So: appearance is ours, behaviour is the native contract. `Select`
  keeps type-ahead, Home/End, arrows, Enter/Space, Esc-cancels-without-commit,
  focus-returns-to-trigger, and full ARIA listbox semantics (all browser-verified).
- **`Select` portals to `document.body` and positions `fixed`.** An absolutely
  positioned menu is clipped by the View strip's overflow and the mobile sheet's
  `overflow-auto`. Do not "simplify" this back to absolute.
- **`ConfirmDialog` is built on the native `<dialog>` + `showModal()`** for focus
  trapping, page inertness, and Esc — not a hand-rolled portal+trap. `window.confirm`
  was rejected for the same reason as the native select: OS chrome in a dark app.
- **The generate gate lives in `useWorldEngine`, not the button**, so every entry
  point inherits it. Fires only when `undoStack.length > 0`. **Auto-update stays
  ungated on purpose** — it fires on every slider change; prompting would be
  unusable.

**Two REAL bugs found by looking at rendered pixels, not code:**

1. **`rg-rise` was clobbering Tailwind transforms** (fixed `f95def3`). The keyframes
   animated `transform`, and `both` fill mode left `transform: translateY(0)` on the
   element permanently, silently overwriting `-translate-x-1/2`. The Do bar sat half
   its own width right of centre, overlapping the Read rail, ever since the docked
   layout landed. Fix: animate the independent **`translate`** property, which
   composes with `transform` instead of replacing it. **If you add a rise/slide
   animation to anything Tailwind also transforms, use `translate`, not `transform`.**
2. **The mobile Make sheet opened on the console** (fixed `e38f22a`). Two causes, and
   the obvious one was the minor one. `ConsoleOutput` called
   `scrollIntoView({behavior:'smooth'})`, which **walks up and scrolls every
   scrollable ancestor**, dragging the sheet down to the log box — and because it is
   smooth, it lands *after* any scroll reset. Resetting `scrollTop` alone still left
   ~150px of drift; that was a **refuted first diagnosis** (I initially blamed only
   the sheet body being one reused DOM node across tabs, which is real but secondary).
   Fix: `ConsoleOutput` sets `scrollTop` on its own box; `NarrowShell` also resets on
   tab change. Verified `scrollTop: 0` with Render Mode in view on first open and
   re-open.

**Still open (Matt scoped these out of this pass):**

- **a11y is the lowest score (1/4) and the biggest available gain.** 58 buttons, only
  5 with semantic ARIA; 44 rely on `title`, which is a tooltip, not an accessible
  name, and never appears on touch. The 15+ biome swatches and faction chips read as
  "button, button, button" to a screen reader.
- **Touch targets ~22px vs the 44px minimum** on the strip chips and EditToolbar
  modes — these are primary controls on the narrow fold.
- **Compounded contrast:** `text-[9px] text-gray-600` on `bg-black/85` is the worst
  case; 30 uses of `text-gray-500`, 14 of `text-gray-600` overall.
- **No token layer** — 18 hard-coded hex, bg opacity spread across /10–/85, informal
  `z-10`/`z-20`. Worth fixing *before* F1b, or the brand pass is a find-and-replace
  across dozens of files instead of one token block.
- **`shellKit`'s stub panels** (`MakePanel`/`ViewPanel`/`DoPanel`/`READ_CARDS`) are
  now used only by `?shell=stub`. Decide whether they ship.

---

## Session 6 (2026-07-25) — F1 docked bucket model SHIPPED

Branch `redesign`, commits `4d804fb`..`5b80764`. Gates: typecheck 0, lint 0
errors / **29** warnings (ratchet was 30 — see below), 138 tests, build OK.
Browser-verified on a fresh `:4180` preview in BOTH folds. Plan doc:
`docs/superpowers/plans/2026-07-25-f1-docked-bucket-model.md`.

**⚠️ DATA-LOSS INCIDENT (read this first).** This session opened with
`HANDOFF.md` clobbered in the working tree: a stale editor buffer had
overwritten it with a pre-Session-5 version, deleting **204 lines** — the whole
F1 record, the full-auto log, the next-pass notes, and the advisor resume notes.
It was uncommitted, so `git checkout` recovered it; Matt's two scratchpad edits
that rode along in the same save were re-applied on top (commit `4d804fb`).
**Lesson: if you have HANDOFF.md open in an editor across a session where an
agent also writes to it, reload the buffer before typing.** The recovery was
committed IMMEDIATELY rather than left in the working tree, because a warning
is not durable and the same buffer was still open.

**The fork, decided (Matt, this session).** The advisor's Session-5 resume note
framed the next step as an explicit fork. Matt chose:

- **(a) migrate ShellApp onto WideShell/NarrowShell**, not (b) grow the bespoke
  frame. Rationale: both branches needed the same two hard pieces (positioning
  surgery + a real `view` slot), so the cost gap was small, and (b) would have
  tangled layout with wiring in one ~450-line file while the approved shells
  rotted unused.
- **Strip each floater's chrome and wrap in the shell `Panel`**, not "render
  bare with its own chrome". This is the "one chrome" consistency win F1 exists
  for.

**Decisions + rationale (the perishable part):**

- **`ViewControls` exports PRIMITIVES, not an orientation-flag component.** The
  Sys tab interleaves render-mode / toggles / view-layer *vertically between
  Make content* (Render Mode at the top, then Seed and points, then toggles,
  then View Layer) while the strip lays them out inline. Layout is the part that
  genuinely differs per host, so each host composes it; only the buttons, toggle
  definitions, and layer list are shared. An `orientation` prop would have been
  the same two-personality smell the advisor flagged for `bare` booleans.
- **`className`-defaulting-to-current-string, NOT a `bare` boolean**, for every
  parameterized component. Classic App keeps working with zero changes.
- **No per-card collapse in the docked shells (explicit assumption, verified).**
  Collapse existed on Legend/MiniMap because they *occlude the globe*; in
  WideShell the Read slot already scrolls and in NarrowShell it is a dismissable
  sheet, so the pressure is gone. Adding collapse to `Panel` would have been
  adding a feature under cover of removing one. If the rail ever reads too tall,
  it goes in `Panel`'s existing `right` slot, implemented once.
- **`MiniMap.isCollapsed` was a PERFORMANCE gate, not presentation** — it
  early-returned the d3 redraw and sat in the effect dep array. The docked card
  therefore conditionally **unmounts** `MiniMapCanvas`; CSS-hiding it would
  redraw 5k paths on every paint stroke. Same reason the Read array is built
  conditionally (MiniMapCanvas returns null without a world → empty titled box).
- **Lint ratchet moved 30 → 29 legitimately.** The extracted `ViewButton` no
  longer needs the `icon: any` the inline version carried. Do not "restore" 30.
  (Mid-flight it spiked to 43: moving the icons out of Controls left 14 dead
  lucide imports. Fixed, not suppressed.)
- **`showHeader` on Controls**: the brand header is *shell* chrome. The shell
  rail draws its own, so docked Controls rendered "RealmGenesis 3D" twice.

**Verification bar was raised deliberately.** This pass restructured
`EditToolbar`, which drives `handlePaint`/`handleUndo` — the paths with zero
test and zero browser coverage. So the interactive smoke was the gate for THIS
pass, not a later one before deleting classic. Verified with **trusted** pointer
events (`page.mouse`, not synthetic dispatch, which R3F raycasting ignores).

Exactly what was covered, so the next session doesn't over-trust this:

- **Wide fold, 3D:** paint stroke in the docked Do bar took undo 0 → 1, undo
  took it 1 → 0, Esc exited edit mode and cleared paint mode. Cell inspection
  populated the docked Inspector with live data and its header tool buttons
  still respond — the `pointer-events-none` wrapper risk did NOT materialize.
- **Wide fold, 2D:** Mercator and Dymaxion both render through the shell's
  `canvas` slot; the Map2D canvas measures 992×720 with CSS size == attribute
  size, so the Dymaxion pick-buffer-mirrors-raster invariant still holds. A
  `ViewStrip` layer chip flips, and the Sys-tab checkbox mirrors it (shared
  state confirmed across both hosts of the extracted primitives).
- **Narrow fold:** Make/View/Do/Read tab bar; Do sheet is edit mode; the Make
  sheet scrolls (340px viewport over 1534px content) with Generate World
  reachable above the tab bar. The `h-full`-inside-a-percentage-capped-sheet
  collapse this was expected to hit did not occur.
- **Classic** re-verified: single brand header, all four collapse chevrons, all
  toggles, no tab bar. **`?shell=stub`** re-verified after the `WideShell`
  `bodyClassName="p-0"` change — still renders A·Tidy correctly.
- **NOT covered:** save/load, lore/`apiKey`, abort-mid-generate, and painting in
  2D projections. Those remain the pre-existing coverage gap.

**Known nits, deliberately not fixed (need Matt's eye / out of scope):**

- **`WorldViewer`'s own pause control (`absolute top-4 right-4`) now sits under
  the Read rail.** Pre-existing element, newly colliding. Fixing it means
  touching 3D presentation chrome — that's F2 territory.
- **Ruler readout (`bottom-6 left-1/2`) overlaps the NarrowShell Do sheet** when
  the ruler is active and the sheet is open. Transient-tool overlap only.
- `?shell=stub` still serves the DesignShell prototype; `?globe=0` still works
  and still only skips WebGL, not generation.

**Next:** classic App and ShellApp are still a fork sharing one hook — mirror
component-wiring changes in both, and retire classic once ShellApp reaches
parity (`index.tsx` routing is the switch). Then F1b, the `/impeccable` visual
identity pass on the settled skeleton.

---

## Session 5 (2026-07-24) — F1 UI redesign: shell prototype BUILT (not merged)

Brainstormed the F-tier + D6, then built a working layout prototype. **On branch
`c3-roads-trade-routes` still (F1 work is NOT yet on its own branch — separate
before committing).** Gates green: typecheck 0, lint 0/30, build OK. Browser-
verified against a fresh `vite preview` build (see gotcha below).

**Decisions + rationale (the perishable part):**

- **Relief hinge resolved (was the crux tying D6/F2/F3 together).** Terrain relief
  can live in *geometry* (displaced mesh, today) or in *texture* (smooth sphere +
  hillshade). Matt's call: **smooth sphere is the DEFAULT across all view layers**
  (better legibility of roads/rivers/borders); displaced-geometry is a *separate
  toggle* for later, where line overlays hug the mesh. **D6 is decoupled** — pure
  gen-algorithm work (realistic plate boundaries, kill seams), NOT a presentation
  decision. F3 = Google-Maps-style vector 2D. This un-blocks treating them
  separately instead of one mega-decision.
- **Sequence: F1 first** (Matt wants the dopamine of a new frontend; it's also the
  most separable — presentation-layer, doesn't care how relief is carried).
- **F1 scope = layout rearchitecture + consistency cleanup** (unify panel chrome,
  one slider-accent color, themed scrollbars). NOT a new visual identity — that's
  **F1b**, a later dedicated `/impeccable` pass on the clean skeleton. Rationale:
  structural work and taste work are different kinds of hard; you can't judge
  type/color until the skeleton is settled and populated.
- **Architecture — content/shell decoupling (advisor-confirmed, preserves the
  no-Context invariant).** I initially feared two shells would force prop-drilling
  through two layout trees or a Context. The advisor corrected the framing: **App
  composes each panel WITH its props into a finished element** (`panels = { make:
  <MakeContent {...}/>, … }`) and hands the map to whichever shell is active as
  named slots; **the shell only POSITIONS pre-built elements, never sees their
  props.** Props still trace one hop to the owner. No Context, no store. Corollary:
  ephemeral presentation state (mobile "which tab is open") lives in the shell's
  own `useState`, NOT App — hoisting it would be the real invariant violation.
- **Layout: A-shell wide, C-shell narrow, same content.** Matt judged A·Tidy best
  on desktop, C·Studio best on mobile. Because content is decoupled from shell,
  "A wide + C narrow" costs ~the same as "C both" (C already needs two shells:
  desktop docks vs mobile tab-bar diverge *behaviorally*, not just by CSS). So no
  compromise needed. Mobile finding: all three wireframe directions converge on
  "globe-hero + sheets" on a phone, so the phone layout is really its own design
  and the desktop pick barely constrains it.
- **Edit mode = single "Edit" toggle, default OFF, summons the contextual Do bar;
  Esc exits.** Inverts today's always-visible EditToolbar (whose first pill is
  "Off"). Keeps Read (click-to-inspect, always on) and Do (paint, modal) from
  stepping on each other. On narrow-C, "Do" is one of the four bottom tabs — same
  metaphor, no special case.
- **F1 STUBS the render modes; it does NOT build F2/F3** (advisor guardrail). The
  placeholder globe + smooth/relief View toggle are stubs so the spec doesn't
  quietly absorb the rendering rework.

**What was built (throwaway-safe prototype, reachable via `?shell=1`):**

- `index.tsx` branches on `?shell` → mounts `DesignShell` instead of `App`.
- `components/shell/shellKit.tsx` — placement-agnostic stub panels (MakePanel,
  ViewPanel, DoPanel, READ_CARDS, PlaceholderGlobe), the `ShellProps` slot
  contract, and the single `PANEL` chrome constant (the "one chrome" the
  consistency cleanup buys — change radius/border/fill in one place).
- `components/shell/WideShell.tsx` (A·Tidy) and `NarrowShell.tsx` (C·Studio).
- `components/shell/DesignShell.tsx` — harness: Auto/Wide/Narrow override toggle
  (preview the fold without resizing), editing state + Esc, composes stubs once.
- **NOT wired to real state.** Panels are dumb stubs; globe is a CSS circle. The
  real F1 = making the actual Inspector/Legend/MiniMap/EditToolbar/Controls
  placement-agnostic (they currently self-position with `absolute …`) and mounting
  them through these shells. That refactor is the next step, not done here.

**Impeccable refinement pass applied (same session, product register).** Removed
three self-inflicted AI tells: per-panel uppercase-mono eyebrow tags, the 4-hue
"rainbow" bucket dots (→ single blue accent, state/selection only), and default
glassmorphism (→ solid `bg-gray-900` panels). Added: state-motion `rg-rise`
(`index.css`, ease-out-quart, `prefers-reduced-motion` fallback) on the Do bar +
mobile sheet; an exported `FOCUS` ring on every interactive control; muted-text
contrast bump. Narrow sheet capped to `max-h-[52%]` (bottom sheet, globe stays
visible). These are chrome/token decisions that carry straight into the real F1.
Skill v3.9.1 installed; v4.0.2 update was offered to Matt (not yet taken). The
full VISUAL identity (type pairing, palette) is still deferred to F1b.

**Gotcha (verified, n=1 but cleanly reproduced):** a LONG-RUNNING `vite` dev
server's Tailwind JIT does **not** pick up brand-new files' unique classes —
`top-3`/`right-3`/`right-[17rem]` were absent from generated CSS (only classes
already used elsewhere rendered), collapsing the layout. `npm run build` +
`vite preview` (or restarting dev) fixes it. **Verify new-file UI against a fresh
build, not the standing dev server.** Also: Tailwind arbitrary `calc()` needs
underscores for spaces — `max-w-[calc(100%_-_18rem)]`, not `-18rem)]`.

**F1 implementation spec written + advisor-reviewed:**
`docs/superpowers/specs/2026-07-24-f1-ui-redesign-design.md`. Covers the
`useWorldEngine()` hook extraction, `Controls` decomposition, the §3.3-A
shared-component positioning decision, the `?globe=0` toggle, and a 5-phase plan
(1 hook extract → 1b positioning extract → 2 shell v1 playable → 3 decompose
Controls → 4 polish). **Two phases (1, 1b) touch the classic app**, each with its
own regression pass; the extraction freezes App's return block (zero-char diff)
so TypeScript + the frozen render carry most of the fidelity proof — manual review
is scoped to effect deps + ref timing. Classic layout is transitional (retired at
shell parity). Second advisor consult DONE for the wiring architecture.

**Phase 1 plan:** `docs/superpowers/plans/2026-07-24-f1-phase1-useworldengine.md`.

### F1 execution log (FULL AUTO, Matt AFK; rate reset ~02:39 2026-07-25)

Matt authorized full-auto execution of F1 Phases 1–4: commit per chunk, update
HANDOFF as I go, delegate only genuinely-parallel well-scoped chunks (most of this
is serial through App/Controls, so mostly self; Phase 3 Controls split is the one
real delegation window). Verify against a FRESH build on `:4180` (dev-server
Tailwind won't hot-scan new files). Matt's `:3000` dev server is his — do not kill.

- **Phase 1 — DONE** (commit `0c373f4`). `useWorldEngine()` extracted verbatim via
  sed (94 exports; return type = `ReturnType<typeof useWorldEngine>`; 4 refs stay
  internal). App = thin consumer, return block byte-identical (diff-verified).
  Gates: typecheck 0, lint 0/30, 138 tests, build OK. Classic app browser-verified
  (generate → full world, labels/colors, zero console errors). Paint/undo drag NOT
  synthetically driven (flaky on 3D) — code moved verbatim; will validate when
  Phase 2 wires EditToolbar into the shell.
- **Phases 1b + 2 RESEQUENCED (full-auto judgment call — decision + rationale):**
  While wiring Phase 2 I found the floaters' fixed viewport anchors (Inspector
  `top-6 left-1/2`, MiniMap `bottom-4 right-4`, etc.) COLLIDE with the shell's View
  strip / Read rail if dropped into the bucket slots. Docking them cleanly needs
  component surgery (each carries its own chrome + centering + `pointer-events`),
  which is aesthetically sensitive — bad to finalize blind while Matt's asleep.
  **So: ship a safe playable increment first** — `ShellApp` = real data in a clean
  reframe with the floaters kept as-is (no collision) + the `?globe=0` toggle.
  **Deferred to a Matt-present pass:** the docked bucket model (bare Read cards,
  View strip, contextual Do via edit toggle) reusing WideShell/NarrowShell, which
  is where the component-docking aesthetics need his eye. `?shell=stub` still
  serves the DesignShell prototype for that layout reference.
- **Phase 2 (v1 reframe) — DONE** (commits `a3e7702`, `a919610`). `ShellApp`
  (`components/shell/ShellApp.tsx`) is the F1 redesign entry, consuming the same
  `useWorldEngine` hook as classic App. Delivered + browser-verified on a fresh
  `:4180` build:
  - **`?shell=1`** → real, playable redesign entry (generate/orbit/inspect all
    work); **`?shell=stub`** → DesignShell prototype; else classic App.
  - **`?globe=0`** → swaps the Three.js globe for `PlaceholderGlobe` + a mode
    banner, hides floaters (fast UI iteration, no WebGL cost). Works via
    `?shell=1&globe=0`.
  - **Contextual Do bucket DONE:** EditToolbar is hidden by default, summoned by
    an "Edit" toggle (top-right, amber when active), dismissed with Esc / toggle
    (which also clears paint mode). Local `editOpen` state — no engine change.
  - v1 keeps the four floaters self-positioned over the canvas (like classic) —
    NO component surgery, NO collisions.

### NEXT PASS (do with Matt present — aesthetically sensitive): docked bucket model

Goal: replace v1's floaters with the approved A-wide/C-narrow docked layout
(`WideShell`/`NarrowShell`). Learnings from tonight that make it fast:

- **The blocker is collisions:** the floaters' fixed viewport anchors (Inspector
  `top-6 left-1/2`, seed HUD `top-4 left-24`, MiniMap `bottom-4 right-4`, Legend
  `bottom-4 left-4`, EditToolbar `bottom-20 left-1/2`) collide with any View strip
  / Read rail. So docking Read (right) and adding the View strip (top) must happen
  together — you can't add the top strip while Inspector floats top-center.
- **Positioning surgery — use the `className`-override-with-default pattern, NOT a
  `bare` boolean** (the boolean is the two-personality smell the advisor flagged;
  a `className` prop defaulting to the current positioning string is idiomatic and
  keeps classic + ShellApp-v1 working with zero changes). Roots to parameterize:
  - Inspector `components/Inspector.tsx:141` — pos: `absolute top-6 left-1/2 -translate-x-1/2 z-10`; internal: `flex flex-col items-center gap-2 pointer-events-none`.
  - Legend `components/Legend.tsx:9` — pos: `absolute bottom-4 left-4 z-10`; internal: `bg-gray-900/80 backdrop-blur border border-gray-700 shadow-xl transition-all duration-300`.
  - MiniMap `components/MiniMap.tsx:50` — pos: `absolute bottom-4 right-4 z-10`; internal: `bg-black/80 border border-gray-700 shadow-2xl overflow-hidden transition-all duration-300`.
  - EditToolbar `components/EditToolbar.tsx:107` — pos: `absolute bottom-20 left-1/2 -translate-x-1/2 z-20`; internal: `flex flex-col items-center gap-1 pointer-events-auto select-none`.
- **Chrome reconciliation:** these components bring their OWN bg/border/shadow. In
  the shell Read slot, either (a) render them bare of the shell `Panel` wrapper
  (accept their own chrome — fastest, slightly inconsistent), or (b) strip their
  chrome too and wrap in `Panel` (cleaner, more surgery). Decide with Matt.
- **View strip content:** split `Controls` (1571 lines) — render-mode + layer
  toggles → a new `ViewControls`; gen params stay in Make. This is Phase 3 and the
  one delegable chunk (tight `sonnet-medium` brief, verify no visual change).
- **Reuse:** `WideShell`/`NarrowShell` already take `make/view/read/doTools/canvas`
  slots; wire `ShellApp` to them. The contextual-Do toggle logic is already built
  in ShellApp — port it to the shell's Edit affordance.

### Full-auto session summary (2026-07-25, ~00:00–early AM)

Shipped on `redesign`, all gates green throughout (typecheck 0, lint 0/30, 138
tests, build OK), each step browser-verified on a fresh `:4180` preview:
`e6dc6ee` spec → `bcddffc` Phase-1 plan → `0c373f4` hook extraction →
`a3e7702` ShellApp+?globe=0 → `a919610` contextual Edit toggle. Preview server may
still be running on `:4180` (static build of this state — kill with
`lsof -ti:4180 | xargs kill`; for live dev restart `:3000`). NOT pushed. Classic
App verified unchanged. Deferred the docked bucket model (above) for Matt's eye.

**Advisor resume notes (read before the next pass):**

1. **ShellApp is a BESPOKE frame — it does NOT use WideShell/NarrowShell yet.**
   The "reuse the shells" line above is aspirational, not current. So the next
   session's FIRST decision is an explicit fork: **(a) migrate ShellApp onto
   WideShell/NarrowShell** (realizes the approved A/C bucket layout, more rework)
   **vs (b) grow the bespoke frame** (keep ShellApp's current structure, add a
   docked Read rail + View strip in place). Pick deliberately; don't drift.
2. **ShellApp duplicates classic App's render JSX** (both consume the one hook).
   Until classic is retired they're a fork — mirror any component-wiring change in
   both, and **retire classic promptly once ShellApp reaches parity** to end the
   fork. `index.tsx` routing is the switch.
3. **Interactive paths are UNVERIFIED** (not by tests, not by browser). The 138
   tests cover the pure engine; the browser smoke only covered *generation*. So
   `handlePaint`, `handleUndo`, abort-mid-generate (`abortControllerRef` never
   actually fired — first generate has no prior controller), lore/`apiKey`, and
   save/load have zero coverage. Verbatim-move keeps them low-risk, but **run a
   ~3-min interactive smoke (paint a stroke, undo, regenerate mid-run, save/load)
   as the GATE before deleting classic** — that deletion is the irreversible step.
4. **`?globe=0` skips WebGL, not generation.** A world still auto-generates on
   mount (seen: full gen ran under globe=0, seed `realmgenesis`); confirm the
   auto-gen source. The "fast UI iteration" benefit is real but partial (no globe
   render/interaction cost; generation still runs).

---

## Session 4 (2026-07-24) — C3 roads & trade routes SHIPPED

Committed on branch `c3-roads-trade-routes` (6 feature commits + spec/plan
docs; NOT pushed, NOT merged to main — Matt fast-forwards when ready).
Final state: typecheck 0, **138 tests**, lint 0 errors / exactly 30 warnings,
build OK. Browser-verified on 3D globe + 2D Mercator + 2D Dymaxion, plus PNG
(4K), SVG (xmllint-clean), and GeoJSON (33 road + 3 searoute features) export.

**Decisions + rationale (the perishable part):**

- **Roads are a FOREST, not one MST per faction.** The advisor caught a real
  bug in the first design: a per-faction great-circle MST can route A–island–B,
  then A* drops both water edges and strands two *mainland* towns that should
  share a road — and it contradicts the connectivity test. Fix: BFS the
  land-connected components first, build one MST per `(faction, land-component)`
  group. The road network is a forest; sea routes bridge the trees. The
  `tests/routes.test.ts` "forest invariant" asserts acyclicity + per-group
  spanning tree (trunk roads excluded — they're cross-faction and can cycle).
- **New `utils/pathfinding.ts`** to avoid a circular import: `MinHeap` (moved
  out of worldGen), `isWaterCell`, extracted `landTerrainStepCost`, and `aStar`.
  Clean DAG: `types ← pathfinding ← {worldGen, routes}`; `worldGen → routes`.
  The civ Dijkstra now consumes `landTerrainStepCost` — **identical by
  construction** (same ops, order, operands; verified by inspection). The
  determinism suite stayed green, but note that only proves self-determinism
  of the new code, NOT equivalence to pre-refactor output — don't trust a green
  suite alone to catch a value-changing refactor.
- **Routes are computed LAZILY (App effect), not at generation.** This is the
  fix for a real regression the advisor caught: `computeRoutes` is O(towns·A*)
  and measured **90ms@20k, 1.8s@60k, 3.6s@120k** cells — and it originally ran
  unconditionally at the tail of `recalculateProvinces`, freezing the main
  thread on *every* generate (even routes-off, the default) AND after the
  progress bar already hit 100%. Now `recalculateProvinces` clears
  `world.routes`; an App `useEffect` computes them only when the toggle is on
  and none are cached (30ms yield + `setIsGenerating` spinner, mutate +
  shallow-copy like paint). Routes-off generations pay zero. Interactive safety
  checked: only the explicit "Update Civs/Provinces" buttons (which already show
  a spinner) reach `recalculateProvinces`; political paint strokes do NOT, so
  no per-stroke route recompute. Tests compute routes explicitly to match.
- **Sea routes use a per-pair goal-scoped `seaStep`** (water cheap, land
  impassable except the destination port) — keeps routes on water and blocks
  cutting across intermediate ports on land. This also sidestepped a `majorSet`
  temporal-dead-zone trap the plan had flagged. Improvement over the plan's
  `majorSet` closure; noted here so it isn't "corrected back" later.
- **Determinism** rides on stable insertion order (same as existing Dijkstra)
  + explicit MST edge tie-break `(weight, minCellId, maxCellId)`. No RNG needed;
  the `civSeed + '_routes'` stream exists only as a hook for future tie-breaks.
- **Dashed sea routes in 3D:** `LineDashedMaterial` needs a `lineDistance`
  attribute, and `LineSegments.computeLineDistances()` isn't reachable through
  the R3F string-element (`'lineSegments' as any`) path — so we build the
  attribute manually in `buildRouteGeometry`. Also: the extra dashed-material
  string const was cast `as typeof LineSegments` (already `any`) rather than a
  fresh `as any`, so the lint ratchet stayed at exactly 30 (no keyword added).

**Finding (n=1, worth knowing):** the raster PNG export (`export.ts`
`renderEquirectangular` / inline `exportMap`) drew **no rivers at all** — the
old HANDOFF note about "rivers in export" was stale. Routes are therefore the
first hydrology-adjacent overlay in PNG (gated on `showRoutes`). SVG export
already had rivers as first-class layers; routes join them there too.

**GLB export omits routes (follow-up, not oversight).** `exportGLB.ts` exports
rivers as line geometry but was scoped out of C3 (PNG/SVG/GeoJSON only), so GLB
is now the one surface where rivers appear and routes don't. Add route line
geometry to the GLB exporter when convenient — small, mirrors the river path.

**Deferred, on purpose:** route connectivity → town importance/population
(ROADMAP wants it; it makes C3 non-additive by feeding civ generation).
`RouteData.fromCellId/toCellId` are stored now so it's a small future step.
**Tuning knob, non-blocking:** "nearest 3 major ports" can draw short sea hops
paralleling a coast road; dedup against road-connected pairs or set a min
crossing distance if it ever reads as clutter (Matt picked the dense web).

---

## Project Snapshot

RealmGenesis 3D is a browser-only procedural fantasy world generator built with
React 19, Three.js/R3F, Canvas2D, d3 geo tooling, and Gemini BYOK lore support.
There is no backend.

- Dev server: `npm run dev` on port 3000
- Quality gates (all enforced in CI): `npm run lint`, `npm run typecheck`,
  `npm test`, `npm run build`
- Vitest engine suite in `tests/`; no formatter is configured
- Deployment target: Netlify static SPA

---

## Work Completed In This Session (Session 3 - 2026-07-18/19)

### "Map identity" milestone — COMPLETE (A1 + A2 + B1 + B3)

All four features of the ROADMAP's recommended first milestone landed this
session, each gate-verified (typecheck 0 / lint 0 errors at the 30-warning
ratchet / tests all green / build OK) and browser-verified via Playwright.
A2/B1/B3 were implemented by delegated Opus subagents and cross-checked by
the orchestrator. Test suite grew 35 → 52 tests (9 files).

### Feature A4: Hillshading & contours — COMPLETE (post-milestone)

- `utils/shading.ts`: `computeShadeMap()` — per-cell Lambert relief factor
  from the tangential height gradient, decoupled from the radial baseline so
  flat land/water sit at exactly 1.0 (no day/night terminator; overlays any
  view). Clamp 0.6–1.15, fixed NW light. `computeContourSegments()` — shared
  Voronoi edges between land cells in different 0.1-elevation bands.
  `drawContourPaths()` shared by Map2D + export.
- Two toggles (default off) through globe (refill-pass color multiply +
  ContourLines segments), both Map2D paths, and PNG export (mirrors
  on-screen toggles). Off = byte-identical rendering. 59 tests green.
- Delegated agent was cut off mid-implementation by the subagent spend
  limit; orchestrator completed Map2D/export/wiring inline.

### Feature A2: Offline namebases — COMPLETE

- `utils/namegen.ts`: order-2 char-level Markov generator, 4 embedded styles
  (fantasy/norse/latin/desert), deterministic from the caller's rng stream.
- Factions/provinces/towns named via dedicated RNG side-streams
  (`civSeed + '_facnames'` / `'_provnames'`) — existing seeds keep
  byte-identical terrain/civ geometry. Gemini is now an optional enhancer.
- `nameStyle` param + Civ-tab select; old saved configs default to 'fantasy'.
- Param-liveness proves nameStyle changes names but not geometry.

### Feature B1: Lakes — COMPLETE

- `generateRivers` returns `{ rivers, lakes }`; contiguous flooded land cells
  (Priority-Flood `waterLevel` above terrain) become `LakeData` entities with
  surface level, outflow pour-point, endorheic + salt flags.
- Salt classification is CLIMATE-driven (mean temp >18°C, moisture <0.3) —
  structural endorheic basins are near-impossible post-Priority-Flood (it
  always finds a spill), so "endorheic → salt" literally would be dead code.
- New LAKE/SALT_LAKE biomes render as water in every view (colors.ts);
  legend auto-picks them up; paint palette excludes them; refreshBiomes
  preserves them. Rivers cut at lake shores and never start in lakes.
- Civs treat lakes as water: no capitals/towns/population, water-cost
  crossing, lake adjacency counts as coast. Heights never mutated.
- Default seed: 29 lakes / 142 cells. Test seed yields 0 lakes so all
  pre-existing signatures stay byte-identical.

### Feature B3: Named geographic features — COMPLETE

- `utils/features.ts`: `detectFeatures(world)` BFS-clusters ranges, deserts,
  forests, oceans/seas, islands; lakes reuse B1 entities 1:1. Read-only, O(n).
- Names via A2 generator on `params.seed + '_geonames'` (TERRAIN seed, not
  civSeed) with kind templates ("X Mountains", "Sea of X", "X Flats" for salt
  lakes). Re-rolling civs never renames terrain — test-enforced.
- Label integration: 7 new LabelKinds through the A1 pipeline; water labels
  italic + blued fill; priorities interleave with civ labels in the shared
  declutter; single `labelVisibility.geography` toggle (default on).
- Inspector shows "Part of: <feature>". Default seed: 5 ranges, 10 deserts,
  7 forests, 1 ocean, 4 islands, 29 lakes.

### Feature A1: Map Labels & Typography — COMPLETE

Multi-tier label system for factions, capitals, provinces, and towns across
3D globe, 2D Mercator, 2D Dymaxion, and PNG export. Committed in 4 chunks
(label engine → 2D/state/export wiring → 3D sprites → docs).

- `utils/labels.ts` (NEW): `MapLabel` model, `collectLabels()` with O(cells)
  bucketing and land-biased centroids, `drawMapLabels()` with greedy priority
  declutter, zoom LOD, and an optional `fontScale` for hi-res exports.
- `utils/geo.ts`: Dymaxion projection promoted from Map2D —
  `projectToDymaxionNet()` returns NET-space coords so each raster pipeline
  applies its own net→canvas fit; `projectDymaxionPoint()` wraps it with the
  screen fit. Export uses its own pad-12/Blender-UV mapping (labels align).
- `types.ts`: `LabelVisibility` + `DEFAULT_LABEL_VISIBILITY` (towns/provinces
  default off). Replaces `showFactionOverlay` everywhere.
- `components/Controls.tsx`: "Map Overlays" checkbox group (borders, faction/
  capital/province/town names); export passes live `labelVisibility` (WYSIWYG).
- `components/Map2D.tsx`: `drawMapLabels()` on Mercator + Dymaxion paths;
  label LOD reads settled zoom via `scaleRef` (no per-wheel-tick cell redraws).
- `components/WorldViewer.tsx`: `PointLabels` canvas-texture sprites at
  r=1.08–1.10 (above 1.05 max terrain + marker pins), camera-distance LOD,
  back-of-globe culling via **world-space** sprite positions (labels spin
  inside the globe group — data-space culling was a real bug, fixed).
- `utils/export.ts`: labels drawn post-raster in `exportMap()` (mirror-
  corrected x, orthographic back-hemisphere skip) and `exportDymaxionRaster()`.

**Key decision — no SDF text in 3D**: drei `<Text>` (troika) spawns a blob
worker and fetches font data from cdn.jsdelivr.net; both violate the strict
CSP (`script-src 'self'`, pinned `connect-src`). Canvas-texture sprites (the
CurvedFactionLabel recipe) keep the CSP untouched and work offline. Do not
reintroduce troika/drei-Text without revisiting the CSP.

**Environment fix**: local `node_modules` was stale (declared `@types/d3`,
`@types/react-dom`, `vitest` missing) — `npm install` fixed it; typecheck is
0 errors, tests 35/35, lint 0 errors/30 warnings (at ratchet), build OK.

**Browser-verified via Playwright**: curved faction labels + capital sprites
on the globe (toggles live), labels on Mercator + Dymaxion, 4K equirect
export PNG contains correctly-placed labels.

---

## Work Completed In This Session (Session 2 - 2026-07-17)

### Audit Hardening Batch (remaining AUDIT.md items)

- **CI added**: `.github/workflows/ci.yml` runs lint + typecheck + test + build
  on pushes to main and all PRs.
- **TypeScript strict mode enabled** (`"strict": true`); zero errors. Added
  `@types/d3`, `@types/react`, `@types/react-dom` and a local `vendor.d.ts`
  shim for `d3-geo-projection` / `d3-geo-voronoi` (no registry types).
  `WorldData.geoJson` is now a typed `GeoJsonCollection` (d3-compatible).
- **Vitest suite added** (`tests/`, 35 tests, ~2 s): RNG determinism, biome
  classification table, generation determinism/structure/abort/progress,
  param liveness (fails if any tunable param stops affecting output), paint
  utils, and import validation. `npm test`.
- **Engine fixes surfaced by the tests**:
  - `civSizeVariance` reworked from expansion budgets (which never bound on
    small maps) to per-faction competitive movement-cost scaling — effective
    at any map resolution.
  - Negative cell populations fixed (high-elevation suitability went negative
    and silently deflated faction totals).
  - `capitalSpacing` threshold was resolution-dependent and almost never
    fired; now a scale-independent squared-chord minimum
    (`spacing^2 * 4 / numFactions`).
- **Tailwind moved from CDN to the build pipeline** (tailwindcss v3 + PostCSS,
  `index.css`, purged ~23 kB output). CSP no longer allows any external script
  host or `unsafe-eval`.
- **`npm audit fix` applied**: 0 vulnerabilities (was 1 critical / 5 high).
- **Lore errors surfaced properly**: `generateWorldLore` throws instead of
  returning sentinel "Error World" lore, and validates the model's JSON field
  by field before mutating civData. `@google/genai` and `GLTFExporter` are now
  dynamically imported (smaller main bundle).
- **Import hardening**: `validateCivData` shape-checks imported civData
  (malformed metadata degrades to terrain-only load); points input capped at
  200k (UI + validation) to match the slider and avoid main-thread freezes.
- **`detailLevel` implemented** as the FBM octave count with a "Detail
  Octaves" Geo slider (default 3 = historical hardcoded value, so default
  worlds are unchanged).
- **Misc**: progress bar reaches exactly 100% (7 ticks / TOTAL_STAGES 7);
  Map2D caption reflects Mercator vs Dymaxion; `ExportResolution` type matches
  the UI (2K/4K/8K); shared Dymaxion geo helpers consolidated into
  `utils/geo.ts`; unused imports removed; lint script has a
  `--max-warnings 30` ratchet (remaining warnings are documented-intentional
  R3F `any`s and hook-dep patterns).

Validation this session: `npm run lint`, `npm run typecheck`, `npm test`
(35/35), `npm run build`, plus a Playwright pass with the bundled Tailwind:
fully styled UI, globe + Dymaxion rendering, new Geo sliders present,
Dymaxion caption correct, no console errors.

---

## Previous Session (Session 1 - 2026-07-17)



### New Features (from AUDIT.md open questions)

- **`civSizeVariance` is now implemented**: each faction draws a per-faction
  expansion budget from `civRng` (base 200, spread `1 ± variance`, clamped to
  0.25x–2x). The Dijkstra frontier stops once a faction's cost exceeds its
  budget. The "Country Size Variance" slider in the Civ tab now does something.
- **`population` view mode implemented**: log-scaled heat gradient on land
  (dark blue → green → bright yellow), uninhabited land dark grey, ocean navy.
  New "Population" View Layer button in the Sys tab.
- **`province` view mode implemented**: faction base color with amplified
  per-province shade variation (variant strength 1.8). New "Provinces" button.
- **`Pangea` land style** now has a real engine branch in `worldGen.ts`
  (`landChance 0.6, landLevel 0.25, oceanLevel -0.45`), on top of the mask the
  preset already set.
- **Cell Jitter slider added to the Geo tab** (was documented but missing).
- **Build-time Gemini key fixed**: `services/gemini.ts` now reads
  `process.env.GEMINI_API_KEY`, matching the `vite.config.ts` define and the
  README's `.env.local` instructions (previously read `API_KEY`, which was
  always undefined). Verified via sentinel-key build.

### 3D Rendering Optimizations (`WorldViewer.tsx`)

- **World mesh geometry is now allocated once per world structure and refilled
  in place** on paint strokes / view changes (`useLayoutEffect`), instead of
  rebuilding + leaking a new ~45k-triangle `BufferGeometry` on every stroke
  event. Attributes use `DynamicDrawUsage`; bounding sphere is fixed (r = 1.1).
- **Dropped `computeVertexNormals` and the normal attribute** from the world
  mesh, selection overlay, and curved labels — the unlit basic material ignores
  normals and the standard material's `flatShading` derives them in-shader.
- **Every `useMemo` geometry now has a disposal effect** (world mesh, rivers,
  faction borders, brush ring, selection overlay, highlight outline, lat/long
  grid, Dymaxion overlay, curved label geometry). GPU memory no longer grows
  across regenerations and paint sessions.
- **Rivers keyed on `world.rivers`** (stable across strokes) so painting never
  re-runs CatmullRom smoothing; **curved faction labels keyed on position
  values** so unchanged centroids don't rebuild patches; the inline cell
  highlight was converted to a memoized `CellHighlightOutline` component.
- **`CityMarkers` no longer allocates Vector3s per instance per update.**
- **Map2D Dymaxion pick buffer keyed on `world.cells`** instead of world
  identity — skips a full-canvas per-pixel reprojection per stroke event.
- **MiniMap redraw debounced (150 ms)** and now passes live faction colors.

### Color Consistency

- New `buildFactionColorMap(civData)` helper in `colors.ts`; `MiniMap`,
  `export.ts` (both render loops), and `exportGLB.ts` now pass live faction
  colors, so user-edited faction colors appear in the minimap, PNG exports,
  and GLB vertex colors (closes the gap invariant 21 warned about).

### Correctness / Type Fixes

- `npx tsc --noEmit` is now clean (was 6 errors): `Controls.tsx` `setParams`
  prop typed as `React.Dispatch<React.SetStateAction<WorldParams>>`; dead
  `dymaxion` switch branch removed from `exportMap` (unreachable after the
  early return); dymaxion projection nodes typed with `children`; dead
  comparison removed in `Inspector.tsx`.
- The `G` keyboard shortcut no longer generates from stale params (and no
  longer reverts recent slider edits): the keydown effect re-registers on
  `params`/lock changes.
- Auto-update now also triggers on `mountainHeight`, `oceanDepth`, and
  `cellJitter` changes (previously only the first change fired, via the
  `landStyle → 'Custom'` side effect).

---

## Files Most Recently Touched

| File | Notes |
|------|-------|
| `utils/labels.ts` | **NEW** — `MapLabel`, `LabelKind`, `collectLabels()`, `drawMapLabels()` with declutter + LOD. |
| `utils/geo.ts` | Promoted Dymaxion helpers from Map2D: `projectDymaxionPoint`, `getDymaxionNetTransform`, `dot3`/`sub3`/`cross3`, `barycentric3D`, `pointInsideSphericalFace`. |
| `types.ts` | Added `LabelVisibility` interface, `DEFAULT_LABEL_VISIBILITY`. |
| `App.tsx` | `showFactionOverlay` → `labelVisibility` state; prop-drills to children. |
| `components/Controls.tsx` | Per-kind Map Overlays toggles replacing single Faction Overlay checkbox. |
| `components/Map2D.tsx` | `drawMapLabels()` replaces `drawFactionLabels()`; `getFactionBorders()` replaces `getFactionOverlayData()`. Removed promoted helpers. |
| `components/WorldViewer.tsx` | `PointLabels` component (drei `Text`+`Billboard`, SDF rendering). `labelVisibility` wired into `CountryLabels` + `FactionBorders`. |
| `utils/export.ts` | Labels drawn in both export paths, honoring visibility toggles. |

---

## Validation

Last successful checks:

- `npm run build`
- `npm run lint` (zero errors; warnings only)
- `npx tsc --noEmit` (zero errors — new since this session)
- Headless-browser pass (Playwright + SwiftShader): generation completes with
  no console errors; Biomes/Provinces/Population layers render on the 3D globe,
  2D Mercator, and minimap; Cell Jitter slider present. (Tailwind CDN is
  blocked in the sandbox, so layout was stubbed — visual styling not covered.)
- `.env.local` sentinel key verified to land in the built bundle.

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

## Feature Roadmap (shortlist)

See `ROADMAP.md` for the full detail. Themes, by leverage:

- **A. Cartographic presentation** — map labels/typography, offline Markov
  namebases (Gemini becomes optional), fantasy map styles, hillshading +
  contours, great-circle ruler + scale bar, geodesic hex grid.
- **B. Physical geography** — lakes as first-class entities (the Priority-
  Flood fill in `generateRivers` already computes them and throws them away),
  river/lake editing, auto-detected + named ranges/seas/deserts/islands.
- **C. Worldbuilding depth** — cultures layer, religions, roads + sea trade
  routes (reuse civ Dijkstra costs), markers/POIs, editor completeness,
  diplomacy later.
- **D. Planet-scale simulation** — seasonal cycle from axial tilt, ocean
  currents feeding climate, ice caps, regional submap generation, planetary
  parameters. These are the sphere-native differentiators vs. flat-map tools.
- **E. Interoperability** — SVG export, GeoJSON export, Azgaar `.map` import
  (stretch).

**"Map identity" milestone (A1, A2, B1, B3) SHIPPED in session 3.** The next
tier per ROADMAP's suggested order: A3 (map styles), A4 (hillshading), A5
(great-circle ruler + scale bar), B2 (river/lake editing), C3 (roads/trade
routes), C4 (markers/POIs), E1/E2 (SVG/GeoJSON export).

---

## Potential Next Tasks

- Manual browser regression pass across 3D, 2D Mercator, and 2D Dymaxion with
  real styling (the sandbox pass could not load Tailwind).
- Remaining items from `AUDIT.md`: CI workflow (lint + typecheck + build),
  Vitest suite over the pure engine, staged TypeScript strict mode,
  `npm audit fix`, Tailwind via build pipeline instead of CDN, code splitting.
- Optional next-level paint optimization: per-cell dirty ranges so strokes
  refill only affected vertices instead of all cells (current in-place refill
  is already allocation-free).
- Tune 3D curved faction label sizing/placement if very long faction names are
  introduced.
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
