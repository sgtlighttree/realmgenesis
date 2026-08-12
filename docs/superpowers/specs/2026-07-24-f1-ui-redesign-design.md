# F1 — UI Redesign: Implementation Design

Status: DRAFT for review · Branch: `redesign` · Date: 2026-07-24 (Session 5)

## 1. Context & goal

The UI is *ungoverned*, not ugly: six self-positioning floaters ring the globe,
the left rail mixes generation params / view toggles / system log in one column,
and contextual tools (the edit bar) are permanent chrome. Worst on mobile.

**F1 delivers:** the app's REAL functionality re-housed into two responsive
shells built on a clean content/shell split, plus a consistency cleanup.
Validated already by the `?shell=1` prototype (stub panels + placeholder globe).

**Explicitly out of scope:**
- **F1b** — the visual *identity* (type pairing, palette, motion language). A
  later dedicated `/impeccable` pass on the finished skeleton.
- **F2 / F3** — the render modes are STUBBED. The smooth/relief globe toggle and
  vector-2D are not built here; the placeholder-globe path (`?globe=0`) is the
  only new render surface.

## 2. Settled decisions (from brainstorm + prototype)

- **Buckets:** Make (regenerates) / View (instant) / Do (contextual) / Read
  (reference). The organizing spine.
- **Layout:** `WideShell` = "A · Tidy" (Make rail, View strip, Read stack right,
  contextual Do bar). `NarrowShell` = "C · Studio" (globe hero, bottom tab bar,
  one sheet at a time). Same content, chosen per breakpoint (768px).
- **Content/shell decoupling (advisor-confirmed, no-Context invariant intact):**
  the owner composes each panel WITH its props into a finished element; the shell
  only POSITIONS named slots and never sees panel props. Ephemeral presentation
  state (which mobile tab is open) lives in the shell's own `useState`.
- **Edit:** single "Edit" toggle, default off, summons the Do bar; Esc exits.
  Narrow: Do is one of the four tabs.
- **Chrome (from the impeccable pass, already in prototype):** one `PANEL`
  constant, single blue accent for state/selection only, solid surfaces,
  `rg-rise` state motion w/ reduced-motion fallback, visible focus rings.

## 3. Architecture

### 3.1 State — `useWorldEngine()` hook

Extract App's ~30 state atoms, refs, and ~12 handlers verbatim into
`hooks/useWorldEngine.ts`, returning one typed object. **Move, don't refactor** —
same declarations, same order, same effect deps and closures. Logic changes are a
separate later commit, never mixed into the extraction.

- `App.tsx` becomes a thin consumer that renders TODAY'S layout from the hook.
- New `ShellApp.tsx` consumes the SAME hook and mounts real components into shell
  slots.
- **Regression gate — a hard textual invariant:** App's return block (lines
  654–810) changes by **zero characters** during the extraction; only the
  declarations above it move into the hook. The git diff of the return block IS
  the proof.
- **Where the real risk is:** TypeScript already guards the structural half — if
  a member isn't threaded through the hook, the render won't compile. The narrow
  human-review surface is **effect dependency arrays and ref timing/identity**:
  `abortControllerRef`, `isPainting`, `currentStrokeSnapshot`,
  `currentPoliticalProvince`, and the `apiKey` effect. Review those; don't re-read
  all 650 lines.

**Classic layout lifespan:** transitional. Once the shell reaches feature parity
it is retired — do not invest in maintaining two layouts through Phases 3–4.

### 3.2 Component decomposition

- **`Controls` (1571 lines) is a Make+View+System monolith.** Split its contents
  by bucket: gen params + Generate + seed/resolution → **Make**; render mode +
  layer toggles + label visibility → **View**; system console → **Read** (or a
  collapsible strip). The GEO/CLIM/CIV tabs are Make; EXP is Make (export lives
  with the world you made). The SYS tab is the mongrel — its contents get
  re-sorted across Make/View/Read.
- **Floaters** (`Inspector`, `BiomeLegend`, `MiniMap`, `EditToolbar`) currently
  self-position with `absolute …`.

### 3.3 The shared-positioning problem (needs a decision — advisor input)

These floaters are used by BOTH the classic layout and the shell. If I strip
their `absolute` positioning so they can sit in shell slots, the CLASSIC layout
(which depends on it) breaks. Two options:

- **(A) Positioning moves to the caller.** Components render bare (fill their
  container); the classic layout wraps each in its existing `absolute` div, the
  shell drops each into a slot. Clean end state, both layouts explicit; touches
  the classic layout.
- **(B) A `bare`/`unstyled` prop** suppresses self-positioning in shell mode
  only. Less churn to the classic layout, but a component with two positioning
  personalities is a smell.

**Decision: (A).** Components render bare; positioning is the caller's job. (B)'s
two-personality component is exactly the smell to avoid. This is a **second
classic-touching change** (the classic layout must add wrapper divs), so it gets
its own phase and its own regression pass — see Phase 1b.

### 3.4 Slot contract & Edit wiring

`ShellApp`: `const e = useWorldEngine(); const panels = { make:…, view:…,
doTools:<EditToolbar …/>, read:[…], canvas: globe }; return effective==='wide' ?
<WideShell {...panels} …/> : <NarrowShell {...panels} …/>`.

- Edit toggle maps to `e.editMode`: "on" sets a sensible default mode (e.g.
  `terrain-raise`), "off"/Esc sets `'off'`. `editing = e.editMode !== 'off'`.

### 3.5 Globe toggle — `?globe=0`

Read once at entry; the canvas slot renders `WorldViewer` normally or
`PlaceholderGlobe` when `globe=0`. Escape hatch for fast UI iteration without the
Three.js cost.

### 3.6 Entry (`index.tsx`)

- `?shell=1` → `ShellApp` (real, playable).
- `?shell=stub` → `DesignShell` (keep the stub harness for isolated layout work).
- otherwise → `App` (classic).

## 4. Phasing

1. **Extract `useWorldEngine`**; App consumes it. Return block frozen (§3.1). No
   UI change. Classic-touching — its own regression pass.
1b. **Positioning extraction (§3.3-A):** `Inspector` / `Legend` / `MiniMap` /
   `EditToolbar` go bare; classic App wraps each in its existing `absolute` div.
   Still no `?shell` involvement. **Classic-touching — its own visual diff.**
2. **`ShellApp` v1 — playable:** real canvas (+`?globe=0`), whole `Controls` in
   the Make slot for now, real bare `Inspector` + `EditToolbar` in slots.
   Generate / orbit / inspect / edit all work in the shell. Behind `?shell`.
3. **Decompose `Controls`** into Make / View; wire the View slot; re-sort SYS.
4. **Polish** spacing/states/rhythm. (Then F1b: visual identity, separate.)

**Two of these touch the classic app (1 and 1b), each with its own regression
pass** — the earlier "only Phase 1 is risky, 2+ are shell-only" framing was
wrong; 1b smuggled a classic layout change into a "shell-only" phase. Phases 2–4
are behind `?shell`.

## 5. Testing

- Determinism + param-liveness suites stay green throughout.
- Phase 1: classic UI unchanged (regression) — the extraction's only job.
- Shell: browser-verify generate, orbit, inspect, paint/edit, both breakpoints,
  and `?globe=0`, against a FRESH build (dev-server Tailwind won't hot-scan new
  files — see HANDOFF Session 5 gotcha).

## 6. Risks

- **Extraction fidelity** — TypeScript covers the structural half (an un-threaded
  member won't compile). The genuine risk is narrow: effect dependency arrays and
  ref timing/identity (§3.1). Mitigation: verbatim move, frozen return block, eyes
  on the effects/refs only.
- **`Controls` decomposition** — 1571 lines; re-sorting the SYS tab is fiddly.
- **Shared-component positioning (§3.3-A)** — a classic-touching change; the
  classic layout must keep working (Phase 1b regression pass).
- Determinism-green ≠ equivalence: it only proves self-consistency, not that the
  refactor preserved old output. The manual regression diff is the real check.
