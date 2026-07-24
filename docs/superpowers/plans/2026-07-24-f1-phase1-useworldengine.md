# F1 Phase 1 — `useWorldEngine` Extraction Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [x]`) syntax for tracking.

**Goal:** Extract all of `App.tsx`'s state, refs, effects, memos, and handlers into a `useWorldEngine()` hook, leaving `App`'s render output byte-identical, so a future `ShellApp` can consume the same state.

**Architecture:** Pure, behavior-preserving refactor. All logic moves verbatim from `App` into `hooks/useWorldEngine.ts`, which returns one typed `WorldEngine` object. `App` becomes a thin consumer: `const e = useWorldEngine()` + destructure + the **unchanged** return block. No new behavior, no new UI, no flag changes — the classic app is the only thing that runs, and it must run identically.

**Tech Stack:** React 19, TypeScript (strict), Vite, Vitest, Tailwind v3.

## Global Constraints

- **Relative imports only** (the `@/` alias is intentionally unused).
- **Return block is frozen:** the JSX inside `App`'s `return ( … )` (currently lines 654–810) changes by **zero characters** in this phase. Only bindings above it move.
- **Lint ratchet:** `npm run lint` = 0 errors, **exactly 30 warnings**. Add none.
- **Determinism is sacred:** the `tests/` suite (determinism, param-liveness, paint) must stay green; a green suite proves self-consistency, NOT old-vs-new equivalence — the frozen return block + manual regression are the real equivalence checks.
- **This is a refactor, not TDD:** there is no new behavior to test-drive. The "test cycle" is the existing suite + typecheck + a return-block diff + a manual browser regression. Do not invent new unit tests for moved code.
- Do NOT touch `HANDOFF.md`, `CLAUDE.md`, `ROADMAP.md`, or the spec in this phase.

---

## File Structure

- **Create `hooks/useWorldEngine.ts`** — owns every `useState`/`useRef`/`useEffect`/`useMemo`/`useCallback` currently in `App`, plus the module-level `DEFAULT_PARAMS` constant and the three pure helpers (`isValidProvinceId`, `resolvePoliticalProvinceId`, `recalculatePoliticalTotals`). Exports `useWorldEngine()` and the `WorldEngine` return type.
- **Modify `App.tsx`** — delete lines 18–652 (the helpers + component body above the return), replace the component body with a single `useWorldEngine()` call + destructuring, keep the return block (654–810) verbatim. Keep the `DEFAULT_PARAMS`/helpers OUT of `App` (they live in the hook now).

---

### Task 1: Extract everything into `useWorldEngine`

This is one atomic relocation: the code only compiles once *all* interdependent bindings live in the same scope, so it moves as a unit. The steps below stage it safely and verify with the regression gate.

**Files:**
- Create: `hooks/useWorldEngine.ts`
- Modify: `App.tsx:18-652` (remove), `App.tsx:1-16` (imports), `App.tsx:105-106` (consume hook)
- Test (regression only): `tests/` (existing suite, unchanged)

**Interfaces:**
- Produces: `useWorldEngine(): WorldEngine`, where `WorldEngine` is an interface exposing **every binding the return block (654–810) references** — state values + their setters, derived memos, and handlers. Exact member list is mechanical: it is precisely the set of identifiers used inside `return ( … )`.

- [x] **Step 1: Snapshot the return block for the frozen-block check**

Before editing, capture the current return block so you can prove it didn't change:

```bash
sed -n '654,810p' App.tsx > /tmp/app-return-before.txt
wc -l /tmp/app-return-before.txt   # expect 157 lines
```

- [x] **Step 2: Create the hook file with all logic moved verbatim**

Create `hooks/useWorldEngine.ts`. Move, **without editing their bodies**:
1. The module-level `DEFAULT_PARAMS` const and the three pure helpers (`isValidProvinceId`, `resolvePoliticalProvinceId`, `recalculatePoliticalTotals`) from `App.tsx:18-103`.
2. Every `useState`, `useRef`, `useEffect`, `useMemo`, and `useCallback` from the `App` component body (`App.tsx:106-652`).
3. All logic imports they need (from `App.tsx:8-15`: `types`, `worldGen`, `routes`, `civEdit`, `paintUtils`, `gemini`, `measure`). Do NOT move component imports (`Controls`, `WorldViewer`, `Map2D`, `MiniMap`, `Inspector`, `Legend`, `EditToolbar`, `lucide-react`) — those stay in `App`.

Structure:

```ts
import { useState, useEffect, useCallback, useRef, useMemo } from 'react';
// ...the logic imports moved from App (types, worldGen, routes, civEdit, paintUtils, gemini, measure)...

const DEFAULT_PARAMS: WorldParams = { /* verbatim from App */ };

// ...the three pure helpers, verbatim from App...

export interface WorldEngine {
  // one member per identifier the return block uses.
  // e.g.: params: WorldParams; setParams: React.Dispatch<React.SetStateAction<WorldParams>>;
  //       world: WorldData | null; handleGenerate: (p?: WorldParams) => Promise<void>; ...
  // Fill from the destructure in Step 3 — the two lists must match exactly.
}

export function useWorldEngine(): WorldEngine {
  // ...ALL state, refs, effects, memos, handlers, verbatim from App:106-652...
  return {
    // every binding the return block references
  };
}
```

Do not change effect dependency arrays, ref initial values, or handler bodies. This is a move, not a rewrite.

- [x] **Step 3: Rewire `App.tsx` to consume the hook, return block untouched**

`App.tsx` collapses to:

```tsx
import React from 'react';
import Controls from './components/Controls';
import WorldViewer from './components/WorldViewer';
import Map2D from './components/Map2D';
import MiniMap from './components/MiniMap';
import Inspector from './components/Inspector';
import { BiomeLegend } from './components/Legend';
import EditToolbar from './components/EditToolbar';
import { Menu, X } from 'lucide-react';
import { useWorldEngine } from './hooks/useWorldEngine';

const App: React.FC = () => {
  const {
    // destructure EVERY name used in the return block below
    params, setParams, world, viewMode, setViewMode, /* …all of them… */
  } = useWorldEngine();

  return (
    // lines 654-810 verbatim — DO NOT EDIT
  );
};

export default App;
```

Keep only the type imports the return block itself needs (e.g. any `types.ts` types referenced in JSX). If the return block references no bare types, drop the `types` import from `App`.

- [x] **Step 4: Typecheck (catches every un-threaded binding)**

Run: `npm run typecheck`
Expected: **0 errors.** Any "Cannot find name 'X'" means `X` is used in the return block but missing from the destructure/`WorldEngine` return — add it to both. This is the compiler doing the fidelity proof for you.

- [x] **Step 5: Prove the return block is unchanged (frozen-block invariant)**

```bash
diff <(sed -n '/^  return (/,/^};$/p' App.tsx) /tmp/app-return-before.txt
```

Expected: **no diff** (the return block content is identical; only its line numbers moved). If `diff` reports changes, you edited JSX — revert those edits.

- [x] **Step 6: Lint + full test suite (regression gate)**

Run: `npm run lint`
Expected: 0 errors, **exactly 30 warnings**.

Run: `npm test`
Expected: all green (138 tests), same count as before.

- [x] **Step 7: Manual browser regression — the real equivalence check**

The dev server picks up new files only on restart (known Tailwind gotcha), so restart it first if it's running.

Run: `npm run dev` (port 3000), open `http://localhost:3000` (NO `?shell`).
Verify against the pre-refactor app:
- Generate the default world (button + `G` shortcut) — completes, same log output, progress reaches 100%.
- Orbit the globe; switch view layers (Biomes/Political/Population); toggle Mercator/Dymaxion.
- Click a cell → Inspector populates; enter an edit mode → paint a stroke → undo.
- No new console errors.

Focus attention on the ref/effect surface the spec flagged: generation abort (regenerate mid-run), paint-stroke undo, and the API-key effect — these are where a mis-moved dep would surface.

- [x] **Step 8: Commit**

```bash
git add hooks/useWorldEngine.ts App.tsx
git commit -m "Extract App state into useWorldEngine hook

Behavior-preserving refactor: all state, refs, effects, memos, and
handlers move verbatim into hooks/useWorldEngine.ts; App consumes the
hook and keeps its return block byte-identical. Prepares a future
ShellApp to share the same engine. No UI or behavior change.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Self-Review

**1. Spec coverage (§3.1):** Task 1 extracts the hook (✓), App consumes it (✓), return block frozen via Step 5 (✓), regression via Steps 4/6/7 (✓), fidelity review pointed at effects/refs via Step 7 (✓). The spec's §3.2–§3.6 (decomposition, positioning, shell, globe toggle, entry) are explicitly **out of scope for Phase 1** and covered by later plans.

**2. Placeholder scan:** The `WorldEngine` member list and the App destructure list are described as "every identifier the return block uses" rather than enumerated — this is intentional and exact (the set is mechanically determined by the frozen return block), not a placeholder. Steps 4–5 make the compiler enforce completeness. No "TODO"/"handle edge cases"/vague steps.

**3. Type consistency:** `useWorldEngine` / `WorldEngine` names are consistent across the hook definition (Step 2) and App's consumption (Step 3). The return object and the destructure are cross-checked by typecheck (Step 4).

---

## Note on later phases

Phase 1b (positioning extraction), Phase 2 (ShellApp v1 + `?globe=0`), Phase 3 (Controls decomposition), and Phase 4 (polish) each get their own plan once Phase 1 lands and the `WorldEngine` shape is real. See `docs/superpowers/specs/2026-07-24-f1-ui-redesign-design.md` §4.
