# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Quick Orientation

RealmGenesis 3D is a browser-only procedural fantasy world engine — spherical Voronoi cells, tectonic plates, erosion, climate, biomes, civilizations, and AI lore. React 19 + Three.js/R3F for 3D, Canvas2D + d3 for 2D projections, Google Gemini for lore. No backend.

Read `HANDOFF.md` first for recent session state, then `docs/` (start at `docs/README.md`) for settled technical context, then `AGENTS.md` for style rules. The old monolithic `ARCHITECTURE.md` is archived at `docs/archive/` and known to drift — don't trust it.

When updating `HANDOFF.md`, add a date and session number.

## Commands

```bash
npm run dev        # Vite dev server on port 3000
npm run build      # Production build → dist/
npm run lint       # ESLint (zero errors, --max-warnings 30 ratchet)
npm run typecheck  # tsc --noEmit (strict mode)
npm test           # Vitest suite (tests/)
```

Run a single test file: `npx vitest run tests/rng.test.ts`

All four gates run in CI (`.github/workflows/ci.yml`) and must pass before merge.

The param-liveness test (`tests/paramLiveness.test.ts`) fails if any `WorldParams` key stops influencing generated output — extend it when adding new params.

## Architecture at a Glance

**All state lives in `hooks/useWorldEngine.ts`** and is prop-drilled by the shell (`components/shell/ShellApp.tsx`, the default route). No Context, Redux, or Zustand. Trace any value up to `useWorldEngine`. (`App.tsx` is the legacy `?shell=classic` route.)

**Generation pipeline** (`utils/worldGen.ts` → `generateWorld()`): runs in a **Web Worker** (`workers/worldGen.worker.ts` via `utils/worldGenClient.ts`), 7 progress stages — Fibonacci+Voronoi → V3 tectonic macro-sim → project to display → erosion → Stage-9b remap → climate/biomes → rivers/lakes → civs. Cancellable via `worker.terminate()`. See `docs/generation-pipeline.md` and `docs/tectonics-v3.md`.

**Rendering**: `WorldViewer.tsx` (3D globe via R3F), `Map2D.tsx` (Mercator/Dymaxion via Canvas2D). Cell colors come from `utils/colors.ts` → `getCellColor(cell, viewMode, seaLevel, factionColors?)`.

**Civilization**: `recalculateCivs()` (Dijkstra expansion), `recalculateProvinces()` (subdivision + towns). Can run independently of terrain generation.

**AI lore**: `services/gemini.ts` → `generateWorldLore()` **mutates `world.civData` in-place**. Caller must `setWorld({ ...world })` to trigger re-render.

## Key Invariants

- **Relative imports only** — `@/` alias is configured but intentionally unused.
- **`seaLevel` must be passed to `getCellColor`** as the third argument (from `world.params.seaLevel`), not hardcoded.
- **`factionColors` map required for political rendering** — any render path calling `getCellColor` for political/province mode must pass a faction-color map (build via `buildFactionColorMap(civData)` from `colors.ts`).
- **R3F element names are strings** (e.g., `"bufferGeometry"`) — intentional pattern to bypass TSX types. `@typescript-eslint/no-explicit-any` is warn, not error.
- **Gemini API key is ephemeral** — never persisted to storage. Set via `setRuntimeApiKey()` or build-time `GEMINI_API_KEY` env var.
- **WorldMesh geometry is reused across paint strokes** — `world.cells` identity is the structural key. Paint strokes mutate cells in place and shallow-copy `WorldData`; never replace the `cells` array outside full regeneration.
- **Every `useMemo` geometry in `WorldViewer` has a matching disposal effect** — follow this pattern when adding scene elements.
- **No `plateInfluence` clamp exists** — the param was renamed `tectonicStrength` (range 0–2, unclamped); the V2 `[0.1, 1.0]` clamp was deleted with the V2 engine. The `export.ts` import validator now bounds `tectonicStrength: [0, 2.0]` directly (the old dead `plateInfluence` key was renamed to it).
- **`mountainHeight`/`oceanDepth`/`seafloorDepth` remap is after normalization, before climate** (Stage 9b) — inserting normalization steps requires adjusting remap placement. `oceanDepth` is a contrast power-curve; `seafloorDepth` is a linear mean-depth datum.
- **Batch `setParams` calls** — use functional updater `setParams(prev => ({ ...prev, ...changes }))` to avoid stale-closure overwrites.
- **Edit undo uses a shared Map reference** — `currentStrokeSnapshot` ref is pushed to `undoStack` at stroke start and mutated in place; never replace it mid-stroke.
- **Dymaxion pick buffer must mirror visible rasterization** — same pipeline, same rotation, same sizing.

## Code Style

- 2-space indent, semicolons, single quotes, trailing commas
- Import order: React → external → local (relative paths), blank line between groups
- `interface` for objects/props, `type` for unions
- Components: `React.FC<Props>`, functional only
- Naming: PascalCase components/types, camelCase functions/vars, SCREAMING_SNAKE enum values
- `handle*` for event handlers, `is*/show*` for booleans
- Tailwind CSS utility classes only — no CSS modules or styled-components
- Dark theme: `bg-gray-950`, `text-gray-200`, `border-gray-800`

## Testing

Vitest suite in `tests/` covers the pure engine (RNG, biomes, generation determinism, param liveness, paint utils, import validation). Rendering behavior is verified manually in the browser.

## Deployment

Netlify static SPA. `public/_redirects` handles SPA routing. Optional `GEMINI_API_KEY` env var at build time or BYOK at runtime.

## Git workflow

Commit locally in small chunks, scoped to one topic or one unit of work each (e.g. "fix screen centering," "add README" — not one giant catch-all commit). Don't batch unrelated changes together just because they landed in the same session. Only push when the user explicitly asks.

## HANDOFF discipline

Update `HANDOFF.md` **as you work**, not only when asked. A long session can be
compacted or die, and unwritten findings die with it — but rationale is the part
that really evaporates. Nobody reconstructs six weeks later *why* a rule was
broken; written at the moment of the decision it costs one paragraph.

What goes where, by type:

- **Decisions + their rationale** — write immediately. Highest value, most
  perishable. Include the alternative that was rejected and why.
- **Findings** — write immediately too, but **at the confidence level you
  actually have.** Say "n=1, unconfirmed" out loud when that's what it is.
  Reserve "verified" for claims you would defend under challenge, naming the
  evidence.
- **Progress narration** — never. That is what `git log` is for.
- **The session entry** — at the end, or when context is about to compact.

The test: *does the next session need this to avoid repeating a mistake or
re-litigating a decision?* If no, it belongs in a commit message, not HANDOFF.

**Record refuted hypotheses, don't delete them.** A wrong idea that looked right
is useful — it stops the next session re-deriving it. Correct in place, state
what refuted it, and keep the original reasoning visible.

This exists because of a real failure (session 2026-07-24): three claims were
written to HANDOFF as settled fact from single observations — "extended thinking
was off," "the agent registry is frozen at session start," "worktrees cut from
session-start HEAD" — and all three were refuted within hours, each needing a
correction commit. The frequency was not the problem. Asserting n=1 findings as
conclusions was. Hedging them at write time would have made every later discovery
an update instead of a contradiction.


## Web Search Workflow

- For fact checking and fan out for research, Claude must use Google Antigravity CLI by way of invoking `agy` in an internal shell and/or using `/antigravity:research`. Claude, or a Claude subagent, is BLOCKED from using its internal WebSearch tool UNLESS instructed to so without needing `agy`, or to run a parallel adversarial search.

## Dev server etiquette

Claude must not kill a preexisting running dev server unless
Claude itself started it (e.g. to drive Playwright screenshots/verification).