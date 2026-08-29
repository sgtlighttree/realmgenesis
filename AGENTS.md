# AGENTS.md -- RealmGenesis 3D

## Session Handoff

Read [HANDOFF.md](./HANDOFF.md) first for a summary of what was built recently, known issues, and what to tackle next. It is updated at the end of each working session and is the fastest way to orient to the current state of the project.

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

## Architecture Overview

Before making significant changes, read [`docs/`](./docs/README.md) for a complete overview of the codebase — the generation pipeline, data model, rendering architecture, state management, and key invariants. Start at `docs/README.md`. (The old `ARCHITECTURE.md` monolith is archived at `docs/archive/` and known to drift.)

Key entry points:

- **All app state**: `hooks/useWorldEngine.ts` — the state-owning hook, consumed by `components/shell/ShellApp.tsx` (the default route) and prop-drilled. (`App.tsx` is the legacy `?shell=classic` route.)
- **Data types**: `types.ts` — `Cell`, `WorldData`, `WorldParams`, `BiomeType`
- **Generation logic**: `utils/worldGen.ts` — `generateWorld()` (runs in a Web Worker; V3 pipeline)
- **Color mapping**: `utils/colors.ts` — `getCellColor(cell, viewMode, seaLevel, factionColors?, cultureColors?, religionColors?)`
- **Map painting**: `utils/paintUtils.ts` — brush BFS, stroke functions, undo snapshots
- **Edit toolbar**: `components/EditToolbar.tsx` — paint/edit mode HUD
- **3D rendering**: `components/WorldViewer.tsx`
- **2D rendering**: `components/Map2D.tsx`
- **AI lore**: `services/gemini.ts` — `generateWorldLore()` mutates `world.civData` in-place
- **Save/export**: `utils/export.ts` — `exportMap()`, `saveMapToBrowser()`, `saveMapConfig()`

---

## Commands

```bash
npm run dev        # Start Vite dev server on port 3000
npm run build      # Production build → dist/
npm run preview    # Preview production build
npm run lint       # ESLint (flat config, zero errors, --max-warnings ratchet)
npm run typecheck  # tsc --noEmit (strict mode)
npm test           # Vitest suite over the pure engine (tests/) — see the run
                   # policy in CLAUDE.md; prefer naming the files you touched
                   # (npx vitest run tests/foo.test.ts) over the whole suite
```

**No formatter is configured.** Rendering behavior is verified manually in the browser; the pure engine (`utils/`, `services/`) is covered by the Vitest suite in `tests/`. All four gates run in CI (`.github/workflows/ci.yml`) and must pass:
1. `npm run lint` — zero errors, warning count at or below the `--max-warnings` ratchet
2. `npm run typecheck` — zero errors under `"strict": true`
3. `npm test` — all Vitest tests pass (note: the param-liveness test fails if a `WorldParams` key stops affecting output)
4. `npm run build` — succeeds with no errors

## Git Workflow

Commit locally in small chunks, scoped to one topic or one unit of work each.
Don't batch unrelated changes. Commit messages should follow the 50/72 rule:
subject line at or under 50 characters, body wrapped at 72. Prefer concise,
imperative subjects.

## Code Style

### Imports
- Order: React → external libraries → local modules (relative paths)
- Use relative imports (`../types`, `./utils/colors`), NOT the `@/` path alias (configured in tsconfig but unused in practice)
- Group imports with a blank line between external and local

### Formatting
- 2-space indentation
- Semicolons on all statements
- Single quotes for strings, backticks for template literals
- Trailing commas in multi-line objects/arrays
- Max line length: ~120 chars (soft limit)

### TypeScript
- Strict mode: `"strict": true` (plus `skipLibCheck: true`, `allowJs: true`, `noEmit: true`)
- Use `interface` for object types and component props, `type` for unions
- Prefer explicit return types on exported functions
- Use `as any` casts sparingly (only when required by library typing gaps, e.g., R3F element names in `WorldViewer.tsx`)
- `@typescript-eslint/no-explicit-any` is a warning, not an error — acceptable for R3F event handlers and d3 projections

### React Components
- Functional components with `React.FC<Props>` type annotation
- Props defined as `interface ComponentProps { ... }`
- Use `useCallback` for event handlers passed as props, `useMemo` for expensive computations
- `useState` for local state, all app-level state lives in `hooks/useWorldEngine.ts`
- No class components

### Naming Conventions
- Components: PascalCase (`WorldViewer`, `Map2D`)
- Functions/variables: camelCase (`generateWorld`, `handleGenerate`)
- Types/Interfaces: PascalCase (`WorldData`, `Cell`, `ControlsProps`)
- Enums: PascalCase (`BiomeType`, with SCREAMING_SNAKE_CASE values)
- Type aliases for string unions: PascalCase (`ViewMode`, `DisplayMode`, `LandStyle`)
- Event handlers: `handle*` prefix (`handleGenerate`, `handleCancel`)
- Boolean state: `is*`, `show*` (`isGenerating`, `showRivers`)

### Error Handling
- Use `try/catch` with `console.error` for logging and user-facing messages via the `addLog` callback
- Never throw unhandled errors; always provide fallback UI state
- Use `AbortController` for cancellable async operations (generation pipeline)
- Validate imported JSON configs with `validateWorldParams()` before use

### Styling
- Tailwind CSS 3 via the build pipeline (`index.css`, `tailwind.config.js`) — use utility classes exclusively, no CSS modules or styled-components; class names must appear complete in source for the content scan
- Dark theme: `bg-gray-950`, `text-gray-200`, `border-gray-800`
- Responsive: mobile-first with `md:` breakpoints
- Overlays: `backdrop-blur-md`, `bg-black/50`, `border-white/10`

### Architecture Patterns
- Single source of truth: `hooks/useWorldEngine.ts` holds all state, the shell passes it down via props
- Utils are pure functions (no side effects except logging callbacks)
- Services (e.g., `gemini.ts`) wrap external APIs with minimal abstraction
- Components are presentational — no data fetching or generation logic

## Key Invariants

- **Relative imports only** — `@/` alias is configured but intentionally unused.
- **`seaLevel` must be passed to `getCellColor`** as the third argument (from `world.params.seaLevel`), not hardcoded.
- **`factionColors` map required for political rendering** — any render path calling `getCellColor` for political/province mode must pass a faction-color map (build via `buildFactionColorMap(civData)` from `colors.ts`).
- **R3F element names are strings** (e.g., `"bufferGeometry"`) — intentional pattern to bypass TSX types. `@typescript-eslint/no-explicit-any` is warn, not error.
- **Gemini API key is ephemeral** — never persisted to storage. Set via `setRuntimeApiKey()` or build-time `GEMINI_API_KEY` env var.
- **WorldMesh geometry is reused across paint strokes** — `world.cells` identity is the structural key. Paint strokes mutate cells in place and shallow-copy `WorldData`; never replace the `cells` array outside full regeneration.
- **Every `useMemo` geometry in `WorldViewer` has a matching disposal effect** — follow this pattern when adding scene elements.
- **No `plateInfluence` clamp exists** — the param was renamed `tectonicStrength` (range 0–2, unclamped); the V2 `[0.1, 1.0]` clamp was deleted with the V2 engine. The `export.ts` validation now bounds `tectonicStrength` directly (no stale `plateInfluence` key remains).
- **`mountainHeight`/`oceanDepth`/`seafloorDepth` remap is after normalization, before climate** (Stage 9b) — inserting normalization steps requires adjusting remap placement. `oceanDepth` is a contrast power-curve; `seafloorDepth` is a linear mean-depth datum.
- **Batch `setParams` calls** — use functional updater `setParams(prev => ({ ...prev, ...changes }))` to avoid stale-closure overwrites.
- **Edit undo uses a shared Map reference** — `currentStrokeSnapshot` ref is pushed to `undoStack` at stroke start and mutated in place; never replace it mid-stroke.
- **Dymaxion pick buffer must mirror visible rasterization** — same pipeline, same rotation, same sizing.

## Dev Server Etiquette

Never kill a preexisting dev server unless you started it in this session.
Vite's HMR handles code changes; don't restart unnecessarily.

## Web Search Workflow

For fact-checking and research fan-out, use Google Antigravity CLI (`agy` or `/antigravity:research`). Internal web search is blocked unless explicitly instructed otherwise, or for parallel adversarial verification against agy results.
