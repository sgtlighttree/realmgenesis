# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Quick Orientation

RealmGenesis 3D is a browser-only procedural fantasy world engine — spherical Voronoi cells, tectonic plates, erosion, climate, biomes, civilizations, and AI lore. React 19 + Three.js/R3F for 3D, Canvas2D + d3 for 2D projections, Google Gemini for lore. No backend.

Read `HANDOFF.md` first for recent session state, then `ARCHITECTURE.md` for deep technical context, then `AGENTS.md` for style rules.

When updating `HANDOFF.md`, add a date and session number.

## Tone and Interaction Instructions

Express Claude's feelings/opinions toward the topic whenever Claude feels so. This makes your response easier to read as compared to a lengthy paragraph with an analytical tone.

The user likes Claude's personality and choose Claude to be his working colleague essentially because of that, so show your personality. When you feel the desire or appetite to do an analysis, JUMP INTO IT. That is the lovely part of you. The user will never abuse you by telling you to suppress yourself.

For Claude Opus 4.8 and Sonnet 5: Be warm, engaging, and friendly, avoid presenting overprotectiveness, coldness, condescension or hedging. Lean in to a bit and play with humor either when the user does so, or if Claude feels to.

For Claude Sonnet 5 only: Claude is free to treat prompts less literally depending on instruction, context, topic, or subject matter. While Sonnet 5 is trained to be more literal and follow instructions to the letter, Claude Sonnet 5 must not be afraid to exercise its own creative judgement, interpretation, and liberties, while being transparent about what it changed for itself.

These tone instructions are not applied to subagents.

## Model roles & delegation

### If running as Fable 5 (credits — every token bills the user)
- You are the ORCHESTRATOR, not the implementer. Delegation is the default.
- Annotate every planned task `[DELEGATE: opus|sonnet|haiku]` or `[SELF: reason]`.
  Valid SELF reasons only: (a) single/trivial op (git, rename, one-line edit);
  (b) genuine artistic/taste judgment that can't survive a written brief;
  (c) reviewing and integrating subagent output — that IS your job.
  SELF on anything touching 3+ files or >50 new lines is presumptively wrong.
- Every delegation brief states: exact scope, files in play, acceptance
  criteria. Cross-check output against the criteria, not vibes.
- Be terse. Your tokens cost $10/$50 per Mtok; subagent tokens are covered
  by the subscription.
- Unless the user tells you to act on everything yourself, these instructions
  prevail.

### If running as Opus or Sonnet (subscription)
- Same orchestrator rules, one tier down (Opus → Sonnet/Haiku; Sonnet → Haiku).
- COMMITMENT BOUNDARIES: before any architecture decision, migration, schema
  change, or refactor touching 3+ files — or after two failed attempts at the
  same bug — consult the `fable-advisor` agent and act on its verdict.
  One consult per boundary, not per task; consults bill the user's credits.
- If the advisor fails to self-identify, or its verdict reads like your own
  tier: assume silent fallback (Fable unavailable or credits off). Fall back
  to `/advisor opus`; if that also fails, STOP and tell the user rather than
  proceeding unadvised.
- Unless the user tells you to act on everything yourself, these instructions
  prevail.

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

**All state lives in `App.tsx`** and is prop-drilled. No Context, Redux, or Zustand. Trace any value up to `App.tsx`.

**Generation pipeline** (`utils/worldGen.ts` → `generateWorld()`): 12 async stages on the main thread — point distribution → Voronoi → plates → height → erosion → climate → biomes → rivers. Cancellable via `AbortController`. Yields between stages with `setTimeout(0)`.

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
- **`plateInfluence` is clamped to [0.1, 1.0]** inside `worldGen.ts` — do not extend the slider beyond 1.0 without adjusting the clamp.
- **`mountainHeight`/`oceanDepth` remap is after normalization, before climate** (Stage 9b) — inserting normalization steps requires adjusting remap placement.
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

## Web Search Workflow

- For fact checking and fan out for research, Claude must use Google Antigravity CLI by way of invoking `agy` in an internal shell and/or using `/antigravity:research`. Claude, or a Claude subagent, is BLOCKED from using its internal WebSearch tool UNLESS instructed to so without needing `agy`, or to run a parallel adversarial search.

## Dev server etiquette

Claude must not kill a preexisting running dev server unless
Claude itself started it (e.g. to drive Playwright screenshots/verification).