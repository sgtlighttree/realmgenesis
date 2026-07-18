# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Quick Orientation

RealmGenesis 3D is a browser-only procedural fantasy world engine — spherical Voronoi cells, tectonic plates, erosion, climate, biomes, civilizations, and AI lore. React 19 + Three.js/R3F for 3D, Canvas2D + d3 for 2D projections, Google Gemini for lore. No backend.

Read `HANDOFF.md` first for recent session state, then `ARCHITECTURE.md` for deep technical context, then `AGENTS.md` for style rules.

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
