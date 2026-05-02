# AGENTS.md -- RealmGenesis 3D

## Session Handoff

Read [HANDOFF.md](./HANDOFF.md) first for a summary of what was built recently, known issues, and what to tackle next. It is updated at the end of each working session and is the fastest way to orient to the current state of the project.

## Architecture Overview

Before making significant changes, read [ARCHITECTURE.md](./ARCHITECTURE.md) for a complete overview of the codebase — the generation pipeline, data model, rendering architecture, state management, and key invariants. It is designed so an LLM (or a new contributor) can navigate the codebase confidently without reading every file.

Key entry points:

- **All app state**: `App.tsx` — 24+ `useState` variables + edit-mode refs, prop-drilled to all children
- **Data types**: `types.ts` — `Cell`, `WorldData`, `WorldParams`, `BiomeType`
- **Generation logic**: `utils/worldGen.ts` — `generateWorld()` (12-stage async pipeline, line 491)
- **Color mapping**: `utils/colors.ts` — `getCellColor(cell, viewMode, seaLevel, factionColors?)`
- **Map painting**: `utils/paintUtils.ts` — brush BFS, stroke functions, undo snapshots
- **Edit toolbar**: `components/EditToolbar.tsx` — paint/edit mode HUD
- **3D rendering**: `components/WorldViewer.tsx`
- **2D rendering**: `components/Map2D.tsx`
- **AI lore**: `services/gemini.ts` — `generateWorldLore()` mutates `world.civData` in-place
- **Save/export**: `utils/export.ts` — `exportMap()`, `saveMapToBrowser()`, `saveMapConfig()`

---

## Commands

```bash
npm run dev       # Start Vite dev server on port 3000
npm run build     # Production build → dist/
npm run preview   # Preview production build
npm run lint      # ESLint (flat config: eslint.config.js)
```

**No test framework or formatter is configured.** This is a visual, interactive app — testing is done by running `npm run dev` and verifying behavior manually in the browser. The quality gates are:
1. `npm run build` must succeed with no errors
2. `npm run lint` must produce zero errors (warnings are acceptable)
3. TypeScript compilation must produce no type errors

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
- Strict-ish mode: `skipLibCheck: true`, `allowJs: true`, `noEmit: true`
- Use `interface` for object types and component props, `type` for unions
- Prefer explicit return types on exported functions
- Use `as any` casts sparingly (only when required by library typing gaps, e.g., R3F element names in `WorldViewer.tsx`)
- `@typescript-eslint/no-explicit-any` is a warning, not an error — acceptable for R3F event handlers and d3 projections

### React Components
- Functional components with `React.FC<Props>` type annotation
- Props defined as `interface ComponentProps { ... }`
- Use `useCallback` for event handlers passed as props, `useMemo` for expensive computations
- `useState` for local state, all app-level state lives in `App.tsx`
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
- Tailwind CSS via CDN — use utility classes exclusively, no CSS modules or styled-components
- Dark theme: `bg-gray-950`, `text-gray-200`, `border-gray-800`
- Responsive: mobile-first with `md:` breakpoints
- Overlays: `backdrop-blur-md`, `bg-black/50`, `border-white/10`

### Architecture Patterns
- Single source of truth: `App.tsx` holds all state, passes down via props
- Utils are pure functions (no side effects except logging callbacks)
- Services (e.g., `gemini.ts`) wrap external APIs with minimal abstraction
- Components are presentational — no data fetching or generation logic
