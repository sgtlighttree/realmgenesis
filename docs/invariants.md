# Invariants & Gotchas

Non-obvious facts that break something if you violate them. Each was
re-verified against current code (not copied from the legacy monolith, which
drifted). Where the legacy `ARCHITECTURE.md` was wrong, the correction is
called out inline.

The short list also lives in the root `CLAUDE.md` "Key Invariants" section for
at-a-glance reference; this file is the reasoned long form.

---

## Generation & threading

### 1. Generation runs in a Web Worker — NOT the main thread
`hooks/useWorldEngine.ts` calls `generateWorldInWorker()`
(`utils/worldGenClient.ts`), which posts to `workers/worldGen.worker.ts`. The
worker runs `generateWorld` and returns a serialized world that the main thread
rehydrates into `Cell[]`.

> **Correction:** legacy invariant #1 said "No Web Workers are used." That was
> true through Session 6 and is now false (D6 Stage 1, Session 7).

Consequences that still bite:
- **Abort is `worker.terminate()`, not a message.** A worker in a synchronous
  generation loop can't drain its message queue, so one worker is spawned per
  generation. Spawn (~1–5 ms) is negligible against a multi-second run.
- **Never route a *partial* recalc (civ-only, province-only) through the worker
  rehydration path.** Rehydration mints a **new `cells` array** every call;
  `WorldViewer`'s mesh is keyed on `world.cells` identity, so a partial recalc
  sent through it forces a full geometry rebuild — a frame-rate regression with
  no obvious cause. Full regeneration is the only correct worker caller today.
- `generateWorldInWorker` cannot load under Vitest (`?worker` import), so its
  production `finish`/settle wiring has **no CI coverage**. Re-run
  `dev/goldenCompare.html` by hand after touching that promise plumbing — it has
  already shipped one hang (a throw inside the `onmessage` handler escaping the
  Promise executor, leaving `isGenerating` stuck true).

### 2. `mountainHeight` / `oceanDepth` / `seafloorDepth` remap is Stage-9b — after all normalization, before climate
`utils/worldGen.ts` remaps below/above-`seaLevel` heights with a power curve
(`mountainHeight`, `oceanDepth`) plus a linear ocean-floor datum
(`seafloorDepth`), holding the coastline fixed. It must stay **after** the final
height normalization (or the curve distorts an already-biased distribution) and
**before** climate/biomes (Stage 10+), because biome assignment uses `seaLevel`
as a hard threshold on the final height values. If you insert a normalization
step after erosion, move the remap after it too.

> **Correction:** legacy invariant #19 named only `mountainHeight`/`oceanDepth`.
> `seafloorDepth` was added to the same block (2026-08-14). See
> [params-reference.md](params-reference.md) for the exact formula and how the
> two ocean knobs differ (contrast curve vs linear datum).

### 3. There is no `plateInfluence` clamp any more
The V3 terrain model uses `tectonicStrength` (renamed from `plateInfluence` in
Session 8) with **no `[0.1, 1.0]` clamp**. A stale validation bound
`plateInfluence: [0, 2.0]` still lingers in `utils/export.ts` but keys nothing
in the engine.

> **Correction:** legacy invariant #18 (and the CLAUDE.md copy) claimed
> `plateInfluence` is hard-clamped to `[0.1, 1.0]` inside `worldGen.ts`. That
> code path was deleted with the V2 engine (Session 11). Do not reintroduce the
> clamp expecting it to exist.

### 4. Determinism is by isolated RNG side-streams
Every subsystem draws from `new RNG(seed + '_<purpose>')` (e.g. `_macro_v3`,
`_crust`, `_plates_v3`, `_names`, `_civs`). Adding a new random draw **must**
use a fresh suffix stream, or you shift the draw order of an existing stream and
silently change every existing seed's output. Geographic names key off the
terrain `seed`; civ names off `civSeed`.

---

## State & React

### 5. All state lives in `hooks/useWorldEngine.ts`
The `useWorldEngine()` hook owns `params`, `world`, view/display/inspect modes,
edit state, and every handler, and returns them as one object. `ShellApp`
(the default route) consumes the hook and prop-drills into `Controls`,
`WorldViewer`/`Map2D`, etc. No Context, Redux, or Zustand.

> **Correction:** legacy docs said "all state lives in `App.tsx`." `App.tsx` is
> now the **legacy classic route** (`?shell=classic`); the live app is
> `components/shell/ShellApp.tsx` on the hook. `DEFAULT_PARAMS` lives in
> `useWorldEngine.ts`, not `App.tsx`.

### 6. Batch multi-field `setParams` with the functional updater
`setParams(prev => ({ ...prev, a, b }))`, never two sequential
`setParams({ ...params, a })` calls — both close over the same stale `params`
and the second overwrites the first.

### 7. Lore mutates `world.civData` in place
`services/gemini.ts` `generateWorldLore()` writes names directly onto existing
`FactionData`/`ProvinceData` objects. The caller must `setWorld({ ...world })`
(shallow copy of `WorldData`) to trigger a re-render. Do **not** replace
`world.civData` with a new object.

### 8. Edit-undo uses a shared, grow-in-place `Map` reference
`currentStrokeSnapshot` (a ref) holds a `Map<cellId, beforeState>` pushed to
`undoStack` at stroke **start**; subsequent stroke events add to the same Map,
so the stack entry grows by reference. The `end` phase only nulls the ref. Never
replace `currentStrokeSnapshot.current` mid-stroke — only add to it. This makes
undo reliable even if the pointer-up is missed off-canvas.

### 9. Paint mutates `world.cells` in place, then shallow-copies `WorldData`
`world.cells` identity is the structural key for `WorldViewer`'s mesh and the
`Map2D` Dymaxion pick buffer. Paint strokes mutate cells in place and
`setWorld({ ...world })` — `cells` stays stable, no geometry reallocation.
Never replace the `cells` array outside a full regeneration.

---

## Rendering

### 10. `seaLevel` must be passed to `getCellColor`
Signature is `getCellColor(cell, viewMode, seaLevel, factionColors?,
cultureColors?, religionColors?)` (`utils/colors.ts`, verified at
`components/Map2D.tsx`). The third arg must be `world.params.seaLevel`, not a
hardcoded value, or ocean/land color boundaries are wrong in every mode.

### 11. Political/province/culture/religion modes need a color map
Pass `factionColors` (build via `buildFactionColorMap(civData)` from
`utils/colors.ts`) to any `getCellColor` call in political or province mode, or
user color edits won't appear. The same pattern holds for culture/religion
modes with their respective color maps (`cultureColors`, `religionColors` from
the hook). Non-React paint paths (`MiniMap`, `export.ts`, `exportGLB.ts`) build
the map themselves.

### 12. R3F element names are strings — intentional
`WorldViewer.tsx` writes Three.js elements as string literals (`"bufferGeometry"`
etc.) to bypass TSX element typing. `@typescript-eslint/no-explicit-any` is
`warn`, not `error`, largely for this file. Don't "fix" the strings.

### 13. WorldMesh geometry is allocated once per `world.cells` identity
The `BufferGeometry` is created once and its position/color buffers are refilled
in place on paint/view changes — no per-stroke allocation. The mesh has **no
normal attribute** (unlit basic material + `flatShading` standard material both
ignore it) and a **fixed bounding sphere (r = 1.1)**. If you add a lit,
smooth-shaded material or geometry past r = 1.1, revisit both.

### 14. Every `useMemo` geometry in `WorldViewer` has a matching disposal effect
`useEffect(() => () => geo.dispose(), [geo])`. Follow this for any new scene
element or GPU memory grows until context loss — geometries are not GC'd.

### 15. Dymaxion pick buffer must mirror the visible rasterization
The `Map2D.tsx` pick buffer encodes cell IDs as RGB and is generated through the
**same** source-canvas flip, `buildDymaxionNet` mapping, rotation constants, and
sizing as the visible raster. It depends on `world.cells` (stable), not `world`,
so it skips per-stroke reprojection. Don't revert to lon/lat nearest-center
picking — it drifts from the rendered map.

### 16. 3D faction labels are curved surface meshes, not sprites
`CountryLabels` renders canvas-text textures on subdivided mesh patches just
above the globe surface, so they rotate with the planet and are occluded by it.
Don't swap them for HTML overlays or billboards unless UI-like floating labels
are explicitly wanted.

### 17. Cell-boundary seams are closed by stroking each polygon with its fill color
Canvas2D anti-aliases path edges, leaving ~1px background seams between adjacent
Voronoi cells (they line up horizontally in equirectangular). Every cell is
filled **and** stroked (`lineWidth = 1`) with the same hex. Applies to `Map2D`
(Mercator + Dymaxion source) and `export.ts`. Don't remove the stroke.

---

## Data & persistence

### 18. Cell `id` == index in `WorldData.cells`
Stable within one generated world (cells are never reordered), not stable across
runs. `WorldData.geoJson` feature index `i` corresponds to `world.cells[i]` and
is cached from Stage 2 for the world's lifetime — don't clear or replace it.

### 19. Gemini API key is ephemeral
Stored only in hook state / `setRuntimeApiKey()`; never written to
`localStorage`, cookies, or any store; resets to empty on reload. A build-time
`GEMINI_API_KEY` bakes into the public bundle — prefer runtime BYOK for shared
deploys. Don't add persistence.

### 20. Relative imports only
The `@/` alias is configured in `tsconfig.json`/`vite.config.ts` but
**intentionally unused**. Use `../types`, `./colors`. Never add `@/` imports.

### 21. GLB vertex colors need a Blender material step
`exportGLB` bakes cell colors into the GLTF `COLOR_0` attribute. Blender doesn't
show them until the material's Base Color is wired to an Attribute node
(`COLOR_0`) or viewport shading is set to Vertex Color.

---

## Quality gates

### 22. Four CI gates, all must pass
`npm run lint` (0 errors, `--max-warnings` ratchet), `npm run typecheck`
(`tsc --noEmit`, `"strict": true`), `npm test` (Vitest), `npm run build`. The
**param-liveness** test fails if any tunable `WorldParams` key stops influencing
output — extend `tests/paramLiveness.test.ts` when adding a param. See
[testing.md](testing.md) for the determinism-instrument model and the known M1
parallel-load flake.
