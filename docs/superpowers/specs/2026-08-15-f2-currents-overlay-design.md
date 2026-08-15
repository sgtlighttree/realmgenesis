# F2 — Screen-space overlay foundation + ocean-current visualization (v1)

**Status:** design approved 2026-08-15 (Session 16). Architectural path.
**Depends on:** D2 ocean currents (shipped S15), which computes the velocity field
this feature persists and draws.
**Advisor:** `fable-advisor` reviewed the architecture 2026-08-15 (verdict folded in
below, marked ▸).

---

## 1. Goal & scope

Visualize the D2 ocean-current velocity field + SST (sea-surface-temperature)
anomaly as a toggleable overlay on both the 3D globe and the 2D map. Beyond the
feature itself, this lays the **foundation of a screen-space 2D dataviz overlay
layer** to replace the current pattern where every globe overlay is a physical R3F
3D object (the ROADMAP F2 critique: 3D overlays "affect visibility and accuracy").

**In scope (v1):**
- Persist the current field on `WorldData` (optional, gen-time).
- A generic `ScreenOverlay` screen-space layer (projection + horizon occlusion +
  DPR/resize lifecycle + tenant registry).
- **Two tenants** on that layer, to prove the abstraction generalizes:
  1. **Currents** — static arrows along the velocity field, warm/cold SST tint.
  2. **Graticule** — migrated off its existing 3D `lineSegments` onto `ScreenOverlay`.
- Currents drawn on the **2D map** via `Map2D`'s existing composite pass.
- A "Currents" overlay toggle in Controls, wired through `useWorldEngine`.

**Explicitly out of scope (deferred, queued in ROADMAP):**
- Animated particle advection (nullschool-style). ▸ Defer — continuous redraw
  thermals the M1 Air; clean later increment once the layer exists.
- Migrating borders / rivers / roads-routes / contours / labels / Dymaxion onto
  `ScreenOverlay`. Documented as a tenant queue in ROADMAP; each is its own increment.
- Any change to the D2 climate coupling or the currents solver itself.

---

## 2. Key decisions (with advisor verdicts)

**2.1 Persist, don't recompute.** ▸ Correct. The field is computed once in the
worker and currently discarded (`utils/worldGen.ts` ~539–543); we stop discarding
it. Recompute-on-toggle would duplicate a 24-pass relaxation on the main thread and
drag `currents.ts`/`seasons.ts`/`planetary.ts` into the render bundle.
- **Two traps the advisor caught:**
  - The field is **optional** on `WorldData` — absent when `currentStrength === 0`
    (preserves the D2 byte-identical escape hatch; nothing to draw anyway).
  - **Exclude it from JSON export** — a `Float32Array` JSON-serializes into
    index-keyed object garbage. It regenerates deterministically, matching the
    existing save philosophy (geometry is never serialized).

**2.2 Globe overlay = screen-space 2D (Option B), analytic occlusion.** ▸ Option B.
- **Occlusion: analytic horizon test, NOT depth-buffer readback.** Deciding fact
  from code: max terrain displacement is **1.05× radius** (`WorldViewer.tsx:377`),
  so relief self-occlusion is a ~5% grazing-angle edge case. The horizon test
  `dot(normalize(camPos − p), n) > ε` with **ε ≈ 0.08** absorbs that edge and is
  visually indistinguishable from correct. `readPixels` per frame stalls the GPU
  pipeline on the M1 — hard no.
- **Camera sync is trivially feasible in R3F:** a component *inside* the `<Canvas>`
  uses `useFrame` + `Vector3.project(camera)` against the live view/projection
  matrices, drawing to a **sibling absolutely-positioned 2D canvas**
  (`pointer-events: none`, DPR-scaled).
- **Redraw only on change** — camera matrix moved OR `world`/overlay identity
  changed — not unconditionally every frame.

**2.3 The half-migration trap → two tenants, not one.** ▸ A single tenant proves
nothing about the abstraction. v1 lands currents **and** migrates the graticule (the
simplest existing overlay, already `depthTest={true}`) onto `ScreenOverlay` in the
same deliverable. Borders/rivers/labels/Dymaxion stay 3D and are a *documented
queue* in ROADMAP, so the two-idiom state is intentional and bounded, not drift.

**2.4 Static arrows for v1.** ▸ Defer animation (see 1, out of scope).

**2.5 Presentation = overlay toggle, not a ViewMode.** Currents composite over any
existing view mode (biome, satellite, temperature, …), matching how the map is
actually read. Overlay-visibility is new state in `useWorldEngine`, not a
`ViewMode` enum value.

---

## 3. Architecture & data flow

```
generateWorld (worker)                     main thread / render
─────────────────────                      ────────────────────
computeOceanCurrents  ─┐
computeSstAnomaly     ─┴─► WorldData.currents?   ──structured clone──►  world.currents
  (only if strength>0)     { vx,vy,vz, sst:                                   │
                             Float32Array }                                   ├─► ScreenOverlay (globe)
                                                                              │     ├─ tenant: currents (arrows + SST tint)
                                                                              │     └─ tenant: graticule (migrated)
                                                                              └─► Map2D composite pass (currents)

useWorldEngine:  showCurrents (bool) overlay state  ──►  Controls toggle
```

### 3.1 Data model — `WorldData.currents`

```ts
// types.ts
export interface OceanCurrentData {
  vx: Float32Array;   // per-cell tangential velocity, unit-sphere frame
  vy: Float32Array;
  vz: Float32Array;
  sst: Float32Array;  // per-cell SST anomaly, °C (ocean cells only; 0 elsewhere)
}
export interface WorldData {
  // …existing fields…
  currents?: OceanCurrentData; // present iff params.currentStrength > 0
}
```
- Populated in `generateWorld` by keeping the already-computed `currentField` +
  `sstAnomaly` instead of dropping them.
- `export.ts`: strip `currents` before `JSON.stringify` on save; never expect it on
  load (regenerates).
- `paramLiveness`: unaffected — `currentStrength` is already exercised via
  temperature/moisture. (No new liveness case needed; the field is derived data.)

### 3.2 `ScreenOverlay` layer (new — the architecture)

A generic component living inside the R3F `<Canvas>` subtree:
- Owns a sibling 2D canvas (absolutely positioned over the WebGL canvas,
  `pointer-events: none`, sized to the drawing buffer × DPR).
- Each qualifying frame: computes camera position, iterates cells, applies the
  horizon test, projects survivors to screen space, and dispatches to registered
  **tenant draw callbacks** with `(ctx, projectedPoints, visibilityMask, world)`.
- Tenants are pure draw functions — they never own scene state or disposal.
- Lifecycle: resize + DPR handling centralized here; no per-tenant geometry, so the
  "every useMemo geometry needs a disposal effect" invariant doesn't apply to
  tenants (there's no geometry) — a deliberate simplification over the 3D pattern.

**Tenant interface (draft):**
```ts
interface OverlayTenant {
  id: string;
  visible: boolean;
  draw(ctx: CanvasRenderingContext2D, proj: ProjectedCells, world: WorldData): void;
}
// ProjectedCells: { x: Float32Array; y: Float32Array; visible: Uint8Array; n: number }
```

### 3.3 Tenants (v1)

- **Currents tenant:** for each visible ocean cell, draw a short arrow from the
  projected center along the screen projection of `(vx,vy,vz)`, length ∝ speed,
  color lerped warm↔cold by `sst`. Decluttering/subsampling by zoom to keep the
  field legible (reuse the label-declutter spirit; simplest first pass = draw every
  Nth ocean cell by speed threshold).
- **Graticule tenant:** replaces the 3D `lineSegments` graticule — draw the 10°
  meridians/parallels as 2D polylines through projected sample points, with the same
  horizon culling. The old 3D graticule component + its geometry/disposal are removed.

### 3.4 2D map (`Map2D`)

Extend the existing composite pass to draw the currents field using `Map2D`'s own
sphere→screen projection (Mercator/Dymaxion). Independent of `ScreenOverlay` (the 2D
map already IS screen-space). Reads `world.currents`.

### 3.5 State & UI

- `useWorldEngine`: `showCurrents: boolean` (+ setter), default off. Overlay state,
  not regen-affecting.
- `Controls`: a "Currents" toggle (near the view/overlay controls). Disabled/greyed
  when `world.currents` is absent (i.e. `currentStrength === 0`), with a hint.

---

## 4. Task graph (parallelization)

```
[0] Lock WorldData.currents contract        in-house, first, tiny
      │
      ├─► [A] Persistence plumbing           DELEGABLE  types.ts, worldGen.ts, export.ts, (paramLiveness check)
      │
      ├─► [B] ScreenOverlay layer            IN-HOUSE   WorldViewer.tsx — the architecture; gates C,D
      │        ├─► [C] currents tenant        in-house   (needs A+B) same file as B
      │        └─► [D] graticule migration    in-house   (needs B)   same file as B
      │
      └─► [E] Map2D currents drawing         DELEGABLE  Map2D.tsx — needs only A's data contract, NOT B
     [F] UI toggle + state                   in-house, small, folds in at end
```

- **Parallel branches** after [0]: **[A]** (persistence) ‖ **[B]** (ScreenOverlay,
  in-house) ‖ **[E]** (Map2D) — disjoint files (A: types/worldGen/export; E: Map2D;
  B: WorldViewer), safe to run concurrently.
- **Serial spine (in-house):** [0] → [B] → [C]+[D], all in `WorldViewer.tsx`.
- **Delegation executor options:** Sonnet, or `opencode-go/deepseek-v4-flash`
  (systems-checked S16 — scratch-dir only, so it generates against a spec in scratch
  and the diff is reviewed/applied by hand). [A] and [E] are the well-specified
  mechanical units that suit either.

---

## 5. Testing

- **Persistence determinism:** same seed → identical `currents` arrays; `currents`
  absent at `currentStrength === 0`; present otherwise.
- **Horizon test unit:** a back-hemisphere cell is culled, a front-facing one kept,
  grazing cells within ε behave.
- **Export round-trip:** save→load drops `currents` cleanly and the world still
  renders (regenerates the field).
- **Gate suite:** typecheck 0, lint within ratchet, full vitest green (watch the
  documented M1 parallel-load flake — re-run suspected failures in isolation).
- **In-browser (Playwright, manual-render convention):** globe overlay occludes the
  far hemisphere correctly; toggle on/off; 2D map draws the field; graticule still
  reads correctly post-migration; 0 console errors.

---

## 6. Risks

- **Two-idiom window** (some overlays 3D, some screen-space) until the queue is
  worked. Bounded by the documented ROADMAP tenant queue; acceptable.
- **Arrow declutter at 200k cells** — naive per-cell draw is too dense; subsample by
  speed/zoom. First pass keeps it simple, refine in-browser.
- **Camera-sync jank** if redraw isn't gated to actual change — gate on camera
  matrix + world identity.
- **DPR / retina sizing** mismatches between the WebGL canvas and the overlay canvas
  — centralize in `ScreenOverlay`, verify on the M1 retina display.
