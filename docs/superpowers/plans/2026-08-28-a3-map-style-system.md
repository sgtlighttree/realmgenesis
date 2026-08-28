# A3 Map Style System (Parchment) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a parchment map style — bare-paper land, hand-drawn relief glyphs, hatched ocean, paper grain — rendering identically in Map2D, PNG export and SVG export.

**Architecture:** A style is an ordered list of draw passes written once against a `Substrate` adapter with two implementations (Canvas2D and SVG-string). Glyph placement is split out as a pure, substrate-independent function of `(cells, projection, widthPx)`, because glyph sizing and collision are screen-space decisions that would otherwise be written twice.

**Tech Stack:** TypeScript (strict), React 19, d3-geo (`d3.GeoProjection` / `d3.geoPath`), Canvas2D, SVG, Vitest.

**Spec:** `docs/superpowers/specs/2026-08-28-a3-map-style-system-design.md`

## Global Constraints

- **Relative imports only** — the `@/` alias is configured but intentionally unused.
- 2-space indent, semicolons, single quotes, trailing commas.
- Import order: React → external → local (relative), blank line between groups.
- `interface` for objects/props, `type` for unions. Components are `React.FC<Props>`, functional only.
- Tailwind utility classes only. Dark theme: `bg-gray-950`, `text-gray-200`, `border-gray-800`.
- Naming: PascalCase components/types, camelCase functions/vars, SCREAMING_SNAKE enum values.
- **Style state must NOT go into `WorldParams`.** It is a render choice, like `viewMode`. Adding it to `WorldParams` would break `tests/paramLiveness.test.ts`, which fails when a param stops influencing *generated* output.
- **The 3D globe (`WorldViewer.tsx`) is out of scope.** It keeps its current rendering.
- Gates that must pass before every commit: `npm run typecheck` (0 errors), `npm run lint` (0 errors, ≤30 warnings), `npm test`.
- `tests/paramLiveness.test.ts` is a **load canary** — if it times out, run `uptime` and re-run it isolated before suspecting the code.

---

### Task 1: Collapse `getCellColor` to a `ColorContext` object

Pure refactor, zero behaviour change, its own commit before any style code. Mixing it with the feature would give a rendering regression two candidate causes.

**Files:**
- Modify: `utils/colors.ts:93` (the `getCellColor` signature)
- Modify: `components/Map2D.tsx:372`, `components/Map2D.tsx:603`
- Modify: `components/WorldViewer.tsx:527`
- Modify: `components/MiniMap.tsx:44`
- Modify: `components/DymaxionPreview2D.tsx:127`
- Modify: `components/DymaxionNetPreview.tsx:64`
- Modify: `utils/export.ts:115`, `utils/export.ts:372`
- Modify: `utils/exportVector.ts:102`
- Modify: `utils/exportGLB.ts:14`
- Modify: `CLAUDE.md` (two invariants keyed to the positional form)

**Interfaces:**
- Consumes: nothing.
- Produces: `ColorContext` (exported from `utils/colors.ts`) and the new
  `getCellColor(cell: Cell, mode: ViewMode, ctx: ColorContext): THREE.Color`.
  Every later task passes a `ColorContext`.

**Note:** `grep -rn "getCellColor" tests/` returns nothing — no test encodes the argument order. The compiler is the verification here, plus the full suite for behaviour.

- [ ] **Step 1: Add the `ColorContext` type and change the signature**

In `utils/colors.ts`, above `getCellColor`:

```ts
/**
 * Everything `getCellColor` needs beyond the cell and the view mode.
 *
 * Replaces a 7-positional-argument signature that had grown two documented
 * footguns (a silently-wrong third argument, and an omitted faction map that
 * rendered political mode blank). An object makes both impossible: `seaLevel`
 * is required and named, and the optional maps are named at every call site.
 */
export interface ColorContext {
  seaLevel: number;
  factionColors?: Map<number, string>;
  cultureColors?: Map<number, string>;
  religionColors?: Map<number, string>;
  seasonalDelta?: number;
}
```

Change the signature to:

```ts
export const getCellColor = (cell: Cell, mode: ViewMode, ctx: ColorContext): THREE.Color => {
  const { seaLevel, factionColors, cultureColors, religionColors, seasonalDelta } = ctx;
  const color = new THREE.Color();
  // ...body unchanged from here down...
```

The body is untouched — destructuring restores every name it already uses.

- [ ] **Step 2: Verify the compiler finds every call site**

Run: `npm run typecheck`
Expected: FAIL, with one error per call site listed in **Files** above (10 errors across 8 files). Confirm the count matches; a missing file means a call site was missed.

- [ ] **Step 3: Update all 10 call sites**

Each becomes the same shape. Example for `components/Map2D.tsx:372`:

```ts
const color = getCellColor(world.cells[i], viewMode, {
  seaLevel: world.params.seaLevel,
  factionColors,
  cultureColors,
  religionColors,
  seasonalDelta: seasonalTemperatureDelta(world.cells[i], world.params),
});
```

`components/DymaxionPreview2D.tsx:127` currently passes only `seaLevel`, so it becomes:

```ts
const color = getCellColor(cell, viewMode, { seaLevel: world.params.seaLevel });
```

- [ ] **Step 4: Rewrite the two CLAUDE.md invariants**

In `CLAUDE.md` under **Key Invariants**, replace:

> - **`seaLevel` must be passed to `getCellColor`** as the third argument (from `world.params.seaLevel`), not hardcoded.
> - **`factionColors` map required for political rendering** — any render path calling `getCellColor` for political/province mode must pass a faction-color map (build via `buildFactionColorMap(civData)` from `colors.ts`).

with:

> - **`getCellColor(cell, mode, ctx)` takes a `ColorContext` object**, not positional arguments. `ctx.seaLevel` is required and comes from `world.params.seaLevel` — never hardcoded. For political/province/culture/religion modes the matching colour map must be present on the context (build via `buildFactionColorMap(civData)` / `buildCultureColorMap` / `buildReligionColorMap` from `colors.ts`); omitting it renders those modes blank.

- [ ] **Step 5: Run the gates**

Run: `npm run typecheck && npm run lint && npm test`
Expected: typecheck 0 errors; lint 0 errors, ≤30 warnings; all tests pass. This is a pure refactor, so **no test count changes and no test changes**.

- [ ] **Step 6: Commit**

```bash
git add utils/colors.ts components/Map2D.tsx components/WorldViewer.tsx components/MiniMap.tsx components/DymaxionPreview2D.tsx components/DymaxionNetPreview.tsx utils/export.ts utils/exportVector.ts utils/exportGLB.ts CLAUDE.md
git commit -m "refactor(A3): getCellColor takes a ColorContext object

Seven positional arguments across ten call sites had grown two documented
footguns, both recorded as CLAUDE.md invariants. An object makes seaLevel
required-and-named and the colour maps named at every call site. Both
invariants rewritten in the same commit. Pure refactor, no behaviour change."
```

---

### Task 2: Style types, registry, state and UI selector

Ships the plumbing with only the existing look ported as `default`. Visually a no-op — that is the point. A reviewer can verify the selector appears and changes nothing.

**Files:**
- Create: `utils/mapStyle/types.ts`
- Create: `utils/mapStyle/styleDefault.ts`
- Create: `utils/mapStyle/index.ts`
- Create: `tests/mapStyle.test.ts`
- Modify: `hooks/useWorldEngine.ts:64` (add state beside `viewMode`), and the returned object at `hooks/useWorldEngine.ts:677`
- Modify: `components/shell/ShellApp.tsx:48` (destructure), and the `Controls` render site
- Modify: `components/Controls.tsx` (selector next to the view-mode control)

**Interfaces:**
- Consumes: `ColorContext` from Task 1.
- Produces:
  - `MapStyleId = 'default' | 'parchment'`
  - `interface MapStyle { id; name; palette; fillPolicy(mode): FillPolicy; passes: StylePass[] }`
  - `type FillPolicy = 'bare' | 'categorical' | 'ramp'`
  - `MAP_STYLES: Record<MapStyleId, MapStyle>` and `getMapStyle(id): MapStyle` from `utils/mapStyle/index.ts`
  - `mapStyleId` / `setMapStyleId` on the `useWorldEngine` return object.

- [ ] **Step 1: Write the failing test**

Create `tests/mapStyle.test.ts`:

```ts
import { describe, it, expect } from 'vitest';

import { getMapStyle, MAP_STYLES } from '../utils/mapStyle';

describe('map style registry', () => {
  it('exposes default and parchment', () => {
    expect(Object.keys(MAP_STYLES).sort()).toEqual(['default', 'parchment']);
  });

  it('returns the default style for an unknown id', () => {
    // @ts-expect-error deliberately probing the runtime fallback
    expect(getMapStyle('nonsense').id).toBe('default');
  });

  it('gives every ViewMode a fill policy in every style', () => {
    const modes = [
      'biome', 'height', 'height_bw', 'temperature', 'moisture', 'plates',
      'political', 'population', 'province', 'satellite', 'culture', 'religion',
    ] as const;
    for (const style of Object.values(MAP_STYLES)) {
      for (const mode of modes) {
        expect(['bare', 'categorical', 'ramp']).toContain(style.fillPolicy(mode));
      }
    }
  });

  it('keeps every mode on the ramp policy in the default style', () => {
    // The default style is the pre-A3 look: every mode paints its own fill.
    expect(getMapStyle('default').fillPolicy('political')).toBe('ramp');
    expect(getMapStyle('default').fillPolicy('satellite')).toBe('ramp');
  });
});
```

- [ ] **Step 2: Run it to verify it fails**

Run: `npx vitest run tests/mapStyle.test.ts`
Expected: FAIL — cannot resolve `../utils/mapStyle`.

- [ ] **Step 3: Create `utils/mapStyle/types.ts`**

```ts
import { ViewMode } from '../../types';

/** How a style paints the land fill for a given view mode. See spec §4. */
export type FillPolicy =
  | 'bare'         // no fill; glyphs, coastline and hillshade carry the terrain
  | 'categorical'  // bare paper plus a muted categorical fill on top
  | 'ramp';        // keep the mode's own continuous fill; suppress glyphs

export type MapStyleId = 'default' | 'parchment';

export type GlyphKind = 'mountain' | 'hill' | 'forest' | 'conifer' | 'dune' | 'marsh';

/** A glyph resolved to output pixels. Both substrates just draw this. */
export interface PlacedGlyph {
  x: number;
  y: number;
  kind: GlyphKind;
  scale: number;   // pixels, tip to base
  seedRot: number; // radians, deterministic per cell
  cellId: number;
}

export interface StylePalette {
  paper: string;
  ink: string;
  inkLight: string;
  sea: string;
  seaHatch: string;
  coast: string;
}

export interface MapStyle {
  id: MapStyleId;
  name: string;
  palette: StylePalette;
  /** Per-view-mode land fill rule. See spec §4. */
  fillPolicy: (mode: ViewMode) => FillPolicy;
}
```

- [ ] **Step 4: Create `utils/mapStyle/styleDefault.ts`**

```ts
import { MapStyle } from './types';

/**
 * The pre-A3 look, expressed as a style so the registry has a neutral member.
 * Every mode is 'ramp': each view mode paints its own fill exactly as before,
 * and no style pass draws anything. Selecting this style is a visual no-op.
 */
export const styleDefault: MapStyle = {
  id: 'default',
  name: 'Default',
  palette: {
    paper: '#000000',
    ink: '#ffffff',
    inkLight: '#999999',
    sea: '#050505',
    seaHatch: '#050505',
    coast: '#ffffff',
  },
  fillPolicy: () => 'ramp',
};
```

- [ ] **Step 5: Create `utils/mapStyle/index.ts`**

```ts
import { MapStyle, MapStyleId } from './types';
import { styleDefault } from './styleDefault';
import { styleParchment } from './styleParchment';

export * from './types';

export const MAP_STYLES: Record<MapStyleId, MapStyle> = {
  default: styleDefault,
  parchment: styleParchment,
};

/** Unknown ids fall back to `default` — a saved/stale id must never blank the map. */
export const getMapStyle = (id: MapStyleId): MapStyle => MAP_STYLES[id] ?? styleDefault;
```

- [ ] **Step 6: Create the parchment stub `utils/mapStyle/styleParchment.ts`**

Passes arrive in Task 7. The fill policy is real from the start — it is the spec §4 table.

```ts
import { ViewMode } from '../../types';
import { FillPolicy, MapStyle } from './types';

const BARE: ViewMode[] = ['satellite', 'biome', 'height_bw'];
const CATEGORICAL: ViewMode[] = ['political', 'province', 'culture', 'religion', 'plates'];

/**
 * Spec §4. Continuous-ramp modes (height, temperature, moisture, population)
 * keep their own fill: their entire information content IS the fill, so bare
 * paper would render them blank. Glyphs are suppressed there too — they would
 * fight the ramp.
 */
const parchmentFillPolicy = (mode: ViewMode): FillPolicy => {
  if (BARE.includes(mode)) return 'bare';
  if (CATEGORICAL.includes(mode)) return 'categorical';
  return 'ramp';
};

export const styleParchment: MapStyle = {
  id: 'parchment',
  name: 'Parchment',
  palette: {
    paper: '#e8d9b5',
    ink: '#3b2f1c',
    inkLight: '#8a7550',
    sea: '#c9bf9a',
    seaHatch: '#a89b74',
    coast: '#2e2414',
  },
  fillPolicy: parchmentFillPolicy,
};
```

- [ ] **Step 7: Run the test to verify it passes**

Run: `npx vitest run tests/mapStyle.test.ts`
Expected: PASS, 4 tests.

- [ ] **Step 8: Add `mapStyleId` state to `useWorldEngine`**

In `hooks/useWorldEngine.ts`, immediately after line 64 (`const [viewMode, setViewMode] = useState<ViewMode>('biome');`):

```ts
  // A3: render-only style axis, orthogonal to viewMode. Deliberately NOT a
  // WorldParam — it never influences generation, so paramLiveness would fail it.
  const [mapStyleId, setMapStyleId] = useState<MapStyleId>('default');
```

Import `MapStyleId` from `../utils/mapStyle`. Add `mapStyleId, setMapStyleId,` to the returned object at line 677, next to `viewMode, setViewMode,`.

- [ ] **Step 9: Drill the prop through `ShellApp`**

In `components/shell/ShellApp.tsx`, add `mapStyleId, setMapStyleId,` to the destructure at line 48, and pass `mapStyleId={mapStyleId} setMapStyleId={setMapStyleId}` to the `Controls` render site alongside the existing `viewMode` props.

- [ ] **Step 10: Add the selector to `Controls.tsx`**

Add these props to the `Controls` props interface:

```ts
  mapStyleId: MapStyleId;
  setMapStyleId: (id: MapStyleId) => void;
```

Render it next to the view-mode control:

```tsx
<div className="space-y-1" title="Map style. Changes how the 2D map and exports are drawn. Does not affect the 3D globe or the generated world.">
  <div className="flex justify-between text-xs text-ink-muted">
    <label>Map Style</label>
  </div>
  <select
    value={mapStyleId}
    onChange={(e) => { setMapStyleId(e.target.value as MapStyleId); }}
    className="w-full bg-gray-950 text-gray-200 border border-gray-800 rounded px-2 py-1 text-xs"
  >
    {Object.values(MAP_STYLES).map(s => (
      <option key={s.id} value={s.id}>{s.name}</option>
    ))}
  </select>
</div>
```

- [ ] **Step 11: Run the gates**

Run: `npm run typecheck && npm run lint && npm test`
Expected: typecheck 0; lint 0 errors, ≤30 warnings; test count rises by 4.

- [ ] **Step 12: Commit**

```bash
git add utils/mapStyle tests/mapStyle.test.ts hooks/useWorldEngine.ts components/shell/ShellApp.tsx components/Controls.tsx
git commit -m "feat(A3): style registry, fill policy and Map Style selector

Ships plumbing only — the default style is a visual no-op. Style state lives
beside viewMode in useWorldEngine, deliberately NOT in WorldParams: it never
influences generation, so paramLiveness would fail it.

fillPolicy encodes spec section 4: continuous-ramp modes keep their own fill
because bare paper would render them blank."
```

---

### Task 3: `placeGlyphs` — the substrate-independent placement stage

The highest-value unit-testable piece, and the one with a wrong answer. Full TDD.

**Files:**
- Create: `utils/mapStyle/placeGlyphs.ts`
- Create: `tests/placeGlyphs.test.ts`

**Interfaces:**
- Consumes: `PlacedGlyph`, `GlyphKind` from Task 2's `types.ts`.
- Produces:
  ```ts
  placeGlyphs(cells: Cell[], projection: d3.GeoProjection, widthPx: number, opts?: GlyphOptions): PlacedGlyph[]
  interface GlyphOptions { seed?: string; minSpacingPx?: number; sizePx?: number; seaLevel: number }
  ```
  Tasks 7 and 8 call this.

- [ ] **Step 1: Write the failing tests**

Create `tests/placeGlyphs.test.ts`:

```ts
import { describe, it, expect } from 'vitest';
import * as d3 from 'd3-geo';

import { generateWorld } from '../utils/worldGen';
import { makeParams } from './helpers';
import { placeGlyphs } from '../utils/mapStyle/placeGlyphs';

const WIDTH = 1024;

const projectionFor = (width: number) =>
  d3.geoEquirectangular().fitSize([width, width / 2], { type: 'Sphere' } as d3.GeoPermissibleObjects);

describe('placeGlyphs', () => {
  it('is deterministic for the same world and width', async () => {
    const world = await generateWorld(makeParams({ points: 2000 }));
    const opts = { seaLevel: world.params.seaLevel, seed: world.params.seed };
    const a = placeGlyphs(world.cells, projectionFor(WIDTH), WIDTH, opts);
    const b = placeGlyphs(world.cells, projectionFor(WIDTH), WIDTH, opts);
    expect(a).toEqual(b);
    expect(a.length).toBeGreaterThan(0);
  }, 30000);

  it('never places two glyphs closer than minSpacingPx', async () => {
    const world = await generateWorld(makeParams({ points: 2000 }));
    const minSpacingPx = 24;
    const placed = placeGlyphs(world.cells, projectionFor(WIDTH), WIDTH, {
      seaLevel: world.params.seaLevel, seed: world.params.seed, minSpacingPx,
    });
    for (let i = 0; i < placed.length; i++) {
      for (let j = i + 1; j < placed.length; j++) {
        const d = Math.hypot(placed[i].x - placed[j].x, placed[i].y - placed[j].y);
        expect(d).toBeGreaterThanOrEqual(minSpacingPx);
      }
    }
  }, 30000);

  it('places no glyph on a water cell', async () => {
    const world = await generateWorld(makeParams({ points: 2000 }));
    const placed = placeGlyphs(world.cells, projectionFor(WIDTH), WIDTH, {
      seaLevel: world.params.seaLevel, seed: world.params.seed,
    });
    for (const g of placed) {
      expect(world.cells[g.cellId].height).toBeGreaterThanOrEqual(world.params.seaLevel);
    }
  }, 30000);

  it('scales glyph size with output width, keeping density roughly constant', async () => {
    const world = await generateWorld(makeParams({ points: 2000 }));
    const opts = { seaLevel: world.params.seaLevel, seed: world.params.seed };
    const small = placeGlyphs(world.cells, projectionFor(512), 512, opts);
    const large = placeGlyphs(world.cells, projectionFor(2048), 2048, opts);
    // Size is proportional to width.
    expect(large[0].scale).toBeCloseTo(small[0].scale * 4, 5);
    // Density (glyphs per unit area) stays within a factor of two, because
    // spacing scales with width too — a big export is the same map, not a
    // denser one.
    const ratio = large.length / Math.max(1, small.length);
    expect(ratio).toBeGreaterThan(0.5);
    expect(ratio).toBeLessThan(2.0);
  }, 30000);

  it('emits only known glyph kinds', async () => {
    const world = await generateWorld(makeParams({ points: 2000 }));
    const placed = placeGlyphs(world.cells, projectionFor(WIDTH), WIDTH, {
      seaLevel: world.params.seaLevel, seed: world.params.seed,
    });
    const kinds = new Set(placed.map(g => g.kind));
    for (const k of kinds) {
      expect(['mountain', 'hill', 'forest', 'conifer', 'dune', 'marsh']).toContain(k);
    }
  }, 30000);
});
```

- [ ] **Step 2: Run to verify it fails**

Run: `npx vitest run tests/placeGlyphs.test.ts`
Expected: FAIL — cannot resolve `../utils/mapStyle/placeGlyphs`.

- [ ] **Step 3: Implement `utils/mapStyle/placeGlyphs.ts`**

```ts
import * as d3 from 'd3-geo';

import { BiomeType, Cell } from '../../types';
import { toLonLat } from '../geo';
import { GlyphKind, PlacedGlyph } from './types';

export interface GlyphOptions {
  seaLevel: number;
  /** World seed, so glyph variation is stable across re-renders. */
  seed?: string;
  /** Minimum centre-to-centre distance, in output pixels, at widthPx = 1024. */
  minSpacingPx?: number;
  /** Glyph height in output pixels, at widthPx = 1024. */
  sizePx?: number;
}

const REFERENCE_WIDTH = 1024;
const DEFAULT_SPACING_PX = 22;
const DEFAULT_SIZE_PX = 16;

/** Height above sea level, as a 0-1 fraction of the land band. */
const landFrac = (cell: Cell, seaLevel: number): number =>
  (cell.height - seaLevel) / Math.max(1e-6, 1 - seaLevel);

/**
 * Which glyph a land cell earns, or null for open ground. Relief wins over
 * vegetation: a forested mountain reads as a mountain, the way a drawn map
 * would show it.
 */
const glyphFor = (cell: Cell, seaLevel: number): GlyphKind | null => {
  const f = landFrac(cell, seaLevel);
  if (f > 0.45) return 'mountain';
  if (f > 0.22) return 'hill';
  switch (cell.biome) {
    case BiomeType.BOREAL_FOREST:
      return 'conifer';
    case BiomeType.TEMPERATE_FOREST:
    case BiomeType.TEMPERATE_RAINFOREST:
    case BiomeType.TROPICAL_RAINFOREST:
      return 'forest';
    case BiomeType.HOT_DESERT:
      return 'dune';
    default:
      return null;
  }
};

/**
 * Prominence orders the greedy thinning pass: the most significant feature in
 * a crowded neighbourhood is the one that survives. Relief sorts by elevation;
 * vegetation sorts below all relief so mountains never lose to a tree.
 */
const prominence = (cell: Cell, kind: GlyphKind, seaLevel: number): number => {
  const f = landFrac(cell, seaLevel);
  if (kind === 'mountain' || kind === 'hill') return 1000 + f * 1000;
  return f * 100;
};

/** Deterministic 0-1 hash from a cell id and the world seed. */
const hash01 = (cellId: number, seed: string): number => {
  let h = 2166136261 >>> 0;
  const str = `${seed}:${cellId}`;
  for (let i = 0; i < str.length; i++) {
    h ^= str.charCodeAt(i);
    h = Math.imul(h, 16777619) >>> 0;
  }
  return (h >>> 8) / 0x1000000;
};

/**
 * Resolve every land cell that earns a glyph into output-pixel coordinates,
 * then thin greedily so no two glyphs collide.
 *
 * Substrate-independent BY DESIGN (spec §3.3): glyph size and collision are
 * screen-space decisions, and Mercator vs Mollweide distort area wildly at high
 * latitude. Both substrates draw the returned list, so this logic exists once.
 *
 * Spacing and size scale with `widthPx`, so an 8192px export is the same map as
 * a 1024px one at higher resolution — not a denser map.
 */
export const placeGlyphs = (
  cells: Cell[],
  projection: d3.GeoProjection,
  widthPx: number,
  opts: GlyphOptions,
): PlacedGlyph[] => {
  const seed = opts.seed ?? '';
  const k = widthPx / REFERENCE_WIDTH;
  const minSpacing = (opts.minSpacingPx ?? DEFAULT_SPACING_PX) * k;
  const size = (opts.sizePx ?? DEFAULT_SIZE_PX) * k;
  const heightPx = widthPx / 2;

  interface Candidate extends PlacedGlyph { rank: number }
  const candidates: Candidate[] = [];

  for (const cell of cells) {
    if (cell.height < opts.seaLevel) continue;
    if (cell.biome === BiomeType.LAKE || cell.biome === BiomeType.SALT_LAKE) continue;
    const kind = glyphFor(cell, opts.seaLevel);
    if (!kind) continue;

    const projected = projection(toLonLat([cell.center.x, cell.center.y, cell.center.z]));
    if (!projected) continue;
    const [x, y] = projected;
    if (!Number.isFinite(x) || !Number.isFinite(y)) continue;
    // Clip to the output box. Projections wrap or blow up past the antimeridian
    // and the poles; an off-canvas glyph is wasted work and can smear.
    if (x < 0 || y < 0 || x > widthPx || y > heightPx * 2) continue;

    const r = hash01(cell.id, seed);
    candidates.push({
      x, y, kind, cellId: cell.id,
      scale: size * (0.85 + r * 0.3),
      seedRot: (r - 0.5) * 0.35,
      rank: prominence(cell, kind, opts.seaLevel),
    });
  }

  // Most prominent first, id as a tiebreak so the order is total and stable.
  candidates.sort((a, b) => (b.rank - a.rank) || (a.cellId - b.cellId));

  // Greedy thinning against a uniform grid: only cells within one spacing can
  // collide, so this stays linear instead of quadratic on large worlds.
  const cellSize = Math.max(1, minSpacing);
  const grid = new Map<string, PlacedGlyph[]>();
  const key = (gx: number, gy: number) => `${gx},${gy}`;
  const accepted: PlacedGlyph[] = [];
  const minSq = minSpacing * minSpacing;

  for (const c of candidates) {
    const gx = Math.floor(c.x / cellSize);
    const gy = Math.floor(c.y / cellSize);
    let collides = false;
    for (let dx = -1; dx <= 1 && !collides; dx++) {
      for (let dy = -1; dy <= 1 && !collides; dy++) {
        for (const other of grid.get(key(gx + dx, gy + dy)) ?? []) {
          const ddx = other.x - c.x;
          const ddy = other.y - c.y;
          if (ddx * ddx + ddy * ddy < minSq) { collides = true; break; }
        }
      }
    }
    if (collides) continue;
    const glyph: PlacedGlyph = {
      x: c.x, y: c.y, kind: c.kind, scale: c.scale, seedRot: c.seedRot, cellId: c.cellId,
    };
    accepted.push(glyph);
    const bucket = grid.get(key(gx, gy));
    if (bucket) bucket.push(glyph); else grid.set(key(gx, gy), [glyph]);
  }

  return accepted;
};
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `npx vitest run tests/placeGlyphs.test.ts`
Expected: PASS, 5 tests.

If the spacing test fails, the grid neighbourhood is wrong — a glyph exactly `minSpacing` away must be *accepted*, so the comparison is strict `<`, not `<=`.

- [ ] **Step 5: Run the gates and commit**

Run: `npm run typecheck && npm run lint && npm test`

```bash
git add utils/mapStyle/placeGlyphs.ts tests/placeGlyphs.test.ts
git commit -m "feat(A3): substrate-independent glyph placement

placeGlyphs(cells, projection, widthPx, opts) resolves land cells to output
pixels, ranks by prominence and thins greedily against a uniform grid.

Split out from drawing deliberately (spec 3.3): glyph size and collision are
screen-space decisions, so placement inside a per-cell draw callback would be
written once for Canvas2D and again for SVG. Both substrates draw this list.

Spacing and size scale with widthPx, so a large export is the same map at
higher resolution rather than a denser one."
```

---

### Task 4: Procedural glyph paths

**Files:**
- Create: `utils/mapStyle/glyphPaths.ts`
- Create: `tests/glyphPaths.test.ts`

**Interfaces:**
- Consumes: `GlyphKind`, `PlacedGlyph` from Task 2.
- Produces: `glyphPathData(g: PlacedGlyph): string` — an SVG path `d` string in **output pixel coordinates**, already translated, scaled and rotated. Both substrates consume this: `SvgSubstrate` embeds it, `Canvas2DSubstrate` feeds it to `new Path2D(d)`.

Using one path string for both substrates is what keeps the shapes identical. `Path2D` accepts SVG path syntax natively, so there is no second implementation.

- [ ] **Step 1: Write the failing test**

Create `tests/glyphPaths.test.ts`:

```ts
import { describe, it, expect } from 'vitest';

import { glyphPathData } from '../utils/mapStyle/glyphPaths';
import { GlyphKind, PlacedGlyph } from '../utils/mapStyle/types';

const KINDS: GlyphKind[] = ['mountain', 'hill', 'forest', 'conifer', 'dune', 'marsh'];

const glyph = (kind: GlyphKind, over: Partial<PlacedGlyph> = {}): PlacedGlyph => ({
  x: 100, y: 200, kind, scale: 16, seedRot: 0, cellId: 1, ...over,
});

describe('glyphPathData', () => {
  it('returns a non-empty path for every kind', () => {
    for (const kind of KINDS) {
      const d = glyphPathData(glyph(kind));
      expect(d.length).toBeGreaterThan(0);
      expect(d.startsWith('M')).toBe(true);
      expect(d).not.toMatch(/NaN|Infinity|undefined/);
    }
  });

  it('is deterministic', () => {
    for (const kind of KINDS) {
      expect(glyphPathData(glyph(kind))).toBe(glyphPathData(glyph(kind)));
    }
  });

  it('scales with the glyph scale', () => {
    const small = glyphPathData(glyph('mountain', { scale: 10 }));
    const large = glyphPathData(glyph('mountain', { scale: 20 }));
    expect(small).not.toBe(large);
  });

  it('keeps every coordinate near the glyph origin', () => {
    // A glyph must not stray more than ~2x its own scale from its anchor, or
    // thinning by centre distance would not prevent visible collisions.
    for (const kind of KINDS) {
      const g = glyph(kind);
      const nums = (glyphPathData(g).match(/-?\d+\.?\d*/g) ?? []).map(Number);
      for (let i = 0; i < nums.length; i += 2) {
        expect(Math.abs(nums[i] - g.x)).toBeLessThan(g.scale * 2);
        expect(Math.abs(nums[i + 1] - g.y)).toBeLessThan(g.scale * 2);
      }
    }
  });
});
```

- [ ] **Step 2: Run to verify it fails**

Run: `npx vitest run tests/glyphPaths.test.ts`
Expected: FAIL — cannot resolve `../utils/mapStyle/glyphPaths`.

- [ ] **Step 3: Implement `utils/mapStyle/glyphPaths.ts`**

```ts
import { GlyphKind, PlacedGlyph } from './types';

/**
 * Procedural glyph shapes as SVG path data, in output pixel coordinates.
 *
 * One path string serves both substrates: `SvgSubstrate` embeds it directly and
 * `Canvas2DSubstrate` feeds it to `new Path2D(d)`, which accepts SVG path
 * syntax. That is what guarantees the raster and vector maps draw the SAME
 * shape rather than two drifting approximations.
 *
 * Shapes are defined in a unit box — x and y in [-1, 1], baseline at y = 0,
 * apex toward -y — then rotated by `seedRot`, scaled and translated.
 */

type UnitPoint = [number, number];
type UnitSubpath = UnitPoint[];

const MOUNTAIN: UnitSubpath[] = [
  [[-1, 0], [-0.35, -1], [0, -0.45], [0.35, -1], [1, 0]],
  [[-0.35, -1], [-0.12, -0.62], [0.06, -0.78]], // snow/shading flick
];

const HILL: UnitSubpath[] = [
  [[-1, 0], [-0.55, -0.55], [0, -0.7], [0.55, -0.55], [1, 0]],
];

const FOREST: UnitSubpath[] = [
  [[-0.55, 0], [-0.55, -0.4]],
  [[-0.9, -0.4], [-0.55, -0.95], [-0.2, -0.4]],
  [[0.55, 0], [0.55, -0.35]],
  [[0.2, -0.35], [0.55, -0.85], [0.9, -0.35]],
];

const CONIFER: UnitSubpath[] = [
  [[0, 0], [0, -0.35]],
  [[-0.6, -0.3], [0, -1], [0.6, -0.3]],
  [[-0.42, -0.62], [0, -1], [0.42, -0.62]],
];

const DUNE: UnitSubpath[] = [
  [[-1, 0], [-0.4, -0.42], [0.35, -0.15], [1, -0.32]],
  [[-0.6, 0.28], [0.1, 0.05], [0.85, 0.22]],
];

const MARSH: UnitSubpath[] = [
  [[-0.8, 0], [0.8, 0]],
  [[-0.45, -0.1], [-0.45, -0.75]],
  [[0.05, -0.1], [0.05, -0.9]],
  [[0.55, -0.1], [0.55, -0.7]],
];

const SHAPES: Record<GlyphKind, UnitSubpath[]> = {
  mountain: MOUNTAIN,
  hill: HILL,
  forest: FOREST,
  conifer: CONIFER,
  dune: DUNE,
  marsh: MARSH,
};

const round = (n: number): string => (Math.round(n * 100) / 100).toString();

/**
 * SVG path `d` for a placed glyph, already positioned in output pixels.
 * Open polylines — these are drawn with `strokePath`, never filled.
 */
export const glyphPathData = (g: PlacedGlyph): string => {
  const cos = Math.cos(g.seedRot);
  const sin = Math.sin(g.seedRot);
  const half = g.scale / 2;

  const parts: string[] = [];
  for (const subpath of SHAPES[g.kind]) {
    const points = subpath.map(([ux, uy]) => {
      const sx = ux * half;
      const sy = uy * half;
      return [
        g.x + sx * cos - sy * sin,
        g.y + sx * sin + sy * cos,
      ] as UnitPoint;
    });
    parts.push(
      `M${round(points[0][0])} ${round(points[0][1])}` +
      points.slice(1).map(p => `L${round(p[0])} ${round(p[1])}`).join(''),
    );
  }
  return parts.join('');
};
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `npx vitest run tests/glyphPaths.test.ts`
Expected: PASS, 4 tests.

- [ ] **Step 5: Run the gates and commit**

Run: `npm run typecheck && npm run lint && npm test`

```bash
git add utils/mapStyle/glyphPaths.ts tests/glyphPaths.test.ts
git commit -m "feat(A3): procedural glyph paths for mountain/hill/forest/conifer/dune/marsh

Shapes are defined in a unit box, then rotated by the glyph's seeded angle,
scaled and translated to output pixels. One SVG path string serves both
substrates — SvgSubstrate embeds it, Canvas2DSubstrate feeds it to Path2D,
which accepts SVG path syntax. That is what keeps raster and vector drawing
the same shape instead of two drifting approximations."
```

---

### Task 5: `Substrate` interface and the Canvas2D implementation

**Files:**
- Create: `utils/mapStyle/substrate.ts` (the interface + shared spec types)
- Create: `utils/mapStyle/substrateCanvas.ts`
- Create: `tests/substrateCanvas.test.ts`

**Interfaces:**
- Consumes: `PlacedGlyph` (Task 2), `glyphPathData` (Task 4).
- Produces:
  ```ts
  interface Substrate {
    fillRect(x, y, w, h, fill: string): void;
    fillFeature(feature: GeoFeatureLike, fill: string, opacity?: number): void;
    strokeFeature(feature: GeoFeatureLike, stroke: string, width: number): void;
    strokeSegments(segments: Array<[Point3, Point3]>, stroke: string, width: number): void;
    hatchRect(x, y, w, h, spec: HatchSpec): void;
    hatchFeature(feature: GeoFeatureLike, spec: HatchSpec): void;
    grain(spec: GrainSpec): void;
    drawGlyph(g: PlacedGlyph, ink: string, width: number): void;
  }
  interface HatchSpec { color: string; spacingPx: number; widthPx: number; angleDeg: number }
  interface GrainSpec { seed: string; opacity: number; scale: number }
  ```
  Task 6 implements the same interface for SVG; Task 7's passes call only these.

- [ ] **Step 1: Write the failing test**

`tests/substrateCanvas.test.ts` uses a recording stub rather than a real canvas, so the test needs no DOM:

```ts
import { describe, it, expect, vi } from 'vitest';

import { Canvas2DSubstrate } from '../utils/mapStyle/substrateCanvas';
import { PlacedGlyph } from '../utils/mapStyle/types';

const makeCtx = () => ({
  save: vi.fn(), restore: vi.fn(), beginPath: vi.fn(), fill: vi.fn(),
  stroke: vi.fn(), fillRect: vi.fn(), clip: vi.fn(), rect: vi.fn(),
  moveTo: vi.fn(), lineTo: vi.fn(), translate: vi.fn(), rotate: vi.fn(),
  fillStyle: '', strokeStyle: '', lineWidth: 0, globalAlpha: 1,
  lineCap: '', lineJoin: '',
}) as unknown as CanvasRenderingContext2D;

const glyph: PlacedGlyph = { x: 10, y: 20, kind: 'mountain', scale: 16, seedRot: 0, cellId: 0 };

describe('Canvas2DSubstrate', () => {
  it('fills a rect with the given colour', () => {
    const ctx = makeCtx();
    new Canvas2DSubstrate(ctx, (() => {}) as never, 100, 50).fillRect(0, 0, 100, 50, '#abcdef');
    expect(ctx.fillRect).toHaveBeenCalledWith(0, 0, 100, 50);
    expect(ctx.fillStyle).toBe('#abcdef');
  });

  it('strokes a glyph without filling it', () => {
    const ctx = makeCtx();
    new Canvas2DSubstrate(ctx, (() => {}) as never, 100, 50).drawGlyph(glyph, '#3b2f1c', 1);
    expect(ctx.stroke).toHaveBeenCalled();
    expect(ctx.fill).not.toHaveBeenCalled();
  });

  it('restores the context for every save', () => {
    const ctx = makeCtx();
    const sub = new Canvas2DSubstrate(ctx, (() => {}) as never, 100, 50);
    sub.hatchRect(0, 0, 100, 50, { color: '#000', spacingPx: 8, widthPx: 1, angleDeg: 45 });
    sub.grain({ seed: 'x', opacity: 0.1, scale: 1 });
    expect((ctx.save as unknown as { mock: { calls: unknown[] } }).mock.calls.length)
      .toBe((ctx.restore as unknown as { mock: { calls: unknown[] } }).mock.calls.length);
  });
});
```

- [ ] **Step 2: Run to verify it fails**

Run: `npx vitest run tests/substrateCanvas.test.ts`
Expected: FAIL — cannot resolve `../utils/mapStyle/substrateCanvas`.

- [ ] **Step 3: Create `utils/mapStyle/substrate.ts`**

```ts
import { Point3 } from '../geo';
import { PlacedGlyph } from './types';

/** A d3 GeoJSON feature, typed structurally so this module needs no d3 import. */
export type GeoFeatureLike = { type: string; geometry: unknown };

export interface HatchSpec {
  color: string;
  spacingPx: number;
  widthPx: number;
  angleDeg: number;
}

export interface GrainSpec {
  seed: string;
  opacity: number;
  /** Feature size multiplier; 1 = fine paper tooth. */
  scale: number;
}

/**
 * The drawing surface a style pass writes to. Two implementations —
 * Canvas2D and SVG — so every parchment pass is authored ONCE.
 *
 * Deliberately narrow: passes may only do things both substrates can express
 * natively. Anything raster-only (a blur, a composite mode) would silently
 * degrade the SVG export, so it does not belong here.
 */
export interface Substrate {
  fillRect(x: number, y: number, w: number, h: number, fill: string): void;
  /**
   * `opacity` is an EXPLICIT parameter, never baked into `fill` as an `rgba()`
   * string. Canvas2D accepts `rgba()`; SVG 1.1 does not — it needs a separate
   * `fill-opacity`, and an `rgba()` fill renders inconsistently in Illustrator
   * and Inkscape. Colour strings stay opaque so both substrates agree.
   */
  fillFeature(feature: GeoFeatureLike, fill: string, opacity?: number): void;
  strokeFeature(feature: GeoFeatureLike, stroke: string, width: number): void;
  strokeSegments(segments: Array<[Point3, Point3]>, stroke: string, width: number): void;
  hatchRect(x: number, y: number, w: number, h: number, spec: HatchSpec): void;
  /** Hatch ONE feature. Full-bleed hatching would cover bare-paper land. */
  hatchFeature(feature: GeoFeatureLike, spec: HatchSpec): void;
  grain(spec: GrainSpec): void;
  drawGlyph(g: PlacedGlyph, ink: string, width: number): void;
}
```

- [ ] **Step 4: Implement `utils/mapStyle/substrateCanvas.ts`**

```ts
import { Point3 } from '../geo';
import { glyphPathData } from './glyphPaths';
import { GeoFeatureLike, GrainSpec, HatchSpec, Substrate } from './substrate';
import { PlacedGlyph } from './types';

type PathGeneratorLike = (object: unknown) => unknown;

/** Deterministic 0-1 noise for the grain pass — no RNG import, no shared state. */
const hash01 = (x: number, y: number, seed: string): number => {
  let h = 2166136261 >>> 0;
  const str = `${seed}:${x}:${y}`;
  for (let i = 0; i < str.length; i++) {
    h ^= str.charCodeAt(i);
    h = Math.imul(h, 16777619) >>> 0;
  }
  return (h >>> 8) / 0x1000000;
};

export class Canvas2DSubstrate implements Substrate {
  constructor(
    private ctx: CanvasRenderingContext2D,
    private path: PathGeneratorLike,
    private width: number,
    private height: number,
  ) {}

  fillRect(x: number, y: number, w: number, h: number, fill: string): void {
    this.ctx.fillStyle = fill;
    this.ctx.fillRect(x, y, w, h);
  }

  fillFeature(feature: GeoFeatureLike, fill: string, opacity = 1): void {
    this.ctx.save();
    this.ctx.globalAlpha = opacity;
    this.ctx.beginPath();
    this.path(feature);
    this.ctx.fillStyle = fill;
    this.ctx.fill();
    this.ctx.restore();
  }

  hatchFeature(feature: GeoFeatureLike, spec: HatchSpec): void {
    this.ctx.save();
    this.ctx.beginPath();
    this.path(feature);
    this.ctx.clip();
    this.hatchRect(0, 0, this.width, this.height, spec);
    this.ctx.restore();
  }

  strokeFeature(feature: GeoFeatureLike, stroke: string, width: number): void {
    this.ctx.beginPath();
    this.path(feature);
    this.ctx.strokeStyle = stroke;
    this.ctx.lineWidth = width;
    this.ctx.stroke();
  }

  strokeSegments(segments: Array<[Point3, Point3]>, stroke: string, width: number): void {
    if (!segments.length) return;
    this.ctx.strokeStyle = stroke;
    this.ctx.lineWidth = width;
    this.ctx.lineCap = 'round';
    this.ctx.lineJoin = 'round';
    this.ctx.beginPath();
    for (const [a, b] of segments) {
      this.path({ type: 'LineString', coordinates: [a, b] });
    }
    this.ctx.stroke();
  }

  hatchRect(x: number, y: number, w: number, h: number, spec: HatchSpec): void {
    this.ctx.save();
    this.ctx.beginPath();
    this.ctx.rect(x, y, w, h);
    this.ctx.clip();
    this.ctx.strokeStyle = spec.color;
    this.ctx.lineWidth = spec.widthPx;
    const rad = (spec.angleDeg * Math.PI) / 180;
    const dx = Math.cos(rad);
    const dy = Math.sin(rad);
    const span = Math.hypot(w, h);
    const steps = Math.ceil((span * 2) / Math.max(1, spec.spacingPx));
    this.ctx.beginPath();
    for (let i = -steps; i <= steps; i++) {
      const off = i * spec.spacingPx;
      const cx = x + w / 2 - dy * off;
      const cy = y + h / 2 + dx * off;
      this.ctx.moveTo(cx - dx * span, cy - dy * span);
      this.ctx.lineTo(cx + dx * span, cy + dy * span);
    }
    this.ctx.stroke();
    this.ctx.restore();
  }

  grain(spec: GrainSpec): void {
    this.ctx.save();
    this.ctx.globalAlpha = spec.opacity;
    const step = Math.max(2, Math.round(3 * spec.scale));
    for (let y = 0; y < this.height; y += step) {
      for (let x = 0; x < this.width; x += step) {
        const n = hash01(x, y, spec.seed);
        if (n < 0.55) continue;
        this.ctx.fillStyle = n > 0.85 ? '#ffffff' : '#000000';
        this.ctx.fillRect(x, y, step, step);
      }
    }
    this.ctx.restore();
  }

  drawGlyph(g: PlacedGlyph, ink: string, width: number): void {
    // Path2D accepts SVG path syntax, so glyphs use the SAME path string the
    // SVG substrate embeds. One shape definition, two surfaces.
    const p = new Path2D(glyphPathData(g));
    this.ctx.strokeStyle = ink;
    this.ctx.lineWidth = width;
    this.ctx.lineCap = 'round';
    this.ctx.lineJoin = 'round';
    this.ctx.stroke(p);
  }
}
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `npx vitest run tests/substrateCanvas.test.ts`
Expected: PASS, 3 tests.

If `Path2D is not defined`, the test environment lacks it — stub it in the test file with `vi.stubGlobal('Path2D', class { constructor(_d: string) {} })` before the `drawGlyph` test.

- [ ] **Step 6: Run the gates and commit**

Run: `npm run typecheck && npm run lint && npm test`

```bash
git add utils/mapStyle/substrate.ts utils/mapStyle/substrateCanvas.ts tests/substrateCanvas.test.ts
git commit -m "feat(A3): Substrate interface and Canvas2D implementation

The interface is deliberately narrow — passes may only do what BOTH substrates
express natively. A raster-only primitive would silently degrade the SVG export.

Glyphs go through Path2D with the same SVG path string the vector substrate
embeds, so both surfaces draw one shape definition."
```

---

### Task 6: SVG substrate

**Files:**
- Create: `utils/mapStyle/substrateSvg.ts`
- Create: `tests/substrateSvg.test.ts`

**Interfaces:**
- Consumes: `Substrate`, `HatchSpec`, `GrainSpec` (Task 5), `glyphPathData` (Task 4).
- Produces: `SvgSubstrate` with `defs(): string` and `body(): string`, both consumed by Task 9.

- [ ] **Step 1: Write the failing test**

`tests/substrateSvg.test.ts`:

```ts
import { describe, it, expect } from 'vitest';

import { SvgSubstrate } from '../utils/mapStyle/substrateSvg';
import { PlacedGlyph } from '../utils/mapStyle/types';

const glyph: PlacedGlyph = { x: 10, y: 20, kind: 'mountain', scale: 16, seedRot: 0, cellId: 0 };

describe('SvgSubstrate', () => {
  it('emits a rect for fillRect', () => {
    const s = new SvgSubstrate((() => '') as never, 100, 50);
    s.fillRect(0, 0, 100, 50, '#abcdef');
    expect(s.body()).toContain('<rect');
    expect(s.body()).toContain('#abcdef');
  });

  it('registers a pattern in defs for hatching, and references it', () => {
    const s = new SvgSubstrate((() => '') as never, 100, 50);
    s.hatchRect(0, 0, 100, 50, { color: '#000', spacingPx: 8, widthPx: 1, angleDeg: 45 });
    expect(s.defs()).toContain('<pattern');
    expect(s.body()).toMatch(/fill="url\(#/);
  });

  it('registers a turbulence filter for grain', () => {
    const s = new SvgSubstrate((() => '') as never, 100, 50);
    s.grain({ seed: 'x', opacity: 0.12, scale: 1 });
    expect(s.defs()).toContain('feTurbulence');
  });

  it('emits a stroked, unfilled path for a glyph', () => {
    const s = new SvgSubstrate((() => '') as never, 100, 50);
    s.drawGlyph(glyph, '#3b2f1c', 1);
    expect(s.body()).toContain('fill="none"');
    expect(s.body()).toContain('stroke="#3b2f1c"');
  });

  it('emits fill-opacity rather than an rgba() fill', () => {
    const s = new SvgSubstrate((() => 'M0 0L1 1') as never, 100, 50);
    s.fillFeature({ type: 'Feature', geometry: {} }, '#000000', 0.25);
    expect(s.body()).toContain('fill-opacity="0.250"');
    expect(s.body()).not.toContain('rgba(');
  });

  it('reuses one pattern id for identical hatch specs', () => {
    const s = new SvgSubstrate((() => '') as never, 100, 50);
    const spec = { color: '#000', spacingPx: 8, widthPx: 1, angleDeg: 45 };
    s.hatchRect(0, 0, 50, 50, spec);
    s.hatchRect(50, 0, 50, 50, spec);
    expect(s.defs().match(/<pattern/g)?.length).toBe(1);
  });
});
```

- [ ] **Step 2: Run to verify it fails**

Run: `npx vitest run tests/substrateSvg.test.ts`
Expected: FAIL — cannot resolve `../utils/mapStyle/substrateSvg`.

- [ ] **Step 3: Implement `utils/mapStyle/substrateSvg.ts`**

```ts
import { Point3 } from '../geo';
import { glyphPathData } from './glyphPaths';
import { GeoFeatureLike, GrainSpec, HatchSpec, Substrate } from './substrate';
import { PlacedGlyph } from './types';

type PathStringGenerator = (object: unknown) => string | null;

const esc = (s: string): string => s.replace(/[<>&"]/g, c =>
  ({ '<': '&lt;', '>': '&gt;', '&': '&amp;', '"': '&quot;' }[c] as string));

/**
 * Accumulates SVG. `defs()` and `body()` are spliced into the document by the
 * export path.
 *
 * Nothing here is a fallback or an approximation of the raster look: SVG has
 * native equivalents for both parchment ingredients that seem raster-only —
 * hatching is `<pattern>` and paper grain is `<filter><feTurbulence>`.
 */
export class SvgSubstrate implements Substrate {
  private defsParts: string[] = [];
  private bodyParts: string[] = [];
  private patternIds = new Map<string, string>();
  private grainId: string | null = null;

  constructor(
    private path: PathStringGenerator,
    private width: number,
    private height: number,
  ) {}

  defs(): string {
    return this.defsParts.length ? `<defs>${this.defsParts.join('')}</defs>` : '';
  }

  body(): string {
    return this.bodyParts.join('');
  }

  fillRect(x: number, y: number, w: number, h: number, fill: string): void {
    this.bodyParts.push(`<rect x="${x}" y="${y}" width="${w}" height="${h}" fill="${esc(fill)}"/>`);
  }

  fillFeature(feature: GeoFeatureLike, fill: string, opacity = 1): void {
    const d = this.path(feature);
    if (!d) return;
    // fill-opacity, NOT an rgba() fill — SVG 1.1 has no rgba() colour syntax.
    const op = opacity < 1 ? ` fill-opacity="${opacity.toFixed(3)}"` : '';
    this.bodyParts.push(
      `<path d="${d}" fill="${esc(fill)}"${op} stroke="${esc(fill)}" stroke-width="0.5"/>`,
    );
  }

  hatchFeature(feature: GeoFeatureLike, spec: HatchSpec): void {
    const d = this.path(feature);
    if (!d) return;
    const id = this.hatchPatternId(spec);
    this.bodyParts.push(`<path d="${d}" fill="url(#${id})"/>`);
  }

  strokeFeature(feature: GeoFeatureLike, stroke: string, width: number): void {
    const d = this.path(feature);
    if (!d) return;
    this.bodyParts.push(`<path d="${d}" fill="none" stroke="${esc(stroke)}" stroke-width="${width}"/>`);
  }

  strokeSegments(segments: Array<[Point3, Point3]>, stroke: string, width: number): void {
    const ds: string[] = [];
    for (const [a, b] of segments) {
      const d = this.path({ type: 'LineString', coordinates: [a, b] });
      if (d) ds.push(d);
    }
    if (!ds.length) return;
    this.bodyParts.push(
      `<path d="${ds.join('')}" fill="none" stroke="${esc(stroke)}" stroke-width="${width}" ` +
      `stroke-linecap="round" stroke-linejoin="round"/>`,
    );
  }

  hatchRect(x: number, y: number, w: number, h: number, spec: HatchSpec): void {
    // One <pattern> per distinct spec — repeating a pattern definition per
    // rect would bloat the file for no visual difference.
    const id = this.hatchPatternId(spec);
    this.bodyParts.push(`<rect x="${x}" y="${y}" width="${w}" height="${h}" fill="url(#${id})"/>`);
  }

  /** One <pattern> per distinct spec, shared by hatchRect and hatchFeature. */
  private hatchPatternId(spec: HatchSpec): string {
    const key = `${spec.color}|${spec.spacingPx}|${spec.widthPx}|${spec.angleDeg}`;
    let id = this.patternIds.get(key);
    if (!id) {
      id = `hatch${this.patternIds.size}`;
      this.patternIds.set(key, id);
      this.defsParts.push(
        `<pattern id="${id}" width="${spec.spacingPx}" height="${spec.spacingPx}" ` +
        `patternUnits="userSpaceOnUse" patternTransform="rotate(${spec.angleDeg})">` +
        `<line x1="0" y1="0" x2="0" y2="${spec.spacingPx}" ` +
        `stroke="${esc(spec.color)}" stroke-width="${spec.widthPx}"/>` +
        `</pattern>`,
      );
    }
    return id;
  }

  grain(spec: GrainSpec): void {
    if (!this.grainId) {
      this.grainId = 'paperGrain';
      this.defsParts.push(
        `<filter id="${this.grainId}" x="0" y="0" width="100%" height="100%">` +
        `<feTurbulence type="fractalNoise" baseFrequency="${0.8 / Math.max(0.1, spec.scale)}" ` +
        `numOctaves="4" seed="${Math.abs(hashSeed(spec.seed)) % 1000}" result="n"/>` +
        `<feColorMatrix type="saturate" values="0" in="n" result="g"/>` +
        `</filter>`,
      );
    }
    this.bodyParts.push(
      `<rect x="0" y="0" width="${this.width}" height="${this.height}" ` +
      `filter="url(#${this.grainId})" opacity="${spec.opacity}"/>`,
    );
  }

  drawGlyph(g: PlacedGlyph, ink: string, width: number): void {
    this.bodyParts.push(
      `<path d="${glyphPathData(g)}" fill="none" stroke="${esc(ink)}" stroke-width="${width}" ` +
      `stroke-linecap="round" stroke-linejoin="round"/>`,
    );
  }
}

const hashSeed = (seed: string): number => {
  let h = 2166136261 >>> 0;
  for (let i = 0; i < seed.length; i++) {
    h ^= seed.charCodeAt(i);
    h = Math.imul(h, 16777619) >>> 0;
  }
  return h | 0;
};
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `npx vitest run tests/substrateSvg.test.ts`
Expected: PASS, 5 tests.

- [ ] **Step 5: Run the gates and commit**

Run: `npm run typecheck && npm run lint && npm test`

```bash
git add utils/mapStyle/substrateSvg.ts tests/substrateSvg.test.ts
git commit -m "feat(A3): SVG substrate with native pattern hatching and turbulence grain

Neither parchment ingredient needs a raster fallback: hatching is <pattern>,
paper grain is <filter><feTurbulence>. Patterns are deduplicated by spec so
repeated hatch fills share one definition."
```

---

### Task 7: Parchment passes

**Files:**
- Modify: `utils/mapStyle/types.ts` (add `StylePass`, `StyleRenderContext` and `passes` to `MapStyle`)
- Modify: `utils/mapStyle/styleDefault.ts` (empty pass list)
- Modify: `utils/mapStyle/styleParchment.ts` (the real pass list)
- Create: `utils/mapStyle/passes.ts`
- Modify: `tests/mapStyle.test.ts` (assert the default style draws nothing)

**Interfaces:**
- Consumes: everything from Tasks 2–6.
- Produces:
  ```ts
  interface StyleRenderContext {
    world: WorldData; viewMode: ViewMode; widthPx: number; heightPx: number;
    projection: d3.GeoProjection; glyphs: PlacedGlyph[];
    shadeMap: Float32Array | null; colorCtx: ColorContext;
    coastlines: Array<[Point3, Point3]>;
  }
  type StylePass = (ctx: StyleRenderContext, sub: Substrate) => void;
  runStyle(style: MapStyle, ctx: StyleRenderContext, sub: Substrate): void
  ```
  Tasks 8 and 9 call `runStyle`.

- [ ] **Step 1: Add the pass types to `utils/mapStyle/types.ts`**

Append:

```ts
export interface StyleRenderContext {
  world: WorldData;
  viewMode: ViewMode;
  widthPx: number;
  heightPx: number;
  projection: d3.GeoProjection;
  /** Pre-computed by placeGlyphs; empty when the fill policy suppresses glyphs. */
  glyphs: PlacedGlyph[];
  shadeMap: Float32Array | null;
  colorCtx: ColorContext;
  coastlines: Array<[Point3, Point3]>;
}

export type StylePass = (ctx: StyleRenderContext, sub: Substrate) => void;
```

Add `passes: StylePass[];` to `MapStyle`. Import `WorldData`, `ViewMode` from `../../types`, `ColorContext` from `../colors`, `Substrate` from `./substrate`, and `d3` from `d3-geo`.

- [ ] **Step 2: Add `passes: []` to `styleDefault.ts`**

```ts
  // Empty by design: the default style is the pre-A3 look, and every view mode
  // already paints its own fill. Selecting it must change nothing.
  passes: [],
```

- [ ] **Step 3: Add the assertion to `tests/mapStyle.test.ts`**

```ts
  it('draws nothing in the default style', () => {
    expect(getMapStyle('default').passes).toEqual([]);
  });
```

- [ ] **Step 4: Implement `utils/mapStyle/passes.ts`**

```ts
import * as THREE from 'three';

import { ViewMode } from '../../types';
import { getCellColor } from '../colors';
import { seasonalTemperatureDelta } from '../seasons';
import { computeShadeMap } from '../shading';
import { isLakeCell } from '../worldGen';
import { Substrate } from './substrate';
import { FillPolicy, MapStyle, StylePalette, StylePass, StyleRenderContext } from './types';

/** Run a style's passes in order. The only entry point render paths need. */
export const runStyle = (style: MapStyle, ctx: StyleRenderContext, sub: Substrate): void => {
  for (const pass of style.passes) pass(ctx, sub);
};

/** Full-bleed paper tone plus grain. Always first — everything sits on it. */
export const paperPass = (palette: StylePalette, seed: string): StylePass =>
  (ctx, sub) => {
    sub.fillRect(0, 0, ctx.widthPx, ctx.heightPx, palette.paper);
    sub.grain({ seed, opacity: 0.10, scale: 1 });
  };

/** Flat sea tone on ocean cells only. */
export const oceanFillPass = (palette: StylePalette): StylePass =>
  (ctx, sub) => {
    const { world, colorCtx } = ctx;
    for (let i = 0; i < world.cells.length; i++) {
      const cell = world.cells[i];
      if (cell.height >= colorCtx.seaLevel) continue;
      const feature = world.geoJson?.features?.[i];
      if (!feature) continue;
      sub.fillFeature(feature, palette.sea);
    }
  };

/**
 * Diagonal hatch on OCEAN CELLS ONLY, never full-bleed.
 *
 * A full-bleed hatchRect would cover bare-paper land: under the 'bare' fill
 * policy — satellite, biome, height_bw, the flagship modes — landPass returns
 * early and paints nothing, so a full-bleed hatch would sit directly on the
 * land. Per-feature hatching is why `hatchFeature` exists on Substrate.
 *
 * Runs AFTER hillshadePass so the hatching sits over the relief shading, per
 * spec §7.4.
 */
export const oceanHatchPass = (palette: StylePalette): StylePass =>
  (ctx, sub) => {
    const { world, colorCtx } = ctx;
    for (let i = 0; i < world.cells.length; i++) {
      const cell = world.cells[i];
      if (cell.height >= colorCtx.seaLevel) continue;
      const feature = world.geoJson?.features?.[i];
      if (!feature) continue;
      sub.hatchFeature(feature, {
        color: palette.seaHatch, spacingPx: 7, widthPx: 0.6, angleDeg: 45,
      });
    }
  };

/**
 * Land fill, obeying the style's fill policy for this view mode (spec §4):
 *   bare        → paper shows through; glyphs carry the terrain
 *   categorical → the mode's own fill, muted toward the paper
 *   ramp        → the mode's own fill at full strength
 *
 * Takes `fillPolicy` and `palette` directly rather than the whole MapStyle:
 * a style holds its passes, so passing the style in would be circular and
 * force a cast.
 */
export const landPass = (
  fillPolicy: (mode: ViewMode) => FillPolicy,
  palette: StylePalette,
): StylePass =>
  (ctx, sub) => {
    const policy = fillPolicy(ctx.viewMode);
    if (policy === 'bare') return;
    const mute = policy === 'categorical' ? 0.45 : 0;
    const paper = new THREE.Color(palette.paper);
    const { world, colorCtx } = ctx;
    for (let i = 0; i < world.cells.length; i++) {
      const cell = world.cells[i];
      if (cell.height < colorCtx.seaLevel && !isLakeCell(cell)) continue;
      const feature = world.geoJson?.features?.[i];
      if (!feature) continue;
      // seasonalDelta is PER CELL — a shared ColorContext would silently kill
      // the D1 seasonal biome re-derivation and D3 sea ice under parchment.
      const color = getCellColor(cell, ctx.viewMode, {
        ...colorCtx,
        seasonalDelta: seasonalTemperatureDelta(cell, world.params),
      });
      if (mute > 0) color.lerp(paper, mute);
      sub.fillFeature(feature, `#${color.getHexString()}`);
    }
  };

/**
 * Existing hillshade, multiplied in softly so relief reads under the glyphs.
 * Opacity is an explicit argument, never an rgba() colour — see Substrate.
 */
export const hillshadePass = (opacityScale: number): StylePass =>
  (ctx, sub) => {
    if (!ctx.shadeMap) return;
    const { world } = ctx;
    for (let i = 0; i < world.cells.length; i++) {
      const s = ctx.shadeMap[world.cells[i].id];
      if (s === 1) continue; // flat ground contributes nothing
      const feature = world.geoJson?.features?.[i];
      if (!feature) continue;
      const a = Math.min(0.5, Math.abs(1 - s) * opacityScale);
      sub.fillFeature(feature, s < 1 ? '#000000' : '#ffffff', a);
    }
  };

/** Heavy ink coastline plus a lighter offset swash line just outside it. */
export const coastlinePass = (palette: StylePalette, widthPx: number): StylePass =>
  (ctx, sub) => {
    const w = Math.max(0.75, (widthPx / 1024) * 1.4);
    sub.strokeSegments(ctx.coastlines, palette.coast, w);
    sub.strokeSegments(ctx.coastlines, palette.inkLight, w * 0.5);
  };

/** Relief and vegetation glyphs. Empty when the fill policy suppressed them. */
export const glyphPass = (palette: StylePalette, widthPx: number): StylePass =>
  (ctx, sub) => {
    const w = Math.max(0.6, (widthPx / 1024) * 1.1);
    for (const g of ctx.glyphs) sub.drawGlyph(g, palette.ink, w);
  };

export { computeShadeMap };
```

- [ ] **Step 5: Wire the pass list into `styleParchment.ts`**

Replace the exported object's tail:

```ts
export const styleParchment: MapStyle = {
  id: 'parchment',
  name: 'Parchment',
  palette: PARCHMENT_PALETTE,
  fillPolicy: parchmentFillPolicy,
  passes: [
    paperPass(PARCHMENT_PALETTE, 'parchment'),
    oceanFillPass(PARCHMENT_PALETTE),
    landPass(parchmentFillPolicy, PARCHMENT_PALETTE),
    hillshadePass(0.5),
    // Hatch AFTER hillshade (spec §7.4: shading sits under the hatching), and
    // per-ocean-cell, so bare-paper land is never covered.
    oceanHatchPass(PARCHMENT_PALETTE),
    coastlinePass(PARCHMENT_PALETTE, 1024),
    glyphPass(PARCHMENT_PALETTE, 1024),
  ],
};
```

Extract the palette to a named `PARCHMENT_PALETTE` const above so both the style and the passes reference one object.

- [ ] **Step 6: Run the gates and commit**

Run: `npm run typecheck && npm run lint && npm test`
Expected: all pass; `tests/mapStyle.test.ts` gains 1 test.

```bash
git add utils/mapStyle tests/mapStyle.test.ts
git commit -m "feat(A3): parchment passes — paper, ocean hatch, land, hillshade, coast, glyphs

Passes are written once against Substrate, so Map2D, PNG and SVG get the same
drawing. landPass obeys the per-view-mode fill policy from spec section 4:
bare paper for glyph-carried modes, muted fill for categorical modes, and the
mode's own fill untouched for continuous ramps."
```

---

### Task 8: Wire parchment into Map2D

**Files:**
- Modify: `components/Map2D.tsx` (props, the main screen render at ~`:598`, the Dymaxion source render at ~`:362`)
- Modify: `components/shell/ShellApp.tsx` (pass `mapStyleId` to `Map2D`)

**Interfaces:**
- Consumes: `getMapStyle`, `runStyle`, `placeGlyphs`, `Canvas2DSubstrate`.
- Produces: nothing new.

- [ ] **Step 1: Add the prop**

Add `mapStyleId?: MapStyleId;` to the `Map2D` props interface, defaulting to `'default'` in the destructure at `components/Map2D.tsx:230`. Pass it from `ShellApp.tsx` at both `Map2D` render sites.

- [ ] **Step 2: Compute glyphs once, memoised**

In `Map2D`, beside the existing `shadeMap` memo at `:260`:

```ts
  const style = useMemo(() => getMapStyle(mapStyleId), [mapStyleId]);

  // Placement is projection- and width-dependent, so it re-runs when either
  // changes — but NOT on every frame.
  const glyphs = useMemo(() => {
    if (!world || style.id === 'default') return [];
    if (style.fillPolicy(viewMode) === 'ramp') return []; // ramps: glyphs would fight the fill
    return placeGlyphs(world.cells, projection, size.width, {
      seaLevel: world.params.seaLevel,
      seed: world.params.seed,
    });
  }, [world, style, viewMode, projection, size.width]);
```

- [ ] **Step 3: Branch the main screen render**

At `components/Map2D.tsx:598`, after `const pathGenerator = d3.geoPath(projection, ctx);`, insert:

```ts
    if (style.id !== 'default') {
      const sub = new Canvas2DSubstrate(ctx, pathGenerator as unknown as (o: unknown) => unknown, size.width, size.height);
      runStyle(style, {
        world, viewMode, widthPx: size.width, heightPx: size.height, projection,
        glyphs, shadeMap, coastlines: coastlineSegments,
        colorCtx: {
          seaLevel: world.params.seaLevel, factionColors, cultureColors, religionColors,
        },
      }, sub);
    } else {
      // ...existing per-cell fill loop, unchanged...
    }
```

`coastlineSegments` comes from the same computation `exportVector.ts:62` uses — extract `computeCoastlineSegments` from `utils/exportVector.ts` into `utils/shading.ts` (which already owns `sharedEdge`) and import it in both places, rather than duplicating it.

- [ ] **Step 4: Do the same for the Dymaxion source buffer**

Repeat at `components/Map2D.tsx:362`. **Critical:** that path wraps drawing in `srcCtx.translate(srcWidth, 0); srcCtx.scale(-1, 1);` (`:367`). Glyphs and any text drawn inside that transform come out **mirrored**. Run the style passes for fills and strokes inside the mirror, then `srcCtx.restore()` and draw the glyph pass **outside** it, placing glyphs with an x-flipped projection.

Simplest correct approach: after the mirrored passes, restore, then draw glyphs at `x' = srcWidth - x`.

- [ ] **Step 5: Verify in the browser**

Start the dev server only if one is not already running on port 3000 — **never kill a preexisting server**.

Check, on a generated world:
1. Style `Default` looks **identical** to before this task.
2. Style `Parchment` in `biome` mode: paper ground, hatched sea, ink coastline, mountain and forest glyphs.
3. `political` mode: muted faction fills on paper.
4. `temperature` mode: the ramp still renders, and **no glyphs**.
5. Switch projection to **Dymaxion** — glyphs must not be mirrored.

- [ ] **Step 6: Run the gates and commit**

Run: `npm run typecheck && npm run lint && npm test`

```bash
git add components/Map2D.tsx components/shell/ShellApp.tsx utils/shading.ts utils/exportVector.ts
git commit -m "feat(A3): render parchment in Map2D

Glyph placement is memoised on [world, style, viewMode, projection, width] so it
does not re-run per frame. Ramp view modes get an empty glyph list — glyphs
would fight a continuous fill.

The Dymaxion source buffer mirrors its canvas, so glyphs draw outside that
transform with an x-flipped coordinate; inside it they came out backwards.

computeCoastlineSegments moved from exportVector to shading so Map2D and the
vector export share one implementation."
```

---

### Task 9: Wire parchment into PNG and SVG export

**Files:**
- Modify: `utils/export.ts:105-120` (the raster cell loop)
- Modify: `utils/exportVector.ts:240` (`exportSVG`)
- Modify: `components/Controls.tsx` (export panel passes the current style)

**Interfaces:**
- Consumes: everything above.
- Produces: nothing new.

- [ ] **Step 1: Thread `mapStyleId` into both export entry points**

Add an optional `mapStyleId: MapStyleId = 'default'` parameter to the PNG export function containing `utils/export.ts:105` and to `exportSVG` / `downloadSVG` in `utils/exportVector.ts`. Pass the live value from the export panel in `Controls.tsx`.

- [ ] **Step 2: Branch the raster export**

In `utils/export.ts`, replace the cell loop at `:107-120` with the same `if (style.id !== 'default')` branch used in Task 8, constructing a `Canvas2DSubstrate` over the export canvas and its `pathGenerator`, and calling `placeGlyphs(world.cells, projection, width, …)` with the **export** width so glyph density matches the output size.

- [ ] **Step 3: Branch the vector export**

In `exportSVG`, when the style is not `default`, build an `SvgSubstrate(pathGenerator as PathStringGenerator, width, height)`, run the passes, and assemble:

```ts
    return (
      `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">` +
      sub.defs() +
      `<g id="style" transform="${mirror}">${sub.body()}</g>` +
      labels +
      `</svg>`
    );
```

**Glyphs must not go inside the mirrored group** — the same flip as Task 8. Emit them into a second, unmirrored group with x-flipped coordinates, exactly as Map2D does.

- [ ] **Step 4: Verify by inspection**

1. Export a **PNG** at 2048 and at 8192 in parchment. They must look like the same map at two resolutions — glyph density roughly equal, not denser at 8192.
2. Export an **SVG** in parchment. Open it in a browser: hatching, grain, coastline and glyphs all present.
3. Open the SVG in a text editor: confirm exactly one `<pattern>` per distinct hatch and one `feTurbulence` filter.
4. Export in `Default` style: byte-comparable to a pre-A3 export.

- [ ] **Step 5: Run the gates and commit**

Run: `npm run typecheck && npm run lint && npm test`

```bash
git add utils/export.ts utils/exportVector.ts components/Controls.tsx
git commit -m "feat(A3): parchment in PNG and SVG export

Both exports run the same pass list Map2D does, so the exported map matches the
screen. placeGlyphs receives the export width, so an 8192px PNG is the same map
at higher resolution rather than a denser one.

Glyphs render outside the mirror transform in SVG, matching Map2D's Dymaxion
handling."
```

---

### Task 10: Documentation

**Files:**
- Create: `docs/map-styles.md`
- Modify: `docs/README.md` (index entry)
- Modify: `ROADMAP.md` (A3 → DONE)
- Modify: `HANDOFF.md` (session entry)

- [ ] **Step 1: Write `docs/map-styles.md`**

Cover: the substrate adapter and why it exists, the pass list, the fill-policy table from spec §4, how to add a style (implement `MapStyle`, register in `index.ts`), the glyph placement contract, and the mirror gotcha in both Map2D-Dymaxion and SVG export.

- [ ] **Step 2: Update `ROADMAP.md`**

Change the A3 heading to `### A3. Map style system (fantasy rendering)  —  ✅ DONE`, summarise what shipped, and copy the preset backlog from spec §9 so the unbuilt ideas live in the roadmap rather than only in a spec.

- [ ] **Step 3: Update `HANDOFF.md`**

Add the session entry: what shipped, the fill-policy decision and why ramp modes were exempted, the mirror gotcha, and that the globe is deliberately excluded.

- [ ] **Step 4: Run the gates and commit**

```bash
git add docs/map-styles.md docs/README.md ROADMAP.md HANDOFF.md
git commit -m "docs(A3): map style system documentation, roadmap and handoff"
```

---

## Self-Review

**Spec coverage:**

| Spec section | Task |
|---|---|
| §3.1 Substrate adapter | 5 (interface + Canvas2D), 6 (SVG) |
| §3.2 Pass list | 7 |
| §3.3 Glyph placement | 3 |
| §3.4 File layout | 2–7 |
| §4 ViewMode fill rule | 2 (`fillPolicy`), 7 (`landPass` honours it), 8 (glyph suppression) |
| §5 `ColorContext` refactor | 1 |
| §5b State location and UI | 2 |
| §6 Data flow | 8, 9 |
| §7.1 Mirror gotcha | 8 step 4, 9 step 3 |
| §7.2 Glyph inside its own cell | 3 (placement uses `cell.center` through the same projection) |
| §7.3 Antimeridian | 3 (pixel-bounds clip) |
| §7.4 D10 water hillshade | 7 (`hillshadePass` skips `s === 1`, so flat water contributes nothing) |
| §8 Testing | 2, 3, 4, 5, 6 unit tests; 8, 9 manual browser checks |
| §9 Preset backlog | 10 step 2 |

**Gap found and closed:** §7.2 originally had no task. Placement in Task 3 uses `cell.center` projected through the *same* `d3.GeoProjection` that draws the cell polygon, so a glyph lands inside its own cell by construction. Noted in the table rather than adding a task.

**Placeholder scan:** clean. Task 9 steps 1–3 describe edits against exact line numbers rather than quoting whole functions, because those functions are long and unchanged apart from one branch — the branch code itself is given in Task 8 step 3 and referenced.

**Type consistency:** `PlacedGlyph`, `GlyphKind`, `FillPolicy`, `MapStyleId`, `StylePass`, `StyleRenderContext`, `Substrate`, `HatchSpec`, `GrainSpec`, `ColorContext` are each defined in exactly one task and used with the same shape thereafter. `placeGlyphs` keeps the signature `(cells, projection, widthPx, opts)` in Tasks 3, 8 and 9. `getMapStyle` / `MAP_STYLES` are defined in Task 2 and used unchanged in 7, 8, 9.

**Corrections applied after review (2026-08-28):**

1. **`oceanPass` split into `oceanFillPass` + `oceanHatchPass`, hatching per ocean
   cell.** The original ended with a full-bleed `hatchRect`, which under the
   `bare` fill policy — `satellite`, `biome`, `height_bw`, the flagship modes —
   would have covered the bare-paper land in sea hatching, because `landPass`
   returns early and paints nothing there. `hatchFeature` was added to
   `Substrate` for this.
2. **`landPass` now spreads `seasonalDelta` per cell.** `ColorContext.seasonalDelta`
   is per-cell — every pre-A3 call site computes it inside its loop. Reusing one
   shared context would have silently disabled D1 seasonal biome re-derivation and
   D3 sea ice under parchment, with no test to catch it.
3. **`fillFeature` takes an explicit `opacity`.** `hillshadePass` previously passed
   an `rgba()` string, which Canvas2D accepts but SVG 1.1 has no syntax for — it
   renders inconsistently in Illustrator and Inkscape. That is exactly the raster
   /vector drift the substrate exists to prevent, so opacity is now a parameter
   and colour strings stay opaque.
4. **Pass order now matches spec §7.4** — hillshade under the ocean hatching.
5. **`landPass` takes `fillPolicy` and `palette` instead of the whole `MapStyle`,**
   removing an `as MapStyle` cast that papered over a circular reference (a style
   holds its passes, so a pass cannot take the style).
