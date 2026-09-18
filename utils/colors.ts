import { BiomeType, Cell, CivData, CultureData, ReligionData, ViewMode, WorldData } from '../types';
import * as THREE from 'three';
import { FACTION_COLORS, CULTURE_COLORS, RELIGION_COLORS } from './palette';
import { determineBiome } from './worldGen';
import { SEAWATER_FREEZE_C } from './seasons';
import { formatElevation } from './datum';

export { FACTION_COLORS, CULTURE_COLORS, RELIGION_COLORS };

// Temperature-view ramp domain, in °C. Module-scoped so the fill (getCellColor)
// and the legend (legendFor) read ONE definition — a shared literal is drift by
// construction, which is the whole reason F5 exists.
export const TEMP_MIN_C = -30;
export const TEMP_MAX_C = 50;

// Builds the live faction-color map from civData so every render/export path
// (viewer, minimap, PNG, GLB) reflects user-edited faction colors identically.
export const buildFactionColorMap = (civData?: CivData): Map<number, string> | undefined =>
  civData ? new Map(civData.factions.map(f => [f.id, f.color])) : undefined;

// Same pattern as buildFactionColorMap, for the culture layer (C1).
export const buildCultureColorMap = (cultures?: CultureData[]): Map<number, string> | undefined =>
  cultures ? new Map(cultures.map(c => [c.id, c.color])) : undefined;

// Same pattern again, for the religion layer (C2).
export const buildReligionColorMap = (religions?: ReligionData[]): Map<number, string> | undefined =>
  religions ? new Map(religions.map(r => [r.id, r.color])) : undefined;

// Earth-like Natural Colors
export const BIOME_COLORS: Record<BiomeType, string> = {
  [BiomeType.DEEP_OCEAN]: '#1a237e', // Deep Blue
  [BiomeType.OCEAN]: '#0277bd',      // Standard Blue
  
  // Polar
  [BiomeType.ICE_CAP]: '#ffffff',    // White
  [BiomeType.TUNDRA]: '#cfd8dc',     // Greyish Cyan/White
  
  // Dry
  [BiomeType.HOT_DESERT]: '#e6c27e', // Sandy Orange
  [BiomeType.COLD_DESERT]: '#bcaaa4', // Greyish Brown
  [BiomeType.STEPPE]: '#c5e1a5',     // Pale dry green
  
  // Tropical
  [BiomeType.TROPICAL_RAINFOREST]: '#004d40', // Deep Jungle Green
  [BiomeType.TROPICAL_SAVANNA]: '#aed581',    // Yellowish Green
  
  // Temperate
  [BiomeType.GRASSLAND]: '#9ccc65',           // Meadow green — temperate prairie/plains
  [BiomeType.MEDITERRANEAN]: '#8d6e63',       // Dry brownish green
  [BiomeType.TEMPERATE_FOREST]: '#2e7d32',    // Standard Forest Green
  [BiomeType.TEMPERATE_RAINFOREST]: '#1b5e20', // Darker Green
  
  // Continental
  [BiomeType.BOREAL_FOREST]: '#00695c',       // Pine Green (Blueish)
  
  // Special
  [BiomeType.BEACH]: '#fff59d',      // Sand
  [BiomeType.VOLCANIC]: '#37474f',   // Dark Grey Rock

  // Hydrology
  [BiomeType.LAKE]: '#3aa0cf',       // Fresh mid blue — lighter/greener than ocean
  [BiomeType.SALT_LAKE]: '#cfe8e4',  // Pale turquoise-white — evaporite basin
};

export const PLATE_COLORS = [
  '#ef5350', '#ab47bc', '#7e57c2', '#5c6bc0', '#42a5f5', '#29b6f6',
  '#26c6da', '#26a69a', '#66bb6a', '#9ccc65', '#d4e157', '#ffee58',
  '#ffca28', '#ffa726', '#ff7043', '#8d6e63', '#bdbdbd', '#78909c'
];

// Folk-religion color: the parent culture's color, darkened and desaturated
// ~30% so folk faith reads as a muted variant of its culture rather than
// competing visually with the vivid RELIGION_COLORS palette.
export const darkenForFolk = (hex: string): string => {
  const c = new THREE.Color(hex);
  const hsl = { h: 0, s: 0, l: 0 };
  c.getHSL(hsl);
  c.setHSL(hsl.h, hsl.s * 0.7, hsl.l * 0.7);
  return '#' + c.getHexString();
};

const getProvinceVariant = (baseColorHex: string, provId: number, strength = 1): THREE.Color => {
  const c = new THREE.Color(baseColorHex);
  const r = Math.sin(provId * 12.9898) * 43758.5453;
  const rnd = r - Math.floor(r);
  const r2 = Math.cos(provId * 78.233) * 43758.5453;
  const rnd2 = r2 - Math.floor(r2);

  const hsl = { h: 0, s: 0, l: 0 };
  c.getHSL(hsl);

  hsl.l = Math.max(0.1, Math.min(0.9, hsl.l + (rnd * 0.5 - 0.25) * strength));
  hsl.s = Math.max(0.1, Math.min(1.0, hsl.s + (rnd2 * 0.6 - 0.3) * strength));
  hsl.h = (hsl.h + (rnd * 0.16 - 0.08) * strength + 1.0) % 1.0;

  c.setHSL(hsl.h, hsl.s, hsl.l);
  return c;
};

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

export const getCellColor = (cell: Cell, mode: ViewMode, ctx: ColorContext): THREE.Color => {
  const { seaLevel, factionColors, cultureColors, religionColors, seasonalDelta } = ctx;
  const color = new THREE.Color();

  // D1: at a non-neutral season the render layer passes a per-cell temperature
  // excursion. Shown temperature = stored (annual-mean) + delta; the DISPLAYED
  // biome is re-derived from it for land cells only (never for water/lakes,
  // which are hydrology-derived and outside determineBiome's remit). cell.biome
  // and cell.temperature themselves are never mutated — civs/export stay canonical.
  //
  // VOLCANIC joins that exclusion list: it is TECTONIC, assigned from
  // `upliftAccum` in its own pass, and `determineBiome` cannot reproduce it from
  // (height, temp, moisture). Without this a volcanic field would silently
  // revert to whatever the climate says the moment the season slider left
  // neutral — and it would do so ONLY in the app, since `upliftAccum` does not
  // survive the worker boundary. Same trap as lakes, same fix.
  const seasonalTemp = cell.temperature + (seasonalDelta ?? 0);
  const displayBiome =
    seasonalDelta && cell.height >= seaLevel &&
    cell.biome !== BiomeType.LAKE && cell.biome !== BiomeType.SALT_LAKE &&
    cell.biome !== BiomeType.VOLCANIC
      ? determineBiome(cell.height, seasonalTemp, cell.moisture, seaLevel)
      : cell.biome;

  // D3: seasonal sea-ice — open-ocean cells whose seasonal temperature is below
  // seawater freezing render as ice, in the physical views only (satellite +
  // biome). Render overlay: cell.biome stays OCEAN (no civ/nav impact). Lakes
  // are excluded (they read as lakes, not sea-ice) and the data/elevation views
  // (height/temperature/…) are untouched so the datum is never obscured.
  if (
    (mode === 'satellite' || mode === 'biome') &&
    cell.height < seaLevel &&
    cell.biome !== BiomeType.LAKE && cell.biome !== BiomeType.SALT_LAKE &&
    seasonalTemp < SEAWATER_FREEZE_C
  ) {
    // Colder ice reads whiter; the thin edge keeps a pale blue cast.
    const iciness = Math.min(1, (SEAWATER_FREEZE_C - seasonalTemp) / 15);
    color.copy(new THREE.Color(0xbcd4e6)).lerp(new THREE.Color(0xffffff), iciness);
    return color;
  }

  switch (mode) {
    case 'satellite':
      if (cell.height < seaLevel) {
         const depth = cell.height / seaLevel;
         const deep = new THREE.Color(0x051e3e);
         const shallow = new THREE.Color(0x006994);
         color.copy(deep).lerp(shallow, Math.pow(depth, 2));
      } else if (cell.biome === BiomeType.LAKE || cell.biome === BiomeType.SALT_LAKE) {
         // Lakes are inland water — render as open water, not terrain.
         color.setHex(cell.biome === BiomeType.SALT_LAKE ? 0xd2e6df : 0x2f7fa6);
      } else {
         const t = (cell.height - seaLevel) / (1 - seaLevel);
         switch(displayBiome) {
             case BiomeType.ICE_CAP: color.setHex(0xffffff); break;
             case BiomeType.TUNDRA: color.setHex(0x78766a); break;
             case BiomeType.HOT_DESERT: color.setHex(0xdabba0); break;
             case BiomeType.COLD_DESERT: color.setHex(0x9e9587); break;
             case BiomeType.TROPICAL_RAINFOREST: color.setHex(0x052e16); break;
             case BiomeType.TEMPERATE_RAINFOREST: color.setHex(0x1a4221); break;
             case BiomeType.TEMPERATE_FOREST: color.setHex(0x285020); break; 
             case BiomeType.BOREAL_FOREST: color.setHex(0x193626); break;
             case BiomeType.TROPICAL_SAVANNA: color.setHex(0x6f7d46); break;
             case BiomeType.STEPPE: color.setHex(0x8a9263); break;
             case BiomeType.GRASSLAND: color.setHex(0x6a8f3c); break;
             case BiomeType.MEDITERRANEAN: color.setHex(0x6b7044); break;
             case BiomeType.BEACH: color.setHex(0xe8ddc5); break;
             case BiomeType.VOLCANIC: color.setHex(0x262626); break;
             default: color.setHex(0x335533);
         }
         let snowThreshold = 0.65;
         if (seasonalTemp > 20) snowThreshold = 0.85;
         if (t > 0.35 && displayBiome !== BiomeType.ICE_CAP) {
             const rockFactor = Math.min(1, (t - 0.35) * 4);
             color.lerp(new THREE.Color(0x524e49), rockFactor);
         }
         if (t > snowThreshold) {
             const snowFactor = Math.min(1, (t - snowThreshold) * 5);
             color.lerp(new THREE.Color(0xffffff), snowFactor);
         }
      }
      break;

    case 'height':
      if (cell.height < seaLevel) {
        const t = cell.height / seaLevel;
        color.setHSL(0.6, 0.7, 0.1 + t * 0.4); 
      } else {
        const t = (cell.height - seaLevel) / (1 - seaLevel);
        if (t < 0.2) {
           color.setHSL(0.25, 0.4, 0.3 + t * 0.5);
        } else if (t < 0.6) {
           color.setHSL(0.1, 0.2, 0.4 + (t - 0.2) * 0.5);
        } else {
           color.setHSL(0, 0, 0.5 + (t - 0.6) * 1.5);
        }
      }
      break;

    case 'height_bw':
      if (cell.height < seaLevel) {
         const t = cell.height / seaLevel;
         color.setScalar(t * 0.15); 
      } else {
         const t = (cell.height - seaLevel) / (1 - seaLevel);
         const val = 0.2 + Math.pow(t, 1.5) * 0.8;
         color.setScalar(val);
      }
      break;

    case 'temperature':
      const tNorm = Math.max(0, Math.min(1, (seasonalTemp - TEMP_MIN_C) / (TEMP_MAX_C - TEMP_MIN_C)));
      color.setHSL(0.65 - (tNorm * 0.65), 0.8, 0.5);
      break;

    case 'moisture':
      if (cell.height < seaLevel) {
         color.setHex(0x004488);
      } else {
         color.setHSL(0.6, cell.moisture, 0.9 - (cell.moisture * 0.5));
      }
      break;

    case 'plates':
      color.set(PLATE_COLORS[cell.plateId % PLATE_COLORS.length]);
      if (cell.height < seaLevel) color.multiplyScalar(0.7);
      break;
      
    case 'political':
      // Lakes read as water regardless of any territorial claim on them.
      if (cell.biome === BiomeType.LAKE || cell.biome === BiomeType.SALT_LAKE) {
        color.setHex(0x1a237e);
        color.multiplyScalar(0.7);
      } else if (cell.regionId !== undefined) {
        const baseColor = factionColors?.get(cell.regionId) ?? FACTION_COLORS[cell.regionId % FACTION_COLORS.length];
        if (cell.provinceId !== undefined) {
          color.copy(getProvinceVariant(baseColor, cell.provinceId));
        } else {
          color.set(baseColor);
        }
      } else if (cell.height < seaLevel) {
        color.setHex(0x1a237e);
        color.multiplyScalar(0.5 + cell.height * 0.5);
      } else {
        color.setHex(0x555555);
      }
      break;

    case 'province':
      // Distinct shade per province: same base hue family as its faction,
      // but with amplified variation so admin borders read at a glance
      if (cell.biome === BiomeType.LAKE || cell.biome === BiomeType.SALT_LAKE) {
        color.setHex(0x1a237e);
        color.multiplyScalar(0.7);
      } else if (cell.regionId !== undefined) {
        const baseColor = factionColors?.get(cell.regionId) ?? FACTION_COLORS[cell.regionId % FACTION_COLORS.length];
        if (cell.provinceId !== undefined) {
          color.copy(getProvinceVariant(baseColor, cell.provinceId, 1.8));
        } else {
          color.set(baseColor);
        }
      } else if (cell.height < seaLevel) {
        color.setHex(0x1a237e);
        color.multiplyScalar(0.5 + cell.height * 0.5);
      } else {
        color.setHex(0x555555);
      }
      break;

    case 'population': {
      if (cell.height < seaLevel) {
        color.setHex(0x0b1026);
        break;
      }
      const pop = cell.population || 0;
      if (pop <= 0) {
        color.setHex(0x263238);
        break;
      }
      // Log-scaled heat gradient: dim blue (~100) -> red -> bright yellow (100k+)
      const p = Math.min(1, Math.log10(1 + pop) / 5);
      color.setHSL(Math.max(0, 0.66 - p * 0.55), 0.85, 0.12 + p * 0.55);
      break;
    }

    case 'culture':
      // Water renders identically to political mode; unassigned land (should
      // only happen transiently, e.g. before recalculateCultures has run)
      // reads as dark grey rather than falling back to a faction color.
      if (cell.biome === BiomeType.LAKE || cell.biome === BiomeType.SALT_LAKE) {
        color.setHex(0x1a237e);
        color.multiplyScalar(0.7);
      } else if (cell.height < seaLevel) {
        color.setHex(0x1a237e);
        color.multiplyScalar(0.5 + cell.height * 0.5);
      } else if (cell.cultureId !== undefined) {
        const cultureColor = cultureColors?.get(cell.cultureId) ?? CULTURE_COLORS[cell.cultureId % CULTURE_COLORS.length];
        color.set(cultureColor);
      } else {
        color.setHex(0x333333);
      }
      break;

    case 'religion':
      // Same water/unassigned handling as 'culture' above — water reads as
      // political water, unassigned land (transient, pre-recalculation) is
      // dark grey rather than falling back to a faction or culture color.
      if (cell.biome === BiomeType.LAKE || cell.biome === BiomeType.SALT_LAKE) {
        color.setHex(0x1a237e);
        color.multiplyScalar(0.7);
      } else if (cell.height < seaLevel) {
        color.setHex(0x1a237e);
        color.multiplyScalar(0.5 + cell.height * 0.5);
      } else if (cell.religionId !== undefined) {
        const religionColor = religionColors?.get(cell.religionId) ?? RELIGION_COLORS[cell.religionId % RELIGION_COLORS.length];
        color.set(religionColor);
      } else {
        color.setHex(0x333333);
      }
      break;

    case 'biome':
    default:
      color.set(BIOME_COLORS[displayBiome] || '#ff00ff');
      break;
  }

  return color;
};

// ---------------------------------------------------------------------------
// F5: layer-aware legend. `legendFor(viewMode, world)` is the single source the
// legend panel reads, beside `getCellColor` so the two cannot drift. Ramp stops
// are SAMPLED from `getCellColor` (not a re-derived gradient), so a change to a
// fill formula moves the legend with it for free. Returns `null` for layers
// with no meaningful key (satellite, plates) — the panel is hidden there.
// ---------------------------------------------------------------------------

export interface LegendSwatch {
  label: string;
  color: string; // #rrggbb
}

export interface LegendRampStop {
  pos: number; // 0..1 along the bar
  color: string; // #rrggbb
}

export interface RampLegend {
  kind: 'ramp';
  title: string;
  stops: LegendRampStop[];
  minLabel: string;
  maxLabel: string;
  // Short-circuit cases the ramp itself never shows (ocean, uninhabited land) —
  // they cover large screen areas, so they must be labelled beside the bar.
  oob?: LegendSwatch[];
  // A marked position on the bar (e.g. the sea-level discontinuity in Elevation).
  tick?: { pos: number; label: string };
}

export interface SwatchesLegend {
  kind: 'swatches';
  title: string;
  items: LegendSwatch[];
}

export type Legend = RampLegend | SwatchesLegend;

const RAMP_STOPS = 16;
const toHex = (c: THREE.Color): string => `#${c.getHexString()}`;

// A synthetic land cell (height above the default sea level) used to probe the
// fill functions. Only the fields a given ramp sweeps are overridden.
const legendCell = (over: Partial<Cell>): Cell => ({
  id: 0,
  center: { x: 0, y: 0, z: 0 },
  vertices: [],
  neighbors: [],
  height: 0.7,
  plateId: 0,
  temperature: 15,
  moisture: 0.5,
  biome: BiomeType.GRASSLAND,
  ...over,
});

const sampleRamp = (
  mode: ViewMode,
  seaLevel: number,
  cellAt: (pos: number) => Cell,
): LegendRampStop[] =>
  Array.from({ length: RAMP_STOPS + 1 }, (_, i) => {
    const pos = i / RAMP_STOPS;
    return { pos, color: toHex(getCellColor(cellAt(pos), mode, { seaLevel })) };
  });

export const legendFor = (mode: ViewMode, world: WorldData): Legend | null => {
  const seaLevel = world.params.seaLevel;
  const maxElevationM = world.params.maxElevationM;
  const ctx: ColorContext = { seaLevel };
  // A representative ocean cell — well below sea level so every water short-circuit fires.
  const oceanColor = () => toHex(getCellColor(legendCell({ height: seaLevel * 0.4 }), mode, ctx));

  switch (mode) {
    // No categorical key and no single ramp — hide the panel entirely.
    case 'satellite':
    case 'plates':
      return null;

    case 'biome': {
      const items: LegendSwatch[] = (Object.keys(BIOME_COLORS) as BiomeType[]).map((b) => ({
        label: b,
        color: BIOME_COLORS[b],
      }));
      // Sea-ice paints open ocean white in biome/satellite at cold seasons, yet
      // it is a render overlay (cell.biome stays OCEAN) with no BIOME_COLORS
      // entry — so it needs its own swatch to name what's on screen.
      items.push({ label: 'Sea Ice', color: '#ffffff' });
      return { kind: 'swatches', title: 'Biomes', items };
    }

    case 'political':
    case 'province': {
      // Province fills are derived shade-variants of each faction's base colour,
      // so a faction swatch is the honest key for both layers (Matt's call, S46).
      const items: LegendSwatch[] = (world.civData?.factions ?? []).map((f) => ({
        label: f.name,
        color: f.color,
      }));
      items.push({ label: 'Unclaimed', color: '#555555' });
      return { kind: 'swatches', title: mode === 'province' ? 'Provinces' : 'Factions', items };
    }

    case 'culture':
      return {
        kind: 'swatches',
        title: 'Cultures',
        items: (world.cultures ?? []).map((c) => ({ label: c.name, color: c.color })),
      };

    case 'religion':
      return {
        kind: 'swatches',
        title: 'Religions',
        items: (world.religions ?? []).map((r) => ({ label: r.name, color: r.color })),
      };

    case 'temperature': {
      const stops = sampleRamp(mode, seaLevel, (pos) =>
        legendCell({ temperature: TEMP_MIN_C + pos * (TEMP_MAX_C - TEMP_MIN_C) }),
      );
      return {
        kind: 'ramp',
        title: 'Temperature',
        stops,
        minLabel: `${TEMP_MIN_C}°C`,
        maxLabel: `${TEMP_MAX_C}°C`,
      };
    }

    case 'moisture': {
      const stops = sampleRamp(mode, seaLevel, (pos) => legendCell({ moisture: pos }));
      return {
        kind: 'ramp',
        title: 'Rainfall',
        stops,
        minLabel: 'Dry',
        maxLabel: 'Wet',
        oob: [{ label: 'Ocean', color: oceanColor() }],
      };
    }

    case 'population': {
      // The fill maps population by log10(1 + pop) / 5, saturating at ~100k. We
      // invert that so each sampled population lands at an even screen position.
      // Floor at 1 so the first stop isn't pop=0 — that hits the uninhabited
      // short-circuit and would open the bar on the grey already shown as OOB.
      const stops = sampleRamp(mode, seaLevel, (pos) =>
        legendCell({ population: Math.max(1, Math.pow(10, pos * 5) - 1) }),
      );
      return {
        kind: 'ramp',
        title: 'Population',
        stops,
        minLabel: '0',
        maxLabel: '100k+',
        oob: [
          { label: 'Ocean', color: oceanColor() },
          {
            label: 'Uninhabited',
            color: toHex(getCellColor(legendCell({ population: 0 }), mode, ctx)),
          },
        ],
      };
    }

    case 'height':
    case 'height_bw': {
      // Full 0->1 field, since the fill is discontinuous at sea level; the tick
      // marks that break. Elevation labels go through datum.ts (frac^2 curve),
      // matching the Inspector rather than a linear map.
      const stops = sampleRamp(mode, seaLevel, (pos) => legendCell({ height: pos }));
      return {
        kind: 'ramp',
        title: mode === 'height_bw' ? 'Elevation (BW)' : 'Elevation',
        stops,
        minLabel: formatElevation(0, seaLevel, maxElevationM),
        maxLabel: formatElevation(1, seaLevel, maxElevationM),
        tick: { pos: seaLevel, label: 'Sea level' },
      };
    }
  }

  return null;
};
