import { WorldData, Cell, LabelVisibility } from '../types';
import { normalizeVec, Point3 } from './geo';

export type LabelKind =
  | 'faction' | 'capital' | 'town' | 'province'
  | 'range' | 'desert' | 'forest' | 'sea' | 'ocean' | 'island' | 'lake'
  | 'marker';

export interface MapLabel {
  kind: LabelKind;
  name: string;
  position: { x: number; y: number; z: number };
  priority: number;
  factionId?: number;
}

export interface LabelStyleConfig {
  fontWeight: number;
  baseFontSize: number;
  alpha: number;
  uppercase: boolean;
  visibilityKey: keyof LabelVisibility;
  italic?: boolean; // cartographic feel for geographic (esp. water) labels
  fill?: string; // override the default light fill (water kinds run slightly blued)
}

// Exported so vector exporters (SVG/GeoJSON) can mirror the exact per-kind
// styling used by the canvas label pass instead of duplicating the table.
export const LABEL_CONFIG: Record<LabelKind, LabelStyleConfig> = {
  faction: { fontWeight: 700, baseFontSize: 14, alpha: 1.0, uppercase: true, visibilityKey: 'factions' },
  capital: { fontWeight: 700, baseFontSize: 11, alpha: 0.95, uppercase: false, visibilityKey: 'capitals' },
  province: { fontWeight: 400, baseFontSize: 10, alpha: 0.85, uppercase: false, visibilityKey: 'provinces' },
  town: { fontWeight: 400, baseFontSize: 9, alpha: 0.75, uppercase: false, visibilityKey: 'towns' },
  // Geographic features (B3) — one shared visibility toggle. Water kinds italic + blued.
  ocean: { fontWeight: 400, baseFontSize: 13, alpha: 0.9, uppercase: false, visibilityKey: 'geography', italic: true, fill: '#dbeafe' },
  sea: { fontWeight: 400, baseFontSize: 10, alpha: 0.85, uppercase: false, visibilityKey: 'geography', italic: true, fill: '#dbeafe' },
  lake: { fontWeight: 400, baseFontSize: 9, alpha: 0.8, uppercase: false, visibilityKey: 'geography', italic: true, fill: '#dbeafe' },
  range: { fontWeight: 400, baseFontSize: 10, alpha: 0.85, uppercase: false, visibilityKey: 'geography' },
  desert: { fontWeight: 400, baseFontSize: 10, alpha: 0.85, uppercase: false, visibilityKey: 'geography' },
  forest: { fontWeight: 400, baseFontSize: 10, alpha: 0.85, uppercase: false, visibilityKey: 'geography' },
  island: { fontWeight: 400, baseFontSize: 9, alpha: 0.8, uppercase: false, visibilityKey: 'geography' },
  // Warm amber, distinct from both civ (blue/white) and geo (light/blued) labels.
  marker: { fontWeight: 600, baseFontSize: 10, alpha: 0.95, uppercase: false, visibilityKey: 'markers', fill: '#fde68a' },
};

const ZOOM_THRESHOLDS: Record<LabelKind, number> = {
  faction: 0,
  capital: 0.7,
  province: 1.2,
  town: 2.0,
  // Oceans always; mid-scale geography from ~1.0; small features only when zoomed in.
  ocean: 0,
  sea: 1.0,
  range: 1.0,
  desert: 1.0,
  forest: 1.0,
  island: 1.8,
  lake: 1.8,
  marker: 0.7,
};

interface LabelRect {
  x: number;
  y: number;
  w: number;
  h: number;
}

const rectsOverlap = (a: LabelRect, b: LabelRect): boolean =>
  a.x < b.x + b.w && a.x + a.w > b.x && a.y < b.y + b.h && a.y + a.h > b.y;

const computeCentroid = (cells: Cell[]): Point3 => {
  let sumX = 0, sumY = 0, sumZ = 0;
  for (const cell of cells) {
    sumX += cell.center.x;
    sumY += cell.center.y;
    sumZ += cell.center.z;
  }
  return normalizeVec([sumX, sumY, sumZ]);
};

// Geographic features (B3) carry their own anchor; priorities interleave with
// the civ labels so the shared declutter pass ranks them sensibly. Oceans sit
// just below faction names; mid kinds around provinces; islands/lakes lowest.
const GEO_PRIORITY: Record<string, number> = {
  ocean: 0.5,
  sea: 2.5,
  range: 2.5,
  desert: 2.5,
  forest: 2.5,
  island: 3.5,
  lake: 3.5,
};

export const collectLabels = (world: WorldData): MapLabel[] => {
  const labels: MapLabel[] = [];
  const seaLevel = world.params.seaLevel;

  // Geographic labels are terrain-derived — emitted even when civData is absent.
  for (const feature of world.features ?? []) {
    labels.push({
      kind: feature.kind,
      name: feature.name,
      position: { x: feature.anchor.x, y: feature.anchor.y, z: feature.anchor.z },
      priority: GEO_PRIORITY[feature.kind] ?? 3,
    });
  }

  // User-placed markers are sphere-anchored, not civ/terrain-derived — emitted
  // independent of civData, same as geo features above.
  for (const marker of world.markers ?? []) {
    labels.push({
      kind: 'marker',
      name: marker.name,
      position: { x: marker.position.x, y: marker.position.y, z: marker.position.z },
      priority: 1.5,
    });
  }

  if (!world.civData) {
    labels.sort((a, b) => a.priority - b.priority);
    return labels;
  }

  const cellsByFaction = new Map<number, Cell[]>();
  const cellsByProvince = new Map<string, Cell[]>();

  for (const cell of world.cells) {
    if (cell.regionId !== undefined) {
      let bucket = cellsByFaction.get(cell.regionId);
      if (!bucket) { bucket = []; cellsByFaction.set(cell.regionId, bucket); }
      bucket.push(cell);

      if (cell.provinceId !== undefined) {
        const key = `${cell.regionId}-${cell.provinceId}`;
        let provBucket = cellsByProvince.get(key);
        if (!provBucket) { provBucket = []; cellsByProvince.set(key, provBucket); }
        provBucket.push(cell);
      }
    }
  }

  for (const faction of world.civData.factions) {
    const factionCells = cellsByFaction.get(faction.id) || [];
    const landCells = factionCells.filter(c => c.height >= seaLevel);
    const labelCells = landCells.length > 0 ? landCells : factionCells;

    if (labelCells.length > 0) {
      const [nx, ny, nz] = computeCentroid(labelCells);
      labels.push({
        kind: 'faction',
        name: faction.name,
        position: { x: nx, y: ny, z: nz },
        priority: 0,
        factionId: faction.id,
      });
    }

    for (const province of faction.provinces) {
      const key = `${faction.id}-${province.id}`;
      const provinceCells = cellsByProvince.get(key) || [];
      const provLandCells = provinceCells.filter(c => c.height >= seaLevel);
      const provLabelCells = provLandCells.length > 0 ? provLandCells : provinceCells;

      if (provLabelCells.length > 0) {
        const [nx, ny, nz] = computeCentroid(provLabelCells);
        labels.push({
          kind: 'province',
          name: province.name,
          position: { x: nx, y: ny, z: nz },
          priority: 2,
          factionId: faction.id,
        });
      }

      for (const town of province.towns) {
        const cell = world.cells[town.cellId];
        if (!cell) continue;
        const [nx, ny, nz] = normalizeVec([cell.center.x, cell.center.y, cell.center.z]);
        labels.push({
          kind: town.isCapital ? 'capital' : 'town',
          name: town.name,
          position: { x: nx, y: ny, z: nz },
          priority: town.isCapital ? 1 : 3,
          factionId: faction.id,
        });
      }
    }
  }

  labels.sort((a, b) => a.priority - b.priority);
  return labels;
};

export const drawMapLabels = (
  ctx: CanvasRenderingContext2D,
  labels: MapLabel[],
  project: (position: { x: number; y: number; z: number }) => [number, number] | null,
  scale: number,
  visibility: LabelVisibility,
  fontScale?: number, // override for high-resolution export canvases
): void => {
  if (labels.length === 0) return;

  const effectiveFontScale = fontScale ?? Math.min(scale, 2);
  const placed: LabelRect[] = [];

  ctx.save();

  for (const label of labels) {
    const config = LABEL_CONFIG[label.kind];
    if (!visibility[config.visibilityKey]) continue;
    if (scale < ZOOM_THRESHOLDS[label.kind]) continue;

    const projected = project(label.position);
    if (!projected) continue;

    const fontSize = Math.max(8, config.baseFontSize * effectiveFontScale);
    const text = config.uppercase ? label.name.toUpperCase() : label.name;
    const padding = Math.max(4, fontSize * 0.3);

    const style = config.italic ? 'italic ' : '';
    ctx.font = `${style}${config.fontWeight} ${fontSize}px Inter, ui-sans-serif, system-ui, sans-serif`;
    const textWidth = ctx.measureText(text).width;
    const textHeight = fontSize * 1.2;

    const rect: LabelRect = {
      x: projected[0] - textWidth / 2 - padding,
      y: projected[1] - textHeight / 2 - padding,
      w: textWidth + padding * 2,
      h: textHeight + padding * 2,
    };

    if (placed.some(r => rectsOverlap(rect, r))) continue;
    placed.push(rect);

    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    ctx.lineJoin = 'round';
    ctx.lineWidth = Math.max(3, fontSize * 0.28);
    ctx.strokeStyle = 'rgba(2, 6, 23, 0.95)';
    ctx.fillStyle = config.fill ?? '#f8fafc';
    ctx.globalAlpha = config.alpha;
    ctx.strokeText(text, projected[0], projected[1]);
    ctx.fillText(text, projected[0], projected[1]);
  }

  ctx.restore();
};
