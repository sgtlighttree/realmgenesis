import { Cell, BiomeType, UndoSnapshot } from '../types';
import { determineBiome } from './worldGen';

export function getCellsInRadius(
    center: Cell,
    rings: number,
    cells: Cell[]
): { cell: Cell; ring: number }[] {
    // Cap by geographic chord distance to prevent polar/antimeridian topology artefacts
    const estimatedRadius = Math.sqrt(4 * Math.PI / cells.length) * 0.9;
    const maxChord = (rings + 1.5) * estimatedRadius;
    const maxChordSq = maxChord * maxChord;
    const cx = center.center.x, cy = center.center.y, cz = center.center.z;

    const visited = new Set<number>([center.id]);
    const result: { cell: Cell; ring: number }[] = [{ cell: center, ring: 0 }];
    let frontier = [center];
    for (let r = 1; r <= rings; r++) {
        const next: Cell[] = [];
        for (const c of frontier) {
            for (const nId of c.neighbors) {
                if (visited.has(nId)) continue;
                const nc = cells[nId];
                const dx = nc.center.x - cx, dy = nc.center.y - cy, dz = nc.center.z - cz;
                if (dx * dx + dy * dy + dz * dz > maxChordSq) continue;
                visited.add(nId);
                next.push(nc);
                result.push({ cell: nc, ring: r });
            }
        }
        frontier = next;
    }
    return result;
}

export function snapshotCells(brush: { cell: Cell }[]): UndoSnapshot {
    return {
        cells: new Map(brush.map(({ cell: c }) => [
            c.id,
            { height: c.height, biome: c.biome, regionId: c.regionId, provinceId: c.provinceId }
        ]))
    };
}

export function applyTerrainStroke(
    brush: { cell: Cell; ring: number }[],
    brushSize: number,
    direction: 'raise' | 'lower',
    style: 'adaptive' | 'freeform',
    cells: Cell[],
    strength: number = 0.5
): void {
    const BASE_DELTA = style === 'freeform' ? 0.06 : 0.025;
    const sign = direction === 'raise' ? 1 : -1;

    brush.forEach(({ cell, ring }) => {
        const weight = 1 - ring / (brushSize + 1);
        cell.height = Math.max(0, Math.min(1, cell.height + sign * BASE_DELTA * strength * weight));
    });

    if (style === 'adaptive') {
        brush.forEach(({ cell, ring }) => {
            if (ring === 0) return;
            const blendWeight = (ring / (brushSize + 1)) * 0.4;
            const avg = cell.neighbors.reduce((s, nId) => s + cells[nId].height, 0) / cell.neighbors.length;
            cell.height = Math.max(0, Math.min(1, cell.height * (1 - blendWeight) + avg * blendWeight));
        });
    }
}

export function applyFlattenStroke(
    brush: { cell: Cell; ring: number }[],
    brushSize: number,
    targetHeight: number,
    strength: number
): void {
    brush.forEach(({ cell, ring }) => {
        const weight = 1 - ring / (brushSize + 1);
        const pull = (targetHeight - cell.height) * strength * weight * 0.4;
        cell.height = Math.max(0, Math.min(1, cell.height + pull));
    });
}

export function applySmoothStroke(
    brush: { cell: Cell; ring: number }[],
    brushSize: number,
    strength: number,
    cells: Cell[]
): void {
    brush.forEach(({ cell, ring }) => {
        const weight = (1 - ring / (brushSize + 1)) * strength * 0.6;
        const avg = cell.neighbors.reduce((s, nId) => s + cells[nId].height, 0) / cell.neighbors.length;
        cell.height = Math.max(0, Math.min(1, cell.height * (1 - weight) + avg * weight));
    });
}

export function applyPoliticalStroke(
    brush: { cell: Cell; ring: number }[],
    factionId: number | undefined,
    provinceId?: number
): void {
    brush.forEach(({ cell }) => {
        cell.regionId = factionId;
        cell.provinceId = provinceId;
    });
}

export function applyBiomeStroke(
    brush: { cell: Cell; ring: number }[],
    biome: BiomeType
): void {
    brush.forEach(({ cell }) => { cell.biome = biome; });
}

export function refreshBiomes(brushCells: Cell[], seaLevel: number): void {
    brushCells.forEach(c => {
        // Lakes are hydrology-derived and determineBiome can never reproduce
        // them; leave them intact so a terrain brush doesn't turn a lake into a
        // stray land biome sitting in a depression. (Lake re-solve is B2.)
        if (c.biome === BiomeType.LAKE || c.biome === BiomeType.SALT_LAKE) return;
        c.biome = determineBiome(c.height, c.temperature, c.moisture, seaLevel);
    });
}
