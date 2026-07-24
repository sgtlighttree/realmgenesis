import { BiomeType, type Cell } from '../types';

// Binary min-heap. Moved verbatim out of worldGen.ts so both the civ Dijkstra
// and the route A* share one implementation (no circular import: this module
// depends only on types.ts).
export class MinHeap<T> {
    private heap: T[];
    private scoreFunction: (t: T) => number;

    constructor(scoreFunction: (t: T) => number) {
        this.heap = [];
        this.scoreFunction = scoreFunction;
    }

    push(node: T) {
        this.heap.push(node);
        this.bubbleUp(this.heap.length - 1);
    }

    pop(): T | undefined {
        if (this.heap.length === 0) return undefined;
        const top = this.heap[0];
        const bottom = this.heap.pop();
        if (this.heap.length > 0 && bottom !== undefined) {
            this.heap[0] = bottom;
            this.sinkDown(0);
        }
        return top;
    }

    size(): number { return this.heap.length; }

    private bubbleUp(index: number) {
        while (index > 0) {
            const parentIndex = Math.floor((index - 1) / 2);
            if (this.scoreFunction(this.heap[index]) >= this.scoreFunction(this.heap[parentIndex])) break;
            [this.heap[index], this.heap[parentIndex]] = [this.heap[parentIndex], this.heap[index]];
            index = parentIndex;
        }
    }

    private sinkDown(index: number) {
        const length = this.heap.length;
        const element = this.heap[index];
        const elemScore = this.scoreFunction(element);

        while (true) {
            let leftChildIdx = 2 * index + 1;
            let rightChildIdx = 2 * index + 2;
            let leftScore, rightScore;
            let swap = null;

            if (leftChildIdx < length) {
                leftScore = this.scoreFunction(this.heap[leftChildIdx]);
                if (leftScore < elemScore) swap = leftChildIdx;
            }
            if (rightChildIdx < length) {
                rightScore = this.scoreFunction(this.heap[rightChildIdx]);
                if (swap === null) {
                    if (rightScore < elemScore) swap = rightChildIdx;
                } else {
                    if (rightScore < leftScore!) swap = rightChildIdx;
                }
            }

            if (swap === null) break;
            [this.heap[index], this.heap[swap]] = [this.heap[swap], this.heap[index]];
            index = swap;
        }
    }
}

export function isWaterCell(cell: Cell, seaLevel: number): boolean {
    return cell.height < seaLevel || cell.biome === BiomeType.LAKE || cell.biome === BiomeType.SALT_LAKE;
}

// Terrain-only step cost, extracted from recalculateCivs. Base land cost 1,
// biome multipliers on the destination, slope penalty. NO water override and
// NO civ-specific multipliers (borderRoughness / costMult) — callers add those.
export function landTerrainStepCost(from: Cell, to: Cell): number {
    let c = 1;
    if (to.biome === BiomeType.ICE_CAP) c *= 4;
    if (to.biome === BiomeType.HOT_DESERT) c *= 2;
    if (to.biome === BiomeType.VOLCANIC) c *= 5;
    c += Math.abs(to.height - from.height) * 20;
    return c;
}

// Directed A* search over the cell graph. stepCost returns Infinity for
// impassable edges. heuristic is a per-cell lower-ish bound on remaining cost
// (accelerator; routes need only be good and deterministic, not provably
// shortest). Returns the cell-id path or null if unreachable / cap exceeded.
export function aStar(
    cells: Cell[],
    startId: number,
    goalId: number,
    stepCost: (fromId: number, toId: number) => number,
    heuristic: (id: number) => number,
    maxExpansions: number,
): number[] | null {
    const open = new MinHeap<{ id: number; f: number }>(x => x.f);
    const g = new Map<number, number>();
    const cameFrom = new Map<number, number>();
    const closed = new Set<number>();
    g.set(startId, 0);
    open.push({ id: startId, f: heuristic(startId) });
    let expansions = 0;

    while (open.size() > 0) {
        const { id } = open.pop()!;
        if (id === goalId) {
            const path: number[] = [id];
            let cur = id;
            while (cameFrom.has(cur)) { cur = cameFrom.get(cur)!; path.push(cur); }
            return path.reverse();
        }
        if (closed.has(id)) continue;
        closed.add(id);
        if (++expansions > maxExpansions) return null;
        // Deterministic neighbor order (neighbors arrays are already stable).
        for (const nId of cells[id].neighbors) {
            const step = stepCost(id, nId);
            if (!isFinite(step)) continue;
            const tentative = g.get(id)! + step;
            if (!g.has(nId) || tentative < g.get(nId)!) {
                g.set(nId, tentative);
                cameFrom.set(nId, id);
                open.push({ id: nId, f: tentative + heuristic(nId) });
            }
        }
    }
    return null;
}
