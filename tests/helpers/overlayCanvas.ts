import { Cell, WorldData } from '../../types';

// Test doubles for the F2 ScreenOverlay tenant seam. Tenants are pure draw
// functions over a CanvasRenderingContext2D and a LocalProjector, so both can
// be recorded without jsdom canvas support.

export interface CanvasOp {
  op: 'beginPath' | 'moveTo' | 'lineTo' | 'stroke' | 'fill' | 'arc' | 'setLineDash';
  args: number[];
}

export interface FakeCtx {
  ops: CanvasOp[];
  dashHistory: number[][]; // every setLineDash argument, in call order
  strokeStyle: string;
  fillStyle: string;
  lineWidth: number;
  lineDash: number[];
}

export function makeFakeCtx(): FakeCtx & CanvasRenderingContext2D {
  const ops: CanvasOp[] = [];
  const dashHistory: number[][] = [];
  const ctx = {
    ops,
    dashHistory,
    strokeStyle: '',
    fillStyle: '',
    lineWidth: 1,
    lineDash: [] as number[],
    beginPath() { ops.push({ op: 'beginPath', args: [] }); },
    moveTo(x: number, y: number) { ops.push({ op: 'moveTo', args: [x, y] }); },
    lineTo(x: number, y: number) { ops.push({ op: 'lineTo', args: [x, y] }); },
    stroke() { ops.push({ op: 'stroke', args: [] }); },
    fill() { ops.push({ op: 'fill', args: [] }); },
    arc(x: number, y: number, r: number) { ops.push({ op: 'arc', args: [x, y, r] }); },
    setLineDash(d: number[]) {
      ctx.lineDash = d;
      dashHistory.push([...d]);
      ops.push({ op: 'setLineDash', args: [...d] });
    },
  };
  return ctx as unknown as FakeCtx & CanvasRenderingContext2D;
}

export interface RecordedProjector {
  project: (x: number, y: number, z: number, out: [number, number]) => boolean;
  radii: number[]; // |P| of every point handed to the projector
}

// `visible` decides culling per point; default = everything visible. Screen
// coords are a trivial stand-in (x, y scaled) — tenants only chain them.
export function makeProjector(
  visible: (x: number, y: number, z: number) => boolean = () => true,
): RecordedProjector {
  const radii: number[] = [];
  return {
    radii,
    project(x, y, z, out) {
      radii.push(Math.hypot(x, y, z));
      if (!visible(x, y, z)) return false;
      out[0] = x * 100;
      out[1] = y * 100;
      return true;
    },
  };
}

// Minimal cell ring on the unit sphere in the XZ plane, each cell linked to its
// two ring neighbors — enough for nearestCellWalk to hill-climb.
export function makeRingWorld(n: number, heights: number[], seaLevel = 0.5): WorldData {
  const cells: Cell[] = [];
  for (let i = 0; i < n; i++) {
    const a = (i / n) * Math.PI * 2;
    cells.push({
      id: i,
      center: { x: Math.cos(a), y: 0, z: Math.sin(a) },
      height: heights[i % heights.length],
      neighbors: [(i + n - 1) % n, (i + 1) % n],
    } as unknown as Cell);
  }
  return { cells, params: { seaLevel } } as unknown as WorldData;
}
