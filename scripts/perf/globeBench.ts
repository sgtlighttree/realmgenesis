// F4 globe-render instrument. Salvaged from the session-28 agent run.
//
// Deliberately NOT a `*.test.ts`: it generates a 30k world and measures for up
// to 15 minutes, which has no business in `npm test`. It lives here and is run
// on purpose, through `scripts/renderGlobeBench.mjs`.
//
// What it measures, and why the shape matters: each overlay tenant is drawn
// against a NULL 2D context, so the number is the JS cost of the tenant — the
// projection, the culling, the polyline building — with rasterisation removed.
// That is the part an optimisation can actually move.
//
// It also pairs every candidate optimisation with a CORRECTNESS number rather
// than a visual check: `max screen error (px)` and `visibility flips` against
// the baseline maths. A faster projection that moves a point half a pixel or
// flips one cell across the horizon is not a win, and neither shows up in a
// frame-time graph.
import { appendFileSync, writeFileSync } from 'node:fs';

const OUT = process.env.F4_OUT ?? '/tmp/f4bench.txt';
const say = (m: string) => { appendFileSync(OUT, m + '\n'); };
import * as THREE from 'three';

import { generateWorld } from '../../utils/worldGen';
import { makeParams } from '../../tests/helpers';
import { WorldData } from '../../types';
import { isVisible, ProjectedCells } from '../../utils/screenProject';
import { displayRadius, CELL_OVERHANG } from '../../utils/displayRadius';
import { getCellColor } from '../../utils/colors';
import { seasonalTemperatureDelta } from '../../utils/seasons';
import { computeShadeMap, computeContourSegments, contourInterval } from '../../utils/shading';
import { computeBorderSegments } from '../../utils/borders';
import { computeRiverPolylines } from '../../utils/riverPaths';
import { collectLabels } from '../../utils/labels';
import { computeLabelAnchors } from '../../utils/labelAnchors';
import { buildFactionColorMap } from '../../utils/colors';
import {
  drawCurrentsTenant, drawGraticuleTenant, drawRoutesTenant, drawContoursTenant,
  drawBordersTenant, drawRiversTenant, drawLabelsTenant, drawDymaxionTenant,
} from '../../components/overlays/tenants';
import { cageEdges } from '../../utils/dymaxionCage';
import { DEFAULT_LABEL_THEME } from '../../utils/mapStyle/labelTheme';
import { DEFAULT_LABEL_VISIBILITY } from '../../types';

const POINTS = Number(process.env.F4_POINTS ?? 30000);
const FRAMES = Number(process.env.F4_FRAMES ?? 120);

// A do-nothing 2D context: measures the JS cost of a tenant draw, not raster.
function nullCtx(): CanvasRenderingContext2D {
  const c = {
    strokeStyle: '', fillStyle: '', lineWidth: 1, lineCap: 'butt', lineJoin: 'miter',
    font: '', textAlign: 'left', textBaseline: 'alphabetic', letterSpacing: '0px',
    globalAlpha: 1, miterLimit: 10,
    beginPath() {}, moveTo() {}, lineTo() {}, stroke() {}, fill() {}, arc() {},
    setLineDash() {}, save() {}, restore() {}, strokeText() {}, fillText() {},
    clearRect() {}, setTransform() {}, rect() {}, closePath() {}, translate() {},
    rotate() {}, scale() {}, quadraticCurveTo() {}, bezierCurveTo() {},
    measureText(t: string) { return { width: t.length * 6 } as TextMetrics; },
  };
  return c as unknown as CanvasRenderingContext2D;
}

function bench(label: string, frames: number, fn: (i: number) => void): number {
  for (let i = 0; i < 10; i++) fn(i); // warm the JIT
  const t0 = performance.now();
  for (let i = 0; i < frames; i++) fn(i);
  const ms = (performance.now() - t0) / frames;
  say(`  ${label.padEnd(46)} ${ms.toFixed(3)} ms/frame`);
  return ms;
}

// ---------------------------------------------------------------------------
// BASELINE transcription of ScreenOverlay's useFrame body (commit fce8df4).
// ---------------------------------------------------------------------------
const IDENT = new THREE.Matrix4().elements;
function matrixKeyBaseline(a: THREE.Matrix4, b: Float32Array | number[]): string {
  const ea = a.elements;
  let s = '';
  for (let i = 0; i < 16; i++) s += ((ea[i] * 1000) | 0) + ',' + ((b[i] * 1000) | 0) + ';';
  return s;
}

export async function run(): Promise<void> {
  {
    writeFileSync(OUT, '');
    const t0 = performance.now();
    const world: WorldData = await generateWorld(
      makeParams({ points: POINTS, plates: 12, numFactions: 8 }),
    );
    say(`world: ${world.cells.length} cells in ${(performance.now() - t0) / 1000 | 0}s`);

    const cssW = 1440, cssH = 900;
    const camera = new THREE.PerspectiveCamera(45, cssW / cssH, 0.1, 1000);
    camera.position.set(0, 0, 2.5);
    camera.updateMatrixWorld();
    camera.updateProjectionMatrix();

    const globe = new THREE.Object3D();
    globe.name = 'globe-mesh';
    const proj: ProjectedCells = {
      x: new Float32Array(world.cells.length),
      y: new Float32Array(world.cells.length),
      visible: new Uint8Array(world.cells.length),
      n: world.cells.length,
    };

    const scratch = new THREE.Vector3();
    const projScratch = new THREE.Vector3();
    const camPos = new THREE.Vector3();
    const smoothGlobe = false;
    const nCells = world.cells.length;

    const advance = (i: number) => {
      globe.rotation.y = i * 0.001;
      globe.updateWorldMatrix(true, false);
    };

    // ---- 1. per-cell projection loop (baseline) ----
    const projectAllBaseline = () => {
      const gm = globe.matrixWorld;
      camera.getWorldPosition(camPos);
      for (let i = 0; i < nCells; i++) {
        const cell = world.cells[i];
        const c = cell.center;
        const r = displayRadius(cell.height, smoothGlobe);
        scratch.set(c.x * r, c.y * r, c.z * r);
        scratch.applyMatrix4(gm);
        if (!isVisible(scratch.x, scratch.y, scratch.z, camPos.x, camPos.y, camPos.z)) {
          proj.visible[i] = 0;
          continue;
        }
        scratch.project(camera);
        proj.x[i] = (scratch.x + 1) / 2 * cssW;
        proj.y[i] = (1 - (scratch.y + 1) / 2) * cssH;
        proj.visible[i] = 1;
      }
    };

    const projectorBaseline = (x: number, y: number, z: number, out: [number, number]) => {
      const gm = globe.matrixWorld;
      projScratch.set(x, y, z);
      projScratch.applyMatrix4(gm);
      if (!isVisible(projScratch.x, projScratch.y, projScratch.z, camPos.x, camPos.y, camPos.z)) return false;
      projScratch.project(camera);
      out[0] = (projScratch.x + 1) / 2 * cssW;
      out[1] = (1 - (projScratch.y + 1) / 2) * cssH;
      return true;
    };

    say('\n--- ScreenOverlay per-frame CPU ---');
    bench('gate: matrixKey string build', FRAMES, (i) => {
      advance(i);
      const active = [{ id: 'rivers' }, { id: 'contours' }, { id: 'graticule' }, { id: 'labels:11111' }];
      const k = (smoothGlobe ? 's|' : 'r|') + active.map((t) => t.id).join(',') + '|'
        + matrixKeyBaseline(camera.matrixWorld, globe.matrixWorld.elements ?? IDENT);
      if (k.length === 0) throw new Error('x');
    });

    bench(`project all ${nCells} cells (baseline)`, FRAMES, (i) => {
      advance(i);
      projectAllBaseline();
    });

    // ---- 2. tenant draws ----
    const ctx = nullCtx();
    const seaLevel = world.params.seaLevel;
    const segs = computeContourSegments(world.cells, seaLevel, contourInterval(world.cells, seaLevel));
    const borders = computeBorderSegments(world.cells);
    const rivers = computeRiverPolylines(world.rivers ?? []);
    const labels = collectLabels(world);
    const anchors = computeLabelAnchors(labels, world.cells);
    const edges = cageEdges({ lon: 0, lat: 0, roll: 0 } as never);
    say(`  [inputs] contours=${segs.length} borders=${borders.length} rivers=${rivers.length} labels=${anchors.length} routes=${world.routes?.length ?? 0}`);

    say('\n--- tenant draw cost (baseline projector, null canvas) ---');
    bench('graticule', FRAMES, (i) => { advance(i); camera.getWorldPosition(camPos); drawGraticuleTenant(ctx, proj, world, projectorBaseline, smoothGlobe); });
    bench('contours', FRAMES, (i) => { advance(i); camera.getWorldPosition(camPos); drawContoursTenant(ctx, segs, projectorBaseline, smoothGlobe); });
    bench('borders', FRAMES, (i) => { advance(i); camera.getWorldPosition(camPos); drawBordersTenant(ctx, borders, projectorBaseline, smoothGlobe); });
    bench('rivers', FRAMES, (i) => { advance(i); camera.getWorldPosition(camPos); drawRiversTenant(ctx, rivers, projectorBaseline, smoothGlobe); });
    bench('routes', FRAMES, (i) => { advance(i); camera.getWorldPosition(camPos); drawRoutesTenant(ctx, proj, world, projectorBaseline, smoothGlobe); });
    bench('dymaxion', FRAMES, (i) => { advance(i); camera.getWorldPosition(camPos); drawDymaxionTenant(ctx, edges, projectorBaseline, smoothGlobe); });
    bench('labels', FRAMES, (i) => { advance(i); camera.getWorldPosition(camPos); drawLabelsTenant(ctx, anchors, projectorBaseline, smoothGlobe, 2.5, DEFAULT_LABEL_VISIBILITY, DEFAULT_LABEL_THEME); });
    if (world.currents) {
      advance(0); camera.getWorldPosition(camPos); projectAllBaseline();
      bench('currents (needs proj)', FRAMES, () => { drawCurrentsTenant(ctx, proj, world, projectorBaseline, smoothGlobe); });
    }

    // ---- 3. WorldViewer geometry refill ----
    let triCount = 0;
    for (const cell of world.cells) triCount += cell.vertices.length;
    const pos = new Float32Array(triCount * 9);
    const col = new Float32Array(triCount * 9);
    const shadeMap = computeShadeMap(world.cells, seaLevel);
    const factionColors = world.civData ? buildFactionColorMap(world.civData) : undefined;

    const refillBaseline = (viewMode: 'biome' | 'political', showHillshade: boolean) => {
      const inflate = CELL_OVERHANG;
      let o = 0;
      for (const cell of world.cells) {
        const c = getCellColor(cell, viewMode, {
          seaLevel, factionColors,
          seasonalDelta: seasonalTemperatureDelta(cell, world.params),
        });
        if (showHillshade) c.multiplyScalar(shadeMap[cell.id]);
        const hMult = displayRadius(cell.height, smoothGlobe);
        const uc = cell.center;
        const cx = uc.x * hMult, cy = uc.y * hMult, cz = uc.z * hMult;
        for (let i = 0; i < cell.vertices.length; i++) {
          const v1 = cell.vertices[i]; const v2 = cell.vertices[(i + 1) % cell.vertices.length];
          pos[o] = cx; pos[o + 1] = cy; pos[o + 2] = cz;
          pos[o + 3] = (uc.x + (v1.x - uc.x) * inflate) * hMult;
          pos[o + 4] = (uc.y + (v1.y - uc.y) * inflate) * hMult;
          pos[o + 5] = (uc.z + (v1.z - uc.z) * inflate) * hMult;
          pos[o + 6] = (uc.x + (v2.x - uc.x) * inflate) * hMult;
          pos[o + 7] = (uc.y + (v2.y - uc.y) * inflate) * hMult;
          pos[o + 8] = (uc.z + (v2.z - uc.z) * inflate) * hMult;
          col[o] = c.r; col[o + 1] = c.g; col[o + 2] = c.b;
          col[o + 3] = c.r; col[o + 4] = c.g; col[o + 5] = c.b;
          col[o + 6] = c.r; col[o + 7] = c.g; col[o + 8] = c.b;
          o += 9;
        }
      }
    };

    say(`\n--- geometry refill (${triCount} tris) ---`);
    bench('full refill (baseline)', 20, () => refillBaseline('biome', false));
    bench('full refill + hillshade (baseline)', 20, () => refillBaseline('biome', true));
    bench('colour math only', 20, () => {
      for (const cell of world.cells) {
        const c = getCellColor(cell, 'biome', { seaLevel, factionColors, seasonalDelta: seasonalTemperatureDelta(cell, world.params) });
        if (c.r < -1) throw new Error('x');
      }
    });
    bench('seasonalDelta only', 20, () => {
      let s = 0;
      for (const cell of world.cells) s += seasonalTemperatureDelta(cell, world.params);
      if (s === 12345.6789) throw new Error('x');
    });
    bench('computeShadeMap', 20, () => { computeShadeMap(world.cells, seaLevel); });

    // -----------------------------------------------------------------------
    // CANDIDATES
    // -----------------------------------------------------------------------
    say('\n--- candidate: fused MVP projector ---');
    const mvp = new THREE.Matrix4();
    const gmInv = new THREE.Matrix4();
    const camLocal = new THREE.Vector3();
    const prepare = () => {
      const gm = globe.matrixWorld;
      mvp.multiplyMatrices(camera.projectionMatrix, camera.matrixWorldInverse).multiply(gm);
      gmInv.copy(gm).invert();
      camera.getWorldPosition(camLocal).applyMatrix4(gmInv);
    };
    const projectorFused = (x: number, y: number, z: number, out: [number, number]) => {
      let dx = camLocal.x - x, dy = camLocal.y - y, dz = camLocal.z - z;
      const len = Math.hypot(dx, dy, dz);
      if (len > 1e-9) {
        dx /= len; dy /= len; dz /= len;
        const pl = Math.hypot(x, y, z);
        if (pl > 1e-9 && (dx * x + dy * y + dz * z) / pl <= 0.005) return false;
      }
      const e = mvp.elements;
      const w = e[3] * x + e[7] * y + e[11] * z + e[15];
      out[0] = ((e[0] * x + e[4] * y + e[8] * z + e[12]) / w + 1) / 2 * cssW;
      out[1] = (1 - ((e[1] * x + e[5] * y + e[9] * z + e[13]) / w + 1) / 2) * cssH;
      return true;
    };

    bench(`project all ${nCells} cells (fused)`, FRAMES, (i) => {
      advance(i); prepare();
      const cx = camLocal.x, cy = camLocal.y, cz = camLocal.z;
      const e = mvp.elements;
      for (let k = 0; k < nCells; k++) {
        const cell = world.cells[k];
        const c = cell.center;
        const r = displayRadius(cell.height, smoothGlobe);
        const x = c.x * r, y = c.y * r, z = c.z * r;
        let dx = cx - x, dy = cy - y, dz = cz - z;
        const len = Math.hypot(dx, dy, dz);
        const pl = Math.hypot(x, y, z);
        if (len > 1e-9 && pl > 1e-9 && (dx * x + dy * y + dz * z) / (len * pl) <= 0.005) {
          proj.visible[k] = 0; continue;
        }
        const w = e[3] * x + e[7] * y + e[11] * z + e[15];
        proj.x[k] = ((e[0] * x + e[4] * y + e[8] * z + e[12]) / w + 1) / 2 * cssW;
        proj.y[k] = (1 - ((e[1] * x + e[5] * y + e[9] * z + e[13]) / w + 1) / 2) * cssH;
        proj.visible[k] = 1;
      }
    });

    const noopProject = (_x: number, _y: number, _z: number, out: [number, number]) => { out[0] = 0; out[1] = 0; return true; };
    say('\n--- graticule cost split ---');
    bench('graticule w/ no-op projector (walk cost)', FRAMES, (i) => { advance(i); drawGraticuleTenant(ctx, proj, world, noopProject, smoothGlobe); });
    bench('graticule w/ fused projector', FRAMES, (i) => { advance(i); prepare(); drawGraticuleTenant(ctx, proj, world, projectorFused, smoothGlobe); });
    bench('graticule smooth (no walk) baseline proj', FRAMES, (i) => { advance(i); camera.getWorldPosition(camPos); drawGraticuleTenant(ctx, proj, world, projectorBaseline, true); });
    bench('graticule smooth (no walk) fused proj', FRAMES, (i) => { advance(i); prepare(); drawGraticuleTenant(ctx, proj, world, projectorFused, true); });

    say('\n--- tenants with the fused projector ---');
    bench('contours (fused)', FRAMES, (i) => { advance(i); prepare(); drawContoursTenant(ctx, segs, projectorFused, smoothGlobe); });
    bench('borders (fused)', FRAMES, (i) => { advance(i); prepare(); drawBordersTenant(ctx, borders, projectorFused, smoothGlobe); });
    bench('rivers (fused)', FRAMES, (i) => { advance(i); prepare(); drawRiversTenant(ctx, rivers, projectorFused, smoothGlobe); });
    advance(0); prepare();
    bench('currents (fused)', FRAMES, () => { drawCurrentsTenant(ctx, proj, world, projectorFused, smoothGlobe); });

    // Accuracy of the fused projector vs the baseline, in screen pixels.
    advance(7); prepare(); camera.getWorldPosition(camPos);
    let maxErr = 0, flips = 0, tested = 0;
    const a: [number, number] = [0, 0]; const b: [number, number] = [0, 0];
    for (let k = 0; k < nCells; k++) {
      const cell = world.cells[k]; const c = cell.center;
      const r = displayRadius(cell.height, smoothGlobe);
      const va = projectorBaseline(c.x * r, c.y * r, c.z * r, a);
      const vb = projectorFused(c.x * r, c.y * r, c.z * r, b);
      if (va !== vb) { flips++; continue; }
      if (!va) continue;
      tested++;
      maxErr = Math.max(maxErr, Math.abs(a[0] - b[0]), Math.abs(a[1] - b[1]));
    }
    say(`\nfused vs baseline: max screen error ${maxErr.toExponential(2)} px over ${tested} pts, ${flips} visibility flips`);

    // ---- candidate: staged typed arrays + sqrt instead of hypot ----
    const pts = new Float32Array(nCells * 3);
    const plen = new Float32Array(nCells);
    const stage = () => {
      for (let k = 0; k < nCells; k++) {
        const c = world.cells[k].center;
        const r = displayRadius(world.cells[k].height, smoothGlobe);
        const x = c.x * r, y = c.y * r, z = c.z * r;
        pts[k * 3] = x; pts[k * 3 + 1] = y; pts[k * 3 + 2] = z;
        plen[k] = Math.sqrt(x * x + y * y + z * z);
      }
    };
    say('\n--- candidate: staged typed arrays ---');
    bench('stage() rebuild (per world, not per frame)', 20, () => stage());
    stage();
    bench(`project all ${nCells} cells (staged+fused+sqrt)`, FRAMES, (i) => {
      advance(i); prepare();
      const cx = camLocal.x, cy = camLocal.y, cz = camLocal.z;
      const e = mvp.elements;
      const e0 = e[0], e4 = e[4], e8 = e[8], e12 = e[12];
      const e1 = e[1], e5 = e[5], e9 = e[9], e13 = e[13];
      const e3 = e[3], e7 = e[7], e11 = e[11], e15 = e[15];
      for (let k = 0; k < nCells; k++) {
        const j = k * 3;
        const x = pts[j], y = pts[j + 1], z = pts[j + 2];
        const dx = cx - x, dy = cy - y, dz = cz - z;
        const len = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if ((dx * x + dy * y + dz * z) / (len * plen[k]) <= 0.005) { proj.visible[k] = 0; continue; }
        const w = e3 * x + e7 * y + e11 * z + e15;
        proj.x[k] = ((e0 * x + e4 * y + e8 * z + e12) / w + 1) / 2 * cssW;
        proj.y[k] = (1 - ((e1 * x + e5 * y + e9 * z + e13) / w + 1) / 2) * cssH;
        proj.visible[k] = 1;
      }
    });

    // ---- candidate: sqrt-based fused projector (no Math.hypot) ----
    const projectorSqrt = (x: number, y: number, z: number, out: [number, number]) => {
      const dx = camLocal.x - x, dy = camLocal.y - y, dz = camLocal.z - z;
      const len = Math.sqrt(dx * dx + dy * dy + dz * dz);
      const pl = Math.sqrt(x * x + y * y + z * z);
      if (len < 1e-9 || pl < 1e-9) { /* degenerate: visible */ }
      else if ((dx * x + dy * y + dz * z) / (len * pl) <= 0.005) return false;
      const e = mvp.elements;
      const w = e[3] * x + e[7] * y + e[11] * z + e[15];
      out[0] = ((e[0] * x + e[4] * y + e[8] * z + e[12]) / w + 1) / 2 * cssW;
      out[1] = (1 - ((e[1] * x + e[5] * y + e[9] * z + e[13]) / w + 1) / 2) * cssH;
      return true;
    };
    say('\n--- candidate: sqrt projector on tenants ---');
    bench('contours (sqrt)', FRAMES, (i) => { advance(i); prepare(); drawContoursTenant(ctx, segs, projectorSqrt, smoothGlobe); });
    bench('borders (sqrt)', FRAMES, (i) => { advance(i); prepare(); drawBordersTenant(ctx, borders, projectorSqrt, smoothGlobe); });
    bench('rivers (sqrt)', FRAMES, (i) => { advance(i); prepare(); drawRiversTenant(ctx, rivers, projectorSqrt, smoothGlobe); });
    bench('graticule smooth (sqrt)', FRAMES, (i) => { advance(i); prepare(); drawGraticuleTenant(ctx, proj, world, projectorSqrt, true); });
    advance(0); prepare();
    bench('currents (sqrt)', FRAMES, () => { drawCurrentsTenant(ctx, proj, world, projectorSqrt, smoothGlobe); });

    // sqrt-vs-hypot accuracy
    advance(11); prepare(); camera.getWorldPosition(camPos);
    let e2 = 0, f2 = 0;
    for (let k = 0; k < nCells; k++) {
      const c = world.cells[k].center; const r = displayRadius(world.cells[k].height, smoothGlobe);
      const va = projectorBaseline(c.x * r, c.y * r, c.z * r, a);
      const vb = projectorSqrt(c.x * r, c.y * r, c.z * r, b);
      if (va !== vb) { f2++; continue; }
      if (!va) continue;
      e2 = Math.max(e2, Math.abs(a[0] - b[0]), Math.abs(a[1] - b[1]));
    }
    say(`sqrt vs baseline: max screen error ${e2.toExponential(2)} px, ${f2} visibility flips`);
  }
}
