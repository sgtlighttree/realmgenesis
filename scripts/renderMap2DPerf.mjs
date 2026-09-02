// Measure what a Map2D redraw costs, and what a pan/zoom frame costs (F4).
//
// Drives `scripts/preview/map2dPerf.html`, which mounts the REAL Map2D and
// scripts a drag and a wheel gesture against it. Prints a JSON blob plus a
// screenshot of the settled frame, so the same run gives both the timing and
// the pixel evidence that nothing changed.
//
// Usage:
//   node scripts/renderMap2DPerf.mjs [--points=30000] [--style=parchment]
//        [--mode=biome] [--seed=S] [--width=1400] [--height=700]
//        [--projection=mercator] [--dpr=2] [--label=before] [--out=tmp]
//        [--port=3000] [--hillshade=1]
//
// --dpr MATTERS. At deviceScaleFactor 1 the settled quality DPR equals
// INTERACTION_DPR, so `setQualityDpr` bails, the offscreen is never redrawn on
// mousedown, and the gesture redraws this script exists to measure do not
// happen at all. 2 is what the machine actually runs at.
//
// It REUSES the dev server on port 3000 like renderGlobePreview does, and only
// starts a private one when nothing answers. Note the page is loaded through a
// path RELATIVE TO THE SERVER ROOT: run from a worktree and the URL carries the
// worktree prefix, so the modules under test are this checkout's.
import { createServer } from 'vite';
import { mkdirSync, existsSync, readdirSync, writeFileSync } from 'fs';
import { relative, resolve } from 'path';
import { execSync } from 'child_process';

const flags = Object.fromEntries(
  process.argv.slice(2).filter(a => a.startsWith('--')).map(a => {
    const [k, v] = a.replace(/^--/, '').split('=');
    return [k, v ?? true];
  }),
);

const outDir = String(flags.out ?? 'tmp');
const label = String(flags.label ?? 'run');
const dpr = Number(flags.dpr ?? 2);
const width = Number(flags.width ?? 1400);
const height = Number(flags.height ?? 700);
const port = Number(flags.port ?? 3000);

if (!existsSync(outDir)) mkdirSync(outDir, { recursive: true });

const npxRoot = `${process.env.HOME}/.npm/_npx`;
let playwrightPath = null;
if (existsSync(npxRoot)) {
  for (const dir of readdirSync(npxRoot)) {
    const candidate = `${npxRoot}/${dir}/node_modules/playwright/index.mjs`;
    if (existsSync(candidate)) { playwrightPath = candidate; break; }
  }
}
if (!playwrightPath) {
  console.error('no playwright in the npx cache. Populate it with:');
  console.error('  npx -y playwright@latest --version');
  process.exit(1);
}
const cacheRoot = `${process.env.HOME}/Library/Caches/ms-playwright`;
const build = readdirSync(cacheRoot)
  .filter(n => /^chromium-\d+$/.test(n))
  .sort((a, b) => Number(a.split('-')[1]) - Number(b.split('-')[1]))
  .pop();
const chromiumDir = `${cacheRoot}/${build}`;
const inner = readdirSync(chromiumDir).find(n => n.startsWith('chrome-mac'));
const app = readdirSync(`${chromiumDir}/${inner}`).find(n => n.endsWith('.app'));
const executablePath =
  `${chromiumDir}/${inner}/${app}/Contents/MacOS/${app.replace(/\.app$/, '')}`;

// Where this checkout sits relative to whatever root the reused server serves.
// A worktree lives under the main repo, so the prefix is what makes the server
// hand back THIS code rather than main's.
const here = process.cwd();
let prefix = '';
let ownServer = null;
let origin = `http://localhost:${port}`;

const alive = await fetch(origin, { signal: AbortSignal.timeout(1500) })
  .then(r => r.ok || r.status < 500)
  .catch(() => false);

if (alive) {
  const mainRoot = execSync('git rev-parse --path-format=absolute --git-common-dir', { cwd: here })
    .toString().trim().replace(/\/\.git$/, '');
  const rel = relative(mainRoot, here);
  prefix = rel && !rel.startsWith('..') ? `/${rel}` : '';
  console.log(`reusing dev server on ${port}${prefix ? ` (worktree prefix ${prefix})` : ''}`);
} else {
  console.log('no dev server on 3000 — starting a private one');
  ownServer = await createServer({ server: { port: 0, host: '127.0.0.1' } });
  await ownServer.listen();
  origin = `http://localhost:${ownServer.config.server.port}`;
}

const params = new URLSearchParams({
  width: String(width),
  height: String(height),
  points: String(flags.points ?? 30000),
  seed: String(flags.seed ?? 'realmgenesis'),
  style: String(flags.style ?? 'parchment'),
  mode: String(flags.mode ?? 'biome'),
  projection: String(flags.projection ?? 'mercator'),
  hillshade: String(flags.hillshade ?? '1'),
});
const url = `${origin}${prefix}/scripts/preview/map2dPerf.html?${params}`;
console.log(url);

const { chromium } = await import(playwrightPath);
const browser = await chromium.launch({ executablePath });
const page = await browser.newPage({
  viewport: { width: width + 40, height: height + 40 },
  deviceScaleFactor: dpr,
});
page.on('console', m => { if (m.type() === 'error') console.log('  [page]', m.text()); });
page.on('pageerror', e => console.log('  [pageerror]', e.message));

await page.goto(url);
await page.waitForFunction('window.__perfDone === true', null, { timeout: 900000 });
const perf = await page.evaluate('window.__perf');

const shot = `${resolve(outDir)}/map2d-${label}.png`;
await page.locator('#host canvas').screenshot({ path: shot });
await browser.close();
if (ownServer) await ownServer.close();

const json = `${resolve(outDir)}/map2d-${label}.json`;
writeFileSync(json, JSON.stringify(perf, null, 2));

console.log('');
console.log(JSON.stringify(perf, null, 2));
console.log('');
console.log(`  ${shot}`);
console.log(`  ${json}`);
process.exit(0);
