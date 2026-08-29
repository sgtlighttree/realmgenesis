// Screenshot the REAL styled 3D globe from fixed camera angles, so an agent can
// LOOK at it.
//
// `renderMapPreview.mjs` renders the flat SVG export and cannot see globe
// defects at all — the polar swirl that survived three fix attempts was
// invisible to it. This drives `scripts/preview/globe.html`, which builds the
// same geometry, the same baked texture and the same material as `WorldViewer`,
// with no app UI and no OrbitControls, so the camera can sit exactly on a pole.
//
// Usage:
//   node scripts/renderGlobePreview.mjs [--style=parchment] [--mode=biome]
//        [--seed=S] [--points=N] [--size=N] [--out=DIR] [--label=NAME]
//        [--views=90,45,0]        (camera latitudes, degrees)
//        [--lon=0] [--dist=2.6]   (camera longitude; lon 180 is the SEAM)
//        [--texture]              (also dump the flat equirectangular bake)
//        [--port=3000]            (dev server to reuse; 0 forces a private one)
//        [--marker]               (paint the canvas top RED / bottom BLUE, to
//                                  see which pole each edge lands on)
//        [--hatch=0]              (ocean hatch angle; 0 = horizontal, which
//                                  makes any offset or skew read off a grid)
//
// --texture is how you tell a TEXTURE defect from a SAMPLING one: if the mark
// is in the flat bake it is a style-pass bug, not a globe bug.
//
// Writes <out>/globe-<label>-lat<NN>.png. Read the PNGs afterwards and look.
//
// It REUSES the dev server on port 3000 if one is already up, and only starts a
// private one when there is none. Vite serves any file under the repo root, so
// the harness page loads straight off it. A second Vite is a second full module
// graph on a laptop that may also be rendering — and starting one while yours
// was running is exactly what made the test suite time out. It never kills a
// server it did not start.
import { createServer } from 'vite';
import { mkdirSync, existsSync, readdirSync, statSync, writeFileSync } from 'fs';

const flags = Object.fromEntries(
  process.argv.slice(2).filter(a => a.startsWith('--')).map(a => {
    const [k, v] = a.replace(/^--/, '').split('=');
    return [k, v ?? true];
  }),
);

const style = String(flags.style ?? 'parchment');
const mode = String(flags.mode ?? 'biome');
const seed = String(flags.seed ?? 'realmgenesis');
const points = Number(flags.points ?? 8000);
const size = Number(flags.size ?? 900);
const outDir = String(flags.out ?? 'tmp');
const label = String(flags.label ?? style);
const views = String(flags.views ?? '90,45,0').split(',').map(Number);
// lon 180 puts the antimeridian dead centre. The fragment shader's atan wraps
// there, so it is the one place the per-fragment mapping has a failure mode the
// per-vertex one did not — check it, do not reason about it.
const lon = Number(flags.lon ?? 0);
const dist = Number(flags.dist ?? 2.6);

if (!existsSync(outDir)) mkdirSync(outDir, { recursive: true });

// Playwright is not a project dependency. Reuse whatever the npx cache holds
// and the shared browsers root from ~/.zshrc, per the house Playwright rules.
// FULL chromium, not chrome-headless-shell: this needs a real WebGL context.
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
// The .app bundle is named differently across Playwright builds and CPU
// architectures ("Chromium.app", "Google Chrome for Testing.app", chrome-mac vs
// chrome-mac-arm64), so find it rather than spelling one layout out.
const executablePath = (() => {
  const root = `${cacheRoot}/${build}`;
  for (const arch of readdirSync(root)) {
    const dir = `${root}/${arch}`;
    if (!statSync(dir).isDirectory()) continue;
    for (const app of readdirSync(dir).filter(n => n.endsWith('.app'))) {
      const macOs = `${dir}/${app}/Contents/MacOS`;
      if (!existsSync(macOs)) continue;
      const bin = readdirSync(macOs).find(n => !n.includes('.'));
      if (bin) return `${macOs}/${bin}`;
    }
  }
  return null;
})();
if (!executablePath) {
  console.error(`no chromium executable under ${cacheRoot}/${build}`);
  console.error('install it with:  npx playwright install chromium');
  process.exit(1);
}
const { chromium } = await import(playwrightPath);

// Reuse the running dev server if there is one. --port=0 forces a private
// server (use it when the dev server is mid-restart or serving another branch).
const wanted = Number(flags.port ?? 3000);
const probe = async (url) => {
  try {
    const res = await fetch(`${url}/scripts/preview/globe.html`, { method: 'HEAD' });
    return res.ok;
  } catch {
    return false;
  }
};

let base = wanted ? `http://localhost:${wanted}` : '';
let server = null;
if (base && await probe(base)) {
  console.log(`reusing the dev server on ${base}`);
} else {
  server = await createServer({ server: { port: 0, strictPort: false } });
  await server.listen();
  const port = server.config.server.port ?? server.httpServer.address().port;
  base = `http://localhost:${port}`;
  console.log(`no dev server on ${wanted}; started a private one on ${port}`);
}

const browser = await chromium.launch({
  executablePath,
  // Headless Chromium falls back to SwiftShader for WebGL; these keep it from
  // refusing the context outright.
  args: ['--use-gl=angle', '--use-angle=swiftshader', '--enable-unsafe-swiftshader'],
});
const page = await browser.newPage({ viewport: { width: size, height: size } });
page.on('console', m => { if (m.type() === 'error') console.error('  page error:', m.text()); });
page.on('pageerror', e => console.error('  page error:', e.message));

const url = `${base}/scripts/preview/globe.html?style=${style}&mode=${mode}`
  + `&seed=${encodeURIComponent(seed)}&points=${points}&size=${size}`
  + (flags.marker ? '&marker=1' : '')
  + (flags.hatch !== undefined ? `&hatch=${flags.hatch}` : '');
console.log(`loading ${url}`);
await page.goto(url, { waitUntil: 'load' });
await page.waitForFunction('window.__ready === true', null, { timeout: 180000 });

for (const lat of views) {
  await page.evaluate(a => window.__look(a[0], a[1], a[2]), [lat, lon, dist]);
  const file = `${outDir}/globe-${label}-lat${lat}.png`;
  await page.locator('#c').screenshot({ path: file });
  console.log(`  ${file}`);
}

if (flags.texture) {
  const dataUrl = await page.evaluate('window.__texture()');
  if (dataUrl) {
    const file = `${outDir}/globe-${label}-texture.png`;
    writeFileSync(file, Buffer.from(dataUrl.split(',')[1], 'base64'));
    console.log(`  ${file}`);
  } else {
    console.log('  (no baked texture — this style draws nothing)');
  }
}

await browser.close();
// Only ever close a server this script started.
if (server) await server.close();
console.log('\nNow LOOK at the PNGs. Test counts cannot see this class of bug.');
