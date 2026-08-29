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
//
// Writes <out>/globe-<label>-lat<NN>.png. Read the PNGs afterwards and look.
//
// Uses its own throwaway Vite server on an ephemeral port; it never touches a
// dev server you already have running.
import { createServer } from 'vite';
import { mkdirSync, existsSync, readdirSync, statSync } from 'fs';

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

const server = await createServer({ server: { port: 0, strictPort: false } });
await server.listen();
const port = server.config.server.port ?? server.httpServer.address().port;
const base = `http://localhost:${port}`;

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
  + `&seed=${encodeURIComponent(seed)}&points=${points}&size=${size}`;
console.log(`loading ${url}`);
await page.goto(url, { waitUntil: 'load' });
await page.waitForFunction('window.__ready === true', null, { timeout: 180000 });

for (const lat of views) {
  await page.evaluate(l => window.__look(l), lat);
  const file = `${outDir}/globe-${label}-lat${lat}.png`;
  await page.locator('#c').screenshot({ path: file });
  console.log(`  ${file}`);
}

await browser.close();
await server.close();
console.log('\nNow LOOK at the PNGs. Test counts cannot see this class of bug.');
