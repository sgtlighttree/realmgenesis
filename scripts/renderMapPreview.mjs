// Render a real map export to SVG (and optionally to PNG) so an agent can LOOK
// at it.
//
// This is the instrument that found the two A3 bugs 352 passing unit tests could
// not: coastlines collapsed to a point (Point3 vectors passed where d3 wanted
// [lon,lat] degrees) and a full-opacity seam stroke that rendered the whole map
// as a dark honeycomb. Every piece was correct in isolation; only the
// composition was broken, and only pixels showed it.
//
// It drives the REAL `exportSVG` path — no dev server, no R3F, no browser
// automation against the app, and none of the auto-rotate CPU trap that makes
// driving the 3D globe expensive on an M1.
//
// Usage:
//   node scripts/renderMapPreview.mjs [--style=parchment] [--mode=biome]
//                                     [--seed=S] [--points=N] [--width=N]
//                                     [--out=DIR] [--png]
//
// --style may be a comma list (e.g. parchment,default) to render several at once.
// --png additionally screenshots each SVG with headless Chromium, if Playwright
//   is available in the npx cache. Read the PNG afterwards and actually look at it.
//
// Defaults: style=parchment mode=biome seed=realmgenesis points=12000 width=1400
//           out=tmp
import { createServer } from 'vite';
import { writeFileSync, mkdirSync, existsSync, readdirSync } from 'fs';

const flags = Object.fromEntries(
  process.argv.slice(2).filter(a => a.startsWith('--')).map(a => {
    const [k, v] = a.replace(/^--/, '').split('=');
    return [k, v ?? true];
  }),
);

const styles = String(flags.style ?? 'parchment').split(',');
const modes = String(flags.mode ?? 'biome').split(',');
const seed = flags.seed ?? 'realmgenesis';
const points = Number(flags.points ?? 12000);
const width = Number(flags.width ?? 1400);
const outDir = String(flags.out ?? 'tmp');

if (!existsSync(outDir)) mkdirSync(outDir, { recursive: true });

const server = await createServer({ server: { middlewareMode: true } });
const { generateWorld } = await server.ssrLoadModule('/utils/worldGen.ts');
const { makeParams } = await server.ssrLoadModule('/tests/helpers.ts');
const { exportSVG } = await server.ssrLoadModule('/utils/exportVector.ts');

console.log(`generating: seed=${seed} points=${points}`);
const world = await generateWorld(makeParams({ seed, points }));

const written = [];
for (const style of styles) {
  for (const mode of modes) {
    const svg = exportSVG(world, mode, 'equirectangular', width, style);
    const file = `${outDir}/${style}-${mode}.svg`;
    writeFileSync(file, svg);
    written.push(file);
    console.log(`  ${file}  ${(svg.length / 1024).toFixed(0)}KB`);
  }
}
await server.close();

if (flags.png) {
  // Playwright is not a project dependency. Reuse whatever the npx cache holds
  // and the shared browsers root from ~/.zshrc, per the house Playwright rules.
  const npxRoot = `${process.env.HOME}/.npm/_npx`;
  let playwrightPath = null;
  if (existsSync(npxRoot)) {
    for (const dir of readdirSync(npxRoot)) {
      const candidate = `${npxRoot}/${dir}/node_modules/playwright/index.mjs`;
      if (existsSync(candidate)) { playwrightPath = candidate; break; }
    }
  }
  if (!playwrightPath) {
    console.log('\n--png skipped: no playwright in the npx cache.');
    console.log('Populate it with:  npx -y playwright@latest --version');
    process.exit(0);
  }

  const cacheRoot = `${process.env.HOME}/Library/Caches/ms-playwright`;
  const shell = readdirSync(cacheRoot)
    .filter(n => n.startsWith('chromium_headless_shell-'))
    .sort()
    .pop();
  const executablePath =
    `${cacheRoot}/${shell}/chrome-headless-shell-mac-arm64/chrome-headless-shell`;

  const { chromium } = await import(playwrightPath);
  const browser = await chromium.launch({ executablePath });
  const page = await browser.newPage({ viewport: { width, height: Math.round(width / 2) } });
  for (const file of written) {
    await page.goto(`file://${process.cwd()}/${file}`);
    await page.waitForTimeout(500);
    const png = file.replace(/\.svg$/, '.png');
    await page.screenshot({ path: png });
    console.log(`  ${png}`);
  }
  await browser.close();
  console.log('\nNow READ the PNGs. Gates do not see a broken map.');
}

process.exit(0);
