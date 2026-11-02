#!/usr/bin/env bash
# Build the web bundle and lay it out for a claude.ai artifact:
#   <out>/index.html  (body-only page; the artifact host adds <html>/<head>)
#   <out>/app.js, <out>/app_bg.wasm  (stable names, published via `files`)
# Usage: scripts/build-web-artifact.sh <out_dir>
set -euo pipefail
OUT="${1:?usage: build-web-artifact.sh <out_dir>}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DIST="$(mktemp -d)"
(cd "$ROOT/crates/rg-app" && trunk build --release=true --dist "$DIST" --public-url ./ index.html)
mkdir -p "$OUT"
cp "$DIST"/realmgenesis-*.js "$OUT/app.js"
cp "$DIST"/realmgenesis-*_bg.wasm "$OUT/app_bg.wasm"
REV="$(git -C "$ROOT" rev-parse --short HEAD)"
cat > "$OUT/index.html" <<HTML
<title>RealmGenesis</title>
<style>
  html, body { margin: 0; height: 100%; overflow: hidden; background: #0e1114; color: #dae0e4; font: 14px system-ui, sans-serif; }
  canvas#realmgenesis { position: absolute; inset: 0; width: 100%; height: 100%; display: block; }
  #loading { position: absolute; inset: 0; display: grid; place-content: center; text-align: center; padding: 16px; }
  #loading code { color: #6eba9a; }
</style>
<canvas id="realmgenesis"></canvas>
<div id="loading"><p>Loading RealmGenesis <code>${REV}</code>…</p></div>
<script type="module">
  try {
    const m = await import('./app.js');
    await m.default({ module_or_path: './app_bg.wasm' });
  } catch (e) {
    document.getElementById('loading').innerHTML =
      '<p>RealmGenesis failed to load.</p><p><code>' + String(e).replace(/</g, '&lt;') + '</code></p>';
  }
</script>
HTML
rm -rf "$DIST"
ls -la "$OUT"
