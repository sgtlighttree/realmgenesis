# RealmGenesis v2 — project rules

A Rust + wgpu fantasy planet generator. Start with `HANDOFF.md`, then the newest file in `docs/session-logs/`, then `docs/ARCHITECTURE.md`.

## Session workflow (Matt, 2026-11-02, overrides the global "push only when asked")

- **Push early and often**, for any purpose. Sessions run in ephemeral cloud containers, and anything unpushed dies with them. Commit in small topical chunks and push after each.
- **HANDOFF.md + session log, written continuously.** After every decision, finding or blocker, append it to the current `docs/session-logs/YYYY-MM-DD-sNN-<slug>.md` and update `HANDOFF.md` (the current state and next steps). Then commit and push. A fresh session must be able to continue from these files alone.
- **Second opinion:** at commitment boundaries (architecture, schema, a refactor across 3+ files, the second failed fix), run a fresh-context reviewer agent as the advisor stand-in, and log its verdict. There is no `advisor` tool in cloud sessions.
- **Dependencies:** assume your knowledge is stale. Check current versions (`cargo search`, docs.rs) before writing against an API. There's no Context7 in cloud sessions.

## Build

```
cargo run -p rg-app --release          # native
trunk serve crates/rg-app/index.html   # web, http://127.0.0.1:8080
trunk build --release crates/rg-app/index.html --dist dist   # web bundle for the artifact preview
cargo test --workspace
cargo clippy --workspace --all-targets
```

## Code rules

- Core crates (`rg-geo`, `rg-noise`, the sim crates, `rg-world`, `rg-render`) never depend on egui or eframe. The UI talks to the core through `Command`s.
- Anything **committed** is stored in the world file, never re-derived (see the determinism contract in ARCHITECTURE §5).
- GPU-first simulations. CPU fallbacks must run single-threaded on wasm32. `rayon` is native-only.
- Positions are f64 on the CPU and camera-relative f32 on the GPU.
- Matt is new to Rust. Keep code idiomatic but plain: no clever trait towers or macro magic without a reason written next to it.

## Preview loop

1. Build the web bundle, then screenshot it with headless Chromium (Playwright is preinstalled; no `playwright install`).
2. Publish the bundle to the project's preview artifact (URL in `HANDOFF.md`). Republish to the same URL.
3. CI (`.github/workflows/`) builds native binaries on push.
