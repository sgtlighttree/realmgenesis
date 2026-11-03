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

## Local (M1) sessions (Matt, 2026-11-03)

Applies only on Matt's MacBook Air M1 (4 P + 4 E cores, 16 GB, fanless; also his After Effects / Blender workstation). Cloud sessions skip this section.

- **Capped is the default** for every build, test, clippy and probe run. The engine sizes its own thread pools from `available_parallelism()`, so `-j` and `--test-threads` alone don't cap it; `taskpolicy -b` (background QoS, E-cores only) does (s02: 2.2 vs 6.8 cores measured).
  ```
  taskpolicy -b cargo build -p rg-app --release -j 4
  taskpolicy -b cargo test --workspace --release -j 4 -- --test-threads=2
  taskpolicy -b cargo clippy --workspace --all-targets --release -j 4
  ```
- **Full speed** (drop `taskpolicy -b` and `-j 4`) only when Matt says the Mac is free. Offer it when After Effects, aerender and Blender aren't using CPU: check `top -l 2 -s 2 -o cpu` (open ≠ busy). If one is busy, ask first or drop to `-j 2 -- --test-threads=1`.
- **Never cap the interactive app:** build capped, then run `target/release/realmgenesis` directly. On E-cores the UI and tile workers would misrepresent performance.
- **Timing runs (`bench_tiles`, probes' ms figures) go uncapped.** E-cores make tiles look ~2–5× slower. Warm up capped once (caches the world in `target/probe-cache/`), then rerun without `taskpolicy`; the timed tile loop is single-threaded, so it takes ~1 P-core. Measured capped (s02): build 11 min, `cargo test --workspace --release` 14 min (`terrain_realism` 232 s), clippy 2 min; peak RSS 1.7 GB.
- **RAM has no hard cap on macOS** (`ulimit -v` fails, no cgroups). Concurrency is the lever; measure peaks with `/usr/bin/time -l`.
- **Release profile for everything**, tests and clippy included: a debug-profile run builds a second dependency tree, and disk is tight.
- **Low data (metered link):** the toolchain is rustup's minimal profile plus `rustup override set 1.97.0` (no wasm32 until a web build; `rustup override unset` restores the toml). Builds pass `--offline`; after a dependency change run `cargo fetch --target aarch64-apple-darwin` first.

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
