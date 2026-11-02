# RealmGenesis v2

A fantasy planet generator, native and fast: a 3D globe and a flat map, plate tectonics, a one-shot climate and ocean run, height detail generated as you zoom, sculpting, and editable political borders. Think Google Earth, but for a world that doesn't exist.

Rust + wgpu + egui. Runs natively on macOS, Windows and Linux, and in WebGPU browsers.

**Status:** early. See `docs/ARCHITECTURE.md` for the plan and `HANDOFF.md` for where things stand.

The V1 browser app (TypeScript, Voronoi cells) lives on `main`.

## Run it

Install Rust once: https://rustup.rs. Then:

```bash
cargo run -p rg-app --release
```
