# 🌏 RealmGenesis 3D

**RealmGenesis** is a powerful, browser-based procedurally generated fantasy world engine. It simulates tectonic plates, hydraulic erosion, moisture transport, and biomes on a spherical 3D globe.

## ✨ Key Features

- **🌋 Core Simulations**: Advanced Voronoi-based world generation with tectonic plate movement, volcanic activity, and realistic hydraulic erosion.
- **❄️ Climate Engine**: Dynamic moisture transport and temperature variance based on axial tilt and latitude, resulting in procedurally accurate biomes.
- **🗺️ Multiple Projections**: 
    - **3D Globe**: Interactive orbital viewer.
    - **2D Mercator**: Classic flat map projection.
    - **Experimental Dymaxion**: High-resolution icosahedral projection with Sharp-DPI support and interactive orientation.
- **🤖 AI Lore (Gemini)**: Detailed world lore, faction backstories, and capital naming powered by Google Gemini Flash.
- **🔑 BYOK (Bring Your Own Key)**: Use your own Gemini API key for AI features. Keys are ephemeral and never stored permanently in the app.
- **⚖️ Political Simulation**: Procedural expansion of factions across provinces and towns with customizable borders and tunable country-size variance.
- **📈 Analytics View Layers**: Population density heatmap and province-level administrative view, alongside biome, satellite, height, climate, plate, and political layers.
- **Interactive Map Editing**: Paint terrain, biomes, and political borders directly in 3D, Mercator, or Dymaxion views, with undo support.
- **Province-Aware Political Brush**: Political painting preserves province ownership, supports an unclaim/eraser brush, and keeps faction population totals in sync.
- **Toggleable Faction Overlay**: Faction borders and names can be shown over any view layer; 3D names curve along the globe surface while 2D labels render over Mercator and Dymaxion maps.
- **📦 3D Export (GLB)**: Export the world as a binary GLB file — world mesh with per-vertex colors, rivers as line geometry, and city markers as low-poly cylinders. Ready to import into Blender, Unity, or Unreal.
- **🔲 Blender UV Net**: Dymaxion export mode that produces a square image matching Blender's default icosphere UV unwrap, compatible with any subdivision level.
- **📊 World Stats Panel**: Live stats in the Sys tab — land coverage %, top biome breakdown, total population, faction/province counts, and river count.
- **⌨️ Keyboard Shortcuts**: `G` to generate, `R` to randomize seed, `Escape` to close the inspector, `1`–`5` to switch tabs. Guards against firing inside text inputs.
- **🗺️ Heightmap PNG Export**: One-click greyscale equirectangular export for use as a displacement map in any 3D application.

## 🚀 Getting Started

### Prerequisites
- [Node.js](https://nodejs.org/) (v18+)
- [npm](https://www.npmjs.com/)

### Running Locally

1. **Clone & Install**:
   ```bash
   npm install
   ```

2. **Configure (Optional)**:
   Create a `.env.local` file and add your `GEMINI_API_KEY` for default AI support, or provide it at runtime in the app settings.

3. **Launch**:
   ```bash
   npm run dev
   ```

## 📖 Documentation

- [**docs/**](./docs/README.md) — The settled technical reference: architecture, generation pipeline, V3 tectonics, data model, params reference, rendering, civilization layers, export, invariants, testing, and engineering notes. Start at `docs/README.md`.
- [**HANDOFF.md**](./HANDOFF.md) — Recent-session state and open work. [**ROADMAP.md**](./ROADMAP.md) — planned features.
- [**AGENTS.md**](./AGENTS.md) — Commands, code style guide, and conventions for contributors and AI agents.

## 🌐 Deployment

RealmGenesis is optimized for deployment on **Netlify**. 
- **SPA Features**: Includes `_redirects` support for deep-linking in browser storage flows.
- **Environmental Safety**: Safe fallback mechanisms for environment variables.

---
*Built with React, Three.js, OpenAI Codex, Google Gemini, and Claude Code.*
