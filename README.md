<!-- <div align="center">
<img width="1200" height="475" alt="GHBanner" src="https://github.com/user-attachments/assets/0aa67016-6eaf-458a-adb2-6e31a0763ed6" />
</div> -->

# 🌏 RealmGenesis 3D

**RealmGenesis** is a powerful, browser-based procedurally generated fantasy world engine. It simulates tectonic plates, hydraulic erosion, moisture transport, and biomes on a spherical 3D globe.

## ✨ Key Features

- **🌋 Core Simulations**: Advanced Voronoi-based world generation with tectonic plate movement, volcanic activity, and realistic hydraulic erosion.
- **❄️ Climate Engine**: Dynamic moisture transport and temperature variance based on axial tilt and latitude, resulting in procedurally accurate biomes.
- **🗺️ Multiple Projections**: 
    - **3D Globe**: Interactive orbital viewer.
    - **2D Mercator**: Classic flat map projection.
    - **Experimental Dymaxion**: High-resolution icosahedral projection with Sharp-DPI support and interactive orientation.
- **🤖 AI Lore (Gemini)**: Detailed world lore, faction backstories, and capital naming powered by Google Gemini 1.5 Flash.
- **🔑 BYOK (Bring Your Own Key)**: Use your own Gemini API key for AI features. Keys are ephemeral and never stored permanently in the app.
- **⚖️ Political Simulation**: Procedural expansion of factions across provinces and towns with customizable borders.

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

## 🌐 Deployment

RealmGenesis is optimized for deployment on **Netlify**. 
- **SPA Features**: Includes `_redirects` support for deep-linking in browser storage flows.
- **Environmental Safety**: Safe fallback mechanisms for environment variables.

---
*Built with React, Three.js, OpenAI Codex, Google Gemini.*
