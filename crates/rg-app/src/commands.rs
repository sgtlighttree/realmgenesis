//! Everything the UI can ask the app to do. The UI never mutates world state directly: it
//! pushes `Command`s, and `RealmApp::apply` executes them. This seam is what lets the UI be
//! redesigned (or swapped for another toolkit) without touching the core.

use rg_pipeline::Layer;
use rg_render::ViewMode;

#[derive(Clone, Debug)]
pub enum Command {
    /// Rebuild the planet from the current parameters (tectonics + climate).
    Regenerate,
    /// Pick a new random seed, then regenerate.
    RandomizeSeed,
    /// Re-run only the climate with the current planet parameters.
    RerunClimate,
    SetLayer(Layer),
    SetView(ViewMode),
    ResetCamera,
    /// Show or hide the controls sidebar.
    TogglePanel,
    /// Stop or resume background terrain-tile generation. Stopping also drops every
    /// queued tile, so the CPU goes quiet within one tile's time.
    ToggleTileStreaming,
    /// Save the world file (native only for now).
    SaveWorld,
    /// Load the most recent world file.
    LoadWorld,
    /// Write a 16-bit equirectangular heightmap PNG.
    ExportHeightmap,
}
