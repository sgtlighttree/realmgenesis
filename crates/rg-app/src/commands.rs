//! Everything the UI can ask the app to do. The UI never mutates world state directly: it
//! pushes `Command`s, and `RealmApp::apply` executes them. This seam is what lets the UI be
//! redesigned (or swapped for another toolkit) without touching the core.

#[derive(Clone, Debug)]
pub enum Command {
    /// Rebuild the planet from the current parameters.
    Regenerate,
    /// Pick a new random seed, then regenerate.
    RandomizeSeed,
    /// Reset the camera to the default orbit.
    ResetCamera,
}
