//! GPU rendering. Takes `&wgpu::Device` / `&wgpu::Queue` and never sees eframe or egui.
//!
//! M0: one preview mesh of the whole planet. M1 replaces it with the streamed tile quadtree.

pub mod camera;
pub mod globe;
pub mod mesh;

pub use camera::OrbitCamera;
pub use globe::{GlobeRenderer, GlobeStyle};
pub use mesh::{PlanetMesh, Vertex};
