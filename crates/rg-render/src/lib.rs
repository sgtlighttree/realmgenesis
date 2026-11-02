//! GPU rendering. Takes `&wgpu::Device` / `&wgpu::Queue` and never sees eframe or egui.
//!
//! M0/MVP: one preview mesh of the whole planet plus equirectangular layer textures.
//! M1 replaces the mesh with the streamed tile quadtree.

pub mod camera;
pub mod mesh;
pub mod planet;

pub use camera::OrbitCamera;
pub use mesh::{PlanetMesh, Vertex};
pub use planet::{HeightImage, LayerImage, MapView, PlanetRenderer, RenderStyle, ViewMode};
