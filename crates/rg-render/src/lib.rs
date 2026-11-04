//! GPU rendering. Takes `&wgpu::Device` / `&wgpu::Queue` and never sees eframe or egui.
//!
//! M0/MVP: one preview mesh of the whole planet plus equirectangular layer textures.
//! M1 replaces the mesh with the streamed tile quadtree.

pub mod camera;
pub mod earth_camera;
pub mod graticule;
pub mod mesh;
pub mod planet;
pub mod terrain;

pub use camera::OrbitCamera;
pub use earth_camera::EarthCamera;
pub use mesh::{PlanetMesh, Vertex};
pub use terrain::{SelectView, Selection, TileCache, TileInstance, TileUpload, TileVertex};
pub use planet::{GradientImage, HeightImage, LayerImage, MapView, PlanetRenderer, RenderStyle, ViewMode};
