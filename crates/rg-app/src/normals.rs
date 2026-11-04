//! Bakes the surface-gradient texture (per-pixel lighting, R2) off the UI thread.
//!
//! A world change starts a bake. Until it lands, the globe is lit flat (the renderer's
//! gradient is reset to zero), so the old world's slopes never light the new one. Natively
//! the bake runs on a utility-QoS thread; the web build bakes synchronously.

use rg_render::GradientImage;
use rg_terrain::TerrainSampler;
use std::sync::mpsc::{Receiver, TryRecvError};

pub struct NormalsJob {
    rx: Receiver<(GradientImage, f64)>,
}

pub enum NormalsPoll {
    Working,
    /// The image and the bake time in seconds.
    Done(GradientImage, f64),
    Failed,
}

/// Bake the gradient for `sampler` (the displayed, macro-only surface).
pub fn bake(sampler: &TerrainSampler) -> GradientImage {
    let res = rg_terrain::gradient::gradient_res(sampler.macro_h.res);
    let field = rg_terrain::gradient::bake_gradient(sampler, res);
    GradientImage { res, data: field.data }
}

impl NormalsJob {
    #[cfg(not(target_arch = "wasm32"))]
    pub fn spawn(sampler: TerrainSampler) -> Result<Self, String> {
        let (tx, rx) = std::sync::mpsc::channel();
        std::thread::Builder::new()
            .name("normals-bake".into())
            .spawn(move || {
                crate::qos::lower_priority();
                let t0 = web_time::Instant::now();
                let img = bake(&sampler);
                let _ = tx.send((img, t0.elapsed().as_secs_f64()));
            })
            .map_err(|e| e.to_string())?;
        Ok(Self { rx })
    }

    pub fn poll(&self) -> NormalsPoll {
        match self.rx.try_recv() {
            Ok((img, s)) => NormalsPoll::Done(img, s),
            Err(TryRecvError::Empty) => NormalsPoll::Working,
            Err(TryRecvError::Disconnected) => NormalsPoll::Failed,
        }
    }
}
