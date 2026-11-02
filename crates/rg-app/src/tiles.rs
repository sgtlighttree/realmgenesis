//! Feeds the renderer's tile cache: generates the tiles selection asks for.
//!
//! - **Native:** a pool of worker threads runs `TileGenerator::generate` (a pure function)
//!   in the background; finished tiles are uploaded on the main thread. Frames never wait.
//! - **Web:** no threads (no SharedArrayBuffer in the preview frame), so tiles are generated
//!   in-frame until a time budget is spent (always at least one per frame).
//!
//! An `epoch` counter tags work, so results for a previous world are dropped.

use rg_geo::{CubeFace, TileId};
use rg_render::{TileCache, TileUpload, TileVertex};
use rg_terrain::{TileData, TileGenerator};
use std::collections::HashSet;

/// Convert generator output into the renderer's vertex layout.
pub fn to_upload(t: TileData) -> TileUpload {
    let vertices = (0..t.positions.len())
        .map(|i| TileVertex {
            pos: t.positions[i],
            height_m: t.heights[i],
            normal: t.normals[i],
            parent_height_m: t.parent_heights[i],
            skirt: t.skirt[i],
        })
        .collect();
    TileUpload {
        id: t.id,
        center: t.center,
        vertices,
        patch_radius_m: t.patch_radius_m,
        max_height_m: t.max_height_m,
        child_error_m: t.child_error_m,
        parent_error_m: t.parent_error_m,
    }
}

pub struct TileFeeder {
    generator: Option<TileGenerator>,
    epoch: u64,
    in_flight: HashSet<TileId>,
    #[cfg(not(target_arch = "wasm32"))]
    pool: native::Pool,
    /// Tiles generated in the last update (for the HUD).
    pub last_generated: usize,
}

impl TileFeeder {
    pub fn new() -> Self {
        Self {
            generator: None,
            epoch: 0,
            in_flight: HashSet::new(),
            #[cfg(not(target_arch = "wasm32"))]
            pool: native::Pool::new(),
            last_generated: 0,
        }
    }

    pub fn in_flight(&self) -> usize {
        self.in_flight.len()
    }

    /// Switch to a new world. Clears the cache and generates the six root tiles right
    /// away, so the globe is never empty.
    pub fn set_world(&mut self, generator: TileGenerator, cache: &mut TileCache, device: &wgpu::Device, now_s: f64) {
        self.epoch += 1;
        self.in_flight.clear();
        cache.clear();
        for face in CubeFace::ALL {
            cache.insert(device, to_upload(generator.generate(TileId::root(face))), now_s);
        }
        self.generator = Some(generator);
    }

    /// Request `wanted` (most important first) and upload whatever has finished.
    pub fn update(&mut self, wanted: &[TileId], cache: &mut TileCache, device: &wgpu::Device, now_s: f64) {
        self.last_generated = 0;
        let Some(generator) = self.generator.as_ref() else { return };

        #[cfg(not(target_arch = "wasm32"))]
        {
            // Collect finished work first, then top the queue up.
            for (epoch, upload) in self.pool.finished() {
                self.in_flight.remove(&upload.id);
                if epoch == self.epoch {
                    cache.insert(device, upload, now_s);
                    self.last_generated += 1;
                }
            }
            let cap = self.pool.workers * 3;
            for id in wanted {
                if self.in_flight.len() >= cap {
                    break;
                }
                if !self.in_flight.contains(id) && !cache.contains(id) {
                    self.in_flight.insert(*id);
                    self.pool.submit(self.epoch, *id, generator.clone());
                }
            }
        }

        #[cfg(target_arch = "wasm32")]
        {
            const BUDGET_S: f64 = 0.006;
            let t0 = web_time::Instant::now();
            for id in wanted {
                if cache.contains(id) {
                    continue;
                }
                cache.insert(device, to_upload(generator.generate(*id)), now_s);
                self.last_generated += 1;
                if t0.elapsed().as_secs_f64() > BUDGET_S {
                    break;
                }
            }
        }
    }
}

#[cfg(not(target_arch = "wasm32"))]
mod native {
    use super::to_upload;
    use rg_geo::TileId;
    use rg_render::TileUpload;
    use rg_terrain::TileGenerator;
    use std::sync::mpsc::{channel, Receiver, Sender};
    use std::sync::{Arc, Mutex};

    type Job = (u64, TileId, TileGenerator);

    /// A small fixed thread pool. Workers pull jobs from a shared queue and send finished
    /// tiles back; dropping the pool closes the queue and the workers exit.
    pub struct Pool {
        jobs: Sender<Job>,
        results: Receiver<(u64, TileUpload)>,
        pub workers: usize,
    }

    impl Pool {
        pub fn new() -> Self {
            let workers = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4).saturating_sub(1).clamp(1, 8);
            let (jobs, job_rx) = channel::<Job>();
            let (result_tx, results) = channel();
            let job_rx = Arc::new(Mutex::new(job_rx));
            for i in 0..workers {
                let rx = Arc::clone(&job_rx);
                let tx = result_tx.clone();
                std::thread::Builder::new()
                    .name(format!("tile-worker-{i}"))
                    .spawn(move || {
                        loop {
                            // Hold the lock only while taking a job, not while working.
                            let job = rx.lock().ok().and_then(|r| r.recv().ok());
                            let Some((epoch, id, generator)) = job else { break };
                            if tx.send((epoch, to_upload(generator.generate(id)))).is_err() {
                                break;
                            }
                        }
                    })
                    .expect("failed to spawn tile worker");
            }
            Self { jobs, results, workers }
        }

        pub fn submit(&self, epoch: u64, id: TileId, generator: TileGenerator) {
            let _ = self.jobs.send((epoch, id, generator));
        }

        pub fn finished(&self) -> Vec<(u64, TileUpload)> {
            self.results.try_iter().collect()
        }
    }
}
