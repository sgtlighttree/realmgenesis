//! Feeds the renderer's tile cache: generates the tiles selection asks for.
//!
//! - **Native:** a pool of worker threads runs `TileGenerator::generate` (a pure function)
//!   in the background; finished tiles are uploaded on the main thread. Frames never wait.
//! - **Web:** no threads (no SharedArrayBuffer in the preview frame), so tiles are generated
//!   in-frame until a time budget is spent (always at least one per frame).
//!
//! An `epoch` counter tags work. Moving it on (a new world, or the stop button) cancels
//! everything queued: workers skip jobs from an old epoch, and late results are dropped.

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
    /// While true, no new tiles are generated ("stop background rendering").
    paused: bool,
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
            paused: false,
            #[cfg(not(target_arch = "wasm32"))]
            pool: native::Pool::new(),
            last_generated: 0,
        }
    }

    pub fn in_flight(&self) -> usize {
        self.in_flight.len()
    }

    pub fn paused(&self) -> bool {
        self.paused
    }

    /// Stop or resume background tile generation. Stopping cancels the queue as well, so
    /// only tiles already being built finish (one tile each, well under a second).
    pub fn set_paused(&mut self, paused: bool) {
        if paused && !self.paused {
            self.new_epoch();
        }
        self.paused = paused;
    }

    /// Start a new epoch: queued jobs and late results from the old one are dropped.
    fn new_epoch(&mut self) {
        self.epoch += 1;
        self.in_flight.clear();
        #[cfg(not(target_arch = "wasm32"))]
        self.pool.set_epoch(self.epoch);
    }

    /// Switch to a new world. Clears the cache and generates the six root tiles right
    /// away, so the globe is never empty.
    /// Pausing survives this: the roots are made, nothing else until resumed.
    pub fn set_world(&mut self, generator: TileGenerator, cache: &mut TileCache, device: &wgpu::Device, now_s: f64) {
        self.new_epoch();
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
            for (epoch, id, upload) in self.pool.finished() {
                // Only current-epoch results count: after a cancel, the same tile may be
                // in flight again under the new epoch.
                if epoch == self.epoch {
                    self.in_flight.remove(&id);
                    if let Some(upload) = upload {
                        cache.insert(device, upload, now_s);
                        self.last_generated += 1;
                    }
                }
            }
            if self.paused {
                return;
            }
            // Workers skip queued tiles the view no longer wants (e.g. after a drag).
            self.pool.set_wanted(wanted);
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
        if !self.paused {
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
    use std::collections::HashSet;
    use std::sync::atomic::{AtomicU64, Ordering};
    use std::sync::mpsc::{channel, Receiver, Sender};
    use std::sync::{Arc, Mutex};

    type Job = (u64, TileId, TileGenerator);

    /// (epoch, tile, the tile or `None` if it was skipped as no longer wanted).
    type Done = (u64, TileId, Option<TileUpload>);

    /// A small fixed thread pool. Workers pull jobs from a shared queue and send finished
    /// tiles back; dropping the pool closes the queue and the workers exit.
    pub struct Pool {
        jobs: Sender<Job>,
        results: Receiver<Done>,
        pub workers: usize,
        /// The feeder's current epoch. Workers skip queued jobs from older ones.
        epoch: Arc<AtomicU64>,
        /// The tiles the view wants this frame. Queued jobs outside it are skipped.
        wanted: Arc<Mutex<HashSet<TileId>>>,
    }

    impl Pool {
        pub fn new() -> Self {
            // Half the cores, at most 4: tiles are background work, and with every core
            // busy the whole machine locked up (Matt, M1, s02). On an 8-core M1 that's 4.
            let workers = (std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4) / 2).clamp(1, 4);
            let (jobs, job_rx) = channel::<Job>();
            let (result_tx, results) = channel();
            let job_rx = Arc::new(Mutex::new(job_rx));
            let epoch = Arc::new(AtomicU64::new(0));
            let wanted = Arc::new(Mutex::new(HashSet::new()));
            for i in 0..workers {
                let rx = Arc::clone(&job_rx);
                let tx = result_tx.clone();
                let current = Arc::clone(&epoch);
                let wanted_now = Arc::clone(&wanted);
                std::thread::Builder::new()
                    .name(format!("tile-worker-{i}"))
                    .spawn(move || {
                        crate::qos::lower_priority();
                        loop {
                            // Hold the lock only while taking a job, not while working.
                            let job = rx.lock().ok().and_then(|r| r.recv().ok());
                            let Some((epoch, id, generator)) = job else { break };
                            // Cancelled while queued (stop button or a new world): skip it.
                            // A stale read costs at most one extra tile, so Relaxed is enough.
                            if epoch != current.load(Ordering::Relaxed) {
                                continue;
                            }
                            // Panned or zoomed away since it was queued: skip, but report
                            // back so the feeder can request it again if it's wanted later.
                            let still_wanted = wanted_now.lock().map_or(true, |w| w.contains(&id));
                            let tile = still_wanted.then(|| to_upload(generator.generate(id)));
                            if tx.send((epoch, id, tile)).is_err() {
                                break;
                            }
                        }
                    })
                    .expect("failed to spawn tile worker");
            }
            Self { jobs, results, workers, epoch, wanted }
        }

        pub fn set_wanted(&self, ids: &[TileId]) {
            if let Ok(mut w) = self.wanted.lock() {
                w.clear();
                w.extend(ids.iter().copied());
            }
        }

        pub fn set_epoch(&self, epoch: u64) {
            self.epoch.store(epoch, Ordering::Relaxed);
        }

        pub fn submit(&self, epoch: u64, id: TileId, generator: TileGenerator) {
            let _ = self.jobs.send((epoch, id, generator));
        }

        pub fn finished(&self) -> Vec<Done> {
            self.results.try_iter().collect()
        }
    }
}
