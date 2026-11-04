//! World generation on a background thread, so the UI keeps drawing and shows live progress
//! (s02). Native only: the web build has no threads and keeps the synchronous path in
//! `app.rs`. (This module compiles everywhere; only the native code calls `spawn`.)

use rg_pipeline::{Progress, World, WorldParams};
use std::sync::mpsc::{channel, Receiver, TryRecvError};

enum Msg {
    Progress(Progress),
    Done(Box<World>),
}

pub enum Poll {
    Working,
    Done(Box<World>),
    /// The thread ended without a world (it panicked).
    Failed,
}

pub struct GenJob {
    rx: Receiver<Msg>,
    /// App clock when the job started (seconds).
    pub started_s: f64,
    /// The latest progress report; `None` until the thread reports in.
    pub latest: Option<Progress>,
}

impl GenJob {
    pub fn spawn(params: WorldParams, now_s: f64) -> Result<Self, String> {
        let (tx, rx) = channel();
        std::thread::Builder::new()
            .name("world-gen".into())
            .spawn(move || {
                crate::qos::lower_priority();
                let report = |p| {
                    let _ = tx.send(Msg::Progress(p));
                };
                let world = rg_pipeline::generate_with_progress(&params, &report);
                let _ = tx.send(Msg::Done(Box::new(world)));
            })
            .map_err(|e| e.to_string())?;
        Ok(Self { rx, started_s: now_s, latest: None })
    }

    /// Take every message waiting. Returns the stage changes seen (for the console) and
    /// the job's state.
    pub fn poll(&mut self) -> (Vec<Progress>, Poll) {
        let mut stage_changes = Vec::new();
        loop {
            match self.rx.try_recv() {
                Ok(Msg::Progress(p)) => {
                    if self.latest.map(|l| l.stage) != Some(p.stage) {
                        stage_changes.push(p);
                    }
                    self.latest = Some(p);
                }
                Ok(Msg::Done(world)) => return (stage_changes, Poll::Done(world)),
                Err(TryRecvError::Empty) => return (stage_changes, Poll::Working),
                Err(TryRecvError::Disconnected) => return (stage_changes, Poll::Failed),
            }
        }
    }
}
