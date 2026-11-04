//! The status-bar console (s02; V1's "System Console" analogue): a short in-memory log of
//! what the app is doing, so long work never looks like a freeze.

use std::collections::VecDeque;

pub struct Console {
    lines: VecDeque<String>,
}

impl Console {
    const MAX_LINES: usize = 100;

    pub fn new() -> Self {
        Self { lines: VecDeque::new() }
    }

    /// Add a line stamped with the time since launch (`mm:ss.s`). Native builds also print
    /// it to stderr, so a terminal launch keeps a full log.
    pub fn log(&mut self, now_s: f64, msg: impl Into<String>) {
        let line = format!("{:02}:{:04.1}  {}", (now_s / 60.0) as u32, now_s % 60.0, msg.into());
        #[cfg(not(target_arch = "wasm32"))]
        eprintln!("{line}");
        self.lines.push_back(line);
        while self.lines.len() > Self::MAX_LINES {
            self.lines.pop_front();
        }
    }

    /// Oldest first.
    pub fn lines(&self) -> impl Iterator<Item = &str> {
        self.lines.iter().map(String::as_str)
    }
}
