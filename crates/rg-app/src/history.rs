//! Edit history (s06, Matt: "edits need to be reversible and non-destructive"): undo/redo of
//! the coast edits (strokes, imports, placement, the land test), one entry per gesture.
//!
//! An entry is a snapshot of the whole `Conditioning`. That's cheap: the rasters sit behind
//! `Arc`s, so a snapshot shares them. What can grow is the set of *distinct* rasters (each
//! crop or channel change decodes a new one, ~250 MB for an 11200×5600 PNG), so the history
//! is capped by those bytes as well as by count.

use rg_pipeline::mask::{Conditioning, RasterMask};
use std::sync::Arc;

const MAX_ENTRIES: usize = 80;
const MAX_RASTER_BYTES: usize = 1 << 30; // 1 GiB of distinct rasters held only by history

/// One undo step: the edits plus the UI's name for the imported file.
#[derive(Clone)]
pub struct Snapshot {
    pub cond: Conditioning,
    pub file_name: String,
}

#[derive(Default)]
pub struct EditHistory {
    undo: Vec<Snapshot>,
    redo: Vec<Snapshot>,
    /// A change is in progress under a held pointer (a slider or patch drag): later changes
    /// in the same gesture don't add entries.
    gesture_open: bool,
    /// Set when history itself changed the edits this frame (undo, redo, a checkpoint), so
    /// the end-of-frame tracker doesn't record that as a new edit.
    skip_track: bool,
    /// The edits the shown world was generated from, to tell when it's out of date.
    pub generated: Option<Conditioning>,
}

/// Same edits, compared cheaply: rasters by pointer, never by content (a full compare of a
/// big mask every frame would cost milliseconds).
pub fn same_edits(a: &Conditioning, b: &Conditioning) -> bool {
    a.strict == b.strict
        && a.outside_procedural == b.outside_procedural
        && a.guard_km == b.guard_km
        && a.islands == b.islands
        && a.islands_guard_km == b.islands_guard_km
        && a.strokes == b.strokes
        && a.masks.len() == b.masks.len()
        && a.masks.iter().zip(&b.masks).all(|(a, b)| same_mask(a, b))
}

fn same_mask(a: &Arc<RasterMask>, b: &Arc<RasterMask>) -> bool {
    Arc::ptr_eq(a, b)
        || (Arc::ptr_eq(&a.data, &b.data)
            && a.place == b.place
            && a.adjust == b.adjust
            && a.decode == b.decode
            && a.enabled == b.enabled
            && a.name == b.name)
}

impl EditHistory {
    /// Record `cur` as an undo point now (before an import, or entering "Move patch").
    /// Returns the entry's index, for [`Self::cancel_to`].
    pub fn checkpoint(&mut self, cur: &Conditioning, file_name: &str) -> usize {
        self.undo.push(Snapshot { cond: cur.clone(), file_name: file_name.to_string() });
        self.redo.clear();
        self.skip_track = true;
        self.gesture_open = false;
        self.trim();
        self.undo.len() - 1
    }

    /// End-of-frame tracking: if the edits changed this frame, the state before the change
    /// becomes an undo entry (once per gesture).
    pub fn track(&mut self, before: &Snapshot, after: &Conditioning, pointer_down: bool) {
        if std::mem::take(&mut self.skip_track) {
            self.gesture_open = false;
            return;
        }
        if !same_edits(&before.cond, after) {
            if !self.gesture_open {
                self.undo.push(before.clone());
                self.redo.clear();
                self.trim();
            }
            self.gesture_open = pointer_down;
        } else if !pointer_down {
            self.gesture_open = false;
        }
    }

    pub fn can_undo(&self) -> bool {
        !self.undo.is_empty()
    }

    pub fn can_redo(&self) -> bool {
        !self.redo.is_empty()
    }

    pub fn undo(&mut self, cur: &mut Snapshot) -> bool {
        let Some(prev) = self.undo.pop() else { return false };
        self.redo.push(std::mem::replace(cur, prev));
        self.skip_track = true;
        true
    }

    pub fn redo(&mut self, cur: &mut Snapshot) -> bool {
        let Some(next) = self.redo.pop() else { return false };
        self.undo.push(std::mem::replace(cur, next));
        self.skip_track = true;
        true
    }

    /// Esc in "Move patch": go back to the checkpoint at `mark`, dropping the steps since
    /// (a cancel, not an undo, so they don't go to redo).
    pub fn cancel_to(&mut self, mark: usize, cur: &mut Snapshot) {
        if mark >= self.undo.len() {
            return;
        }
        self.undo.truncate(mark + 1);
        *cur = self.undo.pop().unwrap();
        self.skip_track = true;
    }

    /// Drop the oldest entries past the count cap, or while the rasters only history holds
    /// exceed the byte cap.
    fn trim(&mut self) {
        while self.undo.len() > MAX_ENTRIES {
            self.undo.remove(0);
        }
        while self.undo.len() > 1 && self.raster_bytes() > MAX_RASTER_BYTES {
            self.undo.remove(0);
        }
    }

    fn raster_bytes(&self) -> usize {
        let mut seen: Vec<*const Vec<f32>> = Vec::new();
        let mut bytes = 0;
        for m in self.undo.iter().chain(&self.redo).flat_map(|s| s.cond.masks.iter()) {
            let p = Arc::as_ptr(&m.data);
            if !seen.contains(&p) {
                seen.push(p);
                bytes += m.data.len() * 4;
            }
        }
        bytes
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rg_pipeline::mask::{Stroke, StrokeOp};

    fn snap(c: &Conditioning) -> Snapshot {
        Snapshot { cond: c.clone(), file_name: String::new() }
    }

    fn stroke() -> Stroke {
        Stroke::new(vec![glam::DVec3::X], 100.0, StrokeOp::Sea)
    }

    #[test]
    fn one_entry_per_gesture_and_undo_redo() {
        let mut h = EditHistory::default();
        let mut c = Conditioning::default();
        // A click (pointer up) that adds a stroke: one entry.
        let before = snap(&c);
        c.strokes.push(stroke());
        h.track(&before, &c, false);
        // A drag over three frames that changes the guard each frame: one entry.
        for k in 0..3 {
            let before = snap(&c);
            c.guard_km = 900.0 + k as f64;
            h.track(&before, &c, true);
        }
        h.track(&snap(&c), &c, false); // release
        assert_eq!(h.undo.len(), 2);
        let mut cur = snap(&c);
        assert!(h.undo(&mut cur));
        assert_eq!(cur.cond.guard_km, 800.0);
        assert_eq!(cur.cond.strokes.len(), 1);
        h.track(&snap(&cur.cond), &cur.cond, false); // the undo frame isn't a new edit
        assert!(h.undo(&mut cur));
        assert!(cur.cond.strokes.is_empty());
        assert!(h.redo(&mut cur) && h.redo(&mut cur));
        assert_eq!(cur.cond.guard_km, 902.0);
        assert!(!h.can_redo());
    }

    #[test]
    fn cancel_returns_to_the_checkpoint() {
        let mut h = EditHistory::default();
        let c0 = Conditioning::default();
        let mark = h.checkpoint(&c0, "");
        let mut cur = snap(&c0);
        cur.cond.guard_km = 1234.0;
        h.track(&snap(&c0), &cur.cond, false);
        h.cancel_to(mark, &mut cur);
        assert_eq!(cur.cond.guard_km, 800.0);
        assert!(!h.can_undo());
    }
}
