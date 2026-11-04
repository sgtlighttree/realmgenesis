//! Thread priority for background work. Tile workers and world generation run at macOS
//! "utility" QoS, so the UI and other apps are scheduled first (s02: with every core busy at
//! normal priority, the app locked up the M1). Threads a utility thread spawns are expected
//! to inherit the class (not yet verified for the scoped pools inside generation).

/// Lower the calling thread to utility QoS on macOS; does nothing elsewhere.
pub fn lower_priority() {
    #[cfg(target_os = "macos")]
    // SAFETY: only changes the calling thread's scheduling class; plain-value arguments.
    // A failure just leaves the default QoS, so the result is ignored.
    unsafe {
        libc::pthread_set_qos_class_self_np(libc::qos_class_t::QOS_CLASS_UTILITY, 0);
    }
}
