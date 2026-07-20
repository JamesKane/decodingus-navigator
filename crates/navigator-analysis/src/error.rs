//! Error type for the analysis layer (plan §6: one `thiserror` enum per layer).

use std::path::PathBuf;

#[derive(Debug, thiserror::Error)]
pub enum AnalysisError {
    #[error("io error on {path}: {source}")]
    Io {
        path: PathBuf,
        #[source]
        source: std::io::Error,
    },

    #[error("{0}")]
    Message(String),

    /// The walk stopped because cancellation was requested (see [`crate::cancel`]).
    ///
    /// A distinct variant, not a `Message`, because callers must be able to tell a user-requested
    /// stop from a real failure: a cancelled walk holds a *partial* result, so its caller has to
    /// skip persisting it, and the UI has to report "cancelled" rather than an error.
    #[error("cancelled")]
    Cancelled,
}

impl AnalysisError {
    pub(crate) fn io(path: impl Into<PathBuf>, source: std::io::Error) -> Self {
        AnalysisError::Io {
            path: path.into(),
            source,
        }
    }
}

/// The text a panic carried, if any — `panic!("…")` payloads are always a `&'static str` or a
/// `String`. Used to surface *what* actually went wrong instead of guessing at a cause.
pub fn panic_text(payload: &(dyn std::any::Any + Send)) -> Option<&str> {
    payload
        .downcast_ref::<&'static str>()
        .copied()
        .or_else(|| payload.downcast_ref::<String>().map(String::as_str))
}

/// Run a BAM/CRAM walk, converting a **panic** into a clean [`AnalysisError`] so one undecodable
/// file fails gracefully instead of unwinding into a cryptic `JoinError`/aborting a worker. The
/// motivating cases are noodles' `todo!()`/`expect()` on inputs it doesn't handle (an unimplemented
/// CRAM data series, or a decode that needs reference bases it wasn't given): without this, such a
/// file panics deep inside the decoder. `what` labels the operation/file for the surfaced message.
///
/// The panic's own text is included rather than a guessed explanation — this is a last-resort net,
/// so it does not know which limitation it caught. Callers that *do* know (see
/// [`crate::index::ensure_index`]) should diagnose the specific case themselves and say what to do
/// about it; anything reaching here is genuinely unclassified.
///
/// `AssertUnwindSafe` is sound here: on a caught panic we discard `f`'s partial state entirely and
/// return an error — no possibly-inconsistent value crosses the boundary. The default panic hook
/// still prints the original message to stderr (useful diagnostics); only the control flow changes.
pub fn guard_walk<T>(what: &str, f: impl FnOnce() -> Result<T, AnalysisError>) -> Result<T, AnalysisError> {
    std::panic::catch_unwind(std::panic::AssertUnwindSafe(f)).unwrap_or_else(|payload| {
        let detail = panic_text(&*payload).unwrap_or("no further detail");
        Err(AnalysisError::Message(format!(
            "{what}: could not decode the alignment — the reader hit a case it does not handle \
             ({detail})"
        )))
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn guard_walk_converts_panic_to_error_and_passes_ok_through() {
        // The default panic hook still prints to stderr; silence it for this test's deliberate panic.
        let prev = std::panic::take_hook();
        std::panic::set_hook(Box::new(|_| {}));
        let caught = guard_walk("decode", || -> Result<u32, AnalysisError> { unimplemented!() });
        std::panic::set_hook(prev);
        assert!(matches!(caught, Err(AnalysisError::Message(_))), "panic → clean Err");

        // A normal Ok/Err result passes straight through (no panic).
        let ok = guard_walk("decode", || Ok::<_, AnalysisError>(7));
        assert!(matches!(ok, Ok(7)));
    }
}
