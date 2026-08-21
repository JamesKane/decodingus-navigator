//! The error type of the analysis layer. See plan §6: one `thiserror` enum in each layer.

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

    /// The walk stopped because somebody asked for a cancel. See [`crate::cancel`].
    ///
    /// This is a variant of its own, and not a `Message`. A caller must be able to separate a stop
    /// that a user asked for from a real failure. A walk that somebody cancelled holds a *partial*
    /// result, so its caller must not put that result into the store. And the UI must report
    /// "cancelled", and not an error.
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

/// The text that a panic carried, when it carried one. The payload of a `panic!("…")` is always a
/// `&'static str` or a `String`. Use it to show *what* went wrong, and do not guess at a cause.
pub fn panic_text(payload: &(dyn std::any::Any + Send)) -> Option<&str> {
    payload
        .downcast_ref::<&'static str>()
        .copied()
        .or_else(|| payload.downcast_ref::<String>().map(String::as_str))
}

/// Run a walk over a BAM or CRAM, and turn a **panic** into a clean [`AnalysisError`]. One file
/// that the code can not decode then fails cleanly. It does not unwind into a `JoinError` that
/// says nothing, and it does not abort a worker.
///
/// The cases that led to this are the `todo!()` and `expect()` calls of noodles, on an input that
/// it does not handle. A CRAM data series that nobody implemented is one. A decode that needs
/// reference bases which nobody gave it is another. Without this net, such a file panics deep
/// inside the decoder. `what` names the operation and the file, for the message that goes out.
///
/// The message holds the own text of the panic, and not an explanation that this code guessed.
/// This is a net of last resort, so it does not know which limit it caught. A caller that *does*
/// know must diagnose its own case and say what to do about it. See
/// [`crate::index::ensure_index`]. This code can put no class on anything that reaches here.
///
/// `AssertUnwindSafe` is sound here. On a panic that this code catches, it throws away the whole
/// partial state of `f`, and it returns an error. No value that could be inconsistent crosses the
/// boundary. The default panic hook still prints the original message to stderr, which is useful,
/// and only the control flow changes.
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
