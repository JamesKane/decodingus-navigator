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
}

impl AnalysisError {
    pub(crate) fn io(path: impl Into<PathBuf>, source: std::io::Error) -> Self {
        AnalysisError::Io {
            path: path.into(),
            source,
        }
    }
}

/// Run a BAM/CRAM walk, converting a **panic** into a clean [`AnalysisError`] so one undecodable
/// file fails gracefully instead of unwinding into a cryptic `JoinError`/aborting a worker. The
/// motivating case is noodles' `todo!()` on CRAM features it hasn't implemented (e.g. a Huffman-
/// coded byte data series, as some FTDNA Big Y CRAMs use): without this, analyzing such a file
/// panics deep inside the decoder. `what` labels the operation/file for the surfaced message.
///
/// `AssertUnwindSafe` is sound here: on a caught panic we discard `f`'s partial state entirely and
/// return an error — no possibly-inconsistent value crosses the boundary. The default panic hook
/// still prints the original message to stderr (useful diagnostics); only the control flow changes.
pub fn guard_walk<T>(what: &str, f: impl FnOnce() -> Result<T, AnalysisError>) -> Result<T, AnalysisError> {
    std::panic::catch_unwind(std::panic::AssertUnwindSafe(f)).unwrap_or_else(|_| {
        Err(AnalysisError::Message(format!(
            "{what}: could not decode the alignment — it may use a BAM/CRAM feature not yet \
             supported by the reader (e.g. a Huffman-coded CRAM data series)"
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
