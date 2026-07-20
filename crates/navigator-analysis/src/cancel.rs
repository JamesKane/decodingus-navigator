//! Cooperative cancellation for the long walks.
//!
//! A whole-genome pass takes minutes, and the UI's Cancel button used to do nothing visible for
//! all of them: the flag it set lived in `navigator-ui` and was only read *between* pipeline steps,
//! while the step itself ran inside a `spawn_blocking` closure that tokio cannot interrupt. Once a
//! walk starts, the only thing that can stop it is the walk itself — so the walkers have to ask.
//!
//! [`CancelToken`] is that question, and the rule for using it is about *where* you ask: often
//! enough that a click feels instant, rarely enough that the check does not show up in a profile.
//! Every place one is checked here sits on a path that already does real per-record or per-contig
//! work, so an atomic load is noise by comparison. Checking inside the innermost per-base loop
//! would not be.
//!
//! A cancelled walk returns [`AnalysisError::Cancelled`] rather than a partial result. Partial
//! coverage is indistinguishable from genuinely low coverage once it is persisted, and silently
//! caching a half-finished walk as if it were complete is a far worse failure than not cancelling
//! at all — so cancellation is an error, and callers skip their persistence step on it.

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

use crate::error::AnalysisError;

/// A shared "stop what you're doing" flag, cheap to clone into worker threads.
///
/// [`CancelToken::none`] is a token that can never be cancelled. It exists so callers with nothing
/// to cancel — tests, CLI one-shots, the non-progress convenience wrappers — pay nothing and read
/// naturally, instead of every signature growing an `Option`.
#[derive(Clone, Debug, Default)]
pub struct CancelToken(Option<Arc<AtomicBool>>);

impl CancelToken {
    /// A live token. Clone it to the canceller and to the work.
    pub fn new() -> Self {
        Self(Some(Arc::new(AtomicBool::new(false))))
    }

    /// A token that is never cancelled.
    pub fn none() -> Self {
        Self(None)
    }

    /// Request cancellation. Idempotent, and safe to call from any thread.
    ///
    /// There is deliberately no way back to the un-cancelled state: a token covers exactly one run,
    /// and reusing one across runs is what let a stale reset clobber a pending cancel before.
    pub fn cancel(&self) {
        if let Some(flag) = &self.0 {
            flag.store(true, Ordering::Relaxed);
        }
    }

    /// Whether cancellation has been requested.
    ///
    /// `Relaxed` is sufficient: this guards no other memory, and the only cost of observing the
    /// store one loop iteration late is one more iteration of work.
    pub fn is_cancelled(&self) -> bool {
        self.0.as_ref().is_some_and(|flag| flag.load(Ordering::Relaxed))
    }

    /// `Err(AnalysisError::Cancelled)` when cancelled, for `?` inside a walk loop.
    pub fn check(&self) -> Result<(), AnalysisError> {
        if self.is_cancelled() {
            return Err(AnalysisError::Cancelled);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_none_token_never_cancels_even_when_asked() {
        let t = CancelToken::none();
        t.cancel();
        assert!(!t.is_cancelled());
        assert!(t.check().is_ok());
    }

    #[test]
    fn cancelling_one_clone_cancels_them_all() {
        let t = CancelToken::new();
        let worker = t.clone();
        assert!(!worker.is_cancelled());
        t.cancel();
        assert!(worker.is_cancelled(), "a clone shares the flag");
        assert!(matches!(worker.check(), Err(AnalysisError::Cancelled)));
    }

    /// `Default` has to be the inert token: a struct that gains a `CancelToken` field by
    /// `..Default::default()` must not silently start out cancellable-but-never-cancelled in a way
    /// that differs from `none()`.
    #[test]
    fn default_is_the_inert_token() {
        let t = CancelToken::default();
        t.cancel();
        assert!(!t.is_cancelled());
    }

    /// The property the whole feature rests on: a token cancelled from *another thread* is observed
    /// by a walk already in progress. This is the case the old design could not express at all —
    /// the flag lived in the UI and the walk had no way to ask.
    #[test]
    fn a_walk_in_progress_observes_a_cancel_from_another_thread() {
        let token = CancelToken::new();
        let canceller = token.clone();
        let handle = std::thread::spawn(move || {
            // Stand in for a record loop: poll on the same cadence the walkers use.
            for i in 0..100_000_000u64 {
                if i % 4096 == 0 && token.check().is_err() {
                    return Err(i);
                }
            }
            Ok(())
        });
        canceller.cancel();
        assert!(handle.join().unwrap().is_err(), "the loop must stop, not run to completion");
    }

    /// Cancellation must be reported as itself, never as a generic failure — the UI branches on
    /// this to avoid telling the user their own click was an error.
    #[test]
    fn cancellation_is_its_own_error_kind() {
        let t = CancelToken::new();
        t.cancel();
        let e = t.check().unwrap_err();
        assert!(matches!(e, AnalysisError::Cancelled));
        assert_eq!(e.to_string(), "cancelled");
    }
}
