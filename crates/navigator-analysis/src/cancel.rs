//! Cooperative cancellation for the long walks.
//!
//! A pass over the whole genome takes minutes. For all of those minutes, the Cancel button of the
//! UI once did nothing that a user could see. The flag that it set lived in `navigator-ui`, and
//! the code read that flag only *between* the steps of the pipeline. The step itself ran inside a
//! `spawn_blocking` closure, and tokio can not interrupt one of those. Once a walk starts, the one
//! thing that can stop it is the walk itself. So a walker has to ask.
//!
//! [`CancelToken`] is that question. The rule for its use is about *where* you ask. Ask often
//! enough that a click feels immediate, and rarely enough that the check does not show in a
//! profile.
//!
//! Every check in this crate sits on a path that already does real work at each record, or at each
//! contig. An atomic load is noise next to that. A check inside the innermost loop over the bases
//! would not be.
//!
//! A walk that somebody cancelled returns [`AnalysisError::Cancelled`], and not a partial result.
//! Once the store holds a partial coverage, nothing can separate it from a coverage that is truly
//! low. To cache a walk that stopped half way, as if it were complete, is a far worse failure than
//! no cancellation at all. So a cancellation is an error, and a caller skips its store step on
//! it.

use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::Arc;

use crate::error::AnalysisError;

/// A shared flag that says "stop what you are doing". A clone of it, into a worker thread, costs
/// almost nothing.
///
/// [`CancelToken::none`] gives a token that nobody can cancel. It exists so that a caller with
/// nothing to cancel pays nothing, and reads naturally. A test, a one-shot CLI command and the
/// wrappers that report no progress are all such callers. Without it, every signature would carry
/// an `Option`.
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
    /// There is no way back to the state before the cancel, and that is deliberate. A token covers
    /// exactly one run. To use one token across two runs is what let a stale reset write over a
    /// cancel that had not yet arrived.
    pub fn cancel(&self) {
        if let Some(flag) = &self.0 {
            flag.store(true, Ordering::Relaxed);
        }
    }

    /// True when somebody has asked for a cancel.
    ///
    /// `Relaxed` is enough. This flag guards no other memory. To see the store one iteration of
    /// the loop late costs one more iteration of work, and nothing else.
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

    /// `Default` must give the token that does nothing. Take a struct that gains a `CancelToken`
    /// field through `..Default::default()`. That struct must not start in a state where somebody
    /// can cancel it, and nobody ever does, and where it differs from `none()`. Nobody would see
    /// that.
    #[test]
    fn default_is_the_inert_token() {
        let t = CancelToken::default();
        t.cancel();
        assert!(!t.is_cancelled());
    }

    /// The property that the whole feature stands on. A walk that already runs sees a token that
    /// *another thread* cancelled. The old design could not express this case at all: the flag
    /// lived in the UI, and the walk had no way to ask.
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
        assert!(
            handle.join().unwrap().is_err(),
            "the loop must stop, not run to completion"
        );
    }

    /// A cancellation must go out as itself, and never as a general failure. The UI branches on
    /// that, so it does not tell the user that their own click was an error.
    #[test]
    fn cancellation_is_its_own_error_kind() {
        let t = CancelToken::new();
        t.cancel();
        let e = t.check().unwrap_err();
        assert!(matches!(e, AnalysisError::Cancelled));
        assert_eq!(e.to_string(), "cancelled");
    }
}
