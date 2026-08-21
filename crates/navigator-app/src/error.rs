//! Application-layer error: store failures plus artifact (de)serialization.

#[derive(Debug, thiserror::Error)]
pub enum AppError {
    #[error(transparent)]
    Store(#[from] navigator_store::StoreError),

    #[error(transparent)]
    Analysis(#[from] navigator_analysis::AnalysisError),

    /// A fault in read mapping, which is realignment stage B.
    ///
    /// This is a separate variant, not part of `Analysis`. `navigator-align` is a separate crate
    /// with its own error type. Also, a mapping fault has a different cause from an analysis
    /// fault.
    #[error("{0}")]
    Align(#[from] navigator_align::AlignError),

    #[error("serialization error: {0}")]
    Serde(#[from] serde_json::Error),

    /// An analysis task on a blocked thread did not join. The task had a panic, or the user
    /// stopped it.
    #[error("analysis task failed: {0}")]
    Join(String),

    #[error("alignment {0} has no BAM/reference path recorded")]
    MissingPaths(i64),

    /// The file of the alignment is no longer on disk. The cause can be a newer vendor download, a
    /// deleted import, or a volume that the user removed.
    ///
    /// This variant is different from [`AppError::MissingPaths`]. That variant means that the
    /// alignment never had a path. Here the alignment has a path, but the path no longer points to
    /// a file.
    ///
    /// This fault has its own variant for two reasons. First, it is the one read fault that a long
    /// life workspace can expect, and no user made a mistake. So a sweep across many alignments
    /// skips the alignment with [`AppError::is_missing_alignment_file`] and does not count a
    /// failure.
    ///
    /// Second, the code raises this error before the large setup that a walk needs. Before this
    /// variant, the reader gave an unclear io error from deep in its code, and a user read that
    /// error as an absent haplotree.
    #[error("alignment {id} file is no longer at {path}")]
    AlignmentFileMissing { id: i64, path: String },

    /// The ancestry reference panel file is absent. Make the file with `navigator-panelbuild` and
    /// install it. As an alternative, set `$NAVIGATOR_ANCESTRY_PANEL`.
    #[error("ancestry panel not found at {0} — build it with navigator-panelbuild")]
    AncestryPanelMissing(std::path::PathBuf),

    /// The bundled panel is for a different reference build than the alignment.
    #[error("ancestry panel build {panel} does not match alignment build {alignment}")]
    AncestryPanelBuildMismatch { panel: String, alignment: String },

    /// Too few sites genotyped for a reliable ancestry estimate.
    #[error("insufficient data for ancestry: {genotyped} SNPs genotyped, {required} required")]
    InsufficientAncestryData { genotyped: usize, required: usize },

    #[error("not signed in — log in to a PDS account first")]
    NotAuthenticated,

    /// An AppView API call failed (e.g. federated IBD). 403 → the device key is not
    /// registered/verified yet; 422 → clock skew; otherwise the server's reason.
    #[error("appview error: {0}")]
    AppView(String),

    #[error("could not read file: {0}")]
    Io(#[from] std::io::Error),

    #[error("import error: {0}")]
    Import(String),

    #[error(transparent)]
    Sync(#[from] navigator_sync::SyncError),

    #[error(transparent)]
    Refgenome(#[from] navigator_refgenome::RefgenomeError),

    /// The import needs one or more reference builds that the cache does not hold. The UI asks the
    /// user, downloads the builds through the gateway, and then tries again. The code wrote nothing
    /// to the database.
    #[error("reference download required: {0:?}")]
    ReferenceNeeded(Vec<crate::BuildNeed>),

    /// The app refuses a change because of the current state. One example is a request to delete a
    /// subject that still has sequence data or a profile.
    #[error("{0}")]
    Conflict(String),

    /// A local-LLM operation failed (server unreachable, bad response, etc.). The message is
    /// plain-language for the Settings "Test connection" UI.
    #[error("{0}")]
    Llm(String),

    /// An installer-update check failed (GitHub Releases unreachable, bad response, etc.). Surfaced
    /// as a plain-language status line; a failed check is non-fatal (the app runs regardless).
    #[error("update check failed: {0}")]
    Update(String),
}

impl AppError {
    /// Shows that the user stopped the job. It does not show a fault.
    ///
    /// A stop moves through the code as an error, because an error unwinds the walk from any point.
    /// But the caller must know the difference between the two. A run that the user stopped holds a
    /// partial result, and the app must not write that result to the store. The UI must also show
    /// the word "cancelled", not an error.
    ///
    /// This method is in this crate, not in the UI. So a layer above `navigator-app` never needs to
    /// look inside this crate.
    pub fn is_cancellation(&self) -> bool {
        matches!(self, AppError::Analysis(navigator_analysis::AnalysisError::Cancelled))
    }

    /// Shows that the only fault is an absent alignment file. A sweep must skip such an alignment
    /// and must not count a failure. See [`AppError::AlignmentFileMissing`].
    pub fn is_missing_alignment_file(&self) -> bool {
        matches!(self, AppError::AlignmentFileMissing { .. })
    }
}

impl From<tokio::task::JoinError> for AppError {
    fn from(e: tokio::task::JoinError) -> Self {
        AppError::Join(e.to_string())
    }
}

impl AppError {
    /// Shows that the user stopped the job. It does not show a fault.
    ///
    /// A long job must know the difference between the two. To report the Cancel click of the user
    /// as an error is incorrect, and it makes the user afraid. This method holds the difference, so
    /// a caller does not compare the text of a message.
    pub fn is_cancelled(&self) -> bool {
        matches!(
            self,
            AppError::Analysis(navigator_analysis::AnalysisError::Cancelled)
                | AppError::Align(navigator_align::AlignError::Cancelled)
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_cancel_is_distinguishable_from_a_failure() {
        assert!(AppError::Analysis(navigator_analysis::AnalysisError::Cancelled).is_cancelled());
        assert!(AppError::Align(navigator_align::AlignError::Cancelled).is_cancelled());
        assert!(!AppError::Import("disk full".into()).is_cancelled());
        assert!(
            !AppError::Analysis(navigator_analysis::AnalysisError::Message("boom".into())).is_cancelled(),
            "a message that happens to be an analysis error is still a failure"
        );
    }
}
