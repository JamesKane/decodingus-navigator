//! Registering a realigned alignment — stage D of the realignment module.
//!
//! Stage C leaves a sorted, duplicate-marked CRAM on disk. That file is not yet an *alignment* as
//! far as the workspace is concerned; this is where it becomes one, so every existing analysis can
//! run against it without knowing it was produced rather than imported.
//!
//! ## Additive, always
//!
//! The realigned row is inserted under the **same `SequenceRun`** as its source — the same physical
//! library, only mapped differently — and the source is left exactly as it was. That mirrors how
//! the sidecar fast path stayed additive, and it is the property that makes realignment safe to
//! offer: nothing a user already had can be lost by trying it, and a realigned alignment can be
//! deleted and rebuilt from a source that never changed.
//!
//! ## The reference is part of the file
//!
//! A CRAM cannot be read without the reference it was compressed against, so `reference_path` is
//! recorded on the row rather than resolved by convention later. An alignment whose reference has
//! moved is unreadable, and the row is the only place that knows which one it was.

use std::path::{Path, PathBuf};

use navigator_domain::workspace::{Alignment, NewAlignment};
use navigator_store::alignment;

use crate::error::AppError;
use crate::{sha256_file_async, App};

/// What produced a derived alignment, recorded in `Alignment::derivation`.
///
/// Spelled `realign:<backend>-<preset>` so the row says both *that* it was realigned and *how* —
/// a re-run under a different backend or preset is a different artifact, and support questions
/// start with which one produced a given file.
pub fn derivation_tag(backend: &str, preset: &str) -> String {
    format!("realign:{backend}-{preset}")
}

impl App {
    /// Register the output of a realignment as a new alignment derived from `source_id`.
    ///
    /// `cram` is stage C's output and `reference` the FASTA it was compressed against — both are
    /// required, because a CRAM without its reference is unreadable rather than merely awkward.
    ///
    /// `backend` and `preset` are taken separately rather than as a ready-made derivation string
    /// so the recorded format stays authoritative here; a caller cannot invent its own spelling
    /// that later queries then fail to recognise.
    pub async fn register_realigned_alignment(
        &self,
        source_id: i64,
        cram: &Path,
        reference: &Path,
        reference_build: &str,
        backend: &str,
        preset: &str,
    ) -> Result<Alignment, AppError> {
        let source = self.alignment_or_err(source_id).await?;

        // Realigning to the build a sample is already on produces a second copy of the same
        // information at hours of cost. The design calls for refusing it; doing so here rather
        // than in the UI means the CLI and any future caller are covered by the same rule.
        if builds_match(&source.reference_build, reference_build) {
            return Err(AppError::Import(format!(
                "alignment #{source_id} is already on {}; realigning it to the same build would \
                 duplicate it",
                source.reference_build
            )));
        }

        if !cram.is_file() {
            return Err(AppError::AlignmentFileMissing {
                id: source_id,
                path: cram.display().to_string(),
            });
        }

        // Hash now rather than lazily. The file was just written, so it is in the page cache and
        // this is nearly free; leaving it for the first analysis would pay the read twice.
        let content_sha256 = sha256_file_async(cram.to_path_buf()).await?;

        let created = self
            .record_alignment(NewAlignment {
                // The same library — only the mapping changed.
                sequence_run_id: source.sequence_run_id,
                reference_build: reference_build.to_string(),
                aligner: backend.to_string(),
                // Realignment maps reads; it calls no variants.
                variant_caller: None,
                bam_path: Some(cram.to_string_lossy().into_owned()),
                reference_path: Some(reference.to_string_lossy().into_owned()),
                content_sha256: Some(content_sha256),
                derived_from_alignment_id: Some(source_id),
                derivation: Some(derivation_tag(backend, preset)),
            })
            .await?;

        Ok(created)
    }

    /// The alignments derived from `source_id`, if any.
    ///
    /// The UI asks this before offering "Realign": a sample that has already been realigned should
    /// be told so rather than silently given a second copy.
    pub async fn derived_alignments(&self, source_id: i64) -> Result<Vec<Alignment>, AppError> {
        let source = self.alignment_or_err(source_id).await?;
        let siblings = alignment::list_for_run(self.store.pool(), source.sequence_run_id).await?;
        Ok(siblings
            .into_iter()
            .filter(|a| a.derived_from_alignment_id == Some(source_id))
            .collect())
    }

    /// The alignment `id` was derived from, or `None` when it is an original.
    pub async fn derivation_source(&self, id: i64) -> Result<Option<Alignment>, AppError> {
        let aln = self.alignment_or_err(id).await?;
        match aln.derived_from_alignment_id {
            Some(parent) => Ok(alignment::get(self.store.pool(), parent).await?),
            None => Ok(None),
        }
    }

    /// Where stage C should write, for a realignment of `source_id` to `build`.
    ///
    /// Derived files live beside the workspace rather than next to the vendor's original: the
    /// source directory may be read-only, on removable media, or somewhere the user does not
    /// expect Navigator to write tens of GB.
    pub fn realigned_output_path(&self, source_id: i64, build: &str) -> PathBuf {
        navigator_domain::paths::decodingus_dir()
            .join("realigned")
            .join(format!("alignment-{source_id}.{build}.cram"))
    }
}

/// Whether two build names refer to the same reference for this purpose.
///
/// Compared case-insensitively on the recorded strings. Deliberately *not* normalised through
/// `canonical_build`: `chm13v2.0` and `chm13v2.0_maskedY_rCRS` share coordinates but differ in
/// chrM and in PAR masking, so realigning between them is a real operation rather than a no-op.
fn builds_match(a: &str, b: &str) -> bool {
    a.trim().eq_ignore_ascii_case(b.trim())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_derivation_tag_names_both_backend_and_preset() {
        assert_eq!(derivation_tag("minimap2", "sr"), "realign:minimap2-sr");
        assert_eq!(derivation_tag("minimap2", "map-hifi"), "realign:minimap2-map-hifi");
    }

    #[test]
    fn build_comparison_ignores_case_and_padding() {
        assert!(builds_match("chm13v2.0", "CHM13v2.0"));
        assert!(builds_match(" GRCh38 ", "grch38"));
        assert!(!builds_match("GRCh38", "chm13v2.0"));
    }

    /// The masked variant shares CHM13's coordinates but differs in chrM and PAR masking, so
    /// moving between them is a real realignment and must not be refused as a no-op.
    #[test]
    fn the_masked_chm13_variant_is_a_different_build() {
        assert!(!builds_match("chm13v2.0", "chm13v2.0_maskedY_rCRS"));
    }
}
