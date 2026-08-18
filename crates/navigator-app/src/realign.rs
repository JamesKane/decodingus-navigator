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
//! A CRAM can not be read without the reference it was compressed against, so `reference_path` is
//! recorded on the row rather than resolved by convention later. An alignment whose reference has
//! moved is unreadable, and the row is the only place that knows which one it was.

use std::path::{Path, PathBuf};

use navigator_domain::du_domain::ids::SampleGuid;
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
    /// so the recorded format stays authoritative here; a caller can not invent its own spelling
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

    /// The subject an alignment belongs to, via its sequencing run.
    ///
    /// The UI needs this to say whose realignment is running. Without it the running card matched on
    /// alignment id alone, and a page showing subject A during a job on subject B told A their
    /// genome was being rebuilt — the ownership question has to be answered where the mapping from
    /// alignment to subject actually lives.
    pub async fn subject_of_alignment(&self, id: i64) -> Result<Option<SampleGuid>, AppError> {
        let Some(aln) = alignment::get(self.store.pool(), id).await? else {
            return Ok(None);
        };
        Ok(
            navigator_store::sequence_run::get(self.store.pool(), aln.sequence_run_id)
                .await?
                .map(|run| run.biosample_guid),
        )
    }

    /// The alignment `id` was derived from, or `None` when it is an original.
    pub async fn derivation_source(&self, id: i64) -> Result<Option<Alignment>, AppError> {
        let aln = self.alignment_or_err(id).await?;
        match aln.derived_from_alignment_id {
            Some(parent) => Ok(alignment::get(self.store.pool(), parent).await?),
            None => Ok(None),
        }
    }

    /// The alignments in `project_id` that a realignment to `target_build` would actually act on.
    ///
    /// Computed up front rather than discovered while running, because the honest thing to tell
    /// someone before a job measured in *days* is how many samples it covers. Skips what would be
    /// refused anyway — anything already on the target build, anything already realigned, and
    /// anything with no file to read — so the count is the real one rather than an upper bound.
    pub async fn realignable_in_project(&self, project_id: i64, target_build: &str) -> Result<Vec<i64>, AppError> {
        // One query for the whole project rather than one per member — the same idiom
        // `project_report` uses on this very tab. Measured on a 2,504-member project: 2.7 ms for the
        // grouped query against 17.7 ms for the per-member loop, and this runs twice per batch.
        let guids: Vec<_> = self
            .list_biosamples(project_id)
            .await?
            .into_iter()
            .map(|s| s.guid)
            .collect();
        let rows = navigator_store::alignment::list_for_biosamples(self.store.pool(), &guids).await?;

        // Grouped by subject, not flattened: the rule's "already realigned" condition asks whether
        // anything *in that subject's own set* was derived from a given alignment, so it has to see
        // one subject's alignments at a time.
        let mut by_subject: std::collections::HashMap<_, Vec<Alignment>> = std::collections::HashMap::new();
        for (guid, alignment) in rows {
            by_subject.entry(guid).or_default().push(alignment);
        }

        let mut out = Vec::new();
        for guid in &guids {
            if let Some(alignments) = by_subject.get(guid) {
                out.extend(realignable_for_subject(alignments, target_build));
            }
        }
        Ok(out)
    }

    /// The cached FASTA for `build`, if it has already been fetched.
    ///
    /// Public because the realignment job needs the reference *before* it starts: mapping to a
    /// build whose FASTA is not cached would otherwise stall at the index stage while gigabytes
    /// download, with nothing on screen explaining the wait.
    pub fn cached_reference_path(&self, build: &str) -> Option<PathBuf> {
        self.gateway.cached_reference(build)
    }

    /// Where stage C should write, for a realignment of `source_id` to `build`.
    ///
    /// Derived files live beside the workspace rather than next to the vendor's original: the
    /// source directory may be read-only, on removable media, or somewhere the user does not
    /// expect Navigator to write tens of GB.
    pub fn realigned_output_path(&self, source_id: i64, build: &str) -> PathBuf {
        navigator_domain::paths::decodingus_dir()
            .join("realigned")
            .join(format!("alignment-{source_id}.{build}.bam"))
    }
}

/// The build realignment targets when nothing says otherwise — the complete assembly, which is the
/// only reason the module exists.
pub const DEFAULT_TARGET_BUILD: &str = "chm13v2.0";

/// Whether `build` is the realignment target — the complete assembly.
///
/// `pub` so the UI can ask the question rather than spelling out its own comparison. Both Advanced
/// realign cards used to do the latter, with `eq_ignore_ascii_case` and no trim, which is a subtly
/// different rule from the one the job enforces.
pub fn is_target_build(build: &str) -> bool {
    builds_match(build, DEFAULT_TARGET_BUILD)
}

/// Which of one subject's `alignments` a realignment to `target_build` would actually act on.
///
/// Takes the whole list rather than one alignment because two of the four conditions are about the
/// *set*: an alignment is skipped if something in the list was already derived from it. Ordered as
/// the input is.
///
/// Shared by the project-wide count and the Simple-mode single-subject offer, so that the two
/// can not drift apart. If they disagreed, a user would be told a batch covers a sample it then
/// silently skips — or be offered four hours of work that the job itself would refuse.
pub(crate) fn realignable_for_subject(alignments: &[Alignment], target_build: &str) -> Vec<i64> {
    alignments
        .iter()
        .filter(|a| a.bam_path.is_some() && !a.is_derived())
        .filter(|a| !builds_match(&a.reference_build, target_build))
        // Already realigned by an earlier run: its output is sitting in this same list.
        .filter(|a| !alignments.iter().any(|d| d.derived_from_alignment_id == Some(a.id)))
        .map(|a| a.id)
        .collect()
}

/// Whether two build names refer to the same reference for this purpose.
///
/// Compared case-insensitively on the recorded strings, after trimming. Deliberately *not*
/// normalised through a `canonical_build`: `chm13v2.0` and `chm13v2.0_maskedY_rCRS` share
/// coordinates but differ in chrM and in PAR masking, so realigning between them is a real
/// operation rather than a no-op.
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

    /// Build one alignment row. Only the five fields the rule looks at are meaningful.
    fn aln(id: i64, build: &str, has_file: bool, derived_from: Option<i64>) -> Alignment {
        Alignment {
            id,
            sequence_run_id: 1,
            reference_build: build.to_string(),
            aligner: "bwa-mem2".into(),
            variant_caller: None,
            bam_path: has_file.then(|| format!("/data/{id}.cram")),
            reference_path: None,
            content_sha256: None,
            derived_from_alignment_id: derived_from,
            derivation: derived_from.map(|_| "realign:minimap2-sr".to_string()),
        }
    }

    #[test]
    fn an_original_on_an_older_build_is_realignable() {
        let list = [aln(1, "GRCh38", true, None)];
        assert_eq!(realignable_for_subject(&list, "chm13v2.0"), vec![1]);
    }

    /// The three ways an alignment disqualifies itself, each on its own.
    #[test]
    fn already_on_target_without_a_file_or_itself_derived_are_all_skipped() {
        assert!(realignable_for_subject(&[aln(1, "chm13v2.0", true, None)], "chm13v2.0").is_empty());
        assert!(realignable_for_subject(&[aln(1, "GRCh38", false, None)], "chm13v2.0").is_empty());
        assert!(realignable_for_subject(&[aln(1, "GRCh38", true, Some(9))], "chm13v2.0").is_empty());
    }

    /// The set-level condition: once its realigned output sits beside it, the source is done. Without
    /// this the offer would reappear after a four-hour job and invite an identical second one.
    #[test]
    fn a_source_that_has_already_been_realigned_is_not_offered_again() {
        let list = [aln(1, "GRCh38", true, None), aln(2, "chm13v2.0", true, Some(1))];
        assert!(realignable_for_subject(&list, "chm13v2.0").is_empty());
    }

    /// A subject holding several originals gets each of them, in input order — the caller picks.
    #[test]
    fn every_qualifying_original_is_returned_in_order() {
        let list = [
            aln(5, "GRCh37", true, None),
            aln(3, "chm13v2.0", true, None),
            aln(8, "GRCh38", true, None),
        ];
        assert_eq!(realignable_for_subject(&list, "chm13v2.0"), vec![5, 8]);
    }

    /// The masked CHM13 variant is a different build, so it is still worth realigning to plain
    /// CHM13 — the same distinction `builds_match` draws above, reaching the rule that uses it.
    #[test]
    fn the_masked_variant_is_still_realignable_to_plain_chm13() {
        let list = [aln(1, "chm13v2.0_maskedY_rCRS", true, None)];
        assert_eq!(realignable_for_subject(&list, "chm13v2.0"), vec![1]);
    }

    #[test]
    fn the_target_build_is_recognised_however_it_is_spelled() {
        assert!(is_target_build("chm13v2.0"));
        assert!(is_target_build(" CHM13v2.0 "));
        assert!(!is_target_build("GRCh38"));
        // The masked variant is a different reference, so a subject holding only that one has not
        // yet had their Y read against plain CHM13.
        assert!(!is_target_build("chm13v2.0_maskedY_rCRS"));
    }
}
