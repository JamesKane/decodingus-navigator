//! This module registers a realigned alignment. It is stage D of the realignment module.
//!
//! Stage C writes a CRAM file to disk. The code sorts that file and marks its duplicates. The
//! workspace does not yet hold it as an *alignment*. This module makes it one. Each analysis can
//! then run against it, and no analysis needs to know that the app made the file.
//!
//! ## This module only adds
//!
//! The code inserts the realigned row under the **same `SequenceRun`** as its source. That run is
//! the same physical library with a different map. The code does not change the source row.
//!
//! The sidecar fast path behaves in the same way, and this behaviour makes a realignment safe to
//! offer. A user can not lose data when they try it. A user can also delete a realigned alignment
//! and build it again from a source that did not change.
//!
//! ## The reference is part of the file
//!
//! A reader can not open a CRAM file without the reference that compressed it. So the row holds
//! `reference_path`, and no later code finds that path by a rule. A reader can not open an
//! alignment whose reference moved, and the row is the only record of that reference.

use std::path::{Path, PathBuf};

use navigator_domain::du_domain::ids::SampleGuid;
use navigator_domain::workspace::{Alignment, NewAlignment};
use navigator_store::alignment;

use crate::error::AppError;
use crate::{sha256_file_async, App};

/// The process that made a derived alignment. `Alignment::derivation` holds this value.
///
/// The format is `realign:<backend>-<preset>`. So the row shows *that* the app realigned the
/// alignment, and it also shows *how*.
///
/// A second run with another backend, or another preset, gives a different file. A support question
/// starts with the process that made the file.
pub fn derivation_tag(backend: &str, preset: &str) -> String {
    format!("realign:{backend}-{preset}")
}

impl App {
    /// Add the output of a realignment as a new alignment that comes from `source_id`.
    ///
    /// `cram` is the output of stage C. `reference` is the FASTA file that compressed it. The
    /// method needs both values. Without its reference, no reader can open a CRAM file.
    ///
    /// The method takes `backend` and `preset` as two values, and not as one derivation string.
    /// This module then controls the format. A caller can not write its own format, which a later
    /// query would not recognize.
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

        // A realignment to the build that a sample already uses gives a second copy of the same
        // data, and it costs many hours. The design refuses that work. This check is in the app
        // layer and not in the UI. So the CLI and each later caller follow the same rule.
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

        // Calculate the hash now, and not at the first use. The code wrote the file a moment ago,
        // so the page cache holds it and the read costs almost nothing. A hash at the first
        // analysis reads the file a second time.
        let content_sha256 = sha256_file_async(cram.to_path_buf()).await?;

        let created = self
            .record_alignment(NewAlignment {
                // The library is the same. Only the mapping changed.
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

    /// The alignments that come from `source_id`, when the workspace holds any.
    ///
    /// The UI calls this method before it offers the "Realign" action. The app must tell a user
    /// that a sample already has a realignment. It must not make a second copy with no message.
    pub async fn derived_alignments(&self, source_id: i64) -> Result<Vec<Alignment>, AppError> {
        let source = self.alignment_or_err(source_id).await?;
        let siblings = alignment::list_for_run(self.store.pool(), source.sequence_run_id).await?;
        Ok(siblings
            .into_iter()
            .filter(|a| a.derived_from_alignment_id == Some(source_id))
            .collect())
    }

    /// The subject of an alignment. The method finds it through the sequence run.
    ///
    /// The UI needs this value to name the subject of a realignment.
    ///
    /// Before this method, the progress card compared the alignment id alone. A page that showed
    /// subject A during a job on subject B then told subject A that the app rebuilt their genome.
    /// The code that maps an alignment to a subject must answer this question.
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

    /// The alignment that gave `id`, or `None` when `id` is an original alignment.
    pub async fn derivation_source(&self, id: i64) -> Result<Option<Alignment>, AppError> {
        let aln = self.alignment_or_err(id).await?;
        match aln.derived_from_alignment_id {
            Some(parent) => Ok(alignment::get(self.store.pool(), parent).await?),
            None => Ok(None),
        }
    }

    /// The alignments in `project_id` that a realignment to `target_build` would act on.
    ///
    /// The method finds them before the job starts. It does not find them during the job. A job of
    /// this size can run for *days*, and the app must tell the user how many samples it covers.
    ///
    /// The method skips each alignment that the job would refuse. Those are an alignment on the
    /// target build, an alignment with a realignment already, and an alignment with no file. So the
    /// count is the true count and not a maximum.
    pub async fn realignable_in_project(&self, project_id: i64, target_build: &str) -> Result<Vec<i64>, AppError> {
        // One query covers the full project. The code does not send one query for each member.
        // `project_report` uses the same method on this tab.
        //
        // A measurement on a project with 2,504 members gave 2.7 ms for the grouped query and
        // 17.7 ms for the loop. This code runs two times in each batch.
        let guids: Vec<_> = self
            .list_biosamples(project_id)
            .await?
            .into_iter()
            .map(|s| s.guid)
            .collect();
        let rows = navigator_store::alignment::list_for_biosamples(self.store.pool(), &guids).await?;

        // The code groups the rows by subject and does not make one flat list. The "already
        // realigned" condition asks whether an alignment *in the set of that subject* comes from a
        // given alignment. So the rule must read the alignments of one subject together.
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

    /// The FASTA file for `build` in the cache, when the app already downloaded it.
    ///
    /// This method is public because the realignment job needs the reference *before* it starts.
    /// Without the check, a map to a build with no FASTA file in the cache stops at the index
    /// stage. The download of some GB then runs, and the screen gives the user no reason for the
    /// wait.
    pub fn cached_reference_path(&self, build: &str) -> Option<PathBuf> {
        self.gateway.cached_reference(build)
    }

    /// The path where stage C writes, for a realignment of `source_id` to `build`.
    ///
    /// A derived file goes beside the workspace. It does not go beside the original file of the
    /// vendor. There are three reasons. A source directory can refuse a write. It can be on
    /// removable media. It can also be in a place where the user does not expect tens of GB from
    /// Navigator.
    pub fn realigned_output_path(&self, source_id: i64, build: &str) -> PathBuf {
        navigator_domain::paths::decodingus_dir()
            .join("realigned")
            .join(format!("alignment-{source_id}.{build}.bam"))
    }
}

/// The default target build of a realignment. It is the complete assembly, and that assembly is the
/// only reason for this module.
pub const DEFAULT_TARGET_BUILD: &str = "chm13v2.0";

/// Shows whether `build` is the target of a realignment, which is the complete assembly.
///
/// The function is `pub`, so the UI calls it and writes no comparison of its own. Both realign
/// cards of Advanced mode wrote their own comparison before. They used `eq_ignore_ascii_case` and
/// removed no space, and that rule is not the rule of the job.
pub fn is_target_build(build: &str) -> bool {
    builds_match(build, DEFAULT_TARGET_BUILD)
}

/// The alignments of one subject that a realignment to `target_build` would act on.
///
/// The function takes the full list and not one alignment. Two of the four conditions are about the
/// *set*. The function skips an alignment when another item in the list comes from it. The output
/// keeps the order of the input.
///
/// The project-wide count and the single-subject offer of Simple mode both call this function. So
/// the two can not become different.
///
/// A difference between them has two effects. The app can tell a user that a batch covers a sample,
/// and then skip that sample with no message. The app can also offer four hours of work that the
/// job then refuses.
pub(crate) fn realignable_for_subject(alignments: &[Alignment], target_build: &str) -> Vec<i64> {
    alignments
        .iter()
        .filter(|a| a.bam_path.is_some() && !a.is_derived())
        .filter(|a| !builds_match(&a.reference_build, target_build))
        // An earlier run already realigned this alignment, and its output is in this same
        // list.
        .filter(|a| !alignments.iter().any(|d| d.derived_from_alignment_id == Some(a.id)))
        .map(|a| a.id)
        .collect()
}

/// Shows whether two build names name the same reference, for this purpose.
///
/// The function removes the space at each end and then compares the two stored strings. It ignores
/// the case of a letter.
///
/// The function does *not* pass the names through `canonical_build`, by design. The builds
/// `chm13v2.0` and `chm13v2.0_maskedY_rCRS` use the same coordinates. But their chrM contigs are
/// different, and their PAR masks are different. So a realignment between the two does real
/// work.
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

    /// The masked build uses the CHM13 coordinates. But its chrM contig and its PAR mask are
    /// different. So a realignment between the two does real work, and the app must not refuse
    /// it.
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

    /// The condition on the set. When the realigned output is beside the source, the work on that
    /// source is complete. Without this rule, the offer returns after a four-hour job, and the user
    /// starts the same job again.
    #[test]
    fn a_source_that_has_already_been_realigned_is_not_offered_again() {
        let list = [aln(1, "GRCh38", true, None), aln(2, "chm13v2.0", true, Some(1))];
        assert!(realignable_for_subject(&list, "chm13v2.0").is_empty());
    }

    /// For a subject with more than one original alignment, the function returns each of them, in
    /// the order of the input. The caller then selects one.
    #[test]
    fn every_qualifying_original_is_returned_in_order() {
        let list = [
            aln(5, "GRCh37", true, None),
            aln(3, "chm13v2.0", true, None),
            aln(8, "GRCh38", true, None),
        ];
        assert_eq!(realignable_for_subject(&list, "chm13v2.0"), vec![5, 8]);
    }

    /// The masked CHM13 build is a different build. So a realignment to plain CHM13 still gives a
    /// result. `builds_match` above makes the same distinction, and this test covers the rule that
    /// calls it.
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
        // The masked build is a different reference. So the app did not yet read the Y chromosome
        // of a subject with only that build against plain CHM13.
        assert!(!is_target_build("chm13v2.0_maskedY_rCRS"));
    }
}
