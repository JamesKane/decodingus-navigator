//! Workspace **chores**. A chore is a batch job that keeps a workspace current, and the user runs
//! it from time to time.
//!
//! Each chore was a CLI subcommand only. They are `private-y --project`,
//! `rebuild-signatures --stale-tree`, and `publish-origins`. The design document of each chore
//! asked for a trigger in the GUI. Each document also asked for *one* answer, not a third separate
//! button.
//!
//! This module is that answer. It names each chore, surveys them, and runs them from one place.
//! The CLI and the GUI then run the same code. A fourth chore is a new row in a table, and not a
//! new screen.
//!
//! **The survey runs only at the request of the user, by design.** Two of the three chores cost
//! real work to measure. One walks each alignment, and another reads and parses a haplotree of many
//! MB. So no code here runs during a paint. The UI asks one time, when the user presses the
//! button.

use super::*;
use crate::fastpath::{chr_m_gvcf_for_alignment, chr_y_gvcf_for_alignment};

/// A workspace chore. The order is the order a fresh workspace wants them in: compute what the
/// cohort views need, re-place anything a new tree invalidated, then publish.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize)]
pub enum Chore {
    /// Calculate the private-Y group of each subject and write it to the cache. A result that
    /// covers more than one subject is not correct before this walk runs one time. Such results are
    /// the unnamed variants that subjects share, and the candidate branches.
    PrivateY,
    /// Place a subject again when an old haplotree scored its calls. Also place a subject again
    /// when its consensus names a branch that this tree no longer holds.
    StaleTree,
    /// Publish MDKA ancestral origins for the subjects the consent predicate allows.
    PublishOrigins,
}

impl Chore {
    pub const ALL: [Chore; 3] = [Chore::PrivateY, Chore::StaleTree, Chore::PublishOrigins];

    /// Stable key for i18n lookups and logs.
    pub fn key(self) -> &'static str {
        match self {
            Chore::PrivateY => "privateY",
            Chore::StaleTree => "staleTree",
            Chore::PublishOrigins => "publishOrigins",
        }
    }
}

/// What a chore would do if run now.
#[derive(Debug, Clone, PartialEq, serde::Serialize)]
pub struct ChoreSurvey {
    pub chore: Chore,
    /// Items the chore would act on.
    pub due: usize,
    /// The items that the chore would examine. The pair `due` of `total` lets the user read "0
    /// due" as "there is no work". Without `total`, that text can also mean "the survey found
    /// nothing".
    pub total: usize,
    /// Why the chore can not run at all (not signed in, no tree). `Some` disables it.
    pub blocked: Option<String>,
}

/// What a chore did.
#[derive(Debug, Clone, Default, PartialEq, serde::Serialize)]
pub struct ChoreOutcome {
    pub done: usize,
    pub skipped: usize,
    pub failed: usize,
    /// One line of text for the status bar. Each chore writes its own text, because "12 done"
    /// refers to a different item in each chore.
    pub summary: String,
}

/// The work that [`App::replace_against_current_tree`] did for one subject.
///
/// The code counts a call failure of one alignment. It does not stop on that failure. An alignment
/// with an absent file is an old vendor download. It is not a reason to leave the subject in its
/// earlier place.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct TreeReplace {
    pub calls_replaced: usize,
    pub calls_failed: usize,
    /// The count of alignments that the chore skipped, because the file of each one is absent.
    ///
    /// This count is separate from `calls_failed`. A user can remove the old vendor downloads of a
    /// workspace. That workspace must not then report many errors for a result that is normal and
    /// harmless.
    pub calls_skipped: usize,
    pub profiles_rebuilt: usize,
}

impl TreeReplace {
    /// Count the call of one alignment. Put an absent file in its own group, and not with the
    /// failures.
    fn record(&mut self, outcome: Result<(), AppError>) {
        match outcome {
            Ok(()) => self.calls_replaced += 1,
            Err(e) if e.is_missing_alignment_file() => self.calls_skipped += 1,
            Err(_) => self.calls_failed += 1,
        }
    }
}

/// The result of a private-Y refresh, for one subject.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct PrivateYRefresh {
    pub computed: usize,
    pub skipped: usize,
    pub failed: usize,
    /// The count of alignments with an absent file. Such a file is an old vendor download. It is
    /// not a fault in the calculation. The count is separate, so an absent file never looks like an
    /// error.
    pub missing_file: usize,
    pub novel: usize,
}

impl App {
    /// Survey each chore. The method runs the real selectors. So it costs the work that a chore
    /// needs to *decide*, and not the work to do that chore. Call this method from a button. Never
    /// call it during a paint.
    pub async fn maintenance_survey(&self) -> Result<Vec<ChoreSurvey>, AppError> {
        let mut out = Vec::with_capacity(Chore::ALL.len());

        // Private-Y: alignments and variant sets with no cached bucket.
        let alignments = self.list_all_alignments().await?;
        let mut due = 0usize;
        for a in &alignments {
            if matches!(self.cached_private_y(a.id).await, Ok(None)) {
                due += 1;
            }
        }
        out.push(ChoreSurvey {
            chore: Chore::PrivateY,
            due,
            total: alignments.len(),
            blocked: None,
        });

        // An old placement has two separate symptoms, and the code joins the two sets. One set
        // does not hold the other. The app derives a consensus and stores it separately. So a
        // consensus becomes old while each call below it
        // stays current.
        let stale = match self.stale_tree_targets(false).await {
            Ok(v) => ChoreSurvey {
                chore: Chore::StaleTree,
                due: v.len(),
                total: self.list_all_biosamples().await?.len(),
                blocked: None,
            },
            // The cache holds no tree, or the parser refused it. Report the reason, not zero.
            Err(e) => ChoreSurvey {
                chore: Chore::StaleTree,
                due: 0,
                total: 0,
                blocked: Some(e.to_string()),
            },
        };
        out.push(stale);

        // A publish needs an account. Without an account there is no destination.
        let publishable = mdka::publishable(self.store.pool(), Lineage::Y.as_str()).await?.len();
        out.push(ChoreSurvey {
            chore: Chore::PublishOrigins,
            due: publishable,
            total: mdka::count_for_lineage(self.store.pool(), Lineage::Y.as_str()).await?,
            blocked: self.current_account().is_none().then(|| "not signed in".to_string()),
        });

        Ok(out)
    }

    /// The subjects that need a new placement. The method joins two sets.
    ///
    /// In the first set, the *source calls* of a subject carry the fingerprint of another tree. In
    /// the second set, the *derived consensus* of a subject names a branch that this tree does not
    /// hold.
    ///
    /// `include_unknown` adds the calls that are older than the fingerprint field. Those calls are
    /// 80% of the total, and most of them are BAM walks. So the user must select this option, and
    /// it is not the default.
    pub async fn stale_tree_targets(&self, include_unknown: bool) -> Result<Vec<SampleGuid>, AppError> {
        let by_fingerprint = self.subjects_placed_against_another_tree(include_unknown).await?;
        let off_tree = self.subjects_labelled_off_tree().await?;
        let mut set: std::collections::HashSet<SampleGuid> = by_fingerprint.into_iter().collect();
        set.extend(off_tree);
        let mut v: Vec<SampleGuid> = set.into_iter().collect();
        v.sort_by_key(|g| g.0);
        Ok(v)
    }

    /// Place one subject against the current tree again. The method scores the Y calls and the mt
    /// calls of each alignment first. It then builds the pooled profiles from those calls.
    ///
    /// The step for each alignment was the half that was absent. A rebuild of only the profiles
    /// refreshes `consensus_profile` and changes no `haplogroup_call` row. So a call with a name
    /// from an older tree stays after each sweep.
    ///
    /// Subject `1087` showed this fault. Its `aln:903` row read `CP086569.2:27785335 G->A`, and its
    /// current `aln:864` row read `R-BY66248`. The Y card then reported "sources diverge below
    /// root".
    ///
    /// The fault was worse than one wrong card. The sweep selects a subject *by* those call
    /// fingerprints, in [`Self::stale_tree_targets`]. So a subject that the sweep never corrected
    /// stayed in the list for all time, and the count never became smaller.
    ///
    /// The order of the two steps matters. The calls are the input of the profile. A new placement
    /// after the build leaves the profile one version behind.
    ///
    /// Each step is independent, and a failure in one step does not stop the others. An alignment
    /// with an absent file must not stop the work on the rest of the subject.
    ///
    /// The fingerprint of the alignment, which is the file hash with the tree hash, guards the new
    /// score. So a subject that is already current costs one comparison and not a walk. When the
    /// tree did change, that walk *is* the chore.
    pub async fn replace_against_current_tree(&self, guid: SampleGuid) -> Result<TreeReplace, AppError> {
        let mut r = TreeReplace::default();
        for aln in self.list_alignments_for_biosample(guid).await.unwrap_or_default() {
            // Do the fast-path calls, which come from a sidecar GVCF, first. Do the CRAM-walk
            // assignment after them.
            //
            // These rows have `external` provenance. Here "external" names the caller of the
            // pipeline. It does not name an authority that the app only relays.
            // `assign_y_from_gvcf` places the calls of the GVCF against *our* tree. So the stored
            // name is only as current as the tree in the cache on the day of the import. Subject
            // `altai363p` carried `chrY:5216846A>C [Node721]` from such an import.
            //
            // The app must derive these rows again. Without this step, the internal assignment
            // defers to those rows, in `has_preferred_external_call`, and they stay the only old
            // rows in the workspace. That case caused this code.
            //
            // Use the recorded paths first. If there are none, look for the GVCF files beside the
            // alignment. This second step is necessary for the corpus that exists today. Each
            // alignment from before the recorded-path change has no path. A key on the record
            // alone corrects only a future import, and it leaves the subjects that caused this
            // change old for all time.
            let recorded = self.recorded_sidecars(aln.id).await.ok().flatten();
            let y_gvcf = recorded
                .as_ref()
                .and_then(|s| s.chr_y_gvcf.clone())
                .or_else(|| chr_y_gvcf_for_alignment(&aln));
            let m_gvcf = recorded
                .as_ref()
                .and_then(|s| s.chr_m_gvcf.clone())
                .or_else(|| chr_m_gvcf_for_alignment(&aln));
            // The code skips a GVCF that is no longer on disk. The cause is an old vendor
            // download or a volume that the user removed. The code does not count it against the
            // subject.
            for outcome in [
                match y_gvcf.filter(|p| p.is_file()) {
                    Some(p) => Some(self.assign_y_from_gvcf(aln.id, &p).await.map(|_| ())),
                    None => None,
                },
                match m_gvcf.filter(|p| p.is_file()) {
                    Some(p) => Some(self.assign_mt_from_gvcf(aln.id, &p).await.map(|_| ())),
                    None => None,
                },
            ]
            .into_iter()
            .flatten()
            {
                r.record(outcome);
            }
            // The CRAM-walk assignments. An alignment with a file that the user removed after
            // the import gives `AlignmentFileMissing` here, and the code skips it. The sidecar
            // calls above can still place the subject correctly without that alignment.
            for outcome in [
                self.assign_y_haplogroup(aln.id).await.map(|_| ()),
                self.assign_mtdna_haplogroup_from_alignment(aln.id).await.map(|_| ()),
            ] {
                r.record(outcome);
            }
        }
        // The code rebuilds both lineages in each case. A lineage with no source gives an empty
        // profile. A subject can also be old on one lineage and current on the other.
        self.build_y_profile(guid).await?;
        self.build_mt_profile(guid).await?;
        r.profiles_rebuilt = 2;
        Ok(r)
    }

    /// Refresh the private-Y of one subject. The method uses each alignment. When the subject has
    /// no alignment, it uses each variant set that is not a chip.
    ///
    /// The variant-set path is necessary, and it is not only for a clean design. Most members of a
    /// real Y project have no alignment. Before this path, those members had no private-Y, and that
    /// gap made each candidate branch inert.
    ///
    /// This code came out of the CLI, so the GUI runs the same code. A second copy would become
    /// different over time.
    pub async fn refresh_private_y(&self, guid: SampleGuid, force: bool) -> Result<PrivateYRefresh, AppError> {
        let mut r = PrivateYRefresh::default();
        let alignments = self.list_alignments_for_biosample(guid).await.unwrap_or_default();
        if alignments.is_empty() {
            for set in self
                .list_variant_sets(guid)
                .await
                .unwrap_or_default()
                .iter()
                .filter(|s| s.source_type != SourceType::Chip)
            {
                match self.private_y_from_variant_set(set).await {
                    Ok(bucket) => {
                        r.computed += 1;
                        r.novel += bucket.novel_in_unique_sequence();
                    }
                    Err(_) => r.failed += 1,
                }
            }
            return Ok(r);
        }
        for a in &alignments {
            // A row with an absent file is not a fault in the calculation. A report of it as a
            // fault hides the real errors.
            if !a.bam_path.as_deref().is_some_and(|p| std::path::Path::new(p).exists()) {
                r.missing_file += 1;
                continue;
            }
            if !force && matches!(self.cached_private_y(a.id).await, Ok(Some(_))) {
                r.skipped += 1;
                continue;
            }
            match self.private_y_variants_self_masked(a.id).await {
                Ok(bucket) => {
                    r.computed += 1;
                    r.novel += bucket.novel_in_unique_sequence();
                }
                Err(_) => r.failed += 1,
            }
        }
        Ok(r)
    }
}
