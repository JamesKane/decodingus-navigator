//! Workspace **chores** — the periodic batch jobs that keep a workspace current.
//!
//! These existed only as CLI subcommands (`private-y --project`, `rebuild-signatures --stale-tree`,
//! `publish-origins`), each noted in its own design doc as wanting a GUI trigger, and each noted as
//! wanting *the same answer* rather than a third bespoke button. This module is that answer: the
//! chores are named, surveyed and driven from one place, so the CLI and the GUI run identical code
//! and a fourth chore is a table entry rather than a new surface.
//!
//! **Surveying is deliberately on demand.** Two of the three cost real work to measure — one walks
//! every alignment, another fetches and parses a multi-MB haplotree — so nothing here runs off a
//! render path. The UI asks once, when a user asks it to.

use super::*;
use crate::fastpath::{chr_m_gvcf_for_alignment, chr_y_gvcf_for_alignment};

/// A workspace chore. The order is the order a fresh workspace wants them in: compute what the
/// cohort views need, re-place anything a new tree invalidated, then publish.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize)]
pub enum Chore {
    /// Compute and cache each subject's private-Y bucket. Nothing cross-subject — shared unnamed
    /// variants, candidate branches — can mean anything until this has been walked once.
    PrivateY,
    /// Re-place subjects whose calls were scored against a superseded haplotree, or whose derived
    /// consensus names a branch this tree no longer carries.
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
    /// Items it would consider — `due` of `total` is what makes "0 due" readable as "nothing to
    /// do" rather than "nothing found".
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
    /// One line for the status bar — chore-specific, since "12 done" means different things.
    pub summary: String,
}

/// What [`App::replace_against_current_tree`] did for one subject. Per-alignment call failures are
/// counted rather than raised: an alignment whose file is gone is a superseded vendor download, not
/// a reason to leave the subject un-replaced.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct TreeReplace {
    pub calls_replaced: usize,
    pub calls_failed: usize,
    /// Alignments skipped because their file is gone. Kept apart from `calls_failed` so a workspace
    /// whose vendor downloads have been cleaned out does not report a wall of errors for the one
    /// outcome that is expected and harmless.
    pub calls_skipped: usize,
    pub profiles_rebuilt: usize,
}

impl TreeReplace {
    /// Tally one per-alignment call, sorting "its file is gone" out of the failures.
    fn record(&mut self, outcome: Result<(), AppError>) {
        match outcome {
            Ok(()) => self.calls_replaced += 1,
            Err(e) if e.is_missing_alignment_file() => self.calls_skipped += 1,
            Err(_) => self.calls_failed += 1,
        }
    }
}

/// Per-subject result of a private-Y refresh.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct PrivateYRefresh {
    pub computed: usize,
    pub skipped: usize,
    pub failed: usize,
    /// Alignments whose file is gone — a superseded vendor download, not a computation failure.
    /// Counted apart so a missing file never reads as an error.
    pub missing_file: usize,
    pub novel: usize,
}

impl App {
    /// Survey every chore. Runs the real selectors, so it costs what the chores cost to *decide*
    /// (not to do) — call it from a button, never from a paint.
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

        // Stale placements: two independent symptoms, unioned. Neither subsumes the other — a
        // consensus is derived and persisted separately, so it rots while every call beneath it
        // stays current.
        let stale = match self.stale_tree_targets(false).await {
            Ok(v) => ChoreSurvey {
                chore: Chore::StaleTree,
                due: v.len(),
                total: self.list_all_biosamples().await?.len(),
                blocked: None,
            },
            // No tree cached yet, or it could not be parsed: report the reason rather than zero.
            Err(e) => ChoreSurvey {
                chore: Chore::StaleTree,
                due: 0,
                total: 0,
                blocked: Some(e.to_string()),
            },
        };
        out.push(stale);

        // Publishing needs an account: without one there is nowhere to publish to.
        let publishable = mdka::publishable(self.store.pool(), Lineage::Y.as_str()).await?.len();
        out.push(ChoreSurvey {
            chore: Chore::PublishOrigins,
            due: publishable,
            total: mdka::count_for_lineage(self.store.pool(), Lineage::Y.as_str()).await?,
            blocked: self.current_account().is_none().then(|| "not signed in".to_string()),
        });

        Ok(out)
    }

    /// Subjects due for re-placement: those whose *source calls* carry another tree's fingerprint,
    /// unioned with those whose *derived consensus* names a branch this tree does not carry.
    ///
    /// `include_unknown` also takes calls that predate the fingerprint field — 80% of them, mostly
    /// BAM re-walks — so it is opt-in rather than the default.
    pub async fn stale_tree_targets(&self, include_unknown: bool) -> Result<Vec<SampleGuid>, AppError> {
        let by_fingerprint = self.subjects_placed_against_another_tree(include_unknown).await?;
        let off_tree = self.subjects_labelled_off_tree().await?;
        let mut set: std::collections::HashSet<SampleGuid> = by_fingerprint.into_iter().collect();
        set.extend(off_tree);
        let mut v: Vec<SampleGuid> = set.into_iter().collect();
        v.sort_by_key(|g| g.0);
        Ok(v)
    }

    /// Re-place one subject against the current tree: its per-alignment Y/mt calls first, then the
    /// pooled profiles built from them.
    ///
    /// The per-alignment step is the half that was missing. Rebuilding only the profiles refreshes
    /// `consensus_profile` while leaving every `haplogroup_call` row untouched, so a call carrying a
    /// name from an older tree survives every sweep — `1087` held `aln:903` reading
    /// `CP086569.2:27785335 G->A` against a current `aln:864` of `R-BY66248`, which is what the Y
    /// card reported as "sources diverge below root". Worse, the sweep selects subjects *by* those
    /// call fingerprints ([`Self::stale_tree_targets`]), so a subject it never corrected stayed due
    /// forever and the count never fell.
    ///
    /// Order matters: the calls are the profile's input, so re-placing them after the build would
    /// leave the profile a version behind. Each step is best-effort and independent — an alignment
    /// whose file is gone must not stop the rest of the subject from being brought current.
    ///
    /// Re-scoring is guarded by the alignment's own fingerprint (file hash + tree hash), so a
    /// subject already current costs a fingerprint comparison rather than a walk. When the tree
    /// genuinely changed, doing that work *is* the chore.
    pub async fn replace_against_current_tree(&self, guid: SampleGuid) -> Result<TreeReplace, AppError> {
        let mut r = TreeReplace::default();
        for aln in self.list_alignments_for_biosample(guid).await.unwrap_or_default() {
            // Fast-path (sidecar GVCF) calls first, and only then the CRAM-walk assignment.
            //
            // These are the `external`-provenance rows, but "external" here means the pipeline's
            // caller, not an authority we merely relay: `assign_y_from_gvcf` places the GVCF's
            // calls against *our* tree, so the stored name is only as current as the tree cached
            // the day it was imported — `altai363p` carried `chrY:5216846A>C [Node721]` from one.
            // They are ours to re-derive. Skipping them would leave the very rows the internal
            // assignment defers to (see `has_preferred_external_call`) as the only stale ones left,
            // which is the case that prompted this.
            // Prefer the recorded paths; fall back to locating the GVCFs beside the alignment.
            // The fallback is what makes this work on the existing corpus at all — every alignment
            // imported before the paths were recorded has none, so keying solely on the record
            // would have fixed only future imports and left the subjects that prompted this
            // permanently stale.
            let recorded = self.recorded_sidecars(aln.id).await.ok().flatten();
            let y_gvcf = recorded
                .as_ref()
                .and_then(|s| s.chr_y_gvcf.clone())
                .or_else(|| chr_y_gvcf_for_alignment(&aln));
            let m_gvcf = recorded
                .as_ref()
                .and_then(|s| s.chr_m_gvcf.clone())
                .or_else(|| chr_m_gvcf_for_alignment(&aln));
            // A GVCF that has since gone (superseded vendor download, unmounted volume) is skipped
            // rather than counted against the subject.
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
            // The CRAM-walk assignments. An alignment whose file has been removed since import
            // reports `AlignmentFileMissing` here and is skipped — the sidecar calls above may still
            // have re-placed the subject perfectly well without it.
            for outcome in [
                self.assign_y_haplogroup(aln.id).await.map(|_| ()),
                self.assign_mtdna_haplogroup_from_alignment(aln.id).await.map(|_| ()),
            ] {
                r.record(outcome);
            }
        }
        // Both arms rebuild regardless: a source-less lineage just yields an empty profile, and a
        // subject can be stale on one arm while current on the other.
        self.build_y_profile(guid).await?;
        self.build_mt_profile(guid).await?;
        r.profiles_rebuilt = 2;
        Ok(r)
    }

    /// Refresh one subject's private-Y: every alignment, or — when there is none — every
    /// non-chip variant set.
    ///
    /// The variant-set arm is not a fallback for tidiness: most members of a real Y project have no
    /// alignment at all, and until it existed they had no private-Y, which is what made candidate
    /// branches inert. Lifted out of the CLI so the GUI runs the same code rather than a second
    /// implementation that drifts.
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
            // A row whose file is gone is not a computation failure — reporting it as one buries
            // the real errors.
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
