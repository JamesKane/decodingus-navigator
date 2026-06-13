//! Y-variant profile — multi-source concordance over a subject's Y-SNP calls.
//!
//! Phase 1 (on-demand, no persistence): each Y-bearing source (a WGS alignment's haplogroup
//! placement, the combined chip/BISDNA placement, the private-Y bucket) contributes per-SNP calls;
//! [`reconcile_y`] groups them **by SNP name** (build-independent — M269 is M269 whether the source
//! aligned to GRCh37 or GRCh38) and weight-votes derived vs ancestral using each source's
//! [`SourceType::snp_weight`]. The result classifies each SNP as confirmed / novel / conflict /
//! single-source with per-source provenance. Mirrors the Scala `YVariantConcordance`, restricted to
//! the SNP axis (STR concordance is served separately by the Y-STR report).

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::variants::SourceType;

/// One source's call state at a Y position.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum YState {
    /// Carries the derived (mutant) allele — positive for the SNP's branch.
    Derived,
    /// Carries the ancestral (reference) allele.
    Ancestral,
    /// No confident call.
    NoCall,
}

/// Cross-source status of a Y SNP after reconciliation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum YVariantStatus {
    /// ≥2 sources agree on the consensus state and the SNP is a known tree variant.
    Confirmed,
    /// Derived but not a known tree variant (private / off-path).
    Novel,
    /// Sources disagree (weighted minority > 30%).
    Conflict,
    /// Only one source reports the SNP.
    SingleSource,
}

/// One source's observation of a SNP (for provenance display).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct YSourceObs {
    pub label: String,
    pub source_type: SourceType,
    pub state: YState,
}

/// A reconciled Y SNP across the subject's sources.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct YProfileVariant {
    /// SNP name (e.g. "M269"); for unnamed/novel calls this is a `@<position>` placeholder.
    pub name: String,
    /// A representative position (from the consensus-side sources; builds may differ).
    pub position: i64,
    pub ancestral: String,
    pub derived: String,
    pub consensus: YState,
    pub status: YVariantStatus,
    /// Sources matching the consensus state.
    pub support: usize,
    /// Sources with any call (excludes NoCall).
    pub total: usize,
    /// Whether the SNP is a known haplotree variant (vs a private/novel call).
    pub in_tree: bool,
    pub sources: Vec<YSourceObs>,
}

/// Per-status counts for the profile header.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct YProfileSummary {
    pub total: usize,
    pub confirmed: usize,
    pub novel: usize,
    pub conflict: usize,
    pub single_source: usize,
}

/// One source's call at a SNP, fed into [`reconcile_y`].
#[derive(Debug, Clone, PartialEq)]
pub struct YObsInput {
    pub name: String,
    pub position: i64,
    pub ancestral: String,
    pub derived: String,
    pub state: YState,
    /// Whether this SNP is a known tree variant (true for placement SNPs, false for private calls).
    pub in_tree: bool,
}

/// Fraction of disagreeing (weighted) support above which a SNP is a conflict.
const CONFLICT_FRACTION: f64 = 0.30;

/// Key a SNP for cross-source/cross-build grouping: by name when present (build-independent), else
/// by position (a novel/unnamed call only ever matches the same build's same position).
fn group_key(obs: &YObsInput) -> String {
    if obs.name.trim().is_empty() {
        format!("@{}", obs.position)
    } else {
        obs.name.trim().to_uppercase()
    }
}

/// Reconcile per-source Y SNP observations into one profile. `sources` is `(label, source_type,
/// observations)`; a source contributes at most one observation per SNP (the last wins on dup).
pub fn reconcile_y(sources: &[(String, SourceType, Vec<YObsInput>)]) -> Vec<YProfileVariant> {
    // Group all observations by SNP key, preserving the source they came from.
    struct Acc {
        repr: YObsInput,
        obs: Vec<(String, SourceType, YState)>, // (label, source_type, state)
    }
    let mut groups: BTreeMap<String, Acc> = BTreeMap::new();

    for (label, source_type, observations) in sources {
        for o in observations {
            let key = group_key(o);
            let acc = groups.entry(key).or_insert_with(|| Acc { repr: o.clone(), obs: Vec::new() });
            // Prefer a named, in-tree representative for display fields.
            if acc.repr.name.trim().is_empty() && !o.name.trim().is_empty() {
                acc.repr = o.clone();
            }
            acc.obs.push((label.clone(), *source_type, o.state));
        }
    }

    let mut out: Vec<YProfileVariant> = groups
        .into_values()
        .map(|acc| {
            let repr = acc.repr;
            // Weighted vote over non-NoCall observations.
            let (mut w_derived, mut w_ancestral) = (0.0f64, 0.0f64);
            let mut total = 0usize;
            for (_, st, state) in &acc.obs {
                match state {
                    YState::Derived => {
                        w_derived += st.snp_weight();
                        total += 1;
                    }
                    YState::Ancestral => {
                        w_ancestral += st.snp_weight();
                        total += 1;
                    }
                    YState::NoCall => {}
                }
            }

            let consensus = if total == 0 {
                YState::NoCall
            } else if w_derived >= w_ancestral {
                YState::Derived
            } else {
                YState::Ancestral
            };

            let support = acc
                .obs
                .iter()
                .filter(|(_, _, s)| *s == consensus && consensus != YState::NoCall)
                .count();

            let total_weight = w_derived + w_ancestral;
            let minority_weight = total_weight - w_derived.max(w_ancestral);
            let minority_fraction = if total_weight > 0.0 { minority_weight / total_weight } else { 0.0 };

            let status = if total == 0 {
                YVariantStatus::SingleSource // degenerate; no real call (shouldn't happen)
            } else if minority_fraction > CONFLICT_FRACTION {
                YVariantStatus::Conflict
            } else if consensus == YState::Derived && !repr.in_tree {
                // Derived off-tree call is novel/private — even from a single source (the common case).
                YVariantStatus::Novel
            } else if total == 1 {
                YVariantStatus::SingleSource
            } else {
                YVariantStatus::Confirmed
            };

            let sources = acc
                .obs
                .iter()
                .map(|(label, st, state)| YSourceObs { label: label.clone(), source_type: *st, state: *state })
                .collect();

            YProfileVariant {
                name: repr.name,
                position: repr.position,
                ancestral: repr.ancestral,
                derived: repr.derived,
                consensus,
                status,
                support,
                total,
                in_tree: repr.in_tree,
                sources,
            }
        })
        .collect();

    // Conflicts first (most actionable), then novel, then by name.
    out.sort_by(|a, b| status_rank(a.status).cmp(&status_rank(b.status)).then_with(|| a.name.cmp(&b.name)));
    out
}

fn status_rank(s: YVariantStatus) -> u8 {
    match s {
        YVariantStatus::Conflict => 0,
        YVariantStatus::Novel => 1,
        YVariantStatus::SingleSource => 2,
        YVariantStatus::Confirmed => 3,
    }
}

/// Per-status counts over a reconciled variant list.
pub fn summarize(variants: &[YProfileVariant]) -> YProfileSummary {
    let mut s = YProfileSummary { total: variants.len(), ..Default::default() };
    for v in variants {
        match v.status {
            YVariantStatus::Confirmed => s.confirmed += 1,
            YVariantStatus::Novel => s.novel += 1,
            YVariantStatus::Conflict => s.conflict += 1,
            YVariantStatus::SingleSource => s.single_source += 1,
        }
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    fn obs(name: &str, pos: i64, state: YState, in_tree: bool) -> YObsInput {
        YObsInput { name: name.into(), position: pos, ancestral: "A".into(), derived: "G".into(), state, in_tree }
    }

    #[test]
    fn two_sources_agree_in_tree_is_confirmed() {
        let v = reconcile_y(&[
            ("aln #1".into(), SourceType::WgsShortRead, vec![obs("M269", 100, YState::Derived, true)]),
            ("consumer".into(), SourceType::Chip, vec![obs("M269", 200, YState::Derived, true)]),
        ]);
        assert_eq!(v.len(), 1); // grouped by name across differing positions/builds
        assert_eq!(v[0].name, "M269");
        assert_eq!(v[0].consensus, YState::Derived);
        assert_eq!(v[0].status, YVariantStatus::Confirmed);
        assert_eq!(v[0].support, 2);
        assert_eq!(v[0].total, 2);
    }

    #[test]
    fn derived_not_in_tree_is_novel() {
        let v = reconcile_y(&[
            ("aln #1".into(), SourceType::WgsShortRead, vec![obs("FT1", 100, YState::Derived, false)]),
            ("aln #2".into(), SourceType::WgsShortRead, vec![obs("FT1", 100, YState::Derived, false)]),
        ]);
        assert_eq!(v[0].status, YVariantStatus::Novel);
    }

    #[test]
    fn comparable_weight_disagreement_is_conflict() {
        // WGS (0.85) derived vs Chip (0.5) ancestral → minority 0.5/1.35 ≈ 0.37 > 0.30 → conflict.
        let v = reconcile_y(&[
            ("aln #1".into(), SourceType::WgsShortRead, vec![obs("M269", 100, YState::Derived, true)]),
            ("consumer".into(), SourceType::Chip, vec![obs("M269", 100, YState::Ancestral, true)]),
        ]);
        assert_eq!(v[0].status, YVariantStatus::Conflict);
        assert_eq!(v[0].consensus, YState::Derived); // higher weight wins the consensus
    }

    #[test]
    fn dominant_weight_disagreement_is_not_conflict() {
        // Sanger (1.0) derived vs Manual (0.3) ancestral → minority 0.3/1.3 ≈ 0.23 ≤ 0.30 → confirmed.
        let v = reconcile_y(&[
            ("sanger".into(), SourceType::Sanger, vec![obs("M269", 100, YState::Derived, true)]),
            ("manual".into(), SourceType::Manual, vec![obs("M269", 100, YState::Ancestral, true)]),
        ]);
        assert_eq!(v[0].consensus, YState::Derived);
        assert_eq!(v[0].status, YVariantStatus::Confirmed);
        assert_eq!(v[0].support, 1); // only the Sanger source matches the derived consensus
    }

    #[test]
    fn single_source_is_single_source() {
        let v = reconcile_y(&[("aln #1".into(), SourceType::WgsShortRead, vec![obs("M269", 100, YState::Derived, true)])]);
        assert_eq!(v[0].status, YVariantStatus::SingleSource);
    }

    #[test]
    fn nocall_excluded_from_vote() {
        let v = reconcile_y(&[
            ("aln #1".into(), SourceType::WgsShortRead, vec![obs("M269", 100, YState::Derived, true)]),
            ("aln #2".into(), SourceType::WgsShortRead, vec![obs("M269", 100, YState::NoCall, true)]),
        ]);
        assert_eq!(v[0].total, 1); // NoCall not counted
        assert_eq!(v[0].status, YVariantStatus::SingleSource);
        assert_eq!(v[0].sources.len(), 2); // but still shown for provenance
    }

    #[test]
    fn summary_counts_by_status() {
        let v = reconcile_y(&[
            ("a".into(), SourceType::WgsShortRead, vec![obs("M269", 1, YState::Derived, true), obs("FT1", 2, YState::Derived, false)]),
            ("b".into(), SourceType::Chip, vec![obs("M269", 1, YState::Derived, true)]),
        ]);
        let s = summarize(&v);
        assert_eq!(s.total, 2);
        assert_eq!(s.confirmed, 1); // M269 (2 sources agree, in tree)
        assert_eq!(s.novel, 1); // FT1 (derived, not in tree → novel even single-source)
        assert_eq!(s.single_source, 0);
    }
}
