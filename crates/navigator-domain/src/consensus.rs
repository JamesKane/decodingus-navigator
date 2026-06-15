//! Multi-source variant **consensus engine** — DNA-type-agnostic.
//!
//! Given a set of sources (a WGS alignment's placement, a chip/BISDNA panel, a private bucket, …),
//! each contributing per-variant calls keyed **by name** (build-independent — M269 is M269 whether
//! the source aligned to GRCh37 or GRCh38), [`reconcile`] groups them and weight-votes the consensus
//! state, classifying each variant as confirmed / novel / conflict / single-source and computing a
//! quality-weighted confidence. Mirrors the Scala `YVariantConcordance`.
//!
//! This engine is the shared foundation for the Y-DNA profile (the [`crate::yprofile`] adapter today)
//! and — by design — the future mtDNA (variants vs rCRS) and autosomal consumers. It carries no
//! DNA-type specifics: callers gather observations and supply the variant identity; the DNA type and
//! consensus label (haplogroup, where applicable) live at the persistence / app layer.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::variants::SourceType;

/// One source's call state at a variant position.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ConsensusState {
    /// Carries the derived (mutant) allele — positive for the variant's branch. For mtDNA this is
    /// "differs from rCRS"; for autosomes a future adapter maps a diploid genotype onto this axis.
    Derived,
    /// Carries the ancestral (reference) allele.
    Ancestral,
    /// No confident call.
    NoCall,
}

/// Cross-source status of a variant after reconciliation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ConsensusStatus {
    /// ≥2 sources agree on the consensus state and the variant is a known tree/reference variant.
    Confirmed,
    /// Derived but not a known tree variant (private / off-path).
    Novel,
    /// Sources disagree (weighted minority > 30%).
    Conflict,
    /// Only one source reports the variant.
    SingleSource,
    /// Has data but the weighted confidence is below the confirmation threshold without crossing the
    /// conflict line (rare — kept for parity with the Scala `YVariantConcordance`).
    Pending,
    /// No source made a confident call (every observation was NoCall).
    NoCoverage,
}

/// Per-position callability of a source's observation — scales its concordance weight (a base in a
/// no-coverage / poor-mapping region carries little confidence). Mirrors the Scala `YCallableState`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CallableState {
    Callable,
    LowCoverage,
    ExcessiveCoverage,
    PoorMappingQuality,
    NoCoverage,
    RefN,
}

impl CallableState {
    /// Confidence multiplier (Scala weights): full for CALLABLE, none for NO_COVERAGE / REF_N.
    pub fn weight(self) -> f64 {
        match self {
            CallableState::Callable => 1.0,
            CallableState::LowCoverage => 0.5,
            CallableState::ExcessiveCoverage | CallableState::PoorMappingQuality => 0.3,
            CallableState::NoCoverage | CallableState::RefN => 0.0,
        }
    }
}

/// One source's observation of a variant (for provenance display).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SourceObs {
    pub label: String,
    pub source_type: SourceType,
    pub state: ConsensusState,
}

/// A reconciled variant across the subject's sources.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ConsensusVariant {
    /// Variant name (e.g. "M269"); for unnamed/novel calls this is a `@<position>` placeholder.
    pub name: String,
    /// A representative position (from the consensus-side sources; builds may differ).
    pub position: i64,
    pub ancestral: String,
    pub derived: String,
    pub consensus: ConsensusState,
    pub status: ConsensusStatus,
    /// Sources matching the consensus state.
    pub support: usize,
    /// Sources with any call (excludes NoCall).
    pub total: usize,
    /// Whether the variant is a known reference/haplotree variant (vs a private/novel call).
    pub in_tree: bool,
    /// Weighted confidence in the consensus = consensusWeight / totalWeight (0 when no call).
    #[serde(default)]
    pub confidence_score: f64,
    pub sources: Vec<SourceObs>,
}

/// Per-status counts for the profile header.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize, Deserialize)]
pub struct ConsensusSummary {
    pub total: usize,
    pub confirmed: usize,
    pub novel: usize,
    pub conflict: usize,
    pub single_source: usize,
    /// Overall profile confidence: `(confirmed + 0.7·novel − 0.5·conflict) / total`, clamped [0,1].
    #[serde(default)]
    pub overall_confidence: f64,
}

/// One source's call at a variant, fed into [`reconcile`]. Quality fields refine the concordance
/// weight (see [`obs_weight`]); sources that don't carry them (chip, tree placement) leave them
/// `None` / `1.0` and fall back to the plain source-type weight.
#[derive(Debug, Clone, PartialEq)]
pub struct ConsensusObs {
    pub name: String,
    pub position: i64,
    pub ancestral: String,
    pub derived: String,
    pub state: ConsensusState,
    /// Whether this variant is a known tree/reference variant (true for placement SNPs, false for
    /// private calls).
    pub in_tree: bool,
    /// Read depth at the call (sequencing sources) — a `√depth/10` bonus, capped at +1.0.
    pub depth: Option<u32>,
    /// Mean mapping quality — an `MQ/60` factor, capped at 1.0.
    pub mapq: Option<f64>,
    /// Callability of the position — scales the weight (`NoCoverage`/`RefN` → 0).
    pub callable: Option<CallableState>,
    /// Region-confidence modifier (e.g. <1 in palindrome/amplicon zones), clamped [0.1, 1.0].
    pub region_modifier: f64,
}

impl ConsensusObs {
    /// A SNP/variant observation with no per-call quality data (weight = the source-type weight).
    /// Quality fields can be set afterward for sources that carry them (e.g. sequencing depth).
    pub fn snp(name: impl Into<String>, position: i64, ancestral: impl Into<String>, derived: impl Into<String>, state: ConsensusState, in_tree: bool) -> Self {
        ConsensusObs {
            name: name.into(),
            position,
            ancestral: ancestral.into(),
            derived: derived.into(),
            state,
            in_tree,
            depth: None,
            mapq: None,
            callable: None,
            region_modifier: 1.0,
        }
    }
}

/// Concordance weight for one observation (Scala `YVariantConcordance.calculateWeight`):
/// `snp_weight · (1 + min(√depth/10, 1)) · min(MQ/60, 1) · callableWeight · clamp(region, 0.1, 1)`.
/// Missing depth → no bonus; missing MQ/callable → factor 1.0.
pub fn obs_weight(source_type: SourceType, depth: Option<u32>, mapq: Option<f64>, callable: Option<CallableState>, region_modifier: f64) -> f64 {
    let method = source_type.snp_weight();
    let depth_bonus = depth.filter(|&d| d > 0).map(|d| ((d as f64).sqrt() / 10.0).min(1.0)).unwrap_or(0.0);
    let mapq_factor = mapq.filter(|&q| q > 0.0).map(|q| (q / 60.0).min(1.0)).unwrap_or(1.0);
    let callable_factor = callable.map(|c| c.weight()).unwrap_or(1.0);
    let region_factor = region_modifier.clamp(0.1, 1.0);
    method * (1.0 + depth_bonus) * mapq_factor * callable_factor * region_factor
}

/// Fraction of disagreeing (weighted) support above which a variant is a conflict.
const CONFLICT_FRACTION: f64 = 0.30;
/// Consensus confidence at or above which a multi-source, non-conflicting variant is confirmed.
const CONFIRMATION_FRACTION: f64 = 0.70;

/// Key a variant for cross-source/cross-build grouping: by name when present (build-independent),
/// else by position (a novel/unnamed call only ever matches the same build's same position).
fn group_key(obs: &ConsensusObs) -> String {
    if obs.name.trim().is_empty() {
        format!("@{}", obs.position)
    } else {
        obs.name.trim().to_uppercase()
    }
}

/// Reconcile per-source variant observations into one profile. `sources` is `(label, source_type,
/// observations)`; a source contributes at most one observation per variant (the last wins on dup).
pub fn reconcile(sources: &[(String, SourceType, Vec<ConsensusObs>)]) -> Vec<ConsensusVariant> {
    // Group all observations by variant key, preserving the source they came from + its weight.
    struct ObsRec {
        label: String,
        source_type: SourceType,
        state: ConsensusState,
        weight: f64,
    }
    struct Acc {
        repr: ConsensusObs,
        obs: Vec<ObsRec>,
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
            let weight = obs_weight(*source_type, o.depth, o.mapq, o.callable, o.region_modifier);
            acc.obs.push(ObsRec { label: label.clone(), source_type: *source_type, state: o.state, weight });
        }
    }

    let mut out: Vec<ConsensusVariant> = groups
        .into_values()
        .map(|acc| {
            let repr = acc.repr;
            // Weighted vote over non-NoCall observations (quality-weighted per observation).
            let (mut w_derived, mut w_ancestral) = (0.0f64, 0.0f64);
            let mut total = 0usize;
            for o in &acc.obs {
                match o.state {
                    ConsensusState::Derived => {
                        w_derived += o.weight;
                        total += 1;
                    }
                    ConsensusState::Ancestral => {
                        w_ancestral += o.weight;
                        total += 1;
                    }
                    ConsensusState::NoCall => {}
                }
            }

            let consensus = if total == 0 {
                ConsensusState::NoCall
            } else if w_derived >= w_ancestral {
                ConsensusState::Derived
            } else {
                ConsensusState::Ancestral
            };

            let support = acc.obs.iter().filter(|o| o.state == consensus && consensus != ConsensusState::NoCall).count();

            let total_weight = w_derived + w_ancestral;
            // Confidence in the consensus = winning weight / total weight (Scala ConfirmationThreshold).
            let confidence_score = if total_weight > 0.0 { w_derived.max(w_ancestral) / total_weight } else { 0.0 };
            let minority_fraction = 1.0 - confidence_score;

            let status = if total == 0 {
                ConsensusStatus::NoCoverage
            } else if minority_fraction > CONFLICT_FRACTION {
                ConsensusStatus::Conflict
            } else if consensus == ConsensusState::Derived && !repr.in_tree {
                // Derived off-tree call is novel/private — even from a single source (the common case).
                ConsensusStatus::Novel
            } else if total == 1 {
                ConsensusStatus::SingleSource
            } else if confidence_score >= CONFIRMATION_FRACTION {
                ConsensusStatus::Confirmed
            } else {
                ConsensusStatus::Pending
            };

            let sources = acc
                .obs
                .iter()
                .map(|o| SourceObs { label: o.label.clone(), source_type: o.source_type, state: o.state })
                .collect();

            ConsensusVariant {
                name: repr.name,
                position: repr.position,
                ancestral: repr.ancestral,
                derived: repr.derived,
                consensus,
                status,
                support,
                total,
                in_tree: repr.in_tree,
                confidence_score,
                sources,
            }
        })
        .collect();

    // Conflicts first (most actionable), then novel, then by name.
    out.sort_by(|a, b| status_rank(a.status).cmp(&status_rank(b.status)).then_with(|| a.name.cmp(&b.name)));
    out
}

fn status_rank(s: ConsensusStatus) -> u8 {
    match s {
        ConsensusStatus::Conflict => 0,
        ConsensusStatus::Novel => 1,
        ConsensusStatus::Pending => 2,
        ConsensusStatus::SingleSource => 3,
        ConsensusStatus::Confirmed => 4,
        ConsensusStatus::NoCoverage => 5,
    }
}

/// Per-status counts + overall confidence over a reconciled variant list.
pub fn summarize(variants: &[ConsensusVariant]) -> ConsensusSummary {
    let mut s = ConsensusSummary { total: variants.len(), ..Default::default() };
    for v in variants {
        match v.status {
            ConsensusStatus::Confirmed => s.confirmed += 1,
            ConsensusStatus::Novel => s.novel += 1,
            ConsensusStatus::Conflict => s.conflict += 1,
            ConsensusStatus::SingleSource => s.single_source += 1,
            // Pending / NoCoverage aren't headline counts; they fold into `total` only.
            ConsensusStatus::Pending | ConsensusStatus::NoCoverage => {}
        }
    }
    // Scala profile confidence: (confirmed + 0.7·novel − 0.5·conflict) / total, clamped [0,1].
    s.overall_confidence = if s.total == 0 {
        0.0
    } else {
        ((s.confirmed as f64 + 0.7 * s.novel as f64 - 0.5 * s.conflict as f64) / s.total as f64).clamp(0.0, 1.0)
    };
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    fn obs(name: &str, pos: i64, state: ConsensusState, in_tree: bool) -> ConsensusObs {
        ConsensusObs::snp(name, pos, "A", "G", state, in_tree)
    }

    #[test]
    fn obs_weight_applies_depth_mapq_callable() {
        // No quality data → bare source-type weight.
        assert!((obs_weight(SourceType::WgsShortRead, None, None, None, 1.0) - 0.85).abs() < 1e-9);
        // depth 100 → bonus min(√100/10,1)=1.0 → ×2; MQ 60 → ×1; callable → ×1.
        assert!((obs_weight(SourceType::WgsShortRead, Some(100), Some(60.0), Some(CallableState::Callable), 1.0) - 1.7).abs() < 1e-9);
        // Low coverage halves; a region modifier <1 scales down further.
        let w = obs_weight(SourceType::WgsShortRead, None, None, Some(CallableState::LowCoverage), 0.5);
        assert!((w - 0.85 * 0.5 * 0.5).abs() < 1e-9);
        // NoCoverage callability zeroes the weight.
        assert_eq!(obs_weight(SourceType::Sanger, Some(50), Some(60.0), Some(CallableState::NoCoverage), 1.0), 0.0);
    }

    #[test]
    fn confidence_score_and_overall() {
        let v = reconcile(&[
            ("a".into(), SourceType::WgsShortRead, vec![obs("M269", 1, ConsensusState::Derived, true)]),
            ("b".into(), SourceType::Chip, vec![obs("M269", 1, ConsensusState::Derived, true)]),
        ]);
        assert!((v[0].confidence_score - 1.0).abs() < 1e-9); // unanimous → full confidence
        let s = summarize(&v);
        assert!((s.overall_confidence - 1.0).abs() < 1e-9); // 1 confirmed / 1 total
    }

    #[test]
    fn two_sources_agree_in_tree_is_confirmed() {
        let v = reconcile(&[
            ("aln #1".into(), SourceType::WgsShortRead, vec![obs("M269", 100, ConsensusState::Derived, true)]),
            ("consumer".into(), SourceType::Chip, vec![obs("M269", 200, ConsensusState::Derived, true)]),
        ]);
        assert_eq!(v.len(), 1); // grouped by name across differing positions/builds
        assert_eq!(v[0].name, "M269");
        assert_eq!(v[0].consensus, ConsensusState::Derived);
        assert_eq!(v[0].status, ConsensusStatus::Confirmed);
        assert_eq!(v[0].support, 2);
        assert_eq!(v[0].total, 2);
    }

    #[test]
    fn derived_not_in_tree_is_novel() {
        let v = reconcile(&[
            ("aln #1".into(), SourceType::WgsShortRead, vec![obs("FT1", 100, ConsensusState::Derived, false)]),
            ("aln #2".into(), SourceType::WgsShortRead, vec![obs("FT1", 100, ConsensusState::Derived, false)]),
        ]);
        assert_eq!(v[0].status, ConsensusStatus::Novel);
    }

    #[test]
    fn comparable_weight_disagreement_is_conflict() {
        // WGS (0.85) derived vs Chip (0.5) ancestral → minority 0.5/1.35 ≈ 0.37 > 0.30 → conflict.
        let v = reconcile(&[
            ("aln #1".into(), SourceType::WgsShortRead, vec![obs("M269", 100, ConsensusState::Derived, true)]),
            ("consumer".into(), SourceType::Chip, vec![obs("M269", 100, ConsensusState::Ancestral, true)]),
        ]);
        assert_eq!(v[0].status, ConsensusStatus::Conflict);
        assert_eq!(v[0].consensus, ConsensusState::Derived); // higher weight wins the consensus
    }

    #[test]
    fn dominant_weight_disagreement_is_not_conflict() {
        // Sanger (1.0) derived vs Manual (0.3) ancestral → minority 0.3/1.3 ≈ 0.23 ≤ 0.30 → confirmed.
        let v = reconcile(&[
            ("sanger".into(), SourceType::Sanger, vec![obs("M269", 100, ConsensusState::Derived, true)]),
            ("manual".into(), SourceType::Manual, vec![obs("M269", 100, ConsensusState::Ancestral, true)]),
        ]);
        assert_eq!(v[0].consensus, ConsensusState::Derived);
        assert_eq!(v[0].status, ConsensusStatus::Confirmed);
        assert_eq!(v[0].support, 1); // only the Sanger source matches the derived consensus
    }

    #[test]
    fn single_source_is_single_source() {
        let v = reconcile(&[("aln #1".into(), SourceType::WgsShortRead, vec![obs("M269", 100, ConsensusState::Derived, true)])]);
        assert_eq!(v[0].status, ConsensusStatus::SingleSource);
    }

    #[test]
    fn nocall_excluded_from_vote() {
        let v = reconcile(&[
            ("aln #1".into(), SourceType::WgsShortRead, vec![obs("M269", 100, ConsensusState::Derived, true)]),
            ("aln #2".into(), SourceType::WgsShortRead, vec![obs("M269", 100, ConsensusState::NoCall, true)]),
        ]);
        assert_eq!(v[0].total, 1); // NoCall not counted
        assert_eq!(v[0].status, ConsensusStatus::SingleSource);
        assert_eq!(v[0].sources.len(), 2); // but still shown for provenance
    }

    #[test]
    fn summary_counts_by_status() {
        let v = reconcile(&[
            ("a".into(), SourceType::WgsShortRead, vec![obs("M269", 1, ConsensusState::Derived, true), obs("FT1", 2, ConsensusState::Derived, false)]),
            ("b".into(), SourceType::Chip, vec![obs("M269", 1, ConsensusState::Derived, true)]),
        ]);
        let s = summarize(&v);
        assert_eq!(s.total, 2);
        assert_eq!(s.confirmed, 1); // M269 (2 sources agree, in tree)
        assert_eq!(s.novel, 1); // FT1 (derived, not in tree → novel even single-source)
        assert_eq!(s.single_source, 0);
    }
}
