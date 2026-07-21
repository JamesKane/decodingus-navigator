//! Donor-level reconciliation of Y/mtDNA haplogroup calls across multiple sources (runs,
//! platforms, Sanger). Phase 1–2 of `documents/design/MultiSource_Reconciliation.md`:
//! per-source [`RunHaplogroupCall`]s combine into a [`Consensus`] by tree topology.
//!
//! Pure types + the consensus algorithm; persistence and the per-source recording live in
//! the app/store. Per-variant concordance (all DNA types) lives in [`crate::consensus`];
//! identity verification and heteroplasmy are later phases.

use serde::{Deserialize, Serialize};

/// Which uniparental lineage a call describes.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum DnaType {
    Y,
    Mt,
}

impl DnaType {
    pub fn as_str(self) -> &'static str {
        match self {
            DnaType::Y => "Y",
            DnaType::Mt => "Mt",
        }
    }
}

/// Where a per-source haplogroup call came from — the precedence tier used when reconciling.
///
/// - `External`: a trusted external caller (a GATK4 GVCF / 1240K call set imported via the sidecar
///   fast path). The user runs an established pipeline and wants these preferred.
/// - `NavigatorWalk`: Navigator's own genotyping — the CRAM walk, and chip/vendor-VCF placement.
/// - `Manual`: a user override (persisted separately today; included for a complete precedence order).
///
/// Higher [`rank`](CallProvenance::rank) wins when the "prefer external caller" policy is on.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum CallProvenance {
    NavigatorWalk,
    External,
    Manual,
}

impl CallProvenance {
    pub fn as_str(self) -> &'static str {
        match self {
            CallProvenance::NavigatorWalk => "navigator-walk",
            CallProvenance::External => "external",
            CallProvenance::Manual => "manual",
        }
    }

    /// Parse a stored provenance token; anything unrecognized — including legacy `NULL`/empty rows
    /// written before the column existed — is the internal `NavigatorWalk` tier.
    pub fn from_token(s: &str) -> CallProvenance {
        match s.trim() {
            "external" => CallProvenance::External,
            "manual" => CallProvenance::Manual,
            _ => CallProvenance::NavigatorWalk,
        }
    }

    /// Precedence tier — higher wins under the prefer-external policy.
    pub fn rank(self) -> u8 {
        match self {
            CallProvenance::NavigatorWalk => 0,
            CallProvenance::External => 1,
            CallProvenance::Manual => 2,
        }
    }
}

/// A haplogroup call from one source (a sequencing run, chip, STR panel, or Sanger entry).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct RunHaplogroupCall {
    /// Human label for the source (e.g. "aln #5 bwa-mem2").
    pub source_label: String,
    /// Terminal haplogroup name.
    pub haplogroup: String,
    /// Root→terminal lineage of haplogroup names.
    pub lineage: Vec<String>,
    /// Assignment score (Kulczynski) — confidence proxy.
    pub score: f64,
    pub matched: i64,
    pub expected: i64,
}

/// How compatible a set of calls is (Scala `CompatibilityLevel`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CompatibilityLevel {
    /// Same branch, differing depths — all calls lie on one root→tip path.
    Compatible,
    /// Diverge near the tips.
    MinorDivergence,
    /// Diverge on the backbone.
    MajorDivergence,
    /// Diverge near the root — likely different individuals.
    Incompatible,
}

/// The reconciled donor-level result for one DNA type.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Consensus {
    pub haplogroup: String,
    pub lineage: Vec<String>,
    pub compatibility: CompatibilityLevel,
    /// The deepest node all sources agree on, when they diverge.
    pub divergence_point: Option<String>,
    pub confidence: f64,
    pub run_count: usize,
    /// True when a user manual override replaced the computed consensus.
    pub overridden: bool,
    pub warnings: Vec<String>,
}

/// An entry in the reconciliation audit log (Scala `AuditEntry`).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AuditEntry {
    /// RFC3339 timestamp.
    pub timestamp: String,
    /// INITIAL / RUN_RECORDED / MANUAL_OVERRIDE / OVERRIDE_CLEARED / RECOMPUTED.
    pub action: String,
    pub note: String,
}

/// Whether multiple sources come from the same individual (Scala `IdentityVerification`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum VerificationStatus {
    VerifiedSame,
    LikelySame,
    Uncertain,
    LikelyDifferent,
    VerifiedDifferent,
}

/// Identity evidence between two sources: autosomal genotype concordance (the primary
/// signal — same individual ≈ 1.0, relatives notably lower) plus Y-STR corroboration.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct IdentityVerification {
    pub status: VerificationStatus,
    pub method: String,
    /// Fraction of shared-called sites with identical dosage (0–1), if any compared.
    pub snp_concordance: Option<f64>,
    pub sites_compared: i64,
    /// Differing Y-STR markers across shared markers, if STR profiles were available.
    pub y_str_distance: Option<i64>,
    pub y_str_markers: i64,
}

/// Classify identity from genotype concordance (primary) + optional Y-STR distance.
/// Concordance thresholds are heuristic: a same-individual pair sits near 1.0 while even
/// parent-child concordance (exact-dosage match) is well below; Y-STR distance corroborates
/// (0 over many markers supports a shared paternal line, a mismatch flags a conflict).
pub fn classify_identity(
    concordance: Option<f64>,
    sites_compared: i64,
    y_str_distance: Option<i64>,
    y_str_markers: i64,
) -> IdentityVerification {
    let (status, method) = match concordance {
        Some(c) if sites_compared > 0 => {
            let s = if c >= 0.95 {
                VerificationStatus::VerifiedSame
            } else if c >= 0.85 {
                VerificationStatus::LikelySame
            } else if c >= 0.65 {
                VerificationStatus::Uncertain
            } else if c >= 0.50 {
                VerificationStatus::LikelyDifferent
            } else {
                VerificationStatus::VerifiedDifferent
            };
            let m = if y_str_markers > 0 {
                "SNP concordance + Y-STR"
            } else {
                "SNP concordance"
            };
            (s, m.to_string())
        }
        // No shared genotypes: Y-STR alone is paternal-line only — never "verified".
        _ if y_str_markers > 0 => {
            let s = match y_str_distance {
                Some(0) => VerificationStatus::LikelySame,
                Some(_) => VerificationStatus::Uncertain,
                None => VerificationStatus::Uncertain,
            };
            (s, "Y-STR only".to_string())
        }
        _ => (VerificationStatus::Uncertain, "no shared data".to_string()),
    };
    IdentityVerification {
        status,
        method,
        snp_concordance: concordance,
        sites_compared,
        y_str_distance,
        y_str_markers,
    }
}

fn is_prefix(short: &[String], long: &[String]) -> bool {
    short.len() <= long.len() && short.iter().zip(long).all(|(a, b)| a == b)
}

/// Longest node-name prefix shared by every lineage.
fn common_prefix(calls: &[RunHaplogroupCall]) -> Vec<String> {
    let Some(first) = calls.first() else { return Vec::new() };
    let mut len = first.lineage.len();
    for c in &calls[1..] {
        len = len.min(c.lineage.len());
        while len > 0 && first.lineage[..len] != c.lineage[..len] {
            len -= 1;
        }
    }
    first.lineage[..len].to_vec()
}

/// Reconcile per-source calls into a donor-level consensus.
///
/// When all calls lie on one root→tip path (compatible), the consensus is the **most
/// confident** call — not blindly the deepest, since a low-coverage source may extend one
/// node further on thin evidence; any strictly-deeper call is reported as a tentative
/// warning. When calls diverge, the consensus is the deepest node they all agree on (the
/// LCA), and the divergence depth sets the compatibility level.
pub fn reconcile(calls: &[RunHaplogroupCall]) -> Option<Consensus> {
    if calls.is_empty() {
        return None;
    }
    let run_count = calls.len();
    let longest = calls.iter().max_by_key(|c| c.lineage.len()).unwrap();
    let all_on_path = calls.iter().all(|c| is_prefix(&c.lineage, &longest.lineage));

    if all_on_path {
        // Most-confident call wins; flag any deeper-but-not-most-confident extension.
        let best = calls
            .iter()
            .max_by(|a, b| a.score.partial_cmp(&b.score).unwrap_or(std::cmp::Ordering::Equal))
            .unwrap();
        let mut warnings = Vec::new();
        for c in calls {
            if c.lineage.len() > best.lineage.len() {
                warnings.push(format!(
                    "{} places deeper at {} (score {:.3}) — tentative",
                    c.source_label, c.haplogroup, c.score
                ));
            }
        }
        return Some(Consensus {
            haplogroup: best.haplogroup.clone(),
            lineage: best.lineage.clone(),
            compatibility: CompatibilityLevel::Compatible,
            divergence_point: None,
            confidence: best.score,
            run_count,
            overridden: false,
            warnings,
        });
    }

    // Divergent: consensus is the LCA; level by how deep the LCA sits.
    let prefix = common_prefix(calls);
    let max_depth = longest.lineage.len().max(1);
    let ratio = prefix.len() as f64 / max_depth as f64;
    // Sharing only the root (≤1 node) means different lineages entirely. Otherwise the
    // LCA's relative depth distinguishes a tip split from a backbone split.
    let compatibility = if prefix.len() <= 1 {
        CompatibilityLevel::Incompatible
    } else if ratio >= 0.66 {
        CompatibilityLevel::MinorDivergence
    } else {
        CompatibilityLevel::MajorDivergence
    };
    let divergence_point = prefix.last().cloned();
    let haplogroup = divergence_point.clone().unwrap_or_else(|| "root".to_string());

    // Distinct terminal calls, for the warning.
    let mut terminals: Vec<&str> = calls.iter().map(|c| c.haplogroup.as_str()).collect();
    terminals.sort_unstable();
    terminals.dedup();
    let warnings = vec![format!(
        "sources diverge below {}: {}",
        divergence_point.as_deref().unwrap_or("root"),
        terminals.join(", ")
    )];
    let confidence = calls.iter().map(|c| c.score).sum::<f64>() / run_count as f64;

    Some(Consensus {
        haplogroup,
        lineage: prefix,
        compatibility,
        divergence_point,
        confidence,
        run_count,
        overridden: false,
        warnings,
    })
}

/// Reconcile per-source calls, honoring call provenance.
///
/// When `prefer_external` is set, only the highest-precedence tier present
/// (Manual > External > NavigatorWalk) is reconciled into the consensus; lower tiers that place a
/// different terminal are surfaced as a warning rather than allowed to drag the call. This is what
/// stops a damaged ancient-DNA CRAM walk from out-scoring — and silently replacing — a clean
/// external GATK4/1240K placement. When `prefer_external` is false, every call is reconciled
/// together by confidence (the source-blind [`reconcile`]), so the two behave identically whenever
/// no external call is present.
pub fn reconcile_with_provenance(
    calls: &[(CallProvenance, RunHaplogroupCall)],
    prefer_external: bool,
) -> Option<Consensus> {
    if calls.is_empty() {
        return None;
    }
    if !prefer_external {
        let flat: Vec<RunHaplogroupCall> = calls.iter().map(|(_, c)| c.clone()).collect();
        return reconcile(&flat);
    }
    let top = calls.iter().map(|(p, _)| p.rank()).max().unwrap();
    let top_calls: Vec<RunHaplogroupCall> = calls
        .iter()
        .filter(|(p, _)| p.rank() == top)
        .map(|(_, c)| c.clone())
        .collect();
    let mut consensus = reconcile(&top_calls)?;
    // Note any lower-precedence source that places a *different* terminal — informative, not
    // authoritative (the external caller was preferred).
    let mut lower: Vec<&str> = calls
        .iter()
        .filter(|(p, _)| p.rank() < top)
        .map(|(_, c)| c.haplogroup.as_str())
        .filter(|h| *h != consensus.haplogroup)
        .collect();
    lower.sort_unstable();
    lower.dedup();
    if !lower.is_empty() {
        consensus
            .warnings
            .push(format!("lower-precedence sources place elsewhere: {} (external caller preferred)", lower.join(", ")));
    }
    Some(consensus)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn call(label: &str, score: f64, lineage: &[&str]) -> RunHaplogroupCall {
        RunHaplogroupCall {
            source_label: label.into(),
            haplogroup: (*lineage.last().unwrap()).into(),
            lineage: lineage.iter().map(|s| s.to_string()).collect(),
            score,
            matched: 0,
            expected: 0,
        }
    }

    #[test]
    fn single_call_is_itself() {
        let c = reconcile(&[call("a", 0.7, &["root", "R", "R-M269"])]).unwrap();
        assert_eq!(c.haplogroup, "R-M269");
        assert_eq!(c.compatibility, CompatibilityLevel::Compatible);
        assert_eq!(c.run_count, 1);
    }

    #[test]
    fn compatible_prefers_the_confident_call_not_the_deepest() {
        // High-confidence short-read at FGC29067; low-confidence HiFi one node deeper.
        let wgs = call("wgs", 0.750, &["root", "R", "R-FGC29067"]);
        let hifi = call("hifi", 0.537, &["root", "R", "R-FGC29067", "R-FGC29071"]);
        let c = reconcile(&[wgs, hifi]).unwrap();
        assert_eq!(c.compatibility, CompatibilityLevel::Compatible);
        assert_eq!(c.haplogroup, "R-FGC29067"); // confident node, not the tentative deeper one
        assert!((c.confidence - 0.750).abs() < 1e-9);
        assert_eq!(c.warnings.len(), 1); // notes R-FGC29071 as tentative
        assert!(c.warnings[0].contains("R-FGC29071"));
    }

    #[test]
    fn compatible_takes_deeper_when_it_is_the_confident_one() {
        let shallow = call("a", 0.5, &["root", "R", "R-M269"]);
        let deep = call("b", 0.9, &["root", "R", "R-M269", "R-L21", "R-DF13"]);
        let c = reconcile(&[shallow, deep]).unwrap();
        assert_eq!(c.haplogroup, "R-DF13");
        assert!(c.warnings.is_empty());
    }

    #[test]
    fn tip_divergence_is_minor() {
        // Agree to R-L21 (depth 3 of 4), then split R-A vs R-B.
        let a = call("a", 0.8, &["root", "R", "R-L21", "R-A"]);
        let b = call("b", 0.8, &["root", "R", "R-L21", "R-B"]);
        let c = reconcile(&[a, b]).unwrap();
        assert_eq!(c.compatibility, CompatibilityLevel::MinorDivergence);
        assert_eq!(c.divergence_point.as_deref(), Some("R-L21"));
        assert_eq!(c.haplogroup, "R-L21");
    }

    #[test]
    fn identity_from_concordance() {
        assert_eq!(
            classify_identity(Some(0.995), 5000, Some(0), 37).status,
            VerificationStatus::VerifiedSame
        );
        assert_eq!(
            classify_identity(Some(0.88), 5000, None, 0).status,
            VerificationStatus::LikelySame
        );
        assert_eq!(
            classify_identity(Some(0.45), 5000, None, 0).status,
            VerificationStatus::VerifiedDifferent
        );
        // no autosomal data -> Y-STR only, never "verified"
        assert_eq!(
            classify_identity(None, 0, Some(0), 111).status,
            VerificationStatus::LikelySame
        );
        assert_eq!(
            classify_identity(None, 0, None, 0).status,
            VerificationStatus::Uncertain
        );
    }

    #[test]
    fn root_divergence_is_incompatible() {
        // Share only "root" — different haplogroups entirely (different individuals?).
        let a = call("a", 0.8, &["root", "R", "R-M269", "R-L21"]);
        let b = call("b", 0.8, &["root", "J", "J-M267"]);
        let c = reconcile(&[a, b]).unwrap();
        assert_eq!(c.compatibility, CompatibilityLevel::Incompatible);
    }

    #[test]
    fn prefer_external_wins_over_a_higher_scoring_walk() {
        // The ancient-DNA case: the CRAM walk out-scores the clean external call (deamination
        // inflates matched sites onto a *different* deep terminal), yet the external call must win.
        let external = call("gatk4 gvcf", 0.60, &["root", "R", "R-M269", "R-L21"]);
        let walk = call("cram walk", 0.95, &["root", "R", "R-M269", "R-L2"]);
        let c = reconcile_with_provenance(
            &[(CallProvenance::External, external), (CallProvenance::NavigatorWalk, walk)],
            true,
        )
        .unwrap();
        assert_eq!(c.haplogroup, "R-L21"); // external terminal, not the higher-scoring walk's R-L2
        assert!(c.warnings.iter().any(|w| w.contains("R-L2")));
    }

    #[test]
    fn without_prefer_external_score_still_wins() {
        // Toggle off → source-blind, identical to plain reconcile: highest score wins.
        let external = call("gatk4 gvcf", 0.60, &["root", "R", "R-M269"]);
        let walk = call("cram walk", 0.95, &["root", "R", "R-M269", "R-L21", "R-DF13"]);
        let c = reconcile_with_provenance(
            &[(CallProvenance::External, external), (CallProvenance::NavigatorWalk, walk)],
            false,
        )
        .unwrap();
        assert_eq!(c.haplogroup, "R-DF13");
    }

    #[test]
    fn prefer_external_with_no_external_is_plain_reconcile() {
        // Two walk calls, prefer_external on → same tier, behaves exactly like reconcile.
        let a = call("a", 0.5, &["root", "R", "R-M269"]);
        let b = call("b", 0.9, &["root", "R", "R-M269", "R-L21"]);
        let c = reconcile_with_provenance(
            &[(CallProvenance::NavigatorWalk, a), (CallProvenance::NavigatorWalk, b)],
            true,
        )
        .unwrap();
        assert_eq!(c.haplogroup, "R-L21");
        assert!(c.warnings.is_empty());
    }
}
