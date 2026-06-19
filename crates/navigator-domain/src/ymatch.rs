//! Cross-subject Y-chromosome matching — the *between-subjects* layer on top of the single-subject
//! Y profile. Given one subject, rank every other by Y relatedness (the FTDNA "Big Y match list"
//! idea): shared derived SNPs, shared private/novel variants, the divergence haplogroup, Y-STR
//! genetic distance, and rough SNP- and STR-based TMRCA estimates.
//!
//! This module is pure (no I/O): the app assembles a [`YMatchProfile`] per subject from cached data
//! — the consensus Y-variant set, the placement-tree lineage, and the imported STR markers — and
//! calls [`rank`]. SNP comparison is keyed by **variant name** (build-independent), matching the
//! consensus engine ([`crate::consensus`]); STR distance reuses [`crate::strprofile::values_match`]
//! so multi-copy markers compare order-independently.

use std::collections::HashSet;

use serde::{Deserialize, Serialize};

use crate::strprofile::{values_match, StrMarker};
use du_domain::ids::SampleGuid;

/// Big-Y-700 convention: ~1 SNP accumulates per this many years on the callable region (FTDNA cites
/// an average ≈ 83 yr/SNP). Used only for the **rough** SNP TMRCA — wide confidence interval.
pub const YEARS_PER_SNP: f64 = 83.0;
/// Years per generation, for converting a year estimate to generations.
pub const YEARS_PER_GEN: f64 = 32.0;
/// Average per-marker, per-generation Y-STR mutation rate (FTDNA-panel order of magnitude). Used only
/// for the **rough** STR TMRCA — wide confidence interval.
pub const MU_PER_MARKER_GEN: f64 = 0.0025;

/// Which evidence backed a pairwise comparison.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum YSignal {
    /// Both subjects have a placed Y-SNP profile AND comparable STR markers.
    SnpStr,
    /// Both have a placed Y-SNP profile (no comparable STR markers).
    Snp,
    /// STR markers only (one or both lack a placed Y-SNP profile).
    Str,
    /// Nothing comparable.
    None,
}

impl YSignal {
    /// Ranking tier — SNP-backed first, then STR-only, then nothing.
    fn tier(self) -> u8 {
        match self {
            YSignal::SnpStr | YSignal::Snp => 0,
            YSignal::Str => 1,
            YSignal::None => 2,
        }
    }
}

/// A rough time-to-most-recent-common-ancestor estimate. Both fields are **approximate** (wide CI).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Tmrca {
    pub generations: f64,
    pub years: f64,
}

/// A lightweight per-subject snapshot, assembled by the app from cached data and fed to [`compare_y`].
#[derive(Debug, Clone)]
pub struct YMatchProfile {
    pub guid: SampleGuid,
    pub donor: String,
    /// Placed terminal haplogroup name (if a Y profile exists).
    pub terminal: Option<String>,
    /// Root→terminal haplogroup lineage (empty when the subject has no placed Y profile).
    pub lineage: Vec<String>,
    /// Derived **tree** (in-tree) SNP names the subject carries.
    pub derived: HashSet<String>,
    /// Derived **off-tree / novel** (private) SNP names.
    pub novel: HashSet<String>,
    /// Imported Y-STR markers (first/preferred panel).
    pub str_markers: Vec<StrMarker>,
}

impl YMatchProfile {
    /// Whether the subject has Y-SNP calls to compare (independent of the tree/lineage being present —
    /// lineage only adds the divergence haplogroup).
    fn has_snp(&self) -> bool {
        !self.derived.is_empty() || !self.novel.is_empty()
    }
}

/// One ranked Y match against a query subject.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct YMatch {
    pub guid: SampleGuid,
    pub donor: String,
    pub terminal: Option<String>,
    /// Count of derived **tree** SNPs both carry.
    pub shared_derived: usize,
    /// Count of **private/novel** SNPs both carry (shared off-tree variants = candidate sub-branch).
    pub shared_novel: usize,
    /// The deepest haplogroup the two lineages share (their LCA), if both are placed.
    pub divergence: Option<String>,
    /// Y-STR genetic distance over markers present in both (None when not comparable).
    pub str_gd: Option<i64>,
    /// Number of STR markers compared.
    pub str_markers: i64,
    /// Rough SNP-based TMRCA (only when both have a placed Y profile).
    pub snp_tmrca: Option<Tmrca>,
    /// Rough STR-based TMRCA (only when STR markers were comparable).
    pub str_tmrca: Option<Tmrca>,
    pub signal: YSignal,
}

/// The deepest haplogroup two lineages share — the longest common prefix of the two root→terminal
/// paths — and its depth (number of shared steps). Returns `(None, 0)` if either lineage is empty.
fn divergence(a: &[String], b: &[String]) -> (Option<String>, usize) {
    let mut last = None;
    let mut depth = 0;
    for (x, y) in a.iter().zip(b.iter()) {
        if x == y {
            last = Some(x.clone());
            depth += 1;
        } else {
            break;
        }
    }
    (last, depth)
}

/// Y-STR genetic distance over markers present in both profiles: count markers whose values differ
/// (multi-copy values compared order-independently). Returns `(differing, compared)`.
fn str_gd(a: &[StrMarker], b: &[StrMarker]) -> (i64, i64) {
    let mut differing = 0;
    let mut compared = 0;
    for ma in a {
        if let Some(mb) = b.iter().find(|m| m.marker.eq_ignore_ascii_case(&ma.marker)) {
            compared += 1;
            if !values_match(&ma.value, &mb.value) {
                differing += 1;
            }
        }
    }
    (differing, compared)
}

/// Rough SNP TMRCA: each lineage accumulates its private (non-shared) variants since divergence at
/// ~[`YEARS_PER_SNP`]; TMRCA years ≈ average private count × yr/SNP. Approximate — depends on equal
/// callable coverage between the two subjects.
fn snp_tmrca(private_a: usize, private_b: usize) -> Tmrca {
    let years = ((private_a + private_b) as f64 / 2.0) * YEARS_PER_SNP;
    Tmrca {
        generations: years / YEARS_PER_GEN,
        years,
    }
}

/// Rough STR TMRCA via a stepwise model: expected differences over two lineages ≈ 2·markers·μ·g, so
/// generations to MRCA ≈ gd / (2·markers·μ). Approximate — single average mutation rate, no per-marker
/// rates or TiP-grade modelling.
fn str_tmrca(gd: i64, markers: i64) -> Option<Tmrca> {
    if markers <= 0 {
        return None;
    }
    let generations = gd as f64 / (2.0 * markers as f64 * MU_PER_MARKER_GEN);
    Some(Tmrca {
        generations,
        years: generations * YEARS_PER_GEN,
    })
}

/// Compare a query subject against one candidate.
pub fn compare_y(query: &YMatchProfile, cand: &YMatchProfile) -> YMatch {
    let snp_possible = query.has_snp() && cand.has_snp();

    let (shared_derived, shared_novel, divergence_hg, snp_tmrca_est) = if snp_possible {
        let shared_derived = query.derived.intersection(&cand.derived).count();
        let shared_novel = query.novel.intersection(&cand.novel).count();
        let (div, _depth) = divergence(&query.lineage, &cand.lineage);
        // Private = variants one side carries that the other does not (tree-derived + novel). A rough
        // proxy for the branch each accumulated since their divergence.
        let union_q = query.derived.len() + query.novel.len();
        let union_c = cand.derived.len() + cand.novel.len();
        let shared_total = shared_derived + shared_novel;
        let private_a = union_q.saturating_sub(shared_total);
        let private_b = union_c.saturating_sub(shared_total);
        (shared_derived, shared_novel, div, Some(snp_tmrca(private_a, private_b)))
    } else {
        (0, 0, None, None)
    };

    let (gd, compared) = str_gd(&query.str_markers, &cand.str_markers);
    let str_possible = compared > 0;
    let (str_gd_val, str_markers, str_tmrca_est) = if str_possible {
        (Some(gd), compared, str_tmrca(gd, compared))
    } else {
        (None, 0, None)
    };

    let signal = match (snp_possible, str_possible) {
        (true, true) => YSignal::SnpStr,
        (true, false) => YSignal::Snp,
        (false, true) => YSignal::Str,
        (false, false) => YSignal::None,
    };

    YMatch {
        guid: cand.guid,
        donor: cand.donor.clone(),
        terminal: cand.terminal.clone(),
        shared_derived,
        shared_novel,
        divergence: divergence_hg,
        str_gd: str_gd_val,
        str_markers,
        snp_tmrca: snp_tmrca_est,
        str_tmrca: str_tmrca_est,
        signal,
    }
}

/// Rank candidates against the query, best match first. SNP-primary: SNP-backed matches first (more
/// shared derived SNPs, then deeper divergence), then STR-only by ascending genetic distance.
/// Candidates with no comparable evidence (`YSignal::None`) are dropped. The query itself is skipped.
pub fn rank(query: &YMatchProfile, candidates: &[YMatchProfile]) -> Vec<YMatch> {
    let mut out: Vec<YMatch> = candidates
        .iter()
        .filter(|c| c.guid != query.guid)
        .map(|c| compare_y(query, c))
        .filter(|m| m.signal != YSignal::None)
        .collect();
    out.sort_by(|a, b| {
        // tier asc; then within tier: shared_derived desc, str_gd asc, donor asc.
        a.signal
            .tier()
            .cmp(&b.signal.tier())
            .then(b.shared_derived.cmp(&a.shared_derived))
            .then(a.str_gd.unwrap_or(i64::MAX).cmp(&b.str_gd.unwrap_or(i64::MAX)))
            .then(a.donor.cmp(&b.donor))
    });
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn marker(m: &str, v: &str) -> StrMarker {
        StrMarker {
            marker: m.into(),
            value: v.into(),
        }
    }

    /// Deterministic distinct guid per donor name (so self-skip works without a real UUID source).
    fn guid_for(donor: &str) -> SampleGuid {
        let mut h: u128 = 0xcbf2_9ce4_8422_2325;
        for b in donor.bytes() {
            h = h.wrapping_mul(0x0100_0000_01b3).wrapping_add(b as u128);
        }
        SampleGuid(uuid::Uuid::from_u128(h))
    }

    fn prof(donor: &str, lineage: &[&str], derived: &[&str], novel: &[&str], strs: &[(&str, &str)]) -> YMatchProfile {
        YMatchProfile {
            guid: guid_for(donor),
            donor: donor.into(),
            terminal: lineage.last().map(|s| s.to_string()),
            lineage: lineage.iter().map(|s| s.to_string()).collect(),
            derived: derived.iter().map(|s| s.to_string()).collect(),
            novel: novel.iter().map(|s| s.to_string()).collect(),
            str_markers: strs.iter().map(|(m, v)| marker(m, v)).collect(),
        }
    }

    #[test]
    fn divergence_is_longest_common_prefix() {
        let a = vec![
            "R".to_string(),
            "R-M269".to_string(),
            "R-L21".to_string(),
            "R-CTS4466".to_string(),
        ];
        let b = vec![
            "R".to_string(),
            "R-M269".to_string(),
            "R-L21".to_string(),
            "R-DF13".to_string(),
        ];
        let (hg, depth) = divergence(&a, &b);
        assert_eq!(hg.as_deref(), Some("R-L21"));
        assert_eq!(depth, 3);
        assert_eq!(divergence(&a, &[]).0, None);
    }

    #[test]
    fn str_gd_is_multicopy_order_independent() {
        let a = [
            marker("DYS393", "13"),
            marker("DYS385", "16-17"),
            marker("DYS390", "24"),
        ];
        let b = [
            marker("DYS393", "13"),
            marker("DYS385", "17-16"),
            marker("DYS390", "25"),
        ];
        // DYS385 16-17 == 17-16 (no diff); DYS390 differs; DYS393 same → gd 1 over 3.
        assert_eq!(str_gd(&a, &b), (1, 3));
    }

    #[test]
    fn identical_subjects_diverge_at_terminal_with_zero_gd() {
        let q = prof(
            "Q",
            &["R", "R-M269", "R-CTS4466"],
            &["M269", "CTS4466"],
            &[],
            &[("DYS393", "13")],
        );
        let c = prof(
            "C",
            &["R", "R-M269", "R-CTS4466"],
            &["M269", "CTS4466"],
            &[],
            &[("DYS393", "13")],
        );
        let m = compare_y(&q, &c);
        assert_eq!(m.signal, YSignal::SnpStr);
        assert_eq!(m.divergence.as_deref(), Some("R-CTS4466"));
        assert_eq!(m.shared_derived, 2);
        assert_eq!(m.str_gd, Some(0));
    }

    #[test]
    fn snp_backed_outranks_str_only_and_more_shared_wins() {
        let q = prof(
            "Q",
            &["R", "R-M269", "R-CTS4466"],
            &["M269", "L21", "CTS4466"],
            &[],
            &[("DYS393", "13")],
        );
        // Close SNP match (shares all 3 derived).
        let close = prof(
            "Close",
            &["R", "R-M269", "R-CTS4466"],
            &["M269", "L21", "CTS4466"],
            &[],
            &[("DYS393", "13")],
        );
        // Distant SNP match (shares only the backbone).
        let distant = prof("Distant", &["R", "R-M269"], &["M269"], &[], &[("DYS393", "14")]);
        // STR-only (no lineage) — must rank below any SNP-backed match.
        let stronly = prof("StrOnly", &[], &[], &[], &[("DYS393", "13")]);
        let ranked = rank(&q, &[distant.clone(), stronly.clone(), close.clone()]);
        assert_eq!(ranked.len(), 3);
        assert_eq!(ranked[0].donor, "Close");
        assert_eq!(ranked[1].donor, "Distant");
        assert_eq!(ranked[2].donor, "StrOnly");
        assert_eq!(ranked[2].signal, YSignal::Str);
    }

    #[test]
    fn no_common_evidence_is_dropped_and_self_skipped() {
        let q = prof("Q", &["R", "R-M269"], &["M269"], &[], &[("DYS393", "13")]);
        // No lineage and no overlapping STR markers → nothing comparable.
        let nothing = prof("Nothing", &[], &[], &[], &[("DYS999", "10")]);
        let ranked = rank(&q, &[q.clone(), nothing]);
        assert!(ranked.is_empty());
    }

    #[test]
    fn snp_tmrca_grows_with_private_variants() {
        let near = snp_tmrca(1, 1);
        let far = snp_tmrca(10, 12);
        assert!(far.years > near.years);
        assert_eq!(near.years, YEARS_PER_SNP); // (1+1)/2 * 83
    }
}
