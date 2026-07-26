//! Archaic (Neanderthal / Denisovan) ancestry — panel types.
//!
//! Design: `documents/design/ArchaicAncestry_Design.md`. This module carries the **Tier A** asset
//! types (the archaic-informative marker panel and the percentile reference); the counting routine
//! and the Tier B segment HMM land on top of them in later milestones.
//!
//! The panel is *computed by us* rather than ingested (design §3a): the Sankararaman 2014 list
//! 23andMe used is not publicly downloadable, and every openly-licensed alternative is either the
//! wrong build or per-individual probabilistic. We therefore derive sites from the EVA archaic
//! VCFs, Ensembl-75 EPO ancestral alleles, and a 1kGP AFR outgroup, redistributing only the
//! derived sites.
//!
//! A consequence worth remembering when validating: **the panel is ours, so its marker count is
//! panel-relative and is not comparable to a vendor's count by equality** — compare the per-site
//! rate on the intersection instead (design §10, M2 validation gate).

use serde::{Deserialize, Serialize};

use crate::error::AnalysisError;

/// The four openly-downloadable archaic reference genomes, in the fixed order every
/// [`ArchaicSite::calls`] array uses.
pub const ARCHAIC_GENOMES: [&str; 4] = ["AltaiNeanderthal", "Vindija33.19", "Chagyrskaya8", "Denisova3"];

/// Index into [`ARCHAIC_GENOMES`] / [`ArchaicSite::calls`] for the sole Denisovan genome. The other
/// three are Neanderthals, which is what [`classify_diagnostic`] keys off.
pub const DENISOVA: usize = 3;

/// One archaic genome's state at a site, expressed **relative to the site's derived allele** rather
/// than to ref/alt — so it survives the ref/alt swap that CHM13 orientation can apply.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ArchaicCall {
    /// Homozygous for the ancestral allele.
    HomAncestral,
    /// One derived copy.
    Het,
    /// Homozygous derived — the introgression donor state (design §4, step 2).
    HomDerived,
    /// Missing / filtered out by the genome's own quality mask.
    NoCall,
}

impl ArchaicCall {
    /// Whether this genome carries the derived allele at all (het or hom).
    pub fn carries_derived(self) -> bool {
        matches!(self, ArchaicCall::Het | ArchaicCall::HomDerived)
    }
}

/// Which archaic lineage a site's derived allele points to.
///
/// The HMM in Tier B cannot itself separate Neanderthal from Denisovan (they coalesce before either
/// meets modern humans, design §3); this classification is what lets called segments be labelled
/// downstream, so it is stored per site at build time.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum DiagnosticClass {
    /// Derived in ≥1 Neanderthal, absent from Denisova.
    Neanderthal,
    /// Derived in Denisova, absent from every Neanderthal.
    Denisovan,
    /// Derived in both lineages — informative for "archaic" but not for attribution.
    SharedArchaic,
}

/// Classify a site from its per-genome calls. `NoCall` is treated as "no evidence of derived",
/// which is deliberately conservative: a masked-out Denisova makes a site Neanderthal-diagnostic
/// only in the sense that we cannot show it is shared, so attribution downstream should lean on
/// [`DiagnosticClass::Neanderthal`] counts in aggregate rather than trusting any single site.
pub fn classify_diagnostic(calls: &[ArchaicCall; 4]) -> DiagnosticClass {
    let nea = calls
        .iter()
        .enumerate()
        .any(|(i, c)| i != DENISOVA && c.carries_derived());
    let den = calls[DENISOVA].carries_derived();
    match (nea, den) {
        (true, false) => DiagnosticClass::Neanderthal,
        (false, true) => DiagnosticClass::Denisovan,
        _ => DiagnosticClass::SharedArchaic,
    }
}

/// One archaic-informative marker.
///
/// Coordinates are CHM13 and **oriented**: `reference_allele` is the actual CHM13 base at
/// `position`. That orientation is not optional — the pipeline lifts GRCh37→CHM13 with `CrossMap
/// bed`, which is not allele-aware, so roughly a third of sites arrive with ref/alt reversed
/// relative to CHM13. Shipping them unoriented is the bug that cost the ancient-ancestry build a
/// full rebuild (that design's §7.16).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ArchaicSite {
    pub contig: String,
    /// 1-based, on the panel's `build`.
    pub position: i64,
    /// The actual base on the panel's build (post-orientation).
    pub reference_allele: char,
    pub alternate_allele: char,
    /// The archaic-derived allele — always one of `reference_allele` / `alternate_allele`. Stored as
    /// a base, not a ref/alt flag, so a later orientation pass cannot silently invert its meaning.
    pub archaic_derived_allele: char,
    /// Per-genome state, indexed by [`ARCHAIC_GENOMES`].
    pub calls: [ArchaicCall; 4],
    pub diagnostic_class: DiagnosticClass,
    /// Derived-allele frequency in the African outgroup, kept for transparency and so the panel can
    /// be re-filtered at a stricter threshold without a rebuild.
    pub afr_freq: f32,
}

impl ArchaicSite {
    /// Copies of the archaic-derived allele (0–2) in a diploid observation.
    ///
    /// Takes the subject's two alleles as bases rather than a dosage, because dosage is defined
    /// against ref/alt whereas the derived allele may be either one.
    pub fn derived_copies(&self, allele_a: char, allele_b: char) -> u8 {
        let d = self.archaic_derived_allele.to_ascii_uppercase();
        u8::from(allele_a.to_ascii_uppercase() == d) + u8::from(allele_b.to_ascii_uppercase() == d)
    }
}

/// The site-selection thresholds a panel was built with, carried in the asset itself.
///
/// Recorded because these are the panel's scientific content: the site list *is* the product, every
/// downstream number inherits it, and a count is meaningless without knowing which filter produced
/// it. Design §10 fixes them by calibration (M1 checkpoint A), not by the illustrative values in §4.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ArchaicPanelThresholds {
    /// Maximum derived-allele frequency in the African outgroup (the introgression signature: rare
    /// or absent in Africans).
    pub max_afr_freq: f32,
    /// Minimum derived-allele frequency outside Africa.
    pub min_non_afr_freq: f32,
}

/// The Tier-A archaic marker panel (`archaic_markers_<build>.bin`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ArchaicMarkerPanel {
    /// Canonical reference build the site coordinates are in (e.g. "chm13v2.0").
    pub build: String,
    pub thresholds: ArchaicPanelThresholds,
    pub sites: Vec<ArchaicSite>,
}

impl ArchaicMarkerPanel {
    /// Deserialize from the built binary (bincode).
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, AnalysisError> {
        bincode::deserialize(bytes).map_err(|e| AnalysisError::Message(format!("archaic panel decode: {e}")))
    }

    /// Serialize to the binary form the builder writes.
    pub fn to_bytes(&self) -> Result<Vec<u8>, AnalysisError> {
        bincode::serialize(self).map_err(|e| AnalysisError::Message(format!("archaic panel encode: {e}")))
    }

    pub fn len(&self) -> usize {
        self.sites.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sites.is_empty()
    }

    /// The denominator a Tier-A count is reported against: two copies per site.
    ///
    /// 23andMe reports exactly this shape — the ground-truth report for the project's validation
    /// sample reads "191 of 7,462", and 7,462 is 2 × 3,731 assayed sites (design §10).
    pub fn possible_copies(&self) -> usize {
        self.sites.len() * 2
    }
}

/// Tier-A count distribution for one reference population (`archaic_marker_dist_<build>.bin`).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct CohortCounts {
    /// Fine population code (e.g. "GBR").
    pub population: String,
    /// Super-population the fine code rolls up to (e.g. "EUR").
    pub super_population: String,
    /// Per-sample archaic-derived copy totals, ascending.
    pub counts: Vec<u32>,
}

/// The percentile reference asset.
///
/// Stored **per population** rather than pre-reduced to super-populations. v1 renders the percentile
/// against a super-population cohort (design §9 Q3 — keying it to the user's inferred fine ancestry
/// would let an ancestry error silently move the archaic headline), but keeping fine-grained counts
/// means a fine-pop percentile is a later re-keying rather than an asset rebuild.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ArchaicCountDistribution {
    pub build: String,
    /// The panel these counts were produced from — a percentile is only meaningful against the same
    /// site list, so a mismatch must be detected rather than silently rendered.
    pub panel_sites: usize,
    pub cohorts: Vec<CohortCounts>,
}

impl ArchaicCountDistribution {
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, AnalysisError> {
        bincode::deserialize(bytes).map_err(|e| AnalysisError::Message(format!("archaic dist decode: {e}")))
    }

    pub fn to_bytes(&self) -> Result<Vec<u8>, AnalysisError> {
        bincode::serialize(self).map_err(|e| AnalysisError::Message(format!("archaic dist encode: {e}")))
    }

    /// Percentile (0–100) of `count` within a super-population cohort, or `None` when that cohort is
    /// absent. Reported as the fraction of reference samples scoring **strictly below** `count`.
    pub fn percentile_in_super(&self, super_population: &str, count: u32) -> Option<f32> {
        let mut below = 0usize;
        let mut total = 0usize;
        for c in self.cohorts.iter().filter(|c| c.super_population == super_population) {
            below += c.counts.iter().filter(|&&v| v < count).count();
            total += c.counts.len();
        }
        (total > 0).then(|| below as f32 * 100.0 / total as f32)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn calls(a: ArchaicCall, v: ArchaicCall, c: ArchaicCall, d: ArchaicCall) -> [ArchaicCall; 4] {
        [a, v, c, d]
    }

    #[test]
    fn classify_splits_the_two_lineages() {
        use ArchaicCall::*;
        // Derived in Neanderthals only.
        assert_eq!(
            classify_diagnostic(&calls(HomDerived, HomDerived, NoCall, HomAncestral)),
            DiagnosticClass::Neanderthal
        );
        // Derived in Denisova only.
        assert_eq!(
            classify_diagnostic(&calls(HomAncestral, HomAncestral, HomAncestral, HomDerived)),
            DiagnosticClass::Denisovan
        );
        // Derived in both → not attributable.
        assert_eq!(
            classify_diagnostic(&calls(HomDerived, HomAncestral, HomAncestral, Het)),
            DiagnosticClass::SharedArchaic
        );
        // A heterozygous Neanderthal still counts as carrying the derived allele.
        assert_eq!(
            classify_diagnostic(&calls(Het, NoCall, NoCall, HomAncestral)),
            DiagnosticClass::Neanderthal
        );
    }

    #[test]
    fn derived_copies_counts_the_derived_allele_not_the_alt() {
        let site = ArchaicSite {
            contig: "chr1".into(),
            position: 100,
            reference_allele: 'A',
            alternate_allele: 'G',
            // The derived allele here is the REFERENCE base, which is exactly the case a
            // dosage-based count would get backwards.
            archaic_derived_allele: 'A',
            calls: [ArchaicCall::HomDerived; 4],
            diagnostic_class: DiagnosticClass::SharedArchaic,
            afr_freq: 0.0,
        };
        assert_eq!(site.derived_copies('A', 'A'), 2);
        assert_eq!(site.derived_copies('A', 'G'), 1);
        assert_eq!(site.derived_copies('G', 'G'), 0);
        // Case-insensitive, since VCF/FASTA sources are not consistent about case.
        assert_eq!(site.derived_copies('a', 'g'), 1);
    }

    #[test]
    fn panel_round_trips_and_reports_the_copy_denominator() {
        let panel = ArchaicMarkerPanel {
            build: "chm13v2.0".into(),
            thresholds: ArchaicPanelThresholds {
                max_afr_freq: 0.01,
                min_non_afr_freq: 0.05,
            },
            sites: vec![ArchaicSite {
                contig: "chr1".into(),
                position: 100,
                reference_allele: 'A',
                alternate_allele: 'G',
                archaic_derived_allele: 'G',
                calls: [ArchaicCall::HomDerived; 4],
                diagnostic_class: DiagnosticClass::SharedArchaic,
                afr_freq: 0.002,
            }],
        };
        let bytes = panel.to_bytes().expect("encode");
        assert_eq!(ArchaicMarkerPanel::from_bytes(&bytes).expect("decode"), panel);
        assert_eq!(panel.possible_copies(), 2);
    }

    #[test]
    fn percentile_is_share_of_cohort_scoring_below() {
        let dist = ArchaicCountDistribution {
            build: "chm13v2.0".into(),
            panel_sites: 10,
            cohorts: vec![
                CohortCounts {
                    population: "GBR".into(),
                    super_population: "EUR".into(),
                    counts: vec![10, 20, 30, 40],
                },
                CohortCounts {
                    population: "YRI".into(),
                    super_population: "AFR".into(),
                    counts: vec![0, 1],
                },
            ],
        };
        // 2 of 4 EUR samples score below 30.
        assert_eq!(dist.percentile_in_super("EUR", 30), Some(50.0));
        assert_eq!(dist.percentile_in_super("EUR", 0), Some(0.0));
        assert_eq!(dist.percentile_in_super("EUR", 100), Some(100.0));
        // Cohorts do not bleed into each other, and an unknown one is None (not a fake 0).
        assert_eq!(dist.percentile_in_super("AFR", 1), Some(50.0));
        assert_eq!(dist.percentile_in_super("EAS", 10), None);
    }
}
