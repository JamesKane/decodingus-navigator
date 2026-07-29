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

use crate::caller::SiteGenotype;
use crate::error::AnalysisError;
use crate::ibd_panel::Locus;

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
    /// GRCh37 locus. Exact, not lifted: these are the archaic VCFs' own hg19 coordinates and
    /// alleles, so there is no liftover and no strand risk on this build.
    #[serde(default)]
    pub grch37: Option<Locus>,
    /// GRCh38 locus, when the site lifted cleanly and could be oriented against an hg38 reference.
    #[serde(default)]
    pub grch38: Option<Locus>,
}

/// The complement of a base, for comparing alleles across builds that may differ in strand.
fn complement(b: char) -> char {
    match b.to_ascii_uppercase() {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        other => other,
    }
}

impl ArchaicSite {
    /// The locus for a build name, mirroring [`crate::ibd_panel::IbdPanelSite::locus`].
    pub fn locus(&self, build: &str) -> Option<&Locus> {
        let b = build.to_ascii_lowercase();
        if b.contains("38") || b == "hg38" {
            self.grch38.as_ref()
        } else if b.contains("37") || b == "hg19" || b == "b37" {
            self.grch37.as_ref()
        } else {
            None
        }
    }

    /// Re-express a dosage measured against some build's `alt` allele as a dosage against this
    /// site's **CHM13** alternate allele.
    ///
    /// Genotyping a GRCh37/38 alignment tallies that build's alleles, which may be ref/alt-swapped
    /// or strand-flipped relative to CHM13. Feeding such a dosage into the count unchanged would
    /// invert those sites silently — the same class of error the CHM13 orientation pass exists to
    /// prevent. Returns `None` when the measured allele corresponds to neither CHM13 allele.
    pub fn rekey_dosage(&self, measured_alt: char, dosage: i32, ploidy: u8) -> Option<i32> {
        let m = measured_alt.to_ascii_uppercase();
        let (r, a) = (
            self.reference_allele.to_ascii_uppercase(),
            self.alternate_allele.to_ascii_uppercase(),
        );
        if m == a || m == complement(a) {
            Some(dosage)
        } else if m == r || m == complement(r) {
            Some(ploidy as i32 - dosage)
        } else {
            None
        }
    }

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

/// The Tier-A result: archaic-derived allele copies carried, over copies assayed.
///
/// Deliberately a **count over what was actually called**, not a "percent Neanderthal" — the same
/// shape a consumer report uses ("191 of 7,462" = copies of 2 × 3,731 assayed sites, design §1).
/// Because the denominator is whatever subset of the panel the subject's data covered, chip and WGS
/// each get an honest headline with no cross-data-type comparison involved.
///
/// [`percentile`](Self::percentile) is the one field that *does* require comparability, so it is an
/// `Option` the caller fills only when the subject's coverage is comparable to the reference cohort
/// (design §10). A chip covers ~3–4 % of the panel and is biased toward its common tail, so ranking
/// a chip count against a WGS-scored cohort is meaningless — it stays `None` there rather than
/// rendering a number that looks authoritative and is not.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ArchaicMarkerResult {
    /// Archaic-derived allele copies carried across called panel sites.
    pub total_copies: u32,
    /// 2 × called sites — the denominator the count is reported against.
    pub possible_copies: u32,
    /// Panel sites with a usable genotype for this subject.
    pub called_sites: usize,
    /// Sites in the panel overall (so the UI can state the coverage honestly).
    pub panel_sites: usize,
    /// `called_sites / panel_sites`.
    pub call_rate: f32,
    /// Copies at Neanderthal-diagnostic sites.
    pub neanderthal_copies: u32,
    /// Copies at Denisovan-diagnostic sites. For a European this should be near the noise floor;
    /// §7 forbids presenting a small value here as a positive Denisovan finding.
    pub denisovan_copies: u32,
    /// Copies at sites derived in both lineages — archaic, but not attributable.
    pub shared_copies: u32,
    /// Percentile within `cohort`, when the comparison is valid (see the type docs).
    pub percentile: Option<f32>,
    /// The cohort `percentile` was computed against (e.g. "EUR").
    pub cohort: Option<String>,
    /// Indices (panel order) of the sites that were actually called for this subject.
    ///
    /// Needed to score the reference cohort over the *same* sites — the whole basis of an honest
    /// percentile for sparse input. Deliberately **not serialized**: it is ~300 k integers on a WGS
    /// subject, and the cached result is stored as JSON.
    #[serde(skip)]
    pub called_indices: Vec<u32>,
}

impl ArchaicMarkerResult {
    /// Copies carried as a fraction of copies assayed — the rate that is comparable *within* one
    /// data type. Returns 0.0 when nothing was called.
    pub fn rate(&self) -> f32 {
        if self.possible_copies == 0 {
            0.0
        } else {
            self.total_copies as f32 / self.possible_copies as f32
        }
    }
}

/// Count a subject's archaic-derived allele copies across the panel.
///
/// Pure dosage arithmetic over the consensus genotypes, so it works identically for chip and WGS
/// input. The one subtlety: `SiteGenotype::dosage` counts the **alternate** allele, while the
/// panel's derived allele may sit on either side, so the dosage is re-expressed against the derived
/// *base*. Reading dosage directly as "archaic copies" would invert every site where CHM13
/// orientation left the derived allele on REF — 3 % of the panel.
pub fn count_archaic_markers(genotypes: &[SiteGenotype], panel: &ArchaicMarkerPanel) -> ArchaicMarkerResult {
    let by_pos: std::collections::HashMap<(&str, i64), &SiteGenotype> = genotypes
        .iter()
        .map(|g| ((g.contig.as_str(), g.position), g))
        .collect();

    let (mut total, mut nea, mut den, mut shared) = (0u32, 0u32, 0u32, 0u32);
    let mut called = 0usize;
    let mut called_indices: Vec<u32> = Vec::new();

    for (idx, site) in panel.sites.iter().enumerate() {
        let Some(g) = by_pos.get(&(site.contig.as_str(), site.position)) else {
            continue;
        };
        // A negative dosage is an explicit no-call; ploidy 0 would make the denominator a lie.
        if g.dosage < 0 || g.ploidy == 0 {
            continue;
        }
        let d = site.archaic_derived_allele.to_ascii_uppercase();
        let matches = |s: &str| s.len() == 1 && s.as_bytes()[0].to_ascii_uppercase() as char == d;
        let copies = if matches(&g.alternate_allele) {
            g.dosage as u32
        } else if matches(&g.reference_allele) {
            (g.ploidy as i32 - g.dosage).max(0) as u32
        } else {
            // The subject's alleles at this position disagree with the panel's — skip rather than
            // guess, and do not count it toward the denominator either.
            continue;
        };
        called += 1;
        called_indices.push(idx as u32);
        total += copies;
        match site.diagnostic_class {
            DiagnosticClass::Neanderthal => nea += copies,
            DiagnosticClass::Denisovan => den += copies,
            DiagnosticClass::SharedArchaic => shared += copies,
        }
    }

    ArchaicMarkerResult {
        total_copies: total,
        possible_copies: (called as u32).saturating_mul(2),
        called_sites: called,
        panel_sites: panel.len(),
        call_rate: if panel.is_empty() {
            0.0
        } else {
            called as f32 / panel.len() as f32
        },
        neanderthal_copies: nea,
        denisovan_copies: den,
        shared_copies: shared,
        percentile: None,
        cohort: None,
        called_indices,
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
    /// Super-population codes indexing [`site_freqs`](Self::site_freqs) and
    /// [`variance_inflation`](Self::variance_inflation).
    #[serde(default)]
    pub populations: Vec<String>,
    /// Per-population derived-allele frequency at each panel site, in **panel order**:
    /// `site_freqs[pop][site_index]`.
    ///
    /// This is what lets a *sparse* subject get an honest percentile. Per-sample totals can only
    /// rank someone scored on the whole panel; frequencies let the cohort's expected count and
    /// variance be computed over exactly the sites a given subject called, whatever those are —
    /// which is the only valid comparison when a chip covers ~3 % of the panel and those sites are
    /// its common tail.
    #[serde(default)]
    pub site_freqs: Vec<Vec<f32>>,
    /// Per-population variance inflation measured at a **ladder of site densities**:
    /// `variance_inflation[pop] = [(density, inflation), …]`, densest first.
    ///
    /// Archaic alleles travel in linked haplotype blocks, so sites are not independent and the
    /// binomial sum understates the spread — but crucially the inflation is **not a constant**. It
    /// was measured at 52.4× on the full panel and 5.3× on a 2.6 % subset of the same panel, because
    /// a sparse subset samples fewer sites per linked block. Applying a single full-panel factor to
    /// a chip would over-widen the deviation ~3× and squash every percentile toward 50.
    ///
    /// No simple block model fits the two measurements (solving for a common block size gives a
    /// negative size), so this is measured empirically at several densities and interpolated in log
    /// space at runtime rather than modelled.
    #[serde(default)]
    pub variance_inflation: Vec<Vec<(f32, f32)>>,
    /// SHA-256 of the panel asset these frequencies were computed against. A percentile is only
    /// rendered when this matches the loaded panel — site_freqs is indexed by panel position, so a
    /// mismatched panel would silently score against the wrong sites.
    #[serde(default)]
    pub panel_fingerprint: String,
}

/// Standard-normal CDF via the Abramowitz & Stegun 7.1.26 error-function approximation
/// (|error| < 1.5e-7) — enough for a percentile rendered to whole numbers.
fn normal_cdf(z: f64) -> f64 {
    let sign = if z < 0.0 { -1.0 } else { 1.0 };
    let x = z.abs() / std::f64::consts::SQRT_2;
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let y = 1.0
        - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592)
            * t
            * (-x * x).exp();
    0.5 * (1.0 + sign * y)
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

    /// Percentile of `observed_copies` for a subject who called exactly `called_sites` (indices into
    /// the panel, in panel order), against `super_population`.
    ///
    /// Works at any coverage — a chip calling 3 % of the panel is compared against the cohort's
    /// expected count *over those same sites*, so the call-rate artefact that would otherwise pin
    /// every chip user near the 0th percentile disappears.
    ///
    /// Modelled as a sum of independent per-site binomials under Hardy-Weinberg (mean `2f`, variance
    /// `2f(1-f)`), normal-approximated — with thousands of sites the CLT is comfortable — then the
    /// variance scaled by the measured LD inflation. Returns `None` when the asset predates the
    /// frequency data, the panel fingerprint disagrees, the population is unknown, or too few sites
    /// were called for the approximation to mean anything.
    pub fn percentile_for_called(
        &self,
        super_population: &str,
        called_sites: &[u32],
        observed_copies: u32,
        panel_fingerprint: &str,
    ) -> Option<f32> {
        if self.site_freqs.is_empty() || self.panel_fingerprint != panel_fingerprint {
            return None;
        }
        if called_sites.len() < MIN_SITES_FOR_PERCENTILE {
            return None;
        }
        let pop = self.populations.iter().position(|p| p == super_population)?;
        let freqs = self.site_freqs.get(pop)?;
        let (mut mean, mut var) = (0.0f64, 0.0f64);
        for &i in called_sites {
            let Some(&f) = freqs.get(i as usize) else { continue };
            let f = f as f64;
            mean += 2.0 * f;
            var += 2.0 * f * (1.0 - f);
        }
        let density = called_sites.len() as f32 / self.panel_sites.max(1) as f32;
        var *= self.inflation_at(pop, density) as f64;
        if var <= 0.0 {
            return None;
        }
        let z = (observed_copies as f64 - mean) / var.sqrt();
        Some((normal_cdf(z) * 100.0).clamp(0.0, 100.0) as f32)
    }
}

impl ArchaicCountDistribution {
    /// Variance inflation for a population at a given site density, log-interpolated between the
    /// measured rungs and clamped to the ends (never extrapolated — an extrapolated inflation is
    /// exactly the kind of confident-but-unfounded number this whole feature avoids).
    fn inflation_at(&self, pop: usize, density: f32) -> f32 {
        let Some(ladder) = self.variance_inflation.get(pop) else {
            return 1.0;
        };
        if ladder.is_empty() {
            return 1.0;
        }
        let mut rungs = ladder.clone();
        rungs.sort_by(|a, b| a.0.total_cmp(&b.0));
        if density <= rungs[0].0 {
            return rungs[0].1.max(1.0);
        }
        if density >= rungs[rungs.len() - 1].0 {
            return rungs[rungs.len() - 1].1.max(1.0);
        }
        for w in rungs.windows(2) {
            let ((d0, i0), (d1, i1)) = (w[0], w[1]);
            if density >= d0 && density <= d1 {
                let (l0, l1, l) = (d0.max(1e-6).ln(), d1.max(1e-6).ln(), density.max(1e-6).ln());
                let t = if (l1 - l0).abs() < f32::EPSILON {
                    0.0
                } else {
                    (l - l0) / (l1 - l0)
                };
                return (i0 + t * (i1 - i0)).max(1.0);
            }
        }
        1.0
    }
}

/// Below this many called sites the normal approximation is not worth rendering as a percentile.
pub const MIN_SITES_FOR_PERCENTILE: usize = 200;

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
            grch37: None,
            grch38: None,
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
                grch37: None,
                grch38: None,
            }],
        };
        let bytes = panel.to_bytes().expect("encode");
        assert_eq!(ArchaicMarkerPanel::from_bytes(&bytes).expect("decode"), panel);
        assert_eq!(panel.possible_copies(), 2);
    }

    fn gt(contig: &str, position: i64, reference_allele: &str, alternate_allele: &str, dosage: i32) -> SiteGenotype {
        SiteGenotype {
            name: String::new(),
            contig: contig.into(),
            position,
            reference_allele: reference_allele.into(),
            alternate_allele: alternate_allele.into(),
            ploidy: 2,
            dosage,
            gq: 0,
            depth: 0,
            ref_depth: 0,
            alt_depth: 0,
            pls: Vec::new(),
            gt: None,
            allele_depths: None,
        }
    }

    fn site(position: i64, reference_allele: char, alternate_allele: char, derived: char, class: DiagnosticClass) -> ArchaicSite {
        ArchaicSite {
            contig: "chr1".into(),
            position,
            reference_allele,
            alternate_allele,
            archaic_derived_allele: derived,
            calls: [ArchaicCall::HomDerived; 4],
            diagnostic_class: class,
            afr_freq: 0.001,
            grch37: None,
            grch38: None,
        }
    }

    #[test]
    fn counting_re_expresses_dosage_against_the_derived_base() {
        // Site 1: derived is the ALT, so dosage counts it directly.
        // Site 2: derived is the REF, so the archaic copies are ploidy - dosage. Reading dosage
        // straight through here would invert the site — the failure this test exists to catch.
        let panel = ArchaicMarkerPanel {
            build: "chm13v2.0".into(),
            thresholds: ArchaicPanelThresholds {
                max_afr_freq: 0.01,
                min_non_afr_freq: 0.0005,
            },
            sites: vec![
                site(100, 'A', 'G', 'G', DiagnosticClass::Neanderthal),
                site(200, 'A', 'G', 'A', DiagnosticClass::Denisovan),
            ],
        };
        let genotypes = vec![gt("chr1", 100, "A", "G", 1), gt("chr1", 200, "A", "G", 1)];
        let r = count_archaic_markers(&genotypes, &panel);
        assert_eq!(r.called_sites, 2);
        assert_eq!(r.possible_copies, 4);
        assert_eq!(r.neanderthal_copies, 1, "derived on ALT: dosage passes through");
        assert_eq!(r.denisovan_copies, 1, "derived on REF: 2 - dosage");
        assert_eq!(r.total_copies, 2);
        assert!((r.rate() - 0.5).abs() < 1e-6);
    }

    #[test]
    fn rekey_dosage_handles_swap_and_strand_flip() {
        // CHM13 alleles A/G, derived G.
        let s = site(100, 'A', 'G', 'G', DiagnosticClass::Neanderthal);
        // Same orientation: the other build also calls G the ALT -> dosage passes through.
        assert_eq!(s.rekey_dosage('G', 1, 2), Some(1));
        assert_eq!(s.rekey_dosage('G', 2, 2), Some(2));
        // Ref/alt swapped on the other build (its ALT is our REF) -> dosage inverts.
        assert_eq!(s.rekey_dosage('A', 2, 2), Some(0));
        assert_eq!(s.rekey_dosage('A', 0, 2), Some(2));
        // Strand-flipped: the other build's ALT is complement(G) = C -> still our ALT.
        assert_eq!(s.rekey_dosage('C', 1, 2), Some(1));
        // Strand-flipped AND swapped: complement(A) = T -> our REF, inverts.
        assert_eq!(s.rekey_dosage('T', 2, 2), Some(0));
    }

    #[test]
    fn locus_lookup_matches_build_aliases() {
        let mut s = site(100, 'A', 'G', 'G', DiagnosticClass::Neanderthal);
        s.grch37 = Some(Locus {
            contig: "1".into(),
            position: 42,
            reference: 'A',
            alternate: 'G',
        });
        assert_eq!(s.locus("GRCh37").map(|l| l.position), Some(42));
        assert_eq!(s.locus("hg19").map(|l| l.position), Some(42));
        assert_eq!(s.locus("b37").map(|l| l.position), Some(42));
        // Not populated / not a per-build lookup.
        assert!(s.locus("GRCh38").is_none());
        assert!(s.locus("chm13v2.0").is_none());
    }

    #[test]
    fn uncalled_and_disagreeing_sites_leave_the_denominator_alone() {
        let panel = ArchaicMarkerPanel {
            build: "chm13v2.0".into(),
            thresholds: ArchaicPanelThresholds {
                max_afr_freq: 0.01,
                min_non_afr_freq: 0.0005,
            },
            sites: vec![
                site(100, 'A', 'G', 'G', DiagnosticClass::Neanderthal),
                site(200, 'A', 'G', 'G', DiagnosticClass::Neanderthal),
                site(300, 'A', 'G', 'G', DiagnosticClass::Neanderthal),
                site(400, 'A', 'G', 'G', DiagnosticClass::Neanderthal),
            ],
        };
        let genotypes = vec![
            gt("chr1", 100, "A", "G", 2),   // counted
            gt("chr1", 200, "A", "G", -1),  // explicit no-call
            gt("chr1", 300, "C", "T", 2),   // alleles disagree with the panel
            // 400 absent entirely
        ];
        let r = count_archaic_markers(&genotypes, &panel);
        assert_eq!(r.called_sites, 1, "only the usable site counts");
        assert_eq!(r.possible_copies, 2);
        assert_eq!(r.total_copies, 2);
        assert_eq!(r.panel_sites, 4);
        assert!((r.call_rate - 0.25).abs() < 1e-6);
        // Never fabricated by the counter — the caller fills it only when comparable.
        assert_eq!(r.percentile, None);
    }

    #[test]
    fn subset_percentile_compares_against_the_same_sites() {
        // Two populations. At every site the "high" cohort carries the derived allele at 50% and the
        // "low" cohort at 5%. A subject calling only 1000 sites is scored against the cohort's
        // expectation OVER THOSE SITES, so sparse coverage no longer drags them to the bottom.
        let n = 2000usize;
        let dist = ArchaicCountDistribution {
            build: "chm13v2.0".into(),
            panel_sites: n,
            cohorts: Vec::new(),
            populations: vec!["HIGH".into(), "LOW".into()],
            site_freqs: vec![vec![0.5; n], vec![0.05; n]],
            variance_inflation: vec![vec![(1.0, 1.0)], vec![(1.0, 1.0)]],
            panel_fingerprint: "fp".into(),
        };
        let called: Vec<u32> = (0..1000).collect();

        // Expected copies over 1000 sites at f=0.5 is 1000; landing exactly there is the median.
        let p = dist.percentile_for_called("HIGH", &called, 1000, "fp").expect("percentile");
        assert!((p - 50.0).abs() < 2.0, "expected ~50th percentile, got {p}");

        // Well above expectation ranks high, well below ranks low.
        assert!(dist.percentile_for_called("HIGH", &called, 1200, "fp").unwrap() > 95.0);
        assert!(dist.percentile_for_called("HIGH", &called, 800, "fp").unwrap() < 5.0);

        // The SAME raw count is unremarkable for HIGH but extraordinary for LOW — which is the whole
        // point of scoring against the right cohort rather than a single pooled distribution.
        assert!(dist.percentile_for_called("LOW", &called, 1000, "fp").unwrap() > 99.0);

        // Guards: wrong panel, unknown population, too few sites -> no number rather than a wrong one.
        assert_eq!(dist.percentile_for_called("HIGH", &called, 1000, "other-panel"), None);
        assert_eq!(dist.percentile_for_called("NOPE", &called, 1000, "fp"), None);
        assert_eq!(dist.percentile_for_called("HIGH", &called[..10], 5, "fp"), None);
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
            populations: Vec::new(),
            site_freqs: Vec::new(),
            variance_inflation: Vec::new(),
            panel_fingerprint: String::new(),
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
