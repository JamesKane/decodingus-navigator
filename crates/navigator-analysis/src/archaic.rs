//! Archaic ancestry, from Neanderthals and Denisovans. This module holds the panel types.
//!
//! The design is in `documents/design/ArchaicAncestry_Design.md`. This module carries the asset
//! types of **Tier A**, which are the panel of archaic-informative markers and the percentile
//! reference. The routine that counts, and the segment HMM of Tier B, go on top of them in later
//! milestones.
//!
//! This project *computes* the panel, and does not take it from elsewhere. See design §3a. The
//! Sankararaman 2014 list that 23andMe used is not available to download. Every alternative with
//! an open licence is either on the wrong build, or probabilistic for each individual. So the
//! code makes the sites from the EVA archaic VCFs, the Ensembl-75 EPO ancestral alleles, and a
//! 1kGP AFR outgroup. It distributes only the derived sites.
//!
//! Keep one consequence in mind during a check. **The panel belongs to this project, so its
//! marker count is relative to the panel.** You can not compare it to the count of a vendor by
//! equality. Compare the rate at each site on the intersection instead. See design §10, the M2
//! gate.

use serde::{Deserialize, Serialize};

use crate::caller::SiteGenotype;
use crate::error::AnalysisError;
use crate::ibd_panel::Locus;

/// The four openly-downloadable archaic reference genomes, in the fixed order every
/// [`ArchaicSite::calls`] array uses.
pub const ARCHAIC_GENOMES: [&str; 4] = ["AltaiNeanderthal", "Vindija33.19", "Chagyrskaya8", "Denisova3"];

/// The index into [`ARCHAIC_GENOMES`] and [`ArchaicSite::calls`] of the one Denisovan genome. The
/// other three are Neanderthals, and [`classify_diagnostic`] uses that fact.
pub const DENISOVA: usize = 3;

/// The state of one archaic genome at a site. It is **against the derived allele of the site**,
/// and not against ref and alt. It stays correct through the exchange of ref and alt that the
/// CHM13 orientation can apply.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ArchaicCall {
    /// Homozygous for the ancestral allele.
    HomAncestral,
    /// One derived copy.
    Het,
    /// Homozygous derived. This is the state of the introgression donor. See design §4, step 2.
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

/// The archaic lineage that the derived allele of a site points to.
///
/// The HMM in Tier B can not separate Neanderthal from Denisovan on its own, because the two
/// coalesce before either one meets modern humans. See design §3. This classification is what
/// lets a later step put a label on a called segment, so the build stores it at each site.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum DiagnosticClass {
    /// Derived in ≥1 Neanderthal, and Denisova positively **called ancestral**.
    Neanderthal,
    /// Derived in Denisova, and ≥1 Neanderthal positively **called ancestral**.
    Denisovan,
    /// The code can not attribute this site to one lineage. Either both lineages are derived, or
    /// the other lineage had no call, so nothing shows that it is absent.
    SharedArchaic,
}

/// Classify a site from the calls of each genome.
///
/// A site belongs to one lineage only with **positive evidence that the allele is absent** in the
/// other lineage. The code must *call* that other lineage homozygous-ancestral. A missing call is
/// not enough.
///
/// An earlier version read `NoCall` as absence, and that was the largest error in the Tier B
/// attribution. A site where a mask had removed the Neanderthals, and where the caller did call
/// Denisova, read as *specific* to Denisova. That raised the Denisovan-diagnostic sites to 18,551
/// against 24,077 Neanderthal ones on chr21 and chr22. It gave about 19% Denisovan for a
/// European, where design §7 expects about zero.
///
/// A site that the code can not attribute falls to [`DiagnosticClass::SharedArchaic`]. So that
/// class means "archaic, but the code can not attribute it". It does not mean strictly "derived
/// in both".
pub fn classify_diagnostic(calls: &[ArchaicCall; 4]) -> DiagnosticClass {
    let nea_derived = calls
        .iter()
        .enumerate()
        .any(|(i, c)| i != DENISOVA && c.carries_derived());
    let nea_ancestral = calls
        .iter()
        .enumerate()
        .any(|(i, c)| i != DENISOVA && *c == ArchaicCall::HomAncestral);
    let den_derived = calls[DENISOVA].carries_derived();
    let den_ancestral = calls[DENISOVA] == ArchaicCall::HomAncestral;

    // A site that is derived in BOTH lineages goes straight to shared. Without this test, the
    // Denisovan branch fires whenever some *other* Neanderthal is ancestral at a site that both
    // lineages carry. The four genomes disagree with each other all the time, so that is common
    // and not a rare case.
    if nea_derived && den_derived {
        DiagnosticClass::SharedArchaic
    } else if nea_derived && den_ancestral {
        DiagnosticClass::Neanderthal
    } else if den_derived && nea_ancestral {
        DiagnosticClass::Denisovan
    } else {
        DiagnosticClass::SharedArchaic
    }
}

/// One archaic-informative marker.
///
/// The coordinates are CHM13, and they are **oriented**: `reference_allele` is the true CHM13
/// base at `position`.
///
/// That orientation is necessary. The pipeline lifts GRCh37 to CHM13 with `CrossMap bed`, which
/// does not know about alleles. So about a third of the sites arrive with ref and alt the wrong
/// way round against CHM13. To send them out without the orientation is the bug that cost the
/// ancient-ancestry build a full rebuild. See §7.16 of that design.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ArchaicSite {
    pub contig: String,
    /// 1-based, on the panel's `build`.
    pub position: i64,
    /// The actual base on the panel's build (post-orientation).
    pub reference_allele: char,
    pub alternate_allele: char,
    /// The archaic-derived allele. It is always `reference_allele` or `alternate_allele`. The
    /// field holds a base and not a ref-or-alt flag. A later orientation pass can then not turn
    /// it around where nobody looks.
    pub archaic_derived_allele: char,
    /// The state of each genome, at the index that [`ARCHAIC_GENOMES`] gives.
    pub calls: [ArchaicCall; 4],
    pub diagnostic_class: DiagnosticClass,
    /// Derived-allele frequency in the African outgroup, kept for transparency and so the panel can
    /// be re-filtered at a stricter threshold without a rebuild.
    pub afr_freq: f32,
    /// GRCh37 locus. Exact, not lifted: these are the archaic VCFs' own hg19 coordinates and
    /// alleles, so there is no liftover and no strand risk on this build.
    #[serde(default)]
    pub grch37: Option<Locus>,
    /// The GRCh38 locus. It is present when the site lifted cleanly, and when the code could
    /// orient it against an hg38 reference.
    #[serde(default)]
    pub grch38: Option<Locus>,
}

/// The complement of a base. Use it to compare alleles across builds that can hold different
/// strands.
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
    /// The locus for a build name. It has the same shape as
    /// [`crate::ibd_panel::IbdPanelSite::locus`].
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
    /// A genotype run on a GRCh37 or GRCh38 alignment counts the alleles of that build. Those
    /// alleles can have ref and alt the other way round against CHM13, or sit on the other
    /// strand. To put such a dosage into the count without a change would turn those sites
    /// around, and nobody would see it. That is the same class of error that the CHM13
    /// orientation pass prevents. Returns `None` when the measured allele matches neither CHM13
    /// allele.
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
    /// It takes the two alleles of the subject as bases, and not as a dosage. A dosage counts
    /// against ref and alt, and the derived allele can be either one of those.
    pub fn derived_copies(&self, allele_a: char, allele_b: char) -> u8 {
        let d = self.archaic_derived_allele.to_ascii_uppercase();
        u8::from(allele_a.to_ascii_uppercase() == d) + u8::from(allele_b.to_ascii_uppercase() == d)
    }
}

/// The site-selection thresholds that made a panel. The asset itself carries them.
///
/// The asset records them because they are the scientific content of the panel. The site list
/// *is* the product, and every number after it inherits that list. A count says nothing unless
/// you know which filter made it. Design §10 sets these thresholds by calibration, at M1
/// checkpoint A. It does not set them from the example values in §4.
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

    /// The denominator of a Tier-A count: two copies at each site.
    ///
    /// 23andMe reports exactly this shape. The ground-truth report of the check sample of this
    /// project reads "191 of 7,462", and 7,462 is 2 × 3,731 assayed sites. See design §10.
    pub fn possible_copies(&self) -> usize {
        self.sites.len() * 2
    }
}

/// The Tier-A result: archaic-derived allele copies carried, over copies assayed.
///
/// This is a **count over the sites that the caller called**, and that choice is deliberate. It
/// is not a "percent Neanderthal". It is the same shape that a consumer report uses: "191 of
/// 7,462" is the copies out of 2 × 3,731 assayed sites. See design §1.
///
/// The denominator is whatever subset of the panel the data of the subject covered. So a chip
/// and a WGS run each get an honest headline, and no comparison between the two data types
/// occurs.
///
/// [`percentile`](Self::percentile) is the one field that *does* need such a comparison, so it is
/// an `Option`. The caller fills it only when the coverage of the subject is comparable to the
/// reference cohort. See design §10.
///
/// A chip covers about 3 to 4% of the panel, and it leans toward the common tail of the panel. To
/// rank a chip count against a cohort that a WGS run scored says nothing. The field stays `None`
/// there, and the report does not show a number that looks authoritative and is not.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ArchaicMarkerResult {
    /// Archaic-derived allele copies carried across called panel sites.
    pub total_copies: u32,
    /// 2 × the count of called sites. This is the denominator of the report.
    pub possible_copies: u32,
    /// Panel sites with a usable genotype for this subject.
    pub called_sites: usize,
    /// Sites in the panel overall (so the UI can state the coverage honestly).
    pub panel_sites: usize,
    /// `called_sites / panel_sites`.
    pub call_rate: f32,
    /// Copies at Neanderthal-diagnostic sites.
    pub neanderthal_copies: u32,
    /// The copies at Denisovan-diagnostic sites. For a European this must lie near the noise
    /// floor. §7 does not let the report show a small value here as Denisovan ancestry.
    pub denisovan_copies: u32,
    /// The copies at sites that are derived in both lineages. They are archaic, but the code can
    /// not attribute them.
    pub shared_copies: u32,
    /// Percentile within `cohort`, when the comparison is valid (see the type docs).
    pub percentile: Option<f32>,
    /// The cohort that `percentile` counts against, for example "EUR".
    pub cohort: Option<String>,
    /// The indices, in panel order, of the sites that the caller called for this subject.
    ///
    /// The code needs these to score the reference cohort over the *same* sites. That is the
    /// whole basis of an honest percentile when the input is sparse. The field is **not
    /// serialized**, and that is deliberate. It holds about 300k integers on a WGS subject, and
    /// the cache holds the result as JSON.
    #[serde(skip)]
    pub called_indices: Vec<u32>,
}

impl ArchaicMarkerResult {
    /// The copies that the subject carries, as a fraction of the copies assayed. You can compare
    /// this rate *inside* one data type. Returns 0.0 when the caller called nothing.
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
/// This is dosage arithmetic over the consensus genotypes and nothing more, so it works the same
/// way for chip input and for WGS input.
///
/// There is one thing to watch. `SiteGenotype::dosage` counts the **alternate** allele, and the
/// derived allele of the panel can sit on either side. So the code expresses the dosage again,
/// against the derived *base*. To read the dosage directly as "archaic copies" would turn around
/// every site where the CHM13 orientation left the derived allele on REF. That is 3% of the
/// panel.
pub fn count_archaic_markers(genotypes: &[SiteGenotype], panel: &ArchaicMarkerPanel) -> ArchaicMarkerResult {
    let by_pos: std::collections::HashMap<(&str, i64), &SiteGenotype> =
        genotypes.iter().map(|g| ((g.contig.as_str(), g.position), g)).collect();

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
            // The alleles of the subject at this position disagree with those of the panel.
            // Skip the site, do not guess, and do not count it in the denominator.
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
    /// The total of archaic-derived copies for each sample, from the lowest up.
    pub counts: Vec<u32>,
}

/// The percentile reference asset.
///
/// The asset stores the counts **for each population**, and does not reduce them to
/// super-populations first. v1 shows the percentile against a super-population cohort. See design
/// §9 Q3. A key on the fine ancestry of the user would let an error in that ancestry move the
/// archaic headline where nobody sees it. The code infers that fine ancestry.
///
/// The fine-grained counts stay in the asset. A percentile for a fine population is then a change
/// of key later, and not a rebuild of the asset.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ArchaicCountDistribution {
    pub build: String,
    /// The panel that made these counts. A percentile is correct only against the same site
    /// list. So the code must find a mismatch, and it must not show a percentile that came from
    /// the wrong list.
    pub panel_sites: usize,
    pub cohorts: Vec<CohortCounts>,
    /// Super-population codes indexing [`site_freqs`](Self::site_freqs) and
    /// [`variance_inflation`](Self::variance_inflation).
    #[serde(default)]
    pub populations: Vec<String>,
    /// The derived-allele frequency of each population at each panel site, in **panel order**:
    /// `site_freqs[pop][site_index]`.
    ///
    /// This is what lets a *sparse* subject get an honest percentile. A total for each sample can
    /// rank only a person whom the code scored on the whole panel. Frequencies instead give the
    /// expected count and the variance of the cohort over exactly the sites that a given subject
    /// called, whatever those sites are. That is the only correct comparison when a chip covers
    /// about 3% of the panel, and those sites are its common tail.
    #[serde(default)]
    pub site_freqs: Vec<Vec<f32>>,
    /// The variance inflation of each population, measured at a **ladder of site densities**:
    /// `variance_inflation[pop] = [(density, inflation), …]`, with the densest first.
    ///
    /// Archaic alleles travel in linked haplotype blocks, so the sites are not independent, and
    /// the binomial sum gives a spread that is too small. But the inflation is **not a constant**,
    /// and that is the important part. A measurement gave 52.4x on the full panel, and 5.3x on a
    /// 2.6% subset of the same panel. A sparse subset takes fewer sites from each linked block. To
    /// apply one full-panel factor to a chip would make the deviation about 3x too wide, and it
    /// would push every percentile toward 50.
    ///
    /// No simple block model agrees with the two measurements. A solution for one common block
    /// size gives a negative size. So the code measures the inflation at some densities, and
    /// interpolates between them in log space at run time. It does not model it.
    #[serde(default)]
    pub variance_inflation: Vec<Vec<(f32, f32)>>,
    /// The SHA-256 of the panel asset that these frequencies count against. The UI shows a
    /// percentile only when this value matches the panel that the app loaded. An index into
    /// site_freqs is a panel position. A panel that does not match would then score against the
    /// wrong sites, and nobody would see it.
    #[serde(default)]
    pub panel_fingerprint: String,
}

/// The standard-normal CDF, through the error-function approximation of Abramowitz and Stegun
/// 7.1.26, whose |error| is less than 1.5e-7. That is accurate enough for a percentile that the
/// UI shows as a whole number.
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

    /// The percentile (0 to 100) of `count` inside a super-population cohort, or `None` when the
    /// asset does not hold that cohort. It is the fraction of the reference samples whose score is
    /// **below** `count`, and not equal to it.
    pub fn percentile_in_super(&self, super_population: &str, count: u32) -> Option<f32> {
        let mut below = 0usize;
        let mut total = 0usize;
        for c in self.cohorts.iter().filter(|c| c.super_population == super_population) {
            below += c.counts.iter().filter(|&&v| v < count).count();
            total += c.counts.len();
        }
        (total > 0).then(|| below as f32 * 100.0 / total as f32)
    }

    /// The percentile of `observed_copies`, against `super_population`, for a subject who called
    /// exactly `called_sites`. Those are indices into the panel, in panel order.
    ///
    /// This works at any coverage. For a chip that calls 3% of the panel, the code compares
    /// against the expected count of the cohort *over those same sites*. The artifact of the call
    /// rate, which would otherwise hold every chip user near the 0th percentile, then goes away.
    ///
    /// The model is a sum of independent binomials, one at each site, under Hardy-Weinberg. Each
    /// has a mean of `2f` and a variance of `2f(1-f)`. The code takes the normal approximation of
    /// that sum, which is safe with thousands of sites. It then scales the variance by the
    /// measured LD inflation.
    ///
    /// Returns `None` in four cases. The asset is older than the frequency data. The fingerprint
    /// of the panel disagrees. The code does not know the population. Or the caller called too
    /// few sites for the approximation to say anything.
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
    /// The variance inflation of a population at a given site density. The code interpolates in
    /// log space between the measured rungs, and it clamps the result to the two ends. It never
    /// goes outside the measured range. A value from outside that range is exactly the kind of
    /// number that looks sure and has no basis, which this whole feature avoids.
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

/// Below this count of called sites, the normal approximation is not good enough to show as a
/// percentile.
pub const MIN_SITES_FOR_PERCENTILE: usize = 200;

#[cfg(test)]
mod tests {
    use super::*;

    fn calls(a: ArchaicCall, v: ArchaicCall, c: ArchaicCall, d: ArchaicCall) -> [ArchaicCall; 4] {
        [a, v, c, d]
    }

    #[test]
    fn classify_requires_positive_evidence_of_absence() {
        use ArchaicCall::*;
        // Derived in Neanderthals, Denisova CALLED ancestral -> Neanderthal-diagnostic.
        assert_eq!(
            classify_diagnostic(&calls(HomDerived, HomDerived, NoCall, HomAncestral)),
            DiagnosticClass::Neanderthal
        );
        // Derived in Denisova, a Neanderthal CALLED ancestral -> Denisovan-diagnostic.
        assert_eq!(
            classify_diagnostic(&calls(HomAncestral, NoCall, NoCall, HomDerived)),
            DiagnosticClass::Denisovan
        );
        // Derived in both -> not attributable.
        assert_eq!(
            classify_diagnostic(&calls(HomDerived, HomAncestral, HomAncestral, Het)),
            DiagnosticClass::SharedArchaic
        );
        // A heterozygous Neanderthal still counts, because it carries the derived allele.
        assert_eq!(
            classify_diagnostic(&calls(Het, NoCall, NoCall, HomAncestral)),
            DiagnosticClass::Neanderthal
        );

        // THE REGRESSION THIS GUARDS. Denisova derived, every Neanderthal MASKED OUT: absence is
        // unknown, not established, so this must NOT read as Denisovan-specific. The old rule
        // called it Denisovan and that single mistake produced ~19% Denisovan for a European.
        assert_eq!(
            classify_diagnostic(&calls(NoCall, NoCall, NoCall, HomDerived)),
            DiagnosticClass::SharedArchaic
        );
        // Mirror case: Neanderthal derived, Denisova masked out -> also not attributable.
        assert_eq!(
            classify_diagnostic(&calls(HomDerived, NoCall, NoCall, NoCall)),
            DiagnosticClass::SharedArchaic
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

    fn site(
        position: i64,
        reference_allele: char,
        alternate_allele: char,
        derived: char,
        class: DiagnosticClass,
    ) -> ArchaicSite {
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
        // At site 1 the derived allele is the ALT, so the dosage counts it directly.
        // At site 2 the derived allele is the REF, so the archaic copies are ploidy - dosage. To
        // read the dosage straight through here would turn the site around. That is the failure
        // that this test catches.
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
        // The field is empty, and it is not a lookup by build.
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
        // Panel site 400 is not in the genotypes at all. That is the fourth way that a site gets
        // no call.
        let genotypes = vec![
            gt("chr1", 100, "A", "G", 2),  // counted
            gt("chr1", 200, "A", "G", -1), // explicit no-call
            gt("chr1", 300, "C", "T", 2),  // alleles disagree with the panel
        ];
        let r = count_archaic_markers(&genotypes, &panel);
        assert_eq!(r.called_sites, 1, "only the usable site counts");
        assert_eq!(r.possible_copies, 2);
        assert_eq!(r.total_copies, 2);
        assert_eq!(r.panel_sites, 4);
        assert!((r.call_rate - 0.25).abs() < 1e-6);
        // The counter never invents this value. The caller fills it in only when a comparison is
        // correct.
        assert_eq!(r.percentile, None);
    }

    #[test]
    fn position_stream_round_trips_and_streams_in_order() {
        let positions = vec![10i64, 11, 300, 70_000, 70_001, 5_000_000];
        let st = PositionStream::encode("chr21", &positions);
        assert_eq!(st.len, positions.len());
        assert_eq!(st.iter().collect::<Vec<_>>(), positions);
        // The point of the encoding: gaps, not positions. A run of far-apart sites must still cost
        // far less than 4 bytes each on average for the dense case.
        let dense: Vec<i64> = (0..10_000).map(|i| i * 40).collect();
        let ds = PositionStream::encode("chr21", &dense);
        assert_eq!(ds.iter().collect::<Vec<_>>(), dense);
        assert!(
            ds.deltas.len() < dense.len() * 2,
            "delta encoding should stay ~1 byte/site here"
        );
    }

    #[test]
    fn retain_private_drops_everything_the_outgroup_carries() {
        let og = ArchaicOutgroup {
            build: "chm13v2.0".into(),
            min_allele_count: 1,
            contigs: vec![PositionStream::encode("chr21", &[100, 200, 300, 400])],
        };
        // Africans also carry 200 and 400, so the code removes those two. The rest are private.
        assert_eq!(
            og.retain_private("chr21", &[50, 200, 250, 400, 500]),
            vec![50, 250, 500]
        );
        // Exact-boundary behaviour: first and last outgroup entries.
        assert_eq!(og.retain_private("chr21", &[100, 400]), Vec::<i64>::new());
        // A contig with no outgroup data gives NOTHING, and not everything. To remove no site at
        // all would call the whole contig archaic.
        assert_eq!(og.retain_private("chr7", &[1, 2, 3]), Vec::<i64>::new());
    }

    #[test]
    fn classify_lookup_returns_sites_in_range_with_their_lineage() {
        let cls = ArchaicClassify {
            build: "chm13v2.0".into(),
            contigs: vec![ClassifyContig {
                positions: PositionStream::encode("chr21", &[100, 200, 300]),
                derived: vec![b'A', b'C', b'G'],
                classes: vec![0, 1, 2],
            }],
        };
        let hits = cls.in_range("chr21", 150, 300);
        assert_eq!(hits.len(), 2);
        assert_eq!(hits[0], (200, 'C', DiagnosticClass::Denisovan));
        assert_eq!(hits[1], (300, 'G', DiagnosticClass::SharedArchaic));
        assert!(cls.in_range("chr21", 1, 50).is_empty());
        assert!(cls.in_range("chr7", 1, 1000).is_empty());
    }

    #[test]
    fn subset_percentile_compares_against_the_same_sites() {
        // Two populations. At every site the "high" cohort carries the derived allele at 50%, and
        // the "low" cohort at 5%. The code scores a subject who calls only 1000 sites against the
        // expectation of the cohort OVER THOSE SITES. Sparse coverage then no longer pulls that
        // subject to the bottom.
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

        // The expected copies over 1000 sites at f=0.5 is 1000. A subject at exactly that value
        // sits at the median.
        let p = dist
            .percentile_for_called("HIGH", &called, 1000, "fp")
            .expect("percentile");
        assert!((p - 50.0).abs() < 2.0, "expected ~50th percentile, got {p}");

        // Well above expectation ranks high, well below ranks low.
        assert!(dist.percentile_for_called("HIGH", &called, 1200, "fp").unwrap() > 95.0);
        assert!(dist.percentile_for_called("HIGH", &called, 800, "fp").unwrap() < 5.0);

        // The SAME raw count is usual for HIGH and very rare for LOW. That is the whole point of a
        // score against the correct cohort, and not against one pooled distribution.
        assert!(dist.percentile_for_called("LOW", &called, 1000, "fp").unwrap() > 99.0);

        // The guards are a wrong panel, a population that the code does not know, and too few
        // sites. Each one gives no number at all, and not a wrong one.
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

// ─────────────────────────── Tier B assets (design §4, assets 2 and 3) ───────────────────────────

/// A variable-length integer encoding for a stream of sorted positions. Each byte holds 7 bits,
/// and the high bit says that more bytes follow.
///
/// Both Tier B assets store the **gap between one position and the next**, and not the positions
/// themselves. The African-outgroup track alone holds about 67M positions across the genome. As
/// raw `u32` values that is about 270 MB, which is too much to hold in memory and too much to
/// send to a user. A usual gap is tens of bases, so a varint delta costs one byte at most
/// sites.
fn push_varint(out: &mut Vec<u8>, mut v: u64) {
    while v >= 0x80 {
        out.push((v as u8) | 0x80);
        v >>= 7;
    }
    out.push(v as u8);
}

/// Decode the varint that starts at `i`. Returns the value and the next index.
fn read_varint(bytes: &[u8], mut i: usize) -> Option<(u64, usize)> {
    let (mut v, mut shift) = (0u64, 0u32);
    loop {
        let b = *bytes.get(i)?;
        i += 1;
        v |= ((b & 0x7f) as u64) << shift;
        if b & 0x80 == 0 {
            return Some((v, i));
        }
        shift += 7;
        if shift > 63 {
            return None;
        }
    }
}

/// One contig's sorted positions, delta-varint encoded.
#[derive(Debug, Clone, PartialEq, Eq, Default, Serialize, Deserialize)]
pub struct PositionStream {
    pub contig: String,
    /// Number of positions encoded (so a reader can size buffers without decoding).
    pub len: usize,
    deltas: Vec<u8>,
}

impl PositionStream {
    /// Encode a **sorted, deduplicated** position list.
    pub fn encode(contig: &str, sorted_positions: &[i64]) -> Self {
        let mut deltas = Vec::with_capacity(sorted_positions.len());
        let mut prev = 0i64;
        for &p in sorted_positions {
            push_varint(&mut deltas, (p - prev).max(0) as u64);
            prev = p;
        }
        PositionStream {
            contig: contig.to_string(),
            len: sorted_positions.len(),
            deltas,
        }
    }

    /// Stream the positions back, from the lowest up.
    ///
    /// This gives an iterator and not a `Vec`, and that is deliberate. The caller does a merge
    /// join of it against the sorted variants of the subject, which are much fewer. Neither side
    /// has to be in memory in full.
    pub fn iter(&self) -> impl Iterator<Item = i64> + '_ {
        let mut i = 0usize;
        let mut pos = 0i64;
        std::iter::from_fn(move || {
            let (d, next) = read_varint(&self.deltas, i)?;
            i = next;
            pos += d as i64;
            Some(pos)
        })
    }
}

/// Asset 2: the positions that are **variable in the African outgroup**. The code uses them to
/// remove shared variants before the segment HMM. See design §5, step 1.
///
/// A variant that a modern African population also carries is not evidence of archaic
/// introgression. What stays after the removal is the "private" derived variation, and the HMM
/// looks in that for an excess of density.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ArchaicOutgroup {
    pub build: String,
    /// The allele count that the outgroup needs at a site before the code calls that site
    /// variable there.
    pub min_allele_count: u32,
    pub contigs: Vec<PositionStream>,
}

impl ArchaicOutgroup {
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, AnalysisError> {
        bincode::deserialize(bytes).map_err(|e| AnalysisError::Message(format!("archaic outgroup decode: {e}")))
    }

    pub fn to_bytes(&self) -> Result<Vec<u8>, AnalysisError> {
        bincode::serialize(self).map_err(|e| AnalysisError::Message(format!("archaic outgroup encode: {e}")))
    }

    pub fn contig(&self, contig: &str) -> Option<&PositionStream> {
        self.contigs.iter().find(|c| c.contig == contig)
    }

    /// Keep only the `sorted_positions` that are **not** in the outgroup. That is the private
    /// set.
    ///
    /// This is a linear merge over both sorted streams. It costs O(subject + outgroup), it needs
    /// no index, and it allocates nothing that grows with the outgroup.
    pub fn retain_private(&self, contig: &str, sorted_positions: &[i64]) -> Vec<i64> {
        let Some(stream) = self.contig(contig) else {
            // There is no outgroup data for this contig. To remove no site at all would call the
            // whole contig archaic. So refuse, and do not guess.
            return Vec::new();
        };
        let mut out = Vec::new();
        let mut og = stream.iter().peekable();
        for &p in sorted_positions {
            while og.peek().is_some_and(|&o| o < p) {
                og.next();
            }
            if og.peek() != Some(&p) {
                out.push(p);
            }
        }
        out
    }
}

/// Asset 3: the archaic diagnostic sites across the genome. They put a Neanderthal or Denisovan
/// label on a called segment. See design §5, step 3.
///
/// The HMM alone can not separate the two lineages, because they coalesce before either one meets
/// modern humans. See §3. The attribution is instead a later count of derived-allele matches
/// against these sites.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ArchaicClassify {
    pub build: String,
    pub contigs: Vec<ClassifyContig>,
}

/// One contig's diagnostic sites: positions delta-encoded, with a parallel derived base and class.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ClassifyContig {
    pub positions: PositionStream,
    /// The derived base at each site, in the same order as `positions`.
    pub derived: Vec<u8>,
    /// The diagnostic class at each site: 0 = Neanderthal, 1 = Denisovan, 2 = shared archaic.
    pub classes: Vec<u8>,
}

impl ArchaicClassify {
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, AnalysisError> {
        bincode::deserialize(bytes).map_err(|e| AnalysisError::Message(format!("archaic classify decode: {e}")))
    }

    pub fn to_bytes(&self) -> Result<Vec<u8>, AnalysisError> {
        bincode::serialize(self).map_err(|e| AnalysisError::Message(format!("archaic classify encode: {e}")))
    }

    /// Diagnostic sites on `contig` within `[start, end]`, as `(position, derived_base, class)`.
    pub fn in_range(&self, contig: &str, start: i64, end: i64) -> Vec<(i64, char, DiagnosticClass)> {
        let Some(c) = self.contigs.iter().find(|c| c.positions.contig == contig) else {
            return Vec::new();
        };
        let mut out = Vec::new();
        for (i, p) in c.positions.iter().enumerate() {
            if p > end {
                break;
            }
            if p >= start {
                let class = match c.classes.get(i).copied().unwrap_or(2) {
                    0 => DiagnosticClass::Neanderthal,
                    1 => DiagnosticClass::Denisovan,
                    _ => DiagnosticClass::SharedArchaic,
                };
                out.push((p, c.derived.get(i).copied().unwrap_or(b'N') as char, class));
            }
        }
        out
    }
}

/// The track of callable regions: the count of **callable bases in each fixed-width window**, for
/// each contig.
///
/// The asset holds a count for each window, and not a list of intervals. That is what the segment
/// HMM needs, and it is much smaller. The archaic masks break into hundreds of thousands of
/// intervals below one kb on each chromosome. A grid of 1 kb windows across the genome is about
/// 3.1M `u16` values, which is about 6 MB.
///
/// This track is necessary. Without it, the HMM finds an excess of private-variant density in
/// repetitive regions and reports those regions as archaic. A measurement gave 4,000 to 9,700
/// variants/Mb there, against 50 to 200/Mb for a true introgressed tract.
///
/// Those regions are exactly where the quality masks of the archaic genomes remove the data. So
/// a region that is callable in all four archaic genomes is the only region where an excess of
/// private variants says anything at all.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ArchaicCallable {
    pub build: String,
    /// The window width that the counts go into.
    pub window_bp: i64,
    pub contigs: Vec<CallableContig>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct CallableContig {
    pub contig: String,
    /// Genomic start of window 0.
    pub start: i64,
    /// The count of callable bases in each window. It stops at `window_bp`.
    pub callable_bp: Vec<u16>,
}

impl ArchaicCallable {
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, AnalysisError> {
        bincode::deserialize(bytes).map_err(|e| AnalysisError::Message(format!("archaic callable decode: {e}")))
    }

    pub fn to_bytes(&self) -> Result<Vec<u8>, AnalysisError> {
        bincode::serialize(self).map_err(|e| AnalysisError::Message(format!("archaic callable encode: {e}")))
    }

    pub fn contig(&self, contig: &str) -> Option<&CallableContig> {
        self.contigs.iter().find(|c| c.contig == contig)
    }

    /// The callable fraction (0.0 to 1.0) of the window that holds `position`. It is 0.0 when the
    /// asset holds neither the contig nor the window. A region that the code does not know counts
    /// as **not** callable. The HMM then skips it, and it does not read a density that it can not
    /// trust.
    pub fn callable_fraction(&self, contig: &str, position: i64) -> f64 {
        let Some(c) = self.contig(contig) else { return 0.0 };
        if position < c.start {
            return 0.0;
        }
        let idx = ((position - c.start) / self.window_bp) as usize;
        match c.callable_bp.get(idx) {
            Some(&bp) => (bp as f64 / self.window_bp as f64).clamp(0.0, 1.0),
            None => 0.0,
        }
    }

    /// The total callable megabases over all of the contigs. It is the honest denominator of a
    /// "% of genome" figure.
    pub fn callable_mb(&self) -> f64 {
        self.contigs
            .iter()
            .flat_map(|c| c.callable_bp.iter())
            .map(|&bp| bp as f64)
            .sum::<f64>()
            / 1_000_000.0
    }
}
