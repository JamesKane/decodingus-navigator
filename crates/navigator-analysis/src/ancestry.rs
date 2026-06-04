//! Ancestry estimation — the genotype → population-proportion path, Navigator-side.
//!
//! Phase 1 is the allele-frequency likelihood (no PCA, no GATK): the bundled [`AncestryPanel`]
//! carries per-(super-)population alt-allele frequencies at a set of ancestry-informative
//! sites; we genotype the sample there with the GL caller ([`crate::caller::genotype_sites`]),
//! then score each population by the binomial likelihood of the observed diploid genotypes
//! under its allele frequencies. The panel is built offline by `navigator-panelbuild` from the
//! 1000G-on-CHM13 VCFs.
//!
//! The result is a [`navigator_domain::ancestry::AncestryResult`]. PCA projection
//! ([`AncestryResult::pca_coordinates`]) is phase 2.

use std::collections::{BTreeMap, HashMap};
use std::path::Path;

use navigator_domain::ancestry::{
    population_color, population_name, population_super, AncestryResult, ConfidenceInterval,
    PopulationComponent, SuperPopulationSummary,
};
use serde::{Deserialize, Serialize};

use crate::caller::{self, HaploidCallerParams, Site, SiteGenotype};
use crate::AnalysisError;

/// One ancestry-informative site with its per-population alt-allele frequencies. `freqs[i]`
/// aligns with [`AncestryPanel::populations`]`[i]`.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PanelSite {
    pub contig: String,
    /// 1-based.
    pub position: i64,
    pub reference_allele: char,
    pub alternate_allele: char,
    pub freqs: Vec<f32>,
}

/// A bundled ancestry reference panel: the populations axis plus the informative sites.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct AncestryPanel {
    /// Canonical reference build the site coordinates are in (e.g. "chm13v2.0").
    pub build: String,
    /// Population codes, defining the axis order of every `PanelSite::freqs`.
    pub populations: Vec<String>,
    pub sites: Vec<PanelSite>,
}

impl AncestryPanel {
    /// Deserialize from the bundled/built binary (bincode).
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, AnalysisError> {
        bincode::deserialize(bytes).map_err(|e| AnalysisError::Message(format!("panel decode: {e}")))
    }

    /// Serialize to the binary form the builder writes and the app bundles.
    pub fn to_bytes(&self) -> Result<Vec<u8>, AnalysisError> {
        bincode::serialize(self).map_err(|e| AnalysisError::Message(format!("panel encode: {e}")))
    }

    pub fn len(&self) -> usize {
        self.sites.len()
    }
    pub fn is_empty(&self) -> bool {
        self.sites.is_empty()
    }
}

/// PCA loadings for projecting a sample onto the reference populations' principal-component
/// space (Phase 2). Built offline by `navigator-panelbuild` from the 1000G genotype matrix:
/// per-SNP loadings + means (for centering), plus each population's centroid and diagonal
/// variance in PC space (for the Mahalanobis/Gaussian assignment and the scatter plot).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PcaLoadings {
    pub build: String,
    /// (contig, 1-based pos) per row, aligned with `means` and the rows of `loadings`.
    pub sites: Vec<(String, i64)>,
    /// Mean dosage per site (reference panel) — used to centre the sample before projecting.
    pub means: Vec<f32>,
    pub n_components: usize,
    /// Row-major `sites.len() × n_components`.
    pub loadings: Vec<f32>,
    pub populations: Vec<String>,
    /// Row-major `populations.len() × n_components`.
    pub centroids: Vec<f32>,
    /// Row-major `populations.len() × n_components` (diagonal covariance).
    pub variances: Vec<f32>,
}

impl PcaLoadings {
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, AnalysisError> {
        bincode::deserialize(bytes).map_err(|e| AnalysisError::Message(format!("pca decode: {e}")))
    }
    pub fn to_bytes(&self) -> Result<Vec<u8>, AnalysisError> {
        bincode::serialize(self).map_err(|e| AnalysisError::Message(format!("pca encode: {e}")))
    }
    pub fn loading(&self, site_idx: usize, component: usize) -> f32 {
        self.loadings[site_idx * self.n_components + component]
    }
    pub fn centroid(&self, pop_idx: usize) -> &[f32] {
        let o = pop_idx * self.n_components;
        &self.centroids[o..o + self.n_components]
    }
    pub fn variance(&self, pop_idx: usize) -> &[f32] {
        let o = pop_idx * self.n_components;
        &self.variances[o..o + self.n_components]
    }
}

/// Project a sample's genotypes onto the reference PCA space: centre each site by its panel
/// mean and accumulate `centered · loading` into each component. A missing genotype contributes
/// 0 (mean-imputed). Returns the sample's coordinate in each principal component.
pub fn project_pca(genotypes: &[SiteGenotype], pca: &PcaLoadings) -> Vec<f64> {
    let dosage: HashMap<(&str, i64), i32> = genotypes
        .iter()
        .filter(|g| g.dosage >= 0)
        .map(|g| ((g.contig.as_str(), g.position), g.dosage))
        .collect();

    let mut coords = vec![0.0f64; pca.n_components];
    for (i, (contig, pos)) in pca.sites.iter().enumerate() {
        let centered = match dosage.get(&(contig.as_str(), *pos)) {
            Some(&d) => d as f64 - pca.means[i] as f64,
            None => continue, // mean-imputed → centred value 0, no contribution
        };
        for (c, coord) in coords.iter_mut().enumerate() {
            *coord += centered * pca.loading(i, c) as f64;
        }
    }
    coords
}

/// Squared Mahalanobis distance (diagonal covariance) from `coords` to a population centroid.
fn mahalanobis_sq(coords: &[f64], centroid: &[f32], variance: &[f32]) -> f64 {
    coords
        .iter()
        .zip(centroid)
        .zip(variance)
        .map(|((&x, &mu), &v)| {
            let d = x - mu as f64;
            let v = v as f64;
            if v > 1e-9 {
                d * d / v
            } else {
                0.0
            }
        })
        .sum()
}

/// Per-population probability that the sample's PCA `coords` belong to each population, by a
/// diagonal-covariance Gaussian (`exp(-½ d²)`), normalized to sum 1. Useful as a PCA-based
/// cross-check of the allele-frequency estimate.
pub fn classify_pca(coords: &[f64], pca: &PcaLoadings) -> Vec<(String, f64)> {
    let raw: Vec<(String, f64)> = pca
        .populations
        .iter()
        .enumerate()
        .map(|(p, code)| {
            let d2 = mahalanobis_sq(coords, pca.centroid(p), pca.variance(p));
            (code.clone(), (-0.5 * d2).exp())
        })
        .collect();
    let total: f64 = raw.iter().map(|(_, p)| p).sum();
    if total > 0.0 {
        raw.into_iter().map(|(c, p)| (c, p / total)).collect()
    } else {
        raw
    }
}

/// Genotype a BAM/CRAM at every panel site (diploid — the panel is autosomal AIMs). Groups
/// sites by contig and runs the GL caller once per contig; returns the per-site genotypes
/// (dosage 0/1/2, or -1 for a no-call). `reference` is required for CRAM.
///
/// Each contig is one full read-scan, so on a whole-genome BAM this is the slow step (minutes);
/// `progress(contigs_done, contigs_total)` is invoked after each contig so the UI can show a
/// bar. Contigs are processed in sorted order so progress is monotonic.
pub fn genotype_panel(
    bam: &Path,
    reference: Option<&Path>,
    panel: &AncestryPanel,
    params: &HaploidCallerParams,
    progress: &mut dyn FnMut(usize, usize),
) -> Result<Vec<SiteGenotype>, AnalysisError> {
    // contig -> caller sites (BTreeMap → deterministic, monotonic progress order)
    let mut by_contig: BTreeMap<String, Vec<Site>> = BTreeMap::new();
    for s in &panel.sites {
        by_contig.entry(s.contig.clone()).or_default().push(Site {
            name: format!("{}:{}", s.contig, s.position),
            contig: s.contig.clone(),
            position: s.position,
            reference_allele: s.reference_allele.to_string(),
            alternate_allele: s.alternate_allele.to_string(),
        });
    }

    let total = by_contig.len();
    let mut out = Vec::with_capacity(panel.sites.len());
    for (done, (contig, sites)) in by_contig.into_iter().enumerate() {
        let calls = caller::genotype_sites(bam, &contig, &sites, 2, params, reference)?;
        out.extend(calls);
        progress(done + 1, total);
    }
    Ok(out)
}

const PIPELINE_VERSION: &str = "1.0.0-af";

/// Estimate ancestry by the per-population binomial allele-frequency likelihood.
///
/// For each population, the log-likelihood sums `ln P(genotype | f)` over genotyped sites,
/// where `f` is that population's alt-allele frequency (clamped to [0.001, 0.999]) and the
/// diploid genotype probability is `(1-f)²` (hom-ref), `2f(1-f)` (het), or `f²` (hom-alt).
/// Likelihoods are exponentiated relative to the best population and normalized to percentages.
pub fn estimate_by_allele_frequency(
    genotypes: &[SiteGenotype],
    panel: &AncestryPanel,
    reference_version: &str,
) -> AncestryResult {
    // (contig, position) -> dosage; missing/no-call dosages (< 0) are dropped.
    let dosage: HashMap<(&str, i64), i32> = genotypes
        .iter()
        .filter(|g| g.dosage >= 0)
        .map(|g| ((g.contig.as_str(), g.position), g.dosage))
        .collect();

    let num_pops = panel.populations.len();
    let mut logl = vec![0.0f64; num_pops];
    let mut snps_with_data = 0usize;

    for site in &panel.sites {
        let Some(&g) = dosage.get(&(site.contig.as_str(), site.position)) else {
            continue;
        };
        if site.freqs.len() != num_pops {
            continue; // malformed site
        }
        snps_with_data += 1;
        for (pop_idx, &f_raw) in site.freqs.iter().enumerate() {
            let f = (f_raw as f64).clamp(0.001, 0.999);
            let p = match g {
                0 => (1.0 - f) * (1.0 - f),
                1 => 2.0 * f * (1.0 - f),
                2 => f * f,
                _ => 1.0,
            };
            logl[pop_idx] += p.max(1e-300).ln();
        }
    }

    // Exponentiate relative to the best population (numerical stability), then normalize.
    let max_ll = logl.iter().cloned().fold(f64::MIN, f64::max);
    let probs: Vec<(String, f64)> = panel
        .populations
        .iter()
        .zip(logl.iter())
        .map(|(code, &ll)| (code.clone(), (ll - max_ll).exp()))
        .collect();

    let confidence = confidence_from_completeness(snps_with_data, panel.sites.len());
    from_probabilities(
        "aims",
        panel.sites.len(),
        snps_with_data,
        &probs,
        confidence,
        reference_version,
    )
}

/// Build an [`AncestryResult`] from raw per-population probabilities (need not be normalized).
/// With the phase-1 super-population panel each component *is* a super-population, so the
/// super-population summary is 1:1 with the components.
fn from_probabilities(
    panel_type: &str,
    snps_analyzed: usize,
    snps_with_genotype: usize,
    population_probs: &[(String, f64)],
    confidence_level: f64,
    reference_version: &str,
) -> AncestryResult {
    let total: f64 = population_probs.iter().map(|(_, p)| p).sum();
    let mut pct: Vec<(String, f64)> = population_probs
        .iter()
        .map(|(code, p)| (code.clone(), if total > 0.0 { p / total * 100.0 } else { 0.0 }))
        .collect();
    pct.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

    let components: Vec<PopulationComponent> = pct
        .iter()
        .enumerate()
        .map(|(idx, (code, p))| {
            let ci = ci_width(*p, snps_with_genotype, snps_analyzed);
            PopulationComponent {
                population_code: code.clone(),
                population_name: population_name(code),
                percentage: *p,
                confidence_interval: ConfidenceInterval {
                    lower: (p - ci).max(0.0),
                    upper: (p + ci).min(100.0),
                },
                rank: idx + 1,
            }
        })
        .collect();

    // Roll components up into super-population summaries. With a super-population panel each
    // component is its own super-population; with a fine-grained (26-pop) panel several
    // components aggregate into one super-population.
    let mut by_super: BTreeMap<String, (f64, Vec<String>)> = BTreeMap::new();
    for (code, p) in &pct {
        let sp = population_super(code).unwrap_or(code.as_str()).to_string();
        let e = by_super.entry(sp).or_insert((0.0, Vec::new()));
        e.0 += *p;
        e.1.push(code.clone());
    }
    let mut super_population_summary: Vec<SuperPopulationSummary> = by_super
        .into_iter()
        .map(|(sp, (pct, members))| SuperPopulationSummary {
            super_population: population_name(&sp),
            percentage: pct,
            populations: members,
        })
        .collect();
    super_population_summary.sort_by(|a, b| b.percentage.partial_cmp(&a.percentage).unwrap_or(std::cmp::Ordering::Equal));

    // Touch the catalog color path so the API stays cohesive; color is consumed by the UI.
    debug_assert!(!population_color("EUR").is_empty());

    AncestryResult {
        panel_type: panel_type.to_string(),
        snps_analyzed,
        snps_with_genotype,
        snps_missing: snps_analyzed.saturating_sub(snps_with_genotype),
        components,
        super_population_summary,
        confidence_level,
        pipeline_version: PIPELINE_VERSION.to_string(),
        reference_version: reference_version.to_string(),
        pca_coordinates: None,
    }
}

/// Binomial-proportion CI half-width (percent), widened for incomplete panels.
fn ci_width(pct: f64, snps_with_data: usize, total_snps: usize) -> f64 {
    let completeness = if total_snps == 0 { 0.0 } else { snps_with_data as f64 / total_snps as f64 };
    let p = pct / 100.0;
    let base = if snps_with_data > 0 {
        1.96 * (p * (1.0 - p) / snps_with_data as f64).sqrt() * 100.0
    } else {
        50.0
    };
    base / completeness.max(0.5)
}

/// Overall confidence from data completeness (Scala `calculateConfidence`).
fn confidence_from_completeness(snps_with_data: usize, total_snps: usize) -> f64 {
    if total_snps == 0 {
        return 0.0;
    }
    let completeness = snps_with_data as f64 / total_snps as f64;
    let adjusted = if completeness < 0.5 {
        completeness * 0.5
    } else {
        0.25 + completeness * 0.75
    };
    adjusted.min(1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sg(contig: &str, pos: i64, dosage: i32) -> SiteGenotype {
        SiteGenotype {
            name: format!("{contig}:{pos}"),
            contig: contig.to_string(),
            position: pos,
            reference_allele: "A".to_string(),
            alternate_allele: "G".to_string(),
            ploidy: 2,
            dosage,
            gq: 50,
            depth: 30,
            ref_depth: 0,
            alt_depth: 0,
            pls: vec![0, 50, 99],
        }
    }

    /// Two populations, A (alt-rich) and B (alt-poor). A sample homozygous-alt at every site
    /// must score overwhelmingly as A.
    #[test]
    fn af_likelihood_picks_the_matching_population() {
        let sites: Vec<PanelSite> = (1..=20)
            .map(|pos| PanelSite {
                contig: "chr1".to_string(),
                position: pos,
                reference_allele: 'A',
                alternate_allele: 'G',
                freqs: vec![0.95, 0.05], // [A, B]
            })
            .collect();
        let panel = AncestryPanel {
            build: "test".to_string(),
            populations: vec!["A".to_string(), "B".to_string()],
            sites,
        };
        let genotypes: Vec<SiteGenotype> = (1..=20).map(|p| sg("chr1", p, 2)).collect();

        let result = estimate_by_allele_frequency(&genotypes, &panel, "test-ref");
        let top = result.primary().unwrap();
        assert_eq!(top.population_code, "A");
        assert!(top.percentage > 99.0, "A% = {}", top.percentage);
        assert_eq!(result.snps_with_genotype, 20);
        assert_eq!(result.snps_analyzed, 20);
    }

    #[test]
    fn missing_genotypes_are_dropped_from_completeness() {
        let sites: Vec<PanelSite> = (1..=10)
            .map(|pos| PanelSite {
                contig: "chr1".to_string(),
                position: pos,
                reference_allele: 'A',
                alternate_allele: 'G',
                freqs: vec![0.9, 0.1],
            })
            .collect();
        let panel = AncestryPanel { build: "t".into(), populations: vec!["A".into(), "B".into()], sites };
        // Half the sites are no-calls (dosage -1).
        let genotypes: Vec<SiteGenotype> =
            (1..=10).map(|p| sg("chr1", p, if p <= 5 { 2 } else { -1 })).collect();

        let result = estimate_by_allele_frequency(&genotypes, &panel, "t");
        assert_eq!(result.snps_with_genotype, 5);
        assert_eq!(result.snps_missing, 5);
        assert!(result.confidence_level < 1.0);
    }

    #[test]
    fn panel_roundtrips_through_bincode() {
        let panel = AncestryPanel {
            build: "chm13v2.0".to_string(),
            populations: vec!["AFR".into(), "EUR".into()],
            sites: vec![PanelSite {
                contig: "chr1".into(),
                position: 12345,
                reference_allele: 'C',
                alternate_allele: 'T',
                freqs: vec![0.3, 0.7],
            }],
        };
        let bytes = panel.to_bytes().unwrap();
        let back = AncestryPanel::from_bytes(&bytes).unwrap();
        assert_eq!(panel, back);
    }

    /// A 1-component PCA where the loading is +1 at every site and the panel mean is 1.0
    /// (a het reference): a hom-alt sample projects to +n_sites, a hom-ref sample to −n_sites.
    #[test]
    fn project_pca_centres_and_accumulates() {
        let sites: Vec<(String, i64)> = (1..=4).map(|p| ("chr1".to_string(), p)).collect();
        let pca = PcaLoadings {
            build: "t".into(),
            sites: sites.clone(),
            means: vec![1.0; 4],
            n_components: 1,
            loadings: vec![1.0; 4],
            populations: vec!["LO".into(), "HI".into()],
            centroids: vec![-4.0, 4.0], // LO at -4, HI at +4 on PC1
            variances: vec![1.0, 1.0],
        };
        let hom_alt: Vec<SiteGenotype> = (1..=4).map(|p| sg("chr1", p, 2)).collect();
        let coords = project_pca(&hom_alt, &pca);
        assert_eq!(coords.len(), 1);
        assert!((coords[0] - 4.0).abs() < 1e-9, "coord = {}", coords[0]); // (2-1)*1 × 4 sites

        // …and the Gaussian classifier places it with the HI population.
        let probs = classify_pca(&coords, &pca);
        let hi = probs.iter().find(|(c, _)| c == "HI").unwrap().1;
        assert!(hi > 0.99, "HI prob = {hi}");
    }

    #[test]
    fn pca_loadings_roundtrip_and_accessors() {
        let pca = PcaLoadings {
            build: "chm13v2.0".into(),
            sites: vec![("chr1".into(), 10), ("chr2".into(), 20)],
            means: vec![0.5, 1.5],
            n_components: 2,
            loadings: vec![0.1, 0.2, 0.3, 0.4],
            populations: vec!["AFR".into(), "EUR".into()],
            centroids: vec![1.0, 2.0, 3.0, 4.0],
            variances: vec![0.5, 0.5, 0.5, 0.5],
        };
        let back = PcaLoadings::from_bytes(&pca.to_bytes().unwrap()).unwrap();
        assert_eq!(pca, back);
        assert_eq!(back.loading(1, 0), 0.3);
        assert_eq!(back.centroid(1), &[3.0, 4.0]);
    }
}
