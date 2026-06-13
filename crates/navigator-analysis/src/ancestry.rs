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
    population_color, population_name, population_super, AncestryResult, AncestrySegment,
    ConfidenceInterval, PopulationComponent, SuperPopulationSummary,
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
/// 0 (mean-imputed), then the projection is rescaled by `total_sites / sites_used` so a sample
/// with missing genotypes isn't shrunk toward the origin (which would pull it off its true
/// cluster). Returns the sample's coordinate in each principal component.
pub fn project_pca(genotypes: &[SiteGenotype], pca: &PcaLoadings) -> Vec<f64> {
    let dosage: HashMap<(&str, i64), i32> = genotypes
        .iter()
        .filter(|g| g.dosage >= 0)
        .map(|g| ((g.contig.as_str(), g.position), g.dosage))
        .collect();

    let mut coords = vec![0.0f64; pca.n_components];
    let mut used = 0usize;
    for (i, (contig, pos)) in pca.sites.iter().enumerate() {
        let centered = match dosage.get(&(contig.as_str(), *pos)) {
            Some(&d) => d as f64 - pca.means[i] as f64,
            None => continue,
        };
        used += 1;
        for (c, coord) in coords.iter_mut().enumerate() {
            *coord += centered * pca.loading(i, c) as f64;
        }
    }
    // Un-shrink: reference coords were built from all sites; scale up for the missing fraction.
    if used > 0 {
        let scale = pca.sites.len() as f64 / used as f64;
        for coord in &mut coords {
            *coord *= scale;
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

/// Parameters for [`paint_local_ancestry`].
#[derive(Debug, Clone)]
pub struct PaintParams {
    /// Per-bp ancestry-switch rate (segment-length knob): switch prob over distance `d` bp is
    /// `1 - exp(-d·rate)`. Smaller → longer segments. Default ≈ one switch per 20 Mb.
    pub rate: f64,
    /// Runs shorter than this many markers are merged into the neighbouring segment.
    pub min_segment_sites: usize,
}

impl Default for PaintParams {
    fn default() -> Self {
        Self { rate: 1.0 / 20_000_000.0, min_segment_sites: 5 }
    }
}

fn logsumexp(xs: &[f64]) -> f64 {
    let m = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if m == f64::NEG_INFINITY {
        return m;
    }
    m + xs.iter().map(|&x| (x - m).exp()).sum::<f64>().ln()
}

/// Diploid genotype log-likelihood under alt-allele frequency `f` (HWE binomial).
fn geno_loglik(g: i32, f: f64) -> f64 {
    let f = f.clamp(1e-4, 1.0 - 1e-4);
    match g {
        0 => 2.0 * (1.0 - f).ln(),
        1 => (2.0 * f * (1.0 - f)).ln(),
        2 => 2.0 * f.ln(),
        _ => 0.0, // missing → uniform (no information)
    }
}

/// Paint each chromosome with local ancestry: an HMM over the panel sites whose hidden states are
/// the super-populations, emissions are the diploid genotype likelihood under each population's
/// allele frequency, and transitions penalise ancestry switches by physical distance. Viterbi
/// gives the segment path; forward-backward gives per-site posteriors (segment confidence).
///
/// `prior` is the genome-wide composition `(population_code, weight)` (rolled to super-populations
/// here) — the HMM's stationary/switch distribution, anchoring the painting to the global estimate.
/// Diploid, single-ancestry-per-locus (not per-haplotype): an unadmixed sample paints one colour.
pub fn paint_local_ancestry(
    genotypes: &[SiteGenotype],
    panel: &AncestryPanel,
    prior: &[(String, f64)],
    params: &PaintParams,
) -> Vec<AncestrySegment> {
    // Super-population states present in the panel (stable order), and each panel pop's state.
    let pop_state: Vec<String> =
        panel.populations.iter().map(|c| population_super(c).unwrap_or(c).to_string()).collect();
    let mut states: Vec<String> = Vec::new();
    for s in &pop_state {
        if !states.contains(s) {
            states.push(s.clone());
        }
    }
    let k = states.len();
    if k == 0 {
        return Vec::new();
    }
    let state_idx = |s: &str| states.iter().position(|x| x == s);

    // Prior π over states (roll the global composition up to super-pops; normalize; uniform fallback).
    let mut pi = vec![0.0f64; k];
    for (code, w) in prior {
        let sp = population_super(code).unwrap_or(code);
        if let Some(j) = state_idx(sp) {
            pi[j] += w.max(0.0);
        }
    }
    let tot: f64 = pi.iter().sum();
    if tot > 0.0 {
        pi.iter_mut().for_each(|p| *p /= tot);
    } else {
        pi.iter_mut().for_each(|p| *p = 1.0 / k as f64);
    }

    let dosage: HashMap<(&str, i64), i32> = genotypes
        .iter()
        .filter(|g| g.dosage >= 0)
        .map(|g| ((g.contig.as_str(), g.position), g.dosage))
        .collect();

    // Per-contig sites with a genotype: (pos, per-state super-pop AF, dosage). Sorted by pos.
    let mut by_contig: BTreeMap<String, Vec<(i64, Vec<f64>, i32)>> = BTreeMap::new();
    for site in &panel.sites {
        if site.freqs.len() != panel.populations.len() {
            continue;
        }
        let Some(&g) = dosage.get(&(site.contig.as_str(), site.position)) else { continue };
        // Mean fine-pop frequency within each super-population state.
        let mut sum = vec![0.0f64; k];
        let mut cnt = vec![0usize; k];
        for (p, &f) in site.freqs.iter().enumerate() {
            if let Some(j) = state_idx(&pop_state[p]) {
                sum[j] += f as f64;
                cnt[j] += 1;
            }
        }
        let af: Vec<f64> = (0..k).map(|j| if cnt[j] > 0 { sum[j] / cnt[j] as f64 } else { 0.5 }).collect();
        by_contig.entry(site.contig.clone()).or_default().push((site.position, af, g));
    }

    let mut segments = Vec::new();
    for (contig, mut sites) in by_contig {
        sites.sort_by_key(|s| s.0);
        if sites.is_empty() {
            continue;
        }
        let path = viterbi(&sites, &pi, params.rate);
        let gamma = posteriors(&sites, &pi, params.rate, k);
        segments.extend(collapse_segments(&contig, &sites, &path, &gamma, &states, params.min_segment_sites));
    }
    segments
}

/// Log transition prob `a → b` given switch probability `sw` and prior `pi`:
/// stay with `1-sw`, else jump to `b` with `pi[b]`.
fn ln_trans(a: usize, b: usize, sw: f64, pi: &[f64]) -> f64 {
    let p = if a == b { (1.0 - sw) + sw * pi[b] } else { sw * pi[b] };
    p.max(1e-300).ln()
}

fn switch_prob(d: i64, rate: f64) -> f64 {
    (1.0 - (-(d.max(0) as f64) * rate).exp()).clamp(0.0, 0.999)
}

/// Viterbi MAP state path over the sites.
fn viterbi(sites: &[(i64, Vec<f64>, i32)], pi: &[f64], rate: f64) -> Vec<usize> {
    let k = pi.len();
    let n = sites.len();
    let mut v = vec![vec![f64::NEG_INFINITY; k]; n];
    let mut bp = vec![vec![0usize; k]; n];
    for s in 0..k {
        v[0][s] = pi[s].max(1e-300).ln() + geno_loglik(sites[0].2, sites[0].1[s]);
    }
    for i in 1..n {
        let sw = switch_prob(sites[i].0 - sites[i - 1].0, rate);
        for s in 0..k {
            let (mut best, mut arg) = (f64::NEG_INFINITY, 0usize);
            for (a, &va) in v[i - 1].iter().enumerate() {
                let val = va + ln_trans(a, s, sw, pi);
                if val > best {
                    best = val;
                    arg = a;
                }
            }
            v[i][s] = best + geno_loglik(sites[i].2, sites[i].1[s]);
            bp[i][s] = arg;
        }
    }
    let mut last = (0..k).max_by(|&a, &b| v[n - 1][a].total_cmp(&v[n - 1][b])).unwrap_or(0);
    let mut path = vec![0usize; n];
    path[n - 1] = last;
    for i in (1..n).rev() {
        last = bp[i][last];
        path[i - 1] = last;
    }
    path
}

/// Per-site posterior of each state (forward-backward, returns γ[i][s]).
fn posteriors(sites: &[(i64, Vec<f64>, i32)], pi: &[f64], rate: f64, k: usize) -> Vec<Vec<f64>> {
    let n = sites.len();
    let mut fwd = vec![vec![0.0f64; k]; n];
    let mut bwd = vec![vec![0.0f64; k]; n];
    for s in 0..k {
        fwd[0][s] = pi[s].max(1e-300).ln() + geno_loglik(sites[0].2, sites[0].1[s]);
    }
    for i in 1..n {
        let sw = switch_prob(sites[i].0 - sites[i - 1].0, rate);
        for s in 0..k {
            let terms: Vec<f64> = (0..k).map(|a| fwd[i - 1][a] + ln_trans(a, s, sw, pi)).collect();
            fwd[i][s] = logsumexp(&terms) + geno_loglik(sites[i].2, sites[i].1[s]);
        }
    }
    // bwd[n-1] is already 0 (log-space) from initialization.
    for i in (0..n - 1).rev() {
        let sw = switch_prob(sites[i + 1].0 - sites[i].0, rate);
        for s in 0..k {
            let terms: Vec<f64> = (0..k)
                .map(|b| ln_trans(s, b, sw, pi) + geno_loglik(sites[i + 1].2, sites[i + 1].1[b]) + bwd[i + 1][b])
                .collect();
            bwd[i][s] = logsumexp(&terms);
        }
    }
    let mut gamma = vec![vec![0.0f64; k]; n];
    for i in 0..n {
        let unn: Vec<f64> = (0..k).map(|s| fwd[i][s] + bwd[i][s]).collect();
        let z = logsumexp(&unn);
        for s in 0..k {
            gamma[i][s] = (unn[s] - z).exp();
        }
    }
    gamma
}

/// Collapse the Viterbi path into segments, merging runs shorter than `min_sites` into the
/// previous segment (keeping its ancestry).
fn collapse_segments(
    contig: &str,
    sites: &[(i64, Vec<f64>, i32)],
    path: &[usize],
    gamma: &[Vec<f64>],
    states: &[String],
    min_sites: usize,
) -> Vec<AncestrySegment> {
    // Runs of equal state: (state, first_idx, last_idx).
    let mut runs: Vec<(usize, usize, usize)> = Vec::new();
    for (i, &s) in path.iter().enumerate() {
        match runs.last_mut() {
            Some(r) if r.0 == s => r.2 = i,
            _ => runs.push((s, i, i)),
        }
    }
    // Merge short runs into the previous run.
    let mut merged: Vec<(usize, usize, usize)> = Vec::new();
    for r in runs {
        if (r.2 - r.1 + 1) < min_sites {
            if let Some(prev) = merged.last_mut() {
                prev.2 = r.2;
                continue;
            }
        }
        merged.push(r);
    }
    merged
        .into_iter()
        .map(|(s, lo, hi)| {
            let post: f64 = (lo..=hi).map(|i| gamma[i][s]).sum::<f64>() / (hi - lo + 1) as f64;
            AncestrySegment {
                contig: contig.to_string(),
                start: sites[lo].0,
                end: sites[hi].0,
                population_code: states[s].clone(),
                posterior: post,
            }
        })
        .collect()
}

/// Reverse-complement a single base (for strand reconciliation; non-ACGT passes through).
fn revcomp_base(b: char) -> char {
    match b.to_ascii_uppercase() {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        other => other,
    }
}

/// Alt-allele dosage (0/1/2) for a chip diploid call `(a1,a2)` against a panel site's
/// `ref_allele`/`alt_allele`. When the call's alleles don't both lie in `{ref,alt}`, retry once on
/// the **reverse-complemented** call (the array reported the other strand); `None` if it still
/// doesn't match (no-call / multi-allelic mismatch). The minimal strand-flip logic chip→panel needs.
pub fn dosage_from_alleles(a1: char, a2: char, ref_allele: char, alt_allele: char) -> Option<i32> {
    let (r, alt) = (ref_allele.to_ascii_uppercase(), alt_allele.to_ascii_uppercase());
    let count = |x: char, y: char| -> Option<i32> {
        let (x, y) = (x.to_ascii_uppercase(), y.to_ascii_uppercase());
        let ok = |b: char| b == r || b == alt;
        (ok(x) && ok(y)).then(|| (x == alt) as i32 + (y == alt) as i32)
    };
    count(a1, a2).or_else(|| count(revcomp_base(a1), revcomp_base(a2)))
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
        "AF_LIKELIHOOD",
        "aims",
        panel.sites.len(),
        snps_with_data,
        &probs,
        confidence,
        reference_version,
    )
}

/// Estimate the sample's **admixture proportions** over the panel populations by supervised
/// ADMIXTURE: the reference allele frequencies `P` (the panel) are fixed and we estimate the
/// mixture vector `Q` (on the simplex, summing to 1) that maximizes the genotype likelihood
/// `∏_j P(g_j | f_j)`, where the mixed alt-allele frequency at site `j` is `f_j = Σ_k q_k·p_{k,j}`
/// and `P(g_j|f_j)` is the diploid binomial under HWE.
///
/// Fitted by the frappe/ADMIXTURE EM: each allele copy has a latent source population; the E-step
/// is its posterior given ref/alt, the M-step re-estimates `q_k` as the mean posterior. Unlike
/// [`estimate_by_allele_frequency`] (which picks the single best-fitting population), this yields
/// a 100%-summing composition — the shape of a consumer ancestry report.
pub fn estimate_admixture(
    genotypes: &[SiteGenotype],
    panel: &AncestryPanel,
    reference_version: &str,
) -> AncestryResult {
    let dosage: HashMap<(&str, i64), i32> = genotypes
        .iter()
        .filter(|g| g.dosage >= 0)
        .map(|g| ((g.contig.as_str(), g.position), g.dosage))
        .collect();

    let k = panel.populations.len();
    // Informative sites: (dosage 0/1/2, clamped per-pop alt frequencies).
    let sites: Vec<(f64, Vec<f64>)> = panel
        .sites
        .iter()
        .filter(|s| s.freqs.len() == k)
        .filter_map(|s| {
            dosage.get(&(s.contig.as_str(), s.position)).map(|&g| {
                let f: Vec<f64> = s.freqs.iter().map(|&p| (p as f64).clamp(0.001, 0.999)).collect();
                (g as f64, f)
            })
        })
        .collect();
    let snps_with_data = sites.len();

    let mut q = vec![1.0 / k.max(1) as f64; k];
    if snps_with_data > 0 {
        // EM to convergence (monotone in the likelihood); cheap — O(sites·k) per iteration.
        for _ in 0..500 {
            let mut acc = vec![0.0f64; k];
            for (g, freqs) in &sites {
                let f: f64 = (0..k).map(|i| q[i] * freqs[i]).sum::<f64>().clamp(1e-9, 1.0 - 1e-9);
                let alt = *g; // expected alt allele copies
                let refc = 2.0 - g; // ref allele copies
                for i in 0..k {
                    acc[i] += alt * (q[i] * freqs[i] / f)
                        + refc * (q[i] * (1.0 - freqs[i]) / (1.0 - f));
                }
            }
            let total: f64 = acc.iter().sum(); // == 2·snps_with_data
            let mut max_delta = 0.0f64;
            if total > 0.0 {
                for i in 0..k {
                    let new = acc[i] / total;
                    max_delta = max_delta.max((new - q[i]).abs());
                    q[i] = new;
                }
            }
            if max_delta < 1e-7 {
                break;
            }
        }
    }

    let probs: Vec<(String, f64)> = panel.populations.iter().cloned().zip(q).collect();
    let confidence = confidence_from_completeness(snps_with_data, panel.sites.len());
    from_probabilities("ADMIXTURE", "genome-wide", panel.sites.len(), snps_with_data, &probs, confidence, reference_version)
}

/// Estimate ancestry by **PCA projection + a diagonal-covariance Gaussian mixture**: project the
/// sample onto the reference PCA space ([`project_pca`]) and assign per-population responsibilities
/// ([`classify_pca`]) — the `PCA_PROJECTION_GMM` method. Unlike the allele-frequency estimators,
/// the composition comes entirely from the sample's position in PC space relative to each
/// population's centroid/variance, so the `pca` asset's `populations` define the components
/// (modern super-pops, or ancient Steppe/EEF/WHG when an ancient asset is supplied). The projected
/// coordinates are attached for the scatter plot.
pub fn estimate_pca_gmm(
    genotypes: &[SiteGenotype],
    pca: &PcaLoadings,
    reference_version: &str,
) -> AncestryResult {
    let coords = project_pca(genotypes, pca);
    let probs = classify_pca(&coords, pca);

    // Coverage: how many of the PCA asset's sites this sample has a genotype for.
    let genotyped: std::collections::HashSet<(&str, i64)> = genotypes
        .iter()
        .filter(|g| g.dosage >= 0)
        .map(|g| (g.contig.as_str(), g.position))
        .collect();
    let snps_with_genotype = pca
        .sites
        .iter()
        .filter(|(c, p)| genotyped.contains(&(c.as_str(), *p)))
        .count();
    let confidence = confidence_from_completeness(snps_with_genotype, pca.sites.len());

    let mut result = from_probabilities(
        "PCA_PROJECTION_GMM",
        "genome-wide",
        pca.sites.len(),
        snps_with_genotype,
        &probs,
        confidence,
        reference_version,
    );
    result.pca_coordinates = Some(coords);
    result
}

/// Squared Euclidean distance between two equal-length coordinate vectors.
fn euclid_sq(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(&x, &y)| (x - y) * (x - y)).sum()
}

/// Non-negative, sum-to-one mixture of `sources` (each a coordinate vector) that best
/// reconstructs `target` in the Euclidean (least-squares) sense — the nMonte/Vahaduo
/// "distance" model. Returns `(weights, fit_distance)` where `fit_distance` is the residual
/// `‖target − Σ wᵢ·sourceᵢ‖`.
///
/// Solved by **Frank–Wolfe** (conditional gradient) on the probability simplex: start at the
/// nearest single source, then repeatedly move toward the vertex (population) that most reduces
/// the residual, with an exact line search. Deterministic, monotone, and projection-free — it
/// stays on the simplex by construction, so weights are always valid proportions.
fn nmonte_fit(target: &[f64], sources: &[Vec<f64>]) -> (Vec<f64>, f64) {
    let k = sources.len();
    if k == 0 {
        return (Vec::new(), f64::INFINITY);
    }
    let dim = target.len();

    // Initialize at the single nearest source vertex.
    let (start, _) = sources
        .iter()
        .enumerate()
        .map(|(i, s)| (i, euclid_sq(target, s)))
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
        .unwrap();
    let mut w = vec![0.0f64; k];
    w[start] = 1.0;
    let mut mix = sources[start].clone(); // current mixture Σ wᵢ·sourceᵢ

    for _ in 0..1000 {
        // residual r = mix − target; gradient of ½‖r‖² wrt wⱼ is r·sourceⱼ.
        let resid: Vec<f64> = (0..dim).map(|c| mix[c] - target[c]).collect();
        // Pick the vertex j minimizing the linear approximation (steepest descent on simplex).
        let (j, _) = sources
            .iter()
            .enumerate()
            .map(|(i, s)| (i, resid.iter().zip(s).map(|(&r, &x)| r * x).sum::<f64>()))
            .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
            .unwrap();
        // Direction d = sourceⱼ − mix; exact step γ = clamp(−r·d / ‖d‖², 0, 1).
        let dvec: Vec<f64> = (0..dim).map(|c| sources[j][c] - mix[c]).collect();
        let num = -resid.iter().zip(&dvec).map(|(&r, &d)| r * d).sum::<f64>();
        let den = dvec.iter().map(|&d| d * d).sum::<f64>();
        if den < 1e-12 {
            break;
        }
        let gamma = (num / den).clamp(0.0, 1.0);
        if gamma < 1e-9 {
            break; // converged (no improving direction)
        }
        for x in w.iter_mut() {
            *x *= 1.0 - gamma;
        }
        w[j] += gamma;
        for c in 0..dim {
            mix[c] += gamma * dvec[c];
        }
    }

    (w, euclid_sq(&mix, target).sqrt())
}

/// Estimate ancestry by **PCA projection + a distance-minimizing mixture fit** (the
/// nMonte/G25-style model) — the `G25_NMONTE` method. Projects the sample into PC space
/// ([`project_pca`]) and fits the non-negative, sum-to-one mixture of the reference populations'
/// centroids that best reconstructs that point ([`nmonte_fit`]). Unlike the GMM classifier (which
/// assigns to the nearest cluster), this *decomposes* an admixed sample into source proportions
/// and reports the fit residual as a quality score (`fit_distance`; lower is better). The `pca`
/// asset's `populations` are the source library, so a richer/global asset yields wider admixtures.
pub fn estimate_nmonte(
    genotypes: &[SiteGenotype],
    pca: &PcaLoadings,
    reference_version: &str,
) -> AncestryResult {
    let coords = project_pca(genotypes, pca);
    let sources: Vec<Vec<f64>> = (0..pca.populations.len())
        .map(|p| pca.centroid(p).iter().map(|&x| x as f64).collect())
        .collect();
    let (weights, distance) = nmonte_fit(&coords, &sources);
    let probs: Vec<(String, f64)> = pca.populations.iter().cloned().zip(weights).collect();

    // Coverage: how many of the PCA asset's sites this sample has a genotype for.
    let genotyped: std::collections::HashSet<(&str, i64)> = genotypes
        .iter()
        .filter(|g| g.dosage >= 0)
        .map(|g| (g.contig.as_str(), g.position))
        .collect();
    let snps_with_genotype = pca
        .sites
        .iter()
        .filter(|(c, p)| genotyped.contains(&(c.as_str(), *p)))
        .count();
    let confidence = confidence_from_completeness(snps_with_genotype, pca.sites.len());

    let mut result = from_probabilities(
        "G25_NMONTE",
        "genome-wide",
        pca.sites.len(),
        snps_with_genotype,
        &probs,
        confidence,
        reference_version,
    );
    result.pca_coordinates = Some(coords);
    result.fit_distance = Some(distance);
    result
}

/// Build an [`AncestryResult`] from raw per-population probabilities (need not be normalized).
/// With the phase-1 super-population panel each component *is* a super-population, so the
/// super-population summary is 1:1 with the components.
fn from_probabilities(
    method: &str,
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
        method: method.to_string(),
        panel_type: panel_type.to_string(),
        snps_analyzed,
        snps_with_genotype,
        snps_missing: snps_analyzed.saturating_sub(snps_with_genotype),
        components,
        super_population_summary,
        confidence_level,
        fit_distance: None,
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

    #[test]
    fn dosage_from_alleles_counts_alt_with_strand_flip() {
        // ref=A alt=G: hom-ref, het, hom-alt.
        assert_eq!(dosage_from_alleles('A', 'A', 'A', 'G'), Some(0));
        assert_eq!(dosage_from_alleles('A', 'G', 'A', 'G'), Some(1));
        assert_eq!(dosage_from_alleles('G', 'G', 'A', 'G'), Some(2));
        // Opposite strand (chip reported C/T for an A/G site) → rev-comp matches: T→A, C→G.
        assert_eq!(dosage_from_alleles('T', 'C', 'A', 'G'), Some(1));
        assert_eq!(dosage_from_alleles('C', 'C', 'A', 'G'), Some(2));
        // A genuine mismatch (neither strand fits) → no-call.
        assert_eq!(dosage_from_alleles('A', 'C', 'A', 'G'), None);
    }

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

    /// estimate_pca_gmm: stamps the method, builds a 100%-summing composition from the GMM
    /// responsibilities, attaches the projected coordinates, and counts covered sites.
    #[test]
    fn pca_gmm_estimate_labels_and_composes() {
        let sites: Vec<(String, i64)> = (1..=4).map(|p| ("chr1".to_string(), p)).collect();
        let pca = PcaLoadings {
            build: "t".into(),
            sites,
            means: vec![1.0; 4],
            n_components: 1,
            loadings: vec![1.0; 4],
            populations: vec!["LO".into(), "HI".into()],
            centroids: vec![-4.0, 4.0],
            variances: vec![1.0, 1.0],
        };
        let hom_alt: Vec<SiteGenotype> = (1..=4).map(|p| sg("chr1", p, 2)).collect();
        let r = estimate_pca_gmm(&hom_alt, &pca, "t");

        assert_eq!(r.method, "PCA_PROJECTION_GMM");
        assert_eq!(r.panel_type, "genome-wide");
        assert_eq!(r.snps_with_genotype, 4); // all PCA sites covered
        assert!(r.pca_coordinates.is_some());
        let hi = r.components.iter().find(|c| c.population_code == "HI").unwrap();
        assert!(hi.percentage > 99.0, "HI% = {}", hi.percentage);
        let sum: f64 = r.components.iter().map(|c| c.percentage).sum();
        assert!((sum - 100.0).abs() < 1e-6, "sum = {sum}");
    }

    /// nMonte fit: a target on a source vertex resolves to ~100% that source (distance ~0);
    /// a target at the midpoint of two sources resolves to ~50/50.
    #[test]
    fn nmonte_fit_recovers_mixtures() {
        let a = vec![0.0, 0.0];
        let b = vec![10.0, 0.0];
        let c = vec![0.0, 10.0];
        let sources = vec![a.clone(), b.clone(), c.clone()];

        // On a vertex → that source, ~zero distance.
        let (w, d) = nmonte_fit(&b, &sources);
        assert!(w[1] > 0.999, "w_b = {}", w[1]);
        assert!(d < 1e-6, "distance = {d}");

        // Midpoint of A and B → ~50/50, ~zero distance (it's in the convex hull).
        let mid = vec![5.0, 0.0];
        let (w, d) = nmonte_fit(&mid, &sources);
        assert!((w[0] - 0.5).abs() < 1e-3 && (w[1] - 0.5).abs() < 1e-3, "w = {w:?}");
        assert!(d < 1e-6, "distance = {d}");

        // A point outside the hull → best projection, with a non-zero residual distance.
        let outside = vec![-3.0, -3.0];
        let (_w, d) = nmonte_fit(&outside, &sources);
        assert!(d > 1.0, "distance = {d}");
    }

    /// estimate_nmonte: labels the method, attaches coords + a fit distance, and (for a sample
    /// projecting onto a population centroid) puts ~all the weight there.
    #[test]
    fn nmonte_estimate_labels_and_fits() {
        let sites: Vec<(String, i64)> = (1..=4).map(|p| ("chr1".to_string(), p)).collect();
        let pca = PcaLoadings {
            build: "t".into(),
            sites,
            means: vec![1.0; 4],
            n_components: 1,
            loadings: vec![1.0; 4],
            populations: vec!["LO".into(), "HI".into()],
            centroids: vec![-4.0, 4.0], // LO at -4, HI at +4 on PC1
            variances: vec![1.0, 1.0],
        };
        // hom-alt projects to +4 → exactly the HI centroid.
        let hom_alt: Vec<SiteGenotype> = (1..=4).map(|p| sg("chr1", p, 2)).collect();
        let r = estimate_nmonte(&hom_alt, &pca, "t");

        assert_eq!(r.method, "G25_NMONTE");
        assert!(r.pca_coordinates.is_some());
        let dist = r.fit_distance.expect("fit_distance set");
        assert!(dist < 1e-6, "distance = {dist}");
        let hi = r.components.iter().find(|c| c.population_code == "HI").unwrap();
        assert!(hi.percentage > 99.0, "HI% = {}", hi.percentage);
    }

    /// Supervised admixture: a sample homozygous-alt where pop A is alt-rich and hom-ref where A
    /// is alt-poor must resolve to ~100% A; the proportions sum to 100%.
    #[test]
    fn admixture_resolves_pure_population() {
        let sites: Vec<PanelSite> = (1..=40)
            .map(|pos| PanelSite {
                contig: "chr1".to_string(),
                position: pos,
                reference_allele: 'A',
                alternate_allele: 'G',
                freqs: if pos % 2 == 0 { vec![0.95, 0.05] } else { vec![0.05, 0.95] },
            })
            .collect();
        let panel = AncestryPanel { build: "t".into(), populations: vec!["A".into(), "B".into()], sites };
        // Genotype to match A: hom-alt (2) at A-rich even sites, hom-ref (0) at A-poor odd sites.
        let genos: Vec<SiteGenotype> =
            (1..=40).map(|p| sg("chr1", p, if p % 2 == 0 { 2 } else { 0 })).collect();

        let r = estimate_admixture(&genos, &panel, "t");
        let a = r.components.iter().find(|c| c.population_code == "A").unwrap();
        assert!(a.percentage > 95.0, "A% = {}", a.percentage);
        let sum: f64 = r.components.iter().map(|c| c.percentage).sum();
        assert!((sum - 100.0).abs() < 1e-6, "sum = {sum}");
    }

    /// A sample that is genotype-wise a 50/50 blend of two divergent populations yields roughly
    /// balanced admixture proportions.
    #[test]
    fn admixture_detects_a_mixture() {
        // Pop A fixed alt, pop B fixed ref. A 50/50 mix → every site heterozygous.
        let sites: Vec<PanelSite> = (1..=60)
            .map(|pos| PanelSite {
                contig: "chr1".to_string(),
                position: pos,
                reference_allele: 'A',
                alternate_allele: 'G',
                freqs: vec![0.99, 0.01],
            })
            .collect();
        let panel = AncestryPanel { build: "t".into(), populations: vec!["A".into(), "B".into()], sites };
        let genos: Vec<SiteGenotype> = (1..=60).map(|p| sg("chr1", p, 1)).collect(); // all het
        let r = estimate_admixture(&genos, &panel, "t");
        let a = r.components.iter().find(|c| c.population_code == "A").unwrap().percentage;
        assert!((40.0..=60.0).contains(&a), "A% = {a} (expected ~50)");
    }

    /// A chromosome whose first half's genotypes match pop A and second half match pop B should
    /// paint two segments (A then B) with a single switch; an all-A chromosome paints one segment.
    #[test]
    fn painting_finds_a_switch() {
        // Two pops, A alt-rich / B alt-poor at every site (so genotype discriminates).
        let n = 80;
        let sites: Vec<PanelSite> = (0..n)
            .map(|i| PanelSite {
                contig: "chr1".to_string(),
                position: 1 + i as i64 * 1_000_000, // 1 Mb spacing
                reference_allele: 'A',
                alternate_allele: 'G',
                freqs: vec![0.95, 0.05],
            })
            .collect();
        let panel = AncestryPanel { build: "t".into(), populations: vec!["A".into(), "B".into()], sites };
        // First half hom-alt (matches A), second half hom-ref (matches B).
        let genos: Vec<SiteGenotype> = (0..n)
            .map(|i| sg("chr1", 1 + i as i64 * 1_000_000, if i < n / 2 { 2 } else { 0 }))
            .collect();
        let prior = vec![("A".to_string(), 0.5), ("B".to_string(), 0.5)];

        let segs = paint_local_ancestry(&genos, &panel, &prior, &PaintParams::default());
        assert_eq!(segs.len(), 2, "expected one switch: {segs:?}");
        assert_eq!(segs[0].population_code, "A");
        assert_eq!(segs[1].population_code, "B");
        assert!(segs[0].end < segs[1].start);

        // All hom-alt → a single A segment.
        let all_a: Vec<SiteGenotype> = (0..n).map(|i| sg("chr1", 1 + i as i64 * 1_000_000, 2)).collect();
        let one = paint_local_ancestry(&all_a, &panel, &prior, &PaintParams::default());
        assert_eq!(one.len(), 1);
        assert_eq!(one[0].population_code, "A");
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
