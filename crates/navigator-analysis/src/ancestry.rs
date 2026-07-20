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

use nalgebra::{DMatrix, DVector};
use navigator_domain::ancestry::{
    fine_population_codes, population_color, population_name, population_super, AncestryResult, AncestrySegment,
    ConfidenceInterval, PopulationComponent, SuperPopulationSummary,
};
use serde::{Deserialize, Serialize};

use crate::caller::SiteGenotype;
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

    /// A panel restricted to `codes` (those present, in `codes` order), projecting each site's
    /// per-population frequencies down to the kept columns. Used to run a well-conditioned
    /// admixture EM over a curated subset of a large fine-frequency panel.
    pub fn subset(&self, codes: &[&str]) -> AncestryPanel {
        let keep: Vec<usize> = codes
            .iter()
            .filter_map(|c| self.populations.iter().position(|p| p == c))
            .collect();
        let populations = keep.iter().map(|&i| self.populations[i].clone()).collect();
        let sites = self
            .sites
            .iter()
            .map(|s| PanelSite {
                contig: s.contig.clone(),
                position: s.position,
                reference_allele: s.reference_allele,
                alternate_allele: s.alternate_allele,
                freqs: keep.iter().map(|&i| s.freqs.get(i).copied().unwrap_or(0.0)).collect(),
            })
            .collect();
        AncestryPanel {
            build: self.build.clone(),
            populations,
            sites,
        }
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
        Self {
            rate: 1.0 / 20_000_000.0,
            min_segment_sites: 5,
        }
    }
}

/// Diploid genotype log-likelihood when the two genome copies draw their alt allele from frequencies
/// `fa` and `fb` (one copy per ancestry): `P(0)=(1-fa)(1-fb)`, `P(1)=fa(1-fb)+(1-fa)fb`, `P(2)=fa·fb`.
/// Missing dosage → uniform. This is the proper diploid (two-copy) emission the pair-state HMM needs.
fn emit_diploid_ln(g: i32, fa: f64, fb: f64) -> f64 {
    let fa = fa.clamp(1e-4, 1.0 - 1e-4);
    let fb = fb.clamp(1e-4, 1.0 - 1e-4);
    let p = match g {
        0 => (1.0 - fa) * (1.0 - fb),
        1 => fa * (1.0 - fb) + (1.0 - fa) * fb,
        2 => fa * fb,
        _ => return 0.0, // missing → uniform
    };
    p.max(1e-300).ln()
}

/// Paint each chromosome with local ancestry: an HMM over the panel sites whose hidden states are
/// the super-populations, emissions are the diploid genotype likelihood under each population's
/// allele frequency, and transitions penalise ancestry switches by physical distance. Viterbi
/// gives the segment path; forward-backward gives per-site posteriors (segment confidence).
///
/// `prior` is the genome-wide composition `(population_code, weight)` (rolled to super-populations
/// here) — the HMM's stationary/switch distribution, anchoring the painting to the global estimate.
/// **Diploid pair-state HMM**: the hidden state is an ancestry pair (both genome copies), so a region
/// where the two copies differ (e.g. EUR/SAS) is shown, not collapsed. Output is two sorted, unphased
/// copies per chromosome (segments tagged `copy` 0/1) — not maternal/paternal (no phasing).
pub fn paint_local_ancestry(
    genotypes: &[SiteGenotype],
    panel: &AncestryPanel,
    prior: &[(String, f64)],
    params: &PaintParams,
) -> Vec<AncestrySegment> {
    // Super-population states present in the panel (stable order), and each panel pop's state.
    let pop_state: Vec<String> = panel
        .populations
        .iter()
        .map(|c| population_super(c).unwrap_or(c).to_string())
        .collect();
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
        let Some(&g) = dosage.get(&(site.contig.as_str(), site.position)) else {
            continue;
        };
        // Mean fine-pop frequency within each super-population state.
        let mut sum = vec![0.0f64; k];
        let mut cnt = vec![0usize; k];
        for (p, &f) in site.freqs.iter().enumerate() {
            if let Some(j) = state_idx(&pop_state[p]) {
                sum[j] += f as f64;
                cnt[j] += 1;
            }
        }
        let af: Vec<f64> = (0..k)
            .map(|j| if cnt[j] > 0 { sum[j] / cnt[j] as f64 } else { 0.5 })
            .collect();
        by_contig
            .entry(site.contig.clone())
            .or_default()
            .push((site.position, af, g));
    }

    let mut segments = Vec::new();
    for (contig, mut sites) in by_contig {
        sites.sort_by_key(|s| s.0);
        if sites.is_empty() {
            continue;
        }
        // Diploid MAP path: one ancestry PAIR per locus, canonicalized (min,max) into two sorted,
        // coherent copies (copy 0 = lower-index ancestry, copy 1 = higher). Unphased.
        let pairs = diploid_viterbi(&sites, &pi, params.rate, k);
        let copy0: Vec<usize> = pairs.iter().map(|&(a, b)| a.min(b)).collect();
        let copy1: Vec<usize> = pairs.iter().map(|&(a, b)| a.max(b)).collect();
        segments.extend(collapse_copy(
            &contig,
            &sites,
            &copy0,
            &states,
            params.min_segment_sites,
            0,
        ));
        segments.extend(collapse_copy(
            &contig,
            &sites,
            &copy1,
            &states,
            params.min_segment_sites,
            1,
        ));
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

/// Diploid Viterbi: the MAP ancestry **pair** `(a1, a2)` per site. Hidden state = an ordered pair of
/// ancestries (the two genome copies, independent Markov chains), so transitions factorize as
/// `ln_trans(a1,b1) + ln_trans(a2,b2)` and the emission is the two-copy [`emit_diploid_ln`]. Returns
/// one `(a1, a2)` per site (state index `a1*k + a2`).
fn diploid_viterbi(sites: &[(i64, Vec<f64>, i32)], pi: &[f64], rate: f64, k: usize) -> Vec<(usize, usize)> {
    let n = sites.len();
    let ns = k * k;
    let lnpi: Vec<f64> = (0..k).map(|s| pi[s].max(1e-300).ln()).collect();
    let mut v = vec![vec![f64::NEG_INFINITY; ns]; n];
    let mut bp = vec![vec![0usize; ns]; n];
    for a1 in 0..k {
        for a2 in 0..k {
            v[0][a1 * k + a2] = lnpi[a1] + lnpi[a2] + emit_diploid_ln(sites[0].2, sites[0].1[a1], sites[0].1[a2]);
        }
    }
    for i in 1..n {
        let sw = switch_prob(sites[i].0 - sites[i - 1].0, rate);
        // Per-chain best predecessor for each target chain-state (factorized, so the pair step is
        // O(k²) not O(k⁴)): for chain value b, max over a of v_chain[a] + ln_trans(a,b).
        for b1 in 0..k {
            for b2 in 0..k {
                let (mut best, mut arg) = (f64::NEG_INFINITY, 0usize);
                for a1 in 0..k {
                    let t1 = ln_trans(a1, b1, sw, pi);
                    for a2 in 0..k {
                        let val = v[i - 1][a1 * k + a2] + t1 + ln_trans(a2, b2, sw, pi);
                        if val > best {
                            best = val;
                            arg = a1 * k + a2;
                        }
                    }
                }
                v[i][b1 * k + b2] = best + emit_diploid_ln(sites[i].2, sites[i].1[b1], sites[i].1[b2]);
                bp[i][b1 * k + b2] = arg;
            }
        }
    }
    let mut last = (0..ns)
        .max_by(|&a, &b| v[n - 1][a].total_cmp(&v[n - 1][b]))
        .unwrap_or(0);
    let mut path = vec![(0usize, 0usize); n];
    path[n - 1] = (last / k, last % k);
    for i in (1..n).rev() {
        last = bp[i][last];
        path[i - 1] = (last / k, last % k);
    }
    path
}

/// Collapse one copy's per-site ancestry path into segments, merging runs shorter than `min_sites`
/// into the previous segment (keeping its ancestry). Each segment is tagged with the `copy` index.
/// `posterior` is set to 1.0 (the MAP path; per-copy posterior shading is a future refinement).
fn collapse_copy(
    contig: &str,
    sites: &[(i64, Vec<f64>, i32)],
    path: &[usize],
    states: &[String],
    min_sites: usize,
    copy: u8,
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
        .map(|(s, lo, hi)| AncestrySegment {
            contig: contig.to_string(),
            start: sites[lo].0,
            end: sites[hi].0,
            population_code: states[s].clone(),
            posterior: 1.0,
            copy,
        })
        .collect()
}

use navigator_domain::seq::complement_base as revcomp_base;

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
                    acc[i] += alt * (q[i] * freqs[i] / f) + refc * (q[i] * (1.0 - freqs[i]) / (1.0 - f));
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
    from_probabilities(
        "ADMIXTURE",
        "genome-wide",
        panel.sites.len(),
        snps_with_data,
        &probs,
        confidence,
        reference_version,
    )
}

/// Fine-population admixture: the same supervised EM as [`estimate_admixture`], run over a curated
/// **modern subset** of a large fine-frequency panel (the `freq_global` asset carries all reference
/// populations incl. ancient; a flat 173-way EM is ill-posed, so we restrict to `modern_codes`).
/// Reuses the *same* genotypes (the fine panel shares the AIM panel's sites). The result is labeled
/// `FINE_ADMIXTURE`; its components roll up to the super-pops via the domain `population_super` map.
pub fn estimate_fine_admixture(
    genotypes: &[SiteGenotype],
    fine_panel: &AncestryPanel,
    reference_version: &str,
) -> AncestryResult {
    let subset = fine_panel.subset(&fine_population_codes());
    let mut result = estimate_admixture(genotypes, &subset, reference_version);
    result.method = "FINE_ADMIXTURE".to_string();
    result.panel_type = "fine".to_string();
    result
}

/// The `ANCIENT_ADMIXTURE` method label — deep (pre-historic) source proportions.
pub const ANCIENT_ADMIXTURE: &str = "ANCIENT_ADMIXTURE";

/// Below this many genotyped panel sites the three-way fit is too noisy to report at all.
const ANCIENT_MIN_SITES: usize = 500;

/// Dispersion above which the sample is **outside the span of the ancient sources** and we report
/// nothing. Under a correct model the dispersion is ≈1 by construction (see [`ancient_dispersion`]).
///
/// Calibrated on simulated reference individuals (`panelbuild validate-ancient`), worst case per
/// population: GBR 1.65 · CEU 1.58 · FIN 1.78 · TSI 2.38 · **IBS 3.65** ‖ CHB 13.1 · JPT 12.4 ·
/// YRI 175 · LWK 158. So 4.0 sits in the wide, empty gap between "every European individual" and
/// "the closest East Asian" — it is not a knob tuned to taste, it is the middle of a real gap.
///
/// It deliberately does **not** try to separate South Asians (PJL 3.3–4.0), who overlap the European
/// tail: no dispersion threshold can do that, which is why there is a second guard,
/// [`ANCIENT_MIN_WEST_EURASIAN`].
const ANCIENT_MAX_DISPERSION: f64 = 4.0;

/// Minimum European share (by the modern super-population admixture) for the deep three-way model to
/// apply at all.
///
/// WHG / Anatolian Farmer / Steppe is a **West-Eurasian** model: those three sources are the ones
/// that actually compose modern Europeans. It has no term for Ancestral South Indian, no term for
/// East Asian, and no term for Sub-Saharan African, so for a person who carries a lot of any of
/// those, a three-way decomposition of their *whole genome* is not an approximation — it is a
/// category error. A Punjabi fits at Steppe 67% here; their real Steppe ancestry is nearer 20–30%,
/// with the rest Iranian-Neolithic and AASI that this model simply cannot see, so it piles the
/// unexplained ancestry onto whichever source is least unlike it.
///
/// Dispersion alone cannot catch that (South Asians overlap the European tail), but the *modern*
/// estimate — which is well validated and independent of this panel — separates them cleanly. So
/// deep ancestry only runs for samples the modern model already calls predominantly European.
const ANCIENT_MIN_WEST_EURASIAN: f64 = 50.0;

/// qpAdm model-fit acceptance: report the deep breakdown only when the model is **not rejected** at
/// this tail probability (docs/design/ancient-ancestry-rebuild.md §7.14). The garbage fits the gate
/// exists to suppress reject at p ≈ 1e-13; a real British WGS/chip accepts at p ≈ 0.15–0.21.
const QPADM_MIN_P: f64 = 0.05;
/// Tolerance for the "weights are valid proportions" check — a fit that needs a source weight outside
/// `[0,1]` is the model failing, not a small numerical overshoot.
const QPADM_WEIGHT_TOL: f64 = 0.02;

/// Estimate **deep ancestral (ancient) source proportions** — the Western Hunter-Gatherer /
/// Anatolian Farmer / Steppe pastoralist decomposition — by the same supervised allele-frequency
/// admixture EM as [`estimate_admixture`], over the dedicated ancient frequency panel
/// (`ancestry_freq_ancient_<build>.bin`, built by `panelbuild ancient-panel` from the AADR).
///
/// This *replaces* an earlier PCA-centroid classifier, which was wrong twice over: it asked "which
/// ancient population **is** this sample?" (a membership posterior) where the question is "what
/// **mixture** of ancient sources is this sample?", and it ran against centroids that had been
/// shrunk on top of the modern European cloud by projection, so they carried no ancient signal at
/// all. A modern European is not a *member* of WHG; they are a *mixture*. Allele frequencies keep
/// the sources genuinely distinct (WHG↔ANF Fst ≈ 0.07) where the projected PCA did not.
///
/// `modern` is the sample's **modern** super-population admixture ([`estimate_admixture`] over the
/// super-pop panel) — an independent, already-validated estimate, used only to decide whether this
/// West-Eurasian model applies to this person at all (see [`ANCIENT_MIN_WEST_EURASIAN`]).
///
/// Returns `None` whenever the model does not apply: too few genotyped sites, too little European
/// ancestry for a WHG/ANF/Steppe decomposition to mean anything, or a fit dispersion above
/// [`ANCIENT_MAX_DISPERSION`] (the sample's ancestry lies outside the span of the three sources — a
/// Yoruba is not *any* mixture of them). Reporting nothing is the entire point: the EM will always
/// return *some* simplex vector, and presenting that vector for a sample the model cannot express is
/// precisely the failure this rebuild exists to prevent.
pub fn estimate_ancient_admixture(
    genotypes: &[SiteGenotype],
    ancient_panel: &AncestryPanel,
    modern: &AncestryResult,
    reference_version: &str,
) -> Option<AncestryResult> {
    if west_eurasian_share(modern) < ANCIENT_MIN_WEST_EURASIAN {
        return None;
    }
    let result = ancient_admixture_fit(genotypes, ancient_panel, reference_version)?;
    let dispersion = result.fit_distance.unwrap_or(f64::INFINITY);
    (dispersion.is_finite() && dispersion <= ANCIENT_MAX_DISPERSION).then_some(result)
}

/// The sample's European share (%) according to a modern super-population estimate — the scope
/// check for the deep three-way model. Reads the `EUR` rollup, so it works whether `modern` came
/// from the 5-way super-pop panel or a finer panel that rolls up to it.
///
/// `SuperPopulationSummary::super_population` carries the *display name*, not the code, so the
/// lookup goes through the catalog rather than hard-coding either spelling.
pub fn west_eurasian_share(modern: &AncestryResult) -> f64 {
    let eur = population_name("EUR");
    modern
        .super_population_summary
        .iter()
        .find(|s| s.super_population == eur || s.super_population == "EUR")
        .map_or(0.0, |s| s.percentage)
}

/// The ancient mixture fit **without** the applicability threshold: the EM result with its
/// dispersion attached as `fit_distance`, for any sample with enough genotyped sites.
///
/// [`estimate_ancient_admixture`] is this plus the [`ANCIENT_MAX_DISPERSION`] gate, and is what the
/// app calls. This variant exists so the offline validator can *report* the dispersion of samples
/// that the gate rejects — the threshold is only defensible if you can see the separation it rests
/// on. Do not use it on a user's data: its components are exactly the numbers the gate exists to
/// suppress.
pub fn ancient_admixture_fit(
    genotypes: &[SiteGenotype],
    ancient_panel: &AncestryPanel,
    reference_version: &str,
) -> Option<AncestryResult> {
    let mut result = estimate_admixture(genotypes, ancient_panel, reference_version);
    if result.snps_with_genotype < ANCIENT_MIN_SITES {
        return None;
    }

    // Recover the fitted mixture on the panel's axis order (`components` are sorted by percentage).
    let q: Vec<f64> = ancient_panel
        .populations
        .iter()
        .map(|code| {
            result
                .components
                .iter()
                .find(|c| &c.population_code == code)
                .map_or(0.0, |c| c.percentage / 100.0)
        })
        .collect();

    result.method = ANCIENT_ADMIXTURE.to_string();
    result.panel_type = "ancient".to_string();
    result.fit_distance = Some(ancient_dispersion(genotypes, ancient_panel, &q));
    Some(result)
}

/// The **app-facing deep-ancestry estimator** (docs/design/ancient-ancestry-rebuild.md §7.14): fit
/// `target = Σ wᵢ · sourcesᵢ` by qpAdm f4 and return it as an [`AncestryResult`] over the source
/// components (WHG / EEF / Steppe), or `None` when the deep model does not apply.
///
/// `sources` and `outgroups` are indices into `panel.populations` — the committed Patterson-2022
/// config: sources first (WHG/EEF/Steppe), then the sister outgroups. Gates, all of which must pass:
/// the sample is West-Eurasian (`modern`, [`ANCIENT_MIN_WEST_EURASIAN`]); enough sites were genotyped
/// ([`ANCIENT_MIN_SITES`]); the qpAdm model is **not rejected** (`p ≥` [`QPADM_MIN_P`]); and the
/// weights are feasible proportions ([`QPADM_WEIGHT_TOL`]). `None` must stay `None` all the way to the
/// UI/PDS — an inapplicable or rejected fit is reported as nothing, never as confident percentages.
///
/// This supersedes [`estimate_ancient_admixture`] (the frequency-mixture EM, which failed the
/// WGS-vs-chip stability gate); the model-fit **p-value** rides out on `fit_distance`.
pub fn estimate_qpadm_ancestry(
    genotypes: &[SiteGenotype],
    panel: &AncestryPanel,
    sources: &[usize],
    outgroups: &[usize],
    modern: &AncestryResult,
    reference_version: &str,
) -> Option<AncestryResult> {
    if west_eurasian_share(modern) < ANCIENT_MIN_WEST_EURASIAN {
        return None;
    }
    let fit = qpadm_fit(genotypes, panel, sources, outgroups, F4_BLOCK_BP)?;
    if fit.n_sites < ANCIENT_MIN_SITES || fit.p_value < QPADM_MIN_P || !fit.weights_feasible(QPADM_WEIGHT_TOL) {
        return None;
    }
    // Report the source weights as an admixture result (clamp the tiny negative overshoots the
    // feasibility gate already bounded to ≥ −tol). `from_probabilities` renormalizes.
    let probs: Vec<(String, f64)> = sources
        .iter()
        .zip(&fit.weights)
        .map(|(&i, &w)| (panel.populations[i].clone(), w.max(0.0)))
        .collect();
    let mut result = from_probabilities(
        ANCIENT_ADMIXTURE,
        "ancient",
        panel.sites.len(),
        fit.n_sites,
        &probs,
        0.9,
        reference_version,
    );
    result.fit_distance = Some(fit.p_value);
    Some(result)
}

/// Goodness of fit of a fitted ancient mixture `q`, as a **variance-ratio dispersion**.
///
/// At each genotyped site the mixture predicts an alt-allele frequency `f = Σ q_k·p_k`, so under
/// the model's own HWE assumption the observed dosage `g` has mean `2f` and variance `2f(1-f)`.
/// Averaging `(g − 2f)² / 2f(1-f)` over sites therefore gives ≈1 **when the model is right**, and
/// grows without bound as the sample's true ancestry moves outside the span of the sources — the
/// mixture is then forced to predict frequencies the genotypes keep contradicting.
///
/// It is a *ratio*, so it does not drift with panel size or the sample's coverage — which is what
/// makes it usable as a fixed applicability threshold rather than a tuned magic number.
fn ancient_dispersion(genotypes: &[SiteGenotype], panel: &AncestryPanel, q: &[f64]) -> f64 {
    let dosage: HashMap<(&str, i64), i32> = genotypes
        .iter()
        .filter(|g| g.dosage >= 0)
        .map(|g| ((g.contig.as_str(), g.position), g.dosage))
        .collect();

    let k = panel.populations.len();
    let mut sum = 0.0f64;
    let mut n = 0usize;
    for site in panel.sites.iter().filter(|s| s.freqs.len() == k) {
        let Some(&g) = dosage.get(&(site.contig.as_str(), site.position)) else {
            continue;
        };
        let f: f64 = (0..k)
            .map(|i| q[i] * (site.freqs[i] as f64).clamp(0.001, 0.999))
            .sum::<f64>()
            .clamp(1e-6, 1.0 - 1e-6);
        let expected_var = 2.0 * f * (1.0 - f);
        let resid = g as f64 - 2.0 * f;
        sum += resid * resid / expected_var;
        n += 1;
    }
    if n == 0 {
        return f64::INFINITY;
    }
    sum / n as f64
}

// ── f-statistics core (Lever 2 / qpAdm) ─────────────────────────────────────────────────────────
//
// `f4(A,B;C,D) = mean_site (a−b)(c−d)` over per-population alt-allele frequencies. It is the
// ascertainment-robust primitive qpAdm is built on (docs/design/ancient-ancestry-rebuild.md §7): a
// difference-of-differences against outgroups that cancels drift shared across the whole set, and is
// **unbiased from *pooled* frequencies** — the estimation noise in each of the four slots is
// independent, so the cross-terms vanish in expectation (no per-sample hzcorr, unlike f2/f3). The
// genotyped sample enters as its own "population" with frequency `dosage/2 ∈ {0, 0.5, 1}`.
//
// This module is the primitive: a jointly-estimated **vector** of f4 statistics with its
// block-jackknife covariance. The qpAdm GLS solve (§7.2) is assembled on top of it in a later step.

/// Genome block size (bp) for the f-statistic block jackknife. ~5 Mb ≫ the LD range, so blocks are
/// effectively independent — the assumption the jackknife variance rests on.
pub const F4_BLOCK_BP: i64 = 5_000_000;

/// A population slot in an f-statistic: either a reference population (index into
/// [`AncestryPanel::populations`]) or the genotyped sample.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Pop {
    /// Reference population `i` (its per-site frequency is `PanelSite::freqs[i]`).
    Ref(usize),
    /// The sample being placed (per-site frequency `dosage/2`).
    Target,
}

/// One f4 quartet `f4(a,b;c,d) = mean_site (a−b)(c−d)`.
#[derive(Clone, Copy, Debug)]
pub struct Quartet {
    pub a: Pop,
    pub b: Pop,
    pub c: Pop,
    pub d: Pop,
}

impl Quartet {
    pub fn new(a: Pop, b: Pop, c: Pop, d: Pop) -> Self {
        Self { a, b, c, d }
    }
}

/// A jointly-estimated vector of f4 statistics with its block-jackknife covariance — the input the
/// qpAdm GLS solve consumes.
#[derive(Clone, Debug)]
pub struct F4Estimate {
    /// Full-sample f4 point estimates, parallel to the requested quartets (ADMIXTOOLS reports the
    /// full-sample estimate as the statistic; the jackknife supplies only the covariance).
    pub values: Vec<f64>,
    /// `d×d` delete-one-block jackknife covariance of `values` (Busing et al. 1999, unequal blocks).
    pub cov: Vec<Vec<f64>>,
    /// Sites contributing (target genotyped ∧ every referenced population present).
    pub n_sites: usize,
    /// Genome blocks with ≥1 contributing site.
    pub n_blocks: usize,
}

impl F4Estimate {
    /// Standard error of statistic `i` from the jackknife covariance diagonal.
    pub fn se(&self, i: usize) -> f64 {
        self.cov.get(i).and_then(|r| r.get(i)).copied().unwrap_or(0.0).max(0.0).sqrt()
    }
}

/// Point estimate of a single `f4(a,b;c,d)` over the genotyped sites (no covariance) — a thin
/// convenience over [`f4_vector`]. `None` if fewer than two genome blocks carry a contributing site.
pub fn f4(genotypes: &[SiteGenotype], panel: &AncestryPanel, q: Quartet, block_bp: i64) -> Option<f64> {
    f4_vector(genotypes, panel, &[q], block_bp).map(|e| e.values[0])
}

/// A jointly-estimated f4 vector over `quartets`, with the Busing et al. (1999) unequal-block
/// jackknife covariance. Every statistic is measured over the **same** informative site set (target
/// genotyped ∧ all referenced populations present at the site), which is what makes the covariance a
/// valid joint covariance for a downstream GLS. `None` if any quartet references a non-existent
/// population, or fewer than two blocks carry a site (a jackknife needs ≥2 blocks).
pub fn f4_vector(
    genotypes: &[SiteGenotype],
    panel: &AncestryPanel,
    quartets: &[Quartet],
    block_bp: i64,
) -> Option<F4Estimate> {
    let d = quartets.len();
    let k = panel.populations.len();
    if d == 0 || block_bp <= 0 {
        return None;
    }
    // Reject out-of-range population indices up front — a mis-built quartet must not panic mid-scan.
    let ref_ok = |p: Pop| matches!(p, Pop::Target) || matches!(p, Pop::Ref(i) if i < k);
    if !quartets.iter().all(|q| ref_ok(q.a) && ref_ok(q.b) && ref_ok(q.c) && ref_ok(q.d)) {
        return None;
    }

    let dosage: HashMap<(&str, i64), i32> = genotypes
        .iter()
        .filter(|g| g.dosage >= 0)
        .map(|g| ((g.contig.as_str(), g.position), g.dosage))
        .collect();

    // Accumulate per genome block over the informative sites: Σ x and site count, plus the totals.
    let mut block_index: HashMap<(&str, i64), usize> = HashMap::new();
    let mut block_sum: Vec<Vec<f64>> = Vec::new();
    let mut block_n: Vec<usize> = Vec::new();
    let mut total = vec![0.0f64; d];
    let mut n_sites = 0usize;

    for site in panel.sites.iter().filter(|s| s.freqs.len() == k) {
        let Some(&g) = dosage.get(&(site.contig.as_str(), site.position)) else {
            continue;
        };
        let tf = g as f64 / 2.0;
        let freq = |p: Pop| -> f64 {
            match p {
                Pop::Ref(i) => site.freqs[i] as f64,
                Pop::Target => tf,
            }
        };
        let bkey = (site.contig.as_str(), site.position / block_bp);
        let bi = *block_index.entry(bkey).or_insert_with(|| {
            block_sum.push(vec![0.0; d]);
            block_n.push(0);
            block_sum.len() - 1
        });
        for (qi, q) in quartets.iter().enumerate() {
            let x = (freq(q.a) - freq(q.b)) * (freq(q.c) - freq(q.d));
            total[qi] += x;
            block_sum[bi][qi] += x;
        }
        block_n[bi] += 1;
        n_sites += 1;
    }

    let g = block_sum.len();
    if g < 2 || n_sites < 2 {
        return None;
    }
    let n = n_sites as f64;
    let theta: Vec<f64> = total.iter().map(|&s| s / n).collect();

    // Delete-one-block estimates θ̂_(j) and per-block weights h_j = n/m_j (Busing et al. 1999, for
    // unequal block sizes). With g ≥ 2 and every block non-empty, n − m_j ≥ 1 and h_j > 1.
    let h: Vec<f64> = block_n.iter().map(|&m| n / m as f64).collect();
    let theta_j: Vec<Vec<f64>> = (0..g)
        .map(|j| {
            let denom = n - block_n[j] as f64;
            (0..d).map(|i| (total[i] - block_sum[j][i]) / denom).collect()
        })
        .collect();

    // Bias-corrected jackknife mean θ̃_J = g·θ̂ − Σ_j (h_j−1)/h_j · θ̂_(j).
    let theta_tilde: Vec<f64> = (0..d)
        .map(|i| g as f64 * theta[i] - (0..g).map(|j| (h[j] - 1.0) / h[j] * theta_j[j][i]).sum::<f64>())
        .collect();

    // Covariance = (1/g) Σ_j d_j d_jᵀ / (h_j − 1), with d_j = h_j·θ̂ − (h_j−1)·θ̂_(j) − θ̃_J.
    let mut cov = vec![vec![0.0f64; d]; d];
    for j in 0..g {
        let dj: Vec<f64> = (0..d)
            .map(|i| h[j] * theta[i] - (h[j] - 1.0) * theta_j[j][i] - theta_tilde[i])
            .collect();
        let w = 1.0 / (h[j] - 1.0);
        for i in 0..d {
            for l in 0..d {
                cov[i][l] += w * dj[i] * dj[l];
            }
        }
    }
    for row in cov.iter_mut() {
        for c in row.iter_mut() {
            *c /= g as f64;
        }
    }

    Some(F4Estimate {
        values: theta,
        cov,
        n_sites,
        n_blocks: g,
    })
}

/// Result of a qpAdm-style f4 fit: the source weights, their standard errors, and the model-fit
/// test. See [`qpadm_fit`] and docs/design/ancient-ancestry-rebuild.md §7.2.
#[derive(Clone, Debug)]
pub struct QpAdmFit {
    /// Weights over the sources, **in the order they were passed** (sums to 1). `weights[0]` is the
    /// base source's weight, recovered as `1 − Σ others`.
    pub weights: Vec<f64>,
    /// Standard error of each weight (from the GLS normal-equations covariance).
    pub std_errors: Vec<f64>,
    /// Model-fit χ² — the minimized GLS objective (residual not explained by the source span).
    pub chi2: f64,
    /// Degrees of freedom `= (#outgroups − 1) − (#sources − 1) = #outgroups − #sources`.
    pub dof: usize,
    /// Tail probability `P(χ²_dof ≥ chi2)`. The model is **rejected** when this is small (< ~0.05):
    /// the sources can't express the target's allele-sharing with the outgroups.
    pub p_value: f64,
    pub n_sites: usize,
    pub n_blocks: usize,
}

impl QpAdmFit {
    /// Whether every weight is a valid proportion (within `tol` of `[0,1]`). qpAdm accepts a model
    /// only when it is not rejected **and** the weights are feasible.
    pub fn weights_feasible(&self, tol: f64) -> bool {
        self.weights.iter().all(|&w| w >= -tol && w <= 1.0 + tol)
    }
}

/// Residual covariance `Σ(w) = Σ_{b,b'} c_b c_{b'} Ω_block(b,b')` with `c = (1, −w₁, …, −w_{n-1})`,
/// plus a tiny ridge for invertibility. `cov` is the joint f4 covariance from [`f4_vector`], laid
/// out as `n` groups of `l` statistics (group 0 = target, groups 1.. = the non-base sources).
fn qpadm_residual_cov(cov: &[Vec<f64>], n: usize, l: usize, w: &[f64]) -> DMatrix<f64> {
    let mut c = vec![0.0f64; n];
    c[0] = 1.0;
    for i in 0..n - 1 {
        c[i + 1] = -w[i];
    }
    let mut sigma = DMatrix::<f64>::zeros(l, l);
    for b in 0..n {
        for bp in 0..n {
            let cc = c[b] * c[bp];
            if cc == 0.0 {
                continue;
            }
            for p in 0..l {
                for q in 0..l {
                    sigma[(p, q)] += cc * cov[b * l + p][bp * l + q];
                }
            }
        }
    }
    let tr: f64 = (0..l).map(|p| sigma[(p, p)].abs()).sum();
    let ridge = (1e-12 * tr / l.max(1) as f64).max(1e-18);
    for p in 0..l {
        sigma[(p, p)] += ridge;
    }
    sigma
}

/// Fit `target = Σ wᵢ · sourcesᵢ` by the qpAdm f4 method (docs §7.2). The weights are estimated from
/// the target's *allele-sharing against outgroups* — differences-of-differences that cancel drift
/// and SNP ascertainment — not from its raw frequencies, which is the property §3's frequency-EM
/// lacked. `sources` and `outgroups` are indices into `panel.populations`; the target enters through
/// `genotypes` (dosage/2 per site).
///
/// Method: for each left population `X ∈ {target, S₂..Sₙ}` form the vector
/// `φ_X = [f4(X, S₁; R₁, Rⱼ)]_{j=2..m}`; the admixture identity is `φ_target = Σ_{i≥2} wᵢ φ_{Sᵢ}`.
/// Solve the weights by iteratively-reweighted GLS against the block-jackknife covariance (the
/// residual covariance depends on the weights, since the sources are themselves estimated), then
/// read the model-fit χ²/p-value from the weighted residual.
///
/// Returns `None` when `sources.len() < 2`, `outgroups.len() < sources.len()`, the f4 vector can't be
/// formed (too few blocks/sites), or the GLS system is singular.
pub fn qpadm_fit(
    genotypes: &[SiteGenotype],
    panel: &AncestryPanel,
    sources: &[usize],
    outgroups: &[usize],
    block_bp: i64,
) -> Option<QpAdmFit> {
    let n = sources.len();
    let m = outgroups.len();
    let k = panel.populations.len();
    if n < 2 || m < n || sources.iter().chain(outgroups).any(|&i| i >= k) {
        return None;
    }
    let l = m - 1; // statistics per left population (outgroups R₂..R_m differenced vs the base R₁)
    let s1 = Pop::Ref(sources[0]);
    let r1 = outgroups[0];

    // Left populations relative to the base source: target, then S₂..Sₙ. Group order in the f4
    // vector is [target, S₂, …, Sₙ], each contributing `l` statistics over the non-base outgroups.
    let lefts: Vec<Pop> = std::iter::once(Pop::Target)
        .chain(sources[1..].iter().map(|&i| Pop::Ref(i)))
        .collect();
    let mut quartets = Vec::with_capacity(n * l);
    for &x in &lefts {
        for &rj in &outgroups[1..] {
            quartets.push(Quartet::new(x, s1, Pop::Ref(r1), Pop::Ref(rj)));
        }
    }
    let est = f4_vector(genotypes, panel, &quartets, block_bp)?;

    // y = φ_target (group 0); A[:, i] = φ_{S_{i+2}} (group i+1). Ω = est.cov, block-structured.
    let y = DVector::from_row_slice(&est.values[0..l]);
    let a = DMatrix::from_fn(l, n - 1, |p, i| est.values[(i + 1) * l + p]);

    // Iteratively-reweighted GLS: recompute Σ(w) and re-solve w = (AᵀΣ⁻¹A)⁻¹ AᵀΣ⁻¹ y until settled.
    let mut w = DVector::from_element(n - 1, 1.0 / n as f64);
    for _ in 0..100 {
        let sigma = qpadm_residual_cov(&est.cov, n, l, w.as_slice());
        let inv = sigma.try_inverse()?;
        let at_si = a.transpose() * &inv;
        let normal = (&at_si * &a).try_inverse()?;
        let new_w = &normal * (&at_si * &y);
        let delta = (&new_w - &w).amax();
        w = new_w;
        if delta < 1e-10 {
            break;
        }
    }

    // Final objective, dof, p-value, and weight SEs at the converged weights.
    let sigma_inv = qpadm_residual_cov(&est.cov, n, l, w.as_slice()).try_inverse()?;
    let r = &y - &a * &w;
    let chi2 = (r.transpose() * &sigma_inv * &r)[(0, 0)];
    let dof = l - (n - 1); // = m − n
    let p_value = chi2_sf(chi2, dof);

    let wcov = (a.transpose() * &sigma_inv * &a).try_inverse()?;
    let mut weights = Vec::with_capacity(n);
    weights.push(1.0 - w.iter().sum::<f64>()); // base source
    weights.extend(w.iter().copied());
    let mut std_errors = vec![0.0f64; n];
    for i in 0..n - 1 {
        std_errors[i + 1] = wcov[(i, i)].max(0.0).sqrt();
    }
    // Var(w_base) = Var(Σ wᵢ) = 1ᵀ Cov(w) 1.
    let ones = DVector::from_element(n - 1, 1.0);
    std_errors[0] = (ones.transpose() * &wcov * &ones)[(0, 0)].max(0.0).sqrt();

    Some(QpAdmFit {
        weights,
        std_errors,
        chi2,
        dof,
        p_value,
        n_sites: est.n_sites,
        n_blocks: est.n_blocks,
    })
}

/// Upper tail of the χ² distribution, `P(χ²_k ≥ x)`, via the regularized upper incomplete gamma
/// `Q(k/2, x/2)`. Used for the qpAdm model-fit p-value.
fn chi2_sf(x: f64, k: usize) -> f64 {
    if k == 0 {
        return if x <= 0.0 { 1.0 } else { 0.0 };
    }
    if x <= 0.0 {
        return 1.0;
    }
    gammq(k as f64 / 2.0, x / 2.0)
}

/// `ln Γ(x)` via the Lanczos approximation (g=7), with the reflection formula for `x < 0.5`.
fn ln_gamma(x: f64) -> f64 {
    const C: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_1,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        std::f64::consts::PI.ln() - (std::f64::consts::PI * x).sin().ln() - ln_gamma(1.0 - x)
    } else {
        let x = x - 1.0;
        let t = x + 7.5;
        let a = C[0] + (1..9).map(|i| C[i] / (x + i as f64)).sum::<f64>();
        0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
    }
}

/// Lower regularized incomplete gamma `P(a,x)` by series expansion (converges fast for `x < a+1`).
fn gser(a: f64, x: f64) -> f64 {
    let gln = ln_gamma(a);
    let mut ap = a;
    let mut del = 1.0 / a;
    let mut sum = del;
    for _ in 0..1000 {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if del.abs() < sum.abs() * 1e-15 {
            break;
        }
    }
    sum * (-x + a * x.ln() - gln).exp()
}

/// Upper regularized incomplete gamma `Q(a,x)` by the Lentz continued fraction (for `x ≥ a+1`).
fn gcf(a: f64, x: f64) -> f64 {
    let gln = ln_gamma(a);
    let tiny = 1e-30;
    let mut b = x + 1.0 - a;
    let mut c = 1.0 / tiny;
    let mut d = 1.0 / b;
    let mut h = d;
    for i in 1..1000 {
        let an = -(i as f64) * (i as f64 - a);
        b += 2.0;
        d = an * d + b;
        if d.abs() < tiny {
            d = tiny;
        }
        c = b + an / c;
        if c.abs() < tiny {
            c = tiny;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 1e-15 {
            break;
        }
    }
    (-x + a * x.ln() - gln).exp() * h
}

/// `Q(a,x) = 1 − P(a,x)`, the regularized upper incomplete gamma.
fn gammq(a: f64, x: f64) -> f64 {
    if x < 0.0 || a <= 0.0 {
        return f64::NAN;
    }
    if x < a + 1.0 {
        1.0 - gser(a, x)
    } else {
        gcf(a, x)
    }
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
    super_population_summary.sort_by(|a, b| {
        b.percentage
            .partial_cmp(&a.percentage)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

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
    let completeness = if total_snps == 0 {
        0.0
    } else {
        snps_with_data as f64 / total_snps as f64
    };
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
            gt: None,
            allele_depths: None,
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
        let panel = AncestryPanel {
            build: "t".into(),
            populations: vec!["A".into(), "B".into()],
            sites,
        };
        // Half the sites are no-calls (dosage -1).
        let genotypes: Vec<SiteGenotype> = (1..=10).map(|p| sg("chr1", p, if p <= 5 { 2 } else { -1 })).collect();

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
    }

    #[test]
    fn admixture_resolves_pure_population() {
        let sites: Vec<PanelSite> = (1..=40)
            .map(|pos| PanelSite {
                contig: "chr1".to_string(),
                position: pos,
                reference_allele: 'A',
                alternate_allele: 'G',
                freqs: if pos % 2 == 0 {
                    vec![0.95, 0.05]
                } else {
                    vec![0.05, 0.95]
                },
            })
            .collect();
        let panel = AncestryPanel {
            build: "t".into(),
            populations: vec!["A".into(), "B".into()],
            sites,
        };
        // Genotype to match A: hom-alt (2) at A-rich even sites, hom-ref (0) at A-poor odd sites.
        let genos: Vec<SiteGenotype> = (1..=40)
            .map(|p| sg("chr1", p, if p % 2 == 0 { 2 } else { 0 }))
            .collect();

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
        let panel = AncestryPanel {
            build: "t".into(),
            populations: vec!["A".into(), "B".into()],
            sites,
        };
        let genos: Vec<SiteGenotype> = (1..=60).map(|p| sg("chr1", p, 1)).collect(); // all het
        let r = estimate_admixture(&genos, &panel, "t");
        let a = r
            .components
            .iter()
            .find(|c| c.population_code == "A")
            .unwrap()
            .percentage;
        assert!((40.0..=60.0).contains(&a), "A% = {a} (expected ~50)");
    }

    #[test]
    fn panel_subset_projects_and_reorders_columns() {
        let sites = vec![PanelSite {
            contig: "chr1".into(),
            position: 1,
            reference_allele: 'A',
            alternate_allele: 'G',
            freqs: vec![0.1, 0.2, 0.3],
        }];
        let p = AncestryPanel {
            build: "t".into(),
            populations: vec!["GBR".into(), "YRI".into(), "Steppe".into()],
            sites,
        };
        let s = p.subset(&["YRI", "GBR"]); // reorder + drop the absent-from-list "Steppe"
        assert_eq!(s.populations, vec!["YRI".to_string(), "GBR".to_string()]);
        assert_eq!(s.sites[0].freqs, vec![0.2, 0.1]); // columns follow the requested order
    }

    #[test]
    fn fine_admixture_restricts_to_modern_subset_and_labels_method() {
        // A fine panel with two modern pops + one ancient (Steppe). The modern subset must drop the
        // ancient column, and the result is labeled FINE_ADMIXTURE rolling up to super-pops.
        let sites: Vec<PanelSite> = (1..=40)
            .map(|pos| PanelSite {
                contig: "chr1".into(),
                position: pos,
                reference_allele: 'A',
                alternate_allele: 'G',
                freqs: vec![0.98, 0.02, 0.5], // GBR alt-rich, YRI alt-poor, Steppe middling
            })
            .collect();
        let fine = AncestryPanel {
            build: "t".into(),
            populations: vec!["GBR".into(), "YRI".into(), "Steppe".into()],
            sites,
        };
        let genos: Vec<SiteGenotype> = (1..=40).map(|p| sg("chr1", p, 2)).collect(); // all hom-alt → GBR
        let r = estimate_fine_admixture(&genos, &fine, "t");
        assert_eq!(r.method, "FINE_ADMIXTURE");
        assert_eq!(r.panel_type, "fine");
        // Ancient component excluded (not in the modern fine-code list).
        assert!(r.components.iter().all(|c| c.population_code != "Steppe"));
        let gbr = r.components.iter().find(|c| c.population_code == "GBR").unwrap();
        assert!(gbr.percentage > 90.0, "GBR% = {}", gbr.percentage);
        // Fine codes roll up to their super-pop (GBR → EUR).
        assert!(r
            .super_population_summary
            .iter()
            .any(|s| s.populations.contains(&"GBR".to_string())));
    }

    // A 2-pop panel (A alt-rich / B alt-poor) for the diploid painting tests.
    fn two_pop_panel(n: usize) -> AncestryPanel {
        let sites: Vec<PanelSite> = (0..n)
            .map(|i| PanelSite {
                contig: "chr1".to_string(),
                position: 1 + i as i64 * 1_000_000, // 1 Mb spacing
                reference_allele: 'A',
                alternate_allele: 'G',
                freqs: vec![0.95, 0.05],
            })
            .collect();
        AncestryPanel {
            build: "t".into(),
            populations: vec!["A".into(), "B".into()],
            sites,
        }
    }

    /// Ancestry-HOMOZYGOUS sample: hom-alt (→ both copies A) first half, hom-ref (→ both copies B)
    /// second half. Diploid painting emits two copies, each switching A→B at the midpoint.
    #[test]
    fn painting_diploid_homozygous_switch() {
        let n = 80;
        let panel = two_pop_panel(n);
        let genos: Vec<SiteGenotype> = (0..n)
            .map(|i| sg("chr1", 1 + i as i64 * 1_000_000, if i < n / 2 { 2 } else { 0 }))
            .collect();
        let prior = vec![("A".to_string(), 0.5), ("B".to_string(), 0.5)];
        let segs = paint_local_ancestry(&genos, &panel, &prior, &PaintParams::default());
        for copy in [0u8, 1u8] {
            let c: Vec<&AncestrySegment> = segs.iter().filter(|s| s.copy == copy).collect();
            assert_eq!(c.len(), 2, "copy {copy}: expected A→B switch, got {c:?}");
            assert_eq!(
                (c[0].population_code.as_str(), c[1].population_code.as_str()),
                ("A", "B")
            );
        }
    }

    /// Ancestry-HETEROZYGOUS sample: every site het (one copy A, one copy B). Diploid painting must
    /// put A on one copy and B on the other across the whole chromosome (the case a single-track
    /// painter cannot express).
    #[test]
    fn painting_diploid_heterozygous_copies_differ() {
        let n = 60;
        let panel = two_pop_panel(n);
        let genos: Vec<SiteGenotype> = (0..n).map(|i| sg("chr1", 1 + i as i64 * 1_000_000, 1)).collect();
        let prior = vec![("A".to_string(), 0.5), ("B".to_string(), 0.5)];
        let segs = paint_local_ancestry(&genos, &panel, &prior, &PaintParams::default());
        let copy0: Vec<&AncestrySegment> = segs.iter().filter(|s| s.copy == 0).collect();
        let copy1: Vec<&AncestrySegment> = segs.iter().filter(|s| s.copy == 1).collect();
        assert_eq!(copy0.len(), 1);
        assert_eq!(copy1.len(), 1);
        // Sorted copies: copy 0 = lower-index ancestry (A), copy 1 = higher (B).
        assert_eq!(copy0[0].population_code, "A");
        assert_eq!(copy1[0].population_code, "B");
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

    // ── deep (ancient) ancestry ─────────────────────────────────────────────────────────────────
    //
    // The three-source model is the one that previously shipped fabricated numbers, so these tests
    // pin the two properties whose absence made that possible: it must recover a mixture it was
    // never told, and it must refuse a sample its sources cannot express.

    /// A deterministic LCG — the simulations below must give the same answer on every run.
    struct Lcg(u64);
    impl Lcg {
        fn next_f64(&mut self) -> f64 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }
        /// A diploid dosage drawn under HWE at alt-frequency `f`.
        fn dosage(&mut self, f: f64) -> i32 {
            (self.next_f64() < f) as i32 + (self.next_f64() < f) as i32
        }
    }

    /// A 3-source panel over `n` sites whose frequencies differ sharply between sources (so the
    /// mixture is well-conditioned), plus an "outsider" frequency track standing in for a sample
    /// from outside the sources' span (a Yoruba against WHG/ANF/Steppe).
    ///
    /// What makes the outsider genuinely unreachable is the last site pattern, where **all three
    /// sources agree** at 0.10: a mixture can only ever predict a value inside the convex hull of
    /// its sources, so at those sites every possible `q` predicts 0.10 while the outsider carries the
    /// allele at 0.95. No mixture can absorb that, which is exactly the situation the applicability
    /// gate exists to detect. Without such sites, a near-pure single source approximates the outsider
    /// well enough to slip under the threshold.
    fn ancient_panel(n: i64) -> (AncestryPanel, Vec<f64>) {
        let mut sites = Vec::new();
        let mut outsider = Vec::new();
        for pos in 1..=n {
            // Cycle through contrasting frequency patterns so every source is identifiable.
            let (a, b, c, out) = match pos % 6 {
                0 => (0.90, 0.10, 0.50, 0.02),
                1 => (0.10, 0.90, 0.50, 0.98),
                2 => (0.50, 0.10, 0.90, 0.02),
                3 => (0.10, 0.50, 0.10, 0.95),
                // All three sources agree → every mixture predicts the same value, and the outsider
                // carries the opposite allele. Unreachable by ANY `q`.
                4 => (0.05, 0.05, 0.05, 0.98),
                _ => (0.95, 0.95, 0.95, 0.02),
            };
            sites.push(PanelSite {
                contig: "chr1".to_string(),
                position: pos,
                reference_allele: 'A',
                alternate_allele: 'G',
                freqs: vec![a, b, c],
            });
            outsider.push(out);
        }
        (
            AncestryPanel {
                build: "t".into(),
                populations: vec!["WHG".into(), "ANF".into(), "Steppe".into()],
                sites,
            },
            outsider,
        )
    }

    fn pct(r: &AncestryResult, code: &str) -> f64 {
        r.components
            .iter()
            .find(|c| c.population_code == code)
            .map_or(0.0, |c| c.percentage)
    }

    /// A stand-in modern super-population estimate that is `eur`% European — the scope input the
    /// deep model gates on.
    fn modern_eur(eur: f64) -> AncestryResult {
        let probs = [("EUR".to_string(), eur), ("AFR".to_string(), 100.0 - eur)];
        from_probabilities("ADMIXTURE", "aims", 1000, 1000, &probs, 0.9, "t")
    }

    /// The estimator recovers a mixture it was never given: simulate a 20/30/50 individual from the
    /// panel's own source frequencies and the EM must return ~20/30/50, with a dispersion near the
    /// model's noise floor of 1. This is the property the PCA-centroid classifier never had — it
    /// answered "which source *is* this?", so a genuine mixture came back as a single population.
    #[test]
    fn ancient_admixture_recovers_a_known_mixture() {
        let (panel, _) = ancient_panel(4000);
        let truth = [0.20, 0.30, 0.50];
        let mut rng = Lcg(12345);
        let genos: Vec<SiteGenotype> = panel
            .sites
            .iter()
            .map(|s| {
                let f: f64 = (0..3).map(|k| truth[k] * s.freqs[k] as f64).sum();
                sg("chr1", s.position, rng.dosage(f))
            })
            .collect();

        let r = estimate_ancient_admixture(&genos, &panel, &modern_eur(95.0), "t")
            .expect("a simulated mixture must be reportable");
        assert_eq!(r.method, ANCIENT_ADMIXTURE);
        assert_eq!(r.panel_type, "ancient");
        for (code, want) in [("WHG", 20.0), ("ANF", 30.0), ("Steppe", 50.0)] {
            let got = pct(&r, code);
            assert!((got - want).abs() < 4.0, "{code}: got {got:.1}, want ~{want}");
        }
        let d = r.fit_distance.expect("dispersion attached");
        assert!(d > 0.5 && d < 1.5, "dispersion of a well-specified sample = {d}");
    }

    /// A sample from outside the sources' span is **rejected**, not decomposed. The EM will always
    /// return *some* simplex vector — that vector is exactly what the old implementation printed as
    /// a result — so the applicability gate, not the EM, is what makes this safe.
    #[test]
    fn ancient_admixture_rejects_a_sample_outside_the_sources() {
        let (panel, outsider) = ancient_panel(4000);
        let mut rng = Lcg(99);
        let genos: Vec<SiteGenotype> = panel
            .sites
            .iter()
            .zip(&outsider)
            .map(|(s, &f)| sg("chr1", s.position, rng.dosage(f)))
            .collect();

        // The raw fit still produces a confident-looking breakdown …
        let raw = ancient_admixture_fit(&genos, &panel, "t").expect("enough sites to fit");
        let total: f64 = raw.components.iter().map(|c| c.percentage).sum();
        assert!((total - 100.0).abs() < 1e-6, "the EM always returns a full simplex");
        assert!(
            raw.fit_distance.unwrap() > ANCIENT_MAX_DISPERSION,
            "an out-of-span sample must be driven above the dispersion threshold"
        );
        // … and the shipping estimator refuses to report it, even for a sample the modern model
        // calls European (so this is the dispersion gate doing the work, not the scope gate).
        assert!(estimate_ancient_admixture(&genos, &panel, &modern_eur(95.0), "t").is_none());
    }

    /// Too few genotyped sites → no estimate at all, rather than a noisy one.
    #[test]
    fn ancient_admixture_needs_enough_sites() {
        let (panel, _) = ancient_panel(4000);
        let mut rng = Lcg(7);
        let genos: Vec<SiteGenotype> = panel
            .sites
            .iter()
            .take(ANCIENT_MIN_SITES - 1)
            .map(|s| sg("chr1", s.position, rng.dosage(s.freqs[0] as f64)))
            .collect();
        assert!(estimate_ancient_admixture(&genos, &panel, &modern_eur(95.0), "t").is_none());
    }

    /// A WHG/ANF/Steppe decomposition is a *West-Eurasian* model. A sample the modern estimate calls
    /// mostly non-European is out of scope and gets nothing — even if its genotypes happen to fit the
    /// three sources well, because "fits the arithmetic" is not the same as "means anything".
    #[test]
    fn ancient_admixture_is_scoped_to_european_samples() {
        let (panel, _) = ancient_panel(4000);
        let truth = [0.20, 0.30, 0.50];
        let mut rng = Lcg(12345);
        let genos: Vec<SiteGenotype> = panel
            .sites
            .iter()
            .map(|s| {
                let f: f64 = (0..3).map(|k| truth[k] * s.freqs[k] as f64).sum();
                sg("chr1", s.position, rng.dosage(f))
            })
            .collect();

        // Same genotypes, same perfect fit — only the scope differs.
        assert!(estimate_ancient_admixture(&genos, &panel, &modern_eur(95.0), "t").is_some());
        assert!(estimate_ancient_admixture(&genos, &panel, &modern_eur(20.0), "t").is_none());
    }

    // ── f-statistics core (Lever 2 / qpAdm) ─────────────────────────────────────────────────────
    //
    // Three properties pin the f4 primitive against known-value graphs: the exact f4-ratio algebra
    // qpAdm rests on, the antisymmetries of f4, and a *calibrated* block-jackknife SE (a symmetric
    // tree reads f4 ≈ 0 within noise, while a real internal edge reads many SE from zero).

    /// A panel over `pops` whose per-site frequency rows are `freqs[site][pop]`, laid out one site
    /// per 100 kb so the 5 Mb block jackknife sees ~50 sites/block. Plus a target genotyped
    /// (dosage 0) at every site, so every site is informative for reference-only quartets.
    fn f4_panel(pops: &[&str], freqs: &[Vec<f32>]) -> (AncestryPanel, Vec<SiteGenotype>) {
        let sites: Vec<PanelSite> = freqs
            .iter()
            .enumerate()
            .map(|(i, row)| PanelSite {
                contig: "chr1".into(),
                position: (i as i64 + 1) * 100_000,
                reference_allele: 'A',
                alternate_allele: 'G',
                freqs: row.clone(),
            })
            .collect();
        let genos = sites.iter().map(|s| sg("chr1", s.position, 0)).collect();
        let panel = AncestryPanel {
            build: "t".into(),
            populations: pops.iter().map(|s| s.to_string()).collect(),
            sites,
        };
        (panel, genos)
    }

    /// The f4-ratio identity qpAdm generalizes: if `X = α·P + (1−α)·Q` (a frequency mixture), then
    /// `f4(X,P;O1,O2) = (1−α)·f4(Q,P;O1,O2)` **exactly** — per site the two are proportional. So the
    /// ratio recovers `1−α` regardless of the outgroups. This is the core arithmetic the whole method
    /// stands on; a wrong sign or a transposed index would blow the recovered α far past f32 noise.
    #[test]
    fn f4_ratio_recovers_the_mixture_weight() {
        let alpha = 0.3_f64;
        let mut rng = Lcg(42);
        let freqs: Vec<Vec<f32>> = (0..4000)
            .map(|_| {
                let (r1, r2) = (rng.next_f64(), rng.next_f64());
                let p = 0.1 + 0.8 * r1;
                let q = 0.1 + 0.8 * r2;
                // Outgroups tied to the sources so f4(Q,P;O1,O2) is large and clean (no small-denom
                // amplification): O1−O2 = 0.7(r1−r2) tracks −(Q−P), giving a firmly non-zero f4.
                let o1 = 0.15 + 0.7 * r1;
                let o2 = 0.15 + 0.7 * r2;
                let x = alpha * p + (1.0 - alpha) * q;
                vec![x as f32, p as f32, q as f32, o1 as f32, o2 as f32]
            })
            .collect();
        let (panel, genos) = f4_panel(&["X", "P", "Q", "O1", "O2"], &freqs);
        let (x, p, q, o1, o2) = (Pop::Ref(0), Pop::Ref(1), Pop::Ref(2), Pop::Ref(3), Pop::Ref(4));
        let est = f4_vector(
            &genos,
            &panel,
            &[Quartet::new(x, p, o1, o2), Quartet::new(q, p, o1, o2)],
            F4_BLOCK_BP,
        )
        .expect("f4 vector");
        assert!(est.values[1].abs() > 0.02, "denominator f4 must be firmly non-degenerate");
        let recovered = 1.0 - est.values[0] / est.values[1];
        assert!((recovered - alpha).abs() < 1e-4, "f4-ratio recovered α={recovered:.6}, want {alpha}");
    }

    /// f4's exact symmetries (pure f64 arithmetic over one fixed site set): swapping either pair
    /// negates it, and swapping the two pairs leaves it unchanged.
    #[test]
    fn f4_obeys_its_antisymmetries() {
        let mut rng = Lcg(7);
        let freqs: Vec<Vec<f32>> = (0..2000)
            .map(|_| (0..4).map(|_| (0.05 + 0.9 * rng.next_f64()) as f32).collect())
            .collect();
        let (panel, genos) = f4_panel(&["A", "B", "C", "D"], &freqs);
        let (a, b, c, d) = (Pop::Ref(0), Pop::Ref(1), Pop::Ref(2), Pop::Ref(3));
        let est = f4_vector(
            &genos,
            &panel,
            &[
                Quartet::new(a, b, c, d),
                Quartet::new(b, a, c, d),
                Quartet::new(a, b, d, c),
                Quartet::new(c, d, a, b),
            ],
            F4_BLOCK_BP,
        )
        .expect("f4 vector");
        let base = est.values[0];
        assert!(base.abs() > 1e-9, "pick a non-degenerate base statistic");
        assert!((est.values[1] + base).abs() < 1e-12, "f4(b,a;c,d) = −f4(a,b;c,d)");
        assert!((est.values[2] + base).abs() < 1e-12, "f4(a,b;d,c) = −f4(a,b;c,d)");
        assert!((est.values[3] - base).abs() < 1e-12, "f4(c,d;a,b) = f4(a,b;c,d)");
        assert!(est.se(0) >= 0.0);
    }

    /// A symmetric tree `((A,B),(C,D))` has `f4(A,B;C,D) = 0` in expectation (the A–B and C–D drift
    /// paths don't overlap), while `f4(A,C;B,D)` sits on the shared internal edge and is non-zero.
    /// Simulate exactly that and require the jackknife SE to *tell them apart*: the null within a few
    /// SE of zero, the real edge many SE away. This is the test that the covariance is calibrated —
    /// the property §5.4 needs and simulation-of-frequencies alone can't fake.
    #[test]
    fn f4_jackknife_se_separates_a_null_from_a_real_edge() {
        let mut rng = Lcg(2024);
        let freqs: Vec<Vec<f32>> = (0..5000)
            .map(|_| {
                let cab = 0.5 + 0.3 * (rng.next_f64() - 0.5); // drift shared by A,B
                let ccd = 0.5 + 0.3 * (rng.next_f64() - 0.5); // drift shared by C,D
                let tip = |rng: &mut Lcg, c: f64| (c + 0.15 * (rng.next_f64() - 0.5)) as f32;
                vec![tip(&mut rng, cab), tip(&mut rng, cab), tip(&mut rng, ccd), tip(&mut rng, ccd)]
            })
            .collect();
        let (panel, genos) = f4_panel(&["A", "B", "C", "D"], &freqs);
        let (a, b, c, d) = (Pop::Ref(0), Pop::Ref(1), Pop::Ref(2), Pop::Ref(3));
        let est = f4_vector(
            &genos,
            &panel,
            &[Quartet::new(a, b, c, d), Quartet::new(a, c, b, d)],
            F4_BLOCK_BP,
        )
        .expect("f4 vector");
        let z_null = est.values[0] / est.se(0);
        let z_edge = est.values[1] / est.se(1);
        assert!(z_null.abs() < 4.0, "symmetric tree: f4(A,B;C,D) must sit near 0, z={z_null:.2}");
        assert!(z_edge.abs() > 8.0, "real internal edge: f4(A,C;B,D) must be many SE from 0, z={z_edge:.2}");
        assert!(
            est.values[0].abs() * 5.0 < est.values[1].abs(),
            "the null statistic must be far smaller than the real edge"
        );
    }

    /// The χ² upper tail against textbook critical points — the model-fit p-value depends on it.
    #[test]
    fn chi2_sf_matches_known_critical_values() {
        assert!((chi2_sf(3.841, 1) - 0.05).abs() < 2e-3);
        assert!((chi2_sf(5.991, 2) - 0.05).abs() < 2e-3);
        assert!((chi2_sf(7.815, 3) - 0.05).abs() < 2e-3);
        assert!((chi2_sf(0.455, 1) - 0.50).abs() < 5e-3);
        assert!((chi2_sf(11.345, 3) - 0.01).abs() < 2e-3);
        assert_eq!(chi2_sf(0.0, 3), 1.0);
        assert!(chi2_sf(100.0, 1) < 1e-10);
    }

    /// A Gaussian increment (Box-Muller) for the drift simulation.
    fn gauss(rng: &mut Lcg, sd: f64) -> f64 {
        let u1 = rng.next_f64().max(1e-12);
        let u2 = rng.next_f64();
        sd * (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }

    /// Simulate frequencies on a small admixture graph: three sources S1/S2/S3, each carrying a
    /// distinct deep component (iA/iB/iC); six outgroups differentially related to those components
    /// (R0 pure near-root; R1,R4→A; R2,R5→B; R3→C); and a target that is an exact per-site frequency
    /// mixture of the three sources, drawn as one diploid genome. Under Brownian drift the f4
    /// tree-identities hold, so this is a graph qpAdm is entitled to decompose.
    fn qpadm_graph(n_sites: usize, weights: [f64; 3], seed: u64) -> (AncestryPanel, Vec<SiteGenotype>) {
        let mut rng = Lcg(seed);
        let clamp = |x: f64| x.clamp(0.02, 0.98);
        let pops = ["S1", "S2", "S3", "R0", "R1", "R2", "R3", "R4", "R5"];
        let mut sites = Vec::with_capacity(n_sites);
        let mut genos = Vec::with_capacity(n_sites);
        for s in 0..n_sites {
            let pos = (s as i64 + 1) * 100_000;
            let p0 = 0.3 + 0.4 * rng.next_f64();
            let (ia, ib, ic) = (gauss(&mut rng, 0.06), gauss(&mut rng, 0.06), gauss(&mut rng, 0.06));
            let pv = |rng: &mut Lcg| gauss(rng, 0.04);
            let f_s1 = clamp(p0 + ia + pv(&mut rng));
            let f_s2 = clamp(p0 + ib + pv(&mut rng));
            let f_s3 = clamp(p0 + ic + pv(&mut rng));
            let f_r0 = clamp(p0 + pv(&mut rng));
            let f_r1 = clamp(p0 + ia + pv(&mut rng));
            let f_r2 = clamp(p0 + ib + pv(&mut rng));
            let f_r3 = clamp(p0 + ic + pv(&mut rng));
            let f_r4 = clamp(p0 + ia + pv(&mut rng));
            let f_r5 = clamp(p0 + ib + pv(&mut rng));
            let row = [f_s1, f_s2, f_s3, f_r0, f_r1, f_r2, f_r3, f_r4, f_r5];
            // Target = exact frequency mixture of the three sources, drawn as a diploid genome.
            let f_t = clamp(weights[0] * f_s1 + weights[1] * f_s2 + weights[2] * f_s3);
            genos.push(sg("chr1", pos, rng.dosage(f_t)));
            sites.push(PanelSite {
                contig: "chr1".into(),
                position: pos,
                reference_allele: 'A',
                alternate_allele: 'G',
                freqs: row.iter().map(|&x| x as f32).collect(),
            });
        }
        let panel = AncestryPanel {
            build: "t".into(),
            populations: pops.iter().map(|s| s.to_string()).collect(),
            sites,
        };
        (panel, genos)
    }

    /// qpAdm recovers the true mixture weights from allele-sharing against the outgroups, accepts the
    /// well-specified model, and **rejects** a model missing a needed source — the property §3's
    /// frequency-EM never had (it always returned a confident simplex).
    #[test]
    fn qpadm_recovers_a_known_mixture_and_rejects_a_deficient_model() {
        let truth = [0.5, 0.3, 0.2];
        let (panel, genos) = qpadm_graph(20_000, truth, 20_240_717);
        let outgroups = [3usize, 4, 5, 6, 7, 8]; // R0 base, then R1..R5

        let fit = qpadm_fit(&genos, &panel, &[0, 1, 2], &outgroups, F4_BLOCK_BP).expect("3-source fit");
        assert_eq!(fit.dof, 3, "dof = #outgroups − #sources = 6 − 3");
        for (i, &want) in truth.iter().enumerate() {
            assert!((fit.weights[i] - want).abs() < 0.08, "w{i} = {:.3}, want {want}", fit.weights[i]);
        }
        assert!(fit.weights_feasible(0.02), "weights must be valid proportions: {:?}", fit.weights);
        assert!(fit.p_value > 0.01, "well-specified model must not be rejected, p = {:.4}", fit.p_value);

        // Drop a needed source (S3): the 2-source model can't express the target's cladeC affinity,
        // so its f4 residual with the cladeC outgroup is large → rejected.
        let deficient = qpadm_fit(&genos, &panel, &[0, 1], &outgroups, F4_BLOCK_BP).expect("2-source fit");
        assert_eq!(deficient.dof, 4);
        assert!(
            deficient.p_value < 0.01,
            "deficient model must be rejected, p = {:.4}",
            deficient.p_value
        );
    }

    /// The app-facing `estimate_qpadm_ancestry`: reports the source weights for a well-specified
    /// European fit, and gates on both scope (non-European → None) and model-fit (deficient → None).
    #[test]
    fn estimate_qpadm_ancestry_reports_european_and_gates_the_rest() {
        let truth = [0.5, 0.3, 0.2];
        let (panel, genos) = qpadm_graph(20_000, truth, 20_240_717);
        let outgroups = [3usize, 4, 5, 6, 7, 8];

        let r = estimate_qpadm_ancestry(&genos, &panel, &[0, 1, 2], &outgroups, &modern_eur(95.0), "t")
            .expect("a well-specified European fit is reported");
        assert_eq!(r.method, ANCIENT_ADMIXTURE);
        assert_eq!(r.panel_type, "ancient");
        // Recovered within the underlying qpAdm test's tolerance (~8 pts), and correctly ordered.
        for (code, want) in [("S1", 50.0), ("S2", 30.0), ("S3", 20.0)] {
            assert!((pct(&r, code) - want).abs() < 9.0, "{code}: {:.1} vs {want}", pct(&r, code));
        }
        assert!(pct(&r, "S1") > pct(&r, "S2") && pct(&r, "S2") > pct(&r, "S3"), "order preserved");
        let p = r.fit_distance.expect("p-value on fit_distance");
        assert!((0.0..=1.0).contains(&p), "p={p}");

        // Scope gate: a mostly-non-European sample gets nothing, even though the arithmetic fits.
        assert!(estimate_qpadm_ancestry(&genos, &panel, &[0, 1, 2], &outgroups, &modern_eur(20.0), "t").is_none());
        // Model-fit gate: a source-deficient model is rejected by the p-value.
        assert!(estimate_qpadm_ancestry(&genos, &panel, &[0, 1], &outgroups, &modern_eur(95.0), "t").is_none());
    }
}
