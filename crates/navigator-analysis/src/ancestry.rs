//! Ancestry estimation. This is the path from genotypes to population proportions, on the
//! Navigator side.
//!
//! Phase 1 is the allele-frequency likelihood. It uses no PCA and no GATK. The bundled
//! [`AncestryPanel`] carries alt-allele frequencies for each (super-)population at a set of
//! ancestry-informative sites. The code genotypes the sample at those sites with the GL caller
//! ([`crate::caller::genotype_sites`]). It then scores each population by the binomial
//! likelihood of the observed diploid genotypes under the allele frequencies of that
//! population. `navigator-panelbuild` builds the panel offline from the 1000G-on-CHM13 VCFs.
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

/// One ancestry-informative site with the alt-allele frequency for each population. `freqs[i]`
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
    /// Population codes. They give the axis order of every `PanelSite::freqs`.
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

    /// A panel that keeps only `codes`, those that are present, in `codes` order. It projects
    /// the frequencies of each site down to the columns that stay. Use it to run a
    /// well-conditioned admixture EM over a curated subset of a large fine-frequency panel.
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

/// PCA loadings. They project a sample onto the principal-component space of the reference
/// populations (Phase 2). `navigator-panelbuild` builds them offline from the 1000G genotype
/// matrix. They hold a loading and a mean for each SNP, and the mean centres the data. They
/// also hold the centroid and the diagonal variance of each population in PC space, for the
/// Mahalanobis/Gaussian assignment and for the scatter plot.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PcaLoadings {
    pub build: String,
    /// A (contig, 1-based pos) for each row, aligned with `means` and the rows of `loadings`.
    pub sites: Vec<(String, i64)>,
    /// The mean dosage at each site, from the reference panel. It centres the sample before
    /// the projection.
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

/// One reference site for [`HaplotypeReference`]: its coordinate and the biallelic ref/alt the
/// packed haplotype bit refers to (`1` = alt allele, `0` = ref).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct HapSite {
    pub contig: String,
    /// 1-based.
    pub position: i64,
    pub reference_allele: char,
    pub alternate_allele: char,
}

/// A reference panel of **phased haplotypes** at the painting loci. It is the substrate for
/// statistical phasing, which uses the Li & Stephens copying model, and, later, for
/// local-ancestry inference with the copying model.
///
/// It is different from [`AncestryPanel`], which stores allele *frequencies* for each
/// population. This panel carries the allele of each individual reference *haplotype* at every
/// site, and the population label of that haplotype. The code can then phase a sample as a
/// mosaic of these haplotypes. `navigator-panelbuild` builds it offline from the **phased**
/// 1000G-on-CHM13 VCFs.
///
/// Only modern, phased references go into it: the 1000G super populations and fine populations.
/// Pseudo-haploid and unphased sources, which are the ancient data and the AADR continental
/// groups, stay out. Those sources go into the frequency panel that the two-tier
/// fine-resolution step uses.
///
/// The alleles are bit-packed row-major. The allele of haplotype `h` at site `s` is bit
/// `h * n_sites + s` of [`alleles`](Self::alleles). This keeps a reference of about 5000
/// haplotypes by about 20k sites near 12 MB.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HaplotypeReference {
    /// Canonical reference build the site coordinates are in (e.g. "chm13v2.0").
    pub build: String,
    /// Reference sites, sorted by `(contig, position)`; the column order of every haplotype row.
    pub sites: Vec<HapSite>,
    /// Distinct population codes; `hap_pop[h]` indexes into this axis.
    pub populations: Vec<String>,
    /// One entry for each haplotype: the index of its label in [`populations`](Self::populations).
    pub hap_pop: Vec<u16>,
    /// `n_haplotypes × n_sites` alleles, bit-packed row-major (see the type doc). `1` = alt.
    pub alleles: Vec<u64>,
    pub n_sites: usize,
    pub n_haplotypes: usize,
}

impl HaplotypeReference {
    /// Deserialize from the bundled/built binary (bincode).
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, AnalysisError> {
        bincode::deserialize(bytes).map_err(|e| AnalysisError::Message(format!("hap reference decode: {e}")))
    }

    /// Serialize to the binary form the builder writes and the app bundles.
    pub fn to_bytes(&self) -> Result<Vec<u8>, AnalysisError> {
        bincode::serialize(self).map_err(|e| AnalysisError::Message(format!("hap reference encode: {e}")))
    }

    /// Pack the allele row of each haplotype (`rows[h][s]` = 0/1) into the bit-packed form.
    /// `hap_pop[h]` is the population index of haplotype `h`. The offline builder and the tests
    /// use this.
    pub fn from_rows(
        build: String,
        sites: Vec<HapSite>,
        populations: Vec<String>,
        hap_pop: Vec<u16>,
        rows: &[Vec<u8>],
    ) -> Self {
        let n_sites = sites.len();
        let n_haplotypes = rows.len();
        let total_bits = n_sites * n_haplotypes;
        let mut alleles = vec![0u64; total_bits.div_ceil(64)];
        for (h, row) in rows.iter().enumerate() {
            let base = h * n_sites;
            for (s, &a) in row.iter().enumerate() {
                if a != 0 {
                    let bit = base + s;
                    alleles[bit / 64] |= 1u64 << (bit % 64);
                }
            }
        }
        HaplotypeReference {
            build,
            sites,
            populations,
            hap_pop,
            alleles,
            n_sites,
            n_haplotypes,
        }
    }

    /// Allele (0 = ref, 1 = alt) of haplotype `hap` at site index `site`.
    #[inline]
    pub fn allele(&self, hap: usize, site: usize) -> u8 {
        let bit = hap * self.n_sites + site;
        ((self.alleles[bit / 64] >> (bit % 64)) & 1) as u8
    }

    /// Population code of haplotype `hap`.
    pub fn population_of(&self, hap: usize) -> &str {
        &self.populations[self.hap_pop[hap] as usize]
    }

    pub fn is_empty(&self) -> bool {
        self.n_sites == 0 || self.n_haplotypes == 0
    }

    /// The same reference with the haplotypes in `drop` removed. The sites, the population axis
    /// and every other haplotype do not change. The leave-one-out check of the copying painter
    /// needs this. A test individual taken *from* the reference would else get its paint from a
    /// copy of itself, which measures nothing. The code ignores an index that it does not know.
    pub fn without_haplotypes(&self, drop: &[usize]) -> Self {
        let dropped: std::collections::HashSet<usize> = drop.iter().copied().collect();
        let keep: Vec<usize> = (0..self.n_haplotypes).filter(|h| !dropped.contains(h)).collect();
        let mut alleles = vec![0u64; (keep.len() * self.n_sites).div_ceil(64)];
        for (new_h, &old_h) in keep.iter().enumerate() {
            let base = new_h * self.n_sites;
            for s in 0..self.n_sites {
                if self.allele(old_h, s) != 0 {
                    let bit = base + s;
                    alleles[bit / 64] |= 1u64 << (bit % 64);
                }
            }
        }
        HaplotypeReference {
            build: self.build.clone(),
            sites: self.sites.clone(),
            populations: self.populations.clone(),
            hap_pop: keep.iter().map(|&h| self.hap_pop[h]).collect(),
            alleles,
            n_sites: self.n_sites,
            n_haplotypes: keep.len(),
        }
    }

    /// The same reference thinned to every `step`-th site. A `step` of 1 or less returns a
    /// clone. Marker density is the constraint that limits the copying model. A shared tract
    /// identifies whose haplotype it is only if enough markers fall inside it. To measure
    /// accuracy against density you need the same panel at more than one density. This
    /// function makes those panels, and it does not build the asset again.
    pub fn thin_sites(&self, step: usize) -> Self {
        if step <= 1 {
            return self.clone();
        }
        let keep: Vec<usize> = (0..self.n_sites).step_by(step).collect();
        let n_sites = keep.len();
        let mut alleles = vec![0u64; (self.n_haplotypes * n_sites).div_ceil(64)];
        for h in 0..self.n_haplotypes {
            let base = h * n_sites;
            for (new_s, &old_s) in keep.iter().enumerate() {
                if self.allele(h, old_s) != 0 {
                    let bit = base + new_s;
                    alleles[bit / 64] |= 1u64 << (bit % 64);
                }
            }
        }
        HaplotypeReference {
            build: self.build.clone(),
            sites: keep.iter().map(|&s| self.sites[s].clone()).collect(),
            populations: self.populations.clone(),
            hap_pop: self.hap_pop.clone(),
            alleles,
            n_sites,
            n_haplotypes: self.n_haplotypes,
        }
    }
}

/// Project the genotypes of a sample onto the reference PCA space. Centre each site by its
/// panel mean, then add `centered · loading` into each component. A missing genotype adds 0,
/// which imputes the mean. The code then scales the projection by `total_sites / sites_used`.
/// Without that scale, a sample with missing genotypes moves toward the origin and away from
/// its true cluster. Returns the coordinate of the sample in each principal component.
pub fn project_pca(genotypes: &[SiteGenotype], pca: &PcaLoadings) -> Vec<f64> {
    let dosage: HashMap<(&str, i64), i32> = genotypes
        .iter()
        .filter(|g| g.dosage >= 0)
        .map(|g| ((g.contig.as_str(), g.position), g.dosage))
        .collect();

    let centered = pca.sites.iter().enumerate().filter_map(|(i, (contig, pos))| {
        let &d = dosage.get(&(contig.as_str(), *pos))?;
        Some((i, d as f64 - pca.means[i] as f64))
    });
    project_centered(pca.sites.len(), pca.n_components, centered, |i, c| {
        pca.loading(i, c) as f64
    })
}

/// The PCA projection kernel. It adds `centered · loading` into each component, over the sites
/// that the sample has. It then scales up by `n_sites / used`, so that a sample with missing
/// genotypes does not move toward the origin. See [`project_pca`].
///
/// `centered` gives `(site index, dosage − site mean)` for each site that is present. `loading`
/// reads the `(site, component)` basis entry. The caller gives both, because the runtime
/// projector and the offline basis builder hold their basis in different layouts: `PcaLoadings`
/// against a `DMatrix`. The scale policy must agree between the two, so it lives here.
pub fn project_centered(
    n_sites: usize,
    n_components: usize,
    centered: impl Iterator<Item = (usize, f64)>,
    loading: impl Fn(usize, usize) -> f64,
) -> Vec<f64> {
    let mut coords = vec![0.0f64; n_components];
    let mut used = 0usize;
    for (i, value) in centered {
        used += 1;
        for (c, coord) in coords.iter_mut().enumerate() {
            *coord += value * loading(i, c);
        }
    }
    // The reference coordinates come from all of the sites. Scale up for the fraction that is
    // missing.
    if used > 0 {
        let scale = n_sites as f64 / used as f64;
        for coord in &mut coords {
            *coord *= scale;
        }
    }
    coords
}

/// Parameters for [`paint_local_ancestry`].
#[derive(Debug, Clone)]
pub struct PaintParams {
    /// The ancestry-switch rate for each bp, which controls the segment length. The switch
    /// probability over a distance of `d` bp is `1 - exp(-d·rate)`. A smaller rate gives longer
    /// segments. The default is about one switch in 20 Mb.
    pub rate: f64,
    /// The code merges a run of fewer markers than this into the segment that is next to it.
    pub min_segment_sites: usize,
    /// The gate on the global composition. The HMM drops a super-population from its state set
    /// if the genome-wide `prior` weight of that population is below this fraction. It always
    /// keeps the dominant ancestry. A donor who is 99% European can then not get a *local*
    /// East-Asian or South-Asian paint from a few noise loci. `0.0` turns the gate off.
    ///
    /// Local ancestry from allele frequencies over coarse super-populations can show a continent
    /// that the genome does not contain at all. To hold the states to the global estimate stops
    /// that.
    pub min_ancestry: f64,
}

impl Default for PaintParams {
    fn default() -> Self {
        Self {
            rate: 1.0 / 20_000_000.0,
            min_segment_sites: 5,
            min_ancestry: 0.02,
        }
    }
}

/// The log-likelihood of a diploid genotype. The two genome copies draw their alt allele from
/// the frequencies `fa` and `fb`, one copy for each ancestry: `P(0)=(1-fa)(1-fb)`,
/// `P(1)=fa(1-fb)+(1-fa)fb`, `P(2)=fa·fb`. A missing dosage gives a uniform value. This is the
/// correct diploid emission, over two copies, that the pair-state HMM needs.
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

/// Paint each chromosome with local ancestry. The model is an HMM over the panel sites. Its
/// hidden states are the super-populations. Its emissions are the diploid genotype likelihood
/// under the allele frequency of each population. Its transitions penalise an ancestry switch
/// by physical distance. Viterbi gives the segment path, and forward-backward gives the
/// posterior at each site, which is the segment confidence.
///
/// `prior` is the genome-wide composition `(population_code, weight)`, which this function rolls
/// up to super-populations. It is the stationary and switch distribution of the HMM, and it
/// holds the painting to the global estimate.
///
/// **The HMM has diploid pair states**: one hidden state is a pair of ancestries, for both
/// genome copies. A region where the two copies are different, for example EUR and SAS, stays
/// visible and does not collapse. The output is two sorted, unphased copies for each chromosome,
/// with the segments tagged `copy` 0 or 1. The copies are not maternal and paternal, because
/// there is no phasing.
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
    let mut all_states: Vec<String> = Vec::new();
    for s in &pop_state {
        if !all_states.contains(s) {
            all_states.push(s.clone());
        }
    }
    if all_states.is_empty() {
        return Vec::new();
    }

    // The prior π over the full state set. Roll the global composition up to the
    // super-populations and normalize it. Fall back to a uniform π when the caller gives no
    // prior.
    let mut full_pi = vec![0.0f64; all_states.len()];
    for (code, w) in prior {
        let sp = population_super(code).unwrap_or(code);
        if let Some(j) = all_states.iter().position(|x| x == sp) {
            full_pi[j] += w.max(0.0);
        }
    }
    let tot: f64 = full_pi.iter().sum();
    if tot > 0.0 {
        full_pi.iter_mut().for_each(|p| *p /= tot);
    } else {
        full_pi.iter_mut().for_each(|p| *p = 1.0 / all_states.len() as f64);
    }

    // The gate on the global composition. Keep only an ancestry that the genome shows
    // (>= min_ancestry), and always keep the dominant one. The local painting can then not show
    // a continent that the genome does not contain at all. With a uniform π, which is the
    // no-prior case, every state clears the default threshold, and the gate does nothing.
    let argmax = full_pi
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.total_cmp(b.1))
        .map(|(i, _)| i)
        .unwrap_or(0);
    let keep: Vec<usize> = (0..all_states.len())
        .filter(|&i| full_pi[i] >= params.min_ancestry || i == argmax)
        .collect();
    let states: Vec<String> = keep.iter().map(|&i| all_states[i].clone()).collect();
    let k = states.len();
    let state_idx = |s: &str| states.iter().position(|x| x == s);

    // π restricted to the kept states, renormalized.
    let mut pi: Vec<f64> = keep.iter().map(|&i| full_pi[i]).collect();
    let ktot: f64 = pi.iter().sum();
    if ktot > 0.0 {
        pi.iter_mut().for_each(|p| *p /= ktot);
    } else {
        pi.iter_mut().for_each(|p| *p = 1.0 / k as f64);
    }

    let dosage: HashMap<(&str, i64), i32> = genotypes
        .iter()
        .filter(|g| g.dosage >= 0)
        .map(|g| ((g.contig.as_str(), g.position), g.dosage))
        .collect();

    // The sites with a genotype, for each contig: (pos, super-pop AF for each state, dosage).
    // Sorted by pos.
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
        // The diploid MAP path. It gives one ancestry PAIR at each locus. The code puts the
        // pair into a canonical (min,max) order. That makes two sorted, coherent copies. Copy 0
        // is the ancestry with the lower index, and copy 1 is the higher. There is no phasing.
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

/// The HMM state scaffold that the two painters share. One painter is diploid and unphased, and
/// the other is haploid and works on a phased side. The scaffold holds the super-population
/// states that stay, the super-population of each panel population, and the prior π over the
/// states that stay.
///
/// The selection and the gate on the global composition are the same as in
/// [`paint_local_ancestry`]. Both painters hold to the same global estimate, and neither
/// can show a continent that the genome does not contain at all.
struct PaintStates {
    states: Vec<String>,
    pop_state: Vec<String>,
    pi: Vec<f64>,
}

fn build_paint_states(panel: &AncestryPanel, prior: &[(String, f64)], params: &PaintParams) -> Option<PaintStates> {
    let pop_state: Vec<String> = panel
        .populations
        .iter()
        .map(|c| population_super(c).unwrap_or(c).to_string())
        .collect();
    let mut all_states: Vec<String> = Vec::new();
    for s in &pop_state {
        if !all_states.contains(s) {
            all_states.push(s.clone());
        }
    }
    if all_states.is_empty() {
        return None;
    }

    let mut full_pi = vec![0.0f64; all_states.len()];
    for (code, w) in prior {
        let sp = population_super(code).unwrap_or(code);
        if let Some(j) = all_states.iter().position(|x| x == sp) {
            full_pi[j] += w.max(0.0);
        }
    }
    let tot: f64 = full_pi.iter().sum();
    if tot > 0.0 {
        full_pi.iter_mut().for_each(|p| *p /= tot);
    } else {
        full_pi.iter_mut().for_each(|p| *p = 1.0 / all_states.len() as f64);
    }

    let argmax = full_pi
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.total_cmp(b.1))
        .map(|(i, _)| i)
        .unwrap_or(0);
    let keep: Vec<usize> = (0..all_states.len())
        .filter(|&i| full_pi[i] >= params.min_ancestry || i == argmax)
        .collect();
    let states: Vec<String> = keep.iter().map(|&i| all_states[i].clone()).collect();

    let mut pi: Vec<f64> = keep.iter().map(|&i| full_pi[i]).collect();
    let ktot: f64 = pi.iter().sum();
    if ktot > 0.0 {
        pi.iter_mut().for_each(|p| *p /= ktot);
    } else {
        pi.iter_mut().for_each(|p| *p = 1.0 / states.len() as f64);
    }
    Some(PaintStates { states, pop_state, pi })
}

/// The alt-allele frequency of each state at one site. It is the mean fine-population frequency
/// in each super-population state that stays. A state with no population in it gets `0.5`. The
/// order is the `states` order.
fn per_state_af(freqs: &[f32], pop_state: &[String], states: &[String]) -> Vec<f64> {
    let k = states.len();
    let mut sum = vec![0.0f64; k];
    let mut cnt = vec![0usize; k];
    for (p, &f) in freqs.iter().enumerate() {
        if let Some(j) = states.iter().position(|x| *x == pop_state[p]) {
            sum[j] += f as f64;
            cnt[j] += 1;
        }
    }
    (0..k)
        .map(|j| if cnt[j] > 0 { sum[j] / cnt[j] as f64 } else { 0.5 })
        .collect()
}

/// Haploid emission: ln P(observe allele `a` (0/1) | copied super-pop alt frequency `f`).
fn emit_haploid_ln(a: u8, f: f64) -> f64 {
    let f = f.clamp(1e-4, 1.0 - 1e-4);
    let p = if a == 1 { f } else { 1.0 - f };
    p.max(1e-300).ln()
}

/// Haploid Viterbi. It gives the MAP super-population at each site for **one** phased
/// haplotype. The hidden state is the super-population that the model copies. The emission is
/// [`emit_haploid_ln`] on the allele of that side. The transitions penalise an ancestry switch
/// by physical distance, with the same `switch_prob` and `ln_trans` as the diploid painter.
/// `sites` are `(pos, AF for each state, allele 0/1)`, in position order.
fn haploid_viterbi(sites: &[(i64, Vec<f64>, u8)], pi: &[f64], rate: f64, k: usize) -> Vec<usize> {
    let n = sites.len();
    let lnpi: Vec<f64> = (0..k).map(|s| pi[s].max(1e-300).ln()).collect();
    let mut v = vec![vec![f64::NEG_INFINITY; k]; n];
    let mut bp = vec![vec![0usize; k]; n];
    for s in 0..k {
        v[0][s] = lnpi[s] + emit_haploid_ln(sites[0].2, sites[0].1[s]);
    }
    for i in 1..n {
        let sw = switch_prob(sites[i].0 - sites[i - 1].0, rate);
        for b in 0..k {
            let (mut best, mut arg) = (f64::NEG_INFINITY, 0usize);
            for (a, &va) in v[i - 1].iter().enumerate() {
                let val = va + ln_trans(a, b, sw, pi);
                if val > best {
                    best = val;
                    arg = a;
                }
            }
            v[i][b] = best + emit_haploid_ln(sites[i].2, sites[i].1[b]);
            bp[i][b] = arg;
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

/// Paint local ancestry from **phased** genotypes. A haploid ancestry HMM runs on each of the
/// two phased sides, and the two runs are independent. The two output tracks are then parental
/// sides that agree with themselves. The segment `copy` is the phased side, 0 or 1, and it keeps
/// that sense across the whole genome. The unphased [`paint_local_ancestry`] can not make that
/// parent split.
///
/// `prior` is the genome-wide composition, which holds the state set in place. `panel` gives the
/// allele frequency of each super-population.
pub fn paint_local_ancestry_phased(
    phased: &crate::phasing::PhasedGenotypes,
    panel: &AncestryPanel,
    prior: &[(String, f64)],
    params: &PaintParams,
) -> Vec<AncestrySegment> {
    let Some(PaintStates { states, pop_state, pi }) = build_paint_states(panel, prior, params) else {
        return Vec::new();
    };
    let k = states.len();

    // The AF of each state at each site, keyed by (contig, pos). The code computes this once,
    // and both sides use it.
    let site_af: HashMap<(&str, i64), Vec<f64>> = panel
        .sites
        .iter()
        .filter(|s| s.freqs.len() == panel.populations.len())
        .map(|s| {
            (
                (s.contig.as_str(), s.position),
                per_state_af(&s.freqs, &pop_state, &states),
            )
        })
        .collect();

    let contigs: std::collections::BTreeSet<&str> = phased.sites.iter().map(|s| s.contig.as_str()).collect();
    let mut segments = Vec::new();
    for contig in contigs {
        for side in [0u8, 1u8] {
            // Build this side's position-sorted (pos, AF, allele) sites for this contig.
            let mut sites: Vec<(i64, Vec<f64>, u8)> = phased
                .sites
                .iter()
                .filter(|s| s.contig == contig)
                .filter_map(|s| {
                    site_af
                        .get(&(contig, s.position))
                        .map(|af| (s.position, af.clone(), if side == 0 { s.side0 } else { s.side1 }))
                })
                .collect();
            sites.sort_by_key(|s| s.0);
            if sites.is_empty() {
                continue;
            }
            let path = haploid_viterbi(&sites, &pi, params.rate, k);
            // collapse_copy needs (pos, _, a dosage-like value). It does not read the AF and
            // allele payload.
            let collapse_sites: Vec<(i64, Vec<f64>, i32)> =
                sites.iter().map(|s| (s.0, Vec::new(), s.2 as i32)).collect();
            segments.extend(collapse_copy(
                contig,
                &collapse_sites,
                &path,
                &states,
                params.min_segment_sites,
                side,
            ));
        }
    }
    segments
}

/// The controls of [`resolve_fine_populations`], which is the two-tier super-to-fine step.
#[derive(Debug, Clone)]
pub struct FineResolveParams {
    /// The count of informative sites that a segment needs before the code tries a fine call.
    pub min_sites: usize,
    /// The mean ln-likelihood advantage at each site that the best fine population needs over
    /// the second one before the code accepts the call. Below this, the segment keeps
    /// `fine_population_code = None`.
    pub min_margin_per_site: f64,
}

impl Default for FineResolveParams {
    fn default() -> Self {
        Self {
            min_sites: 10,
            min_margin_per_site: 0.02,
        }
    }
}

/// Two-tier fine resolution. For each super-population segment that the painter made, this
/// takes the most likely **fine** population *inside that super-population* from the
/// fine-frequency panel. It scores the phased-side alleles of the segment by the haploid
/// likelihood under each candidate fine population.
///
/// It sets [`AncestrySegment::fine_population_code`] in place. It leaves the code `None` when
/// the segment is too short, or when the best fine call is not clearly ahead of the second one.
/// This is the same shape as the super-to-fine admixture hierarchy.
///
/// `fine_panel.populations` are fine-population codes. The `freqs` of each site are the AF of
/// each fine population.
pub fn resolve_fine_populations(
    segments: &mut [AncestrySegment],
    phased: &crate::phasing::PhasedGenotypes,
    fine_panel: &AncestryPanel,
    params: &FineResolveParams,
) {
    if fine_panel.is_empty() || fine_panel.populations.is_empty() {
        return;
    }
    // Fine-panel column index → (code, super-pop code).
    let col_super: Vec<(&str, Option<&'static str>)> = fine_panel
        .populations
        .iter()
        .map(|c| (c.as_str(), population_super(c)))
        .collect();

    // Site lookup: (contig, pos) → fine-panel site index (only well-formed sites).
    let site_idx: HashMap<(&str, i64), usize> = fine_panel
        .sites
        .iter()
        .enumerate()
        .filter(|(_, s)| s.freqs.len() == fine_panel.populations.len())
        .map(|(i, s)| ((s.contig.as_str(), s.position), i))
        .collect();

    // Phased allele lookup: (contig, side, pos) → allele 0/1.
    let mut allele: HashMap<(&str, u8, i64), u8> = HashMap::with_capacity(phased.sites.len() * 2);
    for s in &phased.sites {
        allele.insert((s.contig.as_str(), 0, s.position), s.side0);
        allele.insert((s.contig.as_str(), 1, s.position), s.side1);
    }

    // The sites, in groups by contig and in position order. Each segment can then do a binary
    // search for its own window. A scan of `phased.sites` for each segment instead is
    // O(segments × sites). At genome scale, which is hundreds of segments over about 1M sites,
    // that scan controls the time of the whole fine-resolution pass. The shape of the groups by
    // contig is the same as in `paint_local_ancestry`.
    let mut by_contig: HashMap<&str, Vec<&crate::phasing::PhasedSite>> = HashMap::new();
    for s in &phased.sites {
        by_contig.entry(s.contig.as_str()).or_default().push(s);
    }
    for sites in by_contig.values_mut() {
        sites.sort_by_key(|s| s.position);
    }

    for seg in segments.iter_mut() {
        let sp = seg.population_code.as_str();
        // The candidate fine columns in the super-population of this segment. A fine code that
        // is the same as the super code stays out, because it offers no more resolution. MEA,
        // CAS and OCE are examples.
        let candidates: Vec<usize> = col_super
            .iter()
            .enumerate()
            .filter(|(_, (code, sup))| *sup == Some(sp) && *code != sp)
            .map(|(i, _)| i)
            .collect();
        if candidates.is_empty() {
            continue;
        }

        // The sum of the ln-likelihood of each candidate, over the informative sites of the
        // segment on this side.
        let mut ll = vec![0.0f64; candidates.len()];
        let mut n = 0usize;
        let contig_sites = by_contig.get(seg.contig.as_str()).map(Vec::as_slice).unwrap_or(&[]);
        let lo = contig_sites.partition_point(|s| s.position < seg.start);
        let hi = contig_sites.partition_point(|s| s.position <= seg.end);
        for s in &contig_sites[lo..hi] {
            let Some(&si) = site_idx.get(&(seg.contig.as_str(), s.position)) else {
                continue;
            };
            let Some(&a) = allele.get(&(seg.contig.as_str(), seg.copy, s.position)) else {
                continue;
            };
            let freqs = &fine_panel.sites[si].freqs;
            n += 1;
            for (ci, &col) in candidates.iter().enumerate() {
                ll[ci] += emit_haploid_ln(a, freqs[col] as f64);
            }
        }
        if n < params.min_sites {
            continue;
        }

        // The best and the second. Accept only with a clear margin at each site.
        let mut order: Vec<usize> = (0..candidates.len()).collect();
        order.sort_by(|&a, &b| ll[b].total_cmp(&ll[a]));
        let best = order[0];
        let runner = order.get(1).map(|&i| ll[i]).unwrap_or(f64::NEG_INFINITY);
        let margin = (ll[best] - runner) / n as f64;
        if runner.is_finite() && margin < params.min_margin_per_site {
            continue;
        }
        seg.fine_population_code = Some(fine_panel.populations[candidates[best]].clone());
    }
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

/// Diploid Viterbi. It gives the MAP ancestry **pair** `(a1, a2)` at each site. The hidden state
/// is an ordered pair of ancestries, one for each of the two genome copies, and the two copies
/// are independent Markov chains. The transitions then factorize as
/// `ln_trans(a1,b1) + ln_trans(a2,b2)`, and the emission is the two-copy [`emit_diploid_ln`].
/// Returns one `(a1, a2)` for each site, at the state index `a1*k + a2`.
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
        // The best predecessor in each chain, for each target chain-state. The step factorizes,
        // so the pair step is O(k²) and not O(k⁴). For a chain value b, take the maximum over a
        // of v_chain[a] + ln_trans(a,b).
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

/// Collapse the ancestry path of one copy, which has a value at each site, into segments. A run
/// of fewer sites than `min_sites` merges into the segment before it and takes the ancestry of
/// that segment. Each segment carries the `copy` index. The `posterior` is 1.0, because this is
/// the MAP path. A posterior for each copy is a later improvement.
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
            fine_population_code: None,
        })
        .collect()
}

use navigator_domain::seq::complement_base as revcomp_base;

/// The alt-allele dosage (0/1/2) for a diploid chip call `(a1,a2)`, against the `ref_allele` and
/// `alt_allele` of a panel site. If the two alleles are not both in `{ref,alt}`, the code tries
/// once more on the **reverse-complement** of the call. The array can report the other strand. It returns `None` if the call still does not match, which is a no-call or a
/// multi-allelic mismatch. This is the small amount of strand-flip logic that the path from a
/// chip to the panel needs.
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

/// Estimate ancestry by the binomial allele-frequency likelihood of each population.
///
/// For each population, the log-likelihood is the sum of `ln P(genotype | f)` over the sites
/// with a genotype. `f` is the alt-allele frequency of that population, clamped to
/// [0.001, 0.999]. The diploid genotype probability is `(1-f)²` for hom-ref, `2f(1-f)` for het,
/// or `f²` for hom-alt. The code then takes the exponential of each likelihood against the best
/// population, and normalizes the results to percentages.
pub fn estimate_by_allele_frequency(
    genotypes: &[SiteGenotype],
    panel: &AncestryPanel,
    reference_version: &str,
) -> AncestryResult {
    // A map from (contig, position) to a dosage. The code drops a missing or no-call dosage,
    // which is a value less than 0.
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

/// Estimate the **admixture proportions** of the sample over the panel populations, by
/// supervised ADMIXTURE. The reference allele frequencies `P`, which are the panel, stay fixed.
/// The code estimates the mixture vector `Q`, which lies on the simplex and sums to 1, that
/// maximizes the genotype likelihood `∏_j P(g_j | f_j)`. The mixed alt-allele frequency at site
/// `j` is `f_j = Σ_k q_k·p_{k,j}`, and `P(g_j|f_j)` is the diploid binomial under HWE.
///
/// The frappe/ADMIXTURE EM does the fit. Each allele copy has a latent source population. The
/// E-step gives the posterior of that population, given ref or alt. The M-step estimates `q_k`
/// again, as the mean posterior.
///
/// [`estimate_by_allele_frequency`] takes the one population that fits best. This function
/// instead gives a composition that sums to 100%, which is the shape of a consumer ancestry
/// report.
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
    // The informative sites: (dosage 0/1/2, the clamped alt frequency of each population).
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
        // Run the EM until it converges. It is monotone in the likelihood, and it costs only
        // O(sites·k) in each iteration.
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

/// Fine-population admixture. It is the same supervised EM as [`estimate_admixture`], over a
/// curated **modern subset** of a large fine-frequency panel. The `freq_global` asset holds all
/// of the reference populations, and that includes the ancient ones. A flat 173-way EM is
/// ill-posed, so the code keeps only `modern_codes`.
///
/// It uses the *same* genotypes, because the fine panel and the AIM panel share their sites. The
/// result carries the label `FINE_ADMIXTURE`. Its components roll up to the super-populations
/// through the `population_super` map in the domain crate.
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

/// The `ANCIENT_ADMIXTURE` method label. It marks deep, pre-historic source proportions.
pub const ANCIENT_ADMIXTURE: &str = "ANCIENT_ADMIXTURE";

/// Below this many genotyped panel sites the three-way fit is too noisy to report at all.
const ANCIENT_MIN_SITES: usize = 500;

/// The dispersion above which the sample is **outside the span of the ancient sources**, and the
/// code reports nothing. Under a correct model the dispersion is about 1 by construction. See
/// [`ancient_dispersion`].
///
/// The calibration ran on simulated reference individuals
/// (`panelbuild validate-ancient`). The worst case for each population was:
/// GBR 1.65 · CEU 1.58 · FIN 1.78 · TSI 2.38 · **IBS 3.65** ‖ CHB 13.1 · JPT 12.4 ·
/// YRI 175 · LWK 158. 4.0 sits in the wide, empty gap between "every European individual"
/// and "the closest East Asian". It is the middle of a real gap, and it is not a value that
/// somebody chose because it looked correct.
///
/// This threshold does **not** try to separate South Asians, at PJL 3.3 to 4.0, who lie in the
/// same range as the European tail. No dispersion threshold can do that. That is why there is a
/// second guard, [`ANCIENT_MIN_WEST_EURASIAN`].
const ANCIENT_MAX_DISPERSION: f64 = 4.0;

/// The smallest European share, from the modern super-population admixture, at which the deep
/// three-way model applies at all.
///
/// WHG / Anatolian Farmer / Steppe is a **West-Eurasian** model. Those three sources are the
/// ones that compose modern Europeans. The model has no term for Ancestral South Indian, no term
/// for East Asian, and no term for Sub-Saharan African. For a person who carries much of any of
/// those, a three-way decomposition of the *whole genome* is not an approximation. It is a
/// category error.
///
/// A Punjabi fits at Steppe 67% here. The real Steppe ancestry of that person is nearer to 20 or
/// 30%, and the remainder is Iranian-Neolithic and AASI, which this model can not see. The model
/// then puts all of the ancestry that it can not explain onto the source that is least unlike
/// it.
///
/// The dispersion alone can not catch that, because South Asians lie in the same range as the
/// European tail. But the *modern* estimate separates them cleanly, and that estimate is well
/// checked and independent of this panel. Deep ancestry runs only for a sample that the modern
/// model already calls mostly European.
const ANCIENT_MIN_WEST_EURASIAN: f64 = 50.0;

/// qpAdm model-fit acceptance: report the deep breakdown only when the model is **not rejected** at
/// this tail probability (documents/design/ancient-ancestry-rebuild.md §7.14). The garbage fits the gate
/// exists to suppress reject at p ≈ 1e-13; a real British WGS/chip accepts at p ≈ 0.15–0.21.
const QPADM_MIN_P: f64 = 0.05;
/// The tolerance of the check that the weights are correct proportions. A fit that needs a
/// source weight outside `[0,1]` shows that the model does not hold. It is not a small numerical
/// overshoot.
const QPADM_WEIGHT_TOL: f64 = 0.02;

/// Estimate **deep ancestral (ancient) source proportions**. This is the decomposition into
/// Western Hunter-Gatherer, Anatolian Farmer and Steppe pastoralist. It uses the same supervised
/// allele-frequency admixture EM as [`estimate_admixture`], over the dedicated ancient frequency
/// panel `ancestry_freq_ancient_<build>.bin`, which `panelbuild ancient-panel` builds from the
/// AADR.
///
/// This *replaces* an earlier PCA-centroid classifier. That classifier was wrong in two ways.
/// It asked which ancient population this sample **is**, which is a membership posterior, when
/// the question is what **mixture** of ancient sources the sample is. It also ran against the
/// wrong centroids. The projection had pulled those centroids in on top of the modern European
/// cloud, so they carried no ancient signal at all.
///
/// A modern European is not a *member* of WHG. That person is a *mixture*.
///
/// Allele frequencies keep the sources truly separate, at a WHG-to-ANF Fst of about 0.07, where
/// the projected PCA did not.
///
/// `modern` is the **modern** super-population admixture of the sample, which is
/// [`estimate_admixture`] over the super-population panel. It is an independent estimate that is
/// already checked. It has one use here: to decide if this West-Eurasian model applies to this
/// person at all. See [`ANCIENT_MIN_WEST_EURASIAN`].
///
/// Returns `None` when the model does not apply. There are three such cases. The run genotyped
/// too few sites. Or the sample has too little European ancestry for a WHG/ANF/Steppe
/// decomposition to say anything. Or the fit dispersion is above
/// [`ANCIENT_MAX_DISPERSION`], which puts the ancestry of the sample outside the span of the
/// three sources. A Yoruba is not *any* mixture of them.
///
/// To report nothing is the whole point. The EM always returns *some* simplex vector. To show
/// that vector for a sample that the model can not express is the exact failure that this
/// rebuild prevents.
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

/// The European share (%) of the sample, from a modern super-population estimate. It is the
/// scope check for the deep three-way model. It reads the `EUR` rollup. It works for the
/// 5-way super-population panel, and also for a finer panel that rolls up to that one.
///
/// `SuperPopulationSummary::super_population` carries the *display name* and not the code. The
/// lookup goes through the catalog, and neither form of the name goes into the code.
pub fn west_eurasian_share(modern: &AncestryResult) -> f64 {
    let eur = population_name("EUR");
    modern
        .super_population_summary
        .iter()
        .find(|s| s.super_population == eur || s.super_population == "EUR")
        .map_or(0.0, |s| s.percentage)
}

/// The ancient mixture fit **without** the threshold that decides if the model applies. It is
/// the EM result with its dispersion attached as `fit_distance`, for any sample that has enough
/// genotyped sites.
///
/// [`estimate_ancient_admixture`] is this function plus the [`ANCIENT_MAX_DISPERSION`] gate, and
/// the app calls that one. This one exists so that the offline checker can *report* the
/// dispersion of the samples that the gate refuses. You can defend the threshold only if you can
/// see the separation that it stands on.
///
/// Do not use this function on the data of a user. Its components are the exact numbers that the
/// gate exists to hold back.
pub fn ancient_admixture_fit(
    genotypes: &[SiteGenotype],
    ancient_panel: &AncestryPanel,
    reference_version: &str,
) -> Option<AncestryResult> {
    let mut result = estimate_admixture(genotypes, ancient_panel, reference_version);
    if result.snps_with_genotype < ANCIENT_MIN_SITES {
        return None;
    }

    // Get the fitted mixture again, on the axis order of the panel. The `components` come in
    // the order of their percentage.
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

/// The **deep-ancestry estimator that the app calls**. See
/// documents/design/ancient-ancestry-rebuild.md §7.14. It fits `target = Σ wᵢ · sourcesᵢ` by
/// qpAdm f4. It returns the fit as an [`AncestryResult`] over the source components, which are
/// WHG, EEF and Steppe. It returns `None` when the deep model does not apply.
///
/// `sources` and `outgroups` are indices into `panel.populations`. That order is the committed
/// Patterson-2022 configuration: the sources first, which are WHG, EEF and Steppe, and then the
/// sister outgroups.
///
/// There are four gates, and the fit must pass all of them. The sample is West-Eurasian, which
/// `modern` and [`ANCIENT_MIN_WEST_EURASIAN`] decide. The run genotyped enough sites, which is
/// [`ANCIENT_MIN_SITES`]. The qpAdm model is **not rejected**, at `p ≥` [`QPADM_MIN_P`]. The
/// weights are proportions that can occur, within [`QPADM_WEIGHT_TOL`].
///
/// A `None` must stay a `None` all the way out to the UI and the PDS. A fit that does not apply,
/// or that the test rejects, goes out as nothing. It never goes out as percentages that look
/// sure.
///
/// This function replaces [`estimate_ancient_admixture`], which is the frequency-mixture EM.
/// That EM did not pass the stability gate between WGS and chip data. The **p-value** of the
/// model fit goes out on `fit_distance`.
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

/// How well a fitted ancient mixture `q` agrees with the data, as a **variance-ratio
/// dispersion**.
///
/// At each genotyped site the mixture predicts an alt-allele frequency `f = Σ q_k·p_k`. Under
/// the HWE assumption of the model, the observed dosage `g` then has a mean of `2f` and a
/// variance of `2f(1-f)`. The mean of `(g − 2f)² / 2f(1-f)` over the sites is then about 1
/// **when the model is correct**. It grows without a limit as the true ancestry of the sample
/// moves outside the span of the sources. The mixture must then predict frequencies that the
/// genotypes continue to deny.
///
/// It is a *ratio*, so it does not change with the panel size or with the coverage of the
/// sample. You can use it as a fixed threshold, and it is not a number that somebody tuned.
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
// `f4(A,B;C,D) = mean_site (a−b)(c−d)`, over the alt-allele frequency of each population. It is
// the primitive that qpAdm stands on, and it is robust to ascertainment. See
// documents/design/ancient-ancestry-rebuild.md §7. It is a difference of differences against the
// outgroups, and it cancels the drift that the whole set shares.
//
// It is also **unbiased from *pooled* frequencies**. The estimation noise in each of the four
// slots is independent, so the cross-terms go to zero in expectation. There is no hzcorr for
// each sample, which f2 and f3 need. The genotyped sample goes in as its own "population",
// with the frequency `dosage/2 ∈ {0, 0.5, 1}`.
//
// This module is the primitive. It gives a **vector** of f4 statistics that the code estimates
// together, with its block-jackknife covariance. A later step puts the qpAdm GLS solve (§7.2) on
// top of it.

/// The genome block size, in bp, for the block jackknife of the f-statistics. About 5 Mb is much
/// more than the LD range, so the blocks are independent for this purpose. That is the
/// assumption that the jackknife variance stands on.
pub const F4_BLOCK_BP: i64 = 5_000_000;

/// A population slot in an f-statistic: either a reference population (index into
/// [`AncestryPanel::populations`]) or the genotyped sample.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Pop {
    /// The reference population `i`. Its frequency at each site is `PanelSite::freqs[i]`.
    Ref(usize),
    /// The sample that the code places. Its frequency at each site is `dosage/2`.
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

/// A vector of f4 statistics that the code estimates together, with its block-jackknife
/// covariance. The qpAdm GLS solve reads this.
#[derive(Clone, Debug)]
pub struct F4Estimate {
    /// Full-sample f4 point estimates, parallel to the requested quartets (ADMIXTOOLS reports the
    /// full-sample estimate as the statistic; the jackknife supplies only the covariance).
    pub values: Vec<f64>,
    /// The `d×d` delete-one-block jackknife covariance of `values`. See Busing and others, 1999,
    /// for blocks that are not equal.
    pub cov: Vec<Vec<f64>>,
    /// The count of sites that count. The target has a genotype, and every population that the
    /// statistic names is present.
    pub n_sites: usize,
    /// The count of genome blocks that hold one such site or more.
    pub n_blocks: usize,
}

impl F4Estimate {
    /// Standard error of statistic `i` from the jackknife covariance diagonal.
    pub fn se(&self, i: usize) -> f64 {
        self.cov
            .get(i)
            .and_then(|r| r.get(i))
            .copied()
            .unwrap_or(0.0)
            .max(0.0)
            .sqrt()
    }
}

/// The point estimate of one `f4(a,b;c,d)` over the genotyped sites, with no covariance. It is a
/// thin wrapper over [`f4_vector`]. It returns `None` if fewer than two genome blocks hold a
/// site that counts.
pub fn f4(genotypes: &[SiteGenotype], panel: &AncestryPanel, q: Quartet, block_bp: i64) -> Option<f64> {
    f4_vector(genotypes, panel, &[q], block_bp).map(|e| e.values[0])
}

/// An f4 vector over `quartets` that the code estimates together, with the unequal-block
/// jackknife covariance of Busing and others (1999).
///
/// The code measures every statistic over the **same** set of informative sites. A site counts
/// when the target has a genotype there, and every population that the statistic names is
/// present there. That is what makes the covariance a correct joint covariance for a GLS later
/// in the chain.
///
/// It returns `None` if a quartet names a population that does not exist, or if fewer than two
/// blocks hold a site. A jackknife needs two blocks or more.
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
    // Refuse a population index that is out of range, at the start. A quartet that somebody
    // built wrongly must not panic in the middle of the scan.
    let ref_ok = |p: Pop| matches!(p, Pop::Target) || matches!(p, Pop::Ref(i) if i < k);
    if !quartets
        .iter()
        .all(|q| ref_ok(q.a) && ref_ok(q.b) && ref_ok(q.c) && ref_ok(q.d))
    {
        return None;
    }

    let dosage: HashMap<(&str, i64), i32> = genotypes
        .iter()
        .filter(|g| g.dosage >= 0)
        .map(|g| ((g.contig.as_str(), g.position), g.dosage))
        .collect();

    // Add up over the informative sites, in each genome block: Σ x, the count of sites, and the
    // totals.
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

    // The delete-one-block estimates θ̂_(j), and the weight of each block, h_j = n/m_j. See
    // Busing and others, 1999, for block sizes that are not equal. With g ≥ 2, and with no empty
    // block, n − m_j ≥ 1 and h_j > 1.
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
/// test. See [`qpadm_fit`] and documents/design/ancient-ancestry-rebuild.md §7.2.
#[derive(Clone, Debug)]
pub struct QpAdmFit {
    /// The weights over the sources, **in the order that the caller gave them**. They sum to 1.
    /// `weights[0]` is the weight of the base source, which the code gets as `1 − Σ others`.
    pub weights: Vec<f64>,
    /// Standard error of each weight (from the GLS normal-equations covariance).
    pub std_errors: Vec<f64>,
    /// The χ² of the model fit. It is the minimized GLS objective, which is the residual that
    /// the span of the sources does not explain.
    pub chi2: f64,
    /// Degrees of freedom `= (#outgroups − 1) − (#sources − 1) = #outgroups − #sources`.
    pub dof: usize,
    /// The tail probability `P(χ²_dof ≥ chi2)`. A small value, below about 0.05, **rejects** the
    /// model. The sources can then not express how the target shares alleles with the
    /// outgroups.
    pub p_value: f64,
    pub n_sites: usize,
    pub n_blocks: usize,
}

impl QpAdmFit {
    /// True when every weight is a correct proportion, within `tol` of `[0,1]`. qpAdm accepts a
    /// model only when two things hold. The test does not reject it, **and** the weights can
    /// occur.
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

/// Fit `target = Σ wᵢ · sourcesᵢ` by the qpAdm f4 method. See the design document, §7.2.
///
/// The weights come from how the target *shares alleles against the outgroups*. Those are
/// differences of differences, and they cancel drift and SNP ascertainment. The weights do not
/// come from the raw frequencies of the target. That is the property that the frequency-EM of §3
/// did not have. `sources` and `outgroups` are indices into `panel.populations`. The target comes
/// in through `genotypes`, which holds `dosage/2` at each site.
///
/// The method is this. For each left population `X ∈ {target, S₂..Sₙ}`, make the vector
/// `φ_X = [f4(X, S₁; R₁, Rⱼ)]_{j=2..m}`. The admixture identity is then
/// `φ_target = Σ_{i≥2} wᵢ φ_{Sᵢ}`. Solve for the weights by GLS against the block-jackknife
/// covariance, and reweight at each iteration. That iteration is necessary because the residual
/// covariance depends on the weights, since the code also estimates the sources. Last, read the
/// χ² and the p-value of the model fit from the weighted residual.
///
/// Returns `None` in four cases: `sources.len() < 2`; `outgroups.len() < sources.len()`; the code
/// can not make the f4 vector, because there are too few blocks or sites; or the GLS system is
/// singular.
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

    // The left populations against the base source: the target, then S₂..Sₙ. The group order in
    // the f4 vector is [target, S₂, …, Sₙ]. Each group adds `l` statistics over the outgroups
    // that are not the base.
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

/// The upper tail of the χ² distribution, `P(χ²_k ≥ x)`, through the regularized upper
/// incomplete gamma `Q(k/2, x/2)`. The qpAdm model-fit p-value uses this.
fn chi2_sf(x: f64, k: usize) -> f64 {
    if k == 0 {
        return if x <= 0.0 { 1.0 } else { 0.0 };
    }
    if x <= 0.0 {
        return 1.0;
    }
    gammq(k as f64 / 2.0, x / 2.0)
}

/// `ln Γ(x)` by the Lanczos approximation (g=7), with the reflection formula for `x < 0.5`.
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

/// Build an [`AncestryResult`] from the raw probability of each population. The caller does not
/// have to normalize them. With the phase-1 super-population panel, each component *is* a
/// super-population, so the super-population summary is 1:1 with the components.
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

    // Roll the components up into super-population summaries. With a super-population panel,
    // each component is its own super-population. With a fine-grained 26-population panel, more
    // than one component goes into a single super-population.
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

    // Touch the color path of the catalog, to keep the API coherent. The UI reads the color.
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

    /// When the code drops haplotypes, every row that stays keeps its alleles and its population
    /// label.
    #[test]
    fn without_haplotypes_preserves_the_kept_rows() {
        let sites: Vec<HapSite> = (0..7)
            .map(|i| HapSite {
                contig: "chr1".to_string(),
                position: 100 + i as i64,
                reference_allele: 'A',
                alternate_allele: 'G',
            })
            .collect();
        // 4 haplotypes with distinct patterns (row h has bit set where (s + h) % 3 == 0).
        let rows: Vec<Vec<u8>> = (0..4)
            .map(|h| (0..7).map(|s| ((s + h) % 3 == 0) as u8).collect())
            .collect();
        let full = HaplotypeReference::from_rows(
            "t".to_string(),
            sites,
            vec!["GBR".to_string(), "YRI".to_string()],
            vec![0, 0, 1, 1],
            &rows,
        );

        let reduced = full.without_haplotypes(&[1, 2, 99]); // 99 is not a haplotype — ignored
        assert_eq!(reduced.n_haplotypes, 2);
        assert_eq!(reduced.n_sites, full.n_sites);
        for (new_h, old_h) in [(0usize, 0usize), (1, 3)] {
            assert_eq!(reduced.population_of(new_h), full.population_of(old_h));
            for s in 0..full.n_sites {
                assert_eq!(reduced.allele(new_h, s), full.allele(old_h, s), "hap {new_h} site {s}");
            }
        }
    }

    /// When the code thins the panel, every step-th site keeps its alleles on every haplotype.
    #[test]
    fn thin_sites_keeps_every_nth_column() {
        let sites: Vec<HapSite> = (0..9)
            .map(|i| HapSite {
                contig: "chr1".to_string(),
                position: 100 + i as i64,
                reference_allele: 'A',
                alternate_allele: 'G',
            })
            .collect();
        let rows: Vec<Vec<u8>> = (0..3)
            .map(|h| (0..9).map(|s| ((s + h) % 2 == 0) as u8).collect())
            .collect();
        let full = HaplotypeReference::from_rows("t".to_string(), sites, vec!["GBR".to_string()], vec![0, 0, 0], &rows);
        let thin = full.thin_sites(3);
        assert_eq!(thin.n_sites, 3);
        assert_eq!(thin.n_haplotypes, 3);
        for (new_s, old_s) in [(0usize, 0usize), (1, 3), (2, 6)] {
            assert_eq!(thin.sites[new_s], full.sites[old_s]);
            for h in 0..3 {
                assert_eq!(thin.allele(h, new_s), full.allele(h, old_s), "hap {h} site {new_s}");
            }
        }
        assert_eq!(full.thin_sites(1).n_sites, full.n_sites);
    }

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

    /// A 1-component PCA where the loading is +1 at every site and the panel mean is 1.0, which
    /// is a het reference. A hom-alt sample then projects to +n_sites, and a hom-ref sample to
    /// −n_sites.
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
        // A fine panel with two modern populations and one ancient population, which is Steppe.
        // The modern subset must drop the ancient column. The result carries the label
        // FINE_ADMIXTURE, and it rolls up to the super-populations.
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

    /// A sample that is HOMOZYGOUS in its ancestry. The first half is hom-alt, so both copies are
    /// A. The second half is hom-ref, so both copies are B. The diploid painting gives two
    /// copies, and each one goes from A to B at the middle point.
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
    /// painter can not express).
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

    /// The gate on the global composition. The sample is almost all A, with a short hom-ref run
    /// that paints as B when there is no gate. But B is only a 1% trace across the genome, which
    /// is below the 2% gate. The gate must drop B from the state set, so that the local painting
    /// can not show an ancestry that the genome does not contain. That was the bug that showed a
    /// 99%-European donor as East-Asian on one chromosome arm. This test also covers the k=1
    /// path.
    #[test]
    fn painting_gate_suppresses_globally_absent_ancestry() {
        let n = 80;
        let panel = two_pop_panel(n);
        // Hom-alt (→ A) everywhere except a 15-site hom-ref run (→ B) in the middle.
        let genos: Vec<SiteGenotype> = (0..n)
            .map(|i| {
                sg(
                    "chr1",
                    1 + i as i64 * 1_000_000,
                    if (40..55).contains(&i) { 0 } else { 2 },
                )
            })
            .collect();
        let prior = vec![("A".to_string(), 0.99), ("B".to_string(), 0.01)];

        // Ungated (min_ancestry 0): the hom-ref run surfaces as B.
        let ungated = paint_local_ancestry(
            &genos,
            &panel,
            &prior,
            &PaintParams {
                min_ancestry: 0.0,
                ..PaintParams::default()
            },
        );
        assert!(
            ungated.iter().any(|s| s.population_code == "B"),
            "ungated painting should surface the B run: {ungated:?}"
        );

        // Gated (default 2%): B is globally absent → dropped → the whole chromosome is A.
        let gated = paint_local_ancestry(&genos, &panel, &prior, &PaintParams::default());
        assert!(!gated.is_empty());
        assert!(
            gated.iter().all(|s| s.population_code == "A"),
            "gated painting must not invent globally-absent B: {gated:?}"
        );
    }

    /// Phased painting, with two true parental sides. Side 0 is ancestry A (alt) on the first
    /// half and B (ref) on the second. Side 1 is the mirror of that. Each side must paint as a
    /// coherent two-segment track: A to B on side 0, and B to A on side 1. That is the parent
    /// split that the unphased painter can not express.
    #[test]
    fn painting_phased_two_consistent_sides() {
        use crate::phasing::{PhasedGenotypes, PhasedSite};
        let n = 80usize;
        let panel = two_pop_panel(n);
        let sites: Vec<PhasedSite> = (0..n)
            .map(|i| {
                let first_half = i < n / 2;
                // side0: alt (A) then ref (B); side1: ref (B) then alt (A).
                let (s0, s1) = if first_half { (1u8, 0u8) } else { (0u8, 1u8) };
                PhasedSite {
                    contig: "chr1".to_string(),
                    position: 1 + i as i64 * 1_000_000,
                    side0: s0,
                    side1: s1,
                    confidence: 1.0,
                }
            })
            .collect();
        let phased = PhasedGenotypes { sites };
        let prior = vec![("A".to_string(), 0.5), ("B".to_string(), 0.5)];
        let segs = paint_local_ancestry_phased(&phased, &panel, &prior, &PaintParams::default());

        let side0: Vec<&AncestrySegment> = segs.iter().filter(|s| s.copy == 0).collect();
        let side1: Vec<&AncestrySegment> = segs.iter().filter(|s| s.copy == 1).collect();
        assert_eq!(side0.len(), 2, "side 0: {side0:?}");
        assert_eq!(side1.len(), 2, "side 1: {side1:?}");
        assert_eq!(
            (side0[0].population_code.as_str(), side0[1].population_code.as_str()),
            ("A", "B")
        );
        assert_eq!(
            (side1[0].population_code.as_str(), side1[1].population_code.as_str()),
            ("B", "A")
        );
    }

    /// Two-tier fine resolution. A EUR segment with an alt-rich phased side must resolve to the
    /// alt-rich fine population, GBR, and not to the alt-poor one, TSI. A fine
    /// code that is the same as the super code offers no more resolution, and the code skips
    /// it.
    #[test]
    fn fine_resolution_picks_alt_rich_fine_pop() {
        use crate::phasing::{PhasedGenotypes, PhasedSite};
        let n = 30usize;
        // Fine panel over two EUR fine pops: GBR alt-rich, TSI alt-poor.
        let fine_sites: Vec<PanelSite> = (0..n)
            .map(|i| PanelSite {
                contig: "chr1".to_string(),
                position: 1 + i as i64 * 1_000_000,
                reference_allele: 'A',
                alternate_allele: 'G',
                freqs: vec![0.9, 0.1], // [GBR, TSI]
            })
            .collect();
        let fine_panel = AncestryPanel {
            build: "t".into(),
            populations: vec!["GBR".into(), "TSI".into()],
            sites: fine_sites,
        };
        // Phased side 0 is alt everywhere across a single EUR segment.
        let phased = PhasedGenotypes {
            sites: (0..n)
                .map(|i| PhasedSite {
                    contig: "chr1".to_string(),
                    position: 1 + i as i64 * 1_000_000,
                    side0: 1,
                    side1: 0,
                    confidence: 1.0,
                })
                .collect(),
        };
        let mut segs = vec![AncestrySegment {
            contig: "chr1".to_string(),
            start: 1,
            end: 1 + (n as i64 - 1) * 1_000_000,
            population_code: "EUR".to_string(),
            posterior: 1.0,
            copy: 0,
            fine_population_code: None,
        }];
        resolve_fine_populations(&mut segs, &phased, &fine_panel, &FineResolveParams::default());
        assert_eq!(segs[0].fine_population_code.as_deref(), Some("GBR"));
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
    // The three-source model is the one that once gave invented numbers to a user. These tests
    // hold the two properties whose absence made that possible. The model must find a mixture
    // that nobody told it. And it must refuse a sample that its sources can not express.

    /// A deterministic LCG. The simulations below must give the same answer on every run.
    struct Lcg(u64);
    impl Lcg {
        fn next_f64(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }
        /// A diploid dosage drawn under HWE at alt-frequency `f`.
        fn dosage(&mut self, f: f64) -> i32 {
            (self.next_f64() < f) as i32 + (self.next_f64() < f) as i32
        }
    }

    /// A 3-source panel over `n` sites. The frequencies differ sharply between the sources, so
    /// the mixture is well-conditioned. The panel also holds an "outsider" frequency track, which
    /// represents a sample from outside the span of the sources. A Yoruba against WHG, ANF and
    /// Steppe is such a sample.
    ///
    /// The last site pattern is what makes the outsider truly unreachable. There, **all three
    /// sources agree** at 0.10. A mixture can predict only a value inside the convex hull of its
    /// sources, so at those sites every possible `q` predicts 0.10. The outsider carries the
    /// allele at 0.95. No mixture can absorb that, and that is the exact case that the gate must
    /// detect. Without such sites, one nearly pure source comes close enough to the outsider to
    /// pass below the threshold.
    fn ancient_panel(n: i64) -> (AncestryPanel, Vec<f64>) {
        let mut sites = Vec::new();
        let mut outsider = Vec::new();
        for pos in 1..=n {
            // Go around a set of frequency patterns that differ, so that the code can find
            // every source.
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

    /// A modern super-population estimate for a test, which is `eur`% European. It is the scope
    /// input that the deep model gates on.
    fn modern_eur(eur: f64) -> AncestryResult {
        let probs = [("EUR".to_string(), eur), ("AFR".to_string(), 100.0 - eur)];
        from_probabilities("ADMIXTURE", "aims", 1000, 1000, &probs, 0.9, "t")
    }

    /// The estimator finds a mixture that nobody gave it. The test simulates a 20/30/50
    /// individual from the source frequencies of the panel itself. The EM must then return about
    /// 20/30/50, with a dispersion near 1, which is the noise floor of the model. The
    /// PCA-centroid classifier never had this property. It answered which source the sample
    /// *is*, so a true mixture came back as one population.
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

    /// The code **refuses** a sample from outside the span of the sources. It does not decompose
    /// it. The EM always returns *some* simplex vector, and that vector is exactly what the old
    /// implementation printed as a result. The gate, and not the EM, is what makes this safe.
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

        // The raw fit still gives a breakdown that looks sure …
        let raw = ancient_admixture_fit(&genos, &panel, "t").expect("enough sites to fit");
        let total: f64 = raw.components.iter().map(|c| c.percentage).sum();
        assert!((total - 100.0).abs() < 1e-6, "the EM always returns a full simplex");
        assert!(
            raw.fit_distance.unwrap() > ANCIENT_MAX_DISPERSION,
            "an out-of-span sample must be driven above the dispersion threshold"
        );
        // … and the estimator that the app calls refuses to report it. That holds even for a
        // sample that the modern model calls European. The dispersion gate does the work here,
        // and not the scope gate.
        assert!(estimate_ancient_admixture(&genos, &panel, &modern_eur(95.0), "t").is_none());
    }

    /// With too few genotyped sites, the code gives no estimate at all, and not a noisy one.
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

    /// A WHG/ANF/Steppe decomposition is a *West-Eurasian* model. A sample that the modern
    /// estimate calls mostly non-European is out of scope, and it gets nothing. That holds even
    /// when its genotypes fit the three sources well, because a fit to the arithmetic does not
    /// make the result true.
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

        // The genotypes and the fit are the same. Only the scope is different.
        assert!(estimate_ancient_admixture(&genos, &panel, &modern_eur(95.0), "t").is_some());
        assert!(estimate_ancient_admixture(&genos, &panel, &modern_eur(20.0), "t").is_none());
    }

    // ── f-statistics core (Lever 2 / qpAdm) ─────────────────────────────────────────────────────
    //
    // Three properties test the f4 primitive against graphs with known values. They are the
    // exact f4-ratio algebra that qpAdm stands on, the antisymmetries of f4, and a *calibrated*
    // block-jackknife SE. For that last one, a symmetric tree reads f4 ≈ 0 within the noise,
    // while a real internal edge reads many SE away from zero.

    /// A panel over `pops`. Its frequency row at each site is `freqs[site][pop]`. The sites lie
    /// one in each 100 kb, so the 5 Mb block jackknife sees about 50 sites in each block. The
    /// panel also holds a target with a genotype (dosage 0) at every site, so every site counts
    /// for a quartet that names only references.
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

    /// The f4-ratio identity that qpAdm makes more general. If `X = α·P + (1−α)·Q`, which is a
    /// frequency mixture, then `f4(X,P;O1,O2) = (1−α)·f4(Q,P;O1,O2)` **exactly**. At each site
    /// the two are proportional. The ratio then gives back `1−α` for any outgroups. This is the
    /// core arithmetic that the whole method stands on. A wrong sign, or two indices in the wrong
    /// order, would push the recovered α far outside the f32 noise.
    #[test]
    fn f4_ratio_recovers_the_mixture_weight() {
        let alpha = 0.3_f64;
        let mut rng = Lcg(42);
        let freqs: Vec<Vec<f32>> = (0..4000)
            .map(|_| {
                let (r1, r2) = (rng.next_f64(), rng.next_f64());
                let p = 0.1 + 0.8 * r1;
                let q = 0.1 + 0.8 * r2;
                // Tie the outgroups to the sources, to keep f4(Q,P;O1,O2) large and clean. A
                // small denominator then does not amplify it. O1−O2 = 0.7(r1−r2) follows
                // −(Q−P), which makes f4 clearly different from zero.
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
        assert!(
            est.values[1].abs() > 0.02,
            "denominator f4 must be firmly non-degenerate"
        );
        let recovered = 1.0 - est.values[0] / est.values[1];
        assert!(
            (recovered - alpha).abs() < 1e-4,
            "f4-ratio recovered α={recovered:.6}, want {alpha}"
        );
    }

    /// The exact symmetries of f4, in pure f64 arithmetic over one fixed set of sites. An
    /// exchange of the two members of either pair negates it. An exchange of the two pairs leaves
    /// it the same.
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

    /// A symmetric tree `((A,B),(C,D))` has `f4(A,B;C,D) = 0` in expectation, because the A–B
    /// and C–D drift paths do not overlap. But `f4(A,C;B,D)` lies on the shared internal edge,
    /// and it is not zero. This test simulates that tree, and the jackknife SE must *separate*
    /// the two. The null must lie within a few SE of zero, and the real edge many SE away. Only
    /// this test shows that the covariance carries the correct scale. §5.4 needs that property,
    /// and a simulation of frequencies alone can not produce it.
    #[test]
    fn f4_jackknife_se_separates_a_null_from_a_real_edge() {
        let mut rng = Lcg(2024);
        let freqs: Vec<Vec<f32>> = (0..5000)
            .map(|_| {
                let cab = 0.5 + 0.3 * (rng.next_f64() - 0.5); // drift shared by A,B
                let ccd = 0.5 + 0.3 * (rng.next_f64() - 0.5); // drift shared by C,D
                let tip = |rng: &mut Lcg, c: f64| (c + 0.15 * (rng.next_f64() - 0.5)) as f32;
                vec![
                    tip(&mut rng, cab),
                    tip(&mut rng, cab),
                    tip(&mut rng, ccd),
                    tip(&mut rng, ccd),
                ]
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
        assert!(
            z_null.abs() < 4.0,
            "symmetric tree: f4(A,B;C,D) must sit near 0, z={z_null:.2}"
        );
        assert!(
            z_edge.abs() > 8.0,
            "real internal edge: f4(A,C;B,D) must be many SE from 0, z={z_edge:.2}"
        );
        assert!(
            est.values[0].abs() * 5.0 < est.values[1].abs(),
            "the null statistic must be far smaller than the real edge"
        );
    }

    /// The χ² upper tail against critical points from a textbook. The model-fit p-value depends
    /// on this.
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

    /// Simulate frequencies on a small admixture graph. It has three sources, S1, S2 and S3, and
    /// each one holds a distinct deep component, iA, iB or iC. It has six outgroups, and each one
    /// relates to those components in a different way. R0 is pure and near the root, R1 and R4 go
    /// to A, R2 and R5 go to B, and R3 goes to C.
    ///
    /// The graph also has a target, which is an exact frequency mixture of the three sources at
    /// each site, drawn as one diploid genome. Under Brownian drift the f4 tree identities hold,
    /// so qpAdm may decompose this graph.
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

    /// qpAdm finds the true mixture weights from how the target shares alleles against the
    /// outgroups. It accepts the model that names the correct sources. It **rejects** a model
    /// that leaves out a source that the target needs. The frequency-EM of §3 never had that last
    /// property, because it always returned a simplex that looked sure.
    #[test]
    fn qpadm_recovers_a_known_mixture_and_rejects_a_deficient_model() {
        let truth = [0.5, 0.3, 0.2];
        let (panel, genos) = qpadm_graph(20_000, truth, 20_240_717);
        let outgroups = [3usize, 4, 5, 6, 7, 8]; // R0 base, then R1..R5

        let fit = qpadm_fit(&genos, &panel, &[0, 1, 2], &outgroups, F4_BLOCK_BP).expect("3-source fit");
        assert_eq!(fit.dof, 3, "dof = #outgroups − #sources = 6 − 3");
        for (i, &want) in truth.iter().enumerate() {
            assert!(
                (fit.weights[i] - want).abs() < 0.08,
                "w{i} = {:.3}, want {want}",
                fit.weights[i]
            );
        }
        assert!(
            fit.weights_feasible(0.02),
            "weights must be valid proportions: {:?}",
            fit.weights
        );
        assert!(
            fit.p_value > 0.01,
            "well-specified model must not be rejected, p = {:.4}",
            fit.p_value
        );

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

    /// `estimate_qpadm_ancestry`, which is the function that the app calls. It reports the source
    /// weights for a European fit that names the correct sources. It gates on two things. The
    /// first is the scope, where a non-European sample gives `None`. The second is the model fit,
    /// where a model that lacks a source gives `None`.
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
            assert!(
                (pct(&r, code) - want).abs() < 9.0,
                "{code}: {:.1} vs {want}",
                pct(&r, code)
            );
        }
        assert!(
            pct(&r, "S1") > pct(&r, "S2") && pct(&r, "S2") > pct(&r, "S3"),
            "order preserved"
        );
        let p = r.fit_distance.expect("p-value on fit_distance");
        assert!((0.0..=1.0).contains(&p), "p={p}");

        // Scope gate: a mostly-non-European sample gets nothing, even though the arithmetic fits.
        assert!(estimate_qpadm_ancestry(&genos, &panel, &[0, 1, 2], &outgroups, &modern_eur(20.0), "t").is_none());
        // The gate on the model fit. The p-value rejects a model that lacks a source.
        assert!(estimate_qpadm_ancestry(&genos, &panel, &[0, 1], &outgroups, &modern_eur(95.0), "t").is_none());
    }
}
