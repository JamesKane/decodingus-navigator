//! Haplotype-copying local ancestry inference (RFMix / Li & Stephens copying model) over the phased
//! reference-haplotype panel ([`HaplotypeReference`]).
//!
//! The allele-frequency painter ([`crate::ancestry::paint_local_ancestry_phased`]) scores each site
//! against per-population *allele frequencies*, which throws away the haplotype/linkage structure
//! that separates closely-related populations (e.g. British vs continental NW-European) — so it
//! collapses fine structure to the super-population. This module instead models each phased side as a
//! **mosaic of reference haplotypes**: a Li & Stephens copying HMM whose hidden state is *which
//! reference haplotype we are copying*, with the copied haplotype's population label read off as the
//! local ancestry. Long shared haplotypes (IBD-like tracts) pull the posterior toward the right
//! sub-population, recovering the structure frequency emission can't see.
//!
//! Per phased side and contig: forward–backward over the `K` reference haplotypes (O(N·K) thanks to
//! the rank-1 "stay or uniform-jump" transition), aggregate the per-site copy posterior by population
//! label, then a Viterbi over labels smooths the noisy per-site calls into coherent segments. Output
//! is the same [`AncestrySegment`] type the frequency painter produces, with `fine_population_code`
//! set to the resolved sub-population and `population_code` its super-population roll-up.

use std::collections::{BTreeMap, HashMap};

use navigator_domain::ancestry::{population_super, AncestrySegment};

use crate::ancestry::HaplotypeReference;
use crate::ibd::GeneticMap;
use crate::phasing::{PhasedGenotypes, PhasedSite};

/// Tuning for [`paint_copying_lai`].
#[derive(Debug, Clone)]
pub struct CopyingLaiParams {
    /// Copy mismatch rate μ: probability the copied reference allele is observed flipped (mutation /
    /// divergence since the shared ancestor). Higher tolerates more mismatch before switching copy.
    pub mismatch: f64,
    /// Reference-haplotype switch intensity per centiMorgan (the copying model's recombination): how
    /// readily the mosaic jumps to a different reference haplotype. Too high and the mosaic re-picks
    /// the local allele match almost per-site, discarding the long-haplotype signal that separates
    /// populations; too low and it commits to one reference haplotype for a whole ~10 cM tract, which
    /// at this panel's density (~0.5 markers/Mb) rests each call on a handful of sites — the regime
    /// where drifted isolates (Finnish/Sardinian/Basque) win on chance matches.
    ///
    /// Calibrated on held-out reference individuals (`navigator-panelbuild validate-lai`): over
    /// 0.1 → 1.0 accuracy rises and isolate over-call falls monotonically, flattening around 0.5.
    pub recomb_per_cm: f64,
    /// Ancestry (population-label) switch intensity per cM for the smoothing Viterbi. Lower → longer,
    /// more confident ancestry segments (needs sustained evidence to switch population).
    pub switch_per_cm: f64,
    /// A reference population with fewer than this many haplotypes is not a callable label; its
    /// haplotypes fold into their super-population (suppresses tiny stray-labelled reference groups).
    pub min_ref_haps: usize,
    /// Balance the reference by capping each population at this many haplotypes. Without any cap the
    /// largest 1000G samples (Iberian/Tuscan, ~214 haps) out-vote by count and paint a NW-European
    /// as southern; a per-size *division* over-swings to the tiny HGDP isolates. Capping puts the
    /// populations on a common footing so they compete on match quality, not sample count.
    ///
    /// The cap must stay **well above** the small reference populations (HGDP Orcadian 30, Adygei 32,
    /// Basque 46): set near them it throws away most of the big populations' haplotypes, and a
    /// thinned population copies worse, so the isolates win by default. `validate-lai` measures
    /// exactly that — at 50 a held-out individual scores 12.9% fine with 27.2% of its genome called
    /// into a drifted isolate; at 200, 24.6% and 11.8%. Above ~300 accuracy falls again (the big
    /// populations stop being balanced at all).
    ///
    /// Capping is not the only size correction — see [`Self::size_normalize`], which divides by the
    /// haplotype count and, on a dense panel, does the work capping can not.
    pub max_ref_haps: usize,
    /// Runs shorter than this many **centiMorgans** merge into the neighbouring segment.
    ///
    /// Expressed in genetic distance, not in sites, because the same site count means different
    /// things on different panels: the shipped 15.6k-site panel carries one marker per ~200 kb while
    /// the dense one carries one per ~19 kb, and `validate-lai` finds the same *physical* optimum on
    /// both (~4 cM — 20 sites on the sparse panel, 200 on the dense). A site-count threshold would
    /// silently mean a 40 Mb minimum segment if the panel got 10x denser and the number stayed put.
    pub min_segment_cm: f64,
    /// Correct the copy posterior for how many haplotypes each population contributes: divide its
    /// aggregated copy mass by `count^size_normalize` (`0.0` = off, `1.0` = a full per-haplotype
    /// average). Capping ([`Self::max_ref_haps`]) can only bring populations *down* to a common
    /// size, and this panel still spans 30 (HGDP Orcadian) to 200 haplotypes; dividing lets a small
    /// population keep all its haplotypes and still compete per-haplotype.
    ///
    /// **This knob's answer depends on marker density, and it reversed when the panel got denser.**
    /// On the old 15.6k-site panel every setting made things worse and it was rejected. On the
    /// shipped 165k-site panel `validate-lai` measures the opposite — the small HGDP populations come
    /// back (at 0.25: Sardinian 10→13% / 33→36%, Basque 33→37% / 37→39%, Orcadian 6→10%, Russian
    /// 10→12% / 4→9%) while GBR and TSI also improve and CEU stays put.
    ///
    /// The value is set by the metric the UI actually claims, not by per-site accuracy. Per-site
    /// accuracy keeps rising to 1.0 (32.0→36.4%), but the *largest called label* stops being the
    /// right one for NW-European subjects above 0.25 (top-1 50%→42%), because boosted small
    /// populations displace it — and drifted-isolate over-call climbs 9.6→11→15→18%. 0.25 is the
    /// knee: two thirds of the rescue, top-1 intact, +1.3 pts of isolate noise.
    pub size_normalize: f64,
    /// Global-composition gate: reference haplotypes whose super-population is below this fraction of
    /// the genome-wide `prior` are dropped from the copying set entirely (the dominant super-pop is
    /// always kept). Without it the fine-grained reference invites spurious short copies from every
    /// continent — a 99%-European is painted with scattered East-Asian/African/American specks. `0.0`
    /// disables the gate.
    pub min_ancestry: f64,
}

impl Default for CopyingLaiParams {
    /// Calibrated against known truth by `navigator-panelbuild validate-lai` (held-out reference
    /// individuals + simulated admixture over the shipped `ancestry_haps` panel) rather than by
    /// inspecting one kit's painting. On the dense (165k-site) panel, 41 cases spanning 1000G and
    /// HGDP populations: **32.0% fine-population accuracy against an 11.4% chance level, 68.1% at
    /// the regional level, 98.5% super-population, and 9.6%** of a non-isolate genome called into a
    /// drifted isolate. `recomb_per_cm`, `switch_per_cm` and `mismatch` are flat within
    /// case-to-case noise at this density; `min_segment_cm` is the knob that matters (4 cM beats
    /// 2 cM by 6 points of accuracy and 8 of regional accuracy).
    fn default() -> Self {
        Self {
            mismatch: 0.02,
            recomb_per_cm: 0.5,
            switch_per_cm: 0.05,
            min_ref_haps: 20,
            max_ref_haps: 200,
            min_segment_cm: 4.0,
            min_ancestry: 0.05,
            size_normalize: 0.25,
        }
    }
}

/// Super-populations to keep, given the genome-wide composition `prior` and gate `min_ancestry`:
/// those at or above the threshold, plus the dominant one. `None` (keep all) when the prior is
/// empty/degenerate.
fn kept_super_pops(prior: &[(String, f64)], min_ancestry: f64) -> Option<std::collections::HashSet<String>> {
    let mut comp: BTreeMap<String, f64> = BTreeMap::new();
    for (code, w) in prior {
        let sp = population_super(code).unwrap_or(code);
        *comp.entry(sp.to_string()).or_default() += w.max(0.0);
    }
    let total: f64 = comp.values().sum();
    if total <= 0.0 {
        return None;
    }
    let mut set: std::collections::HashSet<String> = comp
        .iter()
        .filter(|(_, &w)| w / total >= min_ancestry)
        .map(|(k, _)| k.clone())
        .collect();
    if let Some((a, _)) = comp.iter().max_by(|x, y| x.1.total_cmp(y.1)) {
        set.insert(a.clone());
    }
    Some(set)
}

/// Paint local ancestry from **phased** genotypes by haplotype copying against `reference`. Each of
/// the two phased sides is painted independently (segment `copy` = side 0/1); segments carry the
/// resolved fine population in [`AncestrySegment::fine_population_code`] and its super-population in
/// [`AncestrySegment::population_code`]. `prior` is the genome-wide composition `(pop, weight)` — the
/// global-composition gate drops reference haplotypes from super-populations the sample does not have.
/// Returns empty if the reference is empty.
pub fn paint_copying_lai(
    phased: &PhasedGenotypes,
    reference: &HaplotypeReference,
    map: &GeneticMap,
    prior: &[(String, f64)],
    params: &CopyingLaiParams,
) -> Vec<AncestrySegment> {
    let k = reference.n_haplotypes;
    if reference.is_empty() || k == 0 {
        return Vec::new();
    }

    // Global-composition gate (keep only super-populations present genome-wide ≥ min_ancestry; the
    // dominant one is always kept) + per-population capping (≤ max_ref_haps each, balancing the panel
    // so the largest 1000G samples do not out-vote by count). Gating stops spurious continents; capping
    // stops the southern-European / large-sample skew.
    let kept_super = kept_super_pops(prior, params.min_ancestry);
    let mut pop_used = vec![0usize; reference.populations.len()];
    let mut kept_haps: Vec<usize> = Vec::new();
    for h in 0..k {
        let code = reference.population_of(h);
        let sp = population_super(code).unwrap_or(code);
        if !kept_super.as_ref().map_or(true, |set| set.contains(sp)) {
            continue;
        }
        let p = reference.hap_pop[h] as usize;
        if pop_used[p] >= params.max_ref_haps {
            continue;
        }
        pop_used[p] += 1;
        kept_haps.push(h);
    }
    if kept_haps.is_empty() {
        return Vec::new();
    }

    // Callable label per KEPT haplotype (aligned to `kept_haps`): its fine population when that
    // population has enough kept haplotypes, else the super-population it rolls up to (folding tiny
    // stray-labelled groups).
    let mut pop_counts = vec![0usize; reference.populations.len()];
    for &h in &kept_haps {
        pop_counts[reference.hap_pop[h] as usize] += 1;
    }
    let mut label_index: BTreeMap<String, usize> = BTreeMap::new();
    let mut labels: Vec<String> = Vec::new();
    let mut hap_label = vec![0usize; kept_haps.len()];
    for (j, &h) in kept_haps.iter().enumerate() {
        let fine = reference.population_of(h);
        let code = if pop_counts[reference.hap_pop[h] as usize] >= params.min_ref_haps {
            fine.to_string()
        } else {
            population_super(fine).unwrap_or(fine).to_string()
        };
        let next = labels.len();
        hap_label[j] = *label_index.entry(code.clone()).or_insert_with(|| {
            labels.push(code);
            next
        });
    }
    let n_labels = labels.len();
    // Per-label divisor for the size correction: (kept haplotypes carrying the label)^size_normalize.
    let mut label_weight = vec![1.0f64; n_labels];
    if params.size_normalize != 0.0 {
        let mut counts = vec![0usize; n_labels];
        for &l in &hap_label {
            counts[l] += 1;
        }
        for (w, c) in label_weight.iter_mut().zip(counts) {
            *w = 1.0 / (c.max(1) as f64).powf(params.size_normalize);
        }
    }

    // (contig, position) → reference site column.
    let ref_col: HashMap<(&str, i64), usize> = reference
        .sites
        .iter()
        .enumerate()
        .map(|(i, s)| ((s.contig.as_str(), s.position), i))
        .collect();

    // Group phased sites by contig.
    let mut by_contig: BTreeMap<&str, Vec<&PhasedSite>> = BTreeMap::new();
    for s in &phased.sites {
        by_contig.entry(s.contig.as_str()).or_default().push(s);
    }

    let mut segments = Vec::new();
    for (contig, mut psites) in by_contig {
        psites.sort_by_key(|s| s.position);
        for side in [0u8, 1u8] {
            // This side's aligned (ref column, observed allele) and positions, in order.
            let mut cols: Vec<usize> = Vec::new();
            let mut alleles: Vec<u8> = Vec::new();
            let mut positions: Vec<i64> = Vec::new();
            for ps in &psites {
                if let Some(&c) = ref_col.get(&(contig, ps.position)) {
                    cols.push(c);
                    alleles.push(if side == 0 { ps.side0 } else { ps.side1 });
                    positions.push(ps.position);
                }
            }
            if cols.is_empty() {
                continue;
            }
            let post = copying_posteriors(
                &cols,
                &alleles,
                &positions,
                contig,
                reference,
                &kept_haps,
                map,
                &hap_label,
                &label_weight,
                params,
            );
            let path = smooth_viterbi(&post, &positions, contig, map, n_labels, params.switch_per_cm);
            segments.extend(collapse_labels(
                contig,
                &positions,
                &path,
                &post,
                &labels,
                map,
                params.min_segment_cm,
                side,
            ));
        }
    }
    segments
}

/// Forward–backward over the reference haplotypes for one side of one contig. Returns the per-site
/// posterior aggregated by population label: `post[site][label]`. Uses per-site normalization for
/// numerical stability (which cancels in the posterior).
#[allow(clippy::too_many_arguments)]
fn copying_posteriors(
    cols: &[usize],
    alleles: &[u8],
    positions: &[i64],
    contig: &str,
    reference: &HaplotypeReference,
    haps: &[usize],
    map: &GeneticMap,
    hap_label: &[usize],
    label_weight: &[f64],
    params: &CopyingLaiParams,
) -> Vec<Vec<f64>> {
    let n_labels = label_weight.len();
    let n = cols.len();
    let k = haps.len(); // the gated reference-haplotype set
    let (m_match, m_mis) = (1.0 - params.mismatch, params.mismatch);
    // Emission for kept-haplotype `j` at site `i`: match vs mismatch to the observed allele.
    let emit = |i: usize, j: usize| -> f64 {
        if reference.allele(haps[j], cols[i]) == alleles[i] {
            m_match
        } else {
            m_mis
        }
    };
    let rho = |i0: usize, i1: usize| -> f64 {
        let d = map
            .interval_cm(contig, positions[i0] as i32, positions[i1] as i32)
            .unwrap_or(0.0)
            .max(0.0);
        (1.0 - (-d * params.recomb_per_cm).exp()).clamp(1e-6, 0.9999)
    };

    // Forward (store the normalized trellis).
    let mut alpha = vec![vec![0.0f64; k]; n];
    let mut s0 = 0.0;
    for (hap, a) in alpha[0].iter_mut().enumerate() {
        let v = emit(0, hap) / k as f64;
        *a = v;
        s0 += v;
    }
    if s0 > 0.0 {
        alpha[0].iter_mut().for_each(|v| *v /= s0);
    }
    for i in 1..n {
        let r = rho(i - 1, i);
        let jump = r / k as f64; // prev row sums to 1 after normalization
        let stay = 1.0 - r;
        let (prev, cur) = alpha.split_at_mut(i);
        let prev_row = &prev[i - 1];
        let cur_row = &mut cur[0];
        let mut s = 0.0;
        for (hap, (ca, pa)) in cur_row.iter_mut().zip(prev_row.iter()).enumerate() {
            let v = emit(i, hap) * (stay * pa + jump);
            *ca = v;
            s += v;
        }
        if s > 0.0 {
            cur_row.iter_mut().for_each(|v| *v /= s);
        }
    }

    // Backward, accumulating the label-aggregated posterior as we go (no stored β trellis).
    let mut post = vec![vec![0.0f64; n_labels]; n];
    let mut beta = vec![1.0f64 / k as f64; k];
    accumulate_post(&alpha[n - 1], &beta, hap_label, label_weight, &mut post[n - 1]);
    for i in (0..n - 1).rev() {
        let r = rho(i, i + 1);
        let stay = 1.0 - r;
        let mut eb = vec![0.0f64; k];
        let mut eb_sum = 0.0;
        for (hap, (e, b)) in eb.iter_mut().zip(beta.iter()).enumerate() {
            let v = emit(i + 1, hap) * b;
            *e = v;
            eb_sum += v;
        }
        let jump = r * eb_sum / k as f64;
        let mut nb = vec![0.0f64; k];
        let mut s = 0.0;
        for (nbv, ebv) in nb.iter_mut().zip(eb.iter()) {
            let v = stay * ebv + jump;
            *nbv = v;
            s += v;
        }
        if s > 0.0 {
            nb.iter_mut().for_each(|v| *v /= s);
        }
        beta = nb;
        accumulate_post(&alpha[i], &beta, hap_label, label_weight, &mut post[i]);
    }
    post
}

/// `post_out[label] = Σ_{hap: label(hap)=label} normalized(α·β) · label_weight[label]`, renormalized
/// to sum to 1 so the size correction rescales labels against each other without changing the scale
/// the smoothing Viterbi sees.
fn accumulate_post(alpha_i: &[f64], beta_i: &[f64], hap_label: &[usize], label_weight: &[f64], post_out: &mut [f64]) {
    let mut gs = 0.0;
    for (a, b) in alpha_i.iter().zip(beta_i) {
        gs += a * b;
    }
    if gs <= 0.0 {
        return;
    }
    for ((a, b), &l) in alpha_i.iter().zip(beta_i).zip(hap_label) {
        post_out[l] += a * b / gs;
    }
    let mut total = 0.0;
    for (p, w) in post_out.iter_mut().zip(label_weight) {
        *p *= w;
        total += *p;
    }
    if total > 0.0 {
        post_out.iter_mut().for_each(|p| *p /= total);
    }
}

/// Viterbi over population labels with the per-site posterior as (log) emission and a
/// distance-scaled stay/switch transition — smooths noisy per-site calls into coherent segments.
fn smooth_viterbi(
    post: &[Vec<f64>],
    positions: &[i64],
    contig: &str,
    map: &GeneticMap,
    n_labels: usize,
    switch_per_cm: f64,
) -> Vec<usize> {
    let n = post.len();
    if n == 0 {
        return Vec::new();
    }
    let ln = |x: f64| x.max(1e-300).ln();
    let mut v = vec![vec![f64::NEG_INFINITY; n_labels]; n];
    let mut bp = vec![vec![0usize; n_labels]; n];
    for l in 0..n_labels {
        v[0][l] = ln(post[0][l]);
    }
    for i in 1..n {
        let d = map
            .interval_cm(contig, positions[i - 1] as i32, positions[i] as i32)
            .unwrap_or(0.0)
            .max(0.0);
        let sw = (1.0 - (-d * switch_per_cm).exp()).clamp(1e-6, 0.5);
        let stay = (1.0 - sw).ln();
        let jump = (sw / n_labels as f64).max(1e-300).ln();
        for b in 0..n_labels {
            let (mut best, mut arg) = (f64::NEG_INFINITY, 0usize);
            for (a, &va) in v[i - 1].iter().enumerate() {
                let val = va + if a == b { stay } else { jump };
                if val > best {
                    best = val;
                    arg = a;
                }
            }
            v[i][b] = best + ln(post[i][b]);
            bp[i][b] = arg;
        }
    }
    let mut last = (0..n_labels)
        .max_by(|&a, &b| v[n - 1][a].total_cmp(&v[n - 1][b]))
        .unwrap_or(0);
    let mut path = vec![0usize; n];
    path[n - 1] = last;
    for i in (1..n).rev() {
        last = bp[i][last];
        path[i - 1] = last;
    }
    path
}

/// Collapse a per-site label path into segments (merging runs shorter than `min_cm` of genetic
/// distance into the previous run), tagging each with the side `copy`, its super-population, and —
/// when the label is a fine population — the fine code. `posterior` is the mean per-site posterior
/// of the label.
#[allow(clippy::too_many_arguments)]
fn collapse_labels(
    contig: &str,
    positions: &[i64],
    path: &[usize],
    post: &[Vec<f64>],
    labels: &[String],
    map: &GeneticMap,
    min_cm: f64,
    copy: u8,
) -> Vec<AncestrySegment> {
    // Runs of equal label: (label, first_idx, last_idx).
    let mut runs: Vec<(usize, usize, usize)> = Vec::new();
    for (i, &l) in path.iter().enumerate() {
        match runs.last_mut() {
            Some(r) if r.0 == l => r.2 = i,
            _ => runs.push((l, i, i)),
        }
    }
    // A run's genetic span; the map is the same one the copying model used, so a missing interval
    // (contig absent from the map) falls back to 0 and the run merges — the conservative direction.
    let span_cm = |lo: usize, hi: usize| {
        map.interval_cm(contig, positions[lo] as i32, positions[hi] as i32)
            .unwrap_or(0.0)
            .max(0.0)
    };
    let mut merged: Vec<(usize, usize, usize)> = Vec::new();
    for r in runs {
        if span_cm(r.1, r.2) < min_cm {
            if let Some(prev) = merged.last_mut() {
                prev.2 = r.2;
                continue;
            }
        }
        merged.push(r);
    }
    merged
        .into_iter()
        .map(|(l, lo, hi)| {
            let code = labels[l].as_str();
            let super_pop = population_super(code).unwrap_or(code);
            let fine = if super_pop != code {
                Some(code.to_string())
            } else {
                None
            };
            let mean_post =
                (lo..=hi).map(|i| post[i].get(l).copied().unwrap_or(0.0)).sum::<f64>() / (hi - lo + 1) as f64;
            AncestrySegment {
                contig: contig.to_string(),
                start: positions[lo],
                end: positions[hi],
                population_code: super_pop.to_string(),
                posterior: mean_post,
                copy,
                fine_population_code: fine,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ancestry::HapSite;

    /// Deterministic xorshift64* — these are numeric gates, so the same run must give the same
    /// numbers twice.
    struct Rng(u64);

    impl Rng {
        fn next_f64(&mut self) -> f64 {
            let mut x = self.0;
            x ^= x >> 12;
            x ^= x << 25;
            x ^= x >> 27;
            self.0 = x;
            ((x.wrapping_mul(0x2545_F491_4F6C_DD1D)) >> 11) as f64 / (1u64 << 53) as f64
        }
        /// Standard normal (Box–Muller).
        fn normal(&mut self) -> f64 {
            let (u1, u2) = (self.next_f64().max(1e-12), self.next_f64());
            (-2.0 * u1.ln()).sqrt() * (std::f64::consts::TAU * u2).cos()
        }
        /// Allele frequency after drift of magnitude `fst` from `p` (normal approximation to the
        /// Balding–Nichols draw — good enough to give populations distinguishable frequencies).
        fn drift(&mut self, p: f64, fst: f64) -> f64 {
            (p + self.normal() * (fst * p * (1.0 - p)).sqrt()).clamp(0.01, 0.99)
        }
        fn allele(&mut self, f: f64) -> u8 {
            (self.next_f64() < f) as u8
        }
    }

    fn sites(n: usize) -> Vec<HapSite> {
        (0..n)
            .map(|i| HapSite {
                contig: "chr1".to_string(),
                position: 1 + i as i64 * 1_000_000,
                reference_allele: 'A',
                alternate_allele: 'G',
            })
            .collect()
    }

    fn phased_side(alleles: &[u8]) -> PhasedGenotypes {
        PhasedGenotypes {
            sites: alleles
                .iter()
                .enumerate()
                .map(|(i, &a)| PhasedSite {
                    contig: "chr1".to_string(),
                    position: 1 + i as i64 * 1_000_000,
                    side0: a,
                    side1: 0, // side 1 unused in these single-side tests
                    confidence: 1.0,
                })
                .collect(),
        }
    }

    /// Two distinguishable reference populations (GBR→EUR carrying pattern G, YRI→AFR carrying the
    /// opposite pattern). A side that matches GBR haplotypes on the first half and YRI on the second
    /// must paint GBR (EUR) then YRI (AFR) with the fine codes resolved.
    #[test]
    fn copies_the_matching_population_and_resolves_fine() {
        let n = 60usize;
        // Pattern G: 1 at even sites; Pattern Y: its complement.
        let g: Vec<u8> = (0..n).map(|i| (i % 2 == 0) as u8).collect();
        let y: Vec<u8> = (0..n).map(|i| (i % 2 == 1) as u8).collect();
        // 4 GBR haps (pattern G) + 4 YRI haps (pattern Y).
        let mut rows: Vec<Vec<u8>> = Vec::new();
        let mut hap_pop: Vec<u16> = Vec::new();
        for _ in 0..4 {
            rows.push(g.clone());
            hap_pop.push(0);
        }
        for _ in 0..4 {
            rows.push(y.clone());
            hap_pop.push(1);
        }
        let reference = HaplotypeReference::from_rows(
            "t".to_string(),
            sites(n),
            vec!["GBR".to_string(), "YRI".to_string()],
            hap_pop,
            &rows,
        );
        let map = GeneticMap::uniform(1.0, &[("chr1", 250_000_000)]);
        // Side: pattern G first half, pattern Y second half.
        let mut side = g.clone();
        side[n / 2..].copy_from_slice(&y[n / 2..]);
        let phased = phased_side(&side);

        let params = CopyingLaiParams {
            min_ref_haps: 2,
            min_segment_cm: 5.0,
            ..CopyingLaiParams::default()
        };
        // Prior keeps both continents.
        let prior = vec![("EUR".to_string(), 0.5), ("AFR".to_string(), 0.5)];
        let segs = paint_copying_lai(&phased, &reference, &map, &prior, &params);
        let side0: Vec<&AncestrySegment> = segs.iter().filter(|s| s.copy == 0).collect();
        assert_eq!(side0.len(), 2, "expected GBR→YRI switch: {side0:?}");
        assert_eq!(side0[0].fine_population_code.as_deref(), Some("GBR"));
        assert_eq!(side0[0].population_code, "EUR");
        assert_eq!(side0[1].fine_population_code.as_deref(), Some("YRI"));
        assert_eq!(side0[1].population_code, "AFR");

        // Global-composition gate: an EUR-only prior drops the AFR (YRI) reference haplotypes, so even
        // the YRI-matching second half can no longer be painted AFR — the whole side stays European.
        let segs_gated = paint_copying_lai(&phased, &reference, &map, &[("EUR".to_string(), 1.0)], &params);
        assert!(
            segs_gated.iter().all(|s| s.population_code == "EUR"),
            "EUR-only prior must gate out AFR: {segs_gated:?}"
        );
    }

    /// A tiny reference population (below `min_ref_haps`) folds into its super-pop, so it is never a
    /// callable fine label: a 1-hap TSI group is never painted "TSI" (it folds to EUR).
    #[test]
    fn tiny_populations_fold_into_super_pop() {
        let n = 30usize;
        let pat: Vec<u8> = (0..n).map(|i| (i % 2 == 0) as u8).collect();
        // 20 GBR haps + 1 TSI hap, same pattern; TSI is sub-threshold so it folds to its super-pop.
        let mut rows = vec![pat.clone(); 20];
        let mut hap_pop = vec![0u16; 20];
        rows.push(pat.clone());
        hap_pop.push(1);
        let reference = HaplotypeReference::from_rows(
            "t".to_string(),
            sites(n),
            vec!["GBR".to_string(), "TSI".to_string()],
            hap_pop,
            &rows,
        );
        let map = GeneticMap::uniform(1.0, &[("chr1", 250_000_000)]);
        let phased = phased_side(&pat);
        let params = CopyingLaiParams {
            min_segment_cm: 5.0,
            ..CopyingLaiParams::default()
        };
        // Empty prior → gate disabled (keep all haplotypes), so the folding path is exercised.
        let segs = paint_copying_lai(&phased, &reference, &map, &[], &params);
        // The folded tiny pop must never surface as a fine call.
        for s in segs.iter().filter(|s| s.copy == 0) {
            assert_ne!(s.fine_population_code.as_deref(), Some("TSI"));
        }
    }

    /// A synthetic three-population reference with realistic *drift* structure: two close European
    /// populations (GBR, TSI) split off a common European branch, plus FIN — a small, heavily
    /// drifted isolate whose low internal diversity is what made it over-attract the copy in the
    /// field. Haplotypes are drawn from each population's own frequencies.
    ///
    /// Returns `(reference, sides)` where `sides` are two haplotypes drawn from GBR's frequencies
    /// and **not** in the reference (so nothing self-copies).
    fn drifted_reference(n_sites: usize) -> (HaplotypeReference, [Vec<u8>; 2]) {
        let mut rng = Rng(0x5EED_1234_5678_9ABC);
        // Ancestral → European branch → {GBR, TSI}; FIN drifts harder off the same branch. AFR
        // (YRI) splits at the root. The sibling separation here (Fst 0.02) is deliberately wider
        // than real intra-European Fst (~0.005) and the sites are independent: this gate tests the
        // *model*, given information enough to separate the populations. Whether the shipped
        // panel carries that much information is the separate question `validate-lai` answers on
        // the real reference.
        let mut freqs: Vec<[f64; 4]> = Vec::with_capacity(n_sites); // GBR, TSI, FIN, YRI
        for _ in 0..n_sites {
            let anc = 0.05 + rng.next_f64() * 0.9;
            let eur = rng.drift(anc, 0.02);
            freqs.push([
                rng.drift(eur, 0.020),
                rng.drift(eur, 0.020),
                rng.drift(eur, 0.060),
                rng.drift(anc, 0.150),
            ]);
        }
        // Panel sizes deliberately unequal (the balance the capping knob exists to fix): a big GBR
        // and TSI, a small FIN, a big YRI.
        let counts = [60usize, 60, 24, 60];
        let mut rows: Vec<Vec<u8>> = Vec::new();
        let mut hap_pop: Vec<u16> = Vec::new();
        for (p, &count) in counts.iter().enumerate() {
            for _ in 0..count {
                rows.push((0..n_sites).map(|s| rng.allele(freqs[s][p])).collect());
                hap_pop.push(p as u16);
            }
        }
        let reference = HaplotypeReference::from_rows(
            "t".to_string(),
            sites(n_sites),
            vec![
                "GBR".to_string(),
                "TSI".to_string(),
                "FIN".to_string(),
                "YRI".to_string(),
            ],
            hap_pop,
            &rows,
        );
        let side = |rng: &mut Rng| (0..n_sites).map(|s| rng.allele(freqs[s][0])).collect::<Vec<u8>>();
        let sides = [side(&mut rng), side(&mut rng)];
        (reference, sides)
    }

    /// Per-site share of each called fine label across both sides, keyed by label.
    fn call_shares(segs: &[AncestrySegment], positions: &[i64]) -> BTreeMap<String, f64> {
        let mut shares: BTreeMap<String, f64> = BTreeMap::new();
        let mut total = 0.0;
        for seg in segs {
            let n = positions.iter().filter(|p| **p >= seg.start && **p <= seg.end).count() as f64;
            let label = seg
                .fine_population_code
                .clone()
                .unwrap_or_else(|| seg.population_code.clone());
            *shares.entry(label).or_default() += n;
            total += n;
        }
        if total > 0.0 {
            shares.values_mut().for_each(|v| *v /= total);
        }
        shares
    }

    /// **Numeric gate** (the property the calibration commits were chasing, as an assertion rather
    /// than a look at the painted chromosomes): at the shipped defaults a GBR individual painted
    /// against a reference holding a small, heavily drifted FIN must be called GBR — over the
    /// isolate, over its sibling population, above chance, and never outside its continent. The
    /// end-to-end equivalent on the real panel is `navigator-panelbuild validate-lai`, which scores
    /// held-out reference individuals against the same properties.
    #[test]
    fn drifted_isolate_does_not_out_call_the_true_population() {
        let n = 600usize;
        let (reference, sides) = drifted_reference(n);
        let positions: Vec<i64> = (0..n).map(|i| 1 + i as i64 * 1_000_000).collect();
        let phased = PhasedGenotypes {
            sites: (0..n)
                .map(|i| PhasedSite {
                    contig: "chr1".to_string(),
                    position: positions[i],
                    side0: sides[0][i],
                    side1: sides[1][i],
                    confidence: 1.0,
                })
                .collect(),
        };
        let map = GeneticMap::uniform(1.0, &[("chr1", 700_000_000)]);
        let prior = vec![("EUR".to_string(), 0.99), ("AFR".to_string(), 0.01)];
        let segs = paint_copying_lai(&phased, &reference, &map, &prior, &CopyingLaiParams::default());
        let shares = call_shares(&segs, &positions);
        let share = |code: &str| shares.get(code).copied().unwrap_or(0.0);

        // 1. The drifted isolate barely features at all.
        assert!(share("FIN") < 0.05, "drifted FIN over-attracted the copy: {shares:?}");
        // 2. The true population beats its sibling clearly, not by a nose.
        assert!(
            share("GBR") > share("TSI") * 1.5,
            "sibling TSI not separated from the true GBR: {shares:?}"
        );
        // 3. It clears chance (3 callable European labels here → 1/3) by a wide margin.
        assert!(share("GBR") > 0.6, "GBR calls too close to chance: {shares:?}");
        // 4. The continent is never wrong: a 99%-European prior gates AFR out entirely.
        assert!(
            segs.iter().all(|s| s.population_code == "EUR"),
            "painted outside the sample's continent: {shares:?}"
        );
    }

    /// **Calibration regression gate.** The previous defaults (`recomb_per_cm` 0.1, `max_ref_haps`
    /// 50) were tuned by eye and are measurably worse: committing the mosaic to whole ~10 cM tracts
    /// rests each call on a handful of markers, and capping populations down near the size of the
    /// small HGDP isolates thins the large ones until the isolates win by default. `validate-lai`
    /// measures the same gradient on the real panel (fine accuracy 12.9% → 24.6%, isolate over-call
    /// 27.2% → 11.8%). This pins the direction, so reverting the calibration fails here.
    #[test]
    fn shipped_calibration_beats_the_previous_one() {
        let n = 600usize;
        let (reference, sides) = drifted_reference(n);
        let positions: Vec<i64> = (0..n).map(|i| 1 + i as i64 * 1_000_000).collect();
        let phased = PhasedGenotypes {
            sites: (0..n)
                .map(|i| PhasedSite {
                    contig: "chr1".to_string(),
                    position: positions[i],
                    side0: sides[0][i],
                    side1: sides[1][i],
                    confidence: 1.0,
                })
                .collect(),
        };
        let map = GeneticMap::uniform(1.0, &[("chr1", 700_000_000)]);
        let prior = vec![("EUR".to_string(), 1.0)];
        let paint = |params: &CopyingLaiParams| {
            call_shares(
                &paint_copying_lai(&phased, &reference, &map, &prior, params),
                &positions,
            )
        };
        let previous = CopyingLaiParams {
            recomb_per_cm: 0.1,
            max_ref_haps: 50,
            ..CopyingLaiParams::default()
        };
        let (old, new) = (paint(&previous), paint(&CopyingLaiParams::default()));
        let share = |m: &BTreeMap<String, f64>, code: &str| m.get(code).copied().unwrap_or(0.0);
        assert!(
            share(&new, "GBR") > share(&old, "GBR"),
            "shipped calibration calls the true population less often than the old one: {new:?} vs {old:?}"
        );
        assert!(
            share(&new, "FIN") < share(&old, "FIN"),
            "shipped calibration over-calls the drifted isolate more than the old one: {new:?} vs {old:?}"
        );
    }
}
