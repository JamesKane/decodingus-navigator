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
    /// readily the mosaic jumps to a different reference haplotype. Li & Stephens motivates ≈ 4·Ne/K
    /// per Morgan (~0.25/cM for our gated European panel). Too high and the mosaic re-picks the local
    /// allele match almost per-site, discarding the long-haplotype signal that separates NW from SW
    /// European (the "southern smear"); too low and it over-commits to one haplotype's population.
    pub recomb_per_cm: f64,
    /// Ancestry (population-label) switch intensity per cM for the smoothing Viterbi. Lower → longer,
    /// more confident ancestry segments (needs sustained evidence to switch population).
    pub switch_per_cm: f64,
    /// A reference population with fewer than this many haplotypes is not a callable label; its
    /// haplotypes fold into their super-population (suppresses tiny stray-labelled reference groups).
    pub min_ref_haps: usize,
    /// Runs shorter than this many sites merge into the neighbouring segment.
    pub min_segment_sites: usize,
    /// Global-composition gate: reference haplotypes whose super-population is below this fraction of
    /// the genome-wide `prior` are dropped from the copying set entirely (the dominant super-pop is
    /// always kept). Without it the fine-grained reference invites spurious short copies from every
    /// continent — a 99%-European is painted with scattered East-Asian/African/American specks. `0.0`
    /// disables the gate.
    pub min_ancestry: f64,
}

impl Default for CopyingLaiParams {
    fn default() -> Self {
        Self {
            mismatch: 0.02,
            recomb_per_cm: 0.25,
            switch_per_cm: 0.05,
            min_ref_haps: 20,
            min_segment_sites: 20,
            min_ancestry: 0.05,
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
/// global-composition gate drops reference haplotypes from super-populations the sample doesn't have.
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

    // Global-composition gate: keep only reference haplotypes whose super-population is present
    // genome-wide (≥ min_ancestry of the prior; the dominant super-pop is always kept). This is what
    // stops a 99%-European being painted with scattered spurious continents.
    let kept_super = kept_super_pops(prior, params.min_ancestry);
    let kept_haps: Vec<usize> = (0..k)
        .filter(|&h| {
            let code = reference.population_of(h);
            let sp = population_super(code).unwrap_or(code);
            kept_super.as_ref().map_or(true, |set| set.contains(sp))
        })
        .collect();
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

    // (contig, position) → reference site column.
    let ref_col: HashMap<(&str, i64), usize> = reference
        .sites
        .iter()
        .enumerate()
        .map(|(i, s)| ((s.contig.as_str(), s.position), i))
        .collect();

    // Reference-size normalization: a population's posterior is `Σ γ` over its kept haplotypes, so a
    // population with more reference samples would win by sheer count — over-calling the largest 1KG
    // pops (Iberian/Tuscan) and burying smaller ones (British). Divide each label's posterior by its
    // kept-haplotype count so the population whose haplotypes match *best per-capita* wins.
    let mut label_count = vec![0.0f64; n_labels];
    for &l in &hap_label {
        label_count[l] += 1.0;
    }
    let label_weight: Vec<f64> = label_count.iter().map(|&c| if c > 0.0 { 1.0 / c } else { 0.0 }).collect();

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
            let mut post =
                copying_posteriors(&cols, &alleles, &positions, contig, reference, &kept_haps, map, &hap_label, n_labels, params);
            // Size-normalize per label (best per-capita match, not biggest sample).
            for row in &mut post {
                for (l, v) in row.iter_mut().enumerate() {
                    *v *= label_weight[l];
                }
            }
            let path = smooth_viterbi(&post, &positions, contig, map, n_labels, params.switch_per_cm);
            segments.extend(collapse_labels(contig, &positions, &path, &post, &labels, params.min_segment_sites, side));
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
    n_labels: usize,
    params: &CopyingLaiParams,
) -> Vec<Vec<f64>> {
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
    accumulate_post(&alpha[n - 1], &beta, hap_label, &mut post[n - 1]);
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
        accumulate_post(&alpha[i], &beta, hap_label, &mut post[i]);
    }
    post
}

/// `post_out[label] += Σ_{hap: label(hap)=label} normalized(α·β)`.
fn accumulate_post(alpha_i: &[f64], beta_i: &[f64], hap_label: &[usize], post_out: &mut [f64]) {
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
    let mut last = (0..n_labels).max_by(|&a, &b| v[n - 1][a].total_cmp(&v[n - 1][b])).unwrap_or(0);
    let mut path = vec![0usize; n];
    path[n - 1] = last;
    for i in (1..n).rev() {
        last = bp[i][last];
        path[i - 1] = last;
    }
    path
}

/// Collapse a per-site label path into segments (merging runs shorter than `min_sites` into the
/// previous run), tagging each with the side `copy`, its super-population, and — when the label is a
/// fine population — the fine code. `posterior` is the mean per-site posterior of the label.
fn collapse_labels(
    contig: &str,
    positions: &[i64],
    path: &[usize],
    post: &[Vec<f64>],
    labels: &[String],
    min_sites: usize,
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
        .map(|(l, lo, hi)| {
            let code = labels[l].as_str();
            let super_pop = population_super(code).unwrap_or(code);
            let fine = if super_pop != code { Some(code.to_string()) } else { None };
            let mean_post = (lo..=hi).map(|i| post[i].get(l).copied().unwrap_or(0.0)).sum::<f64>() / (hi - lo + 1) as f64;
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
            min_segment_sites: 5,
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

    /// A tiny reference population (below `min_ref_haps`) folds into its super-pop, so it's never a
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
        let params = CopyingLaiParams { min_segment_sites: 5, ..CopyingLaiParams::default() };
        // Empty prior → gate disabled (keep all haplotypes), so the folding path is exercised.
        let segs = paint_copying_lai(&phased, &reference, &map, &[], &params);
        // The folded tiny pop must never surface as a fine call.
        for s in segs.iter().filter(|s| s.copy == 0) {
            assert_ne!(s.fine_population_code.as_deref(), Some("TSI"));
        }
    }
}
