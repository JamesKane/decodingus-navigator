//! Haplotype-copying local ancestry inference (RFMix / Li & Stephens copying model) over the phased
//! reference-haplotype panel ([`HaplotypeReference`]).
//!
//! The allele-frequency painter ([`crate::ancestry::paint_local_ancestry_phased`]) scores each site
//! against the *allele frequency* of each population. That throws away the haplotype and linkage
//! structure. Two populations that are near each other, such as the British and the continental
//! north-west European, need that structure to separate. So the painter collapses the fine
//! structure down to the super-population.
//!
//! This module instead models each phased side as a **mosaic of reference haplotypes**. It is a Li
//! & Stephens copying HMM. Its hidden state is *which reference haplotype the model copies*, and
//! the population label of that haplotype is the local ancestry. A long shared haplotype, which is
//! a tract like an IBD segment, pulls the posterior toward the correct sub-population. That
//! recovers the structure that a frequency emission can not see.
//!
//! For each phased side and each contig, the code does three things. First it runs
//! forward-backward over the `K` reference haplotypes. That costs O(N·K), because the transition
//! has rank 1: it is a "stay, or jump to a uniform choice". Then it adds up the copy posterior at
//! each site, by population label. Last, a Viterbi over the labels turns the noisy call at each
//! site into coherent segments.
//!
//! The output is the same [`AncestrySegment`] type that the frequency painter makes.
//! `fine_population_code` holds the resolved sub-population, and `population_code` holds the
//! super-population that it rolls up to.

use std::collections::{BTreeMap, HashMap};

use navigator_domain::ancestry::{population_super, AncestrySegment};

use crate::ancestry::HaplotypeReference;
use crate::ibd::GeneticMap;
use crate::phasing::{PhasedGenotypes, PhasedSite};

/// The controls of [`paint_copying_lai`].
#[derive(Debug, Clone)]
pub struct CopyingLaiParams {
    /// The copy mismatch rate μ. It is the probability that the copied reference allele reads the
    /// other way round, from a mutation or from divergence since the shared ancestor. A higher
    /// value accepts more mismatch before the model changes the haplotype that it copies.
    pub mismatch: f64,
    /// The switch intensity of the reference haplotype in one centiMorgan. It is the
    /// recombination of the copying model, and it says how easily the mosaic jumps to a different
    /// reference haplotype.
    ///
    /// Set it too high, and the mosaic takes the local allele match again at almost every site.
    /// That throws away the long-haplotype signal which separates the populations. Set it too low,
    /// and the mosaic holds one reference haplotype for a whole tract of about 10 cM. At the
    /// density of this panel, which is about 0.5 markers/Mb, each call then stands on a few sites.
    /// That is where a drifted isolate, such as the Finnish, Sardinian or Basque, wins on a match
    /// by chance.
    ///
    /// The calibration ran on held-out reference individuals, with
    /// `navigator-panelbuild validate-lai`. From 0.1 to 1.0, the accuracy rises and the over-call
    /// of an isolate falls, both without a reverse, and the curve goes flat near 0.5.
    pub recomb_per_cm: f64,
    /// The switch intensity of the ancestry, which is the population label, in one cM, for the
    /// Viterbi that smooths the path. A lower value gives longer ancestry segments with more
    /// confidence, because the model then needs evidence over a longer run to change
    /// population.
    pub switch_per_cm: f64,
    /// A reference population with fewer than this many haplotypes is not a callable label; its
    /// haplotypes fold into their super-population (suppresses tiny stray-labelled reference groups).
    pub min_ref_haps: usize,
    /// Balance the reference: hold each population to this many haplotypes at most.
    ///
    /// With no limit, the largest 1000G samples win by count alone, and they paint a north-west
    /// European as southern. Those are the Iberian and the Tuscan, at about 214 haplotypes. A
    /// *division* by the size instead swings too far, to the very small HGDP isolates. A limit
    /// puts the populations on common ground. They then compete on the quality of the match, and
    /// not on the size of the sample.
    ///
    /// The limit must stay **well above** the small reference populations: HGDP Orcadian has 30,
    /// Adygei 32, and Basque 46. Set near them, it throws away most of the haplotypes of the big
    /// populations. A population that the code thins copies worse, so the isolates then win by
    /// default.
    ///
    /// `validate-lai` measures exactly that. At 50, a held-out individual scores 12.9% fine, and
    /// 27.2% of its genome goes to a drifted isolate. At 200, those become 24.6% and 11.8%. Above
    /// about 300 the accuracy falls again, because the big populations are no longer in balance at
    /// all.
    ///
    /// This limit is not the only correction for size. See [`Self::size_normalize`], which divides
    /// by the haplotype count. On a dense panel that division does the work that this limit can
    /// not.
    pub max_ref_haps: usize,
    /// A run shorter than this many **centiMorgans** merges into the segment beside it.
    ///
    /// This is a genetic distance, and not a count of sites, because the same count of sites means
    /// different things on different panels. The 15.6k-site panel that shipped carries one marker
    /// in about 200 kb, and the dense one carries one in about 19 kb. `validate-lai` finds the same
    /// *physical* best value on both, at about 4 cM. That is 20 sites on the sparse panel and 200
    /// on the dense one.
    ///
    /// Take a threshold in sites, and make the panel 10x denser, and leave the number where it
    /// was. That threshold then means a smallest segment of 40 Mb, and nobody would see it
    /// happen.
    pub min_segment_cm: f64,
    /// Correct the copy posterior for the count of haplotypes that each population gives. The code
    /// divides the copy mass of a population by `count^size_normalize`. `0.0` turns the correction
    /// off, and `1.0` gives a full mean over the haplotypes.
    ///
    /// The limit in [`Self::max_ref_haps`] can bring a population *down* to a common size and
    /// nothing more. This panel still runs from 30 haplotypes (HGDP Orcadian) to 200. The division
    /// instead lets a small population keep all of its haplotypes, and still compete on a
    /// haplotype-for-haplotype basis.
    ///
    /// **The correct value depends on the marker density, and it went the other way when the panel
    /// got denser.** On the old 15.6k-site panel, every setting made the result worse, and this
    /// control stayed off.
    ///
    /// On the 165k-site panel that ships, `validate-lai` measures the opposite. The small HGDP
    /// populations come back. At 0.25: Sardinian 10→13% / 33→36%, Basque 33→37% / 37→39%,
    /// Orcadian 6→10%, Russian 10→12% / 4→9%. GBR and TSI also get better, and CEU does not
    /// move.
    ///
    /// The metric that the UI claims sets this value. The accuracy at each site does not set it.
    ///
    /// That accuracy continues to rise up to 1.0, from 32.0% to 36.4%. But above 0.25 the *largest
    /// called label* stops being the correct one for a north-west European subject. Top-1 falls
    /// from 50% to 42%, because the small populations that the correction lifts displace it. The
    /// over-call of a drifted isolate also climbs, 9.6→11→15→18%.
    ///
    /// 0.25 is the knee of the curve. It gives two thirds of the rescue, it leaves top-1 intact,
    /// and it costs 1.3 points of isolate noise.
    pub size_normalize: f64,
    /// The gate on the global composition. The code drops a reference haplotype from the copying
    /// set when its super-population is below this fraction of the genome-wide `prior`. It always
    /// keeps the dominant super-population.
    ///
    /// Without the gate, the fine-grained reference invites short false copies from every
    /// continent. A subject who is 99% European then gets small East-Asian, African and American
    /// marks over the whole genome. `0.0` turns the gate off.
    pub min_ancestry: f64,
}

impl Default for CopyingLaiParams {
    /// `navigator-panelbuild validate-lai` calibrated these defaults against known truth. It used
    /// held-out reference individuals, and simulated admixture over the `ancestry_haps` panel that
    /// ships. Nobody looked at the painting of one kit and chose a value.
    ///
    /// The run covered the dense 165k-site panel, over 41 cases across the 1000G and HGDP
    /// populations. The result was **32.0% accuracy at the fine population, against a chance level
    /// of 11.4%. At the regional level it was 68.1%, and at the super-population 98.5%.** And
    /// **9.6%** of a genome that is not an isolate went to a drifted isolate.
    ///
    /// At this density, `recomb_per_cm`, `switch_per_cm` and `mismatch` are flat inside the noise
    /// between one case and the next. `min_segment_cm` is the control that matters: 4 cM beats
    /// 2 cM by 6 points of accuracy, and by 8 points of regional accuracy.
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

/// Paint local ancestry from **phased** genotypes, by a copy of haplotypes against `reference`.
///
/// The code paints each of the two phased sides on its own, and the segment `copy` holds the side,
/// 0 or 1. A segment carries the resolved fine population in
/// [`AncestrySegment::fine_population_code`], and the super-population of that fine population in
/// [`AncestrySegment::population_code`].
///
/// `prior` is the genome-wide composition, as `(pop, weight)`. The gate on the global composition
/// uses it to drop a reference haplotype whose super-population the sample does not have. Returns
/// an empty result if the reference is empty.
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

    // Two steps. First comes the gate on the global composition. Keep a super-population only
    // when the genome shows it at min_ancestry or more, and always keep the dominant one. Then
    // comes the limit on each population: max_ref_haps haplotypes at most. That limit balances the
    // panel, so that the largest 1000G samples do not win by count alone. The gate stops a false
    // continent, and the limit stops the skew toward southern Europe and toward a large sample.
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

    // The label of each haplotype that the code KEPT, in line with `kept_haps`. It is the fine
    // population of that haplotype, when the code kept enough haplotypes of that population.
    // Otherwise it is the super-population that the fine one rolls up to. That puts a very small
    // group with a stray label back into a larger one.
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
    // The divisor of each label, for the size correction. It is the count of kept haplotypes with
    // that label, raised to the power size_normalize.
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

/// Forward-backward over the reference haplotypes, for one side of one contig. Returns the
/// posterior at each site, added up by population label, as `post[site][label]`. It normalizes at
/// each site to keep the arithmetic stable, and that normalization cancels in the posterior.
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

    // The backward pass. It adds up the posterior over the labels as it goes, and it stores no β
    // trellis.
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

/// `post_out[label] = Σ_{hap: label(hap)=label} normalized(α·β) · label_weight[label]`. The code
/// then normalizes the result to sum to 1. The size correction then scales the labels against each
/// other, and the scale that the Viterbi sees does not move.
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

/// Viterbi over the population labels. The emission is the log of the posterior at each site, and
/// the transition, which is a stay or a switch, scales with distance. It turns the noisy call at
/// each site into coherent segments.
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

/// Collapse a path that holds a label at each site into segments. A run shorter than `min_cm` of
/// genetic distance merges into the run before it. Each segment carries the side `copy`, the
/// super-population, and, when the label is a fine population, the fine code. `posterior` is the
/// mean posterior of the label over the sites.
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
    // The genetic span of a run. The map is the same one that the copying model used. A missing
    // interval, where the map does not hold the contig, falls back to 0, and the run merges. That
    // is the careful direction.
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

    /// A deterministic xorshift64*. These are numeric gates, so two runs must give the same
    /// numbers.
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
        /// The allele frequency after drift of size `fst` away from `p`. It is the normal
        /// approximation to the draw of Balding and Nichols (1995). That is good enough to give
        /// the populations frequencies that a test can separate.
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

    /// Two reference populations that a test can separate. GBR, which rolls up to EUR, holds
    /// pattern G. YRI, which rolls up to AFR, holds the opposite pattern. Take a side that matches
    /// the GBR haplotypes on the first half and the YRI ones on the second. The painter must give
    /// GBR (EUR) and then YRI (AFR), with both fine codes resolved.
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

        // The gate on the global composition. A prior that holds EUR alone drops the AFR (YRI)
        // reference haplotypes. Even the second half, which matches YRI, can then no longer get
        // an AFR paint. The whole side stays European.
        let segs_gated = paint_copying_lai(&phased, &reference, &map, &[("EUR".to_string(), 1.0)], &params);
        assert!(
            segs_gated.iter().all(|s| s.population_code == "EUR"),
            "EUR-only prior must gate out AFR: {segs_gated:?}"
        );
    }

    /// A very small reference population, below `min_ref_haps`, goes into its super-population.
    /// It is then never a fine label that the painter can call. A TSI group of one haplotype never
    /// gets the paint "TSI", because it goes to EUR.
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
        // An empty prior turns the gate off, and every haplotype stays. This test then covers the
        // path that puts a small group into its super-population.
        let segs = paint_copying_lai(&phased, &reference, &map, &[], &params);
        // The folded tiny pop must never surface as a fine call.
        for s in segs.iter().filter(|s| s.copy == 0) {
            assert_ne!(s.fine_population_code.as_deref(), Some("TSI"));
        }
    }

    /// A synthetic reference of three populations, with a realistic *drift* structure. Two close
    /// European populations, GBR and TSI, come off a common European branch. FIN is the third, and
    /// it is a small isolate with much drift. Its low internal diversity is what made it attract
    /// too much of the copy in the field. The code draws the haplotypes from the frequencies of
    /// each population.
    ///
    /// Returns `(reference, sides)` where `sides` are two haplotypes drawn from GBR's frequencies
    /// and **not** in the reference (so nothing self-copies).
    fn drifted_reference(n_sites: usize) -> (HaplotypeReference, [Vec<u8>; 2]) {
        let mut rng = Rng(0x5EED_1234_5678_9ABC);
        // The ancestral node goes to a European branch, and that branch goes to {GBR, TSI}. FIN
        // drifts further off the same branch. AFR (YRI) splits at the root.
        //
        // The separation between the two close populations here is an Fst of 0.02. That is wider
        // than the real Fst inside Europe, which is about 0.005, and the sites here are
        // independent. Both choices are deliberate. This gate tests the *model*, and it gives the
        // model enough information to separate the populations. Whether the panel that ships
        // carries that much information is a separate question, and `validate-lai` answers it on
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
        // The panel sizes are not equal, and that is deliberate. It is the balance that the limit
        // on each population exists to correct: a big GBR and TSI, a small FIN, and a big YRI.
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

    /// The share of each called fine label, at each site, over both sides. The key is the label.
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

    /// **A numeric gate.** It holds the property that the calibration commits went after, as an
    /// assertion, and not as a look at the painted chromosomes.
    ///
    /// At the defaults that ship, take a GBR individual, and paint it against a reference that
    /// holds a small FIN with much drift. The painter must call it GBR. GBR must beat the isolate,
    /// the population beside it on the tree, and chance. And the painter must never go outside the
    /// continent.
    ///
    /// `navigator-panelbuild validate-lai` is the end-to-end equivalent on the real panel. It
    /// scores held-out reference individuals against the same properties.
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
        // 2. The true population clearly beats the population beside it on the tree. The margin
        //    is not small.
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

    /// **A gate against a return of the old calibration.** Somebody chose the defaults before
    /// these, which were `recomb_per_cm` 0.1 and `max_ref_haps` 50, by eye. A measurement shows
    /// that they are worse.
    ///
    /// There are two reasons. To hold the mosaic on a whole tract of about 10 cM leaves each call
    /// on a few markers alone. And a limit near the size of the small HGDP isolates thins the
    /// large populations until the isolates win by default.
    ///
    /// `validate-lai` measures the same gradient on the real panel: the fine accuracy goes from
    /// 12.9% to 24.6%, and the over-call of an isolate from 27.2% to 11.8%. This test holds the
    /// direction, so a change back to the old calibration fails here.
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
