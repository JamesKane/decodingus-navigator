//! Detection of runs of homozygosity (ROH), which is autozygosity. That is the signal of endogamy
//! and of consanguinity.
//!
//! **This is a spike, from 2026-07-22.** It is a hidden Markov model with two states, Autozygous
//! and Normal, over the autosomal genotypes of a subject.
//!
//! It has the same idiom as the HMM in [`crate::ancestry::paint_local_ancestry`]. That means
//! sorted sites in each contig, and transitions that scale with distance and reset to the prior.
//! It also means a Viterbi and a forward-backward, both in log space.
//!
//! The code joins the runs of the Autozygous state into [`RohSegment`] values. It then rolls those
//! up into an [`RohSummary`], with the genome-wide coefficient F_ROH.
//!
//! **What to give it.** Give it the autosomal consensus genotypes of the subject, from
//! `consensus_genotypes(&DiploidProfile)` in `navigator-app`. The caller calls those at the full
//! 1240k IBD panel, which is a dense set of about 1.15M neutral, biallelic, common SNPs, with full
//! 0/1/2 dosages. That is the density class that a ROH tool for arrays needs, such as PLINK,
//! BCFtools/RoH or detectRUNS. The cM length of a segment, and the denominator of F_ROH, come from
//! the same [`GeneticMap`] that the IBD path already loads.
//!
//! **What this spike leaves simple, deliberately.** See the module tests and the notes that follow
//! them.
//!
//! - One `baseline_het` control holds the heterozygosity that the Normal state expects. A
//!   production version must instead derive it at each site from the panel allele frequencies, as
//!   2·f·(1−f). `AncestryPanel` and `IbdPanel` already carry those, so the emission can know the
//!   frequency.
//! - The [`RohPattern`] that separates endogamy from consanguinity is a heuristic over the
//!   distribution of the ROH length classes. It is not a calibrated model.

use crate::caller::SiteGenotype;
use crate::ibd::{normalize_chromosome, GeneticMap};
use std::collections::BTreeMap;

/// Detector configuration. Defaults target a 1240k-density common-SNP substrate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RohConfig {
    /// P(a heterozygous call | the site is truly autozygous). That is the rate of genotype error
    /// that stays. It lets one stray het live inside a homozygous run, and that run does not
    /// break. The value is small.
    pub het_error: f64,
    /// The heterozygosity that the Normal state expects: P(het | the site is not autozygous).
    ///
    /// `None` is the default, and it estimates the value from the autosomal het fraction of the
    /// sample itself, clamped. It then follows the density of the panel, and its ascertainment,
    /// and it is not a fixed guess. The production upgrade is 2·f·(1−f) at each site, from the
    /// allele frequencies. `Some(v)` fixes the value, mostly for a test, or for somebody who wants
    /// to tune it.
    pub baseline_het: Option<f64>,
    /// The hazard of a state switch in one centimorgan. The switch probability over a gap of `d`
    /// cM is `1 − exp(−d · switch_rate_per_cm)`. A smaller value gives longer runs. The default is
    /// about one switch in 13 cM.
    pub switch_rate_per_cm: f64,
    /// The stationary autozygosity fraction. It is the prior mass on the Autozygous state, and a
    /// switch resets toward it. This is the classic prior of a ROH HMM.
    pub prior_autozygous: f64,
    /// Report a run of this length or more, in **physical Mb**.
    ///
    /// PLINK, detectRUNS and the genealogy field all threshold a ROH on its physical length. ROH
    /// gather in regions of low recombination, near a centromere. There a run of some Mb covers
    /// well under one cM, so a threshold in cM reports too few of them.
    pub min_length_mb: f64,
    /// Report runs with at least this many genotyped sites (guards sparse-coverage false runs).
    pub min_sites: usize,
}

impl Default for RohConfig {
    fn default() -> Self {
        RohConfig {
            het_error: 0.002,
            baseline_het: None,
            switch_rate_per_cm: 1.0 / 13.0,
            prior_autozygous: 0.02,
            min_length_mb: 1.5,
            min_sites: 50,
        }
    }
}

/// A single detected run of homozygosity.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct RohSegment {
    pub chromosome: String,
    /// 1-based inclusive bp span (first/last genotyped site in the run).
    pub start_bp: i64,
    pub end_bp: i64,
    /// Genetic length from the genetic map (cM); falls back to a 1 cM/Mb estimate if the map lacks
    /// the chromosome.
    pub length_cm: f64,
    /// The physical span, in Mb. The report threshold `min_length_mb` applies to this length.
    pub length_mb: f64,
    /// Number of genotyped sites inside the run.
    pub n_sites: usize,
    /// Heterozygous calls inside the run (should be near zero for a clean run).
    pub n_het: usize,
    /// The mean Autozygous posterior over the sites of the run, from the forward-backward pass.
    /// It is a confidence in [0,1].
    pub mean_posterior: f64,
}

/// The length classes of a ROH, in physical Mb. They separate endogamy from consanguinity.
///
/// A short ROH shows distant, background relatedness, which is endogamy. A long ROH shows recent
/// shared ancestry, which is consanguinity. A longer haplotype has had fewer generations of
/// recombination to break it up.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RohClass {
    /// Below 5 Mb. This is deep, background relatedness, which is endogamy.
    Short,
    /// From 5 to 15 Mb. This is the middle class.
    Medium,
    /// 15 Mb or more. This is recent, and it is consanguinity.
    Long,
}

impl RohClass {
    pub fn of(length_mb: f64) -> Self {
        if length_mb < 5.0 {
            RohClass::Short
        } else if length_mb < 15.0 {
            RohClass::Medium
        } else {
            RohClass::Long
        }
    }
}

/// A coarse pattern that the code reads from the distribution of the ROH lengths. It is a
/// heuristic. Use it for narration, and not for a diagnosis.
///
/// It lives in `navigator-domain`, so that the Simple-mode brief can read the answer that
/// [`classify`] reaches here. That brief does not have to derive its own answer from the raw
/// numbers.
pub use navigator_domain::roh::RohPattern;

/// The rollup over the whole genome. Every length is in **physical Mb**.
///
/// The canonical F_ROH, from McQuillan, is a physical ratio. That also keeps it consistent with
/// `min_length_mb`, which filters the runs by physical length. An F_ROH in cM would count far too
/// few of the ROH that sit in a region of low recombination.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct RohSummary {
    pub n_segments: usize,
    pub total_roh_mb: f64,
    /// The autosomal physical length that the input sites cover, in Mb. It is the denominator of
    /// F_ROH.
    pub autosomal_mb: f64,
    /// The coefficient F_ROH. It is the total ROH length divided by the total autosomal length,
    /// and both are in Mb.
    pub f_roh: f64,
    pub longest_mb: f64,
    /// A (count, total Mb) pair for each length class.
    pub short: (usize, f64),
    pub medium: (usize, f64),
    pub long: (usize, f64),
    pub pattern: RohPattern,
}

/// Full result.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct RohResult {
    pub segments: Vec<RohSegment>,
    pub summary: RohSummary,
}

/// True for the autosomes 1 to 22. The code computes a ROH on an autosome alone, and it leaves X,
/// Y and MT out.
fn is_autosome(contig: &str) -> bool {
    matches!(normalize_chromosome(contig).parse::<u8>(), Ok(1..=22))
}

/// Detect runs of homozygosity across the autosomes.
pub fn detect_roh(genotypes: &[SiteGenotype], gmap: &GeneticMap, config: &RohConfig) -> RohResult {
    // The sorted (pos, is_het) pairs of each autosome, over the sites with a *call*, at dosage 0,
    // 1 or 2.
    let mut by_chr: BTreeMap<String, Vec<(i64, bool)>> = BTreeMap::new();
    let (mut called, mut het) = (0u64, 0u64);
    for g in genotypes {
        if !is_autosome(&g.contig) || !(0..=2).contains(&g.dosage) {
            continue;
        }
        called += 1;
        if g.dosage == 1 {
            het += 1;
        }
        by_chr
            .entry(normalize_chromosome(&g.contig))
            .or_default()
            .push((g.position, g.dosage == 1));
    }

    // Effective Normal-state het rate: configured, else the sample's own autosomal het fraction
    // (clamped). Genome-wide het is a slight under-estimate of the non-autozygous rate, but for a
    // mostly-outbred genome the bias is negligible; the clamp guards degenerate inputs.
    let baseline = config.baseline_het.unwrap_or_else(|| {
        if called == 0 {
            0.30
        } else {
            (het as f64 / called as f64).clamp(0.15, 0.45)
        }
    });

    let mut segments: Vec<RohSegment> = Vec::new();
    let mut autosomal_mb = 0.0f64;

    for (chr, mut sites) in by_chr {
        sites.sort_by_key(|(p, _)| *p);
        sites.dedup_by_key(|(p, _)| *p);
        if sites.len() < 2 {
            continue;
        }
        autosomal_mb += (sites.last().unwrap().0 - sites.first().unwrap().0).max(0) as f64 / 1_000_000.0;
        for run in call_chromosome(&chr, &sites, gmap, config, baseline) {
            if run.length_mb >= config.min_length_mb && run.n_sites >= config.min_sites {
                segments.push(run);
            }
        }
    }

    let summary = summarize(&segments, autosomal_mb);
    RohResult { segments, summary }
}

/// cM span between two bp positions on `chr`, with a 1 cM/Mb fallback when the map lacks the contig.
fn span_cm(gmap: &GeneticMap, chr: &str, start_bp: i64, end_bp: i64) -> f64 {
    gmap.interval_cm(chr, start_bp as i32, end_bp as i32)
        .unwrap_or_else(|| (end_bp - start_bp).max(0) as f64 / 1_000_000.0)
}

/// Log-space 2-state HMM (0 = Normal, 1 = Autozygous) over one chromosome's sorted sites; returns
/// the stitched Autozygous runs (unfiltered).
fn call_chromosome(
    chr: &str,
    sites: &[(i64, bool)],
    gmap: &GeneticMap,
    cfg: &RohConfig,
    baseline: f64,
) -> Vec<RohSegment> {
    let n = sites.len();
    let ln = |x: f64| x.max(1e-300).ln();

    // Stationary prior π and its log.
    let pi = [1.0 - cfg.prior_autozygous, cfg.prior_autozygous];
    let ln_pi = [ln(pi[0]), ln(pi[1])];

    // The emission log-likelihood at each site, for each state: [normal, auto].
    let emit = |is_het: bool| -> [f64; 2] {
        if is_het {
            [ln(baseline), ln(cfg.het_error)]
        } else {
            [ln(1.0 - baseline), ln(1.0 - cfg.het_error)]
        }
    };

    // Transition log-prob from state i to j given a cM gap: reset-to-prior with switch prob s.
    // P(j|i) = (1−s)·[i==j] + s·π_j.
    let trans = |i: usize, j: usize, gap_cm: f64| -> f64 {
        let s = 1.0 - (-gap_cm * cfg.switch_rate_per_cm).exp();
        let s = s.clamp(0.0, 1.0);
        let stay = if i == j { 1.0 - s } else { 0.0 };
        ln(stay + s * pi[j])
    };

    // ---- Viterbi (MAP path) ----
    let mut delta = [ln_pi[0] + emit(sites[0].1)[0], ln_pi[1] + emit(sites[0].1)[1]];
    let mut back: Vec<[usize; 2]> = vec![[0, 0]; n];
    for t in 1..n {
        let gap = span_cm(gmap, chr, sites[t - 1].0, sites[t].0);
        let e = emit(sites[t].1);
        let mut next = [f64::NEG_INFINITY; 2];
        for j in 0..2 {
            for (i, &d) in delta.iter().enumerate() {
                let c = d + trans(i, j, gap);
                if c > next[j] {
                    next[j] = c;
                    back[t][j] = i;
                }
            }
            next[j] += e[j];
        }
        delta = next;
    }
    let mut path = vec![0usize; n];
    path[n - 1] = if delta[1] > delta[0] { 1 } else { 0 };
    for t in (0..n - 1).rev() {
        path[t] = back[t + 1][path[t + 1]];
    }

    // ---- The forward-backward posteriors, which give the confidence of each run ----
    let posterior = forward_backward(chr, sites, gmap, cfg, baseline, &ln_pi);

    // ---- Stitch Autozygous runs ----
    let mut runs = Vec::new();
    let mut t = 0;
    while t < n {
        if path[t] != 1 {
            t += 1;
            continue;
        }
        let start = t;
        while t < n && path[t] == 1 {
            t += 1;
        }
        let end = t - 1; // inclusive
        let (s_bp, e_bp) = (sites[start].0, sites[end].0);
        let n_het = sites[start..=end].iter().filter(|(_, h)| *h).count();
        let post: f64 = posterior[start..=end].iter().sum::<f64>() / (end - start + 1) as f64;
        runs.push(RohSegment {
            chromosome: chr.to_string(),
            start_bp: s_bp,
            end_bp: e_bp,
            length_cm: span_cm(gmap, chr, s_bp, e_bp),
            length_mb: (e_bp - s_bp).max(0) as f64 / 1_000_000.0,
            n_sites: end - start + 1,
            n_het,
            mean_posterior: post,
        });
    }
    runs
}

/// The posterior of the Autozygous state at each site, from a scaled forward-backward pass.
fn forward_backward(
    chr: &str,
    sites: &[(i64, bool)],
    gmap: &GeneticMap,
    cfg: &RohConfig,
    baseline: f64,
    ln_pi: &[f64; 2],
) -> Vec<f64> {
    let n = sites.len();
    let pi = [ln_pi[0].exp(), ln_pi[1].exp()];
    let emit = |is_het: bool| -> [f64; 2] {
        if is_het {
            [baseline, cfg.het_error]
        } else {
            [1.0 - baseline, 1.0 - cfg.het_error]
        }
    };
    let trans = |i: usize, j: usize, gap_cm: f64| -> f64 {
        let s = (1.0 - (-gap_cm * cfg.switch_rate_per_cm).exp()).clamp(0.0, 1.0);
        (if i == j { 1.0 - s } else { 0.0 }) + s * pi[j]
    };

    // Forward (scaled).
    let mut alpha = vec![[0.0f64; 2]; n];
    let e0 = emit(sites[0].1);
    let mut a = [pi[0] * e0[0], pi[1] * e0[1]];
    normalize2(&mut a);
    alpha[0] = a;
    for t in 1..n {
        let gap = span_cm(gmap, chr, sites[t - 1].0, sites[t].0);
        let e = emit(sites[t].1);
        let mut nxt = [0.0f64; 2];
        for j in 0..2 {
            let mut acc = 0.0;
            for (i, &ai) in alpha[t - 1].iter().enumerate() {
                acc += ai * trans(i, j, gap);
            }
            nxt[j] = acc * e[j];
        }
        normalize2(&mut nxt);
        alpha[t] = nxt;
    }

    // Backward (scaled).
    let mut beta = vec![[0.0f64; 2]; n];
    beta[n - 1] = [1.0, 1.0];
    for t in (0..n - 1).rev() {
        let gap = span_cm(gmap, chr, sites[t].0, sites[t + 1].0);
        let e = emit(sites[t + 1].1);
        let mut b = [0.0f64; 2];
        for (i, bi) in b.iter_mut().enumerate() {
            let mut acc = 0.0;
            for j in 0..2 {
                acc += trans(i, j, gap) * e[j] * beta[t + 1][j];
            }
            *bi = acc;
        }
        normalize2(&mut b);
        beta[t] = b;
    }

    (0..n)
        .map(|t| {
            let g0 = alpha[t][0] * beta[t][0];
            let g1 = alpha[t][1] * beta[t][1];
            let z = g0 + g1;
            if z > 0.0 {
                g1 / z
            } else {
                0.0
            }
        })
        .collect()
}

fn normalize2(v: &mut [f64; 2]) {
    let z = v[0] + v[1];
    if z > 0.0 {
        v[0] /= z;
        v[1] /= z;
    } else {
        v[0] = 0.5;
        v[1] = 0.5;
    }
}

fn summarize(segments: &[RohSegment], autosomal_mb: f64) -> RohSummary {
    let mut short = (0usize, 0.0f64);
    let mut medium = (0usize, 0.0f64);
    let mut long = (0usize, 0.0f64);
    let mut total = 0.0f64;
    let mut longest = 0.0f64;
    for s in segments {
        total += s.length_mb;
        longest = longest.max(s.length_mb);
        let bucket = match RohClass::of(s.length_mb) {
            RohClass::Short => &mut short,
            RohClass::Medium => &mut medium,
            RohClass::Long => &mut long,
        };
        bucket.0 += 1;
        bucket.1 += s.length_mb;
    }
    let f_roh = if autosomal_mb > 0.0 { total / autosomal_mb } else { 0.0 };
    let pattern = classify(f_roh, total, &short, &long);
    RohSummary {
        n_segments: segments.len(),
        total_roh_mb: total,
        autosomal_mb,
        f_roh,
        longest_mb: longest,
        short,
        medium,
        long,
        pattern,
    }
}

/// The heuristic read of the pattern. It is an illustration, and nobody calibrated it.
///
/// Its key is the normalized F_ROH, so it does not depend on how much of the genome the run
/// covered. The split by length class then separates recent consanguinity, where the long ROH
/// dominate, from endogamy, where the short ones do.
fn classify(f_roh: f64, total: f64, short: &(usize, f64), long: &(usize, f64)) -> RohPattern {
    // Below ~F_ROH 0.02 (roughly a notable-relatedness floor) the sample reads as outbred.
    if f_roh < 0.02 {
        return RohPattern::Outbred;
    }
    let long_frac = if total > 0.0 { long.1 / total } else { 0.0 };
    let short_frac = if total > 0.0 { short.1 / total } else { 0.0 };
    if long_frac > 0.5 {
        RohPattern::RecentConsanguinity
    } else if short_frac > 0.5 {
        RohPattern::Endogamy
    } else {
        RohPattern::Mixed
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A uniform 1 cM/Mb map over one autosome long enough for the runs under test.
    fn map_chr1(len_bp: i32) -> GeneticMap {
        GeneticMap::uniform(1.0, &[("1", len_bp)])
    }

    fn site(pos: i64, het: bool) -> SiteGenotype {
        SiteGenotype {
            name: String::new(),
            contig: "1".into(),
            position: pos,
            reference_allele: "A".into(),
            alternate_allele: "G".into(),
            ploidy: 2,
            dosage: if het { 1 } else { 0 },
            gq: 0,
            depth: 0,
            ref_depth: 0,
            alt_depth: 0,
            pls: Vec::new(),
            gt: None,
            allele_depths: None,
        }
    }

    /// A site at every `step` bp, from 0 to count, and all of them homozygous. That gives one ROH
    /// over the whole chromosome.
    #[test]
    fn all_homozygous_is_one_long_roh() {
        let step = 20_000i64;
        let count = 1000; // 20 Mb → 20 cM
        let genos: Vec<_> = (0..count).map(|i| site(i * step, false)).collect();
        let gmap = map_chr1((count * step) as i32);
        let res = detect_roh(&genos, &gmap, &RohConfig::default());
        assert_eq!(res.segments.len(), 1, "expected a single ROH");
        let seg = &res.segments[0];
        assert!(seg.length_cm > 18.0, "run should span ~20 cM, got {}", seg.length_cm);
        assert_eq!(seg.n_het, 0);
        assert!(seg.mean_posterior > 0.9, "posterior {}", seg.mean_posterior);
        assert!(res.summary.f_roh > 0.9, "F_ROH {}", res.summary.f_roh);
        assert_eq!(res.summary.pattern, RohPattern::RecentConsanguinity);
    }

    /// Heterozygous-rich chromosome → no ROH.
    #[test]
    fn heterozygous_rich_has_no_roh() {
        let step = 20_000i64;
        let count = 1000;
        // Every third site is het. That is dense enough to hold the HMM in the Normal state all
        // the way.
        let genos: Vec<_> = (0..count).map(|i| site(i * step, i % 3 == 0)).collect();
        let gmap = map_chr1((count * step) as i32);
        let res = detect_roh(&genos, &gmap, &RohConfig::default());
        assert!(res.segments.is_empty(), "expected no ROH, got {:?}", res.segments);
        assert!(res.summary.f_roh < 0.05);
        assert_eq!(res.summary.pattern, RohPattern::Outbred);
    }

    /// First half homozygous, second half het → ROH on the first half, boundary near the midpoint.
    #[test]
    fn half_homozygous_calls_only_that_half() {
        let step = 20_000i64;
        let count = 1000;
        let mid = count / 2;
        let genos: Vec<_> = (0..count)
            .map(|i| site(i * step, if i < mid { false } else { i % 2 == 0 }))
            .collect();
        let gmap = map_chr1((count * step) as i32);
        let res = detect_roh(&genos, &gmap, &RohConfig::default());
        assert_eq!(res.segments.len(), 1, "segments: {:?}", res.segments);
        let seg = &res.segments[0];
        assert!(seg.start_bp == 0, "run should start at 0");
        // Boundary within ~1 Mb of the midpoint.
        let mid_bp = mid * step;
        assert!(
            (seg.end_bp - mid_bp).abs() < 1_000_000,
            "end {} vs midpoint {}",
            seg.end_bp,
            mid_bp
        );
    }

    /// The filter removes a short homozygous run that is below `min_length_mb`.
    #[test]
    fn short_run_below_min_length_is_dropped() {
        let step = 20_000i64;
        // 30 hom sites = 0.6 Mb, below the 1.5 Mb floor, embedded in het background.
        let count = 1000;
        let genos: Vec<_> = (0..count)
            .map(|i| site(i * step, !(400..430).contains(&i) && i % 2 == 0))
            .collect();
        let gmap = map_chr1((count * step) as i32);
        let res = detect_roh(&genos, &gmap, &RohConfig::default());
        assert!(
            res.segments.iter().all(|s| s.length_mb >= 1.5),
            "no sub-threshold run should survive: {:?}",
            res.segments
        );
    }
}
