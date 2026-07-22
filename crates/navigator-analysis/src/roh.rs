//! Runs of homozygosity (ROH) / autozygosity detection — the endogamy & consanguinity signal.
//!
//! **Spike (2026-07-22).** A 2-state hidden Markov model (Autozygous / Normal) over a subject's
//! autosomal genotypes, mirroring the [`crate::ancestry::paint_local_ancestry`] HMM idiom
//! (per-contig sorted sites, distance-scaled "reset-to-prior" transitions, log-space Viterbi +
//! forward/backward posteriors). Runs of the Autozygous state are stitched into [`RohSegment`]s and
//! rolled up into an [`RohSummary`] with the genome-wide inbreeding coefficient F_ROH.
//!
//! **Input substrate.** Feed the subject's autosomal-consensus genotypes
//! (`consensus_genotypes(&DiploidProfile)` in `navigator-app`), which are called at the full 1240k
//! IBD panel — a dense (~1.15M), neutral, biallelic common-SNP set with full 0/1/2 dosages. That is
//! the density class array-based ROH tools (PLINK, BCFtools/RoH, detectRUNS) assume. Segment cM
//! lengths and the F_ROH denominator come from the same [`GeneticMap`] the IBD path already loads.
//!
//! **What's deliberately simplified in the spike** (see the module tests + the follow-up notes):
//! - The Normal-state heterozygosity expectation is a single `baseline_het` knob. A production
//!   version should derive it per-site from panel allele frequencies (2·f·(1−f)), which
//!   `AncestryPanel`/`IbdPanel` already carry, so the emission is properly frequency-aware.
//! - The endogamy-vs-consanguinity [`RohPattern`] classification is a heuristic on the ROH
//!   length-class distribution, not a calibrated model.

use crate::caller::SiteGenotype;
use crate::ibd::{normalize_chromosome, GeneticMap};
use std::collections::BTreeMap;

/// Detector configuration. Defaults target a 1240k-density common-SNP substrate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RohConfig {
    /// P(heterozygous call | truly autozygous) — i.e. the residual genotyping-error rate that lets a
    /// stray het survive inside a homozygous run without breaking it. Small.
    pub het_error: f64,
    /// Normal-state expected heterozygosity P(het | not autozygous). `None` (the default) estimates it
    /// from the sample's own autosomal het fraction (clamped), so it tracks the panel's density and
    /// ascertainment instead of a fixed guess. A production upgrade is per-site 2·f·(1−f) from allele
    /// frequencies; `Some(v)` pins it (mainly for tests / advanced tuning).
    pub baseline_het: Option<f64>,
    /// State-switch hazard per centimorgan. Switch probability over a gap of `d` cM is
    /// `1 − exp(−d · switch_rate_per_cm)`. Smaller → longer runs. Default ≈ one switch per ~13 cM.
    pub switch_rate_per_cm: f64,
    /// Stationary autozygosity fraction — the prior mass on the Autozygous state that a switch
    /// resets toward. The classic ROH HMM prior.
    pub prior_autozygous: f64,
    /// Report runs at least this long in **physical Mb**. PLINK/detectRUNS and the genealogy field
    /// threshold ROH on physical length because ROH cluster in low-recombination (pericentromeric)
    /// regions where a multi-Mb run spans well under a cM — a genetic (cM) threshold under-reports them.
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
    /// Physical span in Mb — the length the report threshold (`min_length_mb`) applies to.
    pub length_mb: f64,
    /// Number of genotyped sites inside the run.
    pub n_sites: usize,
    /// Heterozygous calls inside the run (should be near zero for a clean run).
    pub n_het: usize,
    /// Mean Autozygous posterior over the run's sites (forward/backward) — a confidence in [0,1].
    pub mean_posterior: f64,
}

/// ROH length classes (physical Mb), used for the endogamy-vs-consanguinity read. Short ROH reflect
/// distant/background relatedness (endogamy); long ROH reflect recent shared ancestry (consanguinity),
/// because longer haplotypes have had fewer generations of recombination to break them up.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RohClass {
    /// < 5 Mb — deep/background (endogamy).
    Short,
    /// 5–15 Mb — intermediate.
    Medium,
    /// ≥ 15 Mb — recent (consanguinity).
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

/// Coarse pattern read from the ROH length distribution. Heuristic — for narration, not diagnosis.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum RohPattern {
    /// Little total ROH — outbred.
    Outbred,
    /// ROH mass dominated by short segments — background relatedness / endogamous population.
    Endogamy,
    /// ROH mass dominated by long segments — recent consanguinity in the pedigree.
    RecentConsanguinity,
    /// Substantial ROH across all classes.
    Mixed,
}

/// Genome-wide rollup. Lengths are **physical Mb** — the canonical (McQuillan) F_ROH is a physical
/// ratio, and it stays consistent with the physical `min_length_mb` run filter (a genetic/cM F_ROH
/// would badly under-count ROH that sit in low-recombination regions).
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct RohSummary {
    pub n_segments: usize,
    pub total_roh_mb: f64,
    /// Autosomal physical length covered by the input sites (Mb) — the F_ROH denominator.
    pub autosomal_mb: f64,
    /// Inbreeding coefficient F_ROH = total ROH length / total autosomal length (both Mb).
    pub f_roh: f64,
    pub longest_mb: f64,
    /// (count, summed Mb) per length class.
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

/// True for autosomes 1–22 (ROH is computed on autosomes only; X/Y/MT excluded).
fn is_autosome(contig: &str) -> bool {
    matches!(normalize_chromosome(contig).parse::<u8>(), Ok(1..=22))
}

/// Detect runs of homozygosity across the autosomes.
pub fn detect_roh(genotypes: &[SiteGenotype], gmap: &GeneticMap, config: &RohConfig) -> RohResult {
    // Per-autosome sorted (pos, is_het) over *called* sites (dosage 0/1/2).
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
fn call_chromosome(chr: &str, sites: &[(i64, bool)], gmap: &GeneticMap, cfg: &RohConfig, baseline: f64) -> Vec<RohSegment> {
    let n = sites.len();
    let ln = |x: f64| x.max(1e-300).ln();

    // Stationary prior π and its log.
    let pi = [1.0 - cfg.prior_autozygous, cfg.prior_autozygous];
    let ln_pi = [ln(pi[0]), ln(pi[1])];

    // Per-site emission log-likelihoods for each state: [normal, auto].
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

    // ---- Forward/backward posteriors (for per-run confidence) ----
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

/// Autozygous-state posterior per site via scaled forward/backward.
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

/// Heuristic pattern read (illustrative — not calibrated). Keyed on the normalized F_ROH so it is
/// independent of how much genome was analyzed; the length-class split then separates recent
/// consanguinity (long-dominated) from endogamy (short-dominated).
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

    /// Sites every `step` bp from 0..count, all homozygous → one ROH spanning the chromosome.
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
        // Every 3rd site het — dense enough to keep the HMM in the Normal state throughout.
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

    /// A short homozygous run below `min_length_mb` is filtered out.
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
