//! Tier B: the call of archaic **segments**. The design is in
//! `documents/design/ArchaicAncestry_Design.md`, §5.
//!
//! Tier A counts marker copies. This module instead finds the introgressed **tracts** themselves.
//! It is a two-state HMM in the style of hmmix (Skov and others, 2018), over the density of
//! *private* derived variants. A private variant is a variant of the subject that no individual
//! in the African outgroup carries. A variant that Africans also carry is not evidence of
//! introgression, so the removal of those is what makes the density that stays informative.
//!
//! The HMM **can not separate Neanderthal from Denisovan**, because the two lineages coalesce
//! before either one meets modern humans (§3). So it finds the segments, and a later pass puts a
//! label on each one. That pass counts the matches of the derived allele against the archaic
//! genomes, in `ArchaicClassify`.
//!
//! The model uses Viterbi and forward-backward in log space, with transitions that scale in cM.
//! That is the same idiom as [`crate::roh`] and the chromosome painter. Only the emission is
//! different. Here it is a Poisson point process over the count of private variants in each
//! window, and not het against hom.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::archaic::{ArchaicCallable, ArchaicClassify, ArchaicOutgroup, DiagnosticClass};
use crate::caller::SiteGenotype;
use crate::ibd::GeneticMap;

/// The archaic lineage that a called segment belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ArchaicSource {
    Neanderthal,
    Denisovan,
    /// The density says archaic, but the diagnostic sites inside the segment do not point to
    /// either lineage. This is the honest label for a segment that the code can not attribute, and
    /// real data holds a large share of them. Skov 2020 reported about 12 % unknown on
    /// Icelanders.
    Unknown,
}

/// One called archaic tract.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ArchaicSegment {
    pub contig: String,
    /// 1-based inclusive bounds.
    pub start: i64,
    pub end: i64,
    /// Mean posterior probability of the archaic state across the segment's windows.
    pub posterior: f64,
    /// The count of private derived variants inside the segment. That is the evidence that the
    /// call stands on.
    pub n_private: usize,
    pub source: ArchaicSource,
    /// The count of diagnostic-site matches for each lineage. `source` comes from these.
    pub neanderthal_matches: usize,
    pub denisovan_matches: usize,
}

impl ArchaicSegment {
    pub fn length_mb(&self) -> f64 {
        (self.end - self.start).max(0) as f64 / 1_000_000.0
    }
}

/// Genome-level totals.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ArchaicSummary {
    pub total_mb: f64,
    /// The archaic share of the **callable** span, as a percentage. It counts against the callable
    /// length, and not against the nominal length of the genome. A genome with partial coverage
    /// would else read as one with less archaic ancestry, for one reason alone: the run sequenced
    /// less of it.
    pub pct_callable: f64,
    pub callable_mb: f64,
    pub neanderthal_mb: f64,
    pub denisovan_mb: f64,
    pub unknown_mb: f64,
    pub n_segments: usize,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ArchaicSegmentResult {
    pub segments: Vec<ArchaicSegment>,
    pub summary: ArchaicSummary,
}

/// The controls of this module. They have the same shape as `RohConfig` and `PaintParams`.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ArchaicConfig {
    /// Window size in bp for the density process.
    pub window_bp: i64,
    /// The count of private variants that the model expects in one window, outside an archaic
    /// tract. `None` tells the code to estimate it from the genome-wide private rate of the sample
    /// itself. That estimate follows the depth, and it follows how much the outgroup removed.
    pub rate_background: Option<f64>,
    /// How many times the background rate the model expects inside an archaic tract. An
    /// introgressed haplotype carries an excess of derived variants that Africans do not have, and
    /// this is that excess.
    pub archaic_rate_multiple: f64,
    /// Prior probability of the archaic state.
    pub prior_archaic: f64,
    /// The count of state switches that the model expects in one centimorgan. It is the
    /// transition that scales with recombination.
    pub switches_per_cm: f64,
    /// Discard called tracts shorter than this.
    pub min_segment_bp: i64,
    /// Discard tracts whose mean posterior is below this.
    pub min_posterior: f64,
    /// A segment goes to a lineage only when its **enrichment** over the expected carrier rate
    /// points to that lineage by this ratio or more. Below that, the segment is `Unknown`.
    pub min_lineage_ratio: f64,
    /// Expected fraction of Neanderthal-diagnostic sites at which a non-archaic-specific genome
    /// carries the derived allele, and the same for Denisovan-diagnostic sites.
    ///
    /// These base rates are why a raw count of matches can not attribute a lineage. A measurement
    /// on the ground-truth European gave 4.3 % at Neanderthal-diagnostic sites, against 3.9 % at
    /// Denisovan-diagnostic ones. That ratio is 1.10, which separates almost nothing.
    ///
    /// A "Denisovan-diagnostic" allele that a person carries mostly shows ordinary shared
    /// ancestry, and not Denisovan introgression. So the attribution must compare the observed
    /// rate against the expected rate. It must not compare the Neanderthal count against the
    /// Denisovan count.
    pub base_rate_neanderthal: f64,
    pub base_rate_denisovan: f64,
    /// The count of diagnostic matches that the lineage in front needs before the code attributes
    /// anything. At these base rates, a few matches on a segment below one megabase are noise.
    pub min_lineage_matches: usize,
    /// True to try a lineage attribution at each segment. **The default is `false`.**
    ///
    /// The attribution code exists, and unit tests cover it. But nobody has checked it against
    /// truth. On the ground-truth European it gives the opposite of the known pattern: 0.00 Mb
    /// Neanderthal against 0.48 Mb Denisovan. §7 expects almost all Neanderthal, and no Denisovan.
    ///
    /// The base rates show the cause. A European carries the derived allele at 4.3 % of the
    /// Neanderthal-diagnostic sites, and at 3.9 % of the Denisovan-diagnostic ones. That ratio is
    /// 1.10, so at the scale of a segment there is almost no signal that separates the two.
    ///
    /// To send a lineage split out on that basis would make exactly the claim of Denisovan
    /// ancestry in Europeans that §7 forbids. So a segment goes out as archaic, with no lineage,
    /// until somebody checks the method. That is the same discipline that gated the deep
    /// ancestry.
    pub attribute_lineage: bool,
    /// The callable fraction that a window needs before the model uses it at all. The model
    /// **removes** a window below this, and it does not give that window a lower weight.
    ///
    /// Mapping error controls the variant density of a window that is mostly not callable. To keep
    /// such a window is what made the first real run call 3.62 % archaic, out of repetitive
    /// sequence.
    pub min_callable_fraction: f64,
}

impl Default for ArchaicConfig {
    fn default() -> Self {
        ArchaicConfig {
            window_bp: 1_000,
            rate_background: None,
            // The calibration ran against the hmmix 1000G callset (Zenodo, CC BY 4.0), on chr21
            // and chr22 of the ground-truth European. It gave 45 segments over 2.01 Mb, against
            // their EUR target of 43 segments over 2.09 Mb. Both lie inside the p10 to p90 spread,
            // which is 35 to 51 segments and 1.51 to 2.65 Mb.
            archaic_rate_multiple: 6.0,
            prior_archaic: 0.02,
            // 5, and not 1. The transition rate IS the prior on the tract length. At 1.0, with
            // the fallback of 1 cM/Mb, a 1 kb window switches with a probability of about 0.001.
            // That gives tracts of about 1 Mb, against a real median of 31 kb. This one parameter
            // is why the caller gave a third as many segments, and each one was much too long.
            switches_per_cm: 5.0,
            // 5 kb, and not 50 kb. The median European tract of hmmix on these chromosomes is
            // 31 kb, and its p10 is 7 kb. So a floor of 50 kb threw away more than half of all
            // the real segments, by construction.
            min_segment_bp: 5_000,
            min_posterior: 0.70,
            min_lineage_ratio: 2.0,
            base_rate_neanderthal: 0.043,
            base_rate_denisovan: 0.039,
            min_lineage_matches: 10,
            attribute_lineage: false,
            min_callable_fraction: 0.5,
        }
    }
}

/// ln(k!) for the small counts a 1 kb window produces.
fn ln_factorial(k: u32) -> f64 {
    (2..=k).map(|i| (i as f64).ln()).sum()
}

/// Poisson log-pmf.
fn ln_poisson(k: u32, lambda: f64) -> f64 {
    let l = lambda.max(1e-9);
    k as f64 * l.ln() - l - ln_factorial(k)
}

fn ln_sum_exp(a: f64, b: f64) -> f64 {
    let m = a.max(b);
    if m.is_infinite() {
        return m;
    }
    m + ((a - m).exp() + (b - m).exp()).ln()
}

/// The cM between two bp positions. It falls back to 1 cM/Mb when the map does not hold the
/// contig.
fn span_cm(gmap: &GeneticMap, chr: &str, start_bp: i64, end_bp: i64) -> f64 {
    gmap.interval_cm(chr, start_bp as i32, end_bp as i32)
        .unwrap_or_else(|| (end_bp - start_bp).max(0) as f64 / 1_000_000.0)
}

/// Call archaic tracts from a subject's genome-wide diploid calls.
///
/// `calls` must hold the de-novo diploid variant calls of one alignment. Tier B accepts WGS and
/// VCF input alone, because a chip can not give the density that this needs.
pub fn call_archaic_segments(
    calls: &[SiteGenotype],
    outgroup: &ArchaicOutgroup,
    classify: &ArchaicClassify,
    callable: &ArchaicCallable,
    gmap: &GeneticMap,
    cfg: &ArchaicConfig,
) -> ArchaicSegmentResult {
    // For each contig: the variant positions of the subject, which are the positions with a
    // non-reference allele, and the callable extent. A dosage of 0 is a reference site with a
    // call. That is real information about the coverage, but it is not a variant.
    let mut variants: BTreeMap<String, Vec<i64>> = BTreeMap::new();
    let mut extent: BTreeMap<String, (i64, i64)> = BTreeMap::new();
    let mut alleles: BTreeMap<(String, i64), (char, char, i32)> = BTreeMap::new();
    for g in calls {
        if !(0..=g.ploidy as i32).contains(&g.dosage) {
            continue;
        }
        // Key on the EXACT contig name, and never normalize it. The Tier B assets come from the
        // CHM13 VCFs and use the names of those files. The calls of the subject are on the same
        // build and use the same names.
        //
        // A normalization on one side breaks the join, and nobody sees it happen. retain_private
        // then finds no contig, and it returns nothing. The genome reads as 0% archaic, and that
        // answer looks completely reasonable.
        let contig = g.contig.clone();
        let e = extent.entry(contig.clone()).or_insert((g.position, g.position));
        e.0 = e.0.min(g.position);
        e.1 = e.1.max(g.position);
        if g.dosage >= 1 {
            variants.entry(contig.clone()).or_default().push(g.position);
            if let (Some(r), Some(a)) = (g.reference_allele.chars().next(), g.alternate_allele.chars().next()) {
                alleles.insert((contig, g.position), (r, a, g.dosage));
            }
        }
    }

    // Remove everything that the African outgroup also carries. That step turns a raw variant
    // density into a signal of introgression.
    let mut private: BTreeMap<String, Vec<i64>> = BTreeMap::new();
    let mut total_private = 0usize;
    let mut total_windows = 0f64;
    for (contig, mut positions) in variants {
        positions.sort_unstable();
        positions.dedup();
        let kept = outgroup.retain_private(&contig, &positions);
        // Count a callable window toward the background rate, and no other window. Else the
        // regions that the model never looks at make the denominator too large.
        if let Some((lo, hi)) = extent.get(&contig) {
            let mut w = *lo;
            while w <= *hi {
                if callable.callable_fraction(&contig, w) >= cfg.min_callable_fraction {
                    total_windows += 1.0;
                }
                w += cfg.window_bp;
            }
        }
        total_private += kept
            .iter()
            .filter(|&&p| callable.callable_fraction(&contig, p) >= cfg.min_callable_fraction)
            .count();
        private.insert(contig, kept);
    }

    // The background rate. Unless the caller fixes it, it is the genome-wide private density of
    // the sample itself. It follows the depth, and it follows how much the outgroup removed. Both
    // of those change from one sample to the next.
    let background = cfg
        .rate_background
        .unwrap_or_else(|| (total_private as f64 / total_windows.max(1.0)).clamp(0.001, 10.0));
    let archaic_rate = background * cfg.archaic_rate_multiple.max(1.1);

    let mut segments = Vec::new();
    for (contig, positions) in &private {
        let Some(&(lo, hi)) = extent.get(contig) else { continue };
        if positions.is_empty() || hi <= lo {
            continue;
        }
        segments.extend(call_contig(
            contig,
            positions,
            lo,
            hi,
            gmap,
            cfg,
            background,
            archaic_rate,
            classify,
            callable,
            &alleles,
        ));
    }

    // Both sides of the ratio must use the SAME units. The span of a segment holds windows that
    // the mask removed. To count that span against the callable megabases mixes the two, and it
    // makes the percentage too large. The first run read 4.80% of "96.4 Mb callable", while the
    // callable track held 44.6 Mb. So the archaic extent adds up over the callable bases alone.
    let seg_callable_mb = |s: &ArchaicSegment| -> f64 {
        let mut bp = 0.0;
        let mut w = s.start;
        while w <= s.end {
            bp += callable.callable_fraction(&s.contig, w) * cfg.window_bp as f64;
            w += cfg.window_bp;
        }
        bp / 1_000_000.0
    };
    let (mut nea, mut den, mut unk) = (0.0, 0.0, 0.0);
    for s in &segments {
        let mb = seg_callable_mb(s);
        match s.source {
            ArchaicSource::Neanderthal => nea += mb,
            ArchaicSource::Denisovan => den += mb,
            ArchaicSource::Unknown => unk += mb,
        }
    }
    let total_mb = nea + den + unk;
    // Denominator: callable bases within the analysed contigs, from the mask itself.
    let callable_mb: f64 = callable
        .contigs
        .iter()
        .filter(|c| private.contains_key(&c.contig))
        .flat_map(|c| c.callable_bp.iter())
        .map(|&b| b as f64)
        .sum::<f64>()
        / 1_000_000.0;
    let summary = ArchaicSummary {
        total_mb,
        pct_callable: if callable_mb > 0.0 {
            total_mb * 100.0 / callable_mb
        } else {
            0.0
        },
        callable_mb,
        neanderthal_mb: nea,
        denisovan_mb: den,
        unknown_mb: unk,
        n_segments: segments.len(),
    };
    ArchaicSegmentResult { segments, summary }
}

/// The HMM over one contig. It puts the callable span into windows, runs Viterbi and
/// forward-backward, joins the archaic runs, and then attributes each run to a lineage.
#[allow(clippy::too_many_arguments)]
fn call_contig(
    contig: &str,
    private: &[i64],
    lo: i64,
    hi: i64,
    gmap: &GeneticMap,
    cfg: &ArchaicConfig,
    background: f64,
    archaic_rate: f64,
    classify: &ArchaicClassify,
    callable: &ArchaicCallable,
    alleles: &BTreeMap<(String, i64), (char, char, i32)>,
) -> Vec<ArchaicSegment> {
    let n_windows = (((hi - lo) / cfg.window_bp) + 1) as usize;
    if n_windows < 2 {
        return Vec::new();
    }
    let mut counts = vec![0u32; n_windows];
    for &p in private {
        let idx = ((p - lo) / cfg.window_bp) as usize;
        if idx < n_windows {
            counts[idx] += 1;
        }
    }

    // The callable fraction of each window. A window below the floor gives no information. It
    // emits nothing in either state, so it does not support a segment, and it does not break
    // one.
    let frac: Vec<f64> = (0..n_windows)
        .map(|i| callable.callable_fraction(contig, lo + i as i64 * cfg.window_bp))
        .collect();
    let usable: Vec<bool> = frac.iter().map(|f| *f >= cfg.min_callable_fraction).collect();
    if !usable.iter().any(|u| *u) {
        return Vec::new();
    }

    let ln = |x: f64| x.max(1e-300).ln();
    let ln_pi = [ln(1.0 - cfg.prior_archaic), ln(cfg.prior_archaic)];
    // The expected counts scale with how much of the window is callable. A window that is half
    // callable then does not read as a window with few variants.
    let emit = |i: usize, k: u32| -> [f64; 2] {
        if !usable[i] {
            return [0.0, 0.0];
        }
        let f = frac[i].max(1e-3);
        [ln_poisson(k, background * f), ln_poisson(k, archaic_rate * f)]
    };

    // Transition log-probabilities between adjacent windows, scaled by genetic distance: tracts
    // break at recombination, so a wide gap in cM means a switch is likelier.
    let trans = |i: usize| -> [[f64; 2]; 2] {
        let a = lo + i as i64 * cfg.window_bp;
        let b = a + cfg.window_bp;
        let cm = span_cm(gmap, contig, a, b);
        let sw = (1.0 - (-cfg.switches_per_cm * cm.max(0.0)).exp()).clamp(1e-6, 0.5);
        [[ln(1.0 - sw), ln(sw)], [ln(sw), ln(1.0 - sw)]]
    };

    // Forward.
    let mut fwd = vec![[f64::NEG_INFINITY; 2]; n_windows];
    let e0 = emit(0, counts[0]);
    fwd[0] = [ln_pi[0] + e0[0], ln_pi[1] + e0[1]];
    for i in 1..n_windows {
        let t = trans(i - 1);
        let e = emit(i, counts[i]);
        for s in 0..2 {
            fwd[i][s] = ln_sum_exp(fwd[i - 1][0] + t[0][s], fwd[i - 1][1] + t[1][s]) + e[s];
        }
    }
    // Backward.
    let mut bwd = vec![[0.0f64; 2]; n_windows];
    for i in (0..n_windows - 1).rev() {
        let t = trans(i);
        let e = emit(i + 1, counts[i + 1]);
        for s in 0..2 {
            bwd[i][s] = ln_sum_exp(t[s][0] + e[0] + bwd[i + 1][0], t[s][1] + e[1] + bwd[i + 1][1]);
        }
    }
    // The posterior of the archaic state in each window.
    let post: Vec<f64> = (0..n_windows)
        .map(|i| {
            let a = fwd[i][0] + bwd[i][0];
            let b = fwd[i][1] + bwd[i][1];
            let m = a.max(b);
            let (ea, eb) = ((a - m).exp(), (b - m).exp());
            eb / (ea + eb)
        })
        .collect();

    // Stitch runs of archaic-posterior windows, then apply the length and confidence floors.
    let mut out = Vec::new();
    let mut i = 0usize;
    while i < n_windows {
        if post[i] < 0.5 {
            i += 1;
            continue;
        }
        let start_w = i;
        while i < n_windows && post[i] >= 0.5 {
            i += 1;
        }
        let end_w = i - 1;
        let start = lo + start_w as i64 * cfg.window_bp;
        let end = (lo + (end_w as i64 + 1) * cfg.window_bp - 1).min(hi);
        if end - start < cfg.min_segment_bp {
            continue;
        }
        let mean_post = post[start_w..=end_w].iter().sum::<f64>() / (end_w - start_w + 1) as f64;
        if mean_post < cfg.min_posterior {
            continue;
        }
        // Half of the windows of a run must be callable, or more. A run of windows that give no
        // information, and that stand on the transition prior alone, can then not hold up a
        // tract.
        let callable_windows = usable[start_w..=end_w].iter().filter(|u| **u).count();
        if callable_windows * 2 < end_w - start_w + 1 {
            continue;
        }
        let n_private = private.iter().filter(|&&p| p >= start && p <= end).count();
        let (source, nea, den) = attribute(contig, start, end, classify, alleles, cfg);
        out.push(ArchaicSegment {
            contig: contig.to_string(),
            start,
            end,
            posterior: mean_post,
            n_private,
            source,
            neanderthal_matches: nea,
            denisovan_matches: den,
        });
    }
    out
}

/// Attribute a segment by counting the subject's derived-allele matches at diagnostic sites.
///
/// It needs a clear margin, which is `min_lineage_ratio`. A segment whose evidence is equal on the
/// two sides is `Unknown`, and so is a segment with no diagnostic site at all. To guess here would
/// make exactly the claim of Denisovan ancestry in Europeans that §7 forbids.
fn attribute(
    contig: &str,
    start: i64,
    end: i64,
    classify: &ArchaicClassify,
    alleles: &BTreeMap<(String, i64), (char, char, i32)>,
    cfg: &ArchaicConfig,
) -> (ArchaicSource, usize, usize) {
    // Count two things for each lineage: the MATCHES, and the diagnostic sites that are there.
    // The enrichment is matches / (sites x base rate). Without the site counts, the comparison
    // points to whichever lineage happens to have more sites in this segment, and nobody sees
    // that happen.
    let (mut nea, mut den) = (0usize, 0usize);
    let (mut nea_sites, mut den_sites) = (0usize, 0usize);
    for (pos, derived, class) in classify.in_range(contig, start, end) {
        match class {
            DiagnosticClass::Neanderthal => nea_sites += 1,
            DiagnosticClass::Denisovan => den_sites += 1,
            DiagnosticClass::SharedArchaic => {}
        }
        let Some(&(r, a, dosage)) = alleles.get(&(contig.to_string(), pos)) else {
            continue;
        };
        let d = derived.to_ascii_uppercase();
        // Does the subject carry the archaic-derived base here?
        let carries = (a.to_ascii_uppercase() == d && dosage >= 1) || (r.to_ascii_uppercase() == d && dosage <= 1);
        if !carries {
            continue;
        }
        match class {
            DiagnosticClass::Neanderthal => nea += 1,
            DiagnosticClass::Denisovan => den += 1,
            DiagnosticClass::SharedArchaic => {}
        }
    }
    // Enrichment over the expected carrier rate, not raw counts.
    let exp_nea = nea_sites as f64 * cfg.base_rate_neanderthal;
    let exp_den = den_sites as f64 * cfg.base_rate_denisovan;
    let enr_nea = if exp_nea > 0.0 { nea as f64 / exp_nea } else { 0.0 };
    let enr_den = if exp_den > 0.0 { den as f64 / exp_den } else { 0.0 };

    if !cfg.attribute_lineage {
        // Machinery preserved and unit-tested; the label is withheld until it validates.
        return (ArchaicSource::Unknown, nea, den);
    }
    let source = if nea >= cfg.min_lineage_matches && enr_nea >= enr_den * cfg.min_lineage_ratio {
        ArchaicSource::Neanderthal
    } else if den >= cfg.min_lineage_matches && enr_den >= enr_nea * cfg.min_lineage_ratio {
        ArchaicSource::Denisovan
    } else {
        // There is not enough evidence to separate the lineages. For a European this is the
        // expected result, and the honest one. §7 does not let the code invent Denisovan
        // ancestry.
        ArchaicSource::Unknown
    };
    (source, nea, den)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::archaic::{ArchaicCallable, CallableContig, ClassifyContig, PositionStream};

    fn gt(contig: &str, position: i64, dosage: i32) -> SiteGenotype {
        SiteGenotype {
            name: String::new(),
            contig: contig.into(),
            position,
            reference_allele: "A".into(),
            alternate_allele: "G".into(),
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

    /// Background variants at every 5 kb across 2 Mb, plus a dense block at every 200 bp in
    /// [1.0 Mb, 1.3 Mb], which is an introgressed tract. The outgroup holds none of them. The
    /// callable track covers the whole synthetic 2 Mb contig.
    fn callable_all() -> ArchaicCallable {
        ArchaicCallable {
            build: "chm13v2.0".into(),
            window_bp: 1_000,
            contigs: vec![CallableContig {
                contig: "chr21".into(),
                start: 0,
                callable_bp: vec![1_000u16; 2_100],
            }],
        }
    }

    fn sample() -> (Vec<SiteGenotype>, ArchaicOutgroup) {
        let mut calls = Vec::new();
        let mut p = 1i64;
        while p < 2_000_000 {
            calls.push(gt("chr21", p, 1));
            p += 5_000;
        }
        let mut q = 1_000_000i64;
        while q < 1_300_000 {
            calls.push(gt("chr21", q, 1));
            q += 200;
        }
        let og = ArchaicOutgroup {
            build: "chm13v2.0".into(),
            min_allele_count: 1,
            contigs: vec![PositionStream::encode("chr21", &[])],
        };
        (calls, og)
    }

    #[test]
    fn finds_a_dense_private_block_and_ignores_background() {
        let (calls, og) = sample();
        let classify = ArchaicClassify {
            build: "chm13v2.0".into(),
            contigs: Vec::new(),
        };
        let r = call_archaic_segments(
            &calls,
            &og,
            &classify,
            &callable_all(),
            &GeneticMap::from_markers(Vec::new()),
            &ArchaicConfig::default(),
        );
        assert_eq!(r.segments.len(), 1, "exactly the dense block should call");
        let s = &r.segments[0];
        assert!(s.start >= 950_000 && s.start <= 1_050_000, "start {} off", s.start);
        assert!(s.end >= 1_250_000 && s.end <= 1_350_000, "end {} off", s.end);
        assert!(s.posterior > 0.9);
        // No diagnostic sites supplied -> honestly Unknown, not guessed.
        assert_eq!(s.source, ArchaicSource::Unknown);
        assert!(r.summary.total_mb > 0.2 && r.summary.total_mb < 0.4);
    }

    #[test]
    fn stripping_the_outgroup_is_what_makes_the_signal() {
        // The same data, but now the outgroup carries every variant in the dense block. The
        // excess of density then goes away, and the code must call nothing. This is the step that
        // separates introgression from ordinary variation.
        let (calls, _) = sample();
        let mut shared: Vec<i64> = (1_000_000..1_300_000).step_by(200).collect();
        shared.sort_unstable();
        let og = ArchaicOutgroup {
            build: "chm13v2.0".into(),
            min_allele_count: 1,
            contigs: vec![PositionStream::encode("chr21", &shared)],
        };
        let classify = ArchaicClassify {
            build: "chm13v2.0".into(),
            contigs: Vec::new(),
        };
        let r = call_archaic_segments(
            &calls,
            &og,
            &classify,
            &callable_all(),
            &GeneticMap::from_markers(Vec::new()),
            &ArchaicConfig::default(),
        );
        assert!(r.segments.is_empty(), "outgroup-shared density must not call archaic");
    }

    #[test]
    fn an_uncallable_dense_region_is_not_called() {
        // The same dense block, but the mask says that the region is not callable. That is the
        // exact case that made the first real run report 3.62% archaic, out of repetitive
        // sequence.
        let (calls, og) = sample();
        let mut track = callable_all();
        for w in 950..1_350 {
            track.contigs[0].callable_bp[w] = 0;
        }
        let classify = ArchaicClassify {
            build: "chm13v2.0".into(),
            contigs: Vec::new(),
        };
        let r = call_archaic_segments(
            &calls,
            &og,
            &classify,
            &track,
            &GeneticMap::from_markers(Vec::new()),
            &ArchaicConfig::default(),
        );
        assert!(
            r.segments.is_empty(),
            "density in an uncallable region must not be called archaic, got {:?}",
            r.segments
        );
    }

    #[test]
    fn attribution_compares_enrichment_not_raw_counts() {
        // The segment holds 40 Neanderthal-diagnostic sites and 40 Denisovan ones. The subject
        // matches 12 of each, so the raw counts give a tie. Against base rates of 4.3% and 3.9%,
        // the expected counts are 1.7 and 1.6, so the enrichment of the two is almost the same.
        // The result is a true Unknown, which is the honest answer, and not a coin toss.
        let mut positions: Vec<i64> = Vec::new();
        let mut derived: Vec<u8> = Vec::new();
        let mut classes: Vec<u8> = Vec::new();
        for i in 0..80i64 {
            positions.push(1_000_000 + i * 10);
            derived.push(b'G');
            classes.push(if i < 40 { 0 } else { 1 });
        }
        let cls = ArchaicClassify {
            build: "chm13v2.0".into(),
            contigs: vec![ClassifyContig {
                positions: PositionStream::encode("chr21", &positions),
                derived,
                classes,
            }],
        };
        let mut alleles = BTreeMap::new();
        for i in 0..12i64 {
            alleles.insert(("chr21".to_string(), 1_000_000 + i * 10), ('A', 'G', 1));
        }
        for i in 40..52i64 {
            alleles.insert(("chr21".to_string(), 1_000_000 + i * 10), ('A', 'G', 1));
        }
        // Attribution is off by default (see `attribute_lineage`); this test exercises the logic.
        let cfg = ArchaicConfig {
            attribute_lineage: true,
            ..ArchaicConfig::default()
        };
        let (src, nea, den) = attribute("chr21", 1_000_000, 1_001_000, &cls, &alleles, &cfg);
        assert_eq!((nea, den), (12, 12));
        assert_eq!(src, ArchaicSource::Unknown, "equal enrichment must not pick a lineage");

        // Now give Neanderthal a genuine excess and strip the Denisovan matches: it should call.
        let mut alleles2 = BTreeMap::new();
        for i in 0..30i64 {
            alleles2.insert(("chr21".to_string(), 1_000_000 + i * 10), ('A', 'G', 1));
        }
        let (src, nea, den) = attribute("chr21", 1_000_000, 1_001_000, &cls, &alleles2, &cfg);
        assert_eq!((nea, den), (30, 0));
        assert_eq!(src, ArchaicSource::Neanderthal);

        // Too few matches to mean anything at these base rates -> Unknown, not a guess.
        let mut alleles3 = BTreeMap::new();
        for i in 0..3i64 {
            alleles3.insert(("chr21".to_string(), 1_000_000 + i * 10), ('A', 'G', 1));
        }
        let (src, _, _) = attribute("chr21", 1_000_000, 1_001_000, &cls, &alleles3, &cfg);
        assert_eq!(src, ArchaicSource::Unknown);
    }
}
