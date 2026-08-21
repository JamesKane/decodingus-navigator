//! Tier B, the second try. It finds archaic tracts by a **match against the archaic genomes**,
//! and not by a count of mutations.
//!
//! # Why this replaces the density caller
//!
//! [`crate::archaic_segments`] follows Skov 2018 (hmmix). It removes the variants that Africans
//! also carry, and then looks for a region that is dense in what stays. That method exists for a
//! person who does **not** have archaic reference genomes, and who must infer them indirectly.
//! This project has all four, and it already ships [`ArchaicClassify`], which holds 2,031,406
//! sites where the archaics carry a derived allele.
//!
//! A measurement on a real European, against the hmmix calls for the same person, shows a large
//! difference. Both observables carry the same contrast of about 3x. But they differ 30-fold in
//! how much evidence one tract holds, and that is what decides whether the code can call a tract
//! at all:
//!
//! | observable | evidence in a 36 kb tract | sensitivity at 5 % false positives |
//! |---|---|---|
//! | private-variant density | ~1 variant | 14.3 % |
//! | archaic-allele matching (this) | ~30 sites | 95.1 % |
//!
//! Density does not get to 80 % sensitivity at **500 kb**. Matching gets to 95 % at the real
//! median tract of 36 kb. See `documents/design/ArchaicAncestry_Design.md`, § *Why it failed*.
//!
//! # The model
//!
//! An introgressed tract is a haplotype that came down whole from an archaic ancestor. So it
//! carries the archaic allele at a large share of the diagnostic sites that it covers. Elsewhere
//! the subject carries those alleles at the background rate alone.
//!
//! That is a two-state HMM. Its observation is one **bit at each diagnostic site**: the subject
//! carries the allele, or does not. The emissions are Bernoulli. The index runs **over the sites
//! and not over the base pairs**.
//!
//! The index over sites is what makes this robust where the density model was not. The diagnostic
//! sites become the denominator, so their uneven density cancels. The density model needed a map
//! of the mutation rate, and no available proxy gave one. The best proxy explained 38 % of a 14.6x
//! overdispersion. This model does not need that map at all.
//!
//! The transitions scale with recombination between one site and the next, as they do in
//! [`crate::roh`] and in the chromosome painter.
//!
//! # The checks
//!
//! ## Across the genome, which is the configuration that ships
//!
//! Three Europeans, called across all 22 autosomes, and scored against the genome-wide hmmix
//! callset for the same individuals:
//!
//! | | ours | hmmix | ratio | sensitivity | precision | null (max of 400 draws) |
//! |---|---|---|---|---|---|---|
//! | HG00096 | 83.6 Mb | 93.0 | 0.90 | 40.3 % | 44.9 % | 5.5 % |
//! | HG00102 | 83.9 Mb | 89.3 | 0.94 | 42.4 % | 45.1 % | 4.9 % |
//! | HG00112 | 82.1 Mb | 91.0 | 0.90 | 42.9 % | 47.5 % | 5.1 % |
//!
//! All three sit above the *entire* null from random placement. Both the sensitivity and the
//! precision are **better** across the genome than they are on chr21 and chr22: 40 to 43 % against
//! 31.6 %, and about 46 % against 34.9 %.
//!
//! The two-chromosome figures below are careful, and not optimistic. Say that here,
//! because the opposite burned the design of the caller before this one. That design took a
//! chr21+22 target and went outside the measured range, to a value 6 % too low.
//!
//! ## chr21 and chr22, with a split into train and test
//!
//! 60 Europeans on chr21 and chr22, scored against the hmmix calls for the same individuals. A
//! fixed seed split them into 30 for **train** and 30 for **test**. The fit of the thresholds used
//! the train half alone, and every figure below comes from the half that the fit did not see.
//!
//! The split exists because of what happened before. Somebody tuned the caller before this one
//! until a cohort statistic agreed, and then reported that statistic as evidence.
//!
//! | | density caller | this, uncalibrated | this, calibrated |
//! |---|---|---|---|
//! | base-level F1 | n/a | 27.9 % | **34.5 %** |
//! | precision | 1.5 % | 20.2 % | **34.9 %** |
//! | extent ratio ours/theirs | 1.45 | 2.23 | **0.98** |
//! | extent `r` over the individuals | −0.018 (p = 0.94) | +0.520 | **+0.710 (p < 0.0001)** |
//!
//! Look at the extent ratio of 0.98. The caller no longer calls too much everywhere, and the sweep
//! over the emission ratio is what fixed that.
//!
//! Now the locations. All 20 individuals of an earlier cohort scored above their own null from
//! random placement, at a mean sensitivity of 45.3 % against a null of 7.1 %. The density caller
//! scored 2.1 % against a null of 5.0 %, which is below chance.
//!
//! ## Across populations: the detection transfers, and the reported number does not
//!
//! A run on 30 East Asians, with the parameters **frozen** at the European fit. The run fitted
//! nothing again:
//!
//! | | Europe (fitted) | East Asia (new) |
//! |---|---|---|
//! | above own random-placement null | 60/60 | **30/30** |
//! | sensitivity | 31.6 % | **31.6 %** |
//! | precision | 32.2 % | **41.9 %** |
//! | extent `r` over the individuals | +0.620 | **+0.545** |
//!
//! The detection transfers. The sensitivity is the same, and the precision is *better*, on a
//! population that the thresholds never saw. The calibration learned archaic structure, and not
//! European structure.
//!
//! **But the reported extent puts the two populations in the wrong order.** The truth puts the
//! archaic extent of East Asia at **1.217x** that of Europe. The extent that this caller reports
//! is **0.937x**. A user would read that an East Asian carries *less* archaic ancestry than a
//! European. That is the wrong way round, and it is the one reason that this module stays gated.
//!
//! Here is the cause. The reported extent is the true positives *plus* the false positives. The
//! load of false positives depends on the population, at a precision of 32.2 % against 41.9 %.
//! Europeans then collect more extent that is not real.
//!
//! One thing is **not** evidence against this: that the detected sequence reproduces the 1.22x
//! ratio. The detected extent is the sensitivity times the truth, and the sensitivity is equal
//! across the two populations. So that ratio agrees by construction. It repeats the invariance,
//! and it does not test the order.
//!
//! Three causes are out, and a measurement rules out each one. No argument does.
//!
//! Contamination of `p_background`: the rates at which the two populations carry the allele are
//! 11.9 % against 12.2 %, and both states scale together. Tract length: the median is 29 kb in
//! both, and East Asians only have more tracts, at 54 against 46 for each person. Panel
//! ascertainment: the contrast inside a tract is 2.99x against 3.04x, a ratio of 1.014, so the
//! panel carries equal information in both.
//!
//! ## How much of the "false positive" rate belongs to this caller
//!
//! The precision counts against hmmix, but a call that hmmix did not make is not wrong by that
//! fact alone. There is an independent arbiter, and it needs no second caller. The Tier A panel
//! records, at each site, which of the four archaic genomes carries the derived allele. **This
//! caller never sees that**: it reads a derived base and a lineage class alone. The concordance
//! with each genome is evidence that nobody could have fitted it to.
//!
//! Of the sites where a given archaic genome is derived, this is the fraction that the subject
//! carries, for the genome that matches best:
//!
//! | | true positive | false positive | background |
//! |---|---|---|---|
//! | Europe | 93.6 % | **81.3 %** | 59.0 % |
//! | East Asia | 93.5 % | **72.9 %** | 45.5 % |
//!
//! The "false positives" of this caller sit **64 %** and **57 %** of the way from the background
//! to a true positive. They are a mixture: real tracts that hmmix missed, plus true noise, plus
//! calls that sit in the correct place and reach too far. So the precision against hmmix
//! **understates** this caller. That is not enough to dismiss the figure, and F1 stays a usable
//! objective.
//!
//! Note the background rates in the two populations: 59.0 % against 45.5 %. Europeans carry an
//! archaic-derived allele more often *outside* a tract. That is a candidate mechanism for the
//! load of false positives that changes with the population, and so for the inverted order
//! above.
//!
//! ## A concordance filter fixes the precision, and shows a harder limit
//!
//! Score each called segment against the archaic genomes, and drop the poor matches. That raises
//! the **precision from 54 % to 90 %**. The filter is sound. With Denisova held out of it
//! completely, the segments that stay score 74.9 % on Denisova concordance, against 21.5 % for the
//! segments that go. That is a separation of 3.5x on a genome that the filter never saw.
//!
//! It does **not** fix the order of the populations, and a tighter filter makes that order worse.
//! At 90 % precision the reported extent is mostly true positives. It still puts the populations
//! in the wrong order, so false positives are no longer the cause. What remains is a difference in
//! *recovery*: about 46 % of the European truth against 38 % of the East Asian truth.
//!
//! The concordance itself shows the reason. East Asian tracts match our archaic genomes less well
//! than European tracts do, at 83.4 % against 89.2 %. And Denisova is the best match for **32.2 %
//! of the East Asian tracts, against 11.2 % of the European** ones.
//!
//! That 2.9x is the known Denisovan ancestry that East Asians carry and Europeans almost do not,
//! and the data reproduces it. But it also means that our four sequenced archaic genomes
//! **under-represent the archaic diversity of East Asia**. Any filter that uses those references
//! under-calls East Asians. To hold Denisova out, which was the first design here, makes it much
//! worse.
//!
//! That is a limit of the approach, and not a threshold to tune. To fix it would need archaic
//! genomes nearer to the populations that introgressed into East Asia, and those do not exist. **A
//! number that you can compare across populations is not possible this way at present.** The
//! caller is defensible inside one population, and not between two.
//!
//! **This is still not enough to turn the module on.** Beyond the order of the populations, there
//! are three more reasons.
//!
//! The precision is 34.9 % without the filter, on held-out Europeans. The cohort is **chr21 and
//! chr22 alone**. And the reference callset itself has weak support: the tracts of hmmix show an
//! enrichment of only 1.84x for their own archaic SNPs. Agreement with that callset then stops
//! well below 100 %, even for a caller that is correct. F1 alone can not tell you when this work
//! reaches its end.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::archaic::{ArchaicCallable, ArchaicClassify, ArchaicMarkerPanel, DiagnosticClass, ARCHAIC_GENOMES};
use crate::archaic_segments::{ArchaicSegment, ArchaicSegmentResult, ArchaicSource, ArchaicSummary};
use crate::caller::SiteGenotype;
use crate::ibd::GeneticMap;

/// Raise this after a change that would alter the segments that this module makes.
///
/// A stored result carries it as part of its key, in `archaic_segment_sig`. A workspace that holds
/// output from an earlier method derives that output again. The private-variant density
/// caller, which this project withdrew, is the case that matters. Without the key, such a
/// workspace would serve answers that the current code would never make.
pub const METHOD_VERSION: u32 = 1;

/// The concordance that a segment must reach before the code keeps it. This is a measured value,
/// from the sweep over the threshold where the precision stops to rise. It goes from 54 % to 90 %
/// at 0.70. Above 0.70 the precision gets no better, and the recall continues to fall.
pub const MIN_CONCORDANCE: f64 = 0.70;

/// The count of sites that a genome needs in a segment before the code trusts its concordance.
/// Without a floor, a genome with a call at one site alone scores 1.0 and wins every segment.
pub const MIN_CONCORDANCE_SITES: usize = 3;

/// One diagnostic site, reduced to what the HMM consumes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SiteObs {
    pub position: i64,
    /// Whether the subject carries the archaic-derived allele here.
    pub carries: bool,
    pub class: DiagnosticClass,
}

/// The controls of this module.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MatchConfig {
    /// Rate at which a **non-introgressed** genome carries the archaic allele at a diagnostic site.
    /// `None` estimates it from the subject's own genome-wide rate, which absorbs coverage, call
    /// behaviour and ancestry.
    ///
    /// The code estimates this rate, and does not hold it fixed, because it is the denominator of
    /// the whole inference. But it estimates it *directly*, and not by EM.
    ///
    /// Baum-Welch with no constraint, on the caller before this one, went to a degenerate fit. It
    /// gave an emission ratio of 22x, tracts of 9 kb, and 7x the true extent. So a measurement
    /// gives every parameter here, and nobody fits one.
    pub p_background: Option<f64>,
    /// Rate inside an introgressed tract. `None` derives it as `p_background * archaic_ratio`.
    pub p_archaic: Option<f64>,
    /// How many times the background rate the code expects inside a tract, when `p_archaic` is
    /// `None`. A measurement gave 3.04x: 39.5 % inside a real tract, against 13.0 % elsewhere.
    pub archaic_ratio: f64,
    /// The count of state switches that the model expects in one centimorgan.
    pub switches_per_cm: f64,
    /// Discard tracts whose mean posterior is below this.
    pub min_posterior: f64,
    /// The smallest count of diagnostic sites in a tract. A tract that stands on one site or two
    /// is exactly how the density caller failed, in a new observable.
    pub min_sites: usize,
    /// Discard tracts shorter than this.
    pub min_segment_bp: i64,
    /// The callable fraction that the window of a site needs before the code uses that site at
    /// all.
    pub min_callable_fraction: f64,
    /// True to try a Neanderthal or Denisovan attribution at each segment. The default is
    /// `false`, as it was for the density caller. Nobody has shown that the lineage signal works,
    /// and this module alone does not change that.
    pub attribute_lineage: bool,
}

impl Default for MatchConfig {
    fn default() -> Self {
        MatchConfig {
            p_background: None,
            p_archaic: None,
            // A fit gives this value, and no measurement does. The observed enrichment inside a
            // real tract is 3.04x, but the model separates best at 4.5x. The two do not
            // disagree.
            // 3.04x is the *mean* over a tract set from outside this project, and that set itself
            // has weak support. The emission ratio is instead what makes the HMM selective enough
            // to place a boundary. The fit used 30 Europeans, and the report uses 30 held-out
            // ones. This is the parameter that stopped the over-call, and it moved the extent
            // ratio from 2.23 to 0.98.
            archaic_ratio: 4.5,
            switches_per_cm: 1.0,
            // The calibration of all three used the train half, and the report uses the held-out
            // test half. See the module documentation. The objective was the F1 at the base
            // level. Sensitivity alone comes from a call over more sequence: the caller before
            // the calibration called 2.2x too much extent, and it still scored 45 %.
            min_posterior: 0.98,
            min_sites: 16,
            // 5 kb, although the argmax over the grid took 10 kb. Inside the flat part of the
            // curve the two differ by 0.1 F1 points. The 5 kb floor is a little BETTER on the
            // correlation of the extent of each individual, at +0.710 against +0.706. And it
            // throws away half as many real tracts: 8 % of the truth lies below 5 kb, against
            // 16 % below 10 kb.
            //
            // An earlier sweep wanted 40 kb, which would have thrown away 61 %. The design
            // records the same trap once before, at 50 kb. To exclude real tracts by
            // construction is not worth a tenth of a point.
            min_segment_bp: 5_000,
            min_callable_fraction: 0.5,
            attribute_lineage: false,
        }
    }
}

/// Reduce one contig's diagnostic sites to observations.
///
/// `ref_base` gives the reference base at a position. The code drops a site where the
/// archaic-derived allele **is** the reference base.
///
/// There are two reasons. At such a site every genome that matches the reference holds the derived
/// allele, which separates nothing and would weaken the contrast. And the caller emits a variant
/// record alone, so a no-call there means that the subject *does* carry the allele. That is the
/// opposite of what a no-call means everywhere else.
///
/// A site with no variant record is hom-reference, so the subject does **not** carry the allele
/// there. Do not instead keep the sites where the subject happens to have a call. That is the trap
/// that made an early version of this analysis report a rate of 80 %, against a known background
/// of 4.3 %. It takes only the sites where a variant already exists.
pub fn observations_for_contig(
    contig: &str,
    classify: &ArchaicClassify,
    calls_by_pos: &BTreeMap<i64, &SiteGenotype>,
    ref_base: impl Fn(i64) -> Option<u8>,
    callable: &ArchaicCallable,
    min_callable_fraction: f64,
) -> Vec<SiteObs> {
    let Some(c) = classify.contigs.iter().find(|c| c.positions.contig == contig) else {
        return Vec::new();
    };
    let mut out = Vec::new();
    for (i, pos) in c.positions.iter().enumerate() {
        let Some(&derived) = c.derived.get(i) else { continue };
        // Uninformative: the reference already carries the archaic allele.
        if ref_base(pos) == Some(derived) {
            continue;
        }
        if callable.callable_fraction(contig, pos) < min_callable_fraction {
            continue;
        }
        let carries = calls_by_pos
            .get(&pos)
            .is_some_and(|g| g.dosage > 0 && g.alternate_allele.as_bytes().first() == Some(&derived));
        let class = match c.classes.get(i).copied().unwrap_or(2) {
            0 => DiagnosticClass::Neanderthal,
            1 => DiagnosticClass::Denisovan,
            _ => DiagnosticClass::SharedArchaic,
        };
        out.push(SiteObs {
            position: pos,
            carries,
            class,
        });
    }
    out
}

fn ln(x: f64) -> f64 {
    x.max(1e-300).ln()
}

fn ln_sum_exp(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let m = a.max(b);
    m + ((a - m).exp() + (b - m).exp()).ln()
}

/// Posterior probability of the archaic state at each observation.
///
/// This is forward-backward in log space, with the transitions scaled by recombination, as in
/// [`crate::roh`]. It is public so that a test can check the decoding against a posterior computed
/// by hand, and build no asset.
pub fn posteriors(
    obs: &[SiteObs],
    contig: &str,
    gmap: &GeneticMap,
    p_bg: f64,
    p_arch: f64,
    switches_per_cm: f64,
) -> Vec<f64> {
    let n = obs.len();
    if n == 0 {
        return Vec::new();
    }
    let emit = |i: usize| -> [f64; 2] {
        if obs[i].carries {
            [ln(p_bg), ln(p_arch)]
        } else {
            [ln(1.0 - p_bg), ln(1.0 - p_arch)]
        }
    };
    // Switch probability between consecutive sites, from the genetic distance between them.
    let sw = |i: usize| -> f64 {
        let cm = gmap
            .interval_cm(contig, obs[i].position as i32, obs[i + 1].position as i32)
            .unwrap_or_else(|| (obs[i + 1].position - obs[i].position).max(0) as f64 / 1_000_000.0);
        (1.0 - (-switches_per_cm * cm.max(0.0)).exp()).clamp(1e-9, 0.5)
    };

    // The prior is the stationary share of the archaic state. It comes from the rates themselves,
    // and it is not a tuned constant. With p_arch > p_bg, the algebra puts it at a few percent,
    // which agrees with reality.
    let prior_arch = ((p_bg - (1.0 - p_arch) * 0.0) / p_arch).clamp(0.001, 0.5) * 0.1;
    let mut fwd = vec![[f64::NEG_INFINITY; 2]; n];
    let e0 = emit(0);
    fwd[0] = [ln(1.0 - prior_arch) + e0[0], ln(prior_arch) + e0[1]];
    for i in 1..n {
        let s = sw(i - 1);
        let (stay, go) = (ln(1.0 - s), ln(s));
        let e = emit(i);
        for st in 0..2 {
            let from0 = fwd[i - 1][0] + if st == 0 { stay } else { go };
            let from1 = fwd[i - 1][1] + if st == 1 { stay } else { go };
            fwd[i][st] = ln_sum_exp(from0, from1) + e[st];
        }
    }
    let mut bwd = vec![[0.0f64; 2]; n];
    for i in (0..n - 1).rev() {
        let s = sw(i);
        let (stay, go) = (ln(1.0 - s), ln(s));
        let e = emit(i + 1);
        for st in 0..2 {
            let to0 = bwd[i + 1][0] + e[0] + if st == 0 { stay } else { go };
            let to1 = bwd[i + 1][1] + e[1] + if st == 1 { stay } else { go };
            bwd[i][st] = ln_sum_exp(to0, to1);
        }
    }
    let total = ln_sum_exp(fwd[n - 1][0], fwd[n - 1][1]);
    (0..n)
        .map(|i| (fwd[i][1] + bwd[i][1] - total).exp().clamp(0.0, 1.0))
        .collect()
}

/// Call archaic tracts for one subject by matching the archaic genomes.
///
/// `observations` holds one entry for each contig, and [`observations_for_contig`] has already
/// reduced them. This function does no I/O, and it decodes no asset. It is the model alone, and a
/// unit test can cover it as such.
pub fn call_from_observations(
    observations: &BTreeMap<String, Vec<SiteObs>>,
    gmap: &GeneticMap,
    callable: &ArchaicCallable,
    cfg: &MatchConfig,
) -> ArchaicSegmentResult {
    let (carried, total): (usize, usize) = observations
        .values()
        .flatten()
        .fold((0, 0), |(c, t), o| (c + usize::from(o.carries), t + 1));
    if total == 0 {
        return ArchaicSegmentResult {
            segments: Vec::new(),
            summary: ArchaicSummary {
                total_mb: 0.0,
                pct_callable: 0.0,
                callable_mb: 0.0,
                neanderthal_mb: 0.0,
                denisovan_mb: 0.0,
                unknown_mb: 0.0,
                n_segments: 0,
            },
        };
    }
    // Sequence that is not archaic controls the genome-wide rate, because the archaic tracts are
    // a few percent of the genome. That rate gives the background directly.
    let p_bg = cfg
        .p_background
        .unwrap_or((carried as f64 / total as f64).clamp(0.001, 0.5));
    let p_arch = cfg
        .p_archaic
        .unwrap_or((p_bg * cfg.archaic_ratio).clamp(p_bg * 1.1, 0.95));

    let mut segments = Vec::new();
    for (contig, obs) in observations {
        if obs.len() < cfg.min_sites {
            continue;
        }
        let post = posteriors(obs, contig, gmap, p_bg, p_arch, cfg.switches_per_cm);
        let mut i = 0usize;
        while i < post.len() {
            if post[i] < cfg.min_posterior {
                i += 1;
                continue;
            }
            let start = i;
            while i < post.len() && post[i] >= cfg.min_posterior {
                i += 1;
            }
            let end = i - 1;
            let n_sites = end - start + 1;
            let span = obs[end].position - obs[start].position;
            if n_sites < cfg.min_sites || span < cfg.min_segment_bp {
                continue;
            }
            let mean_post = post[start..=end].iter().sum::<f64>() / n_sites as f64;
            let (mut nea, mut den) = (0usize, 0usize);
            for o in &obs[start..=end] {
                if !o.carries {
                    continue;
                }
                match o.class {
                    DiagnosticClass::Neanderthal => nea += 1,
                    DiagnosticClass::Denisovan => den += 1,
                    DiagnosticClass::SharedArchaic => {}
                }
            }
            segments.push(ArchaicSegment {
                contig: contig.clone(),
                start: obs[start].position,
                end: obs[end].position,
                posterior: mean_post,
                n_private: obs[start..=end].iter().filter(|o| o.carries).count(),
                // Attribution stays off by default; the lineage signal is a separate question this
                // module does not answer (see `MatchConfig::attribute_lineage`).
                source: ArchaicSource::Unknown,
                neanderthal_matches: nea,
                denisovan_matches: den,
            });
        }
    }

    let callable_mb: f64 = callable
        .contigs
        .iter()
        .filter(|c| observations.contains_key(&c.contig))
        .flat_map(|c| c.callable_bp.iter())
        .map(|&b| b as f64)
        .sum::<f64>()
        / 1_000_000.0;
    let total_mb: f64 = segments.iter().map(|s| s.length_mb()).sum();
    let summary = ArchaicSummary {
        total_mb,
        pct_callable: if callable_mb > 0.0 {
            total_mb * 100.0 / callable_mb
        } else {
            0.0
        },
        callable_mb,
        neanderthal_mb: 0.0,
        denisovan_mb: 0.0,
        unknown_mb: total_mb,
        n_segments: segments.len(),
    };
    ArchaicSegmentResult { segments, summary }
}

/// How well a segment matches each archaic genome. Take the sites where a given archaic genome
/// carries the derived allele, and ask what fraction of those the subject also carries. The genome
/// with the highest fraction wins.
///
/// The condition is on the **genome**, and not on the subject, and that is what makes the measure
/// separate anything. The other way is to take the sites that the subject carries, and ask how
/// many an archaic genome shares. That scores about 100 % everywhere, and the background too,
/// because at an informative site some archaic holds the derived allele by construction.
///
/// Read the measure the way this function reads it. The background then sits at the genome-wide
/// rate of the subject, and a haplotype that came down from an ancestor sits far above that.
///
/// Returns `None` when no genome has enough called sites in the span to judge.
pub fn segment_concordance(
    panel: &ArchaicMarkerPanel,
    contig: &str,
    start: i64,
    end: i64,
    carried: &BTreeMap<(&str, i64), bool>,
    min_sites: usize,
) -> Option<f64> {
    let mut hits = [0usize; ARCHAIC_GENOMES.len()];
    let mut dens = [0usize; ARCHAIC_GENOMES.len()];
    for s in &panel.sites {
        if s.contig != contig || s.position < start || s.position > end {
            continue;
        }
        let subject_has = carried.get(&(contig, s.position)).copied().unwrap_or(false);
        for (i, call) in s.calls.iter().enumerate() {
            if call.carries_derived() {
                dens[i] += 1;
                hits[i] += usize::from(subject_has);
            }
        }
    }
    (0..ARCHAIC_GENOMES.len())
        .filter(|&i| dens[i] >= min_sites)
        .map(|i| hits[i] as f64 / dens[i] as f64)
        .fold(None, |best: Option<f64>, r| Some(best.map_or(r, |b| b.max(r))))
}

/// Which sites of the Tier A panel the subject carries the archaic-derived allele at.
///
/// A site with no variant record is hom-reference, so the subject does **not** carry the allele
/// there. The orientation of the panel itself removes the sites where the derived allele is the
/// reference base, before the data reaches here.
pub fn carried_panel_sites<'a>(
    panel: &'a ArchaicMarkerPanel,
    calls: &'a [SiteGenotype],
) -> BTreeMap<(&'a str, i64), bool> {
    let by_pos: BTreeMap<(&str, i64), &SiteGenotype> =
        calls.iter().map(|c| ((c.contig.as_str(), c.position), c)).collect();
    let mut out = BTreeMap::new();
    for s in &panel.sites {
        let Some((k, g)) = by_pos.get_key_value(&(s.contig.as_str(), s.position)) else {
            continue;
        };
        let carries = g.dosage > 0 && g.alternate_allele.starts_with(s.archaic_derived_allele);
        out.insert(*k, carries);
    }
    out
}

/// Drop segments that do not look like an inherited archaic haplotype.
///
/// This is the largest measured gain in quality: the **precision goes from 54 % to 90 %** on real
/// data. It works because it reads evidence that the segment caller never sees, which is the
/// question of which archaic genome carries what. That is truly new information, and not a second
/// reading of the same signal.
///
/// The code **keeps** a segment with too few sites to judge. The absence of evidence is not
/// evidence of a bad call, and to drop such a segment would punish a sparse region where nobody
/// looks.
///
/// Read the note about the populations in the module documentation. East Asian tracts match the
/// four sequenced archaic genomes less well than European tracts do, so this filter removes a
/// larger share of them. It raises the precision everywhere, and it makes the extent **less**
/// comparable between two populations.
pub fn filter_by_concordance(
    result: ArchaicSegmentResult,
    panel: &ArchaicMarkerPanel,
    calls: &[SiteGenotype],
    min_concordance: f64,
    min_sites: usize,
) -> ArchaicSegmentResult {
    let carried = carried_panel_sites(panel, calls);
    let kept: Vec<ArchaicSegment> = result
        .segments
        .into_iter()
        .filter(
            |seg| match segment_concordance(panel, &seg.contig, seg.start, seg.end, &carried, min_sites) {
                Some(c) => c >= min_concordance,
                None => true,
            },
        )
        .collect();
    let total_mb: f64 = kept.iter().map(|s| s.length_mb()).sum();
    let callable_mb = result.summary.callable_mb;
    ArchaicSegmentResult {
        summary: ArchaicSummary {
            total_mb,
            pct_callable: if callable_mb > 0.0 {
                total_mb * 100.0 / callable_mb
            } else {
                0.0
            },
            callable_mb,
            neanderthal_mb: 0.0,
            denisovan_mb: 0.0,
            unknown_mb: total_mb,
            n_segments: kept.len(),
        },
        segments: kept,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::archaic::{
        ArchaicCall, ArchaicPanelThresholds, ArchaicSite, CallableContig, ClassifyContig, PositionStream,
    };

    fn gmap(contig: &str, len: i32) -> GeneticMap {
        GeneticMap::uniform(1.0, &[(contig, len)])
    }

    fn callable(contig: &str, windows: usize) -> ArchaicCallable {
        ArchaicCallable {
            build: "chm13v2.0".into(),
            window_bp: 1_000,
            contigs: vec![CallableContig {
                contig: contig.into(),
                start: 0,
                callable_bp: vec![1_000u16; windows],
            }],
        }
    }

    fn obs(positions: &[(i64, bool)]) -> Vec<SiteObs> {
        positions
            .iter()
            .map(|&(position, carries)| SiteObs {
                position,
                carries,
                class: DiagnosticClass::SharedArchaic,
            })
            .collect()
    }

    /// A real tract looks like a run of carried sites, against a background of sites that the
    /// subject does not carry. That is the thing that this model exists to find.
    #[test]
    fn finds_a_run_of_carried_sites() {
        let mut sites: Vec<(i64, bool)> = (0..60).map(|i| (10_000 + i * 500, false)).collect();
        for s in sites.iter_mut().skip(20).take(20) {
            s.1 = true; // a 10 kb tract, 20 diagnostic sites, all carried
        }
        let mut m = BTreeMap::new();
        m.insert("chr21".to_string(), obs(&sites));
        // This test fixes the thresholds, and it does not take them from the defaults. It asks
        // whether the model finds a run at all, and it must not move when a calibrated default
        // moves. It broke once, when `min_posterior` rose to 0.98 and cut the edges off the run.
        // That was correct behaviour, and this test must not react to it.
        let cfg = MatchConfig {
            p_background: Some(0.13),
            p_archaic: Some(0.40),
            min_posterior: 0.80,
            min_sites: 5,
            min_segment_bp: 1_000,
            ..Default::default()
        };
        let r = call_from_observations(&m, &gmap("chr21", 60_000), &callable("chr21", 60), &cfg);
        assert_eq!(r.segments.len(), 1, "one tract expected, got {:?}", r.segments);
        let seg = &r.segments[0];
        assert!(seg.start >= 19_000 && seg.start <= 21_000, "start {} off", seg.start);
        assert!(seg.end >= 29_000 && seg.end <= 31_000, "end {} off", seg.end);
    }

    /// The failure that gated the density caller was calling tracts out of background noise. With
    /// no carried sites at all there is nothing to call, and the model must say so.
    #[test]
    fn calls_nothing_on_a_background_only_contig() {
        let sites: Vec<(i64, bool)> = (0..200).map(|i| (10_000 + i * 500, false)).collect();
        let mut m = BTreeMap::new();
        m.insert("chr21".to_string(), obs(&sites));
        let r = call_from_observations(
            &m,
            &gmap("chr21", 200_000),
            &callable("chr21", 200),
            &MatchConfig {
                p_background: Some(0.13),
                p_archaic: Some(0.40),
                ..Default::default()
            },
        );
        assert!(
            r.segments.is_empty(),
            "background should call nothing, got {:?}",
            r.segments
        );
    }

    /// Carried sites that lie apart, at the background rate, must not add up into a tract. That
    /// is how the density caller failed, in this observable.
    #[test]
    fn scattered_background_carriers_do_not_form_a_tract() {
        // 13 % carried, evenly spread: exactly the background rate, no run.
        let sites: Vec<(i64, bool)> = (0..300).map(|i| (10_000 + i * 500, i % 8 == 0)).collect();
        let mut m = BTreeMap::new();
        m.insert("chr21".to_string(), obs(&sites));
        let r = call_from_observations(
            &m,
            &gmap("chr21", 300_000),
            &callable("chr21", 300),
            &MatchConfig {
                p_background: Some(0.13),
                p_archaic: Some(0.40),
                ..Default::default()
            },
        );
        assert!(
            r.segments.is_empty(),
            "background-rate carriers formed {:?}",
            r.segments
        );
    }

    /// A site whose derived allele IS the reference base separates nothing, and a no-call there
    /// means the opposite of what it means elsewhere. The code must drop such a site, and must not
    /// count it.
    #[test]
    fn observations_drop_sites_where_reference_is_derived() {
        let classify = ArchaicClassify {
            build: "chm13v2.0".into(),
            contigs: vec![ClassifyContig {
                positions: PositionStream::encode("chr21", &[1_000, 2_000, 3_000]),
                derived: vec![b'A', b'C', b'G'],
                classes: vec![0, 1, 2],
            }],
        };
        let calls: BTreeMap<i64, &SiteGenotype> = BTreeMap::new();
        // The reference carries the derived allele at 2_000 only.
        let out = observations_for_contig(
            "chr21",
            &classify,
            &calls,
            |p| if p == 2_000 { Some(b'C') } else { Some(b'T') },
            &callable("chr21", 10),
            0.5,
        );
        assert_eq!(out.len(), 2, "the reference-derived site must be dropped");
        assert!(out.iter().all(|o| o.position != 2_000));
        assert!(out.iter().all(|o| !o.carries), "no calls means nothing carried");
    }

    fn panel_site(pos: i64, derived: char, calls: [ArchaicCall; 4]) -> ArchaicSite {
        ArchaicSite {
            contig: "chr21".into(),
            position: pos,
            reference_allele: 'T',
            alternate_allele: derived,
            archaic_derived_allele: derived,
            calls,
            diagnostic_class: DiagnosticClass::SharedArchaic,
            afr_freq: 0.0,
            grch37: None,
            grch38: None,
        }
    }

    fn genotype(pos: i64, alt: char, dosage: i32) -> SiteGenotype {
        SiteGenotype {
            name: String::new(),
            contig: "chr21".into(),
            position: pos,
            reference_allele: "T".into(),
            alternate_allele: alt.to_string(),
            ploidy: 2,
            dosage,
            gq: 60,
            depth: 30,
            ref_depth: 15,
            alt_depth: 15,
            pls: Vec::new(),
            gt: None,
            allele_depths: None,
        }
    }

    /// Read the concordance over each ARCHAIC GENOME, and not over each carried site. The other
    /// way is to take the sites that the subject carries, and ask how many some archaic
    /// shares. That scores about 100 % everywhere, and the background too, because at an
    /// informative site some archaic is derived by construction. Somebody wrote that version
    /// first, and it separated nothing.
    #[test]
    fn concordance_conditions_on_the_genome_not_the_subject() {
        use ArchaicCall::{HomAncestral as A, HomDerived as D};
        // Altai is derived at all 4 sites; Denisova at only the first.
        let panel = ArchaicMarkerPanel {
            build: "chm13v2.0".into(),
            thresholds: ArchaicPanelThresholds {
                max_afr_freq: 0.01,
                min_non_afr_freq: 0.0005,
            },
            sites: vec![
                panel_site(1_000, 'A', [D, A, A, D]),
                panel_site(2_000, 'A', [D, A, A, A]),
                panel_site(3_000, 'A', [D, A, A, A]),
                panel_site(4_000, 'A', [D, A, A, A]),
            ],
        };
        // The subject carries 3 of Altai's 4 → 0.75 for Altai, 1.0 for Denisova but on ONE site.
        let calls = vec![
            genotype(1_000, 'A', 1),
            genotype(2_000, 'A', 1),
            genotype(3_000, 'A', 1),
        ];
        let carried = carried_panel_sites(&panel, &calls);
        // min_sites = 3 excludes Denisova's single site, so Altai's 0.75 is the answer. Without
        // that floor a 1/1 genome would win every segment.
        let c = segment_concordance(&panel, "chr21", 0, 5_000, &carried, 3).expect("a score");
        assert!((c - 0.75).abs() < 1e-9, "expected Altai's 0.75, got {c}");
    }

    /// The code must KEEP a segment with too few sites to judge. The absence of evidence is not
    /// evidence of a bad call. To drop such a segment would punish a sparse region. Those are
    /// exactly the regions where a caller most needs the doubt to go in its favour.
    #[test]
    fn filter_keeps_segments_it_cannot_judge() {
        let panel = ArchaicMarkerPanel {
            build: "chm13v2.0".into(),
            thresholds: ArchaicPanelThresholds {
                max_afr_freq: 0.01,
                min_non_afr_freq: 0.0005,
            },
            sites: vec![panel_site(1_000, 'A', [ArchaicCall::HomDerived; 4])],
        };
        let seg = ArchaicSegment {
            contig: "chr21".into(),
            start: 500_000,
            end: 600_000, // no panel sites here at all
            posterior: 0.99,
            n_private: 40,
            source: ArchaicSource::Unknown,
            neanderthal_matches: 0,
            denisovan_matches: 0,
        };
        let r = ArchaicSegmentResult {
            summary: ArchaicSummary {
                total_mb: 0.1,
                pct_callable: 1.0,
                callable_mb: 10.0,
                neanderthal_mb: 0.0,
                denisovan_mb: 0.0,
                unknown_mb: 0.1,
                n_segments: 1,
            },
            segments: vec![seg],
        };
        let out = filter_by_concordance(r, &panel, &[], 0.9, 3);
        assert_eq!(out.segments.len(), 1, "an unjudgable segment must survive");
    }

    /// The filter's whole purpose: a segment that does not look like an inherited archaic haplotype
    /// goes, one that does stays. This is the 54 % -> 90 % precision lever.
    #[test]
    fn filter_drops_poorly_matching_segments() {
        use ArchaicCall::{HomAncestral as A, HomDerived as D};
        let mut sites = Vec::new();
        for i in 0..10 {
            sites.push(panel_site(1_000 + i * 100, 'A', [D, A, A, A])); // good segment
        }
        for i in 0..10 {
            sites.push(panel_site(50_000 + i * 100, 'A', [D, A, A, A])); // bad segment
        }
        let panel = ArchaicMarkerPanel {
            build: "chm13v2.0".into(),
            thresholds: ArchaicPanelThresholds {
                max_afr_freq: 0.01,
                min_non_afr_freq: 0.0005,
            },
            sites,
        };
        // Carries 9/10 in the first span, 1/10 in the second.
        let mut calls: Vec<SiteGenotype> = (0..9).map(|i| genotype(1_000 + i * 100, 'A', 1)).collect();
        calls.push(genotype(50_000, 'A', 1));

        let mk = |start: i64, end: i64| ArchaicSegment {
            contig: "chr21".into(),
            start,
            end,
            posterior: 0.99,
            n_private: 20,
            source: ArchaicSource::Unknown,
            neanderthal_matches: 0,
            denisovan_matches: 0,
        };
        let r = ArchaicSegmentResult {
            summary: ArchaicSummary {
                total_mb: 0.002,
                pct_callable: 1.0,
                callable_mb: 10.0,
                neanderthal_mb: 0.0,
                denisovan_mb: 0.0,
                unknown_mb: 0.002,
                n_segments: 2,
            },
            segments: vec![mk(900, 2_000), mk(49_900, 51_000)],
        };
        let out = filter_by_concordance(r, &panel, &calls, 0.7, 3);
        assert_eq!(out.segments.len(), 1, "the poorly-matching segment should go");
        assert_eq!(out.segments[0].start, 900);
        assert_eq!(
            out.summary.n_segments, 1,
            "the summary must be recomputed, not carried over"
        );
    }

    /// A no-call is hom-reference, so the subject does NOT carry the allele. A condition on "has
    /// a call" instead is what made an early version of this analysis report about 80 %, against a
    /// background of 4.3 %.
    #[test]
    fn a_missing_call_is_not_a_carrier() {
        let classify = ArchaicClassify {
            build: "chm13v2.0".into(),
            contigs: vec![ClassifyContig {
                positions: PositionStream::encode("chr21", &[1_000, 2_000]),
                derived: vec![b'A', b'A'],
                classes: vec![0, 0],
            }],
        };
        let carried = SiteGenotype {
            name: String::new(),
            contig: "chr21".into(),
            position: 1_000,
            reference_allele: "T".into(),
            alternate_allele: "A".into(),
            ploidy: 2,
            dosage: 1,
            gq: 60,
            depth: 30,
            ref_depth: 15,
            alt_depth: 15,
            pls: Vec::new(),
            gt: None,
            allele_depths: None,
        };
        let mut calls: BTreeMap<i64, &SiteGenotype> = BTreeMap::new();
        calls.insert(1_000, &carried);
        let out = observations_for_contig("chr21", &classify, &calls, |_| Some(b'T'), &callable("chr21", 10), 0.5);
        assert_eq!(out.len(), 2);
        assert!(out[0].carries, "a called derived allele carries");
        assert!(!out[1].carries, "an absent call is hom-reference, not a carrier");
    }
}
