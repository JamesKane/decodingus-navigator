//! STR genotypes from aligned reads. The model works on a read that encloses the tract.
//!
//! An STR shorter than a read is the case here, and every Y-STR, and every forensic and
//! genealogical marker, qualifies. The reads that carry information are the ones that **enclose**
//! the whole repeat tract, plus clean sequence on both sides of it. GangSTR calls those the
//! "enclosing" class. Its Spanning, Flanking and FRR classes matter only for an expansion that is
//! longer than a read.
//!
//! At each enclosing read the code takes the observed repeat length off the **CIGAR**. It is
//! `tract_bp + (insertions − deletions) inside the tract`, measured against the known reference
//! allele. So it carries no systematic offset, where a count of the motif copies in a loose
//! feature region would.
//!
//! A geometric model of **PCR stutter** then turns the counts of the reads into a
//! maximum-likelihood genotype. That genotype is haploid on chrY, and diploid elsewhere.
//!
//! This is the core of HipSTR and GangSTR that a person can build and defend. It leaves out their
//! stutter EM, their HMM realignment, and their SNP phasing, and those are a later improvement.
//!
//! It trusts the CIGAR of the aligner inside a tight tract. And it lets the genotype, which is the
//! mode over the reads, absorb a read that the aligner put in the wrong place.

use std::path::Path;

use noodles::core::Region;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::RecordBuf;
use serde::{Deserialize, Serialize};

use crate::error::AnalysisError;
use crate::reader;
use crate::strref::StrLocus;

/// Tunables for the STR caller.
#[derive(Debug, Clone, Copy)]
pub struct StrCallerParams {
    /// The mapping quality that a read needs before the code uses it.
    pub min_mapping_quality: u8,
    /// Clean, indel-free reference bases required on each side of the tract to count a read.
    pub flank: i64,
    /// The depth of reads that enclose the tract that the code needs before it emits a genotype.
    pub min_depth: u32,
    /// `P(a read shows the true allele exactly)`. It is the probability of no stutter, and the
    /// HipSTR default is 0.9.
    pub no_stutter: f64,
    /// The geometric decay of the stutter size, at each extra repeat unit. A smaller value makes
    /// ±1 dominate.
    pub stutter_decay: f64,
}

impl Default for StrCallerParams {
    fn default() -> Self {
        Self {
            min_mapping_quality: 20,
            flank: 5,
            min_depth: 5,
            no_stutter: 0.9,
            stutter_decay: 0.5,
        }
    }
}

/// Confidence tier for a called STR genotype.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum StrConfidence {
    High,
    Medium,
    Low,
}

/// A called STR genotype at one locus.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct StrGenotype {
    pub contig: String,
    /// The 0-based start of the tract, in BED form. Use it to join back to the reference, or to a
    /// mapping that a vendor gives.
    pub start: i64,
    pub end: i64,
    pub period: u8,
    pub motif: String,
    /// HipSTR locus id (the result name until a vendor DYS mapping exists).
    pub name: String,
    pub ref_copies: f64,
    /// The called alleles, in **repeat copies**. There is one on a haploid contig, such as chrY,
    /// and one or two on a diploid one.
    pub alleles: Vec<i32>,
    /// The count of enclosing reads that the code used. It is the depth of the genotype.
    pub depth: u32,
    /// The fraction of the enclosing reads that match the called alleles exactly.
    pub concordance: f64,
    pub confidence: StrConfidence,
}

/// `ln P(the observed copies | the true allele copies)`, under the geometric stutter model.
///
/// The read shows the allele exactly, at a probability of `no_stutter`. Otherwise the size of the
/// deviation, in repeat units, follows a geometric distribution, and it is symmetric up and
/// down.
fn obs_lnlik(observed: i32, allele: i32, p: &StrCallerParams) -> f64 {
    if observed == allele {
        p.no_stutter.ln()
    } else {
        let d = (observed - allele).unsigned_abs() as i32;
        // (1-p0) split evenly up/down, geometric over magnitude d>=1.
        let prob = (1.0 - p.no_stutter) * 0.5 * (1.0 - p.stutter_decay) * p.stutter_decay.powi(d - 1);
        prob.max(1e-12).ln()
    }
}

/// Distinct observed copy numbers (the candidate alleles).
fn candidates(observed: &[i32]) -> Vec<i32> {
    let mut c: Vec<i32> = observed.to_vec();
    c.sort_unstable();
    c.dedup();
    c
}

/// Maximum-likelihood **haploid** allele: argmax over candidates of the summed stutter log-lik.
fn call_haploid(observed: &[i32], p: &StrCallerParams) -> Option<i32> {
    candidates(observed)
        .into_iter()
        .map(|a| (a, observed.iter().map(|&o| obs_lnlik(o, a, p)).sum::<f64>()))
        .max_by(|x, y| x.1.partial_cmp(&y.1).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(a, _)| a)
}

/// The maximum-likelihood **diploid** genotype `(A,B)`, where A<=B. It is the argmax over the
/// candidate pairs of the sum of `ln[½P(o|A) + ½P(o|B)]`. Each read comes from either allele with
/// equal probability.
fn call_diploid(observed: &[i32], p: &StrCallerParams) -> Option<(i32, i32)> {
    let cands = candidates(observed);
    let mut best: Option<((i32, i32), f64)> = None;
    for (i, &a) in cands.iter().enumerate() {
        for &b in &cands[i..] {
            let ll: f64 = observed
                .iter()
                .map(|&o| {
                    let pa = obs_lnlik(o, a, p).exp();
                    let pb = obs_lnlik(o, b, p).exp();
                    (0.5 * pa + 0.5 * pb).max(1e-12).ln()
                })
                .sum();
            if best.as_ref().map_or(true, |(_, bll)| ll > *bll) {
                best = Some(((a, b), ll));
            }
        }
    }
    best.map(|(g, _)| g)
}

/// The repeat copies that one enclosing read shows at `locus`. It is `None` when that read does
/// not enclose the tract cleanly. To do so it must sit on `flank` bp of reference with no indel,
/// on both sides.
///
/// The code takes the length off the CIGAR, as `tract_bp + insertions − deletions` inside the
/// tract, divided by the period.
fn observed_copies(ops: &[(Kind, usize)], aln_start: i64, locus: &StrLocus, flank: i64) -> Option<i32> {
    // A HipSTR tract INCLUDES its end. So ref_copies = (end - start + 1)/period. For example,
    // Y:2795644-2795670 at period 4 gives 27/4 = 6.75. The 1-based tract is then [start+1, end+1],
    // and its length is end-start+1. A read that matches the reference then measures exactly
    // ref_copies, with no systematic offset.
    let (ts, te) = (locus.start + 1, locus.end + 1); // 1-based inclusive tract
    let period = locus.period as i64;
    if period == 0 {
        return None;
    }
    let mut ref_pos = aln_start;
    let mut ins_in_tract: i64 = 0;
    let mut del_in_tract: i64 = 0;
    for &(kind, len) in ops {
        let len = len as i64;
        match (kind.consumes_reference(), kind.consumes_read()) {
            (true, true) => ref_pos += len, // M/=/X
            (false, true) => {
                // Insertion at the boundary before `ref_pos`.
                let b = ref_pos;
                if (b > ts - flank && b < ts) || (b > te + 1 && b <= te + flank) {
                    return None; // indel in a flank → not a clean anchor
                }
                if b >= ts && b <= te + 1 {
                    ins_in_tract += len;
                }
            }
            (true, false) => {
                // Deletion spans [ref_pos, ref_pos+len-1].
                let (ds, de) = (ref_pos, ref_pos + len - 1);
                let flank_hit = overlap(ds, de, ts - flank, ts - 1) > 0 || overlap(ds, de, te + 1, te + flank) > 0;
                if flank_hit {
                    return None;
                }
                del_in_tract += overlap(ds, de, ts, te);
                ref_pos += len;
            }
            (false, false) => {} // S/H/P — soft/hard clip, pad: no reference advance
        }
    }
    let observed_len = (te - ts + 1) + ins_in_tract - del_in_tract;
    if observed_len <= 0 {
        return None;
    }
    Some(((observed_len as f64) / period as f64).round() as i32)
}

/// Length of the overlap of `[a0,a1]` and `[b0,b1]` (both inclusive), >= 0.
fn overlap(a0: i64, a1: i64, b0: i64, b1: i64) -> i64 {
    (a1.min(b1) - a0.max(b0) + 1).max(0)
}

fn read_passes(r: &RecordBuf, min_mapq: u8) -> bool {
    let f = r.flags();
    !f.is_unmapped()
        && !f.is_secondary()
        && !f.is_supplementary()
        && !f.is_duplicate()
        && r.mapping_quality().map(|m| m.get()).unwrap_or(0) >= min_mapq
}

/// Build a genotype from the counts of the reads that enclose a locus. A `ploidy` of 1 is haploid,
/// as on chrY. Any other value is diploid.
fn genotype_from_counts(locus: &StrLocus, counts: &[i32], ploidy: u8, p: &StrCallerParams) -> Option<StrGenotype> {
    let depth = counts.len() as u32;
    if depth < p.min_depth {
        return None;
    }
    let alleles: Vec<i32> = if ploidy == 1 {
        vec![call_haploid(counts, p)?]
    } else {
        let (a, b) = call_diploid(counts, p)?;
        if a == b {
            vec![a]
        } else {
            vec![a, b]
        }
    };
    let matching = counts.iter().filter(|&&c| alleles.contains(&c)).count();
    let concordance = matching as f64 / depth as f64;
    let confidence = match (depth, concordance) {
        (d, c) if d >= 10 && c >= 0.7 => StrConfidence::High,
        (d, c) if d >= p.min_depth && c >= 0.5 => StrConfidence::Medium,
        _ => StrConfidence::Low,
    };
    Some(StrGenotype {
        contig: locus.contig.clone(),
        start: locus.start,
        end: locus.end,
        period: locus.period,
        motif: locus.motif.clone(),
        name: locus.name.clone(),
        ref_copies: locus.ref_copies,
        alleles,
        depth,
        concordance,
        confidence,
    })
}

/// Genotype every locus in `loci` from `bam`, in one streaming pass. Every locus must lie on
/// `contig`, and `loci` must come in order of its start.
///
/// Each read gives its observed copy number to every locus that it encloses cleanly. A `ploidy` of
/// 1 is haploid, as on chrY. A CRAM needs `reference`.
pub fn genotype_str_loci(
    bam: &Path,
    contig: &str,
    loci: &[StrLocus],
    ploidy: u8,
    params: &StrCallerParams,
    reference: Option<&Path>,
) -> Result<Vec<StrGenotype>, AnalysisError> {
    if loci.is_empty() {
        return Ok(Vec::new());
    }
    let (header, mut idx) = reader::open_indexed(bam, reference)?;
    let region: Region = contig
        .parse()
        .map_err(|_| AnalysisError::Message(format!("bad region for contig {contig}")))?;
    // The loci come in order of their start. Collect the observed counts at each locus.
    let starts: Vec<i64> = loci.iter().map(|l| l.start).collect();
    let mut counts: Vec<Vec<i32>> = vec![Vec::new(); loci.len()];

    for result in idx.query(&header, &region)? {
        let record = result?;
        if !read_passes(&record, params.min_mapping_quality) {
            continue;
        }
        let Some(start) = record.alignment_start().map(|p| p.get() as i64) else {
            continue;
        };
        let ops: Vec<(Kind, usize)> = record.cigar().as_ref().iter().map(|op| (op.kind(), op.len())).collect();
        let ref_span: i64 = ops
            .iter()
            .filter(|(k, _)| k.consumes_reference())
            .map(|(_, l)| *l as i64)
            .sum();
        let aend = start + ref_span - 1;

        // Loci this read could enclose: tract_start (1-based = start+1) anchored >= flank from the
        // read's left edge, tract_end anchored <= flank from its right edge. Sorted by start, so
        // binary-search the first candidate and iterate while still in range.
        let lo = starts.partition_point(|&s| s + 1 - params.flank < start);
        for (i, locus) in loci.iter().enumerate().skip(lo) {
            if locus.start + 1 - params.flank < start {
                continue;
            }
            if locus.end + 1 + params.flank > aend {
                if locus.start + 1 - params.flank > aend {
                    break; // past this read's reach (sorted) — no later locus can fit either
                }
                continue;
            }
            if let Some(c) = observed_copies(&ops, start, locus, params.flank) {
                counts[i].push(c);
            }
        }
    }

    Ok(loci
        .iter()
        .zip(counts.iter())
        .filter_map(|(locus, c)| genotype_from_counts(locus, c, ploidy, params))
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn locus() -> StrLocus {
        // HipSTR end-inclusive: tract bp = end-start+1 = 20, period 4 → 5.0 ref copies, motif GATA.
        // 1-based tract = [101, 120].
        StrLocus {
            contig: "Y".into(),
            start: 100,
            end: 119,
            period: 4,
            ref_copies: 5.0,
            name: "L".into(),
            motif: "GATA".into(),
        }
    }

    #[test]
    fn haploid_picks_the_modal_allele_through_stutter() {
        let p = StrCallerParams::default();
        // 12 reads at 11 copies, 2 stutter at 10, 1 at 12 → haploid call 11.
        let mut obs = vec![11; 12];
        obs.extend([10, 10, 12]);
        assert_eq!(call_haploid(&obs, &p), Some(11));
    }

    #[test]
    fn diploid_resolves_a_heterozygote() {
        let p = StrCallerParams::default();
        // Balanced 13 and 16 with a little ±1 stutter → (13,16).
        let mut obs = vec![13; 10];
        obs.extend(vec![16; 10]);
        obs.extend([12, 14, 15, 17]);
        assert_eq!(call_diploid(&obs, &p), Some((13, 16)));
    }

    #[test]
    fn observed_copies_reads_indel_off_the_cigar() {
        let l = locus(); // tract 101..120 (20 bp), period 4, flank 5
                         // Read aligned from pos 90, all-match across a wide span → exactly ref (20bp/4 = 5 copies).
        let ops = vec![(Kind::Match, 60)];
        assert_eq!(observed_copies(&ops, 90, &l, 5), Some(5));

        // A 4bp insertion inside the tract (one extra GATA unit): 24bp/4 = 6 copies.
        // CIGAR: 20M (pos90..109) , 4I (insertion before ref 110, inside tract), 20M.
        let ops = vec![(Kind::Match, 20), (Kind::Insertion, 4), (Kind::Match, 20)];
        assert_eq!(observed_copies(&ops, 90, &l, 5), Some(6));

        // A 4bp deletion inside the tract → 16bp/4 = 4 copies.
        let ops = vec![(Kind::Match, 25), (Kind::Deletion, 4), (Kind::Match, 20)];
        assert_eq!(observed_copies(&ops, 90, &l, 5), Some(4));

        // An indel in the left flank → rejected (not a clean anchor).
        let ops = vec![(Kind::Match, 7), (Kind::Deletion, 2), (Kind::Match, 50)];
        assert_eq!(observed_copies(&ops, 90, &l, 5), None);
    }
}
