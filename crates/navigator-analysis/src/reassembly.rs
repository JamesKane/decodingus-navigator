//! The haploid **local reassembly** resolver, which is Option B of the private-Y work, phase 1. It
//! is pure Rust, it needs no external tool, and it is clean on Windows and MSVC.
//!
//! The design is in `documents/design/haploid-reassembly-caller.md`. This module owns **Stages B to
//! E** over one active window. Those four stages are:
//!
//! - the selection of the reads, with a gate on the mapping quality and a dedup of the fragments;
//! - the candidate haplotypes, which in v1 is one for each SNV: the reference, against the
//!   reference with one substitution;
//! - the likelihood of a read against a haplotype, from a **PairHMM that knows the base
//!   qualities** (`bio::stats::pairhmm`);
//! - the haploid genotype, from the log-odds over all of the reads.
//!
//! `caller.rs` owns Stage A and Stage F. Stage A finds the active region, and that code already
//! tallies the counts at each position. Stage F turns a [`ReassemblyCall`] into a `VariantCall`.
//! This module does **no I/O**, and that is deliberate: a unit test can then run it on a synthetic
//! window.
//!
//! Here is why it exists. The pileup caller in `caller.rs` refuses a position whose pileup is near
//! 50/50, because it suspects a paralog artifact. See `is_paralogous`. At a Y locus with a
//! segmental duplication, or an ampliconic one, that throws away a *true* derived SNV. Reads from
//! a paralogous region map to the wrong place and bring the reference base onto the site.
//!
//! GATK resolves those by local reassembly and a PairHMM over the base qualities. This module is
//! the haploid-only equivalent. The POC in `examples/reassembly_probe.rs` showed it on WGS229: the
//! PairHMM over base qualities recovers the misaligned-reference sites where the crude
//! match-against-mismatch pileup gives a tie.
//!
//! v1 works on **one candidate SNV at a time**, with one alternate haplotype for each candidate
//! position. Linked variants and short indels, through a POA assembly of more than one haplotype,
//! are the v2 extension. See the design document. POA still serves here as an optional
//! cross-check for the caller.

use std::collections::HashMap;

use bio::alignment::pairwise::Aligner as PwAligner;
use bio::alignment::AlignmentOperation;
use bio::stats::pairhmm::{EmissionParameters, GapParameters, PairHMM, StartEndGapParameters, XYEmission};
use bio::stats::{LogProb, Prob};

/// Natural-log → Phred scale factor (`10 / ln 10`); `LogProb` is base-*e*.
const PHRED_PER_NAT: f64 = 4.342_944_819_032_518;

/// The controls of the reassembly resolver. The defaults are the start points that the POC
/// checked. The §Open-questions of the design document marks τ and the window size for a
/// calibration against the full truth set.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReassemblyParams {
    /// The code drops a read below this mapping quality, which is the GATK default. That is what
    /// removes a paralog read whose place is in doubt, and which otherwise looks like reference
    /// support at a high base quality.
    pub min_mapping_quality: u8,
    /// Minimum aggregate log-odds (nats) for a haploid DERIVED call; symmetric for ANCESTRAL.
    pub min_log_odds: f64,
    /// A DERIVED call needs this many fragments that support the alt allele, or more, after the
    /// dedup.
    pub min_alt_fragments: u32,
    /// A v2 option: build the alternate haplotype from the reads that support the alt allele. It
    /// is a majority consensus over the reference frame. See [`assemble_alt_haplotype`]. A linked
    /// variant that the true reads carry then does not count against them, in the comparison with
    /// the reference.
    ///
    /// **The default is off.** It helps the synthetic case with a linked variant. But on the real
    /// WGS229 data it moves a site that sits near 50/50, and it broke `chrY:4284195`. There is also
    /// no real truth site with a linked variant yet, so nobody can check the gain.
    ///
    /// Unit tests cover the mechanism, and you turn it on with this flag or with
    /// `NAVIGATOR_REASSEMBLY_ASSEMBLE=1`, until somebody checks it. The floor on the read
    /// likelihood below is the v2 gain that is on by default. See
    /// `haploid-reassembly-caller.md`.
    pub assemble_alt: bool,
    /// A v2 rule: drop a read whose best log-likelihood, over the reference haplotype and the alt
    /// one, is below this. Such a read matches *neither* local haplotype. It is a paralog, or junk
    /// from another locus. One mismatch costs about `-9 nats`, so `-90` accepts real divergence, at
    /// about 9 or 10 mismatches, before the code drops a read.
    pub min_read_loglik: f64,
}

impl Default for ReassemblyParams {
    fn default() -> Self {
        Self {
            min_mapping_quality: 20,
            min_log_odds: 2.0,
            min_alt_fragments: 2,
            assemble_alt: false,
            min_read_loglik: -90.0,
        }
    }
}

/// What a read observes at one candidate site (the base it carries there and that base's quality).
/// `None` in [`WindowRead::site_obs`] means the read does not span that candidate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SiteObs {
    pub base: u8,
    pub qual: u8,
}

/// A read projected onto the reference frame of the active window.
///
/// The caller builds it. That build is the CIGAR walk that gives the sequence in the window frame,
/// the quality of each base, and a [`SiteObs`] at each candidate. This module reads the
/// projection, so it stays free of I/O and a test can cover it.
#[derive(Debug, Clone)]
pub struct WindowRead {
    /// The identity of the fragment, which is the query name. A read and its mate share the name,
    /// and the dedup uses it.
    pub name: Vec<u8>,
    /// Window-frame bases (uppercase), for the whole-read PairHMM realignment.
    pub seq: Vec<u8>,
    /// The Phred quality of each base, in line with `seq`.
    pub quals: Vec<u8>,
    /// Mapping quality of the source record.
    pub mapq: u8,
    /// Observation at each candidate site; parallel to the `candidates` slice given to
    /// [`genotype_window`]. `None` = read does not cover that site.
    pub site_obs: Vec<Option<SiteObs>>,
}

/// A candidate variant position within the window (1-based reference coordinate).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Candidate {
    pub position: i64,
    pub ref_base: u8,
    pub alt_base: u8,
}

/// The haploid genotype the resolver assigns to a candidate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Zygosity {
    /// The alternate haplotype explains the reads (a private variant).
    Derived,
    /// The reference haplotype explains the reads (drop it).
    Ancestral,
    /// Neither side wins by `min_log_odds`. The data does not decide this site, so do not call
    /// it.
    Ambiguous,
}

/// A candidate with a genotype. The caller keeps a [`Zygosity::Derived`] call and turns it into a
/// `VariantCall`, in Stage F. The others come back too, so that a test and a diagnostic can see
/// the decision.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReassemblyCall {
    pub position: i64,
    pub ref_base: u8,
    pub alt_base: u8,
    /// The count of fragments that cover the site, after the MAPQ gate and the dedup of the
    /// mates.
    pub depth: u32,
    /// The count of those fragments whose base at the site is the alternate allele.
    pub alt_depth: u32,
    pub allele_fraction: f64,
    /// Aggregate `Σ ln P(read|alt) − ln P(read|ref)` (nats); >0 favours the alt haplotype.
    pub log_odds: f64,
    /// The confidence of the genotype that won, on the Phred scale. It is like a GQ.
    pub quality: f64,
    pub genotype: Zygosity,
}

/// Genotype every candidate in one active window against the reads projected onto it.
///
/// `ref_window` is the uppercase reference over the window; `window_start` is the 1-based reference
/// coordinate of `ref_window[0]`. Each read's `site_obs` must be parallel to `candidates`.
pub fn genotype_window(
    ref_window: &[u8],
    window_start: i64,
    candidates: &[Candidate],
    reads: &[WindowRead],
    params: &ReassemblyParams,
) -> Vec<ReassemblyCall> {
    let mut hmm = PairHMM::new(&GapParams);
    candidates
        .iter()
        .enumerate()
        .map(|(ci, cand)| genotype_candidate(&mut hmm, ref_window, window_start, ci, cand, reads, params))
        .collect()
}

/// Resolve one candidate. It builds the alt haplotype. It takes the reads that cover the site, and
/// dedups them. It scores each of those against the reference and the alt with the PairHMM. It
/// then gives a genotype from the log-odds over all of them.
fn genotype_candidate(
    hmm: &mut PairHMM,
    ref_window: &[u8],
    window_start: i64,
    ci: usize,
    cand: &Candidate,
    reads: &[WindowRead],
    params: &ReassemblyParams,
) -> ReassemblyCall {
    let off = (cand.position - window_start) as usize;

    // Stage B. Take the reads that clear the MAPQ gate and that cover this candidate. Then put a
    // mate pair that overlaps into one fragment, and keep the record whose base at the site has
    // the higher quality.
    let kept = dedup_spanning_fragments(reads, ci, params);

    // Stage C, the alternate haplotype. In v2, a POA builds it from the reads that support the
    // alt allele. A linked variant that those reads carry then does not count against them in the
    // comparison with the reference. When the assembly is degenerate, the code falls back to the
    // reference plus one substitution. That fallback is the v1 behaviour, so a simple site does
    // not change.
    let mut single_snv = ref_window.to_vec();
    if off < single_snv.len() {
        single_snv[off] = cand.alt_base;
    }
    let alt_hap = if params.assemble_alt {
        assemble_alt_haplotype(reads, &kept, ci, ref_window, off, cand.alt_base).unwrap_or(single_snv)
    } else {
        single_snv
    };

    // Stages D and E. Take the likelihood ratio of each fragment, and the vote of its base at the
    // site. The floor on the absolute likelihood drops a read that matches neither haplotype,
    // which is a paralog, or junk from another locus.
    let mut log_odds = 0.0f64;
    let mut depth = 0u32;
    let mut alt_depth = 0u32;
    for &ri in &kept {
        let read = &reads[ri];
        let lp_ref = hap_likelihood(hmm, &read.seq, &read.quals, ref_window);
        let lp_alt = hap_likelihood(hmm, &read.seq, &read.quals, &alt_hap);
        if (*lp_ref).max(*lp_alt) < params.min_read_loglik {
            continue; // matches neither local haplotype — paralog/junk, don't let it vote
        }
        log_odds += *lp_alt - *lp_ref;
        depth += 1;
        if read.site_obs[ci].map(|o| o.base) == Some(cand.alt_base) {
            alt_depth += 1;
        }
    }

    let allele_fraction = if depth > 0 {
        alt_depth as f64 / depth as f64
    } else {
        0.0
    };
    let genotype = if log_odds > params.min_log_odds && alt_depth >= params.min_alt_fragments {
        Zygosity::Derived
    } else if log_odds < -params.min_log_odds {
        Zygosity::Ancestral
    } else {
        Zygosity::Ambiguous
    };

    ReassemblyCall {
        position: cand.position,
        ref_base: cand.ref_base,
        alt_base: cand.alt_base,
        depth,
        alt_depth,
        allele_fraction,
        log_odds,
        quality: (log_odds.abs() * PHRED_PER_NAT).min(99.0),
        genotype,
    }
}

/// The reads that clear the MAPQ gate and that cover the candidate `ci`. A mate pair that overlaps
/// goes into one fragment, and the code keeps the record whose base at the site has the higher
/// quality.
///
/// It returns the read indices in the order of the fragment names. The assembly that follows is
/// then deterministic, where the order of a `HashMap` is not.
fn dedup_spanning_fragments(reads: &[WindowRead], ci: usize, params: &ReassemblyParams) -> Vec<usize> {
    let mut by_fragment: HashMap<&[u8], usize> = HashMap::new();
    for (ri, read) in reads.iter().enumerate() {
        if read.mapq < params.min_mapping_quality {
            continue;
        }
        let Some(Some(obs)) = read.site_obs.get(ci) else {
            continue; // does not span the site
        };
        by_fragment
            .entry(read.name.as_slice())
            .and_modify(|kept| {
                let kept_q = reads[*kept].site_obs[ci].map(|o| o.qual).unwrap_or(0);
                if obs.qual > kept_q {
                    *kept = ri;
                }
            })
            .or_insert(ri);
    }
    let mut kept: Vec<usize> = by_fragment.into_values().collect();
    kept.sort_by(|&a, &b| reads[a].name.cmp(&reads[b].name));
    kept
}

/// Build the alternate haplotype from the fragments that support the alt allele, which are the
/// ones whose base at the site is `alt_base`. The method is a **majority consensus over the
/// reference frame**.
///
/// It takes three things. The reference. Every position where a strict majority of the alt reads
/// that cover it agree on the same non-reference base. And the candidate substitution at
/// `site_off`. It returns `None` with fewer than two alt reads, and the caller then falls back to
/// the reference plus the SNV.
///
/// This is *not* a raw POA, and that is deliberate. A real read has ragged ends, and it covers the
/// window only in part. A POA over such reads gives a noisy consensus, and that consensus scores a
/// site near 50/50 wrongly. In a test it broke `chrY:4284195`.
///
/// The majority rule comes down to the reference plus the SNV when the alt reads carry no linked
/// variant that they agree on. So it never hurts a site with no linked context. And it still adds
/// a real linked variant, so that the true reads match cleanly. A short indel is v2b, through a
/// POA over the alt reads that the code confirmed.
fn assemble_alt_haplotype(
    reads: &[WindowRead],
    kept: &[usize],
    ci: usize,
    ref_window: &[u8],
    site_off: usize,
    alt_base: u8,
) -> Option<Vec<u8>> {
    let alt_reads: Vec<&WindowRead> = kept
        .iter()
        .map(|&ri| &reads[ri])
        .filter(|r| r.site_obs.get(ci).and_then(|o| *o).map(|o| o.base) == Some(alt_base))
        .collect();
    if alt_reads.len() < 2 {
        return None;
    }

    // Tally the bases of each alt read at each reference position, through a pairwise projection
    // onto the window.
    let mut counts = vec![[0u32; 4]; ref_window.len()];
    let mut cover = vec![0u32; ref_window.len()];
    for r in &alt_reads {
        project_read_onto_ref(&r.seq, ref_window, &mut counts, &mut cover);
    }

    // Reference + concordant (strict-majority, ≥2-read) non-reference substitutions.
    let mut hap = ref_window.to_vec();
    for pos in 0..ref_window.len() {
        if cover[pos] < 2 {
            continue;
        }
        let (bi, cnt) = argmax4(&counts[pos]);
        if cnt * 2 > cover[pos] && BASES[bi] != ref_window[pos].to_ascii_uppercase() {
            hap[pos] = BASES[bi];
        }
    }
    // The candidate substitution is the reason for this call, so force it in. Its column can sit
    // at exactly 50/50.
    if site_off < hap.len() {
        hap[site_off] = alt_base;
    }
    Some(hap)
}

const BASES: [u8; 4] = *b"ACGT";

fn base_index(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// argmax over the four base counts (ties keep the earlier base A<C<G<T).
fn argmax4(counts: &[u32; 4]) -> (usize, u32) {
    let mut bi = 0;
    for i in 1..4 {
        if counts[i] > counts[bi] {
            bi = i;
        }
    }
    (bi, counts[bi])
}

/// Add the bases of `seq` to the `counts` and `cover` tallies at each reference position. It aligns
/// `seq` to `ref_window` in a semiglobal way. Only a column that aligns as a match or a mismatch
/// counts. An insertion and a deletion do not.
fn project_read_onto_ref(seq: &[u8], ref_window: &[u8], counts: &mut [[u32; 4]], cover: &mut [u32]) {
    let score = |a: u8, b: u8| if a == b { 1i32 } else { -4i32 };
    let mut aligner = PwAligner::new(-5, -1, score);
    let aln = aligner.semiglobal(seq, ref_window);
    let mut xi = aln.xstart;
    let mut yi = aln.ystart;
    for op in &aln.operations {
        match op {
            AlignmentOperation::Match | AlignmentOperation::Subst => {
                if let (Some(&b), Some(c), Some(cv)) = (seq.get(xi), counts.get_mut(yi), cover.get_mut(yi)) {
                    if let Some(bi) = base_index(b) {
                        c[bi] += 1;
                        *cv += 1;
                    }
                }
                xi += 1;
                yi += 1;
            }
            AlignmentOperation::Del => yi += 1,
            AlignmentOperation::Ins => xi += 1,
            AlignmentOperation::Xclip(n) => xi += n,
            AlignmentOperation::Yclip(n) => yi += n,
        }
    }
}

/// The log-probability that `hap` gave `read`, which carries `quals`. It marginalises over the
/// alignments.
fn hap_likelihood(hmm: &mut PairHMM, read: &[u8], quals: &[u8], hap: &[u8]) -> LogProb {
    hmm.prob_related(&ReadHapEmission { read, quals, hap }, &Semiglobal, None)
}

// ---- base-quality-aware PairHMM emission model (POC-validated) --------------------------------

/// The error probability that a Phred score gives, clamped to Q2 and Q60. A match or a mismatch is
/// then never sure.
fn phred_err(q: u8) -> f64 {
    let q = q.clamp(2, 60) as f64;
    10f64.powf(-q / 10.0)
}

/// The emission. `x` is the read, which carries a quality at each base. `y` is the candidate
/// haplotype.
struct ReadHapEmission<'a> {
    read: &'a [u8],
    quals: &'a [u8],
    hap: &'a [u8],
}

impl EmissionParameters for ReadHapEmission<'_> {
    fn prob_emit_xy(&self, i: usize, j: usize) -> XYEmission {
        let err = phred_err(self.quals[i]);
        if self.read[i] == self.hap[j] {
            XYEmission::Match(LogProb::from(Prob(1.0 - err)))
        } else {
            XYEmission::Mismatch(LogProb::from(Prob(err / 3.0)))
        }
    }
    fn prob_emit_x(&self, _i: usize) -> LogProb {
        LogProb::ln_one() // insertion in read: base is real; cost is the gap-open prob
    }
    fn prob_emit_y(&self, _j: usize) -> LogProb {
        LogProb::ln_one() // deletion (gap in read): hap base emitted against a gap
    }
    fn len_x(&self) -> usize {
        self.read.len()
    }
    fn len_y(&self) -> usize {
        self.hap.len()
    }
}

/// GATK-ish affine gap model (indels rare relative to substitutions).
struct GapParams;
impl GapParameters for GapParams {
    fn prob_gap_x(&self) -> LogProb {
        LogProb::from(Prob(1e-4))
    }
    fn prob_gap_y(&self) -> LogProb {
        LogProb::from(Prob(1e-4))
    }
    fn prob_gap_x_extend(&self) -> LogProb {
        LogProb::from(Prob(0.1))
    }
    fn prob_gap_y_extend(&self) -> LogProb {
        LogProb::from(Prob(0.1))
    }
}

/// Semiglobal in the read. The offset at the start and at the end is free, so a cut at the edge of
/// the window costs nothing.
struct Semiglobal;
impl StartEndGapParameters for Semiglobal {
    fn free_start_gap_x(&self) -> bool {
        true
    }
    fn free_end_gap_x(&self) -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // A 40 bp reference window; candidate sits at its centre.
    const REF: &[u8] = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
    const WIN_START: i64 = 1000;
    const CAND_OFF: usize = 20;

    fn candidate() -> Candidate {
        Candidate {
            position: WIN_START + CAND_OFF as i64,
            ref_base: REF[CAND_OFF],
            alt_base: b'T', // REF[20] is 'A'; flip A>T
        }
    }

    /// A read whose window sequence carries `site_base` at the candidate offset; uniform quality.
    fn read(name: &str, site_base: u8, qual: u8, mapq: u8) -> WindowRead {
        read_muts(name, site_base, &[], qual, mapq)
    }

    /// The same as [`read`], and it also applies `muts`, which maps an offset to a base, to the
    /// sequence of the window. Use it to make reads that carry a linked variant, or, with many
    /// muts, paralog junk. `site_obs` holds the base at the candidate site alone, as the CIGAR
    /// walk of the caller would give it.
    fn read_muts(name: &str, site_base: u8, muts: &[(usize, u8)], qual: u8, mapq: u8) -> WindowRead {
        let mut seq = REF.to_vec();
        seq[CAND_OFF] = site_base;
        for &(off, b) in muts {
            seq[off] = b;
        }
        WindowRead {
            name: name.as_bytes().to_vec(),
            seq,
            quals: vec![qual; REF.len()],
            mapq,
            site_obs: vec![Some(SiteObs { base: site_base, qual })],
        }
    }

    fn call(reads: &[WindowRead]) -> ReassemblyCall {
        call_with(reads, &ReassemblyParams::default())
    }

    fn call_with(reads: &[WindowRead], params: &ReassemblyParams) -> ReassemblyCall {
        genotype_window(REF, WIN_START, &[candidate()], reads, params)
            .pop()
            .unwrap()
    }

    #[test]
    fn clean_derived_site_is_called() {
        // Twelve fragments, and all of them carry the alt allele. The call is strongly DERIVED.
        let reads: Vec<_> = (0..12).map(|i| read(&format!("r{i}"), b'T', 35, 60)).collect();
        let c = call(&reads);
        assert_eq!(c.genotype, Zygosity::Derived);
        assert_eq!(c.depth, 12);
        assert_eq!(c.alt_depth, 12);
        assert!(c.log_odds > 2.0, "log_odds {}", c.log_odds);
        assert_eq!(c.alt_base, b'T');
    }

    #[test]
    fn low_mapq_paralog_reference_reads_are_dropped_recovering_the_site() {
        // The case where the reference alignment is wrong. There are 8 clean alt fragments at
        // MAPQ 60, and 6 paralog fragments that carry the ref base but whose place is in doubt, at
        // MAPQ 5. The MAPQ gate drops the paralogs, so the site comes back as DERIVED. Without the
        // gate it would go out as a 50/50 rejection.
        let mut reads: Vec<_> = (0..8).map(|i| read(&format!("alt{i}"), b'T', 35, 60)).collect();
        reads.extend((0..6).map(|i| read(&format!("par{i}"), b'A', 35, 5)));
        let c = call(&reads);
        assert_eq!(c.genotype, Zygosity::Derived);
        assert_eq!(c.depth, 8, "paralog low-MAPQ reads must be excluded");
        assert_eq!(c.alt_depth, 8);
    }

    #[test]
    fn genuinely_balanced_high_quality_site_is_not_called() {
        // A test of the specificity. An even split of ref and alt fragments, all of high quality
        // and all placed well, decides nothing. The reassembly must NOT invent a call.
        let mut reads: Vec<_> = (0..6).map(|i| read(&format!("alt{i}"), b'T', 35, 60)).collect();
        reads.extend((0..6).map(|i| read(&format!("ref{i}"), b'A', 35, 60)));
        let c = call(&reads);
        assert_eq!(c.genotype, Zygosity::Ambiguous);
    }

    #[test]
    fn overlapping_mates_are_counted_once() {
        // Four separate alt fragments, plus a read and its mate, which share a name, and which
        // both cover the site. The dedup must put the mate pair into one fragment, so that the
        // depth reads 5 and not 6.
        let mut reads: Vec<_> = (0..4).map(|i| read(&format!("f{i}"), b'T', 35, 60)).collect();
        reads.push(read("pair", b'T', 20, 60)); // read
        reads.push(read("pair", b'T', 35, 60)); // its mate (higher qual → the kept one)
        let c = call(&reads);
        assert_eq!(c.depth, 5, "overlapping mates double-counted");
        assert_eq!(c.alt_depth, 5);
        assert_eq!(c.genotype, Zygosity::Derived);
    }

    #[test]
    fn all_reference_reads_are_ancestral() {
        let reads: Vec<_> = (0..10).map(|i| read(&format!("r{i}"), b'A', 35, 60)).collect();
        let c = call(&reads);
        assert_eq!(c.genotype, Zygosity::Ancestral);
        assert_eq!(c.alt_depth, 0);
        assert!(c.log_odds < -2.0, "log_odds {}", c.log_odds);
    }

    // ---- v2: assembled alt haplotype + read-likelihood floor -----------------------------------

    // Two linked SNVs the true reads carry alongside the derived allele (REF[10]='G', REF[30]='G').
    const LINKED: &[(usize, u8)] = &[(10, b'T'), (30, b'T')];

    #[test]
    fn assembled_alt_haplotype_lifts_confidence_on_linked_variant_site() {
        // The true reads are the majority, and they carry the derived allele PLUS two linked
        // variants. The reference reads are clean.
        //
        // Against a v1 alt haplotype, which is the reference plus one SNV, those linked variants
        // count against the true reads. The v2 haplotype, which a POA assembles, lets them match
        // cleanly. The call is then DERIVED, and it carries more confidence than in v1.
        let mut reads: Vec<_> = (0..10)
            .map(|i| read_muts(&format!("alt{i}"), b'T', LINKED, 35, 60))
            .collect();
        reads.extend((0..4).map(|i| read(&format!("ref{i}"), b'A', 35, 60)));

        let v1 = call_with(&reads, &ReassemblyParams::default()); // assemble_alt: false (ref+SNV)
        let v2 = call_with(
            &reads,
            &ReassemblyParams {
                assemble_alt: true,
                ..Default::default()
            },
        );
        assert_eq!(v1.genotype, Zygosity::Derived);
        assert_eq!(v2.genotype, Zygosity::Derived);
        assert!(
            v2.log_odds > v1.log_odds + 5.0,
            "assembly should raise confidence: v1 {} vs v2 {}",
            v1.log_odds,
            v2.log_odds
        );
    }

    #[test]
    fn paralog_junk_read_matching_neither_haplotype_is_filtered() {
        // Five clean reference reads, plus one "read" that carries the alt base and holds
        // mismatches *across the whole* window. That is a paralog fragment from another locus.
        //
        // The spread of those mismatches matters. The semiglobal PairHMM cuts a clean start and a
        // clean end off a read. So only mismatches that lie across the whole read make it match
        // neither haplotype. The floor on the likelihood must drop it. It must not raise the
        // depth, and it must not move the call away from ANCESTRAL.
        let junk_muts: Vec<(usize, u8)> = (0..REF.len())
            .step_by(2)
            .filter(|&k| k != CAND_OFF)
            .map(|k| (k, if REF[k] == b'A' { b'C' } else { b'A' }))
            .collect();
        let mut reads: Vec<_> = (0..5).map(|i| read(&format!("ref{i}"), b'A', 35, 60)).collect();
        reads.push(read_muts("junk", b'T', &junk_muts, 35, 60));

        let c = call(&reads);
        assert_eq!(c.depth, 5, "the paralog-junk read must be excluded from depth");
        assert_eq!(c.alt_depth, 0);
        assert_eq!(c.genotype, Zygosity::Ancestral);
    }

    #[test]
    fn assembly_falls_back_to_single_snv_when_alt_reads_are_too_few() {
        // One alt read, which is fewer than 2, can not start an assembly. So the code falls back
        // to the reference plus the SNV. With one alt fragment against ten reference reads,
        // the site stays ANCESTRAL, and there is no false call.
        let mut reads: Vec<_> = (0..10).map(|i| read(&format!("ref{i}"), b'A', 35, 60)).collect();
        reads.push(read("lone", b'T', 35, 60));
        let c = call(&reads);
        assert_eq!(c.genotype, Zygosity::Ancestral);
    }
}
