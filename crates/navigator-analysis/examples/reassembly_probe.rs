//! A proof of the reassembly caller, which is Option B. It is pure Rust, over the `bio` crate. It
//! needs no external tool, and it is clean on Windows.
//!
//! It has two purposes. The first is to check the pure-Rust library stack from end to end, on a
//! real CRAM. The second is to find out what the resolver needs.
//!
//! It builds a reference haplotype and an alternate one, over a window around `pos`. It realigns
//! every read that covers the window to both, with `bio::alignment::pairwise`. It drops a read
//! that aligns poorly to both, which is paralog junk that the aligner put in the wrong place. It
//! then assembles the reads that stay with a POA, as a check, through `bio::alignment::poa`.
//!
//! The step that carries the weight comes last. It scores each read that covers the window
//! against both haplotypes, with a **PairHMM that knows the base qualities**, from
//! `bio::stats::pairhmm`. It then takes a genotype from the log-likelihood ratio over all of the
//! reads.
//!
//! **What it found.** The check ran on the real WGS229 CHM13 CRAM, against the own gVCF of GATK,
//! which calls all five sites DERIVED.
//!
//! A plain realignment resolves the *clean* control sites. But it gives a TIE at the marginal
//! misaligned-reference sites that the pileup caller misses. At 4284195 the crude score is 10
//! against 10.
//!
//! The PairHMM over base qualities breaks those ties, and it recovers **4 of the 5**. Those are
//! 3318203 and 16652092, which are the controls; 4284195, where GATK gives AD 9,10, GQ 44 and
//! MQRankSum -3.55, and which is the textbook misaligned-reference case; and 11191589. All four
//! come out DERIVED, and they match GATK.
//!
//! The one that it misses is 20973395. There, paralog reference reads pass the MQ≥20 gate. GATK
//! drops them in its selection of the fragments and reads of the active region, and its own DP
//! falls from 7 to 5 there. That selection is the one part that the full caller still needs.
//!
//! So the Option B caller has four parts. Find the active region. Assemble with a POA. Score each
//! **read against each haplotype with a PairHMM**. Then genotype.
//!
//! This probe shows two things. The whole pure-Rust stack builds and runs on a real CRAM. And the
//! PairHMM is what breaks the tie.
//!
//! ```text
//! cargo run --release --example reassembly_probe -p navigator-analysis -- \
//!     <cram> <ref.fa> chrY <pos[,pos,...]> [window=40]
//! ```

use std::path::Path;

use bio::alignment::pairwise::{Aligner as PwAligner, Scoring};
use bio::alignment::poa::Aligner as PoaAligner;
use bio::alignment::AlignmentOperation;
use bio::stats::pairhmm::{EmissionParameters, GapParameters, PairHMM, StartEndGapParameters, XYEmission};
use bio::stats::{LogProb, Prob};
use navigator_analysis::reader::{open_indexed, read_contig_sequence};
use noodles::core::Region;

fn base_index(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// GATK-default minimum mapping quality: excludes ambiguously-placed paralog reads, which at Y
/// segmental-duplication / ampliconic loci masquerade as high-base-quality reference support.
const MIN_MAPQ: u8 = 20;

/// One read that covers the window. It holds the sequence in the window frame, the Phred quality
/// of each of its bases, its mapping quality, and whether it covers the site.
struct WinRead {
    bases: Vec<u8>,
    quals: Vec<u8>,
    mapq: u8,
    covers_pos: bool,
}

/// The reads inside `[lo, hi]`, as [`WinRead`] values, and the raw A/C/G/T pileup at `pos`.
fn window_reads(cram: &Path, refp: &Path, contig: &str, pos: i64, lo: i64, hi: i64) -> (Vec<WinRead>, [u32; 4]) {
    let (header, mut reader) = open_indexed(cram, Some(refp)).expect("open cram");
    let region: Region = format!("{contig}:{lo}-{hi}").parse().expect("region");
    let mut pile = [0u32; 4];
    let mut reads: Vec<WinRead> = Vec::new();
    for result in reader.query(&header, &region).expect("query") {
        let record = result.expect("record");
        let flags = record.flags();
        if flags.is_secondary() || flags.is_supplementary() || flags.is_duplicate() || flags.is_unmapped() {
            continue;
        }
        let Some(start) = record.alignment_start().map(|p| p.get() as i64) else {
            continue;
        };
        let mapq = record.mapping_quality().map(|m| m.get()).unwrap_or(0);
        let seq = record.sequence();
        let seqb = seq.as_ref();
        let quals = record.quality_scores();
        let qualb = quals.as_ref();
        let mut ref_pos = start;
        let mut qoff = 0usize;
        let mut win: Vec<u8> = Vec::new();
        let mut winq: Vec<u8> = Vec::new();
        let mut covers_pos = false;
        for op in record.cigar().as_ref() {
            let kind = op.kind();
            let len = op.len();
            match (kind.consumes_reference(), kind.consumes_read()) {
                (true, true) => {
                    for i in 0..len {
                        let rp = ref_pos + i as i64;
                        if rp >= lo && rp <= hi {
                            if let Some(&b) = seqb.get(qoff + i) {
                                win.push(b.to_ascii_uppercase());
                                winq.push(qualb.get(qoff + i).copied().unwrap_or(30));
                                if rp == pos {
                                    covers_pos = true;
                                    if let Some(bi) = base_index(b) {
                                        pile[bi] += 1;
                                    }
                                }
                            }
                        }
                    }
                    ref_pos += len as i64;
                    qoff += len;
                }
                (true, false) => ref_pos += len as i64,
                (false, true) => {
                    // Keep the bases of an insertion that lie inside the window, so that an
                    // indel haplotype survives.
                    if kind == noodles::sam::alignment::record::cigar::op::Kind::Insertion
                        && ref_pos > lo
                        && ref_pos <= hi
                    {
                        for i in 0..len {
                            if let Some(&b) = seqb.get(qoff + i) {
                                win.push(b.to_ascii_uppercase());
                                winq.push(qualb.get(qoff + i).copied().unwrap_or(30));
                            }
                        }
                    }
                    qoff += len;
                }
                (false, false) => {}
            }
        }
        // Keep reads that carry enough window sequence to anchor a realignment.
        if win.len() >= 30 {
            reads.push(WinRead {
                bases: win,
                quals: winq,
                mapq,
                covers_pos,
            });
        }
    }
    (reads, pile)
}

/// POA-assemble the window reads into a consensus haplotype (heaviest-bundle).
fn poa_consensus(reads: &[Vec<u8>]) -> Vec<u8> {
    let scoring = Scoring::new(-4, -2, |a: u8, b: u8| if a == b { 2 } else { -4 });
    let mut aligner = PoaAligner::new(scoring, &reads[0]);
    for r in &reads[1..] {
        aligner.global(r).add_to_graph();
    }
    aligner.consensus()
}

/// The base that the consensus carries at the reference coordinate `pos`. It aligns the consensus,
/// as the query, to the reference window `win_ref`, in a semiglobal way. That window starts at the
/// reference coordinate `win_start`.
fn consensus_base_at(consensus: &[u8], win_ref: &[u8], win_start: i64, pos: i64) -> char {
    let score = |a: u8, b: u8| if a == b { 1i32 } else { -4i32 };
    let mut aligner = PwAligner::new(-5, -1, score);
    let aln = aligner.semiglobal(consensus, win_ref);
    if std::env::var("PROBE_DEBUG").is_ok() {
        eprintln!(
            "DEBUG consensus.len={} win_ref.len={} xstart={} ystart={} xend={} yend={} score={}",
            consensus.len(),
            win_ref.len(),
            aln.xstart,
            aln.ystart,
            aln.xend,
            aln.yend,
            aln.score
        );
        eprintln!("  consensus raw[..20]: {:?}", &consensus[..consensus.len().min(20)]);
        eprintln!(
            "  consensus str[..40]: {}",
            String::from_utf8_lossy(&consensus[..consensus.len().min(40)])
        );
        eprintln!(
            "  win_ref   str[..40]: {}",
            String::from_utf8_lossy(&win_ref[..win_ref.len().min(40)])
        );
    }
    let mut xi = aln.xstart; // consensus index
    let mut yi = aln.ystart; // win_ref index (ref coord = win_start + yi)
    for op in &aln.operations {
        match op {
            AlignmentOperation::Match | AlignmentOperation::Subst => {
                if win_start + yi as i64 == pos {
                    return consensus.get(xi).map(|&b| b as char).unwrap_or('?');
                }
                xi += 1;
                yi += 1;
            }
            AlignmentOperation::Del => {
                // gap in consensus: reference has a base the consensus lacks
                if win_start + yi as i64 == pos {
                    return '-';
                }
                yi += 1;
            }
            AlignmentOperation::Ins => xi += 1,
            AlignmentOperation::Xclip(n) => xi += n,
            AlignmentOperation::Yclip(n) => yi += n,
        }
    }
    '?'
}

// ---- base-quality-aware PairHMM: P(read | haplotype) ------------------------------------------
//
// This is what breaks a tie, and the crude alignment score can not. Each base of a read votes for
// the ref allele, or for the alt allele, at `pos`. The Phred quality of that base sets the weight
// of its vote. A Q40 base that does not match a haplotype costs about 10^-4, and a Q10 base costs
// only about 10^-1. A noisy base can then not out-vote a clean one.
//
// The sum of the log-likelihood ratio, over every read that covers the site, is how GATK does
// this. Its HaplotypeCaller resolves a pileup with a wrong reference alignment that way.

/// The error probability that a Phred score gives, clamped to a correct band: Q>0, and never a
/// match or a mismatch that is sure.
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
        LogProb::ln_one() // insertion in read: the base is real; cost is in the gap-open prob
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

/// The log-probability that `hap` gave `read`, which carries `quals`. It marginalises over the
/// alignments.
fn hap_likelihood(hmm: &mut PairHMM, read: &[u8], quals: &[u8], hap: &[u8]) -> LogProb {
    hmm.prob_related(&ReadHapEmission { read, quals, hap }, &Semiglobal, None)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 5 {
        eprintln!("usage: reassembly_probe <cram> <ref.fa> <contig> <pos[,pos,...]> [window=150]");
        std::process::exit(2);
    }
    let cram = Path::new(&args[1]);
    let refp = Path::new(&args[2]);
    let contig = &args[3];
    let positions: Vec<i64> = args[4].split(',').filter_map(|s| s.trim().parse().ok()).collect();
    // Small window so short reads fully span it (homologous POA fragments); ±40 bp fits ~150 bp reads.
    let window: i64 = args.get(5).and_then(|s| s.parse().ok()).unwrap_or(40);

    let mut refseq = read_contig_sequence(refp, contig).expect("read contig");
    refseq.iter_mut().for_each(|b| *b = b.to_ascii_uppercase()); // FASTA may be soft-masked (lowercase)

    let charb = |i: usize| ['A', 'C', 'G', 'T'][i];
    // realign = crude match/mismatch vote (ref/alt/dropped); pairHMM = base-quality vote (ref/alt) +
    // aggregate log-odds. The pairHMM column is the one that resolves the marginal sites.
    println!(
        "{:>10} {:>3} {:>13}  {:>8}   {:>7}   {:>11}   verdict",
        "pos", "ref", "A/C/G/T", "alt@frac", "realign", "pairHMM lo"
    );
    for &pos in &positions {
        let lo = (pos - window).max(1);
        let hi = pos + window;
        let ref_base = refseq[(pos - 1) as usize] as char;
        let (reads, pile) = window_reads(cram, refp, contig, pos, lo, hi);
        let total: u32 = pile.iter().sum();
        // Candidate alt = the most common non-reference base at pos.
        let ref_i = base_index(ref_base as u8).unwrap_or(0);
        let alt_i = (0..4).filter(|&i| i != ref_i).max_by_key(|&i| pile[i]).unwrap_or(ref_i);
        let alt_base = charb(alt_i);
        let alt_frac = if total > 0 {
            (total - pile[ref_i]) as f64 / total as f64
        } else {
            0.0
        };

        // Reference vs alternate haplotype over the window (alt = ref with the SNV at pos).
        let win_ref: Vec<u8> = refseq[(lo - 1) as usize..(hi as usize).min(refseq.len())].to_vec();
        let mut win_alt = win_ref.clone();
        let off = (pos - lo) as usize;
        if off < win_alt.len() {
            win_alt[off] = alt_base as u8;
        }

        // Realign each read that covers the window to both haplotypes. Drop a read that aligns
        // poorly to *both*, because it carries mismatches beyond the site. That read is a paralog
        // that the aligner put in the wrong place, or junk, and it did not come from this locus.
        //
        // For the reads that stay, also score P(read|ref-hap) against P(read|alt-hap), with the
        // PairHMM.
        let mut hmm = PairHMM::new(&GapParams);
        let (mut ref_supp, mut alt_supp, mut dropped) = (0u32, 0u32, 0u32);
        let (mut hmm_ref, mut hmm_alt) = (0u32, 0u32);
        let mut lowmq = 0u32;
        let mut logodds = 0.0f64; // Σ ln P(read|alt) − ln P(read|ref); >0 ⇒ alt haplotype favoured
        let mut kept: Vec<Vec<u8>> = Vec::new();
        for r in &reads {
            if !r.covers_pos {
                continue; // can't distinguish ref from alt without the site
            }
            if r.mapq < MIN_MAPQ {
                lowmq += 1; // ambiguously-placed paralog read — excluded like GATK does
                continue;
            }
            let win = &r.bases;
            let s_ref = realign_score(win, &win_ref);
            let s_alt = realign_score(win, &win_alt);
            let best = s_ref.max(s_alt);
            // A clean read scores ~+len; allow ~3 mismatches (errors/nearby real variants).
            let floor = win.len() as i32 - 15;
            if best < floor {
                dropped += 1;
                continue;
            }
            kept.push(win.clone());
            match s_alt.cmp(&s_ref) {
                std::cmp::Ordering::Greater => alt_supp += 1,
                std::cmp::Ordering::Less => ref_supp += 1,
                std::cmp::Ordering::Equal => {}
            }
            // The likelihood ratio that knows the base qualities. It is what breaks the tie.
            let lp_ref = hap_likelihood(&mut hmm, win, &r.quals, &win_ref);
            let lp_alt = hap_likelihood(&mut hmm, win, &r.quals, &win_alt);
            logodds += *lp_alt - *lp_ref;
            match lp_alt.partial_cmp(&lp_ref) {
                Some(std::cmp::Ordering::Greater) => hmm_alt += 1,
                Some(std::cmp::Ordering::Less) => hmm_ref += 1,
                _ => {}
            }
        }
        // The POA consensus of the reads that stay, as a cross-check. It is for information.
        let cons_note = if kept.len() >= 2 {
            let cons = poa_consensus(&kept);
            let cb = consensus_base_at(&cons, &win_ref, lo, pos);
            format!(" [POA={cb}]")
        } else {
            String::new()
        };

        // Genotype from the aggregate PairHMM log-odds: the haploid site is DERIVED when the alt
        // haplotype is the more likely explanation of the reads overall.
        let verdict = if logodds > 2.0 {
            format!("DERIVED {ref_base}>{alt_base}{cons_note}")
        } else if logodds < -2.0 {
            format!("ancestral{cons_note}")
        } else {
            format!("ambiguous{cons_note}")
        };
        println!(
            "{pos:>10} {ref_base:>3} {:>13}  {alt_base}@{alt_frac:.2}   {ref_supp:>2}/{alt_supp:>2}/{dropped:<2}   {hmm_ref:>2}/{hmm_alt:<2} {logodds:>+7.1}  mq-{lowmq:<2} {verdict}",
            format!("{}/{}/{}/{}", pile[0], pile[1], pile[2], pile[3])
        );
    }
}

/// Semiglobal realignment score of a read window against a haplotype (match +1, mismatch/gap penalised).
fn realign_score(read: &[u8], hap: &[u8]) -> i32 {
    let score = |a: u8, b: u8| if a == b { 1i32 } else { -4i32 };
    let mut aligner = PwAligner::new(-5, -1, score);
    aligner.semiglobal(read, hap).score
}
