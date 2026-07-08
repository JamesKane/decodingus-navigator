//! Prove-out for the reassembly caller (Option B) — pure Rust (bio), no external tools, Windows-clean.
//!
//! Purpose: (1) validate the pure-Rust library stack end-to-end on a real CRAM, and (2) scope what
//! the resolver actually needs. It builds a reference and an alternate haplotype over a window around
//! `pos`, realigns every spanning read to both (`bio::alignment::pairwise`), drops reads that align
//! poorly to both (misaligned paralog junk), POA-assembles the survivors as a check
//! (`bio::alignment::poa`), and reports ref-vs-alt support.
//!
//! Finding: plain realignment resolves the *clean* cases (controls) but NOT the marginal misaligned-ref
//! ~50/50 sites Navigator's pileup caller misses — those stay balanced (e.g. 10/10). GATK's own gVCF
//! calls them only marginally (AD 9,10 / GQ 44), so the tie-breaker is **base-quality-aware PairHMM
//! likelihood** (`bio::stats::pairhmm`), not the crude alignment score. So the full Option B caller =
//! active-region detection + POA assembly + **read-vs-haplotype PairHMM** + genotyping. This probe
//! proves the stack compiles/runs and pinpoints PairHMM as the load-bearing piece to add next.
//!
//!   cargo run --release --example reassembly_probe -p navigator-analysis -- \
//!       <cram> <ref.fa> chrY <pos[,pos,...]> [window=40]

use std::path::Path;

use bio::alignment::pairwise::{Aligner as PwAligner, Scoring};
use bio::alignment::poa::Aligner as PoaAligner;
use bio::alignment::AlignmentOperation;
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

/// Reads overlapping `[lo, hi]` as `(window_bases, covers_pos)`, plus the raw A/C/G/T pileup at `pos`.
fn window_reads(cram: &Path, refp: &Path, contig: &str, pos: i64, lo: i64, hi: i64) -> (Vec<(Vec<u8>, bool)>, [u32; 4]) {
    let (header, mut reader) = open_indexed(cram, Some(refp)).expect("open cram");
    let region: Region = format!("{contig}:{lo}-{hi}").parse().expect("region");
    let mut pile = [0u32; 4];
    let mut reads: Vec<(Vec<u8>, bool)> = Vec::new();
    for result in reader.query(&header, &region).expect("query") {
        let record = result.expect("record");
        let flags = record.flags();
        if flags.is_secondary() || flags.is_supplementary() || flags.is_duplicate() || flags.is_unmapped() {
            continue;
        }
        let Some(start) = record.alignment_start().map(|p| p.get() as i64) else {
            continue;
        };
        let seq = record.sequence();
        let seqb = seq.as_ref();
        let mut ref_pos = start;
        let mut qoff = 0usize;
        let mut win: Vec<u8> = Vec::new();
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
                    // keep insertion bases inside the window so an indel haplotype is preserved
                    if kind == noodles::sam::alignment::record::cigar::op::Kind::Insertion
                        && ref_pos > lo
                        && ref_pos <= hi
                    {
                        for i in 0..len {
                            if let Some(&b) = seqb.get(qoff + i) {
                                win.push(b.to_ascii_uppercase());
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
            reads.push((win, covers_pos));
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

/// Base the consensus carries at reference coordinate `pos`, by semiglobally aligning the consensus
/// (query) to the reference window `win_ref` (which starts at reference coordinate `win_start`).
fn consensus_base_at(consensus: &[u8], win_ref: &[u8], win_start: i64, pos: i64) -> char {
    let score = |a: u8, b: u8| if a == b { 1i32 } else { -4i32 };
    let mut aligner = PwAligner::new(-5, -1, score);
    let aln = aligner.semiglobal(consensus, win_ref);
    if std::env::var("PROBE_DEBUG").is_ok() {
        eprintln!(
            "DEBUG consensus.len={} win_ref.len={} xstart={} ystart={} xend={} yend={} score={}",
            consensus.len(), win_ref.len(), aln.xstart, aln.ystart, aln.xend, aln.yend, aln.score
        );
        eprintln!("  consensus raw[..20]: {:?}", &consensus[..consensus.len().min(20)]);
        eprintln!("  consensus str[..40]: {}", String::from_utf8_lossy(&consensus[..consensus.len().min(40)]));
        eprintln!("  win_ref   str[..40]: {}", String::from_utf8_lossy(&win_ref[..win_ref.len().min(40)]));
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
    println!(
        "{:>10} {:>3} {:>15}   pileup(alt)   -> reassembly (ref-supp/alt-supp/dropped)   verdict",
        "pos", "ref", "A/C/G/T"
    );
    for &pos in &positions {
        let lo = (pos - window).max(1);
        let hi = pos + window;
        let ref_base = refseq[(pos - 1) as usize] as char;
        let (reads, pile) = window_reads(cram, refp, contig, pos, lo, hi);
        let total: u32 = pile.iter().sum();
        // Candidate alt = the most common non-reference base at pos.
        let ref_i = base_index(ref_base as u8).unwrap_or(0);
        let alt_i = (0..4)
            .filter(|&i| i != ref_i)
            .max_by_key(|&i| pile[i])
            .unwrap_or(ref_i);
        let alt_base = charb(alt_i);
        let alt_frac = if total > 0 { (total - pile[ref_i]) as f64 / total as f64 } else { 0.0 };

        // Reference vs alternate haplotype over the window (alt = ref with the SNV at pos).
        let win_ref: Vec<u8> = refseq[(lo - 1) as usize..(hi as usize).min(refseq.len())].to_vec();
        let mut win_alt = win_ref.clone();
        let off = (pos - lo) as usize;
        if off < win_alt.len() {
            win_alt[off] = alt_base as u8;
        }

        // Realign each spanning read to both haplotypes; drop reads that align poorly to *both*
        // (they carry mismatches beyond the site → misaligned paralog / junk, not from this locus).
        let (mut ref_supp, mut alt_supp, mut dropped) = (0u32, 0u32, 0u32);
        let mut kept: Vec<Vec<u8>> = Vec::new();
        for (win, covers) in &reads {
            if !covers {
                continue; // can't distinguish ref from alt without the site
            }
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
        }
        // POA consensus of the surviving reads as a cross-check (informational).
        let cons_note = if kept.len() >= 2 {
            let cons = poa_consensus(&kept);
            let cb = consensus_base_at(&cons, &win_ref, lo, pos);
            format!("  [POA={cb}]")
        } else {
            String::new()
        };

        let verdict = if alt_supp > ref_supp && alt_supp >= 2 {
            format!("DERIVED {ref_base}>{alt_base}{cons_note}")
        } else if ref_supp >= alt_supp && ref_supp > 0 {
            format!("ancestral{cons_note}")
        } else {
            format!("no-call{cons_note}")
        };
        println!(
            "{pos:>10} {ref_base:>3} {:>15}   {alt_base}@{alt_frac:.2}       -> {ref_supp:>2}/{alt_supp:>2}/{dropped:<2}                       {verdict}",
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
