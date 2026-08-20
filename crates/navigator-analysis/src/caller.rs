//! A **haploid** variant caller built for this purpose (plan §4b). It replaces GATK for the Y
//! chromosome and the mtDNA.
//!
//! There is no pure-code Scala caller to port. The legacy app called out to GATK
//! `HaplotypeCaller --sample-ploidy 1`, which force-calls at the tree sites, and to `Mutect2
//! --mitochondria` or a haploid `HaplotypeCaller`, which discovers de-novo. It then subtracted
//! the known tree positions to get the private variants. This module does both of those modes by
//! a **consensus call over a pileup**. That method works here because the Y and the mtDNA are
//! haploid, at ploidy 1, so it needs no diploid local reassembly.
//!
//! There are two modes:
//! 1. [`force_call_sites`] gives a genotype at each known tree `Site`, from alleles that the
//!    caller already has. This is for haplogroup assignment. It makes a pileup, takes the
//!    consensus base, and reports whether that base is the ref allele or the alt allele of the
//!    site.
//! 2. [`call_denovo`] walks the contig and emits each position whose consensus base is different
//!    from the reference. Those are the candidate private variants. [`subtract_known`] then
//!    removes the known tree positions to give the private set.
//!
//! **v1 handles a SNP alone** (plan §4b). An indel or a homopolymer is where a simple pileup call
//! moves away from GATK, and light local realignment is the planned answer. So this module skips
//! an indel allele, and it treats such an allele as advisory until the parity harness of §4c
//! checks it. The defaults are start points, and the harness will tune them.
//!
//! On memory: the de-novo path walks the contig in chunks that overlap (`denovo_chunk`). The
//! chunk, and not the length of the contig, thereby bounds the dense tally at each position. Context overlaps on both sides, so a realignment window that crosses a chunk boundary
//! stays fully visible. The force-call path tallies the target sites alone, which is sparse and
//! costs little at any size.

use std::collections::{HashMap, HashSet};
use std::path::Path;

use noodles::core::Region;
use noodles::fasta;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::RecordBuf;

use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::error::AnalysisError;
use crate::genotype::{self, GenotypeResult};
use crate::reader;
use crate::realign;
use crate::reassembly;

/// The algorithm version of a de-novo caller artifact. Raise it after a change that alters the
/// output. The local realignment raised it to -2, and the local-reassembly resolver to -3.
pub const DENOVO_VERSION: &str = "haploid-denovo-3";

/// Algorithm version for site-genotype (panel) artifacts.
pub const GENOTYPE_VERSION: &str = "genotype-1";

/// The parameters of a haploid call. The defaults are the v1 start points, and §4c gates them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HaploidCallerParams {
    /// The smallest depth of reads that pass, which are the reads that clear the quality
    /// filters, before the caller makes any call.
    pub min_depth: u32,
    /// The caller drops a read below this MAPQ.
    pub min_mapping_quality: u8,
    /// Bases below this quality are not counted.
    pub min_base_quality: u8,
    /// The consensus base must hold at least this fraction of the depth that passes, before the
    /// caller makes a call.
    pub min_allele_fraction: f64,
    /// The allele-balance filter, which finds a paralog, at a haploid site. A true haploid site
    /// on the Y or the mtDNA carries almost one allele alone. A large *second* allele shows a
    /// paralog, or a read that the aligner put in the wrong place, which piles two loci together.
    /// The caller drops such a site.
    ///
    /// The filter fires only when the second most common allele meets two conditions. It has
    /// `min_paralog_minor_reads` reads or more, AND its fraction is above
    /// `max_minor_allele_fraction`. One read that disagrees, which is a sequencing error, does
    /// not fire it. Set the fraction to `1.0` or more to turn the filter off. See
    /// PangenomeExpansion.md, Phase 1.
    pub max_minor_allele_fraction: f64,
    /// Minimum second-allele read count for the paralog filter to engage (guards low depth).
    pub min_paralog_minor_reads: u32,
    /// Run light local realignment around candidate indels before de-novo calling.
    pub local_realign: bool,
    /// Minimum reads with indel evidence at a position to open a realignment window.
    pub realign_min_indel_reads: u32,
    /// Padding (bp) added around indel-evidence runs to form a realignment window.
    pub realign_pad: i64,
    /// The chunk size, in bp, of the de-novo pass. The caller walks the contig in chunks to hold
    /// the memory down. One chunk holds dense arrays for `chunk + 2*overlap` positions.
    pub denovo_chunk: usize,
    /// The context overlap, in bp, that the caller walks on each side of a chunk. A realignment
    /// window that crosses a chunk boundary then stays fully visible. It must be more than
    /// `realign_pad`.
    pub denovo_overlap: usize,
    /// Send a position that the paralog gate would drop to the local-reassembly resolver,
    /// [`crate::reassembly`], and do not throw it away. This recovers a haploid SNV where the
    /// reference alignment is wrong, which is Option B of the private-Y work. When this is off,
    /// the caller keeps its pileup-only behaviour.
    pub reassembly: bool,
    /// Half-width (bp) of the window extracted around a reassembly candidate. Must be ≤
    /// `denovo_overlap` so a boundary window stays inside the processed chunk.
    pub reassembly_window: i64,
}

impl Default for HaploidCallerParams {
    fn default() -> Self {
        HaploidCallerParams {
            min_depth: 4,
            min_mapping_quality: 20,
            min_base_quality: 20,
            min_allele_fraction: 0.5,
            max_minor_allele_fraction: 0.2,
            min_paralog_minor_reads: 2,
            local_realign: true,
            realign_min_indel_reads: 3,
            realign_pad: 15,
            denovo_chunk: 8_000_000,
            denovo_overlap: 500,
            reassembly: true,
            reassembly_window: 40,
        }
    }
}

/// A known tree/ancestry site to genotype (mirrors the Scala `Locus`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Site {
    pub name: String,
    pub contig: String,
    pub position: i64, // 1-based
    pub reference_allele: String,
    pub alternate_allele: String,
}

/// The allele called at a force-call site (haploid → one allele).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CalledAllele {
    Reference,
    Alternate,
    /// Insufficient depth, below-threshold consensus, or consensus is a third allele.
    NoCall,
}

/// Genotype at a known site.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct GenotypeCall {
    pub name: String,
    pub contig: String,
    pub position: i64,
    pub reference_allele: String,
    pub alternate_allele: String,
    pub called: CalledAllele,
    pub depth: u32, // passing depth (all bases)
    pub ref_depth: u32,
    pub alt_depth: u32,
    pub allele_fraction: f64, // alt_depth / depth
}

/// A diploid or haploid genotype at a known site, from the genotype-likelihood model. `dosage` is
/// the count of the alt allele, from 0 to the ploidy, or -1 for a no-call. That is the encoding
/// that the population, ancestry and IBD paths read.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SiteGenotype {
    pub name: String,
    pub contig: String,
    pub position: i64,
    pub reference_allele: String,
    pub alternate_allele: String,
    pub ploidy: u8,
    pub dosage: i32,
    pub gq: u8,
    pub depth: u32,
    pub ref_depth: u32,
    pub alt_depth: u32,
    pub pls: Vec<u8>,
    /// The VCF genotype string, such as `"1/2"`, for a site with more than two alleles. When it
    /// is `None`, the genotype comes from `dosage`, and the site has two alleles. This field is
    /// additive, so an old cached blob decodes to `None`.
    #[serde(default)]
    pub gt: Option<String>,
    /// The read depth of each allele, `[ref, alt1, alt2, …]`, at a site with more than two
    /// alleles. `None` means that the site has two alleles, and you read `ref_depth` and
    /// `alt_depth`.
    #[serde(default)]
    pub allele_depths: Option<Vec<u32>>,
}

/// A de-novo SNP call (consensus base differs from reference).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct VariantCall {
    pub contig: String,
    pub position: i64, // 1-based
    pub reference_allele: char,
    pub alternate_allele: char,
    pub depth: u32,     // passing depth
    pub alt_depth: u32, // reads supporting the consensus alt
    pub allele_fraction: f64,
    /// The confidence, on the Phred scale. The local-reassembly resolver
    /// ([`crate::reassembly`]) sets it on a call that it recovered. It is `None` on the plain
    /// pileup path and the gVCF path. This field is additive, so an old cached blob decodes to
    /// `None`.
    #[serde(default)]
    pub quality: Option<f64>,
}

const BASES: [u8; 4] = *b"ACGT";

pub(crate) fn base_index(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// argmax over the four base counts; ties keep the earlier base (A<C<G<T).
fn consensus(counts: &[u32; 4]) -> (usize, u32) {
    let mut bi = 0;
    let mut best = counts[0];
    for (i, &count) in counts.iter().enumerate().skip(1) {
        if count > best {
            best = count;
            bi = i;
        }
    }
    (bi, best)
}

/// The allele-balance filter, which finds a paralog, over a haploid pileup. A true haploid site
/// carries almost one allele alone. The site looks like it has two alleles when the second most
/// common allele has enough reads (`min_paralog_minor_reads`) and a large enough fraction (above
/// `max_minor_allele_fraction`). That is an artifact of a paralog, or of a read in the wrong
/// place. The caller must drop such a site. One read that disagrees, which is probably a
/// sequencing error, does not fire the filter.
fn is_paralogous(counts: &[u32; 4], depth: u32, params: &HaploidCallerParams) -> bool {
    if depth == 0 {
        return false;
    }
    let (bi, _) = consensus(counts);
    let second = counts
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != bi)
        .map(|(_, &v)| v)
        .max()
        .unwrap_or(0);
    second >= params.min_paralog_minor_reads && (second as f64 / depth as f64) > params.max_minor_allele_fraction
}

/// Read passes the de-novo/force-call filters (primary, not dup/qc-fail, MAPQ ok).
fn passes(record: &RecordBuf, params: &HaploidCallerParams) -> bool {
    let f = record.flags();
    if f.is_secondary() || f.is_supplementary() || f.is_duplicate() || f.is_qc_fail() {
        return false;
    }
    record.mapping_quality().map_or(255u8, |m| m.get()) >= params.min_mapping_quality
}

/// Resolve a contig's length from the BAM header.
pub(crate) fn contig_length(header: &noodles::sam::Header, contig: &str) -> Option<usize> {
    header
        .reference_sequences()
        .iter()
        .find(|(name, _)| {
            let n: &[u8] = name.as_ref();
            n == contig.as_bytes()
        })
        .map(|(_, map)| map.length().get())
}

/// Find the length of a contig. It reads the alignment header at `bam_path`. A CRAM file needs
/// `reference`.
pub(crate) fn read_contig_length(
    bam_path: &Path,
    contig: &str,
    reference: Option<&Path>,
) -> Result<usize, AnalysisError> {
    let header = reader::read_header(bam_path, reference)?;
    contig_length(&header, contig).ok_or_else(|| AnalysisError::Message(format!("contig {contig} not in BAM header")))
}

/// Load the full reference sequence of a contig, in one query against an indexed FASTA. The
/// chunks of the caller share it read-only. Each chunk takes a slice for its own window, and it
/// does not run the query again.
fn load_contig_sequence(reference_path: &Path, contig: &str, length: usize) -> Result<Vec<u8>, AnalysisError> {
    let mut fasta_reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(reference_path)
        .map_err(|e| AnalysisError::io(reference_path, e))?;
    let region: Region = format!("{contig}:1-{length}")
        .parse()
        .map_err(|_| AnalysisError::Message(format!("bad region for contig {contig}")))?;
    Ok(fasta_reader
        .query(&region)
        .map_err(|e| AnalysisError::io(reference_path, e))?
        .sequence()
        .as_ref()
        .to_vec())
}

/// The names of the contigs, which are the reference sequences, in the alignment header. A CRAM
/// file needs `reference`. The caller uses this to skip a lifted position that lands on a contig
/// that the alignment does not hold.
pub fn header_contig_names(bam_path: &Path, reference: Option<&Path>) -> Result<Vec<String>, AnalysisError> {
    let header = reader::read_header(bam_path, reference)?;
    Ok(header
        .reference_sequences()
        .keys()
        .map(|name| String::from_utf8_lossy(name.as_ref()).into_owned())
        .collect())
}

/// Contig name → length from the alignment header (for whole-genome walkers like SV).
pub fn header_contig_lengths(
    bam_path: &Path,
    reference: Option<&Path>,
) -> Result<std::collections::BTreeMap<String, i64>, AnalysisError> {
    let header = reader::read_header(bam_path, reference)?;
    Ok(header
        .reference_sequences()
        .iter()
        .map(|(name, seq)| {
            (
                String::from_utf8_lossy(name.as_ref()).into_owned(),
                seq.length().get() as i64,
            )
        })
        .collect())
}

/// Sparse A/C/G/T tally at the given 1-based target positions (force-call path), keyed
/// by 0-based position. Also returns the contig length.
fn tally_targets(
    bam_path: &Path,
    contig: &str,
    params: &HaploidCallerParams,
    targets: &HashSet<i64>,
    reference: Option<&Path>,
) -> Result<(usize, HashMap<usize, [u32; 4]>), AnalysisError> {
    let (header, mut reader) = reader::open_indexed(bam_path, reference)?;
    let length = contig_length(&header, contig)
        .ok_or_else(|| AnalysisError::Message(format!("contig {contig} not in BAM header")))?;

    let mut counts: HashMap<usize, [u32; 4]> = HashMap::new();
    let region: Region = contig
        .parse()
        .map_err(|_| AnalysisError::Message(format!("bad region for contig {contig}")))?;
    for result in reader.query(&header, &region)? {
        let record = result?;
        if !passes(&record, params) {
            continue;
        }
        let start = match record.alignment_start() {
            Some(p) => p.get(),
            None => continue,
        };
        let seq = record.sequence();
        let quals = record.quality_scores();
        let quals = quals.as_ref();
        let mut ref_pos = start;
        let mut query_off = 0usize;
        for op in record.cigar().as_ref() {
            let kind = op.kind();
            let len = op.len();
            match (kind.consumes_reference(), kind.consumes_read()) {
                (true, true) => {
                    for i in 0..len {
                        let pos = ref_pos + i;
                        if targets.contains(&(pos as i64)) {
                            let base_q = quals.get(query_off + i).copied().unwrap_or(0);
                            if base_q >= params.min_base_quality {
                                if let Some(bi) = seq.get(query_off + i).and_then(base_index) {
                                    counts.entry(pos - 1).or_insert([0; 4])[bi] += 1;
                                }
                            }
                        }
                    }
                    ref_pos += len;
                    query_off += len;
                }
                (true, false) => ref_pos += len,
                (false, true) => query_off += len,
                (false, false) => {}
            }
        }
    }
    Ok((length, counts))
}

/// Call the consensus base at each 1-based `target` position on `contig`. This is the haploid
/// genotype for a haplogroup assignment.
///
/// A position gets a call only when it clears `min_depth` reads that pass, and when the consensus
/// base holds at least `min_allele_fraction` of that depth. A position with no call is absent
/// from the result. Returns a map from a position to an uppercase base.
pub fn call_bases_at(
    bam_path: &Path,
    contig: &str,
    targets: &HashSet<i64>,
    params: &HaploidCallerParams,
    reference: Option<&Path>,
) -> Result<HashMap<i64, char>, AnalysisError> {
    let (_len, counts) = tally_targets(bam_path, contig, params, targets, reference)?;
    const BASES: [char; 4] = ['A', 'C', 'G', 'T'];
    let mut calls = HashMap::new();
    for (pos0, c) in counts {
        let depth: u32 = c.iter().sum();
        if depth < params.min_depth {
            continue;
        }
        let (bi, best) = consensus(&c);
        if (best as f64) < params.min_allele_fraction * depth as f64 {
            continue;
        }
        if is_paralogous(&c, depth, params) {
            continue; // bi-allelic at a haploid site — paralog/mismapping, drop the call
        }
        calls.insert((pos0 + 1) as i64, BASES[bi]);
    }
    Ok(calls)
}

/// A diagnostic. It gives the raw A/C/G/T tally of the reads that pass, at each 1-based
/// `target`. That is the evidence **behind** the consensus that [`call_bases_at`] takes, before
/// the depth filter, the allele-fraction filter and the paralog filter.
///
/// Returns a map from a 1-based position to `[A, C, G, T]` counts. A position is absent when no
/// read that passes covered it. Use it to log what the reads show at a tree SNP.
pub fn tally_at(
    bam_path: &Path,
    contig: &str,
    targets: &HashSet<i64>,
    params: &HaploidCallerParams,
    reference: Option<&Path>,
) -> Result<HashMap<i64, [u32; 4]>, AnalysisError> {
    let (_len, counts) = tally_targets(bam_path, contig, params, targets, reference)?;
    Ok(counts.into_iter().map(|(pos0, c)| ((pos0 + 1) as i64, c)).collect())
}

/// The indel allele that a tree locus expects, from its ancestral and derived alleles in the
/// left-anchored VCF form. It is one of two things. It is an insertion of the bases at the end,
/// where `A`→`ATT` gives Ins("TT"). Or it is a deletion of the difference in length, where
/// `TA`→`T` gives Del(1). It is `None` for a SNP, and for an allele that is complex or not
/// left-anchored.
fn expected_indel_allele(ancestral: &str, derived: &str) -> Option<IndelAllele> {
    let (a, d) = (ancestral.as_bytes(), derived.as_bytes());
    if d.len() > a.len() && d.starts_with(a) {
        Some(IndelAllele::Ins(d[a.len()..].to_ascii_uppercase()))
    } else if a.len() > d.len() && a.starts_with(d) {
        Some(IndelAllele::Del((a.len() - d.len()) as u32))
    } else {
        None
    }
}

/// Walk the CIGAR of one read from `start`, which is 1-based, and collect each indel event as
/// `(anchor 1-based, allele)`. The anchor of a deletion is the first deleted ref base. The anchor
/// of an insertion is the ref base that comes after the insertion. Returns the events, and the
/// inclusive reference end of the read.
fn read_indel_events(record: &RecordBuf, start: i64) -> (Vec<(i64, IndelAllele)>, i64) {
    let seq = record.sequence();
    let mut ref_pos = start;
    let mut query_off = 0usize;
    let mut events = Vec::new();
    for op in record.cigar().as_ref() {
        let (kind, len) = (op.kind(), op.len());
        match (kind.consumes_reference(), kind.consumes_read()) {
            (true, true) => {
                ref_pos += len as i64;
                query_off += len;
            }
            (true, false) => {
                events.push((ref_pos, IndelAllele::Del(len as u32)));
                ref_pos += len as i64;
            }
            (false, true) => {
                if kind == Kind::Insertion {
                    let s: Vec<u8> = (0..len)
                        .filter_map(|i| seq.get(query_off + i).map(|b| b.to_ascii_uppercase()))
                        .collect();
                    events.push((ref_pos, IndelAllele::Ins(s)));
                }
                query_off += len;
            }
            (false, false) => {}
        }
    }
    (events, ref_pos - 1)
}

/// Genotype the **indel** loci of a tree, at given targets. Each target is
/// `(pos, ancestral, derived)`, in the left-anchored VCF form, where `pos` is the anchor base.
/// `A`→`ATT` is an insertion, and `TA`→`T` is a deletion. At each locus the code examines the
/// reads that cover it. A read that carries the matching insertion or deletion, after
/// left-normalization into the reference repeat, supports the derived allele.
///
/// **This function only adds.** A locus with a clear derived majority over `min_depth` comes out
/// as [`haplo::INDEL_DERIVED`] at `pos`. Everything else stays a **no-call**, and it never becomes
/// an ancestral contradiction. That covers a locus with no indel support, a locus with low depth,
/// and reads that only *cover* the site cleanly.
///
/// The reason is noise. Take an indel genotype around a homopolymer or an STR. A read that
/// covers the site cleanly is often the alternate form that the aligner chose for the same indel.
/// To call those ancestral would contradict a thin node for no reason. A node at
/// d == 0 that takes one false ancestral fires the confident-divergence guard, and that vetoes
/// the whole lineage.
///
/// So an indel only ever *confirms* a branch. That matches the intent: cover the many DecodingUs
/// branches that an indel defines, when the sample carries them.
///
/// This needs a `reference`, to left-normalize and to know the deleted bases. It returns an empty
/// result without one, and also when the FASTA does not hold the contig.
pub fn call_indels_at(
    bam_path: &Path,
    contig: &str,
    targets: &[(i64, String, String)],
    params: &HaploidCallerParams,
    reference: Option<&Path>,
) -> Result<HashMap<i64, char>, AnalysisError> {
    let Some(reference) = reference else {
        return Ok(HashMap::new());
    };
    let Ok(refbytes) = reader::read_contig_sequence(reference, contig) else {
        return Ok(HashMap::new()); // contig naming mismatch with the FASTA — skip indels, keep SNPs
    };

    // Parse + left-normalize each target's expected allele. proc_lo = 1 (full-contig reference).
    struct PTarget {
        pos: i64,      // VCF POS (anchor), 1-based
        n_anchor: i64, // normalized CIGAR anchor (= pos+1 canonically)
        n_allele: IndelAllele,
        span_end: i64, // last ref base the ref-spanning read must cover
    }
    let mut ptargets: Vec<PTarget> = Vec::new();
    for (pos, anc, der) in targets {
        let Some(al) = expected_indel_allele(anc, der) else {
            continue;
        };
        let (n_anchor, n_allele) = left_normalize(pos + 1, &al, &refbytes, 1);
        let del_len = match &al {
            IndelAllele::Del(l) => *l as i64,
            IndelAllele::Ins(_) => 0,
        };
        ptargets.push(PTarget {
            pos: *pos,
            n_anchor,
            n_allele,
            span_end: pos + del_len.max(1),
        });
    }
    if ptargets.is_empty() {
        return Ok(HashMap::new());
    }
    ptargets.sort_by_key(|t| t.pos);
    let positions: Vec<i64> = ptargets.iter().map(|t| t.pos).collect();

    // One pass over the whole contig, as the SNP tally does. Walk every read once. At each
    // target that the read covers, add up two kinds of support. One is a match, where the read
    // carries the indel. The other is a clean cover, where the read covers the target and shows
    // the reference there.
    let (header, mut reader) = reader::open_indexed(bam_path, Some(reference))?;
    let region: Region = contig
        .parse()
        .map_err(|_| AnalysisError::Message(format!("bad region for contig {contig}")))?;
    let mut matched = vec![0u32; ptargets.len()];
    let mut refspan = vec![0u32; ptargets.len()];

    for result in reader.query(&header, &region)? {
        let record = result?;
        if !passes(&record, params) {
            continue;
        }
        let Some(start) = record.alignment_start().map(|p| p.get() as i64) else {
            continue;
        };
        let (raw, ref_end) = read_indel_events(&record, start);
        let events: Vec<(i64, IndelAllele)> = raw
            .into_iter()
            .map(|(a, al)| left_normalize(a, &al, &refbytes, 1))
            .collect();
        // Targets whose anchor this read could inform: pos in [start, ref_end].
        let lo = positions.partition_point(|&p| p < start);
        let hi = positions.partition_point(|&p| p <= ref_end);
        for i in lo..hi {
            let t = &ptargets[i];
            if events.iter().any(|(a, al)| *a == t.n_anchor && *al == t.n_allele) {
                matched[i] += 1;
            } else if start <= t.pos && ref_end >= t.span_end && !events.iter().any(|(a, _)| *a == t.n_anchor) {
                refspan[i] += 1;
            }
        }
    }

    let frac = params.min_allele_fraction;
    let mut out = HashMap::new();
    for (i, t) in ptargets.iter().enumerate() {
        let (m, r) = (matched[i], refspan[i]);
        let depth = m + r;
        if depth < params.min_depth {
            continue; // no-call
        }
        if m > r && m as f64 >= frac * depth as f64 {
            out.insert(t.pos, crate::haplo::INDEL_DERIVED);
        }
        // Not a clear derived majority → no-call (additive-only: never an ancestral contradiction).
    }
    Ok(out)
}

/// A dense A/C/G/T tally, and the indel evidence at each position, over the 1-based inclusive
/// region `[lo, hi]`. The index is `pos - lo`. The de-novo path, which works in chunks, uses
/// this.
pub(crate) fn tally_region(
    bam_path: &Path,
    contig: &str,
    params: &HaploidCallerParams,
    lo: usize,
    hi: usize,
    reference: Option<&Path>,
) -> Result<(Vec<[u32; 4]>, Vec<u32>), AnalysisError> {
    let n = hi - lo + 1;
    let mut counts = vec![[0u32; 4]; n];
    let mut indel = vec![0u32; n];

    let (header, mut reader) = reader::open_indexed(bam_path, reference)?;
    let region: Region = format!("{contig}:{lo}-{hi}")
        .parse()
        .map_err(|_| AnalysisError::Message(format!("bad region for {contig}")))?;

    for result in reader.query(&header, &region)? {
        let record = result?;
        if !passes(&record, params) {
            continue;
        }
        let start = match record.alignment_start() {
            Some(p) => p.get(),
            None => continue,
        };
        let seq = record.sequence();
        let quals = record.quality_scores();
        let quals = quals.as_ref();
        let mut ref_pos = start;
        let mut query_off = 0usize;
        for op in record.cigar().as_ref() {
            let kind = op.kind();
            let len = op.len();
            match (kind.consumes_reference(), kind.consumes_read()) {
                (true, true) => {
                    for i in 0..len {
                        let pos = ref_pos + i;
                        if pos >= lo && pos <= hi {
                            let base_q = quals.get(query_off + i).copied().unwrap_or(0);
                            if base_q >= params.min_base_quality {
                                if let Some(bi) = seq.get(query_off + i).and_then(base_index) {
                                    counts[pos - lo][bi] += 1;
                                }
                            }
                        }
                    }
                    ref_pos += len;
                    query_off += len;
                }
                (true, false) => {
                    for k in 0..len {
                        let pos = ref_pos + k;
                        if pos >= lo && pos <= hi {
                            indel[pos - lo] += 1;
                        }
                    }
                    ref_pos += len;
                }
                (false, true) => {
                    if kind == Kind::Insertion && ref_pos >= lo && ref_pos <= hi {
                        indel[ref_pos - lo] += 1;
                    }
                    query_off += len;
                }
                (false, false) => {}
            }
        }
    }
    Ok((counts, indel))
}

/// Record the `(base, qual)` of one read, where it passes, at each of the `targets` that the read
/// covers. `targets` comes in sorted order, and its positions are 1-based. The function does this
/// in a **single** CIGAR walk from the start of the alignment.
///
/// The walk skips a target in a deletion, a ref skip, an insertion or a clip. It also skips a
/// target that is past the read, one below `min_base_quality`, and one whose base is not ACGT.
///
/// This is the many-target form of a probe at one site. One walk feeds many sites. Some nearby
/// panel sites share a long read, and the code decodes and walks that read once, not once for
/// each site.
fn collect_bases(record: &RecordBuf, targets: &[i64], min_base_quality: u8, obs: &mut HashMap<i64, Vec<(u8, u8)>>) {
    let Some(start) = record.alignment_start() else { return };
    let start = start.get() as i64;
    let seq = record.sequence();
    let quals = record.quality_scores();
    let quals = quals.as_ref();

    // First target at/after the read's start; advance through the window as the CIGAR consumes ref.
    let mut ti = targets.partition_point(|&t| t < start);
    let mut ref_pos = start;
    let mut query_off = 0usize;
    for op in record.cigar().as_ref() {
        if ti >= targets.len() {
            break;
        }
        let (cr, cq) = (op.kind().consumes_reference(), op.kind().consumes_read());
        let len = op.len() as i64;
        if cr && cq {
            let end = ref_pos + len; // exclusive
            while ti < targets.len() && targets[ti] < end {
                let t = targets[ti];
                let off = query_off + (t - ref_pos) as usize;
                let base_q = quals.get(off).copied().unwrap_or(0);
                if base_q >= min_base_quality {
                    if let Some(base) = seq.get(off) {
                        if base_index(base).is_some() {
                            obs.entry(t).or_default().push((base, base_q));
                        }
                    }
                }
                ti += 1;
            }
            ref_pos = end;
            query_off += len as usize;
        } else if cr {
            // A deletion or a ref skip. A target inside the gap carries no base.
            let end = ref_pos + len;
            while ti < targets.len() && targets[ti] < end {
                ti += 1;
            }
            ref_pos = end;
        } else if cq {
            query_off += len as usize;
        }
    }
}

/// The `(base, qual)` observations that pass at each target site, keyed by a 1-based position.
/// Those are the ACGT bases that clear the quality filters. This is the input that the
/// genotype-likelihood model needs.
///
/// The code puts the targets into runs that touch each other. It splits a run only where the gap
/// between two adjacent sites is more than one read length. It then gets each run with a
/// **single** streaming query against the index.
///
/// The code then seeks straight to the regions that hold targets, and it never scans the whole
/// contig. It also decodes each read once. A point query at each site would fetch and convert
/// the long HiFi reads again, and one such read covers some nearby sites.
///
/// Inside a run, [`collect_bases`] gives the bases of each read to every target that the read
/// covers, in one CIGAR walk.
fn tally_site_observations(
    bam_path: &Path,
    contig: &str,
    params: &HaploidCallerParams,
    targets: &HashSet<i64>,
    reference: Option<&Path>,
) -> Result<HashMap<i64, Vec<(u8, u8)>>, AnalysisError> {
    let mut positions: Vec<i64> = targets.iter().copied().filter(|&p| p >= 1).collect();
    positions.sort_unstable();
    if positions.is_empty() {
        return Ok(HashMap::new());
    }

    // Split into runs whose consecutive sites lie within MAX_GAP. Past one read length, no read
    // can cover the gap. A split there costs nothing, because no shared read goes away, and it
    // skips the spans that hold no read.
    const MAX_GAP: i64 = 50_000;

    let (header, mut reader) = reader::open_indexed(bam_path, reference)?;
    let mut obs: HashMap<i64, Vec<(u8, u8)>> = HashMap::with_capacity(positions.len());

    let mut i = 0;
    while i < positions.len() {
        let mut j = i + 1;
        while j < positions.len() && positions[j] - positions[j - 1] <= MAX_GAP {
            j += 1;
        }
        let (lo, hi) = (positions[i], positions[j - 1]);
        let run = &positions[i..j];
        let region: Region = format!("{contig}:{lo}-{hi}")
            .parse()
            .map_err(|_| AnalysisError::Message(format!("bad region for {contig}:{lo}-{hi}")))?;
        for result in reader.query(&header, &region)? {
            let record = result?;
            if !passes(&record, params) {
                continue;
            }
            collect_bases(&record, run, params.min_base_quality, &mut obs);
        }
        i = j;
    }
    Ok(obs)
}

/// Genotype the known SNP sites on `contig`, at the given `ploidy`, with the
/// genotype-likelihood model. A ploidy of 1 is a haploid Y, MT or male X. A ploidy of 2 is an
/// autosome or a female X. This is the panel-genotype path that the population, ancestry and IBD
/// analyses read. The code skips a site that is not a SNP.
pub fn genotype_sites(
    bam_path: &Path,
    contig: &str,
    sites: &[Site],
    ploidy: u8,
    params: &HaploidCallerParams,
    reference: Option<&Path>,
) -> Result<Vec<SiteGenotype>, AnalysisError> {
    let targets: HashSet<i64> = sites
        .iter()
        .filter(|s| s.contig == contig && s.reference_allele.len() == 1 && s.alternate_allele.len() == 1)
        .map(|s| s.position)
        .collect();
    if targets.is_empty() {
        return Ok(Vec::new());
    }
    let obs = tally_site_observations(bam_path, contig, params, &targets, reference)?;

    let empty: Vec<(u8, u8)> = Vec::new();
    let mut out = Vec::new();
    for site in sites.iter().filter(|s| s.contig == contig) {
        if site.reference_allele.len() != 1 || site.alternate_allele.len() != 1 {
            continue; // SNP-only
        }
        let site_obs = obs.get(&site.position).unwrap_or(&empty);
        let GenotypeResult {
            dosage,
            pls,
            gq,
            depth,
            ref_depth,
            alt_depth,
        } = genotype::call_genotype(
            site_obs,
            site.reference_allele.as_bytes()[0],
            site.alternate_allele.as_bytes()[0],
            ploidy,
            params.min_depth,
        );
        out.push(SiteGenotype {
            name: site.name.clone(),
            contig: site.contig.clone(),
            position: site.position,
            reference_allele: site.reference_allele.clone(),
            alternate_allele: site.alternate_allele.clone(),
            ploidy,
            dosage,
            gq,
            depth,
            ref_depth,
            alt_depth,
            pls,
            gt: None,
            allele_depths: None,
        });
    }
    Ok(out)
}

/// Genotype `sites` across **every** contig that they cover, with one rayon task for each contig.
/// This is the entry point for a whole-genome panel. [`genotype_sites`] works on one contig, and
/// it is independent and IO-bound on its own index region, so the contigs parallelize cleanly.
///
/// The code joins the results together. The order across the contigs is not defined, because a
/// later step keys on the site and not on the order.
pub fn genotype_sites_all_contigs(
    bam_path: &Path,
    sites: &[Site],
    ploidy: u8,
    params: &HaploidCallerParams,
    reference: Option<&Path>,
    cancel: &crate::cancel::CancelToken,
) -> Result<Vec<SiteGenotype>, AnalysisError> {
    let contigs: Vec<&str> = sites
        .iter()
        .map(|s| s.contig.as_str())
        .collect::<std::collections::BTreeSet<&str>>()
        .into_iter()
        .collect();
    // Run on a pool that is safe for a decode, and not on the global rayon pool, whose stacks are
    // 2 MiB. Each task decodes CRAM records. Those recurse deeply on CRAM 3.1, and the stack would
    // overflow and abort. See [`reader::decode_pool`].
    let pool = crate::reader::decode_pool(contigs.len().max(1).min(crate::unified::analysis_thread_count()))?;
    let per_contig: Result<Vec<Vec<SiteGenotype>>, AnalysisError> = pool.install(|| {
        contigs
            .into_par_iter()
            .map(|contig| {
                cancel.check()?;
                genotype_sites(bam_path, contig, sites, ploidy, params, reference)
            })
            .collect()
    });
    Ok(per_contig?.into_iter().flatten().collect())
}

/// Reconcile the force-call genotypes of each alignment, over a shared set of sites, into one
/// **consensus** diploid genotype at each site. That is the joint genotype of the subject, across
/// the WGS runs of that person.
///
/// Each input holds the [`SiteGenotype`] values of one alignment, at the *union* of the variant
/// sites. All of them are on the same reference build, so `(contig, position, ref, alt)` line
/// up.
///
/// At each site the code holds a vote over the dosage classes {0,1,2}, weighted by depth. An
/// alignment whose depth is below `min_depth` casts its no-call, and the vote leaves it out. A
/// site that one run shows as hom-ref, and that is absent from the list, is then a real vote.
/// That resolves the case where run A is het and run B is hom-ref. A run that truly did not cover
/// the site abstains.
///
/// The result holds the **variant** consensus sites alone, which are het and hom-alt. A hom-ref
/// or no-call consensus is not a variant. The code sums the depth and the AD, and it takes the
/// maximum GQ over the alignments that support the call. It drops the PLs, because the
/// likelihoods of the separate runs do not combine into one PL here.
pub fn reconcile_site_genotypes(per_alignment: &[Vec<SiteGenotype>], min_depth: u32) -> Vec<SiteGenotype> {
    use std::collections::BTreeMap;
    struct Acc {
        repr: SiteGenotype,
        w: [f64; 3],
        counts: [usize; 3],
        depth: u64,
        ref_d: u64,
        alt_d: u64,
        gq: u8,
    }
    let mut groups: BTreeMap<(String, i64, String), Acc> = BTreeMap::new();
    for aln in per_alignment {
        for g in aln {
            let key = (g.contig.clone(), g.position, g.alternate_allele.clone());
            let acc = groups.entry(key).or_insert_with(|| Acc {
                repr: g.clone(),
                w: [0.0; 3],
                counts: [0; 3],
                depth: 0,
                ref_d: 0,
                alt_d: 0,
                gq: 0,
            });
            if g.depth < min_depth {
                continue; // under-covered in this run → abstain (not a hom-ref vote)
            }
            let d = g.dosage;
            if (0..=2).contains(&d) {
                // The weight of the depth bonus. It has the same shape as the WGS term of
                // consensus::obs_weight. The constant method factor cancels in the argmax.
                let weight = 1.0 + ((g.depth as f64).sqrt() / 10.0).min(1.0);
                w_add(&mut acc.w, &mut acc.counts, d as usize, weight);
                acc.depth += g.depth as u64;
                acc.ref_d += g.ref_depth as u64;
                acc.alt_d += g.alt_depth as u64;
                acc.gq = acc.gq.max(g.gq);
            }
        }
    }
    let mut out = Vec::new();
    for (_, acc) in groups {
        // Take the argmax of the weight. A tie goes to the higher count of raw runs that support
        // the call, and then to the lower dosage.
        let mut best = 0usize;
        for d in 1..3 {
            if acc.w[d] > acc.w[best] || (acc.w[d] == acc.w[best] && acc.counts[d] > acc.counts[best]) {
                best = d;
            }
        }
        let total: usize = acc.counts.iter().sum();
        if total == 0 || best == 0 {
            continue; // no-call or hom-ref consensus → not a variant
        }
        let mut g = acc.repr;
        g.name = String::new();
        g.dosage = best as i32;
        g.depth = acc.depth.min(u32::MAX as u64) as u32;
        g.ref_depth = acc.ref_d.min(u32::MAX as u64) as u32;
        g.alt_depth = acc.alt_d.min(u32::MAX as u64) as u32;
        g.gq = acc.gq;
        g.pls = Vec::new();
        g.gt = None;
        g.allele_depths = None;
        out.push(g);
    }
    out
}

#[inline]
fn w_add(w: &mut [f64; 3], counts: &mut [usize; 3], d: usize, weight: f64) {
    w[d] += weight;
    counts[d] += 1;
}

/// Force-call at the known SNP sites on `contig`. The caller already has the alleles. The code
/// skips a site that is not a SNP, because v1 handles a SNP alone. Such a site holds more than
/// one base in its ref allele or in its alt allele.
pub fn force_call_sites(
    bam_path: &Path,
    contig: &str,
    sites: &[Site],
    params: &HaploidCallerParams,
    reference: Option<&Path>,
) -> Result<Vec<GenotypeCall>, AnalysisError> {
    let targets: HashSet<i64> = sites
        .iter()
        .filter(|s| s.contig == contig)
        .map(|s| s.position)
        .collect();
    if targets.is_empty() {
        return Ok(Vec::new());
    }
    let (length, counts) = tally_targets(bam_path, contig, params, &targets, reference)?;

    let mut out = Vec::new();
    for site in sites.iter().filter(|s| s.contig == contig) {
        if site.reference_allele.len() != 1 || site.alternate_allele.len() != 1 {
            continue; // SNP-only
        }
        if site.position < 1 || (site.position as usize) > length {
            continue; // off-contig
        }
        let idx = (site.position - 1) as usize;
        let c = counts.get(&idx).copied().unwrap_or([0; 4]);
        let depth: u32 = c.iter().sum();
        let ref_bi = base_index(site.reference_allele.as_bytes()[0]);
        let alt_bi = base_index(site.alternate_allele.as_bytes()[0]);
        let ref_depth = ref_bi.map_or(0, |i| c[i]);
        let alt_depth = alt_bi.map_or(0, |i| c[i]);

        let (top_bi, top_count) = consensus(&c);
        let called = if depth < params.min_depth
            || top_count == 0
            || (top_count as f64 / depth as f64) < params.min_allele_fraction
            || is_paralogous(&c, depth, params)
        {
            CalledAllele::NoCall // includes the paralog/mismapping (bi-allelic) drop
        } else if Some(top_bi) == alt_bi {
            CalledAllele::Alternate
        } else if Some(top_bi) == ref_bi {
            CalledAllele::Reference
        } else {
            CalledAllele::NoCall // consensus is a third allele
        };

        out.push(GenotypeCall {
            name: site.name.clone(),
            contig: site.contig.clone(),
            position: site.position,
            reference_allele: site.reference_allele.clone(),
            alternate_allele: site.alternate_allele.clone(),
            called,
            depth,
            ref_depth,
            alt_depth,
            allele_fraction: if depth == 0 {
                0.0
            } else {
                alt_depth as f64 / depth as f64
            },
        });
    }
    Ok(out)
}

/// De-novo SNP discovery across `contig`. The code walks it in chunks that overlap, so the chunk
/// bounds the memory, and not the length of the contig. It emits each position whose consensus
/// base passes the depth filter and the fraction filter, and differs from the reference. The
/// context overlaps on both sides, so a realignment window that crosses a chunk boundary stays
/// fully visible.
pub fn call_denovo(
    bam_path: &Path,
    reference_path: &Path,
    contig: &str,
    params: &HaploidCallerParams,
    cancel: &crate::cancel::CancelToken,
) -> Result<Vec<VariantCall>, AnalysisError> {
    let length = read_contig_length(bam_path, contig, Some(reference_path))?;

    // Load the reference of the contig once. The chunks share it read-only, and each chunk takes
    // a slice for its own window. No chunk queries the FASTA again.
    let ref_seq = load_contig_sequence(reference_path, contig, length)?;

    // The emit ranges do not overlap, and they come in order. The code walks each one on its
    // own, with its own region query against the BAM index.
    //
    // The chunks stay large: `denovo_chunk` defaults to 8 MB. A CRAM container covers some MB.
    // Take a chunk smaller than a container. Every chunk over that container decodes it again.
    // The rayon pool limits how many chunks run at one time to its own size, so the peak memory
    // has a bound.
    let threads = crate::unified::analysis_thread_count();
    let chunk = params.denovo_chunk.max(1);
    let mut ranges: Vec<(usize, usize)> = Vec::new();
    let mut emit_lo = 1usize;
    while emit_lo <= length {
        let emit_hi = (emit_lo + chunk - 1).min(length);
        ranges.push((emit_lo, emit_hi));
        emit_lo = emit_hi + 1;
    }

    // Decode-safe worker stack: these tasks decode CRAM records, which recurse deeply on CRAM 3.1
    // (an overflow aborts the process). See [`reader::decode_pool`].
    let pool = crate::reader::decode_pool(threads)?;
    let nested: Vec<Vec<VariantCall>> = pool.install(|| {
        ranges
            .par_iter()
            .map(|&(lo, hi)| {
                // One check at each chunk. The work in a chunk has a bound, at a default of 8 MB
                // of reference. This bounds the delay between a click and a stop, and the chunks
                // do not have to get smaller.
                cancel.check()?;
                denovo_chunk(bam_path, reference_path, contig, params, &ref_seq, length, lo, hi)
            })
            .collect::<Result<Vec<_>, AnalysisError>>()
    })?;
    // The ranges do not overlap, and the code collects them in order. A flat join then keeps the
    // global position order.
    Ok(nested.into_iter().flatten().collect())
}

/// The de-novo SNP calls inside `[region_lo, region_hi]` alone, which is 1-based and inclusive,
/// on `contig`. It is the same tally, realign and reassembly path as [`call_denovo`], over one
/// bounded emit range. It loads the reference of the whole contig once, which costs little, and
/// it queries the region alone.
///
/// Use it in a debug tool or a check tool, for example to look at the recovery at given
/// positions. It does not walk the whole contig.
pub fn call_denovo_region(
    bam_path: &Path,
    reference_path: &Path,
    contig: &str,
    region_lo: usize,
    region_hi: usize,
    params: &HaploidCallerParams,
) -> Result<Vec<VariantCall>, AnalysisError> {
    let length = read_contig_length(bam_path, contig, Some(reference_path))?;
    let ref_seq = load_contig_sequence(reference_path, contig, length)?;
    let lo = region_lo.max(1);
    let hi = region_hi.min(length);
    if lo > hi {
        return Ok(Vec::new());
    }
    denovo_chunk(bam_path, reference_path, contig, params, &ref_seq, length, lo, hi)
}

/// The de-novo SNP calls for one emit range `[emit_lo, emit_hi]`, which is 1-based and inclusive.
/// It tallies a window with `denovo_overlap` of padding, so that a realignment window across the
/// boundary stays fully visible. It emits `[emit_lo, emit_hi]` alone. `ref_seq` is the reference
/// of the whole contig, where index 0 is position 1. Each call opens its own BAM reader, so it is
/// independent and safe across threads.
#[allow(clippy::too_many_arguments)]
fn denovo_chunk(
    bam_path: &Path,
    reference_path: &Path,
    contig: &str,
    params: &HaploidCallerParams,
    ref_seq: &[u8],
    length: usize,
    emit_lo: usize,
    emit_hi: usize,
) -> Result<Vec<VariantCall>, AnalysisError> {
    let overlap = params.denovo_overlap;
    let proc_lo = emit_lo.saturating_sub(overlap).max(1);
    let proc_hi = (emit_hi + overlap).min(length);
    // The reference window [proc_lo, proc_hi], with its index relative to proc_lo. The code
    // clamps it to what the FASTA returned, so a short contig tail reads as 'N', as before.
    let ref_chunk = &ref_seq[(proc_lo - 1).min(ref_seq.len())..proc_hi.min(ref_seq.len())];

    let (mut counts, indel) = tally_region(bam_path, contig, params, proc_lo, proc_hi, Some(reference_path))?;
    if params.local_realign {
        realign_region(
            bam_path,
            contig,
            ref_chunk,
            proc_lo,
            &mut counts,
            &indel,
            params,
            Some(reference_path),
        )?;
    }

    let mut out = Vec::new();
    // Stage A. Take the positions that the paralog gate would drop, but that carry a real
    // non-reference allele. Send those to the local-reassembly resolver, and do not throw them
    // away. That is Option B of the private-Y work.
    let mut active: Vec<reassembly::Candidate> = Vec::new();
    for pos in emit_lo..=emit_hi {
        let r = pos - proc_lo; // index into the chunk arrays
        let c = counts[r];
        let depth: u32 = c.iter().sum();
        if depth < params.min_depth {
            continue;
        }
        let (top_bi, top_count) = consensus(&c);
        if top_count == 0 {
            continue;
        }
        let ref_base = ref_chunk.get(r).copied().unwrap_or(b'N');
        if is_paralogous(&c, depth, params) {
            // The site has two alleles, and it is haploid. The pileup can not separate a true
            // derived SNV from a paralog artifact. Give it to reassembly, in stages B to E
            // below, and do not drop it.
            if params.reassembly {
                if let Some(cand) = active_candidate(pos as i64, &c, ref_base, params) {
                    active.push(cand);
                }
            }
            continue;
        }
        let frac = top_count as f64 / depth as f64;
        if frac < params.min_allele_fraction {
            continue;
        }
        if base_index(ref_base) == Some(top_bi) || base_index(ref_base).is_none() {
            continue; // matches reference, or reference is N/ambiguous
        }
        out.push(VariantCall {
            contig: contig.to_string(),
            position: pos as i64,
            reference_allele: ref_base.to_ascii_uppercase() as char,
            alternate_allele: BASES[top_bi] as char,
            depth,
            alt_depth: top_count,
            allele_fraction: frac,
            quality: None,
        });
    }

    // Stages B to F. Reassemble the windows that stage A sent on, and add the DERIVED calls that
    // they recover.
    if params.reassembly && !active.is_empty() {
        let recovered = resolve_active(bam_path, contig, ref_seq, length, &active, params, Some(reference_path))?;
        out.extend(recovered);
        out.sort_by_key(|v| v.position); // keep the chunk's calls in ascending position order
    }
    Ok(out)
}

/// A reassembly candidate at a position that the paralog gate held. It is the most common
/// **non-reference** base. The code keeps it only when it carries `min_paralog_minor_reads` reads
/// or more, which makes it a real alternate and not one error.
fn active_candidate(
    pos: i64,
    counts: &[u32; 4],
    ref_base: u8,
    params: &HaploidCallerParams,
) -> Option<reassembly::Candidate> {
    let ref_bi = base_index(ref_base)?;
    let (alt_bi, &alt_count) = counts
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != ref_bi)
        .max_by_key(|(_, &v)| v)?;
    if alt_count < params.min_paralog_minor_reads {
        return None;
    }
    Some(reassembly::Candidate {
        position: pos,
        ref_base: ref_base.to_ascii_uppercase(),
        alt_base: BASES[alt_bi],
    })
}

/// Stages B to F, for the candidates that a chunk sent on. It puts them into windows. It takes
/// the reads that cover each window, and projects them onto the reference frame of that window.
/// It genotypes them with the reassembly resolver. It returns the DERIVED recoveries as
/// `VariantCall` values, and a paralog artifact stays dropped.
fn resolve_active(
    bam_path: &Path,
    contig: &str,
    ref_seq: &[u8],
    length: usize,
    active: &[reassembly::Candidate],
    params: &HaploidCallerParams,
    reference: Option<&Path>,
) -> Result<Vec<VariantCall>, AnalysisError> {
    let (header, mut reader) = reader::open_indexed(bam_path, reference)?;
    let rparams = reassembly::ReassemblyParams {
        min_mapping_quality: params.min_mapping_quality,
        // Debug hooks for A/B-ing the v2 levers on real data (unset → the ReassemblyParams default).
        assemble_alt: std::env::var("NAVIGATOR_REASSEMBLY_ASSEMBLE")
            .map(|v| v != "0")
            .unwrap_or(reassembly::ReassemblyParams::default().assemble_alt),
        min_read_loglik: std::env::var("NAVIGATOR_REASSEMBLY_FLOOR")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(reassembly::ReassemblyParams::default().min_read_loglik),
        ..reassembly::ReassemblyParams::default()
    };
    let w = params.reassembly_window.max(1);

    // Merge the candidates whose windows overlap into one extraction window. The candidates come
    // from the emit loop in position order, from the lowest up.
    let mut out = Vec::new();
    let mut i = 0;
    while i < active.len() {
        let mut j = i + 1;
        while j < active.len() && active[j].position - active[j - 1].position <= 2 * w {
            j += 1;
        }
        let group = &active[i..j];
        let win_lo = (group[0].position - w).max(1);
        let win_hi = (group[group.len() - 1].position + w).min(length as i64);

        let ref_window: Vec<u8> = ref_seq[(win_lo - 1) as usize..win_hi as usize]
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect();
        let reads = extract_window_reads(&header, &mut reader, contig, win_lo, win_hi, group, params)?;
        for call in reassembly::genotype_window(&ref_window, win_lo, group, &reads, &rparams) {
            if call.genotype == reassembly::Zygosity::Derived {
                out.push(VariantCall {
                    contig: contig.to_string(),
                    position: call.position,
                    reference_allele: call.ref_base as char,
                    alternate_allele: call.alt_base as char,
                    depth: call.depth,
                    alt_depth: call.alt_depth,
                    allele_fraction: call.allele_fraction,
                    quality: Some(call.quality),
                });
            }
        }
        i = j;
    }
    Ok(out)
}

/// Project every read that covers `[win_lo, win_hi]` onto the reference frame of that window. The
/// result holds the sequence in the window frame, with the quality of each base, which the
/// PairHMM needs. It also holds the base and the quality that each read carries at each
/// candidate, which the depth and the dedup need.
///
/// There is one CIGAR walk for each read. An insertion inside the window stays, so an indel
/// haplotype survives. That matches what the realignment in the pileup intends.
fn extract_window_reads(
    header: &noodles::sam::Header,
    reader: &mut reader::IdxReader,
    contig: &str,
    win_lo: i64,
    win_hi: i64,
    candidates: &[reassembly::Candidate],
    params: &HaploidCallerParams,
) -> Result<Vec<reassembly::WindowRead>, AnalysisError> {
    let region: Region = format!("{contig}:{win_lo}-{win_hi}")
        .parse()
        .map_err(|_| AnalysisError::Message(format!("bad region for {contig}")))?;
    let mut reads = Vec::new();
    for result in reader.query(header, &region)? {
        let record = result?;
        if !passes(&record, params) {
            continue;
        }
        let Some(start) = record.alignment_start().map(|p| p.get() as i64) else {
            continue;
        };
        let mapq = record.mapping_quality().map_or(255u8, |m| m.get());
        let name = record.name().map(|n| n.to_vec()).unwrap_or_default();
        let seq = record.sequence();
        let quals = record.quality_scores();
        let qualb = quals.as_ref();

        let mut ref_pos = start;
        let mut qoff = 0usize;
        let mut wseq: Vec<u8> = Vec::new();
        let mut wq: Vec<u8> = Vec::new();
        let mut site_obs: Vec<Option<reassembly::SiteObs>> = vec![None; candidates.len()];
        for op in record.cigar().as_ref() {
            let kind = op.kind();
            let len = op.len();
            match (kind.consumes_reference(), kind.consumes_read()) {
                (true, true) => {
                    for i in 0..len {
                        let rp = ref_pos + i as i64;
                        if rp >= win_lo && rp <= win_hi {
                            if let Some(b) = seq.get(qoff + i) {
                                let q = qualb.get(qoff + i).copied().unwrap_or(0);
                                wseq.push(b.to_ascii_uppercase());
                                wq.push(q);
                                if q >= params.min_base_quality {
                                    if let Some(ci) = candidates.iter().position(|c| c.position == rp) {
                                        site_obs[ci] = Some(reassembly::SiteObs {
                                            base: b.to_ascii_uppercase(),
                                            qual: q,
                                        });
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
                    if kind == Kind::Insertion && ref_pos > win_lo && ref_pos <= win_hi {
                        for i in 0..len {
                            if let Some(b) = seq.get(qoff + i) {
                                wseq.push(b.to_ascii_uppercase());
                                wq.push(qualb.get(qoff + i).copied().unwrap_or(0));
                            }
                        }
                    }
                    qoff += len;
                }
                (false, false) => {}
            }
        }
        if wseq.len() >= 30 {
            reads.push(reassembly::WindowRead {
                name,
                seq: wseq,
                quals: wq,
                mapq,
                site_obs,
            });
        }
    }
    Ok(reads)
}

/// The nominal base quality of the de-novo **diploid** genotype likelihood.
///
/// The pileup that works in chunks keeps the A/C/G/T counts alone. A quality for each base would
/// use too much memory on a WGS run. Every base that it counted already cleared
/// `min_base_quality`. So the code evaluates the GL at this one representative phred value.
///
/// The genotype that comes out, 0/1 or 1/1 or 0/0, is robust to the exact value. The PL and the
/// GQ are approximate. The [`genotype_sites`] path, which works at one site, keeps the true
/// quality of each read where the exact likelihood matters.
const DENOVO_DIPLOID_Q: u8 = 30;
/// The count of reads that must support the alt allele before a site becomes a candidate variant
/// at all. It holds back a "het" that one sequencing error produced.
const DENOVO_MIN_ALT_READS: u32 = 2;

/// **De-novo diploid** SNV calling over a whole contig. It uses the same parallel pileup in
/// chunks as [`call_denovo`]. But it genotypes each variant site at ploidy 2, with the
/// genotype-likelihood model ([`genotype::call_genotype`]). So it emits heterozygous (0/1) and
/// homozygous-alt (1/1) calls, and not a haploid consensus alone.
///
/// v1 handles two alleles: REF, and the most common non-REF base. It does not call an indel here.
/// The output comes in position order, from the lowest up, as [`SiteGenotype`] at ploidy 2. Give
/// it to [`crate::vcf::write_diploid_vcf`].
pub fn call_denovo_diploid(
    bam_path: &Path,
    reference_path: &Path,
    contig: &str,
    params: &HaploidCallerParams,
    cancel: &crate::cancel::CancelToken,
) -> Result<Vec<SiteGenotype>, AnalysisError> {
    let length = read_contig_length(bam_path, contig, Some(reference_path))?;
    let ref_seq = load_contig_sequence(reference_path, contig, length)?;

    let threads = crate::unified::analysis_thread_count();
    let chunk = params.denovo_chunk.max(1);
    let mut ranges: Vec<(usize, usize)> = Vec::new();
    let mut emit_lo = 1usize;
    while emit_lo <= length {
        let emit_hi = (emit_lo + chunk - 1).min(length);
        ranges.push((emit_lo, emit_hi));
        emit_lo = emit_hi + 1;
    }

    // A worker stack that is safe for a decode. A CRAM 3.1 decode recurses deeply. See
    // [`reader::decode_pool`].
    let pool = crate::reader::decode_pool(threads)?;
    let nested: Vec<Vec<SiteGenotype>> = pool.install(|| {
        ranges
            .par_iter()
            .map(|&(lo, hi)| {
                cancel.check()?;
                denovo_chunk_diploid(bam_path, reference_path, contig, params, &ref_seq, length, lo, hi)
            })
            .collect::<Result<Vec<_>, AnalysisError>>()
    })?;
    Ok(nested.into_iter().flatten().collect())
}

/// De-novo diploid SNV calls for one emit range (mirrors [`denovo_chunk`], but genotypes ploidy 2).
#[allow(clippy::too_many_arguments)]
fn denovo_chunk_diploid(
    bam_path: &Path,
    reference_path: &Path,
    contig: &str,
    params: &HaploidCallerParams,
    ref_seq: &[u8],
    length: usize,
    emit_lo: usize,
    emit_hi: usize,
) -> Result<Vec<SiteGenotype>, AnalysisError> {
    let overlap = params.denovo_overlap;
    let proc_lo = emit_lo.saturating_sub(overlap).max(1);
    let proc_hi = (emit_hi + overlap).min(length);
    let ref_chunk = &ref_seq[(proc_lo - 1).min(ref_seq.len())..proc_hi.min(ref_seq.len())];

    let (mut counts, indel) = tally_region(bam_path, contig, params, proc_lo, proc_hi, Some(reference_path))?;
    if params.local_realign {
        realign_region(
            bam_path,
            contig,
            ref_chunk,
            proc_lo,
            &mut counts,
            &indel,
            params,
            Some(reference_path),
        )?;
    }

    let mut out = Vec::new();
    for pos in emit_lo..=emit_hi {
        let r = pos - proc_lo;
        let c = counts[r];
        let depth: u32 = c.iter().sum();
        if depth < params.min_depth {
            continue;
        }
        let ref_base = ref_chunk.get(r).copied().unwrap_or(b'N');
        let Some(ref_bi) = base_index(ref_base) else { continue }; // reference N/ambiguous
        let ref_byte = BASES[ref_bi];
        let ref_count = c[ref_bi];
        // Every non-reference base that clears the support floor is a candidate alt. The most
        // common one comes first.
        let mut alts: Vec<(usize, u32)> = c
            .iter()
            .enumerate()
            .filter(|&(bi, &n)| bi != ref_bi && n >= DENOVO_MIN_ALT_READS)
            .map(|(bi, &n)| (bi, n))
            .collect();
        if alts.is_empty() {
            continue; // hom-ref (no alt above the floor) — not emitted
        }
        alts.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));

        if alts.len() == 1 {
            // Biallelic: synthesize observations at a nominal quality (see DENOVO_DIPLOID_Q).
            let (alt_bi, alt_count) = alts[0];
            let alt_byte = BASES[alt_bi];
            let mut obs: Vec<(u8, u8)> = Vec::with_capacity((ref_count + alt_count) as usize);
            obs.extend(std::iter::repeat((ref_byte, DENOVO_DIPLOID_Q)).take(ref_count as usize));
            obs.extend(std::iter::repeat((alt_byte, DENOVO_DIPLOID_Q)).take(alt_count as usize));
            let g = genotype::call_genotype(&obs, ref_byte, alt_byte, 2, params.min_depth);
            if g.dosage < 1 {
                continue; // hom-ref or no-call — not a variant record
            }
            out.push(SiteGenotype {
                name: String::new(),
                contig: contig.to_string(),
                position: pos as i64,
                reference_allele: (ref_byte as char).to_string(),
                alternate_allele: (alt_byte as char).to_string(),
                ploidy: 2,
                dosage: g.dosage,
                gq: g.gq,
                depth,
                ref_depth: ref_count,
                alt_depth: alt_count,
                pls: g.pls,
                gt: None,
                allele_depths: None,
            });
            continue;
        }

        // Multiallelic SNV: ref = allele 0, each candidate alt = 1.. (in `alts` order).
        let mut obs: Vec<(usize, u8)> = Vec::with_capacity(depth as usize);
        obs.extend(std::iter::repeat((0usize, DENOVO_DIPLOID_Q)).take(ref_count as usize));
        for (k, &(_, n)) in alts.iter().enumerate() {
            obs.extend(std::iter::repeat((k + 1, DENOVO_DIPLOID_Q)).take(n as usize));
        }
        let mg = genotype::call_genotype_multi(&obs, alts.len() + 1, params.min_depth);
        if mg.gt == (0, 0) {
            continue; // hom-ref — not a variant record
        }
        let alt_depth: u32 = mg.allele_depths.iter().skip(1).sum();
        let dosage = (mg.gt.0 > 0) as i32 + (mg.gt.1 > 0) as i32; // alt-allele count fallback
        out.push(SiteGenotype {
            name: String::new(),
            contig: contig.to_string(),
            position: pos as i64,
            reference_allele: (ref_byte as char).to_string(),
            alternate_allele: alts
                .iter()
                .map(|&(bi, _)| (BASES[bi] as char).to_string())
                .collect::<Vec<_>>()
                .join(","),
            ploidy: 2,
            dosage,
            gq: mg.gq,
            depth,
            ref_depth: *mg.allele_depths.first().unwrap_or(&0),
            alt_depth,
            pls: mg.pls.clone(),
            gt: Some(format!("{}/{}", mg.gt.0, mg.gt.1)),
            allele_depths: Some(mg.allele_depths),
        });
    }

    // Indel pass over the active (indel-evidence) windows; merge into position order.
    let mut indels = indels_in_chunk(
        bam_path,
        contig,
        params,
        proc_lo,
        ref_chunk,
        emit_lo,
        emit_hi,
        &indel,
        Some(reference_path),
    )?;
    out.append(&mut indels);
    out.sort_by_key(|c| c.position);
    Ok(out)
}

/// A candidate indel allele relative to the reference: an insertion of these (uppercased) bases, or
/// a deletion of this many reference bases.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
enum IndelAllele {
    Ins(Vec<u8>),
    Del(u32),
}

/// Left-align an indel inside the repeat structure of the reference. This is the VCF
/// normalization. An aligner may put an indel anywhere inside a homopolymer or an STR run, but
/// the canonical form is the leftmost one. Returns the normalized `anchor`, which is 1-based, and
/// the allele.
///
/// A deletion of `len` bases at `[anchor, anchor+len-1]` moves left while
/// `ref[anchor-1] == ref[anchor+len-1]`. An insertion before `anchor` moves left while
/// `ref[anchor-1]` equals its last base, and the bases turn around the allele. `proc_lo` bounds
/// the move, and it is the 1-based start of the reference window that the code loaded.
fn left_normalize(anchor: i64, allele: &IndelAllele, ref_chunk: &[u8], proc_lo: usize) -> (i64, IndelAllele) {
    let at = |p: i64| -> Option<u8> {
        let i = p - proc_lo as i64;
        (i >= 0 && (i as usize) < ref_chunk.len()).then(|| ref_chunk[i as usize].to_ascii_uppercase())
    };
    match allele {
        IndelAllele::Del(len) => {
            let l = *len as i64;
            let mut a = anchor;
            while a > proc_lo as i64 && at(a - 1).is_some() && at(a - 1) == at(a + l - 1) {
                a -= 1;
            }
            (a, IndelAllele::Del(*len))
        }
        IndelAllele::Ins(seq) => {
            let mut a = anchor;
            let mut s: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();
            while a > proc_lo as i64 && !s.is_empty() && at(a - 1) == s.last().copied() {
                let last = s.pop().unwrap();
                s.insert(0, last); // rotate right: the inserted unit slides left by one ref base
                a -= 1;
            }
            (a, IndelAllele::Ins(s))
        }
    }
}

/// Build a diploid indel [`SiteGenotype`] in the VCF style, left-anchored at `emit_pos`, from the
/// read support for the reference against the indel. It genotypes at ploidy 2, with the GL over
/// sentinel bytes: `b'R'` for a read that covers the reference, and `b'A'` for a read that carries
/// the indel. `ref_byte` is the reference base at `emit_pos`. `deleted` holds the deleted
/// reference bases, and it is empty for an insertion. Returns `None` for a hom-ref call and for a
/// no-call.
#[allow(clippy::too_many_arguments)]
fn indel_site_genotype(
    contig: &str,
    emit_pos: i64,
    ref_byte: u8,
    allele: &IndelAllele,
    deleted: &[u8],
    ref_count: u32,
    alt_count: u32,
    params: &HaploidCallerParams,
) -> Option<SiteGenotype> {
    let r = (ref_byte as char).to_ascii_uppercase();
    let (reference_allele, alternate_allele) = match allele {
        IndelAllele::Ins(seq) => {
            // POS=anchor-1, REF=anchor base, ALT=anchor base + inserted bases.
            (
                r.to_string(),
                format!("{r}{}", String::from_utf8_lossy(seq).to_ascii_uppercase()),
            )
        }
        IndelAllele::Del(_) => {
            // POS=anchor-1, REF=anchor base + deleted bases, ALT=anchor base.
            (
                format!("{r}{}", String::from_utf8_lossy(deleted).to_ascii_uppercase()),
                r.to_string(),
            )
        }
    };
    let mut obs: Vec<(u8, u8)> = Vec::with_capacity((ref_count + alt_count) as usize);
    obs.extend(std::iter::repeat((b'R', DENOVO_DIPLOID_Q)).take(ref_count as usize));
    obs.extend(std::iter::repeat((b'A', DENOVO_DIPLOID_Q)).take(alt_count as usize));
    let g = genotype::call_genotype(&obs, b'R', b'A', 2, params.min_depth);
    if g.dosage < 1 {
        return None; // hom-ref or no-call — not a variant record
    }
    Some(SiteGenotype {
        name: String::new(),
        contig: contig.to_string(),
        position: emit_pos,
        reference_allele,
        alternate_allele,
        ploidy: 2,
        dosage: g.dosage,
        gq: g.gq,
        depth: ref_count + alt_count,
        ref_depth: ref_count,
        alt_depth: alt_count,
        pls: g.pls,
        gt: None,
        allele_depths: None,
    })
}

/// The de-novo diploid **indel** calls for this chunk. An active window is a window with indel
/// evidence. Over each of those, the code takes the indel allele of each read, from the CIGAR I
/// and D operations. It also takes the support of the reads that cover the reference. It tallies
/// the most common allele at each locus, and v1 keeps two alleles. It then genotypes that locus
/// at ploidy 2.
///
/// It emits a locus only when the VCF position of that locus lies in the emit range, which
/// removes a duplicate across a chunk boundary. The anchor is on the left, by the standard VCF
/// convention.
#[allow(clippy::too_many_arguments)]
fn indels_in_chunk(
    bam_path: &Path,
    contig: &str,
    params: &HaploidCallerParams,
    proc_lo: usize,
    ref_chunk: &[u8],
    emit_lo: usize,
    emit_hi: usize,
    indel_evidence: &[u32],
    reference: Option<&Path>,
) -> Result<Vec<SiteGenotype>, AnalysisError> {
    let windows = active_windows(indel_evidence, params.realign_min_indel_reads, params.realign_pad);
    if windows.is_empty() {
        return Ok(Vec::new());
    }
    let (header, mut reader) = reader::open_indexed(bam_path, reference)?;
    let mut out = Vec::new();

    /// One read's reference span + the indel events anchored in this window.
    struct ReadSpan {
        start: i64,
        ref_end: i64,
        events: Vec<(i64, IndelAllele, i64)>, // (anchor 1-based, allele, locus_end 1-based)
    }

    for (w0, w1) in windows {
        let (wlo, whi) = (proc_lo + w0, proc_lo + w1); // 1-based inclusive
        let region: Region = format!("{contig}:{wlo}-{whi}")
            .parse()
            .map_err(|_| AnalysisError::Message(format!("bad region for {contig}")))?;
        let mut reads: Vec<ReadSpan> = Vec::new();
        for result in reader.query(&header, &region)? {
            let record = result?;
            if !passes(&record, params) {
                continue;
            }
            let start = match record.alignment_start() {
                Some(p) => p.get() as i64,
                None => continue,
            };
            let seq = record.sequence();
            let mut ref_pos = start;
            let mut query_off = 0usize;
            let mut events = Vec::new();
            for op in record.cigar().as_ref() {
                let (kind, len) = (op.kind(), op.len());
                match (kind.consumes_reference(), kind.consumes_read()) {
                    (true, true) => {
                        ref_pos += len as i64;
                        query_off += len;
                    }
                    (true, false) => {
                        let anchor = ref_pos; // first deleted ref position (1-based)
                        if (wlo as i64) <= anchor && anchor <= (whi as i64) {
                            events.push((anchor, IndelAllele::Del(len as u32), anchor + len as i64 - 1));
                        }
                        ref_pos += len as i64;
                    }
                    (false, true) => {
                        if kind == Kind::Insertion {
                            let anchor = ref_pos; // insertion precedes this ref position
                            if (wlo as i64) <= anchor && anchor <= (whi as i64) {
                                let s: Vec<u8> = (0..len)
                                    .filter_map(|i| seq.get(query_off + i).map(|b| b.to_ascii_uppercase()))
                                    .collect();
                                events.push((anchor, IndelAllele::Ins(s), anchor));
                            }
                        }
                        query_off += len;
                    }
                    (false, false) => {}
                }
            }
            reads.push(ReadSpan {
                start,
                ref_end: ref_pos - 1,
                events,
            });
        }

        // Tally candidate alleles, normalize each, and group by normalized VCF position so that
        // co-located alleles (compound-het indels) become a single multiallelic record.
        let mut tally: HashMap<(i64, IndelAllele), u32> = HashMap::new();
        for r in &reads {
            for (anchor, al, _) in &r.events {
                *tally.entry((*anchor, al.clone())).or_insert(0) += 1;
            }
        }
        /// A normalized candidate allele grouped at its emit position. `anchor`/`allele` are the
        /// *original* (pre-normalization) key used to match reads; `nal` is the canonical allele.
        struct Cand {
            anchor: i64,
            allele: IndelAllele,
            nal: IndelAllele,
            locus_end: i64, // original locus end — for ref-span support
            count: u32,
        }
        let mut groups: HashMap<i64, Vec<Cand>> = HashMap::new();
        for ((anchor, al), &count) in &tally {
            if count < params.realign_min_indel_reads {
                continue; // sub-threshold noise allele
            }
            let locus_end = match al {
                IndelAllele::Del(l) => anchor + *l as i64 - 1,
                IndelAllele::Ins(_) => *anchor,
            };
            // Left-align within the reference repeat for the canonical VCF position/alleles.
            let (na, nal) = left_normalize(*anchor, al, ref_chunk, proc_lo);
            let emit_pos = na - 1;
            if emit_pos < emit_lo as i64 || emit_pos > emit_hi as i64 {
                continue; // assigned to whichever chunk owns the normalized position (no dup/loss)
            }
            groups.entry(emit_pos).or_default().push(Cand {
                anchor: *anchor,
                allele: al.clone(),
                nal,
                locus_end,
                count,
            });
        }

        for (emit_pos, mut cands) in groups {
            let idx = emit_pos - proc_lo as i64;
            if idx < 0 || idx as usize >= ref_chunk.len() {
                continue;
            }
            let ref_byte = ref_chunk[idx as usize];
            if base_index(ref_byte).is_none() {
                continue; // ambiguous anchor base
            }
            let na = emit_pos + 1;

            if cands.len() == 1 {
                // Biallelic: ref support uses the reads' *actual* indel locus (normalization only
                // changes the VCF representation, not which reads support the allele).
                let c = &cands[0];
                let spanning = reads
                    .iter()
                    .filter(|r| r.start < c.anchor && r.ref_end >= c.locus_end)
                    .count() as u32;
                let ref_count = spanning.saturating_sub(c.count);
                if ref_count + c.count < params.min_depth {
                    continue;
                }
                let deleted: Vec<u8> = if let IndelAllele::Del(len) = &c.nal {
                    let (s, e) = (
                        (na - proc_lo as i64).max(0) as usize,
                        (na - proc_lo as i64 + *len as i64).max(0) as usize,
                    );
                    if s < e && e <= ref_chunk.len() {
                        ref_chunk[s..e].to_vec()
                    } else {
                        continue; // deleted span runs off the loaded reference window
                    }
                } else {
                    Vec::new()
                };
                if let Some(g) =
                    indel_site_genotype(contig, emit_pos, ref_byte, &c.nal, &deleted, ref_count, c.count, params)
                {
                    out.push(g);
                }
                continue;
            }

            // The site has more than two alleles. There is one common REF, which covers the
            // largest deletion, and one ALT for each allele.
            cands.sort_by(|a, b| b.count.cmp(&a.count).then(a.anchor.cmp(&b.anchor))); // dominant first, deterministic
            let maxdel = cands
                .iter()
                .filter_map(|c| {
                    if let IndelAllele::Del(l) = &c.nal {
                        Some(*l as usize)
                    } else {
                        None
                    }
                })
                .max()
                .unwrap_or(0);
            let ref_lo = idx as usize;
            let ref_hi = ref_lo + 1 + maxdel;
            if ref_hi > ref_chunk.len() {
                continue; // REF span runs off the loaded reference window
            }
            let common_ref = ref_chunk[ref_lo..ref_hi].to_ascii_uppercase();
            let tail = &common_ref[1..]; // the `maxdel` reference bases after the anchor
            let anchor_byte = ref_byte.to_ascii_uppercase();
            let alts: Vec<String> = cands
                .iter()
                .map(|c| {
                    let mut v = vec![anchor_byte];
                    match &c.nal {
                        IndelAllele::Ins(seq) => {
                            v.extend(seq.iter().map(|b| b.to_ascii_uppercase()));
                            v.extend_from_slice(tail); // keep the bases a co-located deletion would remove
                        }
                        IndelAllele::Del(l) => v.extend_from_slice(&tail[(*l as usize).min(tail.len())..]),
                    }
                    String::from_utf8_lossy(&v).into_owned()
                })
                .collect();

            // Assign each read to ref (0) or a candidate allele (k+1); synthesize observations.
            let anchor_min = cands.iter().map(|c| c.anchor).min().unwrap();
            let locus_end_max = cands.iter().map(|c| c.locus_end).max().unwrap();
            let mut obs: Vec<(usize, u8)> = Vec::new();
            for r in &reads {
                let carried = r
                    .events
                    .iter()
                    .find_map(|(a, al, _)| cands.iter().position(|c| c.anchor == *a && &c.allele == al));
                match carried {
                    Some(k) => obs.push((k + 1, DENOVO_DIPLOID_Q)),
                    None => {
                        // Ref only if it spans the locus and carries no (other) indel here.
                        let other_indel = r
                            .events
                            .iter()
                            .any(|(a, _, le)| *a <= locus_end_max && *le >= anchor_min);
                        if !other_indel && r.start < anchor_min && r.ref_end >= locus_end_max {
                            obs.push((0, DENOVO_DIPLOID_Q));
                        }
                    }
                }
            }
            if (obs.len() as u32) < params.min_depth {
                continue;
            }
            let mg = genotype::call_genotype_multi(&obs, cands.len() + 1, params.min_depth);
            if mg.gt == (0, 0) {
                continue; // hom-ref — not a variant record
            }
            let alt_depth: u32 = mg.allele_depths.iter().skip(1).sum();
            let dosage = (mg.gt.0 > 0) as i32 + (mg.gt.1 > 0) as i32; // alt-allele count fallback for biallelic consumers
            out.push(SiteGenotype {
                name: String::new(),
                contig: contig.to_string(),
                position: emit_pos,
                reference_allele: String::from_utf8_lossy(&common_ref).into_owned(),
                alternate_allele: alts.join(","),
                ploidy: 2,
                dosage,
                gq: mg.gq,
                depth: obs.len() as u32,
                ref_depth: *mg.allele_depths.first().unwrap_or(&0),
                alt_depth,
                pls: mg.pls.clone(),
                gt: Some(format!("{}/{}", mg.gt.0, mg.gt.1)),
                allele_depths: Some(mg.allele_depths),
            });
        }
    }
    Ok(out)
}

/// Maximal runs of positions with enough indel evidence, each padded by `pad` and
/// merged where they touch. Returns 0-based inclusive `(start, end)` reference windows.
fn active_windows(indel_evidence: &[u32], min_reads: u32, pad: i64) -> Vec<(usize, usize)> {
    let len = indel_evidence.len();
    let mut windows: Vec<(usize, usize)> = Vec::new();
    let mut i = 0usize;
    while i < len {
        if indel_evidence[i] >= min_reads {
            let run_start = i;
            while i < len && indel_evidence[i] >= min_reads {
                i += 1;
            }
            let run_end = i - 1;
            let w0 = (run_start as i64 - pad).max(0) as usize;
            let w1 = ((run_end as i64 + pad) as usize).min(len - 1);
            match windows.last_mut() {
                Some(last) if w0 <= last.1 + 1 => last.1 = last.1.max(w1),
                _ => windows.push((w0, w1)),
            }
        } else {
            i += 1;
        }
    }
    windows
}

/// Fit the reads in each window with indel evidence onto the reference again, and replace the
/// tally over those windows. The index into the arrays is relative to `region_lo`, which is
/// 1-based.
#[allow(clippy::too_many_arguments)]
fn realign_region(
    bam_path: &Path,
    contig: &str,
    ref_chunk: &[u8],
    region_lo: usize,
    counts: &mut [[u32; 4]],
    indel_evidence: &[u32],
    params: &HaploidCallerParams,
    reference: Option<&Path>,
) -> Result<(), AnalysisError> {
    let windows = active_windows(indel_evidence, params.realign_min_indel_reads, params.realign_pad);
    // Keep only in-range windows; they stay sorted + disjoint, so by both `w0` and `w1`.
    let windows: Vec<(usize, usize)> = windows.into_iter().filter(|&(_, w1)| w1 < ref_chunk.len()).collect();
    if windows.is_empty() {
        return Ok(());
    }
    let mut win_counts: Vec<Vec<[u32; 4]>> = windows.iter().map(|&(w0, w1)| vec![[0u32; 4]; w1 - w0 + 1]).collect();

    // ONE query against the index, over all of the active windows. Decode the reads of the
    // region once, and send each read to the window or windows that it overlaps.
    //
    // The code before this one ran a query for each window. A contig with many repeats holds
    // thousands of indel windows. So that code decoded the same CRAM containers again and again,
    // on the hot path of the de-novo pass. The reads are short, so each one overlaps only
    // one window or two, which a binary search over the sorted windows finds.
    let span_lo = region_lo + windows.first().unwrap().0;
    let span_hi = region_lo + windows.last().unwrap().1;
    let region: Region = format!("{contig}:{span_lo}-{span_hi}")
        .parse()
        .map_err(|_| AnalysisError::Message(format!("bad region for {contig}")))?;
    let (header, mut reader) = reader::open_indexed(bam_path, reference)?;
    for result in reader.query(&header, &region)? {
        let record = result?;
        if !passes(&record, params) {
            continue;
        }
        let start = match record.alignment_start() {
            Some(p) => p.get(),
            None => continue,
        };
        let ref_span: usize = record
            .cigar()
            .as_ref()
            .iter()
            .filter(|op| op.kind().consumes_reference())
            .map(|op| op.len())
            .sum();
        let read_end = start + ref_span.saturating_sub(1); // 1-based inclusive

        // First window whose end reaches the read's start, then walk while its start is still
        // within the read (windows are disjoint + sorted).
        let start_rel = start as i64 - region_lo as i64;
        let mut iw = windows.partition_point(|&(_, w1)| (w1 as i64) < start_rel);
        while iw < windows.len() {
            let (w0, w1) = windows[iw];
            let wlo_abs = region_lo + w0;
            if wlo_abs > read_end {
                break;
            }
            let whi_abs = region_lo + w1;
            let target = &ref_chunk[w0..=w1];
            let (qbases, qquals) = window_substring(&record, start, wlo_abs, whi_abs)?;
            if !qbases.is_empty() {
                let (tstart, ops) = realign::fitting_align(&qbases, target);
                for (ref_idx, base, qual) in realign::project(&qbases, &qquals, w0, tstart, &ops) {
                    if qual >= params.min_base_quality {
                        if let Some(bi) = base_index(base) {
                            win_counts[iw][ref_idx - w0][bi] += 1;
                        }
                    }
                }
            }
            iw += 1;
        }
    }

    for (iw, &(w0, _)) in windows.iter().enumerate() {
        for (k, c) in std::mem::take(&mut win_counts[iw]).into_iter().enumerate() {
            counts[w0 + k] = c;
        }
    }
    Ok(())
}

/// Extract a read's bases + qualities over the 1-based reference window `[wlo, whi]`,
/// in reference order, including any inserted bases anchored inside the window.
fn window_substring(
    record: &RecordBuf,
    start: usize,
    wlo: usize,
    whi: usize,
) -> Result<(Vec<u8>, Vec<u8>), AnalysisError> {
    let seq = record.sequence();
    let quals = record.quality_scores();
    let quals = quals.as_ref();
    let mut bases = Vec::new();
    let mut q = Vec::new();
    let mut ref_pos = start; // 1-based
    let mut query_off = 0usize;
    for op in record.cigar().as_ref() {
        let kind = op.kind();
        let len = op.len();
        match (kind.consumes_reference(), kind.consumes_read()) {
            (true, true) => {
                for i in 0..len {
                    let pos = ref_pos + i;
                    if pos >= wlo && pos <= whi {
                        if let Some(b) = seq.get(query_off + i) {
                            bases.push(b);
                            q.push(quals.get(query_off + i).copied().unwrap_or(0));
                        }
                    }
                }
                ref_pos += len;
                query_off += len;
            }
            (true, false) => ref_pos += len,
            (false, true) => {
                // Insertion anchored at ref_pos: include if inside the window.
                if kind == Kind::Insertion && ref_pos >= wlo && ref_pos <= whi {
                    for i in 0..len {
                        if let Some(b) = seq.get(query_off + i) {
                            bases.push(b);
                            q.push(quals.get(query_off + i).copied().unwrap_or(0));
                        }
                    }
                }
                query_off += len;
            }
            (false, false) => {}
        }
    }
    Ok((bases, q))
}

/// Subtract known tree positions from de-novo calls to yield the private variant set
/// (the role `PrivateSnpProcessor` plays after liftover of the tree loci).
pub fn subtract_known(calls: &[VariantCall], known_positions: &HashSet<i64>) -> Vec<VariantCall> {
    calls
        .iter()
        .filter(|v| !known_positions.contains(&v.position))
        .cloned()
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sg(contig: &str, pos: i64, alt: &str, dosage: i32, depth: u32, gq: u8) -> SiteGenotype {
        SiteGenotype {
            name: String::new(),
            contig: contig.into(),
            position: pos,
            reference_allele: "A".into(),
            alternate_allele: alt.into(),
            ploidy: 2,
            dosage,
            gq,
            depth,
            ref_depth: depth.saturating_sub(depth * dosage as u32 / 2),
            alt_depth: depth * dosage as u32 / 2,
            pls: vec![],
            gt: None,
            allele_depths: None,
        }
    }

    #[test]
    fn consensus_reconcile_resolves_homref_and_abstains_on_no_call() {
        // Site 100: run A is het (0/1) and deep, run B is hom-ref (0/0) and deep. That is a real
        // disagreement. In the vote by depth both are deep and their weights are equal, so the
        // lower dosage breaks the tie. The result is hom-ref, and the code does NOT emit it.
        //
        // Site 200: run A is het and deep, run B has no call, because its depth of 1 is below the
        // minimum of 4. B then abstains, A wins, and the result is het.
        //
        // Site 300: both are hom-alt (1/1). The result is hom-alt, with the depths added.
        let a = vec![
            sg("chr1", 100, "G", 1, 30, 50),
            sg("chr1", 200, "G", 1, 30, 50),
            sg("chr1", 300, "T", 2, 20, 60),
        ];
        let b = vec![
            sg("chr1", 100, "G", 0, 30, 50),
            sg("chr1", 200, "G", 0, 1, 0),
            sg("chr1", 300, "T", 2, 25, 55),
        ];
        let out = reconcile_site_genotypes(&[a, b], 4);

        // 100 → hom-ref consensus, not a variant (absent).
        assert!(!out.iter().any(|g| g.position == 100));
        // 200 → het (B abstained, only A's deep het counts).
        let s200 = out.iter().find(|g| g.position == 200).expect("200 emitted");
        assert_eq!(s200.dosage, 1);
        // 300 → hom-alt with summed depth.
        let s300 = out.iter().find(|g| g.position == 300).expect("300 emitted");
        assert_eq!(s300.dosage, 2);
        assert_eq!(s300.depth, 45);
    }

    #[test]
    fn left_normalize_shifts_indels_into_repeats() {
        // ref_chunk starts at proc_lo=1 (1-based). "GAAAAC" → positions 1G 2A 3A 4A 5A 6C.
        let refc = b"GAAAAC";
        // A 1bp deletion reported at the last A (anchor 5) left-aligns to the first A (anchor 2).
        let (a, al) = left_normalize(5, &IndelAllele::Del(1), refc, 1);
        assert_eq!((a, al), (2, IndelAllele::Del(1)));
        // An insertion of "A" before anchor 5 (in the A-run) left-aligns to anchor 2.
        let (a, al) = left_normalize(5, &IndelAllele::Ins(b"A".to_vec()), refc, 1);
        assert_eq!((a, al), (2, IndelAllele::Ins(b"A".to_vec())));
        // A non-repeat deletion does not move: "ACGTC", delete the G (anchor 3).
        let (a, _) = left_normalize(3, &IndelAllele::Del(1), b"ACGTC", 1);
        assert_eq!(a, 3);
    }

    #[test]
    fn indel_site_genotype_builds_left_anchored_alleles() {
        let params = HaploidCallerParams::default();
        // Insertion of "TT" after the anchor base 'C' → REF=C, ALT=CTT; 10/10 → het 0/1.
        let ins = indel_site_genotype(
            "chr1",
            100,
            b'C',
            &IndelAllele::Ins(b"TT".to_vec()),
            &[],
            10,
            10,
            &params,
        )
        .unwrap();
        assert_eq!(
            (ins.reference_allele.as_str(), ins.alternate_allele.as_str()),
            ("C", "CTT")
        );
        assert_eq!((ins.position, ins.dosage), (100, 1));
        // Deletion of "CG" after anchor 'A' → REF=ACG, ALT=A; all-alt → hom-alt 1/1.
        let del = indel_site_genotype("chr1", 200, b'A', &IndelAllele::Del(2), b"CG", 0, 20, &params).unwrap();
        assert_eq!(
            (del.reference_allele.as_str(), del.alternate_allele.as_str()),
            ("ACG", "A")
        );
        assert_eq!(del.dosage, 2);
        // No alt support → hom-ref → not emitted.
        assert!(indel_site_genotype("chr1", 300, b'A', &IndelAllele::Del(1), b"C", 20, 0, &params).is_none());
    }

    #[test]
    fn active_windows_pads_and_merges_indel_runs() {
        // evidence at idx 50 and 52 (>=3 reads); pad 5 -> [45,57] (the two merge).
        let mut ev = vec![0u32; 100];
        ev[50] = 4;
        ev[52] = 3;
        ev[90] = 1; // below threshold -> ignored
        let w = active_windows(&ev, 3, 5);
        assert_eq!(w, vec![(45, 57)]);
        // higher threshold drops everything.
        assert!(active_windows(&ev, 5, 5).is_empty());
    }

    #[test]
    fn consensus_breaks_ties_toward_earlier_base() {
        assert_eq!(consensus(&[3, 3, 0, 0]), (0, 3)); // A wins tie vs C
        assert_eq!(consensus(&[0, 1, 5, 2]), (2, 5)); // G
        assert_eq!(consensus(&[0, 0, 0, 0]), (0, 0)); // empty -> A, count 0
    }

    #[test]
    fn paralog_filter_drops_only_bi_allelic_haploid_sites() {
        let p = HaploidCallerParams::default(); // max_minor 0.20, min_minor_reads 2
        let para = |c: [u32; 4]| is_paralogous(&c, c.iter().sum(), &p);

        // Clean monoallelic call (HiFi-like): not paralogous.
        assert!(!para([11, 0, 0, 0]));
        // One read that disagrees, at low depth. That is a sequencing error, and the site stays.
        assert!(!para([3, 1, 0, 0])); // second=1 (< 2 reads)

        // Errors spread over the other bases, and none of them gets to 2 reads. The site stays.
        assert!(!para([18, 1, 1, 0])); // second=1

        // A true pileup with two alleles, at 7 derived and 4 ancestral. That is a paralog, and
        // the code drops it.
        assert!(para([7, 4, 0, 0])); // second=4, 0.36 > 0.20

        // The boundary. 2/10 = 0.20 is not above the threshold, so the site stays.
        assert!(!para([8, 2, 0, 0]));
        // 3/10 = 0.30, which is more than 0.20, with 3 reads. The code drops the site.
        assert!(para([7, 3, 0, 0]));
        // Empty pileup is never paralogous.
        assert!(!para([0, 0, 0, 0]));
    }

    #[test]
    fn paralog_filter_disables_at_fraction_one() {
        let p = HaploidCallerParams {
            max_minor_allele_fraction: 1.0,
            ..Default::default()
        };
        // With the filter off, even a 50/50 split gets no flag.
        assert!(!is_paralogous(&[5, 5, 0, 0], 10, &p));
    }

    #[test]
    fn subtract_known_removes_listed_positions() {
        let v = |p| VariantCall {
            contig: "chrM".into(),
            position: p,
            reference_allele: 'C',
            alternate_allele: 'A',
            depth: 4,
            alt_depth: 4,
            allele_fraction: 1.0,
            quality: None,
        };
        let calls = vec![v(2), v(3), v(4)];
        let known: HashSet<i64> = [2, 3].into_iter().collect();
        let private = subtract_known(&calls, &known);
        assert_eq!(private.iter().map(|c| c.position).collect::<Vec<_>>(), vec![4]);
    }
}
