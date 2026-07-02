//! Purpose-built **haploid** variant caller (plan §4b) — the GATK replacement for
//! Y/mtDNA. There is no pure-code Scala caller to port: the legacy app shelled out to
//! GATK `HaplotypeCaller --sample-ploidy 1` (force-call at tree sites) and `Mutect2
//! --mitochondria` / haploid `HaplotypeCaller` (de-novo discovery), then subtracted
//! known tree positions to get private variants. This module reproduces both modes by
//! **pileup-consensus calling**, which is tractable precisely because Y and mtDNA are
//! haploid (ploidy 1) — no diploid local reassembly.
//!
//! Two modes:
//! 1. [`force_call_sites`] — genotype-given-alleles at known tree `Site`s (haplogroup
//!    assignment): pileup, take the consensus base, report whether it is the site's
//!    ref or alt allele.
//! 2. [`call_denovo`] — walk the contig, emit positions whose consensus base differs
//!    from the reference (the candidate private variants). [`subtract_known`] removes
//!    known tree positions to yield the private set.
//!
//! **v1 is SNP-only** (plan §4b): indels/homopolymers are where naive pileup calling
//! diverges from GATK (light local realignment is the planned mitigation), so indel
//! alleles are skipped here and treated as advisory until the §4c parity harness
//! validates them. Defaults are starting points the harness will tune.
//!
//! Memory: de-novo processes the contig in overlapping chunks (`denovo_chunk`), so the
//! dense per-position tally is bounded by the chunk, not the contig length. Both-side
//! context overlap keeps realignment windows that straddle a chunk boundary fully
//! visible. Force-call tallies only the target sites (sparse), cheap regardless of size.

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

/// Algorithm version for de-novo caller artifacts; bump on output-affecting changes
/// (e.g. the local-realignment addition bumped this to -2).
pub const DENOVO_VERSION: &str = "haploid-denovo-2";

/// Algorithm version for site-genotype (panel) artifacts.
pub const GENOTYPE_VERSION: &str = "genotype-1";

/// Parameters for haploid calling. Defaults are v1 starting points (gated by §4c).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HaploidCallerParams {
    /// Minimum passing depth (reads clearing the quality filters) to make any call.
    pub min_depth: u32,
    /// Reads below this MAPQ are dropped entirely.
    pub min_mapping_quality: u8,
    /// Bases below this quality are not counted.
    pub min_base_quality: u8,
    /// The consensus base must be at least this fraction of passing depth to call.
    pub min_allele_fraction: f64,
    /// Allele-balance (paralog) filter for haploid sites. A true Y/mt-haploid site is
    /// near-monoallelic; a substantial *second* allele signals paralog/mismapping (two loci
    /// piled together) and the site is dropped. Tripped only when the second-most-common
    /// allele has both at least `min_paralog_minor_reads` reads AND a fraction strictly above
    /// `max_minor_allele_fraction` — a lone discordant read (sequencing error) does not trip
    /// it. Set the fraction `>= 1.0` to disable. See PangenomeExpansion.md (Phase 1).
    pub max_minor_allele_fraction: f64,
    /// Minimum second-allele read count for the paralog filter to engage (guards low depth).
    pub min_paralog_minor_reads: u32,
    /// Run light local realignment around candidate indels before de-novo calling.
    pub local_realign: bool,
    /// Minimum reads with indel evidence at a position to open a realignment window.
    pub realign_min_indel_reads: u32,
    /// Padding (bp) added around indel-evidence runs to form a realignment window.
    pub realign_pad: i64,
    /// De-novo emit chunk size (bp). The contig is processed in chunks so memory is
    /// bounded; a chunk holds dense arrays for `chunk + 2*overlap` positions.
    pub denovo_chunk: usize,
    /// Context overlap (bp) processed on each side of a chunk, so realignment windows
    /// straddling a chunk boundary are still fully seen. Must exceed `realign_pad`.
    pub denovo_overlap: usize,
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

/// A diploid/haploid genotype at a known site (genotype-likelihood model). `dosage` is
/// the alt-allele count (0..=ploidy), or -1 for a no-call — the encoding the
/// population/ancestry/IBD paths consume.
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
    /// Explicit VCF genotype string (e.g. `"1/2"`) for multiallelic sites. When `None`, the
    /// genotype is derived from `dosage` (biallelic). Additive — old cached blobs decode to `None`.
    #[serde(default)]
    pub gt: Option<String>,
    /// Per-allele read depths `[ref, alt1, alt2, …]` for multiallelic sites; `None` → biallelic
    /// (use `ref_depth`/`alt_depth`).
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
}

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

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

/// Allele-balance / paralog filter for a haploid pileup. A true haploid site is near-
/// monoallelic; when the second-most-common allele carries both enough reads
/// (`min_paralog_minor_reads`) and enough fraction (strictly above `max_minor_allele_fraction`)
/// the site looks bi-allelic — a paralog/mismapping artifact — and the caller should drop it.
/// A single discordant read (likely sequencing error) does not trip it.
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

/// Resolve a contig's length by opening the alignment header at `bam_path`. `reference` is
/// required for CRAM.
pub(crate) fn read_contig_length(
    bam_path: &Path,
    contig: &str,
    reference: Option<&Path>,
) -> Result<usize, AnalysisError> {
    let header = reader::read_header(bam_path, reference)?;
    contig_length(&header, contig).ok_or_else(|| AnalysisError::Message(format!("contig {contig} not in BAM header")))
}

/// Load a contig's full reference sequence in one indexed-FASTA query (shared read-only across the
/// caller's chunks — each chunk slices its own window instead of re-querying).
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

/// The contig (reference-sequence) names in the alignment header. `reference` is required
/// for CRAM. Used to skip lifted positions that land on contigs the alignment lacks.
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

/// Call the consensus base at each 1-based `target` position on `contig` (haploid
/// genotyping for haplogroup assignment). A position is called only when it clears
/// `min_depth` passing reads and the consensus base is at least `min_allele_fraction` of
/// that depth; uncalled positions are simply absent. Returns position → uppercase base.
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

/// Diagnostic: the raw passing A/C/G/T read tally at each `target` (1-based) — the evidence
/// **behind** [`call_bases_at`]'s consensus pick, before the depth / allele-fraction / paralog
/// filters. Returns 1-based position → `[A, C, G, T]` counts (absent = no passing read covered it).
/// For "what do the reads actually show at this tree SNP" logging.
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

/// The expected indel allele for a tree locus given its (VCF left-anchored) ancestral/derived
/// alleles: an insertion of the trailing bases (`A`→`ATT` ⇒ Ins("TT")) or a deletion of the length
/// difference (`TA`→`T` ⇒ Del(1)). `None` for a SNP or a complex/non-left-anchored allele.
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

/// Walk one read's CIGAR from `start` (1-based), collecting each indel event as
/// `(anchor 1-based, allele)` — deletion anchor = first deleted ref base; insertion anchor = the ref
/// base the insertion precedes. Returns the events plus the read's inclusive reference end.
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

/// Targeted genotyping of tree **indel** loci. Each target is `(pos, ancestral, derived)` in VCF
/// left-anchored form (`pos` = the anchor base; e.g. `A`→`ATT` insertion, `TA`→`T` deletion). For
/// each locus, reads spanning it are examined: a read carrying the matching insertion/deletion
/// (after left-normalization into the reference repeat) supports the derived allele.
///
/// **Additive-only**: a locus with a clear derived majority over `min_depth` is emitted as
/// [`haplo::INDEL_DERIVED`] at `pos`; everything else — no indel support, low depth, or reads that
/// merely *span* the site cleanly — is left as **no-call**, never an ancestral contradiction. Indel
/// genotyping around homopolymers/STRs is noisy enough that a "clean-spanning" read is often just the
/// aligner's alternate representation of the same indel; calling those ancestral would spuriously
/// contradict sparse nodes (a d==0 node picking up one false ancestral trips the confident-divergence
/// guard and vetoes the whole lineage). So indels only ever *confirm* a branch, matching the intent:
/// cover the many indel-defined DecodingUs branches when the sample carries them. Requires a
/// `reference` (to left-normalize + know deleted bases); returns empty without one, or when the
/// contig isn't in the FASTA.
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
        pos: i64,        // VCF POS (anchor), 1-based
        n_anchor: i64,   // normalized CIGAR anchor (= pos+1 canonically)
        n_allele: IndelAllele,
        span_end: i64,   // last ref base the ref-spanning read must cover
    }
    let mut ptargets: Vec<PTarget> = Vec::new();
    for (pos, anc, der) in targets {
        let Some(al) = expected_indel_allele(anc, der) else { continue };
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

    // Single contig-wide pass (as the SNP tally does): walk every read once, and for each target the
    // read spans, accumulate matched (carries the indel) vs ref-spanning (spans it cleanly) support.
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
        let Some(start) = record.alignment_start().map(|p| p.get() as i64) else { continue };
        let (raw, ref_end) = read_indel_events(&record, start);
        let events: Vec<(i64, IndelAllele)> =
            raw.into_iter().map(|(a, al)| left_normalize(a, &al, &refbytes, 1)).collect();
        // Targets whose anchor this read could inform: pos in [start, ref_end].
        let lo = positions.partition_point(|&p| p < start);
        let hi = positions.partition_point(|&p| p <= ref_end);
        for i in lo..hi {
            let t = &ptargets[i];
            if events.iter().any(|(a, al)| *a == t.n_anchor && *al == t.n_allele) {
                matched[i] += 1;
            } else if start <= t.pos
                && ref_end >= t.span_end
                && !events.iter().any(|(a, _)| *a == t.n_anchor)
            {
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

/// Dense A/C/G/T tally + per-position indel evidence for the 1-based inclusive region
/// `[lo, hi]`, indexed by `pos - lo` (the chunked de-novo path).
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

/// Record one read's passing `(base, qual)` at each of `targets` (sorted, 1-based) that it covers,
/// in a **single** CIGAR walk from the alignment start. Targets in a deletion/skip/insertion/clip,
/// past the read, below `min_base_quality`, or non-ACGT are skipped. This is the multi-target
/// generalization of a per-site probe — one walk feeds many sites, so a long read shared by several
/// nearby panel sites is decoded + walked once instead of once per site.
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
            // Deletion / ref-skip — targets inside the gap carry no base.
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

/// Per-target-site passing `(base, qual)` observations (ACGT bases clearing the quality filters),
/// keyed by 1-based position — the input the genotype-likelihood model needs.
///
/// The targets are grouped into contiguous runs (split only where the gap between adjacent sites
/// exceeds a read length), and each run is fetched with a **single** streaming index query. So we
/// seek straight to the regions that hold targets — never scanning the whole contig — and decode
/// each read once (a point query per site re-fetches + re-converts the long HiFi reads that span
/// several nearby sites). Within a run, [`collect_bases`] distributes each read's bases to every
/// target it covers in one CIGAR walk.
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

    // Split into runs where consecutive sites are within MAX_GAP — beyond a read length no read can
    // span the gap, so splitting there is free (no shared reads lost) and skips read-free spans.
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

/// Genotype known SNP sites on `contig` at the given `ploidy` (1 = haploid Y/MT/male-X,
/// 2 = autosome / female-X) using the genotype-likelihood model — the panel-genotyping
/// path the population / ancestry / IBD analyses consume. Non-SNP sites are skipped.
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

/// Genotype `sites` across **every** contig they span, one contig per rayon task. The panel-genotyping
/// entry point for whole-genome panels (the per-contig [`genotype_sites`] is independent + IO-bound on
/// its own index region, so contigs parallelize cleanly). Results are concatenated (order across
/// contigs is unspecified — downstream consumers key by site, not order).
pub fn genotype_sites_all_contigs(
    bam_path: &Path,
    sites: &[Site],
    ploidy: u8,
    params: &HaploidCallerParams,
    reference: Option<&Path>,
) -> Result<Vec<SiteGenotype>, AnalysisError> {
    let contigs: Vec<&str> = sites
        .iter()
        .map(|s| s.contig.as_str())
        .collect::<std::collections::BTreeSet<&str>>()
        .into_iter()
        .collect();
    // Run on a decode-safe pool rather than rayon's global pool (2 MiB stacks): each task decodes
    // CRAM records, which recurse deeply on CRAM 3.1 and would otherwise overflow + abort. See
    // [`reader::decode_pool`].
    let pool = crate::reader::decode_pool(contigs.len().max(1).min(crate::unified::analysis_thread_count()))?;
    let per_contig: Result<Vec<Vec<SiteGenotype>>, AnalysisError> = pool.install(|| {
        contigs
            .into_par_iter()
            .map(|contig| genotype_sites(bam_path, contig, sites, ploidy, params, reference))
            .collect()
    });
    Ok(per_contig?.into_iter().flatten().collect())
}

/// Reconcile per-alignment force-call genotypes at a shared site set into one **consensus** diploid
/// genotype per site — the subject-level joint genotype across a person's WGS runs. Each input is
/// one alignment's [`SiteGenotype`]s at the *union* of variant sites (all on the same reference
/// build, so `(contig, position, ref, alt)` align). Per site, a depth-weighted vote over the dosage
/// classes {0,1,2}: an alignment whose depth is below `min_depth` is its no-call (excluded), so a
/// site absent-as-hom-ref in one run is a real vote (resolving "run A het vs run B hom-ref") while a
/// genuinely uncovered run abstains. Only **variant** consensus sites (het/hom-alt) are returned;
/// hom-ref / no-call consensus is not a variant. Depth/AD are summed and GQ is the max over the
/// supporting alignments; PLs are dropped (the per-run likelihoods don't compose into one PL here).
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
                // Depth-bonus weight, mirroring consensus::obs_weight's WGS term (constant method
                // factor drops out of the argmax).
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
        // argmax weight; tie → more raw supporting runs, then the lower dosage.
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

/// Force-call (genotype-given-alleles) at known SNP sites on `contig`. Non-SNP sites
/// (multi-base ref/alt) are skipped — v1 is SNP-only.
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

/// De-novo SNP discovery across `contig`, processed in overlapping chunks so memory is
/// bounded by the chunk (not the contig length). Emits positions whose consensus base
/// passes the depth/fraction filters and differs from the reference. Both-side context
/// overlap keeps realignment windows that straddle a chunk boundary fully visible.
pub fn call_denovo(
    bam_path: &Path,
    reference_path: &Path,
    contig: &str,
    params: &HaploidCallerParams,
) -> Result<Vec<VariantCall>, AnalysisError> {
    let length = read_contig_length(bam_path, contig, Some(reference_path))?;

    // Load the contig's reference once, shared read-only across chunks — each chunk slices its own
    // window instead of re-querying the FASTA.
    let ref_seq = load_contig_sequence(reference_path, contig, length)?;

    // Disjoint emit ranges, in order, each processed independently with its own indexed BAM
    // region query. Chunks stay large (`denovo_chunk`, default 8 MB): a CRAM container spans
    // several MB, so chunks smaller than a container would re-decode it in every overlapping
    // chunk. rayon caps in-flight chunks at the pool size, so peak memory is bounded.
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
            .map(|&(lo, hi)| denovo_chunk(bam_path, reference_path, contig, params, &ref_seq, length, lo, hi))
            .collect::<Result<Vec<_>, AnalysisError>>()
    })?;
    // Ranges are disjoint and collected in order, so flattening preserves global position order.
    Ok(nested.into_iter().flatten().collect())
}

/// De-novo SNP calls for one emit range `[emit_lo, emit_hi]` (1-based inclusive). Tallies a
/// `denovo_overlap`-padded window so realignment windows straddling the boundary are fully
/// seen, but emits only `[emit_lo, emit_hi]`. `ref_seq` is the full contig reference (index 0 =
/// position 1). Each call opens its own BAM reader, so it is independent and thread-safe.
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
    // Reference window [proc_lo, proc_hi], indexed relative to proc_lo (clamped to what the
    // FASTA actually returned, so a short contig tail reads as 'N' like before).
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
        let frac = top_count as f64 / depth as f64;
        if frac < params.min_allele_fraction {
            continue;
        }
        if is_paralogous(&c, depth, params) {
            continue; // bi-allelic at a haploid site — paralog/mismapping, not a private call
        }
        let ref_base = ref_chunk.get(r).copied().unwrap_or(b'N');
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
        });
    }
    Ok(out)
}

/// Nominal base quality for the de-novo **diploid** genotype likelihood. The chunked pileup keeps
/// only A/C/G/T counts (per-base quals would blow up WGS memory), and every counted base already
/// cleared `min_base_quality`, so the GL is evaluated at this representative phred. The resulting
/// genotype (0/1 vs 1/1 vs 0/0) is robust to the exact value; PL/GQ are approximate (the per-site
/// [`genotype_sites`] path keeps true per-read quals when exact likelihoods matter).
const DENOVO_DIPLOID_Q: u8 = 30;
/// Minimum reads supporting the alt allele before a site is even considered a candidate variant —
/// suppresses singleton sequencing-error "hets".
const DENOVO_MIN_ALT_READS: u32 = 2;

/// Whole-contig **de-novo diploid** SNV calling: the same chunked, parallel pileup as
/// [`call_denovo`], but each variant site is genotyped at ploidy 2 via the genotype-likelihood model
/// ([`genotype::call_genotype`]) — emitting heterozygous (0/1) and homozygous-alt (1/1) calls, not
/// just a haploid consensus. Biallelic (REF + the top non-REF base) for v1; indels are not called
/// here. Output is in ascending position order, as [`SiteGenotype`] (ploidy 2) — feed it to
/// [`crate::vcf::write_diploid_vcf`].
pub fn call_denovo_diploid(
    bam_path: &Path,
    reference_path: &Path,
    contig: &str,
    params: &HaploidCallerParams,
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

    // Decode-safe worker stack (CRAM 3.1 decode recursion — see [`reader::decode_pool`]).
    let pool = crate::reader::decode_pool(threads)?;
    let nested: Vec<Vec<SiteGenotype>> = pool.install(|| {
        ranges
            .par_iter()
            .map(|&(lo, hi)| denovo_chunk_diploid(bam_path, reference_path, contig, params, &ref_seq, length, lo, hi))
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
        // All non-reference bases clearing the support floor are candidate alts (dominant first).
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

/// Left-align an indel within the reference repeat structure (VCF normalization): an aligner may
/// place an indel anywhere within a homopolymer/STR run, but the canonical representation is the
/// leftmost. Returns the normalized `anchor` (1-based) and allele. A deletion of `len` bases at
/// `[anchor, anchor+len-1]` shifts left while `ref[anchor-1] == ref[anchor+len-1]`; an insertion
/// before `anchor` shifts left while `ref[anchor-1]` equals its last base (rotating the bases).
/// Bounded by `proc_lo` (the loaded reference window start, 1-based).
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

/// Build a diploid indel [`SiteGenotype`] (VCF-style, left-anchored at `emit_pos`) from ref-vs-indel
/// read support, genotyped at ploidy 2 via the sentinel-byte GL (`b'R'` ref-spanning, `b'A'`
/// indel-carrying). `ref_byte` is the reference base at `emit_pos`; `deleted` is the deleted
/// reference bases (empty for an insertion). `None` for a hom-ref / no-call.
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

/// De-novo diploid **indel** calls for this chunk: over each active (indel-evidence) window, extract
/// per-read indel alleles (CIGAR I/D) + ref-spanning support, tally the dominant allele per locus
/// (biallelic v1), and genotype it at ploidy 2. Emits only loci whose VCF position is in the emit
/// range (dedup across chunk boundaries). Left-anchored at the standard VCF convention.
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

            // Multiallelic: one common REF spanning the largest deletion, one ALT per allele.
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

/// Re-fit reads in each indel-active window onto the reference and replace the tally
/// over those windows. Arrays are indexed relative to `region_lo` (1-based).
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

    // ONE indexed query spanning all active windows: decode the region's reads once and route
    // each read to the window(s) it overlaps. The previous code re-queried per window, which on
    // a repeat-rich contig (thousands of indel windows) re-decoded the same CRAM containers over
    // and over — the de-novo hot path. Reads are short, so each overlaps only a window or two,
    // found by binary search over the sorted windows.
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
        // Site 100: run A het (0/1, deep), run B hom-ref (0/0, deep) → real disagreement; depth-
        // weighted vote, both deep, equal weight → tie broken by lower dosage = hom-ref → NOT emitted.
        // Site 200: run A het (deep), run B no-call (depth 1 < min 4) → B abstains, A wins → het.
        // Site 300: both hom-alt (1/1) → hom-alt, depths summed.
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
        // A non-repeat deletion doesn't move: "ACGTC", delete the G (anchor 3).
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
        // One discordant read at low depth — a sequencing error, kept.
        assert!(!para([3, 1, 0, 0])); // second=1 (< 2 reads)
                                      // Scattered errors across other bases, none reaching 2 reads — kept.
        assert!(!para([18, 1, 1, 0])); // second=1
                                       // Genuine bi-allelic pileup (7 derived / 4 ancestral) — paralog, dropped.
        assert!(para([7, 4, 0, 0])); // second=4, 0.36 > 0.20
                                     // Boundary: 2/10 = 0.20 is not strictly above the threshold — kept.
        assert!(!para([8, 2, 0, 0]));
        // 3/10 = 0.30 > 0.20 with 3 reads — dropped.
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
        // Even a 50/50 split is not flagged when the filter is disabled.
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
        };
        let calls = vec![v(2), v(3), v(4)];
        let known: HashSet<i64> = [2, 3].into_iter().collect();
        let private = subtract_known(&calls, &known);
        assert_eq!(private.iter().map(|c| c.position).collect::<Vec<_>>(), vec![4]);
    }
}
