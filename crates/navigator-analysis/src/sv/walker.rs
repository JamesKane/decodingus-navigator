//! SV evidence walker — port of the Scala `SvEvidenceWalker`. One pass over the alignment
//! collecting per-bin read depth (CNV), discordant read pairs (BreakDancer-style), and
//! split reads from the SA tag (Pindel-style).
//!
//! Two walks share one per-record body ([`EvidenceSink::accept_read`]), so they can not drift:
//! [`collect_evidence_parallel`] fans one region query per contig across a decode-safe rayon pool,
//! and [`collect_evidence`] makes a single sequential pass for files with no coordinate index.
//! Prefer the parallel entry point — it falls back to the sequential one on its own.
//!
//! It walks records as [`AlnRead`] views rather than a concrete record type: the BAM path stays on
//! the lazy, zero-copy `bam::Record` (this is a whole-genome pass, so a per-read owned copy would
//! be costly), while the CRAM path gets the decoded record it has no cheaper form of.

use std::collections::BTreeMap;
use std::path::Path;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::Arc;

use noodles::core::Region;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::data::field::Tag;
use rayon::prelude::*;

use super::evidence::{DiscordantPair, DiscordantReason, SplitRead, SvEvidenceCollection};
use super::types::SvCallerConfig;
use crate::error::AnalysisError;
use crate::reader::{self, RecordSink};
use crate::readview::AlnRead;

const SA_TAG: Tag = Tag::new(b'S', b'A');

/// How often the record loop polls the cancel token, matching the indexed reader's own cadence.
const CANCEL_CHECK_RECORDS: u32 = 4096;

/// The running evidence tally. One instance covers the whole file in the sequential walk and one
/// contig in the parallel fan-out — the difference is only which records get fed to it and how
/// many entries `depth_bins` starts with, so both walks run the identical per-record body.
struct EvidenceSink<'a> {
    /// Reference id -> interned name, header order. Needed whole even per contig: a discordant pair
    /// names its *mate's* contig, which is routinely a different one.
    names: &'a [Arc<str>],
    config: &'a SvCallerConfig,
    budget: &'a EvidenceBudget,
    /// The empty name a pair falls back to when its mate's reference id resolves to nothing,
    /// interned so that path allocates no more than the normal one.
    unknown_contig: Arc<str>,
    insert_min: f64,
    insert_max: f64,
    depth_bins: BTreeMap<String, Vec<u32>>,
    discordant_pairs: Vec<DiscordantPair>,
    split_reads: Vec<SplitRead>,
    discordant_dropped: u64,
    split_dropped: u64,
}

impl<'a> EvidenceSink<'a> {
    fn new(
        names: &'a [Arc<str>],
        config: &'a SvCallerConfig,
        budget: &'a EvidenceBudget,
        insert_min: f64,
        insert_max: f64,
        depth_bins: BTreeMap<String, Vec<u32>>,
    ) -> Self {
        Self {
            names,
            config,
            budget,
            unknown_contig: Arc::from(""),
            insert_min,
            insert_max,
            depth_bins,
            discordant_pairs: Vec::new(),
            split_reads: Vec::new(),
            discordant_dropped: 0,
            split_dropped: 0,
        }
    }

    fn accept_read(&mut self, record: &impl AlnRead) {
        let config = self.config;
        let flags = record.flags();
        if flags.is_unmapped() {
            return;
        }
        let Some(ref_id) = record.reference_sequence_id() else {
            return;
        };
        // `names` is a shared slice, so this borrows the interned name for `'a` — not from `self`,
        // which leaves the tallies below free to take `&mut self`.
        let Some(contig) = self.names.get(ref_id) else { return };
        let Some(start) = record.alignment_start().map(|p| p as i64) else {
            return;
        };
        let mapq = record.mapping_quality().unwrap_or(255);
        let secondary_or_supp = flags.is_secondary() || flags.is_supplementary();

        // 1. Depth (primary, non-supplementary only).
        if !secondary_or_supp {
            if let Some(bins) = self.depth_bins.get_mut(&**contig) {
                let bin = (start / config.bin_size) as usize;
                if bin < bins.len() {
                    bins[bin] += 1;
                }
            }
        }

        // 2. Discordant pairs (primary, paired only).
        if !secondary_or_supp && flags.is_segmented() {
            if let Some(dp) = self.detect_discordant_pair(record, contig, mapq) {
                if self.budget.claim_discordant() {
                    self.discordant_pairs.push(dp);
                } else {
                    self.discordant_dropped += 1;
                }
            }
        }

        // 3. Split reads (SA tag).
        if mapq >= config.min_mapq {
            if let Some(sr) = self.extract_split_read(record, contig, mapq) {
                if self.budget.claim_split() {
                    self.split_reads.push(sr);
                } else {
                    self.split_dropped += 1;
                }
            }
        }
    }

    /// Drop the borrowed lookup tables, keeping the tally. The fan-out returns this rather than the
    /// sink itself so the per-contig results carry no lifetime.
    fn into_parts(self) -> ContigEvidence {
        ContigEvidence {
            depth_bins: self.depth_bins,
            discordant_pairs: self.discordant_pairs,
            split_reads: self.split_reads,
            discordant_dropped: self.discordant_dropped,
            split_dropped: self.split_dropped,
        }
    }
}

/// A ceiling on retained evidence, shared across the fan-out's contig workers so the bound is
/// genome-wide rather than per contig (24 contigs each allowed the full cap would bound nothing).
///
/// Counts only what is *kept*. Evidence past the cap is still detected and counted as dropped, so
/// the reported totals stay honest and a truncated run is visible rather than silent.
struct EvidenceBudget {
    cap: u64,
    discordant_kept: AtomicU64,
    split_kept: AtomicU64,
}

impl EvidenceBudget {
    fn new(cap: u64) -> Self {
        Self {
            cap,
            discordant_kept: AtomicU64::new(0),
            split_kept: AtomicU64::new(0),
        }
    }

    /// Claim one slot in `counter`, or return false once the cap is met. Relaxed ordering: the
    /// counters guard memory growth, and nothing is ordered against them.
    fn claim(&self, counter: &AtomicU64) -> bool {
        counter
            .fetch_update(Ordering::Relaxed, Ordering::Relaxed, |n| {
                (n < self.cap).then_some(n + 1)
            })
            .is_ok()
    }

    fn claim_discordant(&self) -> bool {
        self.claim(&self.discordant_kept)
    }

    fn claim_split(&self) -> bool {
        self.claim(&self.split_kept)
    }
}

/// An evidence tally with no borrows: one contig's slice of it in the parallel walk, the whole
/// file's in the sequential one.
#[derive(Default)]
struct ContigEvidence {
    depth_bins: BTreeMap<String, Vec<u32>>,
    discordant_pairs: Vec<DiscordantPair>,
    split_reads: Vec<SplitRead>,
    discordant_dropped: u64,
    split_dropped: u64,
}

impl ContigEvidence {
    fn into_collection(self, expected_insert_size: f64, insert_size_sd: f64, cap: u64) -> SvEvidenceCollection {
        // Loud, because the alternative reading of a short call set is "this sample is clean".
        if self.discordant_dropped > 0 || self.split_dropped > 0 {
            eprintln!(
                "warning: SV evidence hit the {cap}-record cap — dropped {} discordant pair(s) and \
                 {} split read(s). Calls in the truncated regions may be missing; raise \
                 SvCallerConfig::max_evidence_records to keep them.",
                self.discordant_dropped, self.split_dropped
            );
        }
        SvEvidenceCollection {
            discordant_pairs: self.discordant_pairs,
            split_reads: self.split_reads,
            depth_bins: self.depth_bins,
            sample_name: "unknown".to_string(), // RG SM only feeds deferred VCF/summary output
            expected_insert_size,
            insert_size_sd,
            discordant_pairs_dropped: self.discordant_dropped,
            split_reads_dropped: self.split_dropped,
        }
    }
}

impl RecordSink for EvidenceSink<'_> {
    fn accept(&mut self, record: &impl AlnRead) {
        self.accept_read(record);
    }
}

/// Reference id -> name (header order), interned once per walk. Every retained pair and split read
/// holds `Arc`s from this table, so a contig name is stored once for the whole run instead of once
/// per record.
fn contig_names(header: &noodles::sam::Header) -> Vec<Arc<str>> {
    header
        .reference_sequences()
        .keys()
        .map(|n| Arc::from(String::from_utf8_lossy(n.as_ref()).as_ref()))
        .collect()
}

/// A zeroed depth-bin vector for a contig of `len` bases.
fn zeroed_bins(len: i64, bin_size: i64) -> Vec<u32> {
    let num_bins = ((len + bin_size - 1) / bin_size).max(0) as usize;
    vec![0u32; num_bins]
}

/// Zeroed depth bins for every requested contig. Contigs absent from the alignment header keep
/// their zeroed vector, so the segmenter sees the same key set either way.
fn all_zeroed_bins(contig_lengths: &BTreeMap<String, i64>, bin_size: i64) -> BTreeMap<String, Vec<u32>> {
    contig_lengths
        .iter()
        .map(|(c, &len)| (c.clone(), zeroed_bins(len, bin_size)))
        .collect()
}

/// The insert-size window outside which a pair counts as an outlier.
fn insert_bounds(expected_insert_size: f64, insert_size_sd: f64, config: &SvCallerConfig) -> (f64, f64) {
    let max = expected_insert_size + config.insert_size_z_threshold * insert_size_sd;
    let min = (expected_insert_size - config.insert_size_z_threshold * insert_size_sd).max(0.0);
    (min, max)
}

/// Collect SV evidence one contig at a time, in parallel. `contig_lengths` selects which contigs
/// get depth bins (and their sizes); `expected_insert_size`/`insert_size_sd` come from
/// read-metrics. `reference` is required for CRAM (ignored for BAM) — SV evidence never consults
/// reference *bases*, but decoding a CRAM record at all does.
///
/// SV was the last whole-genome analysis still walking the file on one thread, which on a 30x WGS
/// CRAM is hours of single-core decode: 2–5 h per sample, against ~55 min for the per-contig
/// [`crate::unified`] walk over the same files. Region-querying each contig separately spends the
/// same total decode across every core instead of one.
///
/// Evidence is concatenated in header order, which for a coordinate-sorted file is exactly the
/// order the sequential walk emits — and the clusterer sorts by position regardless, so the calls
/// are identical either way. Falls back to [`collect_evidence`] when there is no `.bai`/`.crai`.
pub fn collect_evidence_parallel(
    bam_path: &Path,
    reference: Option<&Path>,
    contig_lengths: &BTreeMap<String, i64>,
    expected_insert_size: f64,
    insert_size_sd: f64,
    config: &SvCallerConfig,
    cancel: &crate::cancel::CancelToken,
) -> Result<SvEvidenceCollection, AnalysisError> {
    // Per-contig region queries need a coordinate index. Without one the only way to reach the
    // records is a sequential pass, so take it rather than failing.
    if !reader::has_region_index(bam_path) {
        return collect_evidence(
            bam_path,
            reference,
            contig_lengths,
            expected_insert_size,
            insert_size_sd,
            config,
            cancel,
        );
    }

    let header = reader::read_header(bam_path, reference)?;
    let names = contig_names(&header);
    let (insert_min, insert_max) = insert_bounds(expected_insert_size, insert_size_sd, config);
    // Shared across the workers so the ceiling is on the genome, not on each contig.
    let budget = EvidenceBudget::new(config.max_evidence_records);

    // One work item per *header* contig, not per `contig_lengths` entry: depth bins are limited to
    // the requested contigs, but discordant pairs and split reads are collected genome-wide (the
    // sequential walk sees every record in the file), so every contig has to be visited.
    //
    // Records with no reference position are skipped by both walks — the sequential one via the
    // `is_unmapped` guard, this one by never querying for them — so nothing is lost by not sweeping
    // the unmapped tail here.
    let process_contig = |name: &Arc<str>| -> Result<ContigEvidence, AnalysisError> {
        // Bail before paying for this contig's reader.
        cancel.check()?;
        let (h, mut idx) = reader::open_indexed(bam_path, reference)?;
        let region = Region::new(name.as_bytes().to_vec(), ..); // whole contig
        let bins = match contig_lengths.get(&**name) {
            Some(&len) => BTreeMap::from([(name.to_string(), zeroed_bins(len, config.bin_size))]),
            None => BTreeMap::new(),
        };
        let mut sink = EvidenceSink::new(&names, config, &budget, insert_min, insert_max, bins);
        idx.for_each(&h, &region, &mut sink, cancel)?;
        Ok(sink.into_parts())
    };

    // noodles' CRAM decoder can recurse deeply enough to blow rayon's default 2 MiB worker stack
    // (the main thread's larger stack handles the same file in the sequential walker) — and an
    // overflow aborts the whole process, so the workers get a decode-safe stack.
    let pool = reader::decode_pool(crate::unified::analysis_thread_count())?;
    let per_contig: Vec<ContigEvidence> = pool.install(|| {
        names
            .par_iter()
            .map(&process_contig)
            .collect::<Result<Vec<_>, AnalysisError>>()
    })?;

    // Merge. Each contig owns a disjoint depth-bin key, so inserting over the pre-zeroed map is a
    // fill rather than a sum; a requested contig with no records keeps its zeros.
    let mut merged = ContigEvidence {
        depth_bins: all_zeroed_bins(contig_lengths, config.bin_size),
        ..ContigEvidence::default()
    };
    for part in per_contig {
        merged.depth_bins.extend(part.depth_bins);
        merged.discordant_pairs.extend(part.discordant_pairs);
        merged.split_reads.extend(part.split_reads);
        merged.discordant_dropped += part.discordant_dropped;
        merged.split_dropped += part.split_dropped;
    }
    Ok(merged.into_collection(expected_insert_size, insert_size_sd, config.max_evidence_records))
}

/// Collect SV evidence in a single sequential pass — the parity reference for
/// [`collect_evidence_parallel`], and the walk used when the alignment has no coordinate index.
/// Arguments are as documented there.
pub fn collect_evidence(
    bam_path: &Path,
    reference: Option<&Path>,
    contig_lengths: &BTreeMap<String, i64>,
    expected_insert_size: f64,
    insert_size_sd: f64,
    config: &SvCallerConfig,
    cancel: &crate::cancel::CancelToken,
) -> Result<SvEvidenceCollection, AnalysisError> {
    let (header, mut reader) = reader::open_seq(bam_path, reference)?;
    let names = contig_names(&header);
    let (insert_min, insert_max) = insert_bounds(expected_insert_size, insert_size_sd, config);
    let budget = EvidenceBudget::new(config.max_evidence_records);
    let mut sink = EvidenceSink::new(
        &names,
        config,
        &budget,
        insert_min,
        insert_max,
        all_zeroed_bins(contig_lengths, config.bin_size),
    );

    let mut seen = 0u32;
    for result in reader.records_lazy(&header) {
        seen += 1;
        if seen % CANCEL_CHECK_RECORDS == 0 {
            cancel.check()?;
        }
        sink.accept_read(&result?);
    }

    Ok(sink
        .into_parts()
        .into_collection(expected_insert_size, insert_size_sd, config.max_evidence_records))
}

impl EvidenceSink<'_> {
    /// Classify one primary paired read, building a [`DiscordantPair`] only if it is discordant.
    ///
    /// This runs for essentially every read in the file and all but a fraction of a percent are
    /// concordant, so nothing is built before the verdict is known. It used to allocate the read
    /// name and clone the mate contig name up front — two mallocs per read, thrown away almost
    /// every time — which is what made a 30x WGS walk malloc-bound. The pair that does get built is
    /// now allocation-free: both contigs are `Arc` clones from the interned table.
    fn detect_discordant_pair(&self, record: &impl AlnRead, contig: &Arc<str>, mapq: u8) -> Option<DiscordantPair> {
        let flags = record.flags();
        if flags.is_mate_unmapped() || mapq < self.config.min_mapq {
            return None;
        }
        let ref_id = record.reference_sequence_id();
        let mate_ref_id = record.mate_reference_sequence_id();
        let pos1 = record.alignment_start().map_or(0, |p| p as i64);
        let mate_pos = record.mate_alignment_start().map_or(0, |p| p as i64);

        // Same precedence as before: inter-chromosomal, then insert size, then orientation.
        let template_len = record.template_length().abs();
        let (insert_size, reason) = if ref_id != mate_ref_id {
            (0, DiscordantReason::InterChromosomal)
        } else if template_len as f64 > self.insert_max || (template_len > 0 && (template_len as f64) < self.insert_min)
        {
            (template_len, DiscordantReason::InsertSizeOutlier)
        } else if !is_expected_orientation(record, pos1, mate_pos) {
            (template_len, DiscordantReason::WrongOrientation)
        } else {
            return None;
        };

        Some(DiscordantPair {
            chrom1: Arc::clone(contig),
            pos1,
            strand1: if flags.is_reverse_complemented() { '-' } else { '+' },
            // An unresolvable mate contig kept the empty string before; `unknown_contig` is the
            // same value interned once, so the fallback allocates nothing either.
            chrom2: mate_ref_id
                .and_then(|i| self.names.get(i))
                .unwrap_or(&self.unknown_contig)
                .clone(),
            pos2: mate_pos,
            strand2: if flags.is_mate_reverse_complemented() { '-' } else { '+' },
            insert_size,
            mapq,
            reason,
        })
    }
}

/// Standard Illumina FR orientation check (mirrors the Scala logic).
fn is_expected_orientation(record: &impl AlnRead, pos1: i64, mate_pos: i64) -> bool {
    let flags = record.flags();
    let read_neg = flags.is_reverse_complemented();
    let mate_neg = flags.is_mate_reverse_complemented();
    if read_neg == mate_neg {
        false // tandem
    } else if pos1 < mate_pos {
        !read_neg && mate_neg // read upstream -> read on +
    } else {
        read_neg && !mate_neg // read downstream -> read on -
    }
}

impl EvidenceSink<'_> {
    /// Parse the first SA-tag alignment into a [`SplitRead`]; clip length is the read's own
    /// soft/hard-clip total. Only reads that actually carry an `SA` tag get this far, so unlike
    /// the discordant-pair path this one is not hot — a genome-wide walk yields ~0.1–1 M of these
    /// against billions of reads, which is why interning the supplementary contig off the tag text
    /// is not worth a lookup table.
    fn extract_split_read(&self, record: &impl AlnRead, contig: &Arc<str>, mapq: u8) -> Option<SplitRead> {
        let sa = record.string_tag(SA_TAG)?;
        if sa.is_empty() {
            return None;
        }
        let first = sa.split(';').next().unwrap_or("");
        let parts: Vec<&str> = first.split(',').collect();
        if parts.len() < 5 {
            return None;
        }
        let (Ok(supp_pos), Some(supp_strand), Ok(supp_mapq)) =
            (parts[1].parse::<i64>(), parts[2].chars().next(), parts[4].parse::<u8>())
        else {
            return None;
        };

        // Clip length from this read's CIGAR (sum of S/H ops).
        let clip_length: i32 = record.cigar_with(|ops| {
            ops.filter(|(kind, _)| matches!(kind, Kind::SoftClip | Kind::HardClip))
                .map(|(_, len)| len as i32)
                .sum()
        });

        if supp_mapq >= self.config.min_mapq && clip_length >= 10 {
            Some(SplitRead {
                primary_chrom: Arc::clone(contig),
                primary_pos: record.alignment_start().map_or(0, |p| p as i64),
                primary_strand: if record.flags().is_reverse_complemented() {
                    '-'
                } else {
                    '+'
                },
                supp_chrom: Arc::from(parts[0]),
                supp_pos,
                supp_strand,
                clip_length,
                mapq: mapq.min(supp_mapq),
            })
        } else {
            None
        }
    }
}
