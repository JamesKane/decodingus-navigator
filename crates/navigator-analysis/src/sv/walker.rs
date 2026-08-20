//! The walker over the SV evidence. It is the port of the Scala `SvEvidenceWalker`. It makes one
//! pass over the alignment, and in that pass it collects three things:
//!
//! - the read depth in each bin, for a CNV call;
//! - the discordant read pairs, in the style of BreakDancer;
//! - the split reads from the SA tag, in the style of Pindel.
//!
//! Two walks share one body at the record level, which is [`EvidenceSink::accept_read`]. The two
//! can then never come apart. [`collect_evidence_parallel`] fans one region query for each contig
//! across a rayon pool whose stacks are safe for a decode. [`collect_evidence`] makes one
//! sequential pass, for a file with no coordinate index. Call the parallel one: it falls back to
//! the sequential one by itself.
//!
//! It walks the records as [`AlnRead`] views, and not as one concrete record type. So the BAM path
//! stays on the lazy, zero-copy `bam::Record`. This is a pass over the whole genome, so an
//! owned copy at each read would cost much. The CRAM path gets the decoded record, because it has
//! no cheaper form.

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

/// The tally of the evidence as the walk goes on. One instance covers the whole file in the
/// sequential walk, and one contig in the parallel fan-out. The only difference is which records
/// go into it, and how many entries `depth_bins` holds at the start. So both walks run the same
/// body at each record.
struct EvidenceSink<'a> {
    /// A map from a reference id to an interned name, in header order. Even a walk over one contig
    /// needs the whole map. A discordant pair names the contig of its *mate*, and that is often a
    /// different contig.
    names: &'a [Arc<str>],
    config: &'a SvCallerConfig,
    budget: &'a EvidenceBudget,
    /// The empty name that a pair takes when the reference id of its mate resolves to nothing. The
    /// code interns it, so that path allocates no more than the usual one.
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
        // `names` is a shared slice, so this borrows the interned name for `'a`. It does not
        // borrow from `self`, which leaves the tallies below free to take `&mut self`.
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

    /// Drop the borrowed lookup tables, and keep the tally. The fan-out returns this, and not the
    /// sink itself, so that the result of each contig carries no lifetime.
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

/// An upper limit on the evidence that the code keeps. The contig workers of the fan-out share it,
/// so the limit covers the whole genome, and not one contig. With 24 contigs, and the full limit
/// for each, the limit would hold nothing back.
///
/// It counts what the code *keeps*, and nothing else. The walk still finds the evidence past the
/// limit, and it counts that evidence as dropped. So the totals in the report stay honest, and a
/// run that the limit cut short is visible to the user.
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

    /// Take one slot in `counter`, or return false once the count reaches the limit. It uses a
    /// relaxed ordering. These counters hold the memory down, and nothing orders itself against
    /// them.
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

/// A map from a reference id to a name, in header order, interned once for each walk. Every pair
/// and split read that the code keeps holds an `Arc` from this table. A contig name goes into
/// memory once for the whole run, and not once at each record.
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

/// Collect the SV evidence one contig at a time, in parallel. `contig_lengths` says which contigs
/// get depth bins, and how large those contigs are. `expected_insert_size` and `insert_size_sd`
/// come from read-metrics. A CRAM needs `reference`, and a BAM ignores it. The SV evidence never
/// reads a reference *base*, but the decode of a CRAM record does.
///
/// SV was the last analysis over the whole genome that still walked the file on one thread. On a
/// 30x WGS CRAM that is hours of decode on one core. It took 2 to 5 h for each sample, against
/// about 55 min for the [`crate::unified`] walk over the contigs, on the same files. A separate
/// region query for each contig spends the same total decode across every core, and not on one.
///
/// The code joins the evidence together in header order. For a file in coordinate order that is
/// exactly the order that the sequential walk gives. And the clusterer sorts by position in any
/// case, so the calls are the same either way. This function falls back to [`collect_evidence`]
/// when there is no `.bai` and no `.crai`.
pub fn collect_evidence_parallel(
    bam_path: &Path,
    reference: Option<&Path>,
    contig_lengths: &BTreeMap<String, i64>,
    expected_insert_size: f64,
    insert_size_sd: f64,
    config: &SvCallerConfig,
    cancel: &crate::cancel::CancelToken,
) -> Result<SvEvidenceCollection, AnalysisError> {
    // A region query on one contig needs a coordinate index. Without one, a sequential pass is
    // the only way to reach the records. So take that pass, and do not fail.
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
    // The workers share this, so the limit covers the genome, and not each contig.
    let budget = EvidenceBudget::new(config.max_evidence_records);

    // One work item for each *header* contig, and not for each `contig_lengths` entry. The depth
    // bins cover the contigs that the caller asked for. But the discordant pairs and the split
    // reads cover the whole genome, because the sequential walk sees every record in the file. So
    // the code must visit every contig.
    //
    // Both walks skip a record with no reference position. The sequential one does that with its
    // `is_unmapped` guard, and this one never queries for such a record. Nothing goes missing when
    // this code does not sweep the unmapped tail.
    let process_contig = |name: &Arc<str>| -> Result<ContigEvidence, AnalysisError> {
        // Stop before the cost of the reader of this contig.
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

    // The CRAM decoder of noodles can recurse deep enough to overflow the default 2 MiB worker
    // stack of rayon. In the sequential walker, the larger stack of the main thread holds the same
    // file. An overflow aborts the whole process, so the workers get a stack that is safe for a
    // decode.
    let pool = reader::decode_pool(crate::unified::analysis_thread_count())?;
    let per_contig: Vec<ContigEvidence> = pool.install(|| {
        names
            .par_iter()
            .map(&process_contig)
            .collect::<Result<Vec<_>, AnalysisError>>()
    })?;

    // The merge. The depth-bin keys of two contigs never meet. A write over the map, which starts
    // at zero, then fills a slot and does not add to one. A contig that the caller asked for, and
    // that holds no record, keeps its zeros.
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

/// Collect the SV evidence in one sequential pass. It is the reference that
/// [`collect_evidence_parallel`] must agree with, and it is the walk for an alignment with no
/// coordinate index. The arguments are the same, and that function documents them.
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
    /// Classify one primary read that has a pair. It builds a [`DiscordantPair`] only when that
    /// read is discordant.
    ///
    /// This runs at almost every read in the file, and all but a small fraction of a percent are
    /// concordant. So the code builds nothing before it knows the answer.
    ///
    /// An earlier version allocated the read name, and cloned the name of the mate contig, at the
    /// start. That was two mallocs at each read, and it threw away almost all of them. It is what
    /// made a 30x WGS walk spend its time in malloc. The pair that the code does build now
    /// allocates nothing: both contigs are `Arc` clones from the interned table.
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
    /// Parse the first alignment of the SA tag into a [`SplitRead`]. The clip length is the total
    /// soft clip and hard clip of the read itself.
    ///
    /// Only a read that carries an `SA` tag reaches this code, so this path is not hot, and the
    /// discordant-pair path is. A walk over the whole genome gives about 0.1M to 1M of these,
    /// against billions of reads. That is why the supplementary contig comes straight off the tag
    /// text, and a lookup table to intern it is not worth the code.
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
