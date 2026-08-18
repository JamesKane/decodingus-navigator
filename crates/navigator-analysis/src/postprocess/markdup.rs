//! Mark PCR/optical duplicates on a coordinate-sorted alignment.
//!
//! Two fragments that start and end at the same place, on the same strands, are almost certainly
//! copies of one original molecule rather than two independent observations. Counting both inflates
//! coverage and, worse, makes a sequencing error present as a confidently-supported variant. So one
//! member of each such group keeps its flag clear and the rest get `0x400`.
//!
//! Nothing is removed. `0x400` is advice, and every consumer decides for itself — the coverage walk
//! honours it, a structural-variant caller may not. Deleting reads would take that choice away.
//!
//! ## Short reads only
//!
//! Long-read libraries (HiFi, ONT) are typically PCR-free, and long reads genuinely share endpoints
//! far less often, so the inference "same endpoints therefore same molecule" does not hold. Marking
//! them would discard real coverage, which is why [`MarkDupParams::enabled`] exists and why stage C
//! turns this off for the long-read presets.
//!
//! ## Unclipped positions
//!
//! Grouping uses each end's **unclipped** 5′ position, not its alignment start. Two copies of one
//! molecule can be soft-clipped differently — a mismatch near one copy's end is enough — which
//! moves the alignment start without moving the fragment. Adding the clipped bases back recovers
//! where the molecule actually began, which is the thing being compared.
//!
//! ## Both ends of a template must agree
//!
//! A template with one end marked and the other not is a real corruption: consumers that filter on
//! `0x400` would see half a pair. This marks each end independently, which is safe only because a
//! template's signature is symmetric — it combines this end's 5′ position and strand with its
//! mate's, so both ends of two duplicate templates land in groups with identical membership, and
//! first-seen-in-coordinate-order picks the same template at both. The sort upstream is
//! deterministic precisely so that "first" means the same thing on every run. There is a test that
//! holds this property directly rather than trusting the argument.

use std::collections::{HashMap, VecDeque};
use std::path::Path;

use noodles::sam;
use noodles::sam::alignment::io::Write as _;
use noodles::sam::alignment::record::cigar::op::Kind;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::RecordBuf;

use super::bamio;
use crate::cancel::CancelToken;
use crate::error::AnalysisError;

const CANCEL_CHECK_INTERVAL: u64 = 4096;

/// How far back a candidate duplicate can sit and still be compared, in bases.
///
/// Only clipping separates two copies' alignment starts, so this needs to exceed the largest
/// plausible clip — a read length, comfortably. It also bounds memory: only signatures within the
/// window are held.
const DEFAULT_WINDOW: usize = 1_000;

/// Tuning for [`mark_duplicates`].
#[derive(Debug, Clone)]
pub struct MarkDupParams {
    /// Off for long-read data. See the module docs.
    pub enabled: bool,
    /// Lookback in bases; see [`DEFAULT_WINDOW`].
    pub window: usize,
}

impl Default for MarkDupParams {
    fn default() -> Self {
        Self {
            enabled: true,
            window: DEFAULT_WINDOW,
        }
    }
}

/// What the marking did.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MarkDupStats {
    /// Records read, and written — marking never drops one.
    pub records: u64,
    /// Records flagged `0x400`.
    pub duplicates: u64,
    /// Records not eligible: unmapped, secondary, or supplementary.
    pub ineligible: u64,
    /// True when [`MarkDupParams::enabled`] was false and records were copied through unmarked.
    pub skipped: bool,
}

/// Mark duplicates in `input`, writing to `output`.
///
/// `input` must be coordinate-sorted — grouping relies on copies of a molecule being near each
/// other in the file, which is only true after the sort.
pub fn mark_duplicates(
    input: &Path,
    output: &Path,
    params: &MarkDupParams,
    cancel: &CancelToken,
    progress: &mut dyn FnMut(u64),
) -> Result<MarkDupStats, AnalysisError> {
    if let Some(parent) = output.parent() {
        std::fs::create_dir_all(parent).map_err(|e| AnalysisError::io(parent, e))?;
    }

    let mut reader = bamio::open(input)?;
    let header = reader.read_header().map_err(|e| AnalysisError::io(input, e))?;

    let mut writer = bamio::create(output)?;
    writer.write_header(&header).map_err(|e| AnalysisError::io(output, e))?;

    let mut stats = MarkDupStats {
        skipped: !params.enabled,
        ..Default::default()
    };
    let mut seen = SeenSignatures::new(params.window);

    for (i, result) in reader.record_bufs(&header).enumerate() {
        if i as u64 % CANCEL_CHECK_INTERVAL == 0 {
            cancel.check()?;
            progress(stats.records);
        }
        let mut record = result.map_err(|e| AnalysisError::io(input, e))?;
        stats.records += 1;

        if params.enabled {
            match signature(&record) {
                Some(sig) => {
                    // Recompute from scratch: an input that was marked before (a re-run, or a
                    // vendor file) must not inherit a verdict this pass did not reach.
                    set_duplicate(&mut record, false);
                    if seen.is_duplicate(sig) {
                        set_duplicate(&mut record, true);
                        stats.duplicates += 1;
                    }
                }
                None => stats.ineligible += 1,
            }
        }

        writer
            .write_alignment_record(&header, &record)
            .map_err(|e| AnalysisError::io(output, e))?;
    }

    bamio::finish(writer, output)?;
    progress(stats.records);
    Ok(stats)
}

// ---- signatures -----------------------------------------------------------

/// What makes two records copies of one molecule.
///
/// Symmetric across a pair by construction: it carries this end's 5′ position and strand *and* the
/// mate's, so both ends of two duplicate templates produce groups with the same membership.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct Signature {
    reference: usize,
    five_prime: i64,
    reverse: bool,
    mate_reference: i64,
    mate_position: i64,
    mate_reverse: bool,
}

/// The signature of an eligible record, or `None` when the record can not be marked.
///
/// Unmapped records have no position to compare. Secondary and supplementary records describe an
/// alignment of a read whose primary is elsewhere; marking them would double-count a template that
/// the primary already represents.
fn signature(record: &RecordBuf) -> Option<Signature> {
    let flags = record.flags();
    if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
        return None;
    }
    let reference = record.reference_sequence_id()?;
    let start = record.alignment_start()?.get() as i64;

    let (leading, trailing, span) = clip_and_span(record);
    let reverse = flags.is_reverse_complemented();
    // The 5′ end of the molecule: the alignment start for a forward read, the alignment end for a
    // reverse one — in both cases with the clipped bases added back.
    let five_prime = if reverse {
        start + span - 1 + trailing
    } else {
        start - leading
    };

    let (mate_reference, mate_position, mate_reverse) = if flags.is_segmented() {
        (
            record.mate_reference_sequence_id().map(|r| r as i64).unwrap_or(-1),
            record.mate_alignment_start().map(|p| p.get() as i64).unwrap_or(-1),
            flags.is_mate_reverse_complemented(),
        )
    } else {
        // Unpaired: nothing to combine with, so the read's own end is the whole signature.
        (-1, -1, false)
    };

    Some(Signature {
        reference,
        five_prime,
        reverse,
        mate_reference,
        mate_position,
        mate_reverse,
    })
}

/// Leading clip, trailing clip, and reference span, from the CIGAR.
///
/// Both soft (`S`) and hard (`H`) clips count: a hard clip means bases were removed from the
/// record, but the molecule still started that far back.
fn clip_and_span(record: &RecordBuf) -> (i64, i64, i64) {
    use noodles::sam::alignment::record::Cigar as _;

    let cigar = record.cigar();
    let ops: Vec<(Kind, usize)> = cigar
        .iter()
        .filter_map(|op| op.ok().map(|o| (o.kind(), o.len())))
        .collect();

    let is_clip = |k: Kind| matches!(k, Kind::SoftClip | Kind::HardClip);
    let leading: i64 = ops
        .iter()
        .take_while(|(k, _)| is_clip(*k))
        .map(|(_, n)| *n as i64)
        .sum();
    let trailing: i64 = ops
        .iter()
        .rev()
        .take_while(|(k, _)| is_clip(*k))
        .map(|(_, n)| *n as i64)
        .sum();
    let span: i64 = ops
        .iter()
        .filter(|(k, _)| {
            matches!(
                k,
                Kind::Match | Kind::Deletion | Kind::Skip | Kind::SequenceMatch | Kind::SequenceMismatch
            )
        })
        .map(|(_, n)| *n as i64)
        .sum();

    (leading, trailing, span.max(1))
}

fn set_duplicate(record: &mut RecordBuf, duplicate: bool) {
    let mut flags = record.flags();
    flags.set(Flags::DUPLICATE, duplicate);
    *record.flags_mut() = flags;
}

// ---- the sliding window ---------------------------------------------------

/// Signatures seen recently, evicted once the file has moved past them.
///
/// Bounded by the window rather than the file, so a WGS costs the same as an exome.
struct SeenSignatures {
    window: i64,
    live: HashMap<Signature, ()>,
    order: VecDeque<(i64, Signature)>,
}

impl SeenSignatures {
    fn new(window: usize) -> Self {
        Self {
            window: window.max(1) as i64,
            live: HashMap::new(),
            order: VecDeque::new(),
        }
    }

    /// Whether this signature has been seen inside the window. The first record of a group is the
    /// one kept; everything matching it afterwards is a duplicate.
    fn is_duplicate(&mut self, sig: Signature) -> bool {
        self.evict(sig.five_prime);
        if self.live.contains_key(&sig) {
            return true;
        }
        self.live.insert(sig, ());
        self.order.push_back((sig.five_prime, sig));
        false
    }

    /// Drop signatures the file has moved past. Anything further back than the window can not be a
    /// duplicate of the current record, because only clipping separates two copies' positions.
    fn evict(&mut self, position: i64) {
        while let Some((pos, sig)) = self.order.front().copied() {
            if position - pos <= self.window {
                break;
            }
            self.order.pop_front();
            self.live.remove(&sig);
        }
    }
}

/// Copy `header` and every record through unchanged. Used when marking is disabled, so callers get
/// the same output path either way rather than branching on whether stage C ran.
pub fn copy_through(input: &Path, output: &Path) -> Result<u64, AnalysisError> {
    let mut reader = bamio::open(input)?;
    let header: sam::Header = reader.read_header().map_err(|e| AnalysisError::io(input, e))?;

    let mut writer = bamio::create(output)?;
    writer.write_header(&header).map_err(|e| AnalysisError::io(output, e))?;

    let mut n = 0;
    for result in reader.record_bufs(&header) {
        let record = result.map_err(|e| AnalysisError::io(input, e))?;
        writer
            .write_alignment_record(&header, &record)
            .map_err(|e| AnalysisError::io(output, e))?;
        n += 1;
    }
    bamio::finish(writer, output)?;
    Ok(n)
}
