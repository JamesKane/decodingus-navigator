//! Mark PCR/optical duplicates on a coordinate-sorted alignment.
//!
//! Take two fragments that start and end at the same place, on the same strands. They are almost
//! surely copies of one original molecule, and not two independent observations. To count both
//! makes the coverage too high. Worse, it makes a sequencing error look like a variant with strong
//! support. So one member of each such group keeps its flag clear, and the rest get `0x400`.
//!
//! The code removes nothing. `0x400` is advice, and each consumer decides for itself. The coverage
//! walk follows it, and a structural-variant caller may not. To delete a read would take that
//! choice away.
//!
//! ## Short reads alone
//!
//! A long-read library, from HiFi or ONT, usually has no PCR step, and two long reads share their
//! endpoints far less often. So the inference from the same endpoints to the same molecule does
//! not hold. A mark on them would throw away real coverage. That is why
//! [`MarkDupParams::enabled`] exists, and why stage C turns this off for the long-read presets.
//!
//! ## The position before the clip
//!
//! The groups use the 5′ position of each end **before its clip**, and not its alignment start.
//! Two copies of one molecule can carry different soft clips, and a mismatch near the end of one
//! copy is enough to cause that. The alignment start then moves, and the fragment does not. To add
//! the clipped bases back gives the place where the molecule began, and that is what the code must
//! compare.
//!
//! ## Both ends of a template must agree
//!
//! A template with a mark on one end, and none on the other, is real damage. A consumer that
//! filters on `0x400` would then see half of a pair.
//!
//! This code marks each end on its own. That is safe for one reason only: the signature of a
//! template is symmetric. It puts the 5′ position and strand of this end together with those of
//! its mate. Both ends of two duplicate templates land in groups with the same membership. The
//! rule "the first one in coordinate order" then takes the same template at both ends.
//!
//! The sort before this is deterministic for exactly that reason: "first" must mean the same thing
//! on every run. A test holds this property directly, and it does not depend on the argument.

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

/// How far back a candidate duplicate can sit, in bases, and still enter a comparison.
///
/// A clip is the only thing that separates the alignment starts of two copies. So this value must
/// be more than the largest clip that can occur, and one read length covers that with room to
/// spare. It also bounds the memory, because the code holds a signature only while it is inside
/// the window.
const DEFAULT_WINDOW: usize = 1_000;

/// The controls of [`mark_duplicates`].
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

/// What this pass did.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct MarkDupStats {
    /// The count of records that the code read, and wrote. This pass never drops one.
    pub records: u64,
    /// Records flagged `0x400`.
    pub duplicates: u64,
    /// Records not eligible: unmapped, secondary, or supplementary.
    pub ineligible: u64,
    /// True when [`MarkDupParams::enabled`] was false, and the code copied every record through
    /// with no mark.
    pub skipped: bool,
}

/// Mark duplicates in `input`, writing to `output`.
///
/// `input` must be in coordinate order. The groups depend on the copies of a molecule that lie
/// near each other in the file, and that holds only after the sort.
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
                    // Compute this again from nothing. Take an input that already carries a mark,
                    // from a second run or from a vendor. It must not keep an answer that this
                    // pass did not reach.
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
/// It is symmetric across a pair by construction. It carries the 5′ position and strand of this
/// end, *and* those of the mate. Both ends of two duplicate templates then make groups with the
/// same membership.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct Signature {
    reference: usize,
    five_prime: i64,
    reverse: bool,
    mate_reference: i64,
    mate_position: i64,
    mate_reverse: bool,
}

/// The signature of a record that the code may mark, or `None` when it may not mark that
/// record.
///
/// A record with no mapping has no position to compare. A secondary record, and a supplementary
/// one, each describe an alignment of a read whose primary record sits elsewhere. To mark one of
/// those would count a template twice, because its primary record already stands for it.
fn signature(record: &RecordBuf) -> Option<Signature> {
    let flags = record.flags();
    if flags.is_unmapped() || flags.is_secondary() || flags.is_supplementary() {
        return None;
    }
    let reference = record.reference_sequence_id()?;
    let start = record.alignment_start()?.get() as i64;

    let (leading, trailing, span) = clip_and_span(record);
    let reverse = flags.is_reverse_complemented();
    // The 5′ end of the molecule. For a forward read that is the alignment start. For a reverse
    // read it is the alignment end. In both cases the code adds the clipped bases back.
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

/// The clip at the start, the clip at the end, and the reference span, from the CIGAR.
///
/// Both a soft clip (`S`) and a hard clip (`H`) count. A hard clip means that somebody removed
/// bases from the record. But the molecule still began that far back.
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

// ---- the window that slides -----------------------------------------------

/// Signatures seen recently, evicted once the file has moved past them.
///
/// The window bounds this, and the file does not, so a WGS costs the same as an exome.
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

    /// True when this signature already occurred inside the window. The code keeps the first
    /// record of a group. Everything after it that matches is a duplicate.
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

/// Copy `header`, and every record, through unchanged. The code calls this when the mark is off. A
/// caller then gets the same output path either way, and it does not have to ask whether stage C
/// ran.
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
