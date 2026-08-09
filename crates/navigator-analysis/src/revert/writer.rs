//! Emit collated reads as paired FASTQ.
//!
//! The invariant the mapper depends on: `_1.fastq` and `_2.fastq` must stay in lockstep, record
//! for record. A template only reaches those files if it has exactly one R1 and exactly one R2;
//! anything else — an unpaired library, a mate lost to hard clipping, flags that don't say which
//! end a read is — goes to singletons. Silently writing an unmatched read into `_1` would shift
//! every later pair by one and mis-pair the entire rest of the file, so the check is per-template
//! rather than a trailing reconciliation.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use super::collate::Merged;
use super::transform::{Mate, RevertedRead};
use super::RevertStats;
use crate::cancel::CancelToken;
use crate::error::AnalysisError;

/// Output buffer per FASTQ file. FASTQ is small records and many of them, so the write path is
/// syscall-bound without a generous buffer.
const FASTQ_BUFFER: usize = 1024 * 1024;

/// Phred offset for FASTQ's ASCII quality encoding (Sanger / Illumina 1.8+).
const PHRED_OFFSET: u8 = 33;

/// Same cadence as the record loop; the merge is long enough to need its own cancellation point.
const CANCEL_CHECK_INTERVAL: u64 = 4096;

/// Drain `merged` into `_1.fastq` / `_2.fastq` / `_singletons.fastq` under `out_dir`.
///
/// Returns the three paths in that order.
pub fn write_fastq(
    mut merged: Merged,
    out_dir: &Path,
    stats: &mut RevertStats,
    cancel: &CancelToken,
) -> Result<(PathBuf, PathBuf, PathBuf), AnalysisError> {
    let p1 = out_dir.join("reverted_1.fastq");
    let p2 = out_dir.join("reverted_2.fastq");
    let ps = out_dir.join("reverted_singletons.fastq");

    let mut w1 = open(&p1)?;
    let mut w2 = open(&p2)?;
    let mut ws = open(&ps)?;

    let mut group: Vec<RevertedRead> = Vec::new();
    // Reused across every record: the ASCII quality line is the only part that needs building
    // rather than copying, and allocating one per read would dominate the write path.
    let mut qual_scratch: Vec<u8> = Vec::new();
    let mut templates = 0u64;

    while merged.next_group(&mut group)? {
        if templates % CANCEL_CHECK_INTERVAL == 0 {
            cancel.check()?;
        }
        templates += 1;

        match pair_of(&group) {
            Some((one, two)) => {
                write_record(&mut w1, &group[one], &p1, &mut qual_scratch)?;
                write_record(&mut w2, &group[two], &p2, &mut qual_scratch)?;
                stats.pairs += 1;
                stats.reads_emitted += 2;
            }
            None => {
                for read in &group {
                    write_record(&mut ws, read, &ps, &mut qual_scratch)?;
                    stats.singletons += 1;
                    stats.reads_emitted += 1;
                }
            }
        }
    }

    flush(&mut w1, &p1)?;
    flush(&mut w2, &p2)?;
    flush(&mut ws, &ps)?;

    Ok((p1, p2, ps))
}

/// The indices of the R1 and R2 of a complete pair, or `None` if this template is not exactly one
/// of each. Duplicated segment bits (two R1s under one name) are treated as unpaired: which of the
/// two is "the" R1 is unanswerable, and guessing would silently corrupt the pairing.
fn pair_of(group: &[RevertedRead]) -> Option<(usize, usize)> {
    let mut one = None;
    let mut two = None;
    for (i, read) in group.iter().enumerate() {
        match read.mate {
            Mate::One if one.is_none() => one = Some(i),
            Mate::Two if two.is_none() => two = Some(i),
            Mate::Unpaired => {}
            // A second R1 or R2 under the same name.
            _ => return None,
        }
    }
    match (one, two) {
        // Exactly one of each, and nothing else tagging along under the same name.
        (Some(a), Some(b)) if group.len() == 2 => Some((a, b)),
        _ => None,
    }
}

fn open(path: &Path) -> Result<BufWriter<File>, AnalysisError> {
    let file = File::create(path).map_err(|e| AnalysisError::io(path, e))?;
    Ok(BufWriter::with_capacity(FASTQ_BUFFER, file))
}

fn flush(w: &mut BufWriter<File>, path: &Path) -> Result<(), AnalysisError> {
    w.flush().map_err(|e| AnalysisError::io(path, e))
}

/// One FASTQ record. Names are written bare — no `/1` or `/2` — so R1/R2 pair by position; see the
/// module docs on why. Qualities are shifted into ASCII here, the inverse of the decode in
/// [`super::transform`], through `scratch` so the shift costs no allocation per read.
fn write_record(
    w: &mut BufWriter<File>,
    read: &RevertedRead,
    path: &Path,
    scratch: &mut Vec<u8>,
) -> Result<(), AnalysisError> {
    scratch.clear();
    scratch.extend(read.qualities.iter().map(|q| q.saturating_add(PHRED_OFFSET)));

    let write = |w: &mut BufWriter<File>| -> std::io::Result<()> {
        w.write_all(b"@")?;
        w.write_all(&read.name)?;
        w.write_all(b"\n")?;
        w.write_all(&read.sequence)?;
        w.write_all(b"\n+\n")?;
        w.write_all(scratch)?;
        w.write_all(b"\n")
    };
    write(w).map_err(|e| AnalysisError::io(path, e))
}
