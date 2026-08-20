//! Emit collated reads as paired FASTQ.
//!
//! Here is the invariant that the mapper depends on. `_1.fastq` and `_2.fastq` must stay in step,
//! record for record. A template reaches those two files only when it holds exactly one R1 and
//! exactly one R2.
//!
//! Everything else goes to the singletons file. That covers a library with no pairs, a mate that
//! a hard clip removed, and flags that do not say which end a read is.
//!
//! To write a read with no partner into `_1` would move every later pair by one. The whole rest of
//! the file would then hold the wrong pairs, and nobody would see it happen. So the check runs at
//! each template, and not as a reconciliation at the end.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use flate2::write::GzEncoder;
use flate2::Compression;

use super::collate::Merged;
use super::transform::{Mate, RevertedRead};
use super::RevertStats;
use crate::cancel::CancelToken;
use crate::error::AnalysisError;

/// The size of the output buffer of each FASTQ file. FASTQ holds many small records, so without a
/// large buffer the syscalls control the time of the write path.
const FASTQ_BUFFER: usize = 1024 * 1024;

/// The reverted reads go out **through gzip**, and that is not an improvement for its own sake.
///
/// A 30x WGS reverted to plain FASTQ is about 200 GB. That is more than the free space on a usual
/// machine, and some times the size of the alignment that it came from. FASTQ stores one ASCII
/// quality byte at each base, and a BAM packs its data where FASTQ does not. Through gzip the same
/// data is nearer to 55 GB. That is the difference between a pipeline that runs and one that fills
/// the disk in its second stage.
///
/// The level is `fast`, and not `default`. This is a temporary file, and the mapper reads it once.
/// A few percent of the ratio buys much less CPU, in a job that already takes hours. That is the
/// correct side of the curve. The mapper finds gzip from the first bytes of the file, so
/// nothing after this changes.
const FASTQ_COMPRESSION: Compression = Compression::fast();

/// A FASTQ sink that puts its output through gzip.
type FastqWriter = GzEncoder<BufWriter<navigator_resource::PacedFile>>;

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
    let p1 = out_dir.join("reverted_1.fastq.gz");
    let p2 = out_dir.join("reverted_2.fastq.gz");
    let ps = out_dir.join("reverted_singletons.fastq.gz");

    let mut w1 = open(&p1)?;
    let mut w2 = open(&p2)?;
    let mut ws = open(&ps)?;

    let mut group: Vec<RevertedRead> = Vec::new();
    // Every record uses this again. The ASCII quality line is the one part that the code must
    // build, and not copy. To allocate one at each read would control the time of the write
    // path.
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

    finish(w1, &p1)?;
    finish(w2, &p2)?;
    finish(ws, &ps)?;

    Ok((p1, p2, ps))
}

/// The indices of the R1 and the R2 of a complete pair, or `None` when this template does not hold
/// exactly one of each.
///
/// A template with the same segment bit twice, which is two R1 records under one name, counts as
/// unpaired. Nobody can say which of the two is "the" R1. A guess would put the wrong reads
/// together, and nobody would see it happen.
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
        // Exactly one of each, and nothing else under the same name.
        (Some(a), Some(b)) if group.len() == 2 => Some((a, b)),
        _ => None,
    }
}

fn open(path: &Path) -> Result<FastqWriter, AnalysisError> {
    let file = File::create(path).map_err(|e| AnalysisError::io(path, e))?;
    Ok(GzEncoder::new(
        BufWriter::with_capacity(FASTQ_BUFFER, navigator_resource::PacedFile::new(file)),
        FASTQ_COMPRESSION,
    ))
}

/// Finish the gzip stream. Do not only flush it. A gzip member with no trailer is a file that
/// stops early. A reader would then stop early too, and it would give no error. That is the same
/// class of bug as a missing BGZF EOF block, and it is as hard to see.
fn finish(w: FastqWriter, path: &Path) -> Result<(), AnalysisError> {
    let mut inner = w.finish().map_err(|e| AnalysisError::io(path, e))?;
    inner.flush().map_err(|e| AnalysisError::io(path, e))
}

/// One FASTQ record. The names go out bare, with no `/1` and no `/2`, so an R1 and an R2 pair by
/// their position. The module documentation says why. The code shifts the qualities into ASCII
/// here, which is the inverse of the decode in [`super::transform`]. It does that through
/// `scratch`, so the shift allocates nothing at each read.
///
/// The code builds the whole record in `scratch`, and it gives that to **one** `write_all`. There
/// were seven of those: `@`, the name, a newline, the sequence, `\n+\n`, the qualities, and a
/// newline. Each one went into the state machine of the gzip encoder on its own. Each read then
/// paid that cost seven times, and not once.
///
/// A measurement on this exact stack, at 151 bp, gave **1,882 ns for each record, against
/// 539 ns**. That is a difference of 3.5x, on the write path of the stage that already holds the
/// peak of the scratch space. At about 600M reads for a 30x WGS, that is about thirteen minutes of
/// CPU on one thread, for each realignment. The bytes that come out are the same.
fn write_record(
    w: &mut FastqWriter,
    read: &RevertedRead,
    path: &Path,
    scratch: &mut Vec<u8>,
) -> Result<(), AnalysisError> {
    scratch.clear();
    scratch.push(b'@');
    scratch.extend_from_slice(&read.name);
    scratch.push(b'\n');
    scratch.extend_from_slice(&read.sequence);
    scratch.extend_from_slice(b"\n+\n");
    scratch.extend(read.qualities.iter().map(|q| q.saturating_add(PHRED_OFFSET)));
    scratch.push(b'\n');

    w.write_all(scratch).map_err(|e| AnalysisError::io(path, e))
}
