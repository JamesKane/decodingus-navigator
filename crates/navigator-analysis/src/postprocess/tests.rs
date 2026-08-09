//! Tests for the coordinate sort.
//!
//! The property that matters most is that the sort is *lossless* — a sort that quietly drops
//! records would show up downstream only as coverage reading a bit low, which is close to
//! undetectable. The spill/no-spill equivalence test is the second: the result must not depend on
//! how much of the input happened to fit in memory.

use std::path::{Path, PathBuf};

use noodles::sam::alignment::io::Write as _;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record_buf::{QualityScores, Sequence};
use noodles::sam::alignment::RecordBuf;
use noodles::sam::header::record::value::{map, Map};
use noodles::{bam, sam};

use super::sort::*;
use crate::cancel::CancelToken;
use crate::error::AnalysisError;

fn scratch(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("dun-sort-{}-{tag}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn header(contigs: &[(&str, usize)]) -> sam::Header {
    let mut builder = sam::Header::builder();
    for (name, len) in contigs {
        builder = builder.add_reference_sequence(
            name.as_bytes(),
            Map::<map::ReferenceSequence>::new(std::num::NonZeroUsize::new(*len).unwrap()),
        );
    }
    builder.build()
}

/// A record at `(ref_id, pos)`, or unplaced when `ref_id` is `None`.
fn record(name: &str, ref_id: Option<usize>, pos: usize) -> RecordBuf {
    let mut builder = RecordBuf::builder()
        .set_name(name)
        .set_sequence(Sequence::from(b"ACGTACGTAC".to_vec()))
        .set_quality_scores(QualityScores::from(vec![30; 10]));
    match ref_id {
        Some(id) => {
            builder = builder
                .set_flags(Flags::from(0u16))
                .set_reference_sequence_id(id)
                .set_alignment_start(noodles::core::Position::new(pos).unwrap());
        }
        None => builder = builder.set_flags(Flags::UNMAPPED),
    }
    builder.build()
}

fn write_bam(path: &Path, header: &sam::Header, records: &[RecordBuf]) {
    let file = std::fs::File::create(path).unwrap();
    let mut writer = bam::io::Writer::new(file);
    writer.write_header(header).unwrap();
    for record in records {
        writer.write_alignment_record(header, record).unwrap();
    }
    writer.try_finish().unwrap();
}

/// `(name, ref_id, pos)` for every record, in file order.
fn read_back(path: &Path) -> (sam::Header, Vec<(String, Option<usize>, usize)>) {
    let file = std::fs::File::open(path).unwrap();
    let mut reader = bam::io::Reader::new(file);
    let header = reader.read_header().unwrap();
    let records = reader
        .record_bufs(&header)
        .map(|r| {
            let r = r.unwrap();
            (
                String::from_utf8_lossy(r.name().unwrap_or_default()).to_string(),
                r.reference_sequence_id(),
                r.alignment_start().map(|p| p.get()).unwrap_or(0),
            )
        })
        .collect();
    (header, records)
}

fn run(input: &Path, dir: &Path, params: SortParams) -> (PathBuf, SortStats) {
    let output = dir.join("sorted.bam");
    let stats = sort_alignment(
        input,
        &output,
        &dir.join("scratch"),
        &params,
        &CancelToken::none(),
        &mut |_| {},
    )
    .expect("sort should succeed");
    (output, stats)
}

/// Shuffled input across two contigs, plus unplaced reads.
fn unsorted_fixture(dir: &Path) -> (PathBuf, sam::Header, usize) {
    let hdr = header(&[("chr1", 100_000), ("chr2", 100_000)]);
    let mut records = Vec::new();
    for i in 0..60 {
        // Deliberately out of order, and interleaved between the two contigs.
        let pos = ((i * 4703) % 90_000) + 1;
        let ref_id = i % 2;
        records.push(record(&format!("r{i:03}"), Some(ref_id), pos));
    }
    for i in 0..5 {
        records.push(record(&format!("u{i}"), None, 0));
    }
    let input = dir.join("unsorted.bam");
    let total = records.len();
    write_bam(&input, &hdr, &records);
    (input, hdr, total)
}

// ---- the properties that matter -------------------------------------------

/// Coordinate order, with unplaced reads at the end where SAM puts them.
#[test]
fn records_come_out_in_coordinate_order_with_unplaced_reads_last() {
    let dir = scratch("order");
    let (input, _, total) = unsorted_fixture(&dir);
    let (output, stats) = run(&input, &dir, SortParams::default());

    assert_eq!(stats.records as usize, total);
    assert_eq!(stats.unplaced, 5);

    let (_, records) = read_back(&output);
    assert_eq!(records.len(), total);

    let placed: Vec<_> = records.iter().filter(|r| r.1.is_some()).collect();
    let unplaced: Vec<_> = records.iter().filter(|r| r.1.is_none()).collect();
    assert_eq!(unplaced.len(), 5);
    assert_eq!(
        records[records.len() - 5..].iter().filter(|r| r.1.is_none()).count(),
        5,
        "the unplaced reads must be the final records, not scattered"
    );

    for pair in placed.windows(2) {
        let (a, b) = (&pair[0], &pair[1]);
        assert!((a.1, a.2) <= (b.1, b.2), "out of order: {:?} before {:?}", a, b);
    }
}

/// A sort that loses records would show up downstream only as slightly low coverage — which is
/// close to undetectable. Nothing may be dropped, including the unplaced reads that realignment
/// exists to recover.
#[test]
fn the_sort_is_lossless() {
    let dir = scratch("lossless");
    let (input, _, total) = unsorted_fixture(&dir);
    let (output, stats) = run(&input, &dir, SortParams { buffer_bytes: 1 });

    let (_, records) = read_back(&output);
    assert_eq!(records.len(), total);
    assert_eq!(stats.records as usize, total);

    let mut names: Vec<&str> = records.iter().map(|r| r.0.as_str()).collect();
    names.sort_unstable();
    let mut expected: Vec<String> = (0..60).map(|i| format!("r{i:03}")).collect();
    expected.extend((0..5).map(|i| format!("u{i}")));
    expected.sort_unstable();
    assert_eq!(names, expected, "every input record must appear exactly once");
}

/// The result must not depend on how much of the input fit in memory — the same claim the revert
/// stage's external sort makes, and the reason the disk path is safe to rely on.
#[test]
fn spilling_produces_the_same_order_as_sorting_in_memory() {
    let dir = scratch("spill");
    let (input, _, _) = unsorted_fixture(&dir);

    let in_memory_dir = dir.join("mem");
    std::fs::create_dir_all(&in_memory_dir).unwrap();
    let (a, a_stats) = run(
        &input,
        &in_memory_dir,
        SortParams {
            buffer_bytes: 64 * 1024 * 1024,
        },
    );

    let spilled_dir = dir.join("spill");
    std::fs::create_dir_all(&spilled_dir).unwrap();
    let (b, b_stats) = run(&input, &spilled_dir, SortParams { buffer_bytes: 1 });

    assert_eq!(a_stats.runs, 1, "everything fit");
    assert!(b_stats.runs > 1, "the budget forced multiple runs");
    assert_eq!(read_back(&a).1, read_back(&b).1);
}

/// An index is only valid for a coordinate-sorted file, and readers decide whether they may query
/// a region by reading this. Sorting correctly but failing to say so means every reader rescans.
#[test]
fn the_output_header_declares_coordinate_order() {
    let dir = scratch("header");
    let (input, _, _) = unsorted_fixture(&dir);
    let (output, _) = run(&input, &dir, SortParams::default());

    let (header, _) = read_back(&output);
    let hd = header.header().expect("@HD must be present");
    use noodles::sam::header::record::value::map::header::tag;
    let so = hd.other_fields().get(&tag::SORT_ORDER).expect("@HD must carry SO");
    assert_eq!(so.as_slice(), b"coordinate");
}

/// Runs are the size of the alignment itself, so they must not outlive the sort.
#[test]
fn spilled_runs_are_cleaned_up() {
    let dir = scratch("cleanup");
    let (input, _, _) = unsorted_fixture(&dir);
    let _ = run(&input, &dir, SortParams { buffer_bytes: 1 });

    let leftovers: Vec<String> = std::fs::read_dir(dir.join("scratch"))
        .unwrap()
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().to_string())
        .filter(|n| n.starts_with("sort-run-"))
        .collect();
    assert!(leftovers.is_empty(), "run files left behind: {leftovers:?}");
}

/// Sorting a WGS is long enough to need a cancel that works, reported as itself.
#[test]
fn cancellation_stops_the_sort() {
    let dir = scratch("cancel");
    let (input, _, _) = unsorted_fixture(&dir);
    let cancel = CancelToken::new();
    cancel.cancel();

    let err = sort_alignment(
        &input,
        &dir.join("out.bam"),
        &dir.join("scratch"),
        &SortParams::default(),
        &cancel,
        &mut |_| {},
    )
    .unwrap_err();
    assert!(matches!(err, AnalysisError::Cancelled));
}

/// An already-sorted file is a normal input (a re-run, or a mapper that happened to emit in
/// order), and must come out unchanged rather than subtly reordered.
#[test]
fn an_already_sorted_file_is_unchanged() {
    let dir = scratch("idempotent");
    let (input, _, _) = unsorted_fixture(&dir);
    let (once, _) = run(&input, &dir, SortParams::default());

    let twice_dir = dir.join("again");
    std::fs::create_dir_all(&twice_dir).unwrap();
    let (twice, _) = run(&once, &twice_dir, SortParams::default());

    assert_eq!(read_back(&once).1, read_back(&twice).1);
}
