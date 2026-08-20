//! Tests for the revert stage.
//!
//! These tests build their records by hand, and they do not read a BAM fixture. The container
//! layer belongs to [`crate::reader`], and the tests there cover it. Everything that matters here
//! is independent of how the bytes arrived. That is the flags, the orientation, how two mates
//! pair, and the spill and merge boundary of the external sort.
//!
//! A record that a test builds directly also lets that test make the inputs that are not correct,
//! and those are the ones that matter. A real aligner rarely gives flags that contradict
//! themselves, or a record with no qualities.

use noodles::sam::alignment::record::cigar::op::{Kind, Op};
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record_buf::data::field::Value;
use noodles::sam::alignment::record_buf::{Cigar, Data, QualityScores, Sequence};
use noodles::sam::alignment::RecordBuf;

use super::*;

const OQ: Tag = Tag::new(b'O', b'Q');

/// A scratch directory of its own for each test, under the temp directory of the system. That is
/// the convention elsewhere in this crate, and it needs no `tempfile` dependency.
fn scratch(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("dun-revert-{}-{tag}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// A small record. It holds a name, the flags, a sequence, and qualities whose length matches
/// that sequence.
fn record(name: &str, flags: u16, seq: &str, quals: &[u8]) -> RecordBuf {
    RecordBuf::builder()
        .set_name(name)
        .set_flags(Flags::from(flags))
        .set_sequence(Sequence::from(seq.as_bytes().to_vec()))
        .set_quality_scores(QualityScores::from(quals.to_vec()))
        .build()
}

fn paired(name: &str, first: bool, reverse: bool, seq: &str, quals: &[u8]) -> RecordBuf {
    let mut flags = 0x1u16 | if first { 0x40 } else { 0x80 };
    if reverse {
        flags |= 0x10;
    }
    record(name, flags, seq, quals)
}

fn run(records: Vec<RecordBuf>, dir: &Path, params: &RevertParams) -> RevertOutput {
    revert_records(records.into_iter().map(Ok), dir, params, &CancelToken::none()).unwrap()
}

fn read_lines(path: &Path) -> Vec<String> {
    use std::io::BufRead;
    // The reverted FASTQ goes through gzip. writer.rs says why. So a test reads it the same way
    // that the mapper does, and it does not expect plain text.
    crate::gzio::open_maybe_gz(path)
        .unwrap()
        .lines()
        .map(|l| l.unwrap())
        .collect()
}

// ---- the transform --------------------------------------------------------

/// The core property. A reverse-strand alignment stores the read the other way round. The FASTQ
/// must carry what the sequencer gave, and not what the aligner stored.
#[test]
fn a_reverse_strand_read_is_restored_to_sequencer_orientation() {
    let dir = scratch("revcomp");
    let out = run(
        vec![record("r", 0x10, "AACCG", &[10, 20, 30, 40, 50])],
        &dir,
        &RevertParams::default(),
    );

    let lines = read_lines(&out.singletons);
    assert_eq!(lines[1], "CGGTT", "reverse-complemented back");
    // Qualities travel with the bases they describe, so they reverse too (but are not complemented).
    assert_eq!(
        lines[3].bytes().map(|b| b - 33).collect::<Vec<_>>(),
        vec![50, 40, 30, 20, 10]
    );
}

/// A forward-strand read must go through unchanged. This is the mirror of the test above. A bug
/// that takes the reverse complement of every read can then not pass both.
#[test]
fn a_forward_strand_read_is_left_alone() {
    let dir = scratch("forward");
    let out = run(
        vec![record("r", 0x0, "AACCG", &[10, 20, 30, 40, 50])],
        &dir,
        &RevertParams::default(),
    );
    let lines = read_lines(&out.singletons);
    assert_eq!(lines[1], "AACCG");
    assert_eq!(
        lines[3].bytes().map(|b| b - 33).collect::<Vec<_>>(),
        vec![10, 20, 30, 40, 50]
    );
}

/// A secondary record and a supplementary one each repeat a read whose full sequence lives on the
/// primary. To keep them would emit the same read more than once, and a supplementary one would
/// come out short.
#[test]
fn secondary_and_supplementary_records_are_dropped() {
    let dir = scratch("nonprimary");
    let out = run(
        vec![
            record("r", 0x0, "ACGT", &[30; 4]),
            record("r", 0x100, "ACG", &[30; 3]),
            record("r", 0x800, "CGT", &[30; 3]),
        ],
        &dir,
        &RevertParams::default(),
    );

    assert_eq!(out.stats.secondary_dropped, 1);
    assert_eq!(out.stats.supplementary_dropped, 1);
    assert_eq!(out.stats.reads_emitted, 1, "only the primary survives");
}

/// The reads with no mapping are the gain of the realignment. They must come through, and the code
/// must count them, so that somebody can measure the gain.
#[test]
fn unmapped_reads_are_kept_and_counted() {
    let dir = scratch("unmapped");
    let out = run(vec![record("u", 0x4, "ACGT", &[30; 4])], &dir, &RevertParams::default());
    assert_eq!(out.stats.unmapped_reads, 1);
    assert_eq!(out.stats.reads_emitted, 1);
}

/// `OQ` holds the qualities from before a recalibration wrote over `QUAL`, in ASCII. To take `OQ`
/// first is the difference between a revert to the original read and a revert to one that a
/// pipeline already changed.
#[test]
fn original_qualities_are_preferred_over_recalibrated_ones() {
    let dir = scratch("oq");
    let rec = RecordBuf::builder()
        .set_name("r")
        .set_flags(Flags::from(0x0u16))
        .set_sequence(Sequence::from(b"ACGT".to_vec()))
        .set_quality_scores(QualityScores::from(vec![2, 2, 2, 2]))
        // Phred 40 as FASTQ ASCII ('I' == 40 + 33).
        .set_data([(OQ, Value::String("IIII".into()))].into_iter().collect::<Data>())
        .build();

    let out = run(vec![rec], &dir, &RevertParams::default());
    assert_eq!(out.stats.original_qualities_used, 1);
    assert_eq!(read_lines(&out.singletons)[3], "IIII", "OQ won, decoded and re-encoded");
}

/// The option that turns this off must turn it off. Else the flag says nothing.
#[test]
fn original_qualities_can_be_declined() {
    let dir = scratch("oq-off");
    let rec = RecordBuf::builder()
        .set_name("r")
        .set_flags(Flags::from(0x0u16))
        .set_sequence(Sequence::from(b"ACGT".to_vec()))
        .set_quality_scores(QualityScores::from(vec![2, 2, 2, 2]))
        .set_data([(OQ, Value::String("IIII".into()))].into_iter().collect::<Data>())
        .build();

    let params = RevertParams {
        prefer_original_qualities: false,
        ..Default::default()
    };
    let out = run(vec![rec], &dir, &params);
    assert_eq!(out.stats.original_qualities_used, 0);
    assert_eq!(read_lines(&out.singletons)[3], "####", "QUAL 2 == '#', one per base");
}

/// A primary record with a hard clip has already lost sequence. To skip it is the default. A read
/// that the code drops shows in the statistics, and a read that comes out short does not.
#[test]
fn hard_clipped_primaries_are_skipped_by_default_and_emittable_on_request() {
    let rec = || {
        RecordBuf::builder()
            .set_name("h")
            .set_flags(Flags::from(0x0u16))
            .set_sequence(Sequence::from(b"ACGT".to_vec()))
            .set_quality_scores(QualityScores::from(vec![30; 4]))
            .set_cigar(Cigar::from(vec![Op::new(Kind::HardClip, 10), Op::new(Kind::Match, 4)]))
            .build()
    };

    let skipped = run(vec![rec()], &scratch("hardclip-skip"), &RevertParams::default());
    assert_eq!(skipped.stats.hard_clipped, 1);
    assert_eq!(skipped.stats.reads_emitted, 0);

    let params = RevertParams {
        hard_clipped: HardClipPolicy::Emit,
        ..Default::default()
    };
    let emitted = run(vec![rec()], &scratch("hardclip-emit"), &params);
    assert_eq!(emitted.stats.hard_clipped, 1, "counted under either policy");
    assert_eq!(emitted.stats.reads_emitted, 1);
}

/// A `QUAL` of `*` is legal. A mapper can still map the read, so the code keeps it. But the code
/// invents the qualities that come out, and the statistic is the one thing that says so.
#[test]
fn missing_qualities_are_synthesized_and_counted() {
    let dir = scratch("noqual");
    let out = run(vec![record("r", 0x0, "ACGT", &[])], &dir, &RevertParams::default());
    assert_eq!(out.stats.qualities_synthesized, 1);
    let lines = read_lines(&out.singletons);
    assert_eq!(lines[3].len(), 4, "one quality per base, as FASTQ requires");
}

// ---- how two mates pair ---------------------------------------------------

/// The main case. In coordinate order the two mates arrive far apart, and they must come back
/// together.
#[test]
fn mates_separated_in_the_input_are_paired_in_the_output() {
    let dir = scratch("pairing");
    // Interleave three templates so no pair is adjacent on input.
    let out = run(
        vec![
            paired("a", true, false, "AAAA", &[30; 4]),
            paired("b", true, false, "BBBB", &[30; 4]),
            paired("c", true, false, "CCCC", &[30; 4]),
            paired("b", false, false, "TTTT", &[30; 4]),
            paired("c", false, false, "GGGG", &[30; 4]),
            paired("a", false, false, "CCCC", &[30; 4]),
        ],
        &dir,
        &RevertParams::default(),
    );

    assert_eq!(out.stats.pairs, 3);
    assert_eq!(out.stats.singletons, 0);

    let r1 = read_lines(&out.read1);
    let r2 = read_lines(&out.read2);
    // The names come in order. And, as the invariant that the mapper depends on says, R1 and R2
    // stay in step.
    assert_eq!([r1[0].as_str(), r1[4].as_str(), r1[8].as_str()], ["@a", "@b", "@c"]);
    assert_eq!(r1.len(), r2.len(), "files stay the same length");
    for i in (0..r1.len()).step_by(4) {
        assert_eq!(r1[i], r2[i], "record {} pairs the same template", i / 4);
    }
}

/// A mate whose partner the code dropped must not go into `_1`. That would move every later pair
/// by one, and the rest of the file would hold the wrong pairs.
#[test]
fn a_read_whose_mate_was_dropped_becomes_a_singleton() {
    let dir = scratch("orphan");
    let out = run(
        vec![
            paired("a", true, false, "AAAA", &[30; 4]),
            paired("a", false, false, "TTTT", &[30; 4]),
            // 'b' has only its first segment; its mate never appears.
            paired("b", true, false, "GGGG", &[30; 4]),
        ],
        &dir,
        &RevertParams::default(),
    );

    assert_eq!(out.stats.pairs, 1);
    assert_eq!(out.stats.singletons, 1);
    assert_eq!(read_lines(&out.read1).len(), 4, "only the complete pair");
    assert_eq!(read_lines(&out.singletons)[0], "@b");
}

/// Flags that say "part of a pair", and that do not say which end, have no place in a file that
/// must stay in step.
#[test]
fn a_paired_record_with_contradictory_segment_flags_is_a_singleton() {
    let dir = scratch("contradictory");
    // 0x1 segmented, with both 0x40 and 0x80 set.
    let out = run(
        vec![record("x", 0x1 | 0x40 | 0x80, "ACGT", &[30; 4])],
        &dir,
        &RevertParams::default(),
    );
    assert_eq!(out.stats.pairs, 0);
    assert_eq!(out.stats.singletons, 1);
}

/// Two records with the same name *and* the same end have no answer. To take one of them would
/// drop a read where nobody sees it, and it could put the wrong reads together.
#[test]
fn duplicate_segment_bits_under_one_name_do_not_pair() {
    let dir = scratch("dupe-segment");
    let out = run(
        vec![
            paired("d", true, false, "AAAA", &[30; 4]),
            paired("d", true, false, "CCCC", &[30; 4]),
            paired("d", false, false, "GGGG", &[30; 4]),
        ],
        &dir,
        &RevertParams::default(),
    );
    assert_eq!(out.stats.pairs, 0);
    assert_eq!(out.stats.singletons, 3, "all three go to singletons");
}

// ---- the external sort ----------------------------------------------------

/// The property that the whole design stands on: the result must not depend on whether the input
/// fit in memory. A budget of 1 byte makes the code spill at every read, and that covers the k-way
/// merge. The output must match the output of the one-run case exactly.
#[test]
fn spilling_to_disk_produces_the_same_output_as_sorting_in_memory() {
    let records = || {
        let mut v = Vec::new();
        for i in 0..50 {
            // Names deliberately out of lexical order relative to emission order.
            let name = format!("read{:03}", (i * 37) % 50);
            v.push(paired(&name, true, false, "ACGT", &[30; 4]));
            v.push(paired(&name, false, true, "ACGT", &[30; 4]));
        }
        v
    };

    let in_memory = run(
        records(),
        &scratch("sort-memory"),
        &RevertParams {
            sort_buffer_bytes: 64 * 1024 * 1024,
            ..Default::default()
        },
    );
    let spilled = run(
        records(),
        &scratch("sort-spill"),
        &RevertParams {
            sort_buffer_bytes: 1,
            ..Default::default()
        },
    );

    assert_eq!(in_memory.stats.runs_spilled, 1, "everything fit");
    assert!(spilled.stats.runs_spilled > 1, "actually spilled multiple runs");
    assert_eq!(in_memory.stats.pairs, 50);
    assert_eq!(spilled.stats.pairs, in_memory.stats.pairs);
    assert_eq!(
        read_lines(&spilled.read1),
        read_lines(&in_memory.read1),
        "identical output regardless of how much fit in memory"
    );
    assert_eq!(read_lines(&spilled.read2), read_lines(&in_memory.read2));
}

/// The names must come out in sorted order after a merge across many runs. The code that makes the
/// groups reads a run of equal names. A merge that does not sort would then split a template in
/// two, and nobody would see it.
#[test]
fn the_merge_emits_names_in_sorted_order() {
    let dir = scratch("sorted");
    let mut records = Vec::new();
    for i in (0..30).rev() {
        records.push(record(&format!("r{i:03}"), 0x0, "ACGT", &[30; 4]));
    }
    let out = run(
        records,
        &dir,
        &RevertParams {
            sort_buffer_bytes: 1,
            ..Default::default()
        },
    );

    let names: Vec<String> = read_lines(&out.singletons).into_iter().step_by(4).collect();
    let mut sorted = names.clone();
    sorted.sort();
    assert_eq!(names, sorted);
    assert_eq!(names.len(), 30);
}

/// The scratch files are large, at tens of GB for a WGS, so none of them must stay after the run
/// ends.
#[test]
fn run_files_are_cleaned_up() {
    let dir = scratch("cleanup");
    let mut records = Vec::new();
    for i in 0..20 {
        records.push(record(&format!("r{i:03}"), 0x0, "ACGT", &[30; 4]));
    }
    let _ = run(
        records,
        &dir,
        &RevertParams {
            sort_buffer_bytes: 1,
            ..Default::default()
        },
    );

    let leftovers: Vec<_> = std::fs::read_dir(&dir)
        .unwrap()
        .filter_map(|e| e.ok())
        .map(|e| e.file_name().to_string_lossy().to_string())
        .filter(|n| n.starts_with("revert-run-"))
        .collect();
    assert!(leftovers.is_empty(), "spill files left behind: {leftovers:?}");
}

// ---- cancellation ---------------------------------------------------------

/// A revert takes hours. A token that somebody already cancelled must stop it, and the job must
/// not run to its end. The job must also report itself as cancelled, and not as a failure.
#[test]
fn an_already_cancelled_token_stops_the_revert() {
    let dir = scratch("cancel");
    let cancel = CancelToken::new();
    cancel.cancel();

    let records = (0..10).map(|i| Ok(record(&format!("r{i}"), 0x0, "ACGT", &[30; 4])));
    let err = revert_records(records, &dir, &RevertParams::default(), &cancel).unwrap_err();
    assert!(matches!(err, AnalysisError::Cancelled));
}
