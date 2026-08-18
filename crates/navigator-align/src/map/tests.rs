//! Tests for the mapping pass.
//!
//! The one that matters is [`a_split_index_places_reads_exactly_where_a_whole_index_does`]. Every
//! other property here is ordinary; that one is the design's central claim, and if it fails the
//! memory bound that justifies this whole module is not free after all.

use std::path::PathBuf;

use super::*;
use crate::batch::BatchSize;
use crate::index::build_index;

fn scratch(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("dun-map-{}-{tag}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// A deterministic pseudo-random reference. Real-ish base composition matters: a low-complexity
/// sequence collapses into a few minimizer buckets and maps ambiguously everywhere, which would
/// make these tests measure the fixture rather than the code.
fn reference_bases(contig: usize, len: usize) -> Vec<u8> {
    let mut state = 0x9E3779B97F4A7C15u64 ^ (contig as u64).wrapping_mul(0xD1B54A32D192ED03);
    (0..len)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[(state >> 33) as usize % 4]
        })
        .collect()
}

fn write_reference(dir: &Path, contigs: usize, len: usize) -> PathBuf {
    let path = dir.join("ref.fa");
    let mut text = Vec::new();
    for c in 0..contigs {
        text.extend_from_slice(format!(">chr{c}\n").as_bytes());
        for chunk in reference_bases(c, len).chunks(60) {
            text.extend_from_slice(chunk);
            text.push(b'\n');
        }
    }
    std::fs::write(&path, text).unwrap();
    path
}

/// Reads lifted straight out of the reference, so the true origin of each is known from its name.
fn write_reads(dir: &Path, contigs: usize, len: usize, per_contig: usize, read_len: usize) -> PathBuf {
    let path = dir.join("reads.fq");
    let mut text = Vec::new();
    for c in 0..contigs {
        let bases = reference_bases(c, len);
        for i in 0..per_contig {
            let start = (i * 7919) % (len - read_len);
            let read = &bases[start..start + read_len];
            text.extend_from_slice(format!("@chr{c}_{start}\n").as_bytes());
            text.extend_from_slice(read);
            text.extend_from_slice(b"\n+\n");
            text.extend_from_slice(&vec![b'I'; read_len]);
            text.push(b'\n');
        }
    }
    std::fs::write(&path, text).unwrap();
    path
}

fn no_cancel() -> impl Fn() -> bool {
    || false
}

fn run(index: &Path, reads: &Path, out: &Path, dir: &Path, preset: Preset) -> MapStats {
    map_reads(
        index,
        reads,
        out,
        &dir.join("scratch"),
        &MapParams {
            preset,
            threads: 1,
            read_group: None,
            // These assertions read columns out of SAM text, so keep the text container here;
            // BAM output has its own round-trip test.
            format: OutputFormat::Sam,
            reference: None,
        },
        &no_cancel(),
        &mut |_, _, _| {},
    )
    .expect("mapping should succeed")
}

/// SAM alignment lines as `(qname, rname, pos, flag)`, headers dropped.
fn placements(sam: &Path) -> Vec<(String, String, String, String)> {
    std::fs::read_to_string(sam)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[0].to_string(), f[2].to_string(), f[3].to_string(), f[1].to_string())
        })
        .collect()
}

fn mapq(sam: &Path) -> Vec<(String, String)> {
    std::fs::read_to_string(sam)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            (f[0].to_string(), f[4].to_string())
        })
        .collect()
}

// ---- the claim ------------------------------------------------------------

/// **The central property.** Splitting the index is a memory optimization; it must not be an
/// accuracy decision. A read mapped against a 3-part index has to land exactly where it lands
/// against a whole one — same contig, same position, same flags — because the per-part hits are
/// merged and re-ranked before anything is emitted.
///
/// If this ever fails, the memory table in `crate::batch` stops being free and the design has to
/// be re-argued.
#[test]
fn a_split_index_places_reads_exactly_where_a_whole_index_does() {
    let dir = scratch("equivalence");
    let (contigs, len) = (6, 500_000);
    let reference = write_reference(&dir, contigs, len);
    let reads = write_reads(&dir, contigs, len, 40, 150);

    // Whole index: one part, no merge.
    let whole = dir.join("whole.mmi");
    build_index(
        &reference,
        &whole,
        Preset::ShortRead,
        BatchSize::new(8_000_000_000),
        &mut |_, _| {},
    )
    .unwrap();
    let whole_sam = dir.join("whole.sam");
    let whole_stats = run(&whole, &reads, &whole_sam, &dir, Preset::ShortRead);

    // Split index: several parts, exercising the write-merge path.
    let split_idx = dir.join("split.mmi");
    build_index(
        &reference,
        &split_idx,
        Preset::ShortRead,
        BatchSize::new(1_000_000),
        &mut |_, _| {},
    )
    .unwrap();
    let split_sam = dir.join("split.sam");
    let split_stats = run(&split_idx, &reads, &split_sam, &dir, Preset::ShortRead);

    assert_eq!(whole_stats.parts, 1, "the whole index must be one part");
    assert!(split_stats.parts > 1, "the split index must actually split");
    assert_eq!(split_stats.queries, whole_stats.queries);
    assert_eq!(
        split_stats.mapped, whole_stats.mapped,
        "the same reads must map either way"
    );
    assert_eq!(
        placements(&split_sam),
        placements(&whole_sam),
        "a split index must place every read identically to a whole one"
    );
}

/// MAPQ is the part of the answer a split index is *allowed* to differ on — a read's second-best
/// hit can fall in another part — but the merge re-runs the MAPQ calculation across all parts
/// precisely so it does not. On a reference with no duplicated sequence there is nothing to be
/// ambiguous about, so it must agree exactly.
#[test]
fn merging_restores_mapq_across_parts() {
    let dir = scratch("mapq");
    let (contigs, len) = (6, 500_000);
    let reference = write_reference(&dir, contigs, len);
    let reads = write_reads(&dir, contigs, len, 30, 150);

    let whole = dir.join("whole.mmi");
    build_index(
        &reference,
        &whole,
        Preset::ShortRead,
        BatchSize::new(8_000_000_000),
        &mut |_, _| {},
    )
    .unwrap();
    let whole_sam = dir.join("whole.sam");
    run(&whole, &reads, &whole_sam, &dir, Preset::ShortRead);

    let split_idx = dir.join("split.mmi");
    build_index(
        &reference,
        &split_idx,
        Preset::ShortRead,
        BatchSize::new(1_000_000),
        &mut |_, _| {},
    )
    .unwrap();
    let split_sam = dir.join("split.sam");
    run(&split_idx, &reads, &split_sam, &dir, Preset::ShortRead);

    assert_eq!(mapq(&split_sam), mapq(&whole_sam));
}

// ---- ordinary properties --------------------------------------------------

/// Reads were lifted out of the reference, so they must come back to the contig they came from.
/// Without this the equivalence test above could pass with both sides equally wrong.
#[test]
fn reads_map_back_to_where_they_came_from() {
    let dir = scratch("placement");
    let (contigs, len) = (3, 300_000);
    let reference = write_reference(&dir, contigs, len);
    let reads = write_reads(&dir, contigs, len, 20, 150);

    let index = dir.join("ref.mmi");
    build_index(
        &reference,
        &index,
        Preset::ShortRead,
        BatchSize::default(),
        &mut |_, _| {},
    )
    .unwrap();
    let sam = dir.join("out.sam");
    let stats = run(&index, &reads, &sam, &dir, Preset::ShortRead);

    assert_eq!(stats.queries, 60);
    assert!(stats.mapped >= 58, "expected almost all to map, got {}", stats.mapped);

    for (qname, rname, pos, _) in placements(&sam) {
        if rname == "*" {
            continue;
        }
        // Name is `chr<c>_<start>`; SAM POS is 1-based.
        let (expect_chr, expect_start) = qname.split_once('_').unwrap();
        assert_eq!(rname, expect_chr, "{qname} landed on the wrong contig");
        let pos: i64 = pos.parse().unwrap();
        let expected: i64 = expect_start.parse::<i64>().unwrap() + 1;
        assert!(
            (pos - expected).abs() <= 5,
            "{qname} placed at {pos}, expected ~{expected}"
        );
    }
}

/// A SAM record without a CIGAR is not an alignment, it is a coordinate guess — and nothing
/// downstream (coverage, callable, SV, the variant caller) can consume it.
///
/// This exists because its absence was invisible: the placement and equivalence tests all passed
/// while every record carried `*`, since chaining alone yields coordinates. It also pins the
/// primary/supplementary flag, which failed the same way and for the same underlying reason —
/// without base-level alignment, `map_query` never ranks the regions it returns.
#[test]
fn records_carry_a_real_cigar_and_a_primary_flag() {
    let dir = scratch("cigar");
    let (contigs, len) = (2, 300_000);
    let reference = write_reference(&dir, contigs, len);
    let reads = write_reads(&dir, contigs, len, 10, 150);

    let index = dir.join("ref.mmi");
    build_index(
        &reference,
        &index,
        Preset::ShortRead,
        BatchSize::default(),
        &mut |_, _| {},
    )
    .unwrap();
    let sam = dir.join("out.sam");
    run(&index, &reads, &sam, &dir, Preset::ShortRead);

    let mut primaries = 0;
    for line in std::fs::read_to_string(&sam).unwrap().lines() {
        if line.starts_with('@') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let (flag, cigar) = (f[1].parse::<u16>().unwrap(), f[5]);
        if flag & 0x4 != 0 {
            continue; // unmapped records legitimately carry `*`
        }
        assert_ne!(cigar, "*", "mapped record has no CIGAR: {line}");
        assert!(
            cigar.contains('M') || cigar.contains('='),
            "CIGAR has no aligned bases: {cigar}"
        );
        if flag & 0x900 == 0 {
            primaries += 1;
        }
    }
    assert_eq!(primaries, 20, "every read should have exactly one primary record");
}

/// The SAM has to be readable by everything downstream, which starts with a header naming every
/// contig — including, in the split case, contigs from parts that were never resident together.
#[test]
fn the_header_names_every_contig_across_every_part() {
    let dir = scratch("header");
    let (contigs, len) = (6, 500_000);
    let reference = write_reference(&dir, contigs, len);
    let reads = write_reads(&dir, contigs, len, 2, 150);

    let index = dir.join("split.mmi");
    build_index(
        &reference,
        &index,
        Preset::ShortRead,
        BatchSize::new(1_000_000),
        &mut |_, _| {},
    )
    .unwrap();
    let sam = dir.join("out.sam");
    let stats = run(&index, &reads, &sam, &dir, Preset::ShortRead);
    assert!(stats.parts > 1);

    let text = std::fs::read_to_string(&sam).unwrap();
    for c in 0..contigs {
        assert!(
            text.contains(&format!("SN:chr{c}\t")),
            "@SQ missing chr{c}:\n{}",
            text.lines().take(10).collect::<Vec<_>>().join("\n")
        );
    }
    assert!(text.contains("@PG"), "a realigned header must carry a @PG record");
}

/// An unmappable read gets a record rather than vanishing. Which reads failed is information —
/// realignment exists partly to recover reads a previous reference could not place.
#[test]
fn unmappable_reads_are_written_as_unmapped_records() {
    let dir = scratch("unmapped");
    let reference = write_reference(&dir, 2, 200_000);
    let index = dir.join("ref.mmi");
    build_index(
        &reference,
        &index,
        Preset::ShortRead,
        BatchSize::default(),
        &mut |_, _| {},
    )
    .unwrap();

    // Poly-A shares nothing with the pseudo-random reference.
    let reads = dir.join("junk.fq");
    std::fs::write(
        &reads,
        "@junk\n".to_string() + &"A".repeat(150) + "\n+\n" + &"I".repeat(150) + "\n",
    )
    .unwrap();

    let sam = dir.join("out.sam");
    let stats = run(&index, &reads, &sam, &dir, Preset::ShortRead);

    assert_eq!(stats.queries, 1);
    assert_eq!(stats.unmapped, 1);
    let records = placements(&sam);
    assert_eq!(records.len(), 1, "the read is present, not dropped");
    assert_eq!(records[0].0, "junk");
    assert_eq!(records[0].3, "4", "SAM flag 4 == unmapped");
}

/// Per-part hit blocks are one file per part for the whole read set — at genome scale that is
/// large, so it must not survive the call.
#[test]
fn split_scratch_is_cleaned_up() {
    let dir = scratch("cleanup");
    let (contigs, len) = (6, 500_000);
    let reference = write_reference(&dir, contigs, len);
    let reads = write_reads(&dir, contigs, len, 5, 150);

    let index = dir.join("split.mmi");
    build_index(
        &reference,
        &index,
        Preset::ShortRead,
        BatchSize::new(1_000_000),
        &mut |_, _| {},
    )
    .unwrap();
    let stats = run(&index, &reads, &dir.join("out.sam"), &dir, Preset::ShortRead);
    assert!(stats.parts > 1);

    let leftovers: Vec<String> = std::fs::read_dir(dir.join("scratch"))
        .map(|entries| {
            entries
                .filter_map(|e| e.ok())
                .map(|e| e.file_name().to_string_lossy().to_string())
                .collect()
        })
        .unwrap_or_default();
    assert!(leftovers.is_empty(), "scratch left behind: {leftovers:?}");
}

/// A multi-hour job has to stop when asked, and report the stop as itself rather than as a
/// failure.
#[test]
fn cancellation_stops_the_mapping_pass() {
    let dir = scratch("cancel");
    let reference = write_reference(&dir, 2, 200_000);
    let reads = write_reads(&dir, 2, 200_000, 10, 150);
    let index = dir.join("ref.mmi");
    build_index(
        &reference,
        &index,
        Preset::ShortRead,
        BatchSize::default(),
        &mut |_, _| {},
    )
    .unwrap();

    let err = map_reads(
        &index,
        &reads,
        &dir.join("out.sam"),
        &dir.join("scratch"),
        &MapParams {
            preset: Preset::ShortRead,
            threads: 1,
            read_group: None,
            // These assertions read columns out of SAM text, so keep the text container here;
            // BAM output has its own round-trip test.
            format: OutputFormat::Sam,
            reference: None,
        },
        &|| true,
        &mut |_, _, _| {},
    )
    .unwrap_err();

    assert!(matches!(err, AlignError::Cancelled), "got {err:?}");
}

/// An index file with nothing in it must fail loudly rather than produce an empty, valid-looking
/// SAM that a later stage would treat as "this sample simply has no reads".
#[test]
fn an_empty_index_is_an_error() {
    let dir = scratch("emptyidx");
    let index = dir.join("empty.mmi");
    std::fs::write(&index, b"").unwrap();
    let reads = write_reads(&dir, 1, 100_000, 2, 150);

    let err = map_reads(
        &index,
        &reads,
        &dir.join("out.sam"),
        &dir.join("scratch"),
        &MapParams::default(),
        &no_cancel(),
        &mut |_, _, _| {},
    );
    assert!(err.is_err());
}

// ---- output containers ----------------------------------------------------

/// BAM is the default container, and it has to hold exactly what SAM did. Read back through
/// noodles as *typed* records, so this asserts on fields rather than on column positions — the
/// same reason the writer exists.
#[test]
fn bam_output_round_trips_the_same_records_as_sam() {
    let dir = scratch("bam");
    let (contigs, len) = (2, 300_000);
    let reference = write_reference(&dir, contigs, len);
    let reads = write_reads(&dir, contigs, len, 20, 150);
    let index = dir.join("ref.mmi");
    build_index(
        &reference,
        &index,
        Preset::ShortRead,
        BatchSize::default(),
        &mut |_, _| {},
    )
    .unwrap();

    let params = |format| MapParams {
        preset: Preset::ShortRead,
        threads: 1,
        read_group: None,
        format,
        reference: None,
    };

    let sam = dir.join("out.sam");
    map_reads(
        &index,
        &reads,
        &sam,
        &dir.join("s1"),
        &params(OutputFormat::Sam),
        &no_cancel(),
        &mut |_, _, _| {},
    )
    .unwrap();

    let bam = dir.join("out.bam");
    map_reads(
        &index,
        &reads,
        &bam,
        &dir.join("s2"),
        &params(OutputFormat::Bam),
        &no_cancel(),
        &mut |_, _, _| {},
    )
    .unwrap();

    let (sam_header, sam_records) = crate::output::read_all(&sam).unwrap();
    let (bam_header, bam_records) = crate::output::read_all_bam(&bam).unwrap();

    assert_eq!(
        sam_header.reference_sequences().len(),
        bam_header.reference_sequences().len(),
        "both headers describe the same references"
    );
    assert_eq!(sam_records.len(), bam_records.len());
    assert!(!bam_records.is_empty());
    for (s, b) in sam_records.iter().zip(&bam_records) {
        assert_eq!(s.name(), b.name());
        assert_eq!(s.flags(), b.flags());
        assert_eq!(s.reference_sequence_id(), b.reference_sequence_id());
        assert_eq!(s.alignment_start(), b.alignment_start());
        assert_eq!(s.mapping_quality(), b.mapping_quality());
        assert_eq!(s.cigar(), b.cigar(), "CIGAR must survive the BAM encode");
        assert_eq!(s.sequence(), b.sequence());
    }
}

/// A BAM whose BGZF end-of-file block is missing reads as truncated. Writing it is the writer's
/// job on `finish`, and nothing else in the test suite would notice if it stopped happening.
#[test]
fn bam_output_is_a_complete_bgzf_stream() {
    let dir = scratch("bgzf");
    let (contigs, len) = (1, 200_000);
    let reference = write_reference(&dir, contigs, len);
    let reads = write_reads(&dir, contigs, len, 5, 150);
    let index = dir.join("ref.mmi");
    build_index(
        &reference,
        &index,
        Preset::ShortRead,
        BatchSize::default(),
        &mut |_, _| {},
    )
    .unwrap();

    let bam = dir.join("out.bam");
    map_reads(
        &index,
        &reads,
        &bam,
        &dir.join("scratch"),
        &MapParams {
            preset: Preset::ShortRead,
            threads: 1,
            read_group: None,
            format: OutputFormat::Bam,
            reference: None,
        },
        &no_cancel(),
        &mut |_, _, _| {},
    )
    .unwrap();

    // The 28-byte BGZF EOF marker, which every reader uses to tell "done" from "truncated".
    const EOF: &[u8] = &[
        0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00, 0x1b, 0x00,
        0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
    ];
    let bytes = std::fs::read(&bam).unwrap();
    assert!(bytes.ends_with(EOF), "BAM is missing its BGZF EOF block");
}

/// The container follows the filename when a caller does not say otherwise.
#[test]
fn the_output_format_can_be_read_off_the_path() {
    use crate::output::OutputFormat as F;
    assert_eq!(F::from_path(Path::new("x.sam")), F::Sam);
    assert_eq!(F::from_path(Path::new("x.bam")), F::Bam);
    assert_eq!(F::from_path(Path::new("x.cram")), F::Cram);
    assert_eq!(F::from_path(Path::new("x.CRAM")), F::Cram, "case-insensitive");
    assert_eq!(F::from_path(Path::new("x")), F::Bam, "BAM is the default");
}

/// CRAM can not be written without the reference it is compressed against, and saying so up front
/// beats failing partway through a multi-hour job.
#[test]
fn cram_without_a_reference_is_refused_before_any_work() {
    let dir = scratch("cramref");
    let err = crate::output::AlignmentWriter::create(&dir.join("out.cram"), OutputFormat::Cram, "@HD\tVN:1.6\n", None);
    assert!(err.is_err());
}
