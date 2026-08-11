//! Tests for paired-end mapping.
//!
//! Two things get the most attention because they are what this module adds over the single-end
//! path and what nothing else checks: the paired SAM fields (flags, `RNEXT`/`PNEXT`, `TLEN`), and
//! that pairing survives a split index — pairing decided per part would rest on a fraction of the
//! genome.

use std::path::{Path, PathBuf};

use super::*;
use crate::batch::BatchSize;
use crate::index::build_index;
use crate::output::OutputFormat;
use crate::preset::Preset;

fn scratch(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("dun-pe-{}-{tag}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

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

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => *other,
        })
        .collect()
}

/// Proper FR pairs lifted out of the reference: R1 forward at `start`, R2 reverse-complemented
/// from the far end of a `fragment`-length template. Names encode the truth.
fn write_pairs(
    dir: &Path,
    contigs: usize,
    len: usize,
    per_contig: usize,
    read_len: usize,
    fragment: usize,
) -> (PathBuf, PathBuf) {
    let p1 = dir.join("r1.fq");
    let p2 = dir.join("r2.fq");
    let (mut t1, mut t2) = (Vec::new(), Vec::new());
    for c in 0..contigs {
        let bases = reference_bases(c, len);
        for i in 0..per_contig {
            let start = (i * 7919) % (len - fragment - 1);
            let mate_start = start + fragment - read_len;
            let r1 = &bases[start..start + read_len];
            let r2 = revcomp(&bases[mate_start..mate_start + read_len]);

            t1.extend_from_slice(format!("@chr{c}_{start}\n").as_bytes());
            t1.extend_from_slice(r1);
            t1.extend_from_slice(b"\n+\n");
            t1.extend_from_slice(&vec![b'I'; read_len]);
            t1.push(b'\n');

            t2.extend_from_slice(format!("@chr{c}_{start}\n").as_bytes());
            t2.extend_from_slice(&r2);
            t2.extend_from_slice(b"\n+\n");
            t2.extend_from_slice(&vec![b'I'; read_len]);
            t2.push(b'\n');
        }
    }
    std::fs::write(&p1, t1).unwrap();
    std::fs::write(&p2, t2).unwrap();
    (p1, p2)
}

fn run(index: &Path, r1: &Path, r2: &Path, out: &Path, dir: &Path) -> MapStats {
    map_pairs(
        index,
        r1,
        r2,
        out,
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
        &|| false,
        &mut |_, _, _| {},
    )
    .expect("paired mapping should succeed")
}

struct Rec {
    qname: String,
    flag: u16,
    rname: String,
    pos: i64,
    rnext: String,
    pnext: i64,
    tlen: i64,
}

fn records(sam: &Path) -> Vec<Rec> {
    std::fs::read_to_string(sam)
        .unwrap()
        .lines()
        .filter(|l| !l.starts_with('@'))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            Rec {
                qname: f[0].into(),
                flag: f[1].parse().unwrap(),
                rname: f[2].into(),
                pos: f[3].parse().unwrap(),
                rnext: f[6].into(),
                pnext: f[7].parse().unwrap(),
                tlen: f[8].parse().unwrap(),
            }
        })
        .collect()
}

/// Primary records only — supplementary/secondary ones repeat a template and would double-count.
fn primaries(sam: &Path) -> Vec<Rec> {
    records(sam).into_iter().filter(|r| r.flag & 0x900 == 0).collect()
}

fn setup(dir: &Path, contigs: usize, len: usize, batch: BatchSize) -> (PathBuf, PathBuf, PathBuf) {
    let reference = write_reference(dir, contigs, len);
    let (r1, r2) = write_pairs(dir, contigs, len, 20, 150, 500);
    let index = dir.join("ref.mmi");
    build_index(&reference, &index, Preset::ShortRead, batch, &mut |_, _| {}).unwrap();
    (index, r1, r2)
}

// ---- the paired SAM fields ------------------------------------------------

/// The fields this module exists to add. Every one of them is something a single-end record does
/// not have and that downstream tools rely on to reconstruct the template.
#[test]
fn proper_pairs_get_paired_flags_mate_fields_and_opposing_tlen() {
    let dir = scratch("fields");
    let (index, r1, r2) = setup(&dir, 2, 300_000, BatchSize::default());
    let sam = dir.join("out.sam");
    let stats = run(&index, &r1, &r2, &sam, &dir);

    assert_eq!(stats.queries, 80, "40 templates, both ends");

    let recs = primaries(&sam);
    assert_eq!(recs.len(), 80);

    let mut proper = 0;
    for pair in recs.chunks(2) {
        let (a, b) = (&pair[0], &pair[1]);
        assert_eq!(a.qname, b.qname, "both ends share a QNAME");

        // Paired, and exactly one of first/last on each end.
        assert_ne!(a.flag & 0x1, 0, "0x1 paired");
        assert_ne!(b.flag & 0x1, 0);
        assert_ne!(a.flag & 0x40, 0, "R1 is the first segment");
        assert_ne!(b.flag & 0x80, 0, "R2 is the last segment");
        assert_eq!(a.flag & 0x80, 0, "R1 must not also claim last");
        assert_eq!(b.flag & 0x40, 0);

        if a.flag & 0x4 != 0 || b.flag & 0x4 != 0 {
            continue;
        }

        // Mate fields point at each other.
        assert_eq!(a.rnext, "=", "same contig is written as '='");
        assert_eq!(b.rnext, "=");
        assert_eq!(a.pnext, b.pos, "PNEXT is the mate's POS");
        assert_eq!(b.pnext, a.pos);

        // The strand bits are consistent between the two ends.
        assert_eq!(
            (a.flag & 0x10 != 0),
            (b.flag & 0x20 != 0),
            "mate-reverse mirrors reverse"
        );
        assert_eq!((b.flag & 0x10 != 0), (a.flag & 0x20 != 0));

        // Observed template length: equal magnitude, opposite sign.
        assert_eq!(a.tlen, -b.tlen, "TLEN must oppose across the pair");
        assert!(a.tlen.abs() > 0, "a mapped pair has a template length");
        if a.flag & 0x2 != 0 {
            proper += 1;
            assert_ne!(b.flag & 0x2, 0, "0x2 must agree across the pair");
            assert!(
                (a.tlen.abs() - 500).abs() <= 20,
                "TLEN {} should be near the 500 bp simulated fragment",
                a.tlen
            );
        }
    }
    assert!(proper >= 35, "expected most FR pairs to be proper, got {proper}");
}

/// An FR pair is one forward and one reverse read; if both came out on the same strand the
/// orientation handling is wrong and `proper_frag` would be meaningless.
#[test]
fn the_two_ends_map_to_opposite_strands() {
    let dir = scratch("strand");
    let (index, r1, r2) = setup(&dir, 2, 300_000, BatchSize::default());
    let sam = dir.join("out.sam");
    run(&index, &r1, &r2, &sam, &dir);

    let recs = primaries(&sam);
    let mut opposite = 0;
    for pair in recs.chunks(2) {
        if pair[0].flag & 0x4 != 0 || pair[1].flag & 0x4 != 0 {
            continue;
        }
        if (pair[0].flag & 0x10 != 0) != (pair[1].flag & 0x10 != 0) {
            opposite += 1;
        }
    }
    assert!(opposite >= 35, "expected FR orientation on most pairs, got {opposite}");
}

/// A read whose mate did not map still has to say so, and both records must stay at the same
/// locus so a coordinate sort keeps the template together.
#[test]
fn an_unmappable_mate_is_flagged_and_placed_with_its_partner() {
    let dir = scratch("halfmapped");
    let reference = write_reference(&dir, 1, 200_000);
    let index = dir.join("ref.mmi");
    build_index(
        &reference,
        &index,
        Preset::ShortRead,
        BatchSize::default(),
        &mut |_, _| {},
    )
    .unwrap();

    // R1 comes from the reference; R2 is poly-A and shares nothing with it.
    let bases = reference_bases(0, 200_000);
    let r1 = dir.join("r1.fq");
    let r2 = dir.join("r2.fq");
    let read: Vec<u8> = bases[5000..5150].to_vec();
    let mut t1 = b"@solo\n".to_vec();
    t1.extend_from_slice(&read);
    t1.extend_from_slice(b"\n+\n");
    t1.extend_from_slice(&[b'I'; 150]);
    t1.push(b'\n');
    std::fs::write(&r1, t1).unwrap();
    std::fs::write(&r2, format!("@solo\n{}\n+\n{}\n", "A".repeat(150), "I".repeat(150))).unwrap();

    let sam = dir.join("out.sam");
    let stats = run(&index, &r1, &r2, &sam, &dir);
    assert_eq!(stats.queries, 2);
    assert_eq!(stats.unmapped, 1, "the poly-A end does not map");

    let recs = records(&sam);
    let mapped = recs.iter().find(|r| r.flag & 0x4 == 0).expect("R1 maps");
    let unmapped = recs.iter().find(|r| r.flag & 0x4 != 0).expect("R2 does not");

    assert_ne!(mapped.flag & 0x8, 0, "the mapped end reports its mate unmapped");
    assert_eq!(
        unmapped.rname, mapped.rname,
        "the unmapped read is placed with its mate"
    );
    assert_eq!(unmapped.pos, mapped.pos);
    assert_eq!(mapped.flag & 0x2, 0, "a half-mapped template is not a proper pair");
}

/// Vendor FASTQ commonly carries `/1` and `/2`. Both ends must end up with the same QNAME or
/// every downstream tool loses the pairing.
#[test]
fn mate_suffixes_are_stripped_so_both_ends_share_a_qname() {
    assert_eq!(strip_mate_suffix("read/1"), "read");
    assert_eq!(strip_mate_suffix("read/2"), "read");
    assert_eq!(strip_mate_suffix("read"), "read");
    // Not a mate suffix — a name that merely ends in a digit, or in /3.
    assert_eq!(strip_mate_suffix("read1"), "read1");
    assert_eq!(strip_mate_suffix("read/3"), "read/3");
}

/// TLEN is signed by which end is leftmost, and zero when neither is.
#[test]
fn tlen_is_signed_by_position_and_zero_on_a_tie() {
    let reg = |rs: i32, re: i32| AlignReg {
        rs,
        re,
        ..Default::default()
    };
    assert_eq!(tlen(&reg(100, 250), &reg(400, 550)), 450, "leftmost is positive");
    assert_eq!(tlen(&reg(400, 550), &reg(100, 250)), -450, "rightmost is negative");
    assert_eq!(tlen(&reg(100, 250), &reg(100, 250)), 0, "no leftmost end");
}

// ---- the split index ------------------------------------------------------

/// The same claim the single-end path makes, for pairs: splitting the index is a memory decision,
/// not an accuracy one. Pairing in particular must be re-derived after the merge — done per part
/// it would rest on whichever fraction of the genome was resident.
#[test]
fn a_split_index_pairs_reads_exactly_as_a_whole_index_does() {
    let dir = scratch("equivalence");
    let (contigs, len) = (6, 500_000);
    let reference = write_reference(&dir, contigs, len);
    let (r1, r2) = write_pairs(&dir, contigs, len, 20, 150, 500);

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
    let whole_stats = run(&whole, &r1, &r2, &whole_sam, &dir);

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
    let split_stats = run(&split_idx, &r1, &r2, &split_sam, &dir);

    assert_eq!(whole_stats.parts, 1);
    assert!(split_stats.parts > 1, "the split index must actually split");
    assert_eq!(split_stats.mapped, whole_stats.mapped);

    let (a, b) = (records(&split_sam), records(&whole_sam));
    assert_eq!(a.len(), b.len(), "same number of records");
    for (x, y) in a.iter().zip(&b) {
        assert_eq!(
            (&x.qname, x.flag, &x.rname, x.pos, &x.rnext, x.pnext, x.tlen),
            (&y.qname, y.flag, &y.rname, y.pos, &y.rnext, y.pnext, y.tlen),
            "split and whole disagree for {}",
            x.qname
        );
    }
}

/// R1/R2 that have drifted out of step would pair every later read with the wrong mate — a
/// corruption that produces confident, wrong alignments. It must refuse rather than truncate.
#[test]
fn mismatched_read_counts_are_refused() {
    let dir = scratch("lockstep");
    let reference = write_reference(&dir, 1, 200_000);
    let index = dir.join("ref.mmi");
    build_index(
        &reference,
        &index,
        Preset::ShortRead,
        BatchSize::default(),
        &mut |_, _| {},
    )
    .unwrap();

    let bases = reference_bases(0, 200_000);
    let mut t1 = Vec::new();
    for i in 0..3 {
        let s = 1000 + i * 500;
        t1.extend_from_slice(format!("@p{i}\n").as_bytes());
        t1.extend_from_slice(&bases[s..s + 150]);
        t1.extend_from_slice(b"\n+\n");
        t1.extend_from_slice(&[b'I'; 150]);
        t1.push(b'\n');
    }
    let r1 = dir.join("r1.fq");
    let r2 = dir.join("r2.fq");
    std::fs::write(&r1, t1).unwrap();
    // Only one mate for three R1 reads.
    let mut t2 = b"@p0\n".to_vec();
    t2.extend_from_slice(&revcomp(&bases[1350..1500]));
    t2.extend_from_slice(b"\n+\n");
    t2.extend_from_slice(&[b'I'; 150]);
    t2.push(b'\n');
    std::fs::write(&r2, t2).unwrap();

    let err = map_pairs(
        &index,
        &r1,
        &r2,
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
        &|| false,
        &mut |_, _, _| {},
    );
    assert!(err.is_err(), "out-of-step mates must not be paired silently");
}

#[test]
fn cancellation_stops_paired_mapping() {
    let dir = scratch("cancel");
    let (index, r1, r2) = setup(&dir, 1, 200_000, BatchSize::default());
    let err = map_pairs(
        &index,
        &r1,
        &r2,
        &dir.join("out.sam"),
        &dir.join("scratch"),
        &MapParams::default(),
        &|| true,
        &mut |_, _, _| {},
    )
    .unwrap_err();
    assert!(matches!(err, AlignError::Cancelled));
}

/// **Regression.** Every fixture above uses fixed-length reads, and that is exactly why they all
/// passed while a real WGS failed 74 minutes into a run.
///
/// The paired reader batches by *bases*. If the two files are read independently with the same
/// base budget, files whose reads differ in length return different record counts — and real data
/// always has a tail of shorter reads from adapter and quality trimming. On WGS229 that surfaced
/// as 332,653 against 332,722 in one batch, tripping the lockstep guard on files that were
/// perfectly in step.
///
/// So this fixture deliberately gives R1 and R2 *different* length distributions, and asserts
/// every pair still comes back matched.
#[test]
fn pairs_stay_in_step_when_reads_have_different_lengths() {
    let dir = scratch("varlen");
    let contig_len = 200_000;
    let reference = write_reference(&dir, 1, contig_len);
    let bases = reference_bases(0, contig_len);

    let (p1, p2) = (dir.join("r1.fq"), dir.join("r2.fq"));
    let (mut t1, mut t2) = (Vec::new(), Vec::new());
    for i in 0..400usize {
        let start = 1000 + i * 200;
        // R1 keeps full length; R2 is trimmed by a varying amount, as a trimmer would leave it.
        let len1 = 150;
        let len2 = 150 - (i % 37);
        let r1 = &bases[start..start + len1];
        let r2 = revcomp(&bases[start + 300 - len2..start + 300]);

        t1.extend_from_slice(format!("@p{i}\n").as_bytes());
        t1.extend_from_slice(r1);
        t1.extend_from_slice(b"\n+\n");
        t1.extend_from_slice(&vec![b'I'; len1]);
        t1.push(b'\n');

        t2.extend_from_slice(format!("@p{i}\n").as_bytes());
        t2.extend_from_slice(&r2);
        t2.extend_from_slice(b"\n+\n");
        t2.extend_from_slice(&vec![b'I'; len2]);
        t2.push(b'\n');
    }
    std::fs::write(&p1, t1).unwrap();
    std::fs::write(&p2, t2).unwrap();

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
    let stats = run(&index, &p1, &p2, &sam, &dir);

    assert_eq!(stats.queries, 800, "400 templates, both ends, none dropped");

    // Both ends of every template must be present and share a QNAME.
    let recs = primaries(&sam);
    assert_eq!(recs.len(), 800);
    for pair in recs.chunks(2) {
        assert_eq!(pair[0].qname, pair[1].qname, "a template's ends were separated");
        assert_ne!(pair[0].flag & 0x40, 0);
        assert_ne!(pair[1].flag & 0x80, 0);
    }
}

/// The guard still has to fire when the files really are mismatched, which is the case it exists
/// for — one file ending before the other.
#[test]
fn a_truncated_mate_file_is_still_refused() {
    let dir = scratch("truncated");
    let contig_len = 200_000;
    let reference = write_reference(&dir, 1, contig_len);
    let bases = reference_bases(0, contig_len);

    let (p1, p2) = (dir.join("r1.fq"), dir.join("r2.fq"));
    let mut t1 = Vec::new();
    for i in 0..5usize {
        let start = 1000 + i * 400;
        t1.extend_from_slice(format!("@p{i}\n").as_bytes());
        t1.extend_from_slice(&bases[start..start + 150]);
        t1.extend_from_slice(b"\n+\n");
        t1.extend_from_slice(&[b'I'; 150]);
        t1.push(b'\n');
    }
    std::fs::write(&p1, t1).unwrap();
    // Only two mates for five R1 reads.
    let mut t2 = Vec::new();
    for i in 0..2usize {
        let start = 1000 + i * 400;
        t2.extend_from_slice(format!("@p{i}\n").as_bytes());
        t2.extend_from_slice(&revcomp(&bases[start + 250..start + 400]));
        t2.extend_from_slice(b"\n+\n");
        t2.extend_from_slice(&[b'I'; 150]);
        t2.push(b'\n');
    }
    std::fs::write(&p2, t2).unwrap();

    let index = dir.join("ref.mmi");
    build_index(
        &reference,
        &index,
        Preset::ShortRead,
        BatchSize::default(),
        &mut |_, _| {},
    )
    .unwrap();
    let err = map_pairs(
        &index,
        &p1,
        &p2,
        &dir.join("out.sam"),
        &dir.join("scratch"),
        &MapParams {
            preset: Preset::ShortRead,
            threads: 1,
            read_group: None,
            format: OutputFormat::Sam,
            reference: None,
        },
        &|| false,
        &mut |_, _, _| {},
    );
    assert!(err.is_err(), "a truncated mate file must not be paired silently");
}
