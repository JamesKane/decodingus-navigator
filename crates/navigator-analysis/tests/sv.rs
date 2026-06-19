//! Structural variant tests: the walker against sv.bam, plus the pure segmenter /
//! clusterer / confidence logic against hand-built evidence.

use std::collections::BTreeMap;
use std::path::PathBuf;

use navigator_analysis::sv::evidence::{DepthSegment, DiscordantPair, DiscordantReason, SvEvidenceCollection};
use navigator_analysis::sv::types::SvCall;
use navigator_analysis::sv::{calculate_confidence, clusterer, segmenter, walker, SvCallerConfig, SvType};

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

// ---- walker (sv.bam) -------------------------------------------------------

#[test]
fn walker_extracts_discordant_pairs_split_reads_and_depth() {
    let lengths = BTreeMap::from([("chr1".to_string(), 5000i64), ("chr2".to_string(), 5000)]);
    let ev = walker::collect_evidence(
        &fixtures().join("sv.bam"),
        &lengths,
        400.0, // expected insert
        50.0,  // sd  -> outlier if insert > 600 or < 200
        &SvCallerConfig::default(),
    )
    .expect("walker should succeed");

    // 2 inter-chromosomal (one per mate) + 2 insert-size outliers.
    assert_eq!(ev.total_discordant_pairs(), 4);
    let inter = ev.inter_chromosomal_pairs();
    assert_eq!(inter.len(), 2);
    let outliers = ev
        .discordant_pairs
        .iter()
        .filter(|p| p.reason == DiscordantReason::InsertSizeOutlier)
        .count();
    assert_eq!(outliers, 2);
    // One inter pair is chr1->chr2 from r_inter.
    assert!(inter
        .iter()
        .any(|p| p.chrom1 == "chr1" && p.chrom2 == "chr2" && p.pos1 == 100));

    // One split read with 20 bp clip, supplementary on chr1:2000.
    assert_eq!(ev.total_split_reads(), 1);
    let sr = &ev.split_reads[0];
    assert_eq!(sr.clip_length, 20);
    assert_eq!(sr.supp_chrom, "chr1");
    assert_eq!(sr.supp_pos, 2000);

    // Depth bins.
    assert_eq!(ev.depth_bins["chr1"], vec![2, 1, 0, 0, 1]);
    assert_eq!(ev.depth_bins["chr2"], vec![1, 0, 0, 0, 0]);
}

// ---- segmenter (pure) ------------------------------------------------------

#[test]
fn segmenter_calls_del_and_dup_and_applies_size_filter() {
    // expected reads/bin = 30 * 1000 / 150 = 200.
    let mut bins = vec![200u32; 60];
    for b in bins.iter_mut().take(30).skip(10) {
        *b = 20; // bins 10..29: deletion (20 bins = 20 kb)
    }
    for b in bins.iter_mut().take(50).skip(40) {
        *b = 400; // bins 40..49: duplication (10 bins = 10 kb)
    }
    bins[55] = 0; // single aberrant bin (1 kb) -> filtered by min_cnv_size

    let depth_bins = BTreeMap::from([("chr1".to_string(), bins)]);
    let lengths = BTreeMap::from([("chr1".to_string(), 60_000i64)]);
    let segs = segmenter::segment(&depth_bins, &lengths, 30.0, 150.0, &SvCallerConfig::default());

    assert_eq!(segs.len(), 2, "got {segs:?}");
    assert_eq!(segs[0].sv_type, SvType::Del);
    assert_eq!((segs[0].start, segs[0].end), (10_000, 30_000));
    assert_eq!(segs[0].num_bins, 20);
    assert!(segs[0].z_score < 0.0);
    assert_eq!(segs[1].sv_type, SvType::Dup);
    assert_eq!((segs[1].start, segs[1].end), (40_000, 50_000));
    assert_eq!(segs[1].num_bins, 10);
    assert!((segs[1].log2_ratio - 1.0).abs() < 1e-9); // 400/200 = 2x -> log2 = 1
}

#[test]
fn merge_nearby_segments_joins_same_type_within_gap() {
    let del = |start, end, bins| DepthSegment {
        chrom: "chr1".into(),
        start,
        end,
        mean_depth: 20.0,
        log2_ratio: -3.0,
        z_score: -10.0,
        num_bins: bins,
        sv_type: SvType::Del,
    };
    let segs = vec![del(0, 10_000, 10), del(20_000, 30_000, 10)]; // gap 10kb <= 50kb
    let merged = segmenter::merge_nearby_segments(&segs, 50_000);
    assert_eq!(merged.len(), 1);
    assert_eq!((merged[0].start, merged[0].end), (0, 30_000));
    assert_eq!(merged[0].num_bins, 20);
}

// ---- clusterer (pure) ------------------------------------------------------

fn pair(pos1: i64, pos2: i64, s1: char, s2: char, reason: DiscordantReason) -> DiscordantPair {
    DiscordantPair {
        read_name: "r".into(),
        chrom1: "chr1".into(),
        pos1,
        strand1: s1,
        chrom2: if reason == DiscordantReason::InterChromosomal {
            "chr2".into()
        } else {
            "chr1".into()
        },
        pos2,
        strand2: s2,
        insert_size: 6000,
        mapq: 60,
        reason,
    }
}

fn collection(pairs: Vec<DiscordantPair>) -> SvEvidenceCollection {
    SvEvidenceCollection {
        discordant_pairs: pairs,
        split_reads: Vec::new(),
        depth_bins: BTreeMap::new(),
        sample_name: "test".into(),
        expected_insert_size: 400.0,
        insert_size_sd: 50.0,
    }
}

#[test]
fn clusterer_calls_deletion_from_fr_insert_outliers() {
    let pairs = vec![
        pair(1000, 2000, '+', '-', DiscordantReason::InsertSizeOutlier),
        pair(1100, 2100, '+', '-', DiscordantReason::InsertSizeOutlier),
        pair(1200, 2200, '+', '-', DiscordantReason::InsertSizeOutlier),
    ];
    let calls = clusterer::cluster(&collection(pairs), &[], &SvCallerConfig::default());
    assert_eq!(calls.len(), 1);
    let c = &calls[0];
    assert_eq!(c.sv_type, SvType::Del);
    assert_eq!(c.start, 1100); // mean of pos1
    assert_eq!(c.end, 2100); // start + |meanMate - start|
    assert_eq!(c.sv_len, -1000);
    assert_eq!(c.paired_end_support, 3);
    assert_eq!(c.filter, "PASS");
}

#[test]
fn clusterer_calls_inversion_for_same_strand_pairs() {
    let pairs = vec![
        pair(1000, 2000, '+', '+', DiscordantReason::WrongOrientation),
        pair(1100, 2100, '+', '+', DiscordantReason::WrongOrientation),
        pair(1200, 2200, '+', '+', DiscordantReason::WrongOrientation),
    ];
    let calls = clusterer::cluster(&collection(pairs), &[], &SvCallerConfig::default());
    assert_eq!(calls.len(), 1);
    assert_eq!(calls[0].sv_type, SvType::Inv);
}

#[test]
fn clusterer_calls_translocation_for_interchromosomal_pairs() {
    let pairs = vec![
        pair(500, 5000, '+', '-', DiscordantReason::InterChromosomal),
        pair(600, 5100, '+', '-', DiscordantReason::InterChromosomal),
    ];
    let calls = clusterer::cluster(&collection(pairs), &[], &SvCallerConfig::default());
    assert_eq!(calls.len(), 1);
    let c = &calls[0];
    assert_eq!(c.sv_type, SvType::Bnd);
    assert_eq!(c.start, 550); // mean pos1
    assert_eq!(c.mate_chrom.as_deref(), Some("chr2"));
    assert_eq!(c.mate_pos, Some(5050)); // mean pos2
}

#[test]
fn confidence_weights_pe_sr_and_depth() {
    let call = SvCall {
        id: "x".into(),
        chrom: "chr1".into(),
        start: 1,
        end: 2,
        sv_type: SvType::Del,
        sv_len: -1,
        ci_pos: (0, 0),
        ci_end: (0, 0),
        quality: 50.0,
        paired_end_support: 10,    // -> 1.0
        split_read_support: 5,     // -> 1.0
        relative_depth: Some(0.5), // deviation 0.5 -> 1.0
        mate_chrom: None,
        mate_pos: None,
        filter: "PASS".into(),
        genotype: "0/1".into(),
    };
    assert!((calculate_confidence(&call) - 1.0).abs() < 1e-9);

    let mut weak = call.clone();
    weak.paired_end_support = 0;
    weak.split_read_support = 0;
    weak.relative_depth = None;
    assert!((calculate_confidence(&weak) - 0.0).abs() < 1e-9);
}
