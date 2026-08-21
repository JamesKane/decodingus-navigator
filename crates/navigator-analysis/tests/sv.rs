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
        None, // BAM needs no reference
        &lengths,
        400.0, // expected insert
        50.0,  // sd  -> outlier if insert > 600 or < 200
        &SvCallerConfig::default(),
        &navigator_analysis::CancelToken::none(),
    )
    .expect("walker should succeed");

    // 2 pairs across two chromosomes, one for each mate, and 2 pairs whose insert size is an
    // outlier.
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
        .any(|p| &*p.chrom1 == "chr1" && &*p.chrom2 == "chr2" && p.pos1 == 100));

    // One split read with 20 bp clip, supplementary on chr1:2000.
    assert_eq!(ev.total_split_reads(), 1);
    let sr = &ev.split_reads[0];
    assert_eq!(sr.clip_length, 20);
    assert_eq!(&*sr.supp_chrom, "chr1");
    assert_eq!(sr.supp_pos, 2000);

    // Depth bins.
    assert_eq!(ev.depth_bins["chr1"], vec![2, 1, 0, 0, 1]);
    assert_eq!(ev.depth_bins["chr2"], vec![1, 0, 0, 0, 0]);
}

/// The same reads, stored as a CRAM, must give evidence that matches to the last byte. The walker
/// once opened every file with the BAM reader, which reads BGZF. So a CRAM failed at the open,
/// with `invalid BGZF header`. Every CRAM in a workspace then had no SV call at all, and nobody saw
/// it.
#[test]
fn walker_reads_cram_with_the_same_result_as_bam() {
    let lengths = BTreeMap::from([("chr1".to_string(), 5000i64), ("chr2".to_string(), 5000)]);
    let config = SvCallerConfig::default();
    let cfg = (400.0, 50.0);

    let from_bam = walker::collect_evidence(
        &fixtures().join("sv.bam"),
        None,
        &lengths,
        cfg.0,
        cfg.1,
        &config,
        &navigator_analysis::CancelToken::none(),
    )
    .expect("BAM walk should succeed");
    let from_cram = walker::collect_evidence(
        &fixtures().join("sv.cram"),
        Some(&fixtures().join("svref.fa")),
        &lengths,
        cfg.0,
        cfg.1,
        &config,
        &navigator_analysis::CancelToken::none(),
    )
    .expect("CRAM walk should succeed");

    assert_eq!(from_cram.depth_bins, from_bam.depth_bins);
    assert_eq!(from_cram.total_discordant_pairs(), from_bam.total_discordant_pairs());
    assert_eq!(from_cram.total_split_reads(), from_bam.total_split_reads());

    // Compare the evidence itself, and not the counts alone. The split read carries the fields
    // that come from the accessors that a CRAM implements in a different way. Those are the SA tag
    // and the clip length of the CIGAR.
    let (b, c) = (&from_bam.split_reads[0], &from_cram.split_reads[0]);
    assert_eq!(
        (c.clip_length, &c.supp_chrom, c.supp_pos, c.primary_pos),
        (b.clip_length, &b.supp_chrom, b.supp_pos, b.primary_pos)
    );
    let placed = |e: &SvEvidenceCollection| {
        let mut v: Vec<_> = e
            .discordant_pairs
            .iter()
            .map(|p| (p.chrom1.clone(), p.pos1, p.pos2, format!("{:?}", p.reason)))
            .collect();
        v.sort();
        v
    };
    assert_eq!(placed(&from_cram), placed(&from_bam));
}

/// The parallel fan-out over the contigs must give speed and nothing else. The depth bins, the
/// discordant pairs and the split reads must all match, and they must come in the same order.
///
/// SV was the last analysis over the whole genome that still decoded on one thread, at 2 to 5 h for
/// each 30x CRAM in a batch. So the fan-out is worth its code only when it changes the wall time
/// and nothing more.
///
/// This runs over a BAM *and* over a CRAM. The two take different decode paths into the same
/// sink.
#[test]
fn parallel_walk_matches_sequential_on_bam_and_cram() {
    let lengths = BTreeMap::from([("chr1".to_string(), 5000i64), ("chr2".to_string(), 5000)]);
    let config = SvCallerConfig::default();
    let cases = [("sv.bam", None), ("sv.cram", Some(fixtures().join("svref.fa")))];

    for (file, reference) in cases {
        let path = fixtures().join(file);
        // With no index, the parallel entry point falls back to the sequential walk. This test
        // would then compare a walk against itself, and it would pass whatever the fan-out
        // does.
        assert!(
            navigator_analysis::reader::has_region_index(&path),
            "{file} needs its .bai/.crai for this test to exercise the parallel path"
        );

        let seq = walker::collect_evidence(
            &path,
            reference.as_deref(),
            &lengths,
            400.0,
            50.0,
            &config,
            &navigator_analysis::CancelToken::none(),
        )
        .unwrap_or_else(|e| panic!("{file} sequential walk: {e}"));
        let par = walker::collect_evidence_parallel(
            &path,
            reference.as_deref(),
            &lengths,
            400.0,
            50.0,
            &config,
            &navigator_analysis::CancelToken::none(),
        )
        .unwrap_or_else(|e| panic!("{file} parallel walk: {e}"));

        assert_eq!(par.depth_bins, seq.depth_bins, "{file} depth bins");
        assert_eq!(par.discordant_pairs, seq.discordant_pairs, "{file} discordant pairs");
        assert_eq!(par.split_reads, seq.split_reads, "{file} split reads");
    }
}

/// The evidence cap must cut short what the code *keeps*, and it must not change what the code
/// *found*. A cap that lowered `total_discordant_pairs`, where nobody saw it, would make a run that
/// the cap cut short read as a clean sample. That is the one failure that a safety valve must not
/// have. This test covers both walks, because each one claims against the shared budget on its
/// own.
#[test]
fn evidence_cap_truncates_retained_evidence_but_not_the_reported_totals() {
    let lengths = BTreeMap::from([("chr1".to_string(), 5000i64), ("chr2".to_string(), 5000)]);
    let capped = SvCallerConfig {
        max_evidence_records: 1,
        ..SvCallerConfig::default()
    };

    for parallel in [false, true] {
        let walk = if parallel {
            walker::collect_evidence_parallel
        } else {
            walker::collect_evidence
        };
        let ev = walk(
            &fixtures().join("sv.bam"),
            None,
            &lengths,
            400.0,
            50.0,
            &capped,
            &navigator_analysis::CancelToken::none(),
        )
        .expect("capped walk should succeed");

        // 4 discordant pairs exist (see the uncapped test); the cap keeps 1 and counts the rest.
        assert_eq!(ev.discordant_pairs.len(), 1, "parallel={parallel} retained");
        assert_eq!(ev.discordant_pairs_dropped, 3, "parallel={parallel} dropped");
        assert_eq!(ev.total_discordant_pairs(), 4, "parallel={parallel} reported total");
        // There is only 1 split read. So its count reaches the cap exactly, and the code drops
        // nothing.
        assert_eq!(ev.split_reads.len(), 1, "parallel={parallel} split retained");
        assert_eq!(ev.split_reads_dropped, 0, "parallel={parallel} split dropped");
        // Depth bins are not evidence records and are never capped.
        assert_eq!(ev.depth_bins["chr1"], vec![2, 1, 0, 0, 1], "parallel={parallel} depth");
    }
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
        discordant_pairs_dropped: 0,
        split_reads_dropped: 0,
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
