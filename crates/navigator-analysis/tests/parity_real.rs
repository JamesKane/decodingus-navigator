//! A smoke test on real data. It is the seed of the parity harness of phase 3, §4c. The test
//! carries `#[ignore]`, and it runs only when the environment variables point it at a local BAM
//! and a reference:
//!
//! ```text
//! HG002_CHRM_BAM=/tmp/hg002.chrM.bam CHM13_REF=/Users/.../chm13v2.0.fa \
//!   cargo test -p navigator-analysis --test parity_real -- --ignored --nocapture
//! ```
//!
//! This test checks two things. noodles handles a real BAM, which carries CIGARs and MAPQ values
//! of every shape. And the chrM numbers look reasonable. It is NOT a strict parity check. Phase 3
//! measures that against the Scala walker and against GATK.

use std::collections::HashSet;
use std::path::PathBuf;

use navigator_analysis::caller::{call_denovo, HaploidCallerParams};
use navigator_analysis::coverage::{collect_coverage_callable, CallableLociParams};
use navigator_analysis::parity::{compare_denovo_snps, parse_truth_vcf};
use navigator_analysis::read_metrics::collect_read_metrics;
use navigator_analysis::sex::{infer_from_bam, InferredSex};
use navigator_analysis::unified::{collect_unified_metrics, collect_unified_metrics_parallel};

fn real_data() -> Option<(String, String)> {
    match (std::env::var("HG002_CHRM_BAM"), std::env::var("CHM13_REF")) {
        (Ok(bam), Ok(reference)) => Some((bam, reference)),
        _ => {
            eprintln!("set HG002_CHRM_BAM and CHM13_REF to run this test");
            None
        }
    }
}

#[test]
#[ignore = "requires local HG002_CHRM_BAM + CHM13_REF env vars"]
fn hg002_chrm_smoke() {
    let Some((bam, reference)) = real_data() else { return };

    let allow: HashSet<String> = ["chrM".to_string()].into_iter().collect();
    let result = collect_coverage_callable(
        &PathBuf::from(bam),
        &PathBuf::from(reference),
        &CallableLociParams::default(),
        Some(&allow),
    )
    .expect("walker should succeed on real data");

    eprintln!("genome_territory = {}", result.genome_territory);
    eprintln!("mean_coverage    = {:.3}", result.mean_coverage);
    eprintln!("median_coverage  = {}", result.median_coverage);
    eprintln!("sd_coverage      = {:.3}", result.sd_coverage);
    eprintln!(
        "pct_10x/20x/30x  = {:.4} / {:.4} / {:.4}",
        result.pct_10x, result.pct_20x, result.pct_30x
    );
    eprintln!("callable_bases   = {}", result.callable_bases);
    eprintln!("callable metrics = {:?}", result.contig_callable);
    eprintln!("coverage stats   = {:?}", result.contig_coverage_stats);

    // The reads must cover the whole of chrM, at a high depth.
    assert_eq!(result.genome_territory, 16569);
    let cs = &result.contig_coverage_stats[0];
    assert_eq!(cs.contig, "chrM");
    assert_eq!(cs.cov_bases, 16569);
    assert!((cs.coverage - 100.0).abs() < 1e-9);
    assert!(cs.mean_depth > 50.0, "mean depth {} unexpectedly low", cs.mean_depth);

    let cm = &result.contig_callable[0];
    let total =
        cm.callable + cm.low_coverage + cm.no_coverage + cm.poor_mapping_quality + cm.ref_n + cm.excessive_coverage;
    assert_eq!(total, 16569);
}

#[test]
#[ignore = "requires local HG002_CHRM_BAM + CHM13_REF env vars"]
fn hg002_chrm_denovo_smoke() {
    let Some((bam, reference)) = real_data() else { return };

    let calls = call_denovo(
        &PathBuf::from(bam),
        &PathBuf::from(reference),
        "chrM",
        &HaploidCallerParams::default(),
        &navigator_analysis::CancelToken::none(),
    )
    .expect("de-novo should succeed on real data");

    eprintln!("chrM de-novo SNP calls: {}", calls.len());
    for c in calls.iter().take(15) {
        eprintln!(
            "  chrM:{} {}>{} depth={} af={:.3}",
            c.position, c.reference_allele, c.alternate_allele, c.depth, c.allele_fraction
        );
    }

    // The mtDNA of HG002 against the chrM of CHM13. There are a few real differences, from a
    // handful up to a few dozen, at a high depth. There are never thousands. Thousands would mean
    // that a gate no longer works: either the one on the consensus, or the one on the
    // fraction.
    assert!(!calls.is_empty(), "expected some mtDNA variants");
    assert!(calls.len() < 1000, "implausibly many calls: {}", calls.len());
    for c in &calls {
        assert!(c.allele_fraction >= 0.5);
        assert!(c.depth >= 4);
        assert_ne!(c.reference_allele, c.alternate_allele);
    }
}

/// A de-novo call over chrY, at 57 Mb, against the full BAM. It covers the tally that works in
/// chunks, on a large contig. The chunk bounds the memory, and the length of chrY does not.
#[test]
#[ignore = "requires HG002_BAM + CHM13_REF (chrY de-novo, chunked)"]
fn hg002_chry_denovo_streams() {
    let (Ok(bam), Ok(reference)) = (std::env::var("HG002_BAM"), std::env::var("CHM13_REF")) else {
        eprintln!("set HG002_BAM and CHM13_REF to run this test");
        return;
    };
    let calls = call_denovo(
        &PathBuf::from(bam),
        &PathBuf::from(reference),
        "chrY",
        &HaploidCallerParams::default(),
        &navigator_analysis::CancelToken::none(),
    )
    .expect("chrY de-novo should succeed");
    eprintln!("chrY de-novo calls: {}", calls.len());
    // Shallow single lane -> few high-AF haploid calls; just assert it completed sanely.
    for c in &calls {
        assert!(c.allele_fraction >= 0.5 && c.depth >= 4);
    }
}

/// The coverage over the whole genome, from the full BAM, with no allowlist. This is possible only
/// because the walker streams a window that slides. The dense version before it allocated an array
/// at each position. It held one for every contig of the main assembly at one time, and that came
/// to about 84 GB.
#[test]
#[ignore = "requires HG002_BAM + CHM13_REF (whole-genome streaming coverage)"]
fn hg002_wgs_coverage_streams_all_contigs() {
    let (Ok(bam), Ok(reference)) = (std::env::var("HG002_BAM"), std::env::var("CHM13_REF")) else {
        eprintln!("set HG002_BAM and CHM13_REF to run this test");
        return;
    };
    let result = collect_coverage_callable(
        &PathBuf::from(bam),
        &PathBuf::from(reference),
        &CallableLociParams::default(),
        None,
    )
    .expect("whole-genome coverage should succeed");

    eprintln!(
        "genome_territory={} mean_coverage={:.4} contigs={}",
        result.genome_territory,
        result.mean_coverage,
        result.contig_coverage_stats.len()
    );
    // CHM13 main assembly: chr1-22, X, Y, M = 25 contigs, ~3.05 Gb.
    assert_eq!(result.contig_coverage_stats.len(), 25);
    assert!(
        result.genome_territory > 3_000_000_000,
        "territory {}",
        result.genome_territory
    );
    assert_eq!(result.genome_territory, result.coverage_histogram.iter().sum::<u64>());
    assert!(result.mean_coverage > 0.0);
}

/// Whole-BAM smoke tests for sex + read_metrics. Driven by HG002_BAM.
#[test]
#[ignore = "requires local HG002_BAM env var (full whole-genome BAM)"]
fn hg002_sex_smoke() {
    let Ok(bam) = std::env::var("HG002_BAM") else {
        eprintln!("set HG002_BAM to run this test");
        return;
    };
    let r = infer_from_bam(&PathBuf::from(bam), None).expect("sex inference should succeed");
    eprintln!(
        "inferred {:?} ({:?}); ratio={:.3} autosome={:.1} chrX={:.1}",
        r.inferred_sex, r.confidence, r.x_autosome_ratio, r.autosome_mean_coverage, r.x_coverage
    );
    // HG002 / NA24385 is male.
    assert_eq!(r.inferred_sex, InferredSex::Male);
}

#[test]
#[ignore = "requires local HG002_BAM env var (full whole-genome BAM scan)"]
fn hg002_read_metrics_smoke() {
    let Ok(bam) = std::env::var("HG002_BAM") else {
        eprintln!("set HG002_BAM to run this test");
        return;
    };
    let m = collect_read_metrics(&PathBuf::from(bam), None).expect("read metrics should succeed");
    eprintln!(
        "total={} pf_aligned={} ({:.3}) proper={:.3} chimera={:.4} orient={} \
         read_len mean={:.1} median={} insert mean={:.1} median={} mapq={:.1}",
        m.total_reads,
        m.pf_reads_aligned,
        m.pct_pf_reads_aligned,
        m.pct_proper_pairs,
        m.pct_chimeras,
        m.pair_orientation.as_str(),
        m.mean_read_length,
        m.median_read_length,
        m.mean_insert_size,
        m.median_insert_size,
        m.mean_mapping_quality
    );
    assert!(m.total_reads > 0);
    assert!(m.pct_pf_reads_aligned > 0.5);
    assert!(m.proper_pairs > 0);
}

/// In one pass over the whole genome, the fused unified walker must give exactly what the three
/// separate walkers give in three passes. That covers the coverage, the read-metrics and the sex,
/// field for field. This test guards the parity of the one-pass change, on real data.
#[test]
#[ignore = "requires HG002_BAM + CHM13_REF (whole-genome single-pass parity)"]
fn hg002_unified_matches_standalone() {
    let (Ok(bam), Ok(reference)) = (std::env::var("HG002_BAM"), std::env::var("CHM13_REF")) else {
        eprintln!("set HG002_BAM and CHM13_REF to run this test");
        return;
    };
    let bam = PathBuf::from(bam);
    let reference = PathBuf::from(reference);
    let params = CallableLociParams::default();

    let unified = collect_unified_metrics(&bam, &reference, &params, None).expect("unified walker should succeed");

    let cov = collect_coverage_callable(&bam, &reference, &params, None).unwrap();
    let rm = collect_read_metrics(&bam, Some(&reference)).unwrap();
    let sex = infer_from_bam(&bam, None).unwrap();

    assert_eq!(unified.coverage, cov, "coverage diverged from standalone");
    assert_eq!(unified.read_metrics, rm, "read metrics diverged from standalone");
    // The separate sex path on a BAM uses the fast path over the BAI. The fused path tallies from
    // the record stream. The two put the sample into the same class. But the read-density floats
    // can differ in the last ULP. That happens when the mapped-record count of the BAI differs
    // from the mapped-read count of the stream. So compare the call itself.
    let fused_sex = unified.sex.expect("HG002 has autosomes + chrX");
    assert_eq!(fused_sex.inferred_sex, sex.inferred_sex, "sex call diverged");
    assert_eq!(fused_sex.inferred_sex, InferredSex::Male);
}

/// The parallel walker over the contigs must give a result that matches the sequential one to the
/// last byte. It needs an indexed BAM. Without one it falls back to the sequential path, and the
/// test then proves nothing.
#[test]
#[ignore = "requires indexed HG002_BAM + CHM13_REF (parallel == sequential parity)"]
fn hg002_unified_parallel_matches_sequential() {
    let (Ok(bam), Ok(reference)) = (std::env::var("HG002_BAM"), std::env::var("CHM13_REF")) else {
        eprintln!("set HG002_BAM and CHM13_REF to run this test");
        return;
    };
    let bam = PathBuf::from(bam);
    let reference = PathBuf::from(reference);
    let params = CallableLociParams::default();

    let seq = collect_unified_metrics(&bam, &reference, &params, None).unwrap();
    let par = collect_unified_metrics_parallel(&bam, &reference, &params, None).unwrap();

    // The whole struct must be equal: the coverage, which holds the values of each contig and the
    // histogram, the read-metrics, and the sex.
    assert_eq!(par.coverage, seq.coverage, "parallel coverage diverged");
    assert_eq!(par.read_metrics, seq.read_metrics, "parallel read-metrics diverged");
    assert_eq!(par.sex, seq.sex, "parallel sex diverged");
}

/// A comparison of the wall time of the sequential unified metrics against the parallel ones,
/// which run one task at each contig. It asserts nothing about correctness, and the parity test
/// covers that. It prints both times, and the speedup.
#[test]
#[ignore = "perf smoke: requires indexed HG002_BAM + CHM13_REF"]
fn hg002_unified_parallel_timing() {
    let (Ok(bam), Ok(reference)) = (std::env::var("HG002_BAM"), std::env::var("CHM13_REF")) else {
        eprintln!("set HG002_BAM and CHM13_REF to run this test");
        return;
    };
    let bam = PathBuf::from(bam);
    let reference = PathBuf::from(reference);
    let params = CallableLociParams::default();

    let t0 = std::time::Instant::now();
    let seq = collect_unified_metrics(&bam, &reference, &params, None).unwrap();
    let seq_secs = t0.elapsed().as_secs_f64();

    let t1 = std::time::Instant::now();
    let par = collect_unified_metrics_parallel(&bam, &reference, &params, None).unwrap();
    let par_secs = t1.elapsed().as_secs_f64();

    assert_eq!(par.coverage.genome_territory, seq.coverage.genome_territory);
    eprintln!(
        "unified sequential = {seq_secs:.1}s, parallel = {par_secs:.1}s, speedup = {:.2}x",
        seq_secs / par_secs
    );
}

/// The parity gate of §4c. It puts the de-novo SNP calls of the Rust caller against a truth VCF
/// from GATK, on the chrM of HG002. Make that truth file with:
///
/// ```text
/// gatk HaplotypeCaller -I hg002.chrM.bam -R chm13v2.0.fa -L chrM \
///   --sample-ploidy 1 -O hg002.chrM.gatk.vcf.gz && bgzip -d hg002.chrM.gatk.vcf.gz
/// ```
#[test]
#[ignore = "requires GATK_CHRM_VCF + HG002_CHRM_BAM + CHM13_REF env vars"]
fn hg002_chrm_gatk_parity() {
    let (Ok(vcf), Ok(bam), Ok(reference)) = (
        std::env::var("GATK_CHRM_VCF"),
        std::env::var("HG002_CHRM_BAM"),
        std::env::var("CHM13_REF"),
    ) else {
        eprintln!("set GATK_CHRM_VCF, HG002_CHRM_BAM, CHM13_REF to run this test");
        return;
    };

    let truth = parse_truth_vcf(&PathBuf::from(vcf)).expect("parse GATK VCF");
    let calls = call_denovo(
        &PathBuf::from(&bam),
        &PathBuf::from(&reference),
        "chrM",
        &HaploidCallerParams::default(),
        &navigator_analysis::CancelToken::none(),
    )
    .expect("de-novo should succeed");

    let report = compare_denovo_snps(&truth, &calls);
    eprintln!(
        "SNP parity: matched={} rust_only(FP)={} truth_only(FN)={} truth_indels_skipped={}",
        report.matched_count(),
        report.rust_only.len(),
        report.truth_only.len(),
        report.truth_non_snp_alleles
    );
    eprintln!(
        "precision={:.3} recall={:.3} f1={:.3}",
        report.precision(),
        report.recall(),
        report.f1()
    );
    for fp in &report.rust_only {
        eprintln!("  FP {}:{} {}>{}", fp.chrom, fp.pos, fp.reference, fp.alternate);
    }
    for f_n in &report.truth_only {
        eprintln!("  FN {}:{} {}>{}", f_n.chrom, f_n.pos, f_n.reference, f_n.alternate);
    }

    // Gate: catch every GATK SNP; allow a few homopolymer-adjacent FPs (the §4b indel
    // risk) until local realignment lands.
    assert!(report.recall() >= 0.95, "recall {:.3} below gate", report.recall());
    assert!(
        report.precision() >= 0.80,
        "precision {:.3} below gate",
        report.precision()
    );

    // Chunked: a chunk boundary at 16300 splits the 16294 insertion (chunk 1) from the
    // smeared 16302 position (chunk 2). Both-side overlap must keep realignment correct,
    // so the call set is identical to the unchunked run.
    let chunked = HaploidCallerParams {
        denovo_chunk: 16_300,
        denovo_overlap: 500,
        ..HaploidCallerParams::default()
    };
    let chunked_calls = call_denovo(
        &PathBuf::from(&bam),
        &PathBuf::from(&reference),
        "chrM",
        &chunked,
        &navigator_analysis::CancelToken::none(),
    )
    .unwrap();
    let chunked_report = compare_denovo_snps(&truth, &chunked_calls);
    eprintln!(
        "chunked SNP parity: precision={:.3} recall={:.3} FP={}",
        chunked_report.precision(),
        chunked_report.recall(),
        chunked_report.rust_only.len()
    );
    assert_eq!(
        chunked_calls.iter().map(|c| c.position).collect::<Vec<_>>(),
        calls.iter().map(|c| c.position).collect::<Vec<_>>(),
        "chunk boundary changed the call set"
    );
}
