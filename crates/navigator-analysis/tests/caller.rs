//! Tests of the haploid caller, against the coverage fixture that the tests share, in
//! tests/fixtures.
//!
//! Every read of that fixture is `A`, over the reference chrM, which is `ACGTACGT...` with an N at
//! 25:
//!
//! ```text
//! pos 1-10   depth 4  MAPQ 60  -> passes the filters, consensus A
//! pos 11-20  depth 2  MAPQ 60  -> below the min_depth of 4
//! pos 26-30  depth 5  MAPQ 0   -> the min_mapping_quality of 20 drops it
//! ```
//!
//! With the default parameters, a de-novo call lands only at 1 to 10, where the reference base is
//! not `A`. Those bases are A C G T A C G T A C, so the SNPs are at {2,3,4,6,7,8,10}.

use std::collections::HashSet;
use std::path::PathBuf;

use navigator_analysis::caller::{
    call_denovo, force_call_sites, subtract_known, CalledAllele, HaploidCallerParams, Site,
};

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

fn site(name: &str, pos: i64, r: &str, a: &str) -> Site {
    Site {
        name: name.into(),
        contig: "chrM".into(),
        position: pos,
        reference_allele: r.into(),
        alternate_allele: a.into(),
    }
}

#[test]
fn denovo_calls_snps_only_where_consensus_differs_from_ref() {
    let dir = fixtures();
    let calls = call_denovo(
        &dir.join("coverage.bam"),
        &dir.join("ref.fa"),
        "chrM",
        &HaploidCallerParams::default(),
        &navigator_analysis::CancelToken::none(),
    )
    .expect("de-novo should succeed");

    let positions: Vec<i64> = calls.iter().map(|c| c.position).collect();
    assert_eq!(positions, vec![2, 3, 4, 6, 7, 8, 10]);

    // Every call is ref -> A at depth 4, fraction 1.0.
    for c in &calls {
        assert_eq!(c.alternate_allele, 'A');
        assert_eq!(c.depth, 4);
        assert_eq!(c.alt_depth, 4);
        assert!((c.allele_fraction - 1.0).abs() < 1e-9);
    }
    // Reference bases at the called positions.
    let by_pos: std::collections::HashMap<i64, char> = calls.iter().map(|c| (c.position, c.reference_allele)).collect();
    assert_eq!(by_pos[&2], 'C');
    assert_eq!(by_pos[&3], 'G');
    assert_eq!(by_pos[&4], 'T');
}

#[test]
fn force_call_genotypes_known_sites() {
    let dir = fixtures();
    let sites = vec![
        site("ref_match", 1, "A", "G"),    // consensus A == ref      -> Reference
        site("alt_match", 2, "C", "A"),    // consensus A == alt      -> Alternate
        site("third", 3, "G", "T"),        // consensus A is neither  -> NoCall
        site("shallow", 11, "G", "A"),     // depth 2 < min 4         -> NoCall
        site("filtered", 26, "A", "C"),    // MAPQ-0 reads dropped    -> NoCall
        site("indel", 5, "AT", "A"),       // not a SNP               -> skipped
        site("offcontig", 9999, "A", "C"), // beyond contig         -> skipped
    ];

    let calls = force_call_sites(
        &dir.join("coverage.bam"),
        "chrM",
        &sites,
        &HaploidCallerParams::default(),
        None,
    )
    .expect("force-call should succeed");

    // The code drops the indel and the site off the contig. 5 SNP sites stay.
    assert_eq!(calls.len(), 5);
    let by_name: std::collections::HashMap<&str, &_> = calls.iter().map(|c| (c.name.as_str(), c)).collect();

    assert_eq!(by_name["ref_match"].called, CalledAllele::Reference);
    assert_eq!(by_name["ref_match"].ref_depth, 4);
    assert_eq!(by_name["ref_match"].alt_depth, 0);

    assert_eq!(by_name["alt_match"].called, CalledAllele::Alternate);
    assert_eq!(by_name["alt_match"].alt_depth, 4);
    assert!((by_name["alt_match"].allele_fraction - 1.0).abs() < 1e-9);

    assert_eq!(by_name["third"].called, CalledAllele::NoCall);
    assert_eq!(by_name["shallow"].called, CalledAllele::NoCall);
    assert_eq!(by_name["shallow"].depth, 2);
    assert_eq!(by_name["filtered"].called, CalledAllele::NoCall);
    assert_eq!(by_name["filtered"].depth, 0);
}

#[test]
fn denovo_chunking_matches_unchunked() {
    let dir = fixtures();
    let base = call_denovo(
        &dir.join("coverage.bam"),
        &dir.join("ref.fa"),
        "chrM",
        &HaploidCallerParams::default(),
        &navigator_analysis::CancelToken::none(),
    )
    .unwrap();
    // Force many tiny chunks over the 50 bp fixture; result must be identical.
    let chunked = HaploidCallerParams {
        denovo_chunk: 8,
        denovo_overlap: 3,
        ..HaploidCallerParams::default()
    };
    let got = call_denovo(
        &dir.join("coverage.bam"),
        &dir.join("ref.fa"),
        "chrM",
        &chunked,
        &navigator_analysis::CancelToken::none(),
    )
    .unwrap();
    assert_eq!(got, base);
    assert_eq!(
        got.iter().map(|c| c.position).collect::<Vec<_>>(),
        vec![2, 3, 4, 6, 7, 8, 10]
    );
}

#[test]
fn private_set_subtracts_known_tree_positions() {
    let dir = fixtures();
    let calls = call_denovo(
        &dir.join("coverage.bam"),
        &dir.join("ref.fa"),
        "chrM",
        &HaploidCallerParams::default(),
        &navigator_analysis::CancelToken::none(),
    )
    .unwrap();

    // This test treats positions 2 and 3 as known tree sites.
    let known: HashSet<i64> = [2, 3].into_iter().collect();
    let private = subtract_known(&calls, &known);
    assert_eq!(
        private.iter().map(|c| c.position).collect::<Vec<_>>(),
        vec![4, 6, 7, 8, 10]
    );
}
