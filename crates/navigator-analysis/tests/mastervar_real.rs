//! A check of the parse of a real CompleteGenomics masterVar file. It carries `#[ignore]`, and it
//! needs a local dump.
//!
//! Run it against a real file, compressed or plain:
//!
//! ```text
//! MASTERVAR_TSV=/path/to/var-GS00253-DNA_A01_200_37-ASM.tsv.bz2 \
//!   cargo test -p navigator-analysis --test mastervar_real -- --ignored --nocapture
//! ```
//!
//! It streams the whole genome. It then prints the sample id and the reference build, the tallies
//! of the loci and the SNPs, and the call count of each contig. Last, it prints a breakdown of the
//! genotype shapes.
//!
//! That is enough to confirm two things. The rows that hold two alleles come together correctly.
//! And the haploid contigs, chrY and chrM, come through as hemizygous.

use std::collections::BTreeMap;

#[test]
#[ignore]
fn parse_real_master_var() {
    let Ok(path) = std::env::var("MASTERVAR_TSV") else {
        eprintln!("set MASTERVAR_TSV to run");
        return;
    };
    let t0 = std::time::Instant::now();
    let out = navigator_analysis::mastervar::parse_file(std::path::Path::new(&path)).expect("parse");
    let elapsed = t0.elapsed();

    println!("file: {path}");
    println!("sample: {:?}  build: {}", out.sample_id, out.reference_build);
    println!(
        "loci_seen: {}  snp_loci: {}  calls: {}  ({:.1}s)",
        out.loci_seen,
        out.snp_loci,
        out.calls.len(),
        elapsed.as_secs_f64()
    );
    assert_eq!(out.snp_loci as usize, out.calls.len());

    let mut per_contig: BTreeMap<&str, usize> = BTreeMap::new();
    let mut gt: BTreeMap<&str, usize> = BTreeMap::new();
    for c in &out.calls {
        *per_contig.entry(c.contig.as_str()).or_default() += 1;
        *gt.entry(c.genotype.as_deref().unwrap_or("?")).or_default() += 1;
    }
    println!("per-contig: {per_contig:?}");
    println!("genotype shapes: {gt:?}");

    // chrY and chrM must both be hemizygous, at the genotype "1". Neither may be diploid.
    for c in out.calls.iter().filter(|c| c.contig == "chrY" || c.contig == "chrM") {
        assert_eq!(
            c.genotype.as_deref(),
            Some("1"),
            "{}:{} should be hemizygous",
            c.contig,
            c.position
        );
    }
    // Every call is a clean single-base biallelic SNP.
    for c in &out.calls {
        assert_eq!(c.reference.len(), 1);
        assert_eq!(c.alternate.len(), 1);
    }
}
