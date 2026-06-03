//! Diploid genotype-likelihood genotyping against the diploid fixture (chr1, depth 20):
//!   pos1 hom-ref(A), pos2 het(C/G), pos5 het(A/T), pos8 hom-alt(T->A).

use std::path::PathBuf;

use navigator_analysis::caller::{genotype_sites, Site, SiteGenotype};
use navigator_analysis::caller::HaploidCallerParams;

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

fn site(name: &str, pos: i64, r: &str, a: &str) -> Site {
    Site { name: name.into(), contig: "chr1".into(), position: pos, reference_allele: r.into(), alternate_allele: a.into() }
}

fn by_name(calls: &[SiteGenotype], name: &str) -> SiteGenotype {
    calls.iter().find(|c| c.name == name).unwrap().clone()
}

#[test]
fn diploid_panel_genotyping_yields_dosages() {
    let dir = fixtures();
    let sites = vec![
        site("homref", 1, "A", "G"),
        site("het_cg", 2, "C", "G"),
        site("het_at", 5, "A", "T"),
        site("homalt", 8, "T", "A"),
        site("missing", 9999, "A", "C"), // off-contig -> no coverage -> no-call
    ];
    let calls = genotype_sites(&dir.join("diploid.bam"), "chr1", &sites, 2, &HaploidCallerParams::default(), None)
        .expect("genotyping should succeed");

    let homref = by_name(&calls, "homref");
    assert_eq!(homref.dosage, 0);
    assert_eq!(homref.depth, 20);
    assert_eq!((homref.ref_depth, homref.alt_depth), (20, 0));
    assert!(homref.gq > 50);

    let het = by_name(&calls, "het_cg");
    assert_eq!(het.dosage, 1);
    assert_eq!((het.ref_depth, het.alt_depth), (10, 10));
    assert!(het.gq > 50);

    assert_eq!(by_name(&calls, "het_at").dosage, 1);

    let homalt = by_name(&calls, "homalt");
    assert_eq!(homalt.dosage, 2);
    assert_eq!((homalt.ref_depth, homalt.alt_depth), (0, 20));

    let missing = by_name(&calls, "missing");
    assert_eq!(missing.dosage, -1); // no-call
    assert_eq!(missing.depth, 0);
}

#[test]
fn haploid_genotyping_calls_zero_or_one() {
    let dir = fixtures();
    // ploidy 1: the hom-alt site reads as the alt allele (dosage 1); hom-ref as 0.
    let sites = vec![site("homref", 1, "A", "G"), site("homalt", 8, "T", "A")];
    let calls = genotype_sites(&dir.join("diploid.bam"), "chr1", &sites, 1, &HaploidCallerParams::default(), None).unwrap();
    assert_eq!(by_name(&calls, "homref").dosage, 0);
    assert_eq!(by_name(&calls, "homalt").dosage, 1);
    assert!(by_name(&calls, "homalt").pls.len() == 2); // ploidy 1 -> 2 genotypes
}
