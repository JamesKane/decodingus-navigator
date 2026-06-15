//! Diploid genotype-likelihood genotyping against the diploid fixture (chr1, depth 20):
//!   pos1 hom-ref(A), pos2 het(C/G), pos5 het(A/T), pos8 hom-alt(T->A).

use std::path::PathBuf;

use navigator_analysis::caller::{call_denovo_diploid, genotype_sites, Site, SiteGenotype};
use navigator_analysis::caller::HaploidCallerParams;
use navigator_analysis::vcf::write_diploid_vcf;

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/fixtures")
}

/// The diploid fixture is chr1 over reference `ACGTACGTAC` (10 bp). Write that reference + its `.fai`
/// to a temp dir (the bundled `ref.fa` is chrM-only) so the de-novo caller can read it.
fn chr1_reference() -> PathBuf {
    let dir = std::env::temp_dir().join(format!("dun-diploid-ref-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let fa = dir.join("chr1.fa");
    std::fs::write(&fa, b">chr1\nACGTACGTAC\n").unwrap();
    // .fai: name, length, offset-of-first-base, bases-per-line, bytes-per-line.
    std::fs::write(dir.join("chr1.fa.fai"), b"chr1\t10\t6\t10\t11\n").unwrap();
    fa
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
fn denovo_diploid_calls_het_and_hom_alt_then_writes_vcf() {
    let dir = fixtures();
    let reference = chr1_reference();
    let calls =
        call_denovo_diploid(&dir.join("diploid.bam"), &reference, "chr1", &HaploidCallerParams::default()).unwrap();

    // Only the variant sites are emitted (the 7 hom-ref positions are not), in position order.
    let by_pos = |p: i64| calls.iter().find(|c| c.position == p).cloned();
    assert_eq!(calls.len(), 3, "expected pos 2/5/8 only, got {:?}", calls.iter().map(|c| c.position).collect::<Vec<_>>());

    let p2 = by_pos(2).expect("het at pos 2");
    assert_eq!((p2.reference_allele.as_str(), p2.alternate_allele.as_str()), ("C", "G"));
    assert_eq!(p2.dosage, 1); // 10 C / 10 G → 0/1
    assert_eq!((p2.ref_depth, p2.alt_depth), (10, 10));

    assert_eq!(by_pos(5).expect("het at pos 5").dosage, 1); // A/T het

    let p8 = by_pos(8).expect("hom-alt at pos 8");
    assert_eq!((p8.reference_allele.as_str(), p8.alternate_allele.as_str()), ("T", "A"));
    assert_eq!(p8.dosage, 2); // 20 A → 1/1

    // VCF round-trips the genotypes.
    let vcf = write_diploid_vcf("FIX", &calls);
    assert!(vcf.contains("chr1\t2\t.\tC\tG\t"));
    assert!(vcf.contains("\t0/1:10,10:20:"));
    assert!(vcf.contains("\t1/1:0,20:20:"));
    let _ = std::fs::remove_dir_all(reference.parent().unwrap());
}

#[test]
fn denovo_diploid_calls_a_heterozygous_deletion() {
    // indel.bam (chrM): 10 reads 50M (ref) + 10 reads 5M2D43M (2bp deletion of ref pos 6-7 = C,G).
    // The bundled ref.fa chrM is ACGTACGTAC… so pos5=A, pos6=C, pos7=G → REF=ACG, ALT=A, het 0/1.
    let dir = fixtures();
    let calls = call_denovo_diploid(&dir.join("indel.bam"), &dir.join("ref.fa"), "chrM", &HaploidCallerParams::default())
        .unwrap();
    let del = calls
        .iter()
        .find(|c| c.position == 5 && c.reference_allele.len() > c.alternate_allele.len())
        .expect("a deletion call at pos 5");
    assert_eq!((del.reference_allele.as_str(), del.alternate_allele.as_str()), ("ACG", "A"));
    assert_eq!(del.dosage, 1); // 10 deletion-reads / 10 ref-reads → heterozygous
    assert_eq!(del.alt_depth, 10);
    // It renders as a standard indel VCF record.
    let vcf = write_diploid_vcf("FIX", &calls);
    assert!(vcf.contains("chrM\t5\t.\tACG\tA\t"));
    assert!(vcf.contains("\t0/1:"));
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
