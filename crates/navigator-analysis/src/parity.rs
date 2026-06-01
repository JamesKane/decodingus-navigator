//! §4c golden-truth parity harness — compares the Rust haploid caller's de-novo SNP
//! calls against a GATK truth VCF and reports concordance (precision/recall).
//!
//! The v1 caller is SNP-only, so GATK indel/MNP alleles are counted separately and
//! excluded from the SNP concordance rather than scored as misses — an honest fair
//! comparison until local realignment lands (plan §4b). This is the cutover gate and
//! the regression guard for the analysis layer; the comparison logic is pure and
//! unit-tested, and an ignored test drives it against real GATK output.

use std::collections::BTreeSet;
use std::io::BufReader;
use std::path::Path;

pub use du_bio::vcf::VcfVariant;

use crate::caller::VariantCall;
use crate::error::AnalysisError;

/// Parse a (plain-text) truth VCF — e.g. GATK output decompressed with `bgzip -d`.
/// Reuses the shared `du-bio` variant-column parser.
pub fn parse_truth_vcf(path: &Path) -> Result<Vec<VcfVariant>, AnalysisError> {
    let file = std::fs::File::open(path).map_err(|e| AnalysisError::io(path, e))?;
    du_bio::vcf::parse(BufReader::new(file))
        .map_err(|e| AnalysisError::Message(format!("parsing {}: {e}", path.display())))
}

/// A normalized single-nucleotide variant for set comparison.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct SnpKey {
    pub chrom: String,
    pub pos: i64,
    pub reference: char,
    pub alternate: char,
}

fn snp_key(chrom: &str, pos: i64, reference: &str, alternate: &str) -> Option<SnpKey> {
    let r = single_base(reference)?;
    let a = single_base(alternate)?;
    Some(SnpKey { chrom: chrom.to_string(), pos, reference: r, alternate: a })
}

/// Accepts a one-character A/C/G/T allele (case-insensitive); rejects indels/symbolic.
fn single_base(allele: &str) -> Option<char> {
    if allele.len() != 1 {
        return None;
    }
    let c = allele.as_bytes()[0].to_ascii_uppercase();
    matches!(c, b'A' | b'C' | b'G' | b'T').then(|| c as char)
}

/// Concordance between Rust de-novo SNP calls and a GATK truth set.
#[derive(Debug, Clone)]
pub struct ParityReport {
    /// SNPs called by both (position + ref + alt agree).
    pub matched: Vec<SnpKey>,
    /// Called only by the Rust caller (candidate false positives).
    pub rust_only: Vec<SnpKey>,
    /// In the truth set but not called by Rust (candidate false negatives).
    pub truth_only: Vec<SnpKey>,
    /// Truth alt alleles excluded from the SNP comparison (indels/MNPs/symbolic).
    pub truth_non_snp_alleles: usize,
}

impl ParityReport {
    pub fn matched_count(&self) -> usize {
        self.matched.len()
    }

    /// matched / (matched + rust_only). 1.0 when there are no calls at all.
    pub fn precision(&self) -> f64 {
        let denom = self.matched.len() + self.rust_only.len();
        if denom == 0 { 1.0 } else { self.matched.len() as f64 / denom as f64 }
    }

    /// matched / (matched + truth_only). 1.0 when the truth set is empty.
    pub fn recall(&self) -> f64 {
        let denom = self.matched.len() + self.truth_only.len();
        if denom == 0 { 1.0 } else { self.matched.len() as f64 / denom as f64 }
    }

    pub fn f1(&self) -> f64 {
        let (p, r) = (self.precision(), self.recall());
        if p + r == 0.0 { 0.0 } else { 2.0 * p * r / (p + r) }
    }
}

/// Compare Rust de-novo SNP calls against a parsed GATK truth VCF.
pub fn compare_denovo_snps(truth: &[VcfVariant], rust_calls: &[VariantCall]) -> ParityReport {
    let mut truth_snps: BTreeSet<SnpKey> = BTreeSet::new();
    let mut truth_non_snp_alleles = 0usize;
    for v in truth {
        for alt in &v.alternate {
            match snp_key(&v.chrom, v.pos, &v.reference, alt) {
                Some(k) => {
                    truth_snps.insert(k);
                }
                None => truth_non_snp_alleles += 1,
            }
        }
    }

    let rust_snps: BTreeSet<SnpKey> = rust_calls
        .iter()
        .filter_map(|c| {
            snp_key(
                &c.contig,
                c.position,
                &c.reference_allele.to_string(),
                &c.alternate_allele.to_string(),
            )
        })
        .collect();

    ParityReport {
        matched: truth_snps.intersection(&rust_snps).cloned().collect(),
        rust_only: rust_snps.difference(&truth_snps).cloned().collect(),
        truth_only: truth_snps.difference(&rust_snps).cloned().collect(),
        truth_non_snp_alleles,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rust_call(pos: i64, r: char, a: char) -> VariantCall {
        VariantCall {
            contig: "chrM".into(),
            position: pos,
            reference_allele: r,
            alternate_allele: a,
            depth: 100,
            alt_depth: 100,
            allele_fraction: 1.0,
        }
    }

    #[test]
    fn compares_snp_sets_and_excludes_truth_indels() {
        let truth = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chrM\t100\t.\tC\tT\t.\t.\t.
chrM\t200\t.\tA\tG\t.\t.\t.
chrM\t300\t.\tG\tC\t.\t.\t.
chrM\t400\t.\tAT\tA\t.\t.\t.
";
        let truth = du_bio::vcf::parse(truth.as_bytes()).unwrap();
        let calls = vec![
            rust_call(100, 'C', 'T'), // match
            rust_call(200, 'A', 'G'), // match
            rust_call(500, 'T', 'A'), // rust-only (FP)
        ];
        // 300 G>C is truth-only (FN); 400 AT>A is an indel, excluded.
        let r = compare_denovo_snps(&truth, &calls);
        assert_eq!(r.matched_count(), 2);
        assert_eq!(r.rust_only.len(), 1);
        assert_eq!(r.rust_only[0].pos, 500);
        assert_eq!(r.truth_only.len(), 1);
        assert_eq!(r.truth_only[0].pos, 300);
        assert_eq!(r.truth_non_snp_alleles, 1);
        assert!((r.precision() - 2.0 / 3.0).abs() < 1e-9);
        assert!((r.recall() - 2.0 / 3.0).abs() < 1e-9);
    }

    #[test]
    fn empty_sets_are_fully_concordant() {
        let r = compare_denovo_snps(&[], &[]);
        assert_eq!(r.matched_count(), 0);
        assert!((r.precision() - 1.0).abs() < 1e-9);
        assert!((r.recall() - 1.0).abs() < 1e-9);
    }
}
