//! Minimal VCF 4.2 writer for diploid genotype calls (the de-novo diploid caller + known-site
//! genotyping output). Emits one variant record per [`SiteGenotype`] with `FORMAT GT:AD:DP:GQ:PL`.
//! Records are written in the order given (the caller returns ascending position per contig).

use crate::caller::SiteGenotype;

/// `dosage` → diploid `GT`. `-1` (no-call) → `./.`, else `0/0` | `0/1` | `1/1`.
fn genotype_field(dosage: i32) -> &'static str {
    match dosage {
        0 => "0/0",
        1 => "0/1",
        2 => "1/1",
        _ => "./.",
    }
}

/// Write a diploid VCF (single sample) from `calls`. `sample` names the sample column. The QUAL of
/// each record is the hom-ref PL (how unlikely the no-variant genotype is), i.e. `pls[0]`.
pub fn write_diploid_vcf(sample: &str, calls: &[SiteGenotype]) -> String {
    let mut out = String::new();
    out.push_str("##fileformat=VCFv4.2\n");
    out.push_str("##source=navigator-diploid-caller\n");
    out.push_str("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    out.push_str("##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths (ref,alt)\">\n");
    out.push_str("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth\">\n");
    out.push_str("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n");
    out.push_str("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods\">\n");
    out.push_str(&format!(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
    ));

    for c in calls {
        let qual = c.pls.first().map(|p| p.to_string()).unwrap_or_else(|| ".".to_string());
        let pl = if c.pls.is_empty() {
            ".".to_string()
        } else {
            c.pls.iter().map(|p| p.to_string()).collect::<Vec<_>>().join(",")
        };
        // Multiallelic sites carry an explicit GT string + per-allele AD; biallelic sites derive
        // both from `dosage` / `ref_depth,alt_depth`.
        let gt = c.gt.clone().unwrap_or_else(|| genotype_field(c.dosage).to_string());
        let ad = match &c.allele_depths {
            Some(d) => d.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(","),
            None => format!("{},{}", c.ref_depth, c.alt_depth),
        };
        out.push_str(&format!(
            "{}\t{}\t.\t{}\t{}\t{}\t.\t.\tGT:AD:DP:GQ:PL\t{}:{}:{}:{}:{}\n",
            c.contig, c.position, c.reference_allele, c.alternate_allele, qual, gt, ad, c.depth, c.gq, pl,
        ));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn call(pos: i64, r: &str, a: &str, dosage: i32, ref_d: u32, alt_d: u32, pls: Vec<u8>) -> SiteGenotype {
        SiteGenotype {
            name: String::new(),
            contig: "chr1".into(),
            position: pos,
            reference_allele: r.into(),
            alternate_allele: a.into(),
            ploidy: 2,
            dosage,
            gq: 99,
            depth: ref_d + alt_d,
            ref_depth: ref_d,
            alt_depth: alt_d,
            pls,
            gt: None,
            allele_depths: None,
        }
    }

    #[test]
    fn writes_header_and_diploid_records() {
        let calls = vec![
            call(2, "C", "G", 1, 10, 10, vec![120, 0, 130]), // het 0/1
            call(8, "T", "A", 2, 0, 20, vec![200, 90, 0]),   // hom-alt 1/1
        ];
        let vcf = write_diploid_vcf("KANE-0001", &calls);
        assert!(vcf.starts_with("##fileformat=VCFv4.2"));
        assert!(vcf.contains("\tFORMAT\tKANE-0001\n"));
        // het: QUAL = pls[0] = 120; GT 0/1; AD 10,10; PL joined.
        assert!(vcf.contains("chr1\t2\t.\tC\tG\t120\t.\t.\tGT:AD:DP:GQ:PL\t0/1:10,10:20:99:120,0,130\n"));
        // hom-alt: GT 1/1.
        assert!(vcf.contains("chr1\t8\t.\tT\tA\t200\t.\t.\tGT:AD:DP:GQ:PL\t1/1:0,20:20:99:200,90,0\n"));
    }

    #[test]
    fn writes_a_multiallelic_record() {
        // A 1/2 site (two indel alleles): explicit GT + 3-value AD + PL over the 6 diploid genotypes.
        let mut c = call(5, "ACG", "A,AT", -1, 0, 0, vec![200, 90, 0, 95, 5, 0]);
        c.gt = Some("1/2".into());
        c.allele_depths = Some(vec![0, 9, 8]);
        c.depth = 17;
        let vcf = write_diploid_vcf("S", &[c]);
        assert!(
            vcf.contains("chr1\t5\t.\tACG\tA,AT\t200\t.\t.\tGT:AD:DP:GQ:PL\t1/2:0,9,8:17:99:200,90,0,95,5,0\n"),
            "{vcf}"
        );
    }

    #[test]
    fn no_call_genotype_field() {
        assert_eq!(genotype_field(-1), "./.");
        assert_eq!(genotype_field(0), "0/0");
        assert_eq!(genotype_field(1), "0/1");
        assert_eq!(genotype_field(2), "1/1");
    }
}
