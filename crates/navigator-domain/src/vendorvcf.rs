//! Vendor-VCF classification — recognize FTDNA Big Y, Full Genomes Y Elite, YSEQ, etc. from a
//! `.vcf` so the import can tag it (vendor label + a meaningful `SourceType`) instead of treating
//! every VCF as a generic `IMPORTED` set.
//!
//! Signals (from real exports): the `##source` meta line — FTDNA Big Y stamps `##source=aengine`
//! (its Arpeggi caller) — plus the contig set (chrY-only ⇒ Y-targeted, chrM-only ⇒ mtDNA), the file
//! name, and the sibling `readme.txt` FTDNA ships ("…BigY raw data…"). Mirrors the Scala
//! `VcfCache.VcfVendor`.

use crate::variants::SourceType;

/// A recognized sequencing vendor behind a VCF (or `Generic` when nothing matches).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VendorVcf {
    FtdnaBigY,
    FtdnaMtFull,
    Yseq,
    FullGenomes,
    Nebula,
    Dante,
    Generic,
}

impl VendorVcf {
    /// Human label for the variant-set `source_label`.
    pub fn display(self) -> &'static str {
        match self {
            VendorVcf::FtdnaBigY => "FTDNA Big Y",
            VendorVcf::FtdnaMtFull => "FTDNA mtFull Sequence",
            VendorVcf::Yseq => "YSEQ",
            VendorVcf::FullGenomes => "Full Genomes",
            VendorVcf::Nebula => "Nebula Genomics",
            VendorVcf::Dante => "Dante Labs",
            VendorVcf::Generic => "Imported VCF",
        }
    }

    /// Concordance weighting for the calls: vendor-grade targeted Y/mt sequencing is `TargetedNgs`;
    /// consumer WGS vendors are short-read WGS; an unrecognized VCF stays generic `Imported`.
    pub fn source_type(self) -> SourceType {
        match self {
            VendorVcf::FtdnaBigY | VendorVcf::FtdnaMtFull | VendorVcf::Yseq | VendorVcf::FullGenomes => {
                SourceType::TargetedNgs
            }
            VendorVcf::Nebula | VendorVcf::Dante => SourceType::WgsShortRead,
            VendorVcf::Generic => SourceType::Imported,
        }
    }

    pub fn is_recognized(self) -> bool {
        self != VendorVcf::Generic
    }
}

/// Classify a VCF from its header `meta` (the `##` lines, lower-casing handled here), the set of
/// contig names it declares, its `filename`, and an optional sibling `readme` text.
pub fn classify(meta: &str, contigs: &[String], filename: &str, readme: Option<&str>) -> VendorVcf {
    let hay = format!("{} {} {}", meta, filename, readme.unwrap_or("")).to_lowercase();
    let only = |pred: fn(&str) -> bool| !contigs.is_empty() && contigs.iter().all(|c| pred(c));
    let mt_only = only(is_mt_contig);
    let y_only = only(is_y_contig);

    // FTDNA: the Arpeggi caller (`aengine`) or an explicit "big y" mention. mtFull if it's chrM-only.
    if hay.contains("aengine") || hay.contains("big y") || hay.contains("bigy") || hay.contains("mtfull") {
        return if mt_only {
            VendorVcf::FtdnaMtFull
        } else {
            VendorVcf::FtdnaBigY
        };
    }
    if hay.contains("yseq") {
        return VendorVcf::Yseq;
    }
    if hay.contains("full genomes") || hay.contains("fullgenomes") || hay.contains("y elite") || hay.contains("yelite")
    {
        return VendorVcf::FullGenomes;
    }
    if hay.contains("nebula") {
        return VendorVcf::Nebula;
    }
    if hay.contains("dante") {
        return VendorVcf::Dante;
    }
    // No vendor token: a chrM-only or chrY-only VCF is still a recognizable mtFull / targeted-Y shape.
    if mt_only {
        return VendorVcf::FtdnaMtFull;
    }
    let _ = y_only;
    VendorVcf::Generic
}

fn core(name: &str) -> &str {
    name.strip_prefix("chr")
        .unwrap_or(name)
        .split('_')
        .next()
        .unwrap_or(name)
}
fn is_y_contig(name: &str) -> bool {
    core(name) == "Y"
}
fn is_mt_contig(name: &str) -> bool {
    matches!(core(name), "M" | "MT")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ftdna_big_y_from_source_meta() {
        // Real Big Y export: ##source=aengine, chrY-only contigs.
        let meta = "##fileformat=VCFv4.1\n##reference=ucsc.hg38.fasta\n##source=aengine\n##ae_version=v2.1.3";
        let contigs = vec!["chrY".to_string(), "chrY_KI270740v1_random".to_string()];
        let v = classify(meta, &contigs, "variants.vcf", Some("…the BigY raw data files…"));
        assert_eq!(v, VendorVcf::FtdnaBigY);
        assert_eq!(v.source_type(), SourceType::TargetedNgs);
        assert!(v.is_recognized());
    }

    #[test]
    fn ftdna_mtfull_when_chrm_only() {
        let contigs = vec!["chrM".to_string()];
        assert_eq!(
            classify("##source=aengine", &contigs, "variants.vcf", None),
            VendorVcf::FtdnaMtFull
        );
        // chrM-only with no vendor token is still recognizably mtFull-shaped.
        assert_eq!(
            classify("##source=other", &contigs, "mt.vcf", None),
            VendorVcf::FtdnaMtFull
        );
    }

    #[test]
    fn other_vendors_by_token() {
        let c: Vec<String> = vec![];
        assert_eq!(classify("", &c, "yseq_results.vcf", None), VendorVcf::Yseq);
        assert_eq!(
            classify("", &c, "x.vcf", Some("Full Genomes Y Elite")),
            VendorVcf::FullGenomes
        );
        assert_eq!(classify("##source=nebula", &c, "x.vcf", None), VendorVcf::Nebula);
        assert_eq!(classify("", &c, "dante_labs.vcf", None), VendorVcf::Dante);
    }

    #[test]
    fn generic_when_nothing_matches() {
        let contigs = vec!["chr1".to_string(), "chr2".to_string(), "chrX".to_string()];
        let v = classify("##source=GATK", &contigs, "sample.vcf", None);
        assert_eq!(v, VendorVcf::Generic);
        assert_eq!(v.source_type(), SourceType::Imported);
        assert!(!v.is_recognized());
    }
}
