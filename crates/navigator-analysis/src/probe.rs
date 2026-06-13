//! Probe a BAM/CRAM **header** for the metadata a user would otherwise type in: the reference
//! build (from `@SQ`), the aligner (from `@PG`), and the sequencing platform / instrument /
//! test type (from `@RG`). Only the SAM header is read — no records, and (for CRAM) no
//! reference FASTA — so it is cheap and runs before the reference is resolved.

use std::path::Path;

use noodles::sam::header::record::value::map::{program, read_group, reference_sequence};
use noodles::{bam, cram, sam};

use crate::error::AnalysisError;
use crate::reader::{detect_format, Format};

/// What we could infer from an alignment file's header. Every field is best-effort.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct AlignmentProbe {
    /// Canonical reference build, e.g. `"chm13v2.0"`, `"GRCh38"`, `"GRCh37"`.
    pub reference_build: Option<String>,
    /// Alignment program, e.g. `"pbmm2"`, `"minimap2"`, `"bwa-mem2"`.
    pub aligner: Option<String>,
    /// Sequencing platform (`@RG PL`), upper-cased: `ILLUMINA`, `PACBIO`, `ONT`, …
    pub platform: Option<String>,
    /// Instrument model (`@RG PM`), e.g. `"NovaSeq 6000"`, `"Sequel II"`.
    pub instrument_model: Option<String>,
    /// Best-guess test type code from the catalog, e.g. `"WGS_HIFI"`, `"WGS"`.
    pub test_type: Option<String>,
    /// Vendor hint scraped from the header (`@RG CN` center / `@PG` / `@CO`), e.g. `"FamilyTreeDNA"`,
    /// `"Full Genomes"`, `"YSEQ"` — refines a coverage-detected targeted-Y test into the specific
    /// vendor product. `None` when no recognizable vendor token appears.
    pub vendor_hint: Option<String>,
}

/// Probe `path`'s header. Reads only the SAM header (no records / no reference needed).
pub fn probe_alignment(path: &Path) -> Result<AlignmentProbe, AnalysisError> {
    let header = read_header_only(path)?;
    let reference_build = detect_build(&header);
    let aligner = detect_aligner(&header);
    let (platform, instrument_model) = detect_platform(&header);
    let test_type = detect_test_type(platform.as_deref(), aligner.as_deref(), &header);
    let vendor_hint = detect_vendor_hint(&header);
    Ok(AlignmentProbe { reference_build, aligner, platform, instrument_model, test_type, vendor_hint })
}

/// Recognizable sequencing-vendor tokens, mapped to a canonical display the test-type catalog
/// keys on (`testtype::targeted_y_for_vendor`).
const VENDOR_TOKENS: &[(&str, &str)] = &[
    ("familytreedna", "FamilyTreeDNA"),
    ("family tree dna", "FamilyTreeDNA"),
    ("ftdna", "FamilyTreeDNA"),
    ("big y", "FamilyTreeDNA"),
    ("full genomes", "Full Genomes"),
    ("fullgenomes", "Full Genomes"),
    ("y elite", "Full Genomes"),
    ("yseq", "YSEQ"),
];

/// Scrape a vendor hint from the `@RG CN` (sequencing center), `@PG` (program id/name/command),
/// and `@CO` (free-text comment) header lines. First recognizable token wins.
fn detect_vendor_hint(header: &sam::Header) -> Option<String> {
    let mut hay = String::new();
    for map in header.read_groups().values() {
        if let Some(cn) = map.other_fields().get(&read_group::tag::SEQUENCING_CENTER) {
            hay.push_str(&s(cn));
            hay.push(' ');
        }
    }
    for (id, map) in header.programs().as_ref().iter() {
        hay.push_str(&s(id));
        hay.push(' ');
        if let Some(pn) = map.other_fields().get(&program::tag::NAME) {
            hay.push_str(&s(pn));
            hay.push(' ');
        }
        if let Some(cl) = map.other_fields().get(&program::tag::COMMAND_LINE) {
            hay.push_str(&s(cl));
            hay.push(' ');
        }
    }
    for comment in header.comments() {
        hay.push_str(&s(comment));
        hay.push(' ');
    }
    let hay = hay.to_lowercase();
    VENDOR_TOKENS.iter().find(|(tok, _)| hay.contains(tok)).map(|(_, canon)| (*canon).to_string())
}

/// Read just the SAM header from a BAM or CRAM (CRAM's header doesn't need the reference).
fn read_header_only(path: &Path) -> Result<sam::Header, AnalysisError> {
    match detect_format(path) {
        Format::Bam => {
            let mut inner = bam::io::reader::Builder.build_from_path(path).map_err(|e| AnalysisError::io(path, e))?;
            inner.read_header().map_err(|e| AnalysisError::io(path, e))
        }
        Format::Cram => {
            let mut inner = cram::io::reader::Builder::default().build_from_path(path).map_err(|e| AnalysisError::io(path, e))?;
            inner.read_header().map_err(|e| AnalysisError::io(path, e))
        }
    }
}

/// Lossy UTF-8 of a header field value (a `bstr::BString` or sequence name — any `AsRef<[u8]>`).
fn s<T: AsRef<[u8]>>(v: &T) -> String {
    String::from_utf8_lossy(v.as_ref()).into_owned()
}

/// Reference build from `@SQ`: prefer the assembly (`AS`) tag, else the chr1 length signature.
fn detect_build(header: &sam::Header) -> Option<String> {
    // The Y-PAR-masked + rCRS CHM13 analysis set is indistinguishable from plain chm13v2.0
    // by `@SQ` (same contig names/lengths), so check the reference filename the aligner
    // recorded (`@PG CL` / `@SQ UR`) first. Plain chm13 won't match this signature.
    if header_mentions_masked_rcrs(header) {
        return Some("chm13v2.0_maskedY_rCRS".into());
    }
    let seqs = header.reference_sequences();
    for map in seqs.values() {
        if let Some(asm) = map.other_fields().get(&reference_sequence::tag::ASSEMBLY_ID) {
            if let Some(b) = normalize_build(&s(asm)) {
                return Some(b);
            }
        }
    }
    for (name, map) in seqs.iter() {
        let n = s(name);
        if n == "chr1" || n == "1" {
            return match map.length().get() {
                248_956_422 => Some("GRCh38".into()),
                249_250_621 => Some("GRCh37".into()),
                248_387_328 => Some("chm13v2.0".into()),
                _ => None,
            };
        }
    }
    None
}

/// Does the header name the masked+rCRS analysis-set FASTA (`chm13v2.0_maskedY_rCRS`)? The
/// aligner records the reference path in `@PG CL`, and some pipelines stamp it as the `@SQ`
/// URI (`UR`). Case-insensitive match on the distinctive `maskedY_rCRS` token.
fn header_mentions_masked_rcrs(header: &sam::Header) -> bool {
    const NEEDLE: &str = "maskedy_rcrs";
    for (_id, map) in header.programs().as_ref().iter() {
        if let Some(cl) = map.other_fields().get(&program::tag::COMMAND_LINE) {
            if s(cl).to_lowercase().contains(NEEDLE) {
                return true;
            }
        }
    }
    for map in header.reference_sequences().values() {
        if let Some(uri) = map.other_fields().get(&reference_sequence::tag::URI) {
            if s(uri).to_lowercase().contains(NEEDLE) {
                return true;
            }
        }
    }
    false
}

/// Map an assembly label (`AS` tag) onto a canonical build, or `None` if unrecognized.
fn normalize_build(asm: &str) -> Option<String> {
    let l = asm.to_lowercase();
    if l.contains("chm13") || l.contains("t2t") {
        Some("chm13v2.0".into())
    } else if l.contains("grch38") || l.contains("hg38") {
        Some("GRCh38".into())
    } else if l.contains("grch37") || l.contains("hg19") {
        Some("GRCh37".into())
    } else {
        None
    }
}

/// Known aligners, longest/most-specific first (so `bwa-mem2` wins over `bwa`).
const ALIGNERS: &[(&str, &str)] = &[
    ("pbmm2", "pbmm2"),
    ("bwa-mem2", "bwa-mem2"),
    ("minimap2", "minimap2"),
    ("winnowmap", "winnowmap"),
    ("bowtie2", "bowtie2"),
    ("novoalign", "novoalign"),
    ("ngmlr", "ngmlr"),
    ("dragen", "dragen"),
    ("hisat2", "hisat2"),
    ("bwa", "bwa"),
];

/// Aligner from `@PG`: match the program id / name (`PN`) / command line (`CL`) against the
/// known-aligner list, so non-aligner programs (samtools, gatk, …) are ignored.
fn detect_aligner(header: &sam::Header) -> Option<String> {
    for (id, map) in header.programs().as_ref().iter() {
        let pn = map.other_fields().get(&program::tag::NAME).map(s).unwrap_or_default();
        let cl = map.other_fields().get(&program::tag::COMMAND_LINE).map(s).unwrap_or_default();
        let hay = format!("{} {} {}", s(id), pn, cl).to_lowercase();
        for (needle, canon) in ALIGNERS {
            if hay.contains(needle) {
                return Some((*canon).to_string());
            }
        }
    }
    None
}

/// Platform (`PL`, upper-cased) + instrument model (`PM`) from the first informative `@RG`.
fn detect_platform(header: &sam::Header) -> (Option<String>, Option<String>) {
    for map in header.read_groups().values() {
        let pl = map.other_fields().get(&read_group::tag::PLATFORM).map(|v| s(v).to_uppercase());
        let pm = map.other_fields().get(&read_group::tag::PLATFORM_MODEL).map(s);
        if pl.is_some() || pm.is_some() {
            return (pl, pm);
        }
    }
    (None, None)
}

/// Best-guess test-type code: by platform (PacBio→HiFi, Nanopore→ONT WGS, Illumina→WGS); a
/// recognized aligner with no platform still implies a whole-genome run.
fn detect_test_type(platform: Option<&str>, aligner: Option<&str>, _header: &sam::Header) -> Option<String> {
    match platform {
        Some(p) if p.contains("PACBIO") => Some("WGS_HIFI".into()),
        Some(p) if p.contains("NANOPORE") || p == "ONT" => Some("WGS_NANOPORE".into()),
        Some(p) if p.contains("ILLUMINA") => Some("WGS".into()),
        _ => aligner.map(|_| "WGS".into()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalize_build_variants() {
        assert_eq!(normalize_build("GRCh38"), Some("GRCh38".into()));
        assert_eq!(normalize_build("hg38"), Some("GRCh38".into()));
        assert_eq!(normalize_build("T2T-CHM13v2.0"), Some("chm13v2.0".into()));
        assert_eq!(normalize_build("hg19"), Some("GRCh37".into()));
        assert_eq!(normalize_build("mystery"), None);
    }

    fn header_from_sam(text: &str) -> sam::Header {
        sam::io::Reader::new(text.as_bytes()).read_header().unwrap()
    }

    #[test]
    fn detects_masked_rcrs_from_the_reference_filename() {
        // Same @SQ as plain chm13 (chr1 length), but the aligner's @PG names the masked FASTA.
        let masked = header_from_sam(
            "@HD\tVN:1.6\n\
             @SQ\tSN:chr1\tLN:248387328\n\
             @PG\tID:bwa-mem2\tPN:bwa-mem2\tCL:bwa-mem2 mem /refs/chm13v2.0_maskedY_rCRS.fa r.fq\n",
        );
        assert_eq!(detect_build(&masked), Some("chm13v2.0_maskedY_rCRS".into()));

        // Without that signature the same @SQ falls back to plain chm13 (chr1 length).
        let plain = header_from_sam(
            "@HD\tVN:1.6\n\
             @SQ\tSN:chr1\tLN:248387328\n\
             @PG\tID:bwa-mem2\tPN:bwa-mem2\tCL:bwa-mem2 mem /refs/chm13v2.0.fa r.fq\n",
        );
        assert_eq!(detect_build(&plain), Some("chm13v2.0".into()));
    }

    #[test]
    fn vendor_hint_from_center_program_or_comment() {
        // @RG CN (sequencing center).
        let rg = header_from_sam("@HD\tVN:1.6\n@SQ\tSN:chrY\tLN:57227415\n@RG\tID:r1\tCN:FamilyTreeDNA\tPL:ILLUMINA\n");
        assert_eq!(detect_vendor_hint(&rg).as_deref(), Some("FamilyTreeDNA"));
        // @CO free-text comment naming the product.
        let co = header_from_sam("@HD\tVN:1.6\n@SQ\tSN:chrY\tLN:57227415\n@CO\tFull Genomes Y Elite v2\n");
        assert_eq!(detect_vendor_hint(&co).as_deref(), Some("Full Genomes"));
        // No vendor token.
        let plain = header_from_sam("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:248387328\n@RG\tID:r1\tPL:PACBIO\n");
        assert_eq!(detect_vendor_hint(&plain), None);
    }

    #[test]
    fn test_type_by_platform() {
        assert_eq!(detect_test_type(Some("PACBIO"), Some("pbmm2"), &sam::Header::default()), Some("WGS_HIFI".into()));
        assert_eq!(detect_test_type(Some("ILLUMINA"), Some("bwa-mem2"), &sam::Header::default()), Some("WGS".into()));
        assert_eq!(detect_test_type(None, Some("minimap2"), &sam::Header::default()), Some("WGS".into()));
        assert_eq!(detect_test_type(None, None, &sam::Header::default()), None);
    }
}
