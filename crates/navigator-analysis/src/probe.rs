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
}

/// Probe `path`'s header. Reads only the SAM header (no records / no reference needed).
pub fn probe_alignment(path: &Path) -> Result<AlignmentProbe, AnalysisError> {
    let header = read_header_only(path)?;
    let reference_build = detect_build(&header);
    let aligner = detect_aligner(&header);
    let (platform, instrument_model) = detect_platform(&header);
    let test_type = detect_test_type(platform.as_deref(), aligner.as_deref(), &header);
    Ok(AlignmentProbe { reference_build, aligner, platform, instrument_model, test_type })
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

    #[test]
    fn test_type_by_platform() {
        assert_eq!(detect_test_type(Some("PACBIO"), Some("pbmm2"), &sam::Header::default()), Some("WGS_HIFI".into()));
        assert_eq!(detect_test_type(Some("ILLUMINA"), Some("bwa-mem2"), &sam::Header::default()), Some("WGS".into()));
        assert_eq!(detect_test_type(None, Some("minimap2"), &sam::Header::default()), Some("WGS".into()));
        assert_eq!(detect_test_type(None, None, &sam::Header::default()), None);
    }
}
