//! Short-tandem-repeat reference loci, parsed from a HipSTR-format reference BED.
//!
//! The HipSTR reference gives **tight repeat tracts**, and not loose feature regions. That is the
//! coordinate precision that a repeat counter over enclosing reads needs. Each line holds tabs
//! between its fields:
//!
//! ```text
//! chrom  start(0-based)  end  period  ref_copies  locus_id  motif
//! Y      10001           10038  6      6.33333     Human_STR_1604566  AACCCT
//! ```
//!
//! Sometimes `motif` holds a set of alternatives, with a `/` between them, as in `CCTT/CCCT`. The
//! code takes the first one as canonical. A contig name is bare, as `1` or `Y`, and a caller
//! normalizes it against the names in the BAM.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use flate2::read::MultiGzDecoder;

use crate::error::AnalysisError;

/// One STR locus: a tight repeat tract with its period (motif length) and reference copy number.
#[derive(Debug, Clone, PartialEq)]
pub struct StrLocus {
    /// The contig, as the reference BED writes it. It is bare: `1`, `X` or `Y`.
    pub contig: String,
    /// 0-based, half-open tract start (BED convention).
    pub start: i64,
    /// Exclusive tract end. Tract length in bp = `end - start`.
    pub end: i64,
    /// Repeat-unit length (motif size) in bp.
    pub period: u8,
    /// The copy number of the reference allele. It can hold a fraction, when the last unit is
    /// partial.
    pub ref_copies: f64,
    /// Locus id (HipSTR `Human_STR_N`), used as the result name until a vendor mapping exists.
    pub name: String,
    /// Canonical repeat motif (the first when the BED lists `A/B` alternatives).
    pub motif: String,
}

impl StrLocus {
    /// True when `name`, which is the contig of the BED, matches `query`. The comparison removes a
    /// `chr` prefix from either side first. The BAM can say `chrY` where the BED says `Y`. See the
    /// convention for the names of the contigs.
    pub fn contig_matches(&self, query: &str) -> bool {
        crate::contig::bare(&self.contig).eq_ignore_ascii_case(crate::contig::bare(query))
    }
}

/// Parse one BED line into a locus. Returns `None` for blank/comment/short lines.
fn parse_line(line: &str) -> Option<StrLocus> {
    if line.is_empty() || line.starts_with('#') {
        return None;
    }
    let mut f = line.split('\t');
    let contig = f.next()?.to_string();
    let start: i64 = f.next()?.parse().ok()?;
    let end: i64 = f.next()?.parse().ok()?;
    let period: u8 = f.next()?.parse().ok()?;
    let ref_copies: f64 = f.next()?.parse().ok()?;
    let name = f.next().unwrap_or("").to_string();
    // Motif is optional in the spec; HipSTR's reference includes it. Take the first alternative.
    let motif = f.next().unwrap_or("").split('/').next().unwrap_or("").to_string();
    Some(StrLocus {
        contig,
        start,
        end,
        period,
        ref_copies,
        name,
        motif,
    })
}

/// Read the STR loci from a HipSTR reference BED. That BED may come through gzip. This keeps a
/// locus only when that locus is on `contig`, and when its `period` is `min_period` or more. The
/// match on the contig ignores a `chr` prefix.
///
/// The filter runs as the code streams the file. It thereby never holds the genome-wide set of
/// about 1.6M loci in memory, when the caller needs one chromosome. The results come back in order
/// of their start.
///
/// A `min_period` of 2 drops the homopolymers, which have a period of 1. Those are noisy, and they
/// are not genealogical markers.
pub fn load_hipstr_contig(bed_gz: &Path, contig: &str, min_period: u8) -> Result<Vec<StrLocus>, AnalysisError> {
    let file = File::open(bed_gz).map_err(|e| AnalysisError::io(bed_gz, e))?;
    let reader = BufReader::new(MultiGzDecoder::new(file));
    let mut out = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(|e| AnalysisError::io(bed_gz, e))?;
        if let Some(locus) = parse_line(&line) {
            if locus.period >= min_period && locus.contig_matches(contig) {
                out.push(locus);
            }
        }
    }
    out.sort_by_key(|l| l.start);
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_a_hipstr_line() {
        let l = parse_line("Y\t10001\t10038\t6\t6.33333\tHuman_STR_1604566\tAACCCT").unwrap();
        assert_eq!(l.contig, "Y");
        assert_eq!((l.start, l.end), (10001, 10038));
        assert_eq!(l.period, 6);
        assert!((l.ref_copies - 6.33333).abs() < 1e-4);
        assert_eq!(l.name, "Human_STR_1604566");
        assert_eq!(l.motif, "AACCCT");
    }

    #[test]
    fn takes_first_motif_alternative_and_matches_contig_prefix_insensitively() {
        let l = parse_line("Y\t12946\t13016\t4\t17.75\tHuman_STR_1604569\tCCTT/CCCT").unwrap();
        assert_eq!(l.motif, "CCTT");
        assert!(l.contig_matches("chrY"));
        assert!(l.contig_matches("Y"));
        assert!(!l.contig_matches("chr1"));
    }

    #[test]
    fn skips_comments_and_short_lines() {
        assert!(parse_line("# header").is_none());
        assert!(parse_line("").is_none());
        assert!(parse_line("Y\t1").is_none());
    }
}
