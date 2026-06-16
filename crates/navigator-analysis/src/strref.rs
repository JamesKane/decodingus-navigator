//! Short-tandem-repeat reference loci, parsed from a HipSTR-format reference BED.
//!
//! The HipSTR reference defines **tight repeat tracts** (not loose feature regions) — the
//! coordinate precision an enclosing-read repeat counter needs. Each line is tab-delimited:
//!
//! ```text
//! chrom  start(0-based)  end  period  ref_copies  locus_id  motif
//! Y      10001           10038  6      6.33333     Human_STR_1604566  AACCCT
//! ```
//!
//! `motif` is occasionally a `/`-separated alternative set (e.g. `CCTT/CCCT`) — the first is taken
//! as canonical. Contig names are bare (`1`, `Y`); callers normalize against the BAM's naming.

use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use flate2::read::MultiGzDecoder;

use crate::error::AnalysisError;

/// One STR locus: a tight repeat tract with its period (motif length) and reference copy number.
#[derive(Debug, Clone, PartialEq)]
pub struct StrLocus {
    /// Contig as written in the reference BED (bare — `1`, `X`, `Y`).
    pub contig: String,
    /// 0-based, half-open tract start (BED convention).
    pub start: i64,
    /// Exclusive tract end. Tract length in bp = `end - start`.
    pub end: i64,
    /// Repeat-unit length (motif size) in bp.
    pub period: u8,
    /// Copy number in the reference allele (can be fractional — a partial final unit).
    pub ref_copies: f64,
    /// Locus id (HipSTR `Human_STR_N`), used as the result name until a vendor mapping exists.
    pub name: String,
    /// Canonical repeat motif (the first when the BED lists `A/B` alternatives).
    pub motif: String,
}

impl StrLocus {
    /// Whether `name`/the BED contig matches `query` after stripping an optional `chr` prefix on
    /// either side (the BAM may be `chrY`, the BED `Y`; see the contig-naming convention).
    pub fn contig_matches(&self, query: &str) -> bool {
        let strip = |s: &str| s.strip_prefix("chr").unwrap_or(s).to_string();
        strip(&self.contig).eq_ignore_ascii_case(&strip(query))
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
    Some(StrLocus { contig, start, end, period, ref_copies, name, motif })
}

/// Read STR loci from a (gzipped) HipSTR reference BED, keeping only those on `contig` (matched
/// prefix-insensitively) with `period >= min_period`. Filtering while streaming avoids holding the
/// genome-wide ~1.6M-locus set in memory when only one chromosome is needed. Results are sorted by
/// start. `min_period` of 2 drops homopolymers (period 1) — noisy and not genealogical markers.
pub fn load_hipstr_contig(
    bed_gz: &Path,
    contig: &str,
    min_period: u8,
) -> Result<Vec<StrLocus>, AnalysisError> {
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
