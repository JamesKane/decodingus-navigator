//! FTDNA Big Y CSV variant reports — the "lesser access" substitute for the BAM/CRAM/VCF, for
//! project admins whose access tier exposes only the browser CSV exports. Two report flavors,
//! both **GRCh38 chrY** derived-allele calls:
//!
//!   Named Variants:   `SNP_Name,Position,On_Haplotree,Ancestral,Derived`
//!   Private Variants: `Position,Ancestral,Derived`
//!
//! Each row is a position where the sample carries the **Derived** allele (a positive call), so we
//! emit a SNP [`VariantCall`] with `reference = ancestral`, `alternate = derived`, `genotype = "1"`
//! (derived), and `rs_id` = the SNP name when present. SNP-only (single-base ACGT); other rows are
//! skipped. The calls land in GRCh38 space, which is FTDNA's native Y-tree build — so Y placement
//! matches positions directly, no liftover.

use serde::{Deserialize, Serialize};

use crate::variants::{self, VariantCall};

/// Which FTDNA Big Y CSV report a file is.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum FtdnaReport {
    /// Named (haplotree/known) SNPs the sample is derived for.
    Named,
    /// Novel/unnamed variants (position only) the sample is derived for.
    Private,
}

impl FtdnaReport {
    pub fn label(self) -> &'static str {
        match self {
            FtdnaReport::Named => "FTDNA Big Y Named Variants",
            FtdnaReport::Private => "FTDNA Big Y Private Variants",
        }
    }
}

/// chrY contig label for the emitted calls — matches the FTDNA GRCh38 Y-tree's contig.
const CONTIG: &str = "chrY";

/// Split a CSV line into trimmed, unquoted cells.
fn cells(line: &str) -> Vec<String> {
    line.split(',').map(|s| s.trim().trim_matches('"').to_string()).collect()
}

/// Recognize the report flavor from a header row's columns, or `None` if it isn't an FTDNA Big Y
/// Named/Private Variants header. Case-insensitive; quotes already stripped by [`cells`].
pub fn report_of_header(cols: &[String]) -> Option<FtdnaReport> {
    let norm: Vec<String> = cols.iter().map(|c| c.to_ascii_lowercase()).collect();
    let has = |n: &str| norm.iter().any(|c| c == n);
    if has("snp_name") && has("position") && has("ancestral") && has("derived") {
        Some(FtdnaReport::Named)
    } else if norm.len() == 3 && norm[0] == "position" && norm[1] == "ancestral" && norm[2] == "derived" {
        Some(FtdnaReport::Private)
    } else {
        None
    }
}

/// True if `text`'s first non-empty line is an FTDNA Big Y Named/Private Variants header.
pub fn looks_like_ftdna_variant_csv(text: &str) -> bool {
    text.lines()
        .map(str::trim)
        .find(|l| !l.is_empty())
        .map(|h| report_of_header(&cells(h)).is_some())
        .unwrap_or(false)
}

/// Parse an FTDNA Big Y Named/Private Variants CSV into chrY derived-allele SNP calls, returning
/// the report flavor alongside. Errors if the header isn't a recognized FTDNA report or no SNP
/// rows parse.
pub fn parse(text: &str) -> Result<(FtdnaReport, Vec<VariantCall>), String> {
    let mut lines = text.lines().map(str::trim).filter(|l| !l.is_empty());
    let header = lines.next().ok_or("empty FTDNA variant CSV")?;
    let hcols = cells(header);
    let report = report_of_header(&hcols)
        .ok_or("not an FTDNA Big Y Named/Private Variants CSV (unrecognized header)")?;

    let col = |name: &str| hcols.iter().position(|c| c.eq_ignore_ascii_case(name));
    let i_name = col("SNP_Name");
    let i_pos = col("Position").ok_or("missing Position column")?;
    let i_anc = col("Ancestral").ok_or("missing Ancestral column")?;
    let i_der = col("Derived").ok_or("missing Derived column")?;

    let mut calls = Vec::new();
    for line in lines {
        let c = cells(line);
        let get = |i: usize| c.get(i).map(String::as_str).unwrap_or("");
        let Ok(position) = get(i_pos).parse::<i64>() else { continue };
        let name = i_name.map(|i| get(i).to_string()).filter(|s| !s.is_empty());
        // Each row is a derived (positive) call: ref = ancestral, alt = derived, gt = "1".
        if let Some(call) = variants::snp_call(CONTIG, position, get(i_anc), get(i_der), name, Some("1".into())) {
            calls.push(call);
        }
    }
    if calls.is_empty() {
        return Err("FTDNA variant CSV header recognized but no SNP rows parsed".into());
    }
    Ok((report, calls))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_named_variants() {
        let csv = "\"SNP_Name\",\"Position\",\"On_Haplotree\",\"Ancestral\",\"Derived\"\n\
                   \"M9056\",\"10002698\",\"No\",\"A\",\"G\"\n\
                   \"PF5885\",\"10004091\",\"Yes\",\"G\",\"C\"\n";
        let (report, calls) = parse(csv).unwrap();
        assert_eq!(report, FtdnaReport::Named);
        assert_eq!(calls.len(), 2);
        let pf = &calls[1];
        assert_eq!(pf.contig, "chrY");
        assert_eq!(pf.position, 10004091);
        assert_eq!(pf.reference, "G"); // ancestral
        assert_eq!(pf.alternate, "C"); // derived
        assert_eq!(pf.rs_id.as_deref(), Some("PF5885"));
        assert_eq!(pf.genotype.as_deref(), Some("1"));
    }

    #[test]
    fn parses_private_variants() {
        let csv = "\"Position\",\"Ancestral\",\"Derived\"\n\
                   \"13603685\",\"C\",\"T\"\n\
                   \"15487200\",\"C\",\"A\"\n";
        let (report, calls) = parse(csv).unwrap();
        assert_eq!(report, FtdnaReport::Private);
        assert_eq!(calls.len(), 2);
        assert_eq!(calls[0].position, 13603685);
        assert_eq!(calls[0].reference, "C");
        assert_eq!(calls[0].alternate, "T");
        assert!(calls[0].rs_id.is_none());
    }

    #[test]
    fn skips_non_snp_rows_and_rejects_foreign_headers() {
        // An indel row is dropped; the SNP row survives.
        let csv = "SNP_Name,Position,On_Haplotree,Ancestral,Derived\n\
                   ins1,100,No,A,AT\n\
                   M1,200,Yes,C,T\n";
        let (_, calls) = parse(csv).unwrap();
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].position, 200);

        // A generic variant table (has contig) is not an FTDNA report.
        assert!(!looks_like_ftdna_variant_csv("contig,position,ref,alt\nchrY,100,A,G\n"));
        assert!(parse("contig,position,ref,alt\nchrY,100,A,G\n").is_err());
    }
}
