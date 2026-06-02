//! Genotyping-array (chip) profiles — the QC summary of a vendor raw-data export
//! (23andMe, AncestryDNA, MyHeritage, …), a pragmatic port of the Scala `ChipProfile`.
//! We don't keep every genotype (a chip is ~600–700k markers); we keep the call/no-call/
//! het summary and per-region counts that drive quality and downstream eligibility.
//! [`summarize`] is a pure pass over the file text (no IO) that also guesses the vendor.

use du_domain::ids::SampleGuid;
use serde::{Deserialize, Serialize};

/// Quality summary computed from a chip's genotype calls.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct ChipSummary {
    pub total_markers_possible: i64,
    pub total_markers_called: i64,
    /// Fraction of markers with no call (0.0–1.0).
    pub no_call_rate: f64,
    /// Heterozygosity rate over *called autosomal* markers, if any.
    pub het_rate: Option<f64>,
    pub y_markers_called: i64,
    pub mt_markers_called: i64,
    pub autosomal_markers_called: i64,
}

/// A subject's chip profile (QC summary + provenance).
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ChipProfile {
    pub id: i64,
    pub biosample_guid: SampleGuid,
    /// Vendor (one of [`KNOWN_PROVIDERS`]).
    pub provider: String,
    pub chip_version: Option<String>,
    pub summary: ChipSummary,
    pub source_file_name: Option<String>,
}

/// Fields for creating a chip profile (the store assigns the id).
#[derive(Debug, Clone, PartialEq)]
pub struct NewChipProfile {
    pub biosample_guid: SampleGuid,
    pub provider: String,
    pub chip_version: Option<String>,
    pub summary: ChipSummary,
    pub source_file_name: Option<String>,
}

/// Known array vendors (for the import form's dropdown).
pub const KNOWN_PROVIDERS: &[&str] =
    &["23andMe", "AncestryDNA", "MyHeritage", "FTDNA", "LivingDNA", "OTHER"];

/// A genotype call's zygosity, derived only from the called bases.
enum Zygosity {
    NoCall,
    Hom,
    Het,
}

/// Classify a genotype token (e.g. "AA", "AG", "--", "00", "DI"). Non-A/C/G/T characters
/// are ignored, so any all-symbol token (no-call, indel) classifies as a no-call.
fn classify(genotype: &str) -> Zygosity {
    let bases: Vec<u8> = genotype
        .trim()
        .trim_matches('"')
        .bytes()
        .map(|b| b.to_ascii_uppercase())
        .filter(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
        .collect();
    match bases.len() {
        1 => Zygosity::Hom, // haploid call (Y/MT)
        2 if bases[0] == bases[1] => Zygosity::Hom,
        2 => Zygosity::Het,
        _ => Zygosity::NoCall,
    }
}

enum Region {
    Autosomal,
    Y,
    Mt,
    Other,
}

fn region(chrom: &str) -> Region {
    let lc = chrom.trim().trim_matches('"').to_ascii_lowercase();
    let core = lc.strip_prefix("chr").unwrap_or(&lc);
    match core {
        "y" | "24" => Region::Y,
        "mt" | "m" | "26" => Region::Mt,
        c if c.parse::<u32>().map(|n| (1..=22).contains(&n)).unwrap_or(false) => Region::Autosomal,
        _ => Region::Other, // X (23) and anything else
    }
}

/// Is this row a header (`rsid …`) rather than data?
fn is_header(first_field: &str) -> bool {
    let f = first_field.trim().trim_matches('"').to_ascii_lowercase();
    f == "rsid" || f == "rs_id" || f == "snp" || f == "#rsid"
}

/// Guess the vendor from the file's header/comment text.
fn detect_provider(lower_header: &str) -> Option<String> {
    let has = |s: &str| lower_header.contains(s);
    if has("23andme") {
        Some("23andMe".into())
    } else if has("ancestrydna") || has("ancestry.com") {
        Some("AncestryDNA".into())
    } else if has("myheritage") {
        Some("MyHeritage".into())
    } else if has("familytreedna") || has("ftdna") {
        Some("FTDNA".into())
    } else if has("livingdna") || has("living dna") {
        Some("LivingDNA".into())
    } else if lower_header.contains("allele1") && lower_header.contains("allele2") {
        Some("AncestryDNA".into()) // allele1/allele2 layout
    } else {
        None
    }
}

/// Summarize a vendor raw-data export into QC metrics and a guessed provider. Accepts the
/// common layouts: tab- or comma-separated, optional `#` comment header, then either
/// `rsid,chrom,pos,genotype` (23andMe/MyHeritage) or `rsid,chrom,pos,allele1,allele2`
/// (AncestryDNA). Errors only if no marker rows are found.
pub fn summarize(text: &str) -> Result<(ChipSummary, Option<String>), String> {
    let mut possible = 0i64;
    let mut called = 0i64;
    let (mut y, mut mt, mut auto) = (0i64, 0i64, 0i64);
    let (mut auto_called, mut auto_het) = (0i64, 0i64);
    let mut header_lower = String::new();

    for raw in text.lines() {
        let line = raw.trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('#') {
            if header_lower.len() < 4096 {
                header_lower.push_str(&line.to_ascii_lowercase());
                header_lower.push('\n');
            }
            continue;
        }
        let sep = if line.contains('\t') { '\t' } else { ',' };
        let cols: Vec<&str> = line.split(sep).map(|s| s.trim().trim_matches('"')).collect();
        if cols.len() < 4 {
            continue;
        }
        if is_header(cols[0]) {
            header_lower.push_str(&line.to_ascii_lowercase());
            header_lower.push('\n');
            continue;
        }
        let chrom = cols[1];
        // genotype is either one column (4-col) or two allele columns (>=5-col).
        let genotype = if cols.len() >= 5 {
            format!("{}{}", cols[3], cols[4])
        } else {
            cols[3].to_string()
        };

        possible += 1;
        let call = classify(&genotype);
        let is_called = !matches!(call, Zygosity::NoCall);
        if is_called {
            called += 1;
        }
        match region(chrom) {
            Region::Y => {
                if is_called {
                    y += 1;
                }
            }
            Region::Mt => {
                if is_called {
                    mt += 1;
                }
            }
            Region::Autosomal => {
                if is_called {
                    auto += 1;
                    auto_called += 1;
                    if matches!(call, Zygosity::Het) {
                        auto_het += 1;
                    }
                }
            }
            Region::Other => {}
        }
    }

    if possible == 0 {
        return Err("no marker rows found (expected rsid,chrom,pos,genotype)".into());
    }
    let summary = ChipSummary {
        total_markers_possible: possible,
        total_markers_called: called,
        no_call_rate: (possible - called) as f64 / possible as f64,
        het_rate: (auto_called > 0).then(|| auto_het as f64 / auto_called as f64),
        y_markers_called: y,
        mt_markers_called: mt,
        autosomal_markers_called: auto,
    };
    Ok((summary, detect_provider(&header_lower)))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn summarizes_23andme_tsv() {
        let f = "# This data file generated by 23andMe\nrsid\tchromosome\tposition\tgenotype\n\
                 rs1\t1\t100\tAA\nrs2\t1\t200\tAG\nrs3\t1\t300\t--\nrs4\tY\t400\tG\nrs5\tMT\t500\tT\n";
        let (s, provider) = summarize(f).unwrap();
        assert_eq!(provider.as_deref(), Some("23andMe"));
        assert_eq!(s.total_markers_possible, 5);
        assert_eq!(s.total_markers_called, 4); // the "--" is a no-call
        assert_eq!(s.autosomal_markers_called, 2);
        assert_eq!(s.y_markers_called, 1);
        assert_eq!(s.mt_markers_called, 1);
        assert_eq!(s.het_rate, Some(0.5)); // 1 het of 2 autosomal
        assert!((s.no_call_rate - 0.2).abs() < 1e-9);
    }

    #[test]
    fn summarizes_ancestry_allele_columns() {
        let f = "#AncestryDNA raw data download\nrsid\tchromosome\tposition\tallele1\tallele2\n\
                 rs1\t1\t100\tA\tA\nrs2\t2\t200\tA\tG\nrs3\t23\t300\t0\t0\n";
        let (s, provider) = summarize(f).unwrap();
        assert_eq!(provider.as_deref(), Some("AncestryDNA"));
        assert_eq!(s.total_markers_possible, 3);
        assert_eq!(s.total_markers_called, 2); // 0/0 is a no-call; X not counted in autosomal
        assert_eq!(s.autosomal_markers_called, 2);
    }

    #[test]
    fn empty_errors() {
        assert!(summarize("# only comments\n\n").is_err());
    }
}
