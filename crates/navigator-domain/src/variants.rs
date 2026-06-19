//! A subject's SNP variant calls, imported from a VCF or a CSV/TSV table and grouped into
//! a named [`VariantSet`] (Scala's `DataType.Variants`). A pragmatic port of the Scala
//! `VariantCall`: the columns a basic VCF/CSV carries — contig, position, ref/alt, rsID,
//! and (CSV only) a genotype — without QUAL/depth, which the shared VCF parser doesn't
//! surface yet. Types are pure; [`parse_csv`] turns a marker table into calls with no IO.

use du_domain::ids::SampleGuid;
use serde::{Deserialize, Serialize};

/// A single biallelic SNP call.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct VariantCall {
    pub contig: String,
    /// 1-based position.
    pub position: i64,
    pub reference: String,
    pub alternate: String,
    /// dbSNP rsID (or other variant id), if known.
    pub rs_id: Option<String>,
    /// Genotype call (e.g. "0/1", "1/1"), if the source provides one.
    pub genotype: Option<String>,
}

/// The kind of source a variant set came from — carries the SNP-concordance weight used
/// when reconciling across sources (Scala `YProfileSourceType`). Sanger is the gold
/// standard (1.0); a low-confidence manual entry is 0.3.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SourceType {
    Sanger,
    WgsLongRead,
    WgsShortRead,
    TargetedNgs,
    Chip,
    Manual,
    /// VCF/CSV of unknown provenance.
    Imported,
}

impl SourceType {
    /// SNP-concordance weight (0–1): reliability for SNP calls.
    pub fn snp_weight(self) -> f64 {
        match self {
            SourceType::Sanger => 1.0,
            SourceType::WgsLongRead => 0.95,
            SourceType::WgsShortRead => 0.85,
            SourceType::TargetedNgs => 0.75,
            SourceType::Chip => 0.5,
            SourceType::Manual => 0.3,
            SourceType::Imported => 0.7,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            SourceType::Sanger => "SANGER",
            SourceType::WgsLongRead => "WGS_LONG_READ",
            SourceType::WgsShortRead => "WGS_SHORT_READ",
            SourceType::TargetedNgs => "TARGETED_NGS",
            SourceType::Chip => "CHIP",
            SourceType::Manual => "MANUAL",
            SourceType::Imported => "IMPORTED",
        }
    }

    pub fn from_code(s: &str) -> SourceType {
        match s {
            "SANGER" => SourceType::Sanger,
            "WGS_LONG_READ" => SourceType::WgsLongRead,
            "WGS_SHORT_READ" => SourceType::WgsShortRead,
            "TARGETED_NGS" => SourceType::TargetedNgs,
            "CHIP" => SourceType::Chip,
            "MANUAL" => SourceType::Manual,
            _ => SourceType::Imported,
        }
    }

    /// All types, for a UI dropdown.
    pub const ALL: &'static [SourceType] = &[
        SourceType::Sanger,
        SourceType::WgsLongRead,
        SourceType::WgsShortRead,
        SourceType::TargetedNgs,
        SourceType::Chip,
        SourceType::Manual,
        SourceType::Imported,
    ];
}

/// A subject's variant calls from one import (a VCF, CSV export, YSEQ/Sanger panel, …).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct VariantSet {
    pub id: i64,
    pub biosample_guid: SampleGuid,
    /// A label for the source (typically the file name).
    pub source_label: String,
    pub source_type: SourceType,
    /// Reference build the call positions are on (`"hs1"`, `"GRCh38"`, …), when known. `None`
    /// for sources of unknown build (a generic VCF/CSV import). Lets build-specific consumers
    /// (e.g. Y-SNP-panel placement) read the build directly instead of re-deriving it.
    pub reference_build: Option<String>,
    pub calls: Vec<VariantCall>,
}

/// Fields for creating a variant set (the store assigns the id).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NewVariantSet {
    pub biosample_guid: SampleGuid,
    pub source_label: String,
    pub source_type: SourceType,
    /// Reference build the calls are on, when known (see [`VariantSet::reference_build`]).
    pub reference_build: Option<String>,
    pub calls: Vec<VariantCall>,
}

/// True for a one-base A/C/G/T allele (case-insensitive) — used to keep SNP rows only.
fn is_snp_allele(a: &str) -> bool {
    a.len() == 1 && matches!(a.as_bytes()[0].to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

/// Build a SNP `VariantCall`, returning `None` for indels/symbolic alleles.
pub fn snp_call(
    contig: &str,
    position: i64,
    reference: &str,
    alternate: &str,
    rs_id: Option<String>,
    genotype: Option<String>,
) -> Option<VariantCall> {
    (is_snp_allele(reference) && is_snp_allele(alternate)).then(|| VariantCall {
        contig: contig.to_string(),
        position,
        reference: reference.to_ascii_uppercase(),
        alternate: alternate.to_ascii_uppercase(),
        rs_id: rs_id.filter(|s| !s.is_empty() && s != "."),
        genotype: genotype.filter(|s| !s.is_empty() && s != "."),
    })
}

/// Which CSV column holds each field. Positions are 0-based column indices.
struct Layout {
    contig: usize,
    position: usize,
    reference: usize,
    alternate: usize,
    rs_id: Option<usize>,
    genotype: Option<usize>,
}

impl Layout {
    /// Default positional layout: contig, position, ref, alt, [rsid], [genotype].
    fn positional() -> Self {
        Layout {
            contig: 0,
            position: 1,
            reference: 2,
            alternate: 3,
            rs_id: Some(4),
            genotype: Some(5),
        }
    }

    /// Map columns by a recognized header row, or `None` if the row isn't a header.
    fn from_header(cols: &[&str]) -> Option<Self> {
        let find = |names: &[&str]| {
            cols.iter().position(|c| {
                let c = c.trim().to_ascii_lowercase();
                names.contains(&c.as_str())
            })
        };
        let contig = find(&["contig", "chrom", "chromosome", "#chrom"])?;
        let position = find(&["position", "pos"])?;
        let reference = find(&["reference", "ref"])?;
        let alternate = find(&["alternate", "alt"])?;
        Some(Layout {
            contig,
            position,
            reference,
            alternate,
            rs_id: find(&["rsid", "rs_id", "id"]),
            genotype: find(&["genotype", "gt"]),
        })
    }
}

/// Parse a CSV/TSV variant table into SNP calls. The first non-comment row is treated as a
/// header when it names known columns (contig/pos/ref/alt[/rsid/genotype], any order),
/// otherwise columns are read positionally as contig,position,ref,alt[,rsid][,genotype].
/// Non-SNP rows and rows with an unparseable position are skipped. Errors if none parse.
pub fn parse_csv(text: &str) -> Result<Vec<VariantCall>, String> {
    let mut rows = text
        .lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with("##"));

    let Some(first) = rows.next() else {
        return Err("no variant rows found".into());
    };
    let sep = if first.contains('\t') { '\t' } else { ',' };
    let split = |line: &str| line.split(sep).map(|s| s.trim().to_string()).collect::<Vec<_>>();

    let first_cols: Vec<&str> = first.split(sep).map(str::trim).collect();
    let layout = Layout::from_header(&first_cols);
    let mut calls = Vec::new();
    // If the first row wasn't a header, it's data — parse it positionally too.
    let header_layout = match layout {
        Some(l) => l,
        None => {
            push_row(&mut calls, &split(first), &Layout::positional());
            Layout::positional()
        }
    };
    for line in rows {
        push_row(&mut calls, &split(line), &header_layout);
    }
    if calls.is_empty() {
        return Err("no SNP variant rows parsed (expected contig,position,ref,alt)".into());
    }
    Ok(calls)
}

fn push_row(out: &mut Vec<VariantCall>, cols: &[String], l: &Layout) {
    let get = |i: usize| cols.get(i).map(String::as_str).unwrap_or("");
    let Ok(position) = get(l.position).parse::<i64>() else {
        return;
    };
    let opt_at = |idx: Option<usize>| idx.map(|i| get(i).to_string()).filter(|s| !s.is_empty());
    if let Some(call) = snp_call(
        get(l.contig),
        position,
        get(l.reference),
        get(l.alternate),
        opt_at(l.rs_id),
        opt_at(l.genotype),
    ) {
        out.push(call);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_header_csv_and_keeps_only_snps() {
        let csv = "contig,position,ref,alt,rsid\nchr1,1000,A,G,rs1\nchr1,2000,A,AT,rs2\nchrM,73,G,A,.\n";
        let v = parse_csv(csv).unwrap();
        assert_eq!(v.len(), 2); // the A>AT indel is dropped
        assert_eq!(
            v[0],
            VariantCall {
                contig: "chr1".into(),
                position: 1000,
                reference: "A".into(),
                alternate: "G".into(),
                rs_id: Some("rs1".into()),
                genotype: None
            }
        );
        assert_eq!(v[1].rs_id, None); // "." normalized away
    }

    #[test]
    fn parses_positional_tsv_with_genotype() {
        let tsv = "chr1\t100\tC\tT\trs9\t0/1\nchr2\t200\tG\tA\t.\t1/1\n";
        let v = parse_csv(tsv).unwrap();
        assert_eq!(v.len(), 2);
        assert_eq!(v[0].genotype.as_deref(), Some("0/1"));
        assert_eq!(v[1].genotype.as_deref(), Some("1/1"));
    }

    #[test]
    fn errors_when_no_snps() {
        assert!(parse_csv("contig,position,ref,alt\nchr1,10,A,ACGT\n").is_err());
    }
}
