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
    /// Absolute path of the imported raw-data file, for re-reading the autosomal genotypes on
    /// demand (ancestry). `None` for older rows imported before this was tracked.
    pub source_path: Option<String>,
}

/// Fields for creating a chip profile (the store assigns the id).
#[derive(Debug, Clone, PartialEq)]
pub struct NewChipProfile {
    pub biosample_guid: SampleGuid,
    pub provider: String,
    pub chip_version: Option<String>,
    pub summary: ChipSummary,
    pub source_file_name: Option<String>,
    pub source_path: Option<String>,
}

/// Known array vendors (for the import form's dropdown).
pub const KNOWN_PROVIDERS: &[&str] = &["23andMe", "AncestryDNA", "MyHeritage", "FTDNA", "LivingDNA", "OTHER"];

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

/// Which haploid lineage a [`ChipHaploCall`] belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChipDna {
    Y,
    Mt,
}

/// A single haploid Y or mtDNA genotype pulled from a chip export — the raw observed allele on
/// the vendor's reference build, for on-import haplogroup placement. (Consumer arrays report
/// Y/MT as a single haploid base; we keep only unambiguous single-base calls.)
#[derive(Debug, Clone, PartialEq)]
pub struct ChipHaploCall {
    pub dna: ChipDna,
    pub rsid: String,
    pub position: i64,
    /// Observed allele, uppercase A/C/G/T.
    pub base: char,
}

/// The single haploid base of a genotype token, or `None` if it's a no-call, an indel
/// (`I`/`D`), or heterozygous (two different bases — on a true haploid Y/MT that's
/// contamination, so we drop it rather than guess).
fn haploid_base(genotype: &str) -> Option<char> {
    let mut bases = genotype
        .bytes()
        .map(|b| b.to_ascii_uppercase())
        .filter(|b| matches!(b, b'A' | b'C' | b'G' | b'T'));
    let first = bases.next()?;
    bases.all(|b| b == first).then_some(first as char)
}

/// Extract the Y and mtDNA haploid calls from a vendor raw-data export, for on-import
/// haplogroup placement. Skips autosomal/X rows, no-calls, indels, and heterozygous calls.
/// Positions are on the vendor build (consumer arrays are GRCh37 — see [`detect_build`]).
/// Pairs with [`summarize`]: same row layouts (tab/comma, optional `#` header, then
/// `rsid,chrom,pos,genotype` or `rsid,chrom,pos,allele1,allele2`).
pub fn haplo_calls(text: &str) -> Vec<ChipHaploCall> {
    let mut out = Vec::new();
    for raw in text.lines() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let sep = if line.contains('\t') { '\t' } else { ',' };
        let cols: Vec<&str> = line.split(sep).map(|s| s.trim().trim_matches('"')).collect();
        if cols.len() < 4 || is_header(cols[0]) {
            continue;
        }
        let dna = match region(cols[1]) {
            Region::Y => ChipDna::Y,
            Region::Mt => ChipDna::Mt,
            _ => continue,
        };
        let Ok(position) = cols[2].parse::<i64>() else { continue };
        let genotype = if cols.len() >= 5 {
            format!("{}{}", cols[3], cols[4])
        } else {
            cols[3].to_string()
        };
        let Some(base) = haploid_base(&genotype) else { continue };
        out.push(ChipHaploCall {
            dna,
            rsid: cols[0].to_string(),
            position,
            base,
        });
    }
    out
}

/// A single autosomal diploid genotype from a chip export — the two observed alleles at a SNP, on
/// the vendor build (GRCh37). Fed (after liftover to the panel build) into the ancestry estimators.
#[derive(Debug, Clone, PartialEq)]
pub struct ChipAutosomalCall {
    /// Chromosome, normalized to `chr1`..`chr22`.
    pub contig: String,
    pub position: i64,
    pub a1: char,
    pub a2: char,
}

/// The two A/C/G/T bases of a diploid genotype token (`"AG"`, or two allele columns joined), or
/// `None` for a no-call / indel / not-exactly-two-bases token.
fn diploid_bases(genotype: &str) -> Option<(char, char)> {
    let bases: Vec<char> = genotype
        .bytes()
        .map(|b| b.to_ascii_uppercase())
        .filter(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
        .map(|b| b as char)
        .collect();
    (bases.len() == 2).then(|| (bases[0], bases[1]))
}

/// Extract the **autosomal** diploid SNP calls from a vendor raw-data export, for ancestry. Keeps
/// only chr1–22, called, biallelic-SNP rows (drops Y/MT/X, no-calls, indels). Same row layouts as
/// [`summarize`]/[`haplo_calls`]. Positions are on the vendor build (GRCh37 — see [`detect_build`]).
pub fn autosomal_calls(text: &str) -> Vec<ChipAutosomalCall> {
    let mut out = Vec::new();
    for raw in text.lines() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let sep = if line.contains('\t') { '\t' } else { ',' };
        let cols: Vec<&str> = line.split(sep).map(|s| s.trim().trim_matches('"')).collect();
        if cols.len() < 4 || is_header(cols[0]) {
            continue;
        }
        if !matches!(region(cols[1]), Region::Autosomal) {
            continue;
        }
        // Normalize the chromosome to chrN (1..22) — matches the CHM13 panel + liftover contig naming.
        let core = cols[1].trim().trim_matches('"').to_ascii_lowercase();
        let core = core.strip_prefix("chr").unwrap_or(&core);
        let Ok(position) = cols[2].parse::<i64>() else { continue };
        let genotype = if cols.len() >= 5 {
            format!("{}{}", cols[3], cols[4])
        } else {
            cols[3].to_string()
        };
        let Some((a1, a2)) = diploid_bases(&genotype) else {
            continue;
        };
        out.push(ChipAutosomalCall {
            contig: format!("chr{core}"),
            position,
            a1,
            a2,
        });
    }
    out
}

/// The reference build a vendor export is reported on. Consumer arrays (23andMe v4/v5,
/// AncestryDNA v1/v2) are GRCh37, so that's the default; a header naming build 38 / GRCh38 /
/// hg38 overrides it. Scans only the comment header.
pub fn detect_build(text: &str) -> String {
    for raw in text.lines() {
        let line = raw.trim();
        if line.starts_with('#') {
            let lc = line.to_ascii_lowercase();
            if lc.contains("build 38") || lc.contains("grch38") || lc.contains("hg38") {
                return "GRCh38".into();
            }
        } else if !line.is_empty() {
            break; // past the header block
        }
    }
    "GRCh37".into()
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

    #[test]
    fn autosomal_calls_keeps_only_called_autosomal_snps() {
        // 23andMe 4-col + AncestryDNA 5-col rows, mixed with Y/MT/X/no-call/indel to drop.
        let f = "rsid\tchromosome\tposition\tgenotype\n\
                 rs1\t1\t100\tAG\n\
                 rs2\t22\t200\tCC\n\
                 rs3\t1\t300\t--\n\
                 rs4\tY\t400\tG\n\
                 rs5\tMT\t500\tT\n\
                 rs6\t23\t600\tAA\n\
                 rsI\t2\t700\tII\n";
        let calls = autosomal_calls(f);
        assert_eq!(calls.len(), 2);
        assert_eq!(
            calls[0],
            ChipAutosomalCall {
                contig: "chr1".into(),
                position: 100,
                a1: 'A',
                a2: 'G'
            }
        );
        assert_eq!(
            calls[1],
            ChipAutosomalCall {
                contig: "chr22".into(),
                position: 200,
                a1: 'C',
                a2: 'C'
            }
        );
    }

    #[test]
    fn autosomal_calls_ancestry_allele_columns() {
        let f = "#AncestryDNA\nrsid\tchromosome\tposition\tallele1\tallele2\n\
                 rs1\t5\t1000\tA\tG\nrs2\t5\t2000\t0\t0\n";
        let calls = autosomal_calls(f);
        assert_eq!(calls.len(), 1);
        assert_eq!(
            calls[0],
            ChipAutosomalCall {
                contig: "chr5".into(),
                position: 1000,
                a1: 'A',
                a2: 'G'
            }
        );
    }

    #[test]
    fn haplo_calls_extracts_y_and_mt_haploid_bases() {
        // autosomal + X rows ignored; Y/MT haploid kept; het Y dropped; no-call/indel dropped.
        let f = "# 23andMe\nrsid\tchromosome\tposition\tgenotype\n\
                 rs1\t1\t100\tAG\n\
                 rs2\tX\t200\tA\n\
                 rsY1\tY\t2800000\tG\n\
                 rsY2\t24\t2900000\tCC\n\
                 rsYhet\tY\t3000000\tAT\n\
                 rsYnc\tY\t3100000\t--\n\
                 rsM1\tMT\t263\tG\n\
                 rsM2\t26\t750\tA\n\
                 rsMii\tMT\t900\tII\n";
        let calls = haplo_calls(f);
        let ys: Vec<_> = calls.iter().filter(|c| c.dna == ChipDna::Y).collect();
        let mts: Vec<_> = calls.iter().filter(|c| c.dna == ChipDna::Mt).collect();
        assert_eq!(ys.len(), 2, "two valid Y haploid calls (chr Y + chr 24)");
        assert_eq!(
            ys[0],
            &ChipHaploCall {
                dna: ChipDna::Y,
                rsid: "rsY1".into(),
                position: 2_800_000,
                base: 'G'
            }
        );
        assert_eq!(ys[1].base, 'C'); // "CC" homozygous → C
        assert_eq!(mts.len(), 2, "two valid MT haploid calls (chr MT + chr 26)");
        assert_eq!(
            mts[0],
            &ChipHaploCall {
                dna: ChipDna::Mt,
                rsid: "rsM1".into(),
                position: 263,
                base: 'G'
            }
        );
    }

    #[test]
    fn haplo_calls_handles_ancestry_allele_columns() {
        let f = "#AncestryDNA\nrsid\tchromosome\tposition\tallele1\tallele2\n\
                 rsY\t24\t2800000\tA\tA\nrsX\t23\t100\tC\tT\n";
        let calls = haplo_calls(f);
        assert_eq!(
            calls,
            vec![ChipHaploCall {
                dna: ChipDna::Y,
                rsid: "rsY".into(),
                position: 2_800_000,
                base: 'A'
            }]
        );
    }

    #[test]
    fn detect_build_defaults_grch37_and_honors_header() {
        assert_eq!(
            detect_build("# This data file generated by 23andMe\nrsid\t..\n"),
            "GRCh37"
        );
        assert_eq!(
            detect_build("# reference human assembly build 38\nrsid\t..\n"),
            "GRCh38"
        );
        assert_eq!(detect_build("#GRCh38 export\n"), "GRCh38");
    }
}
