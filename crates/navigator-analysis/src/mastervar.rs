//! CompleteGenomics **masterVar** reader — the whole-genome variant table `cgatools` emitted
//! for the old CG sequencing service (`var-<slide>-ASM.tsv[.bz2]`, `FORMAT_VERSION 2.0`).
//!
//! The file is a tab-separated per-allele call over the *entire* genome: `ref` / `no-call`
//! blocks interleaved with the actual `snp` / `ins` / `del` / `sub` calls. Each locus is one
//! or more consecutive rows sharing the leading `locus` id — a diploid heterozygous SNP is two
//! rows (allele `1` and allele `2`), a homozygous or haploid call one row (allele `all` / `1`).
//! Coordinates are **0-based half-open** (`begin`/`end`); the 1-based SNP position is `begin + 1`.
//!
//! We extract SNPs only (matching the VCF/CSV variant importer, which is SNP-only): each locus
//! becomes at most one biallelic [`VariantCall`] with a reconstructed genotype (`1`, `0/1`,
//! `1/1`, `1/2`, or `1/.` when the partner allele is a no-call). Indels/substitutions and the
//! `ref`/`no-call` spans are skipped. NCBI build 37 ⇒ the calls are on **GRCh37** (chrM = rCRS).
//!
//! The whole file is streamed (a genome masterVar is multi-GB uncompressed), transparently
//! decoding `.bz2` / `.gz` via [`crate::gzio::open_maybe_compressed`].

use std::io::BufRead;
use std::path::Path;

use navigator_domain::variants::{snp_call, VariantCall};

use crate::gzio;

#[derive(Debug, thiserror::Error)]
pub enum MasterVarError {
    #[error(transparent)]
    Io(#[from] std::io::Error),
    #[error("not a CompleteGenomics masterVar file: {0}")]
    Format(String),
}

/// The result of reading a masterVar file: the SNP calls plus the sample id / reference build
/// parsed from the `#`-comment header, and a couple of tallies for the import summary.
#[derive(Debug, Default, Clone)]
pub struct MasterVarImport {
    /// The `#SAMPLE` header value, if present (e.g. `GS00253-DNA_A01`).
    pub sample_id: Option<String>,
    /// Reference build the calls are on — `"GRCh37"` for `NCBI build 37` (the only build the CG
    /// service shipped), else a best-effort echo of the `#GENOME_REFERENCE` header.
    pub reference_build: String,
    /// The extracted biallelic SNP calls (one per SNP locus).
    pub calls: Vec<VariantCall>,
    /// Total loci examined (distinct `locus` ids).
    pub loci_seen: u64,
    /// Loci that yielded a SNP call.
    pub snp_loci: u64,
}

/// Column indices for the fields we read, resolved from the `>locus ploidy allele …` header.
struct Columns {
    locus: usize,
    ploidy: usize,
    allele: usize,
    chromosome: usize,
    begin: usize,
    var_type: usize,
    reference: usize,
    allele_seq: usize,
    x_ref: Option<usize>,
}

impl Columns {
    /// Map the masterVar column header (the `>`-prefixed line) to field indices by name. Returns
    /// `None` if a required column is missing (so a look-alike table isn't parsed as masterVar).
    fn from_header(line: &str) -> Option<Columns> {
        let header = line.strip_prefix('>').unwrap_or(line);
        let names: Vec<&str> = header.split('\t').map(str::trim).collect();
        let idx = |want: &str| names.iter().position(|n| n.eq_ignore_ascii_case(want));
        Some(Columns {
            locus: idx("locus")?,
            ploidy: idx("ploidy")?,
            allele: idx("allele")?,
            chromosome: idx("chromosome")?,
            begin: idx("begin")?,
            var_type: idx("varType")?,
            reference: idx("reference")?,
            allele_seq: idx("alleleSeq")?,
            x_ref: idx("xRef"),
        })
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Allele {
    All,
    One,
    Two,
    Other,
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum VarType {
    Snp,
    Ref,
    Other,
}

/// One parsed data row (only the fields we need). Borrows nothing — a locus group is buffered
/// before it is resolved, so the owning line has already been consumed.
struct Row {
    locus: u64,
    ploidy: u8,
    allele: Allele,
    var_type: VarType,
    chromosome: String,
    begin: i64,
    reference: String,
    allele_seq: String,
    rs_id: Option<String>,
}

impl Row {
    fn parse(fields: &[&str], c: &Columns) -> Option<Row> {
        let get = |i: usize| fields.get(i).copied().unwrap_or("");
        let locus = get(c.locus).parse::<u64>().ok()?;
        let begin = get(c.begin).parse::<i64>().ok()?;
        let allele = match get(c.allele) {
            "all" => Allele::All,
            "1" => Allele::One,
            "2" => Allele::Two,
            _ => Allele::Other,
        };
        let var_type = match get(c.var_type) {
            "snp" => VarType::Snp,
            "ref" => VarType::Ref,
            _ => VarType::Other,
        };
        Some(Row {
            locus,
            ploidy: get(c.ploidy).parse::<u8>().unwrap_or(0),
            allele,
            var_type,
            chromosome: get(c.chromosome).to_string(),
            begin,
            reference: get(c.reference).to_string(),
            allele_seq: get(c.allele_seq).to_string(),
            rs_id: c.x_ref.and_then(|i| first_rs_id(get(i))),
        })
    }
}

/// Pull the first `rs<digits>` accession out of a CG `xRef` cell like
/// `dbsnp.100:rs2748067;dbsnp.131:rs76046194`. Returns `None` when the cell has no rsID.
fn first_rs_id(x_ref: &str) -> Option<String> {
    let bytes = x_ref.as_bytes();
    let mut i = 0;
    while i + 2 < bytes.len() {
        // Match "rs" not preceded by another letter/digit (so it starts an accession token).
        let starts = (bytes[i] == b'r' || bytes[i] == b'R')
            && (bytes[i + 1] == b's' || bytes[i + 1] == b'S')
            && bytes[i + 2].is_ascii_digit();
        if starts {
            let mut end = i + 2;
            while end < bytes.len() && bytes[end].is_ascii_digit() {
                end += 1;
            }
            return Some(x_ref[i..end].to_string());
        }
        i += 1;
    }
    None
}

/// The reconstructed base on one haplotype of a diploid locus.
enum Hap {
    Ref,
    Alt(String),
    Missing,
}

/// Resolve one locus's rows into at most one biallelic SNP [`VariantCall`]. Returns `None` for
/// a locus with no `snp` row (a `ref`/`no-call`/indel span) or one whose SNP alleles aren't a
/// clean single-base substitution.
fn locus_call(rows: &[Row]) -> Option<VariantCall> {
    // Site anchor: the first snp row carries the reference base + coordinates.
    let snp = rows.iter().find(|r| r.var_type == VarType::Snp)?;
    let contig = snp.chromosome.as_str();
    let position = snp.begin + 1; // 0-based begin → 1-based SNP position
    let reference = snp.reference.as_str();

    // A dbSNP rsID from any of the locus's rows.
    let rs_id = rows.iter().find_map(|r| r.rs_id.clone());

    // Haploid contigs (chrY, chrM, male non-PAR chrX) come through as ploidy 1 with a single
    // called allele — emit the alt with a hemizygous genotype.
    if snp.ploidy == 1 {
        return snp_call(contig, position, reference, &snp.allele_seq, rs_id, Some("1".into()));
    }

    // Diploid: reconstruct the two haplotypes. An `all` snp row means both alleles carry it.
    if let Some(all) = rows.iter().find(|r| r.allele == Allele::All && r.var_type == VarType::Snp) {
        return snp_call(contig, position, reference, &all.allele_seq, rs_id, Some("1/1".into()));
    }
    let hap = |which: Allele| -> Hap {
        // A compound locus can list several rows for one allele (e.g. a `ref` segment beside the
        // `snp` segment). Prefer the SNP call for that allele; fall back to ref, else missing —
        // so a snp isn't hidden behind a same-allele ref/no-call row and the locus lost.
        let mut result = Hap::Missing;
        for r in rows.iter().filter(|r| r.allele == which) {
            match r.var_type {
                VarType::Snp => return Hap::Alt(r.allele_seq.clone()),
                VarType::Ref => result = Hap::Ref,
                VarType::Other => {}
            }
        }
        result
    };
    let (h1, h2) = (hap(Allele::One), hap(Allele::Two));

    let (alt, genotype) = match (&h1, &h2) {
        (Hap::Alt(a), Hap::Alt(b)) if a == b => (a.clone(), "1/1"),
        (Hap::Alt(a), Hap::Alt(_)) => (a.clone(), "1/2"), // tri-allelic het; keep allele 1's alt
        (Hap::Alt(a), Hap::Ref) | (Hap::Ref, Hap::Alt(a)) => (a.clone(), "0/1"),
        (Hap::Alt(a), Hap::Missing) | (Hap::Missing, Hap::Alt(a)) => (a.clone(), "1/."),
        // No alt on either haplotype — not a variant (shouldn't occur given the snp row above).
        _ => return None,
    };
    snp_call(contig, position, reference, &alt, rs_id, Some(genotype.into()))
}

/// Reference build for a `#GENOME_REFERENCE` header value. CG only ever shipped NCBI build 37.
fn build_for_reference(genome_reference: &str) -> String {
    let g = genome_reference.to_ascii_lowercase();
    if g.contains("build 37") || g.contains("grch37") || g.contains("b37") || g.contains("hg19") {
        "GRCh37".to_string()
    } else if g.contains("build 38") || g.contains("grch38") || g.contains("hg38") {
        "GRCh38".to_string()
    } else {
        // Unknown/blank header: CG service data is build 37, so default there rather than to None.
        "GRCh37".to_string()
    }
}

/// Read and parse a masterVar file (transparently decompressing `.bz2` / `.gz`).
pub fn parse_file(path: &Path) -> Result<MasterVarImport, MasterVarError> {
    let reader = gzio::open_maybe_compressed(path)?;
    parse_reader(reader)
}

/// Parse masterVar rows from any buffered reader (the streaming core; testable without IO).
pub fn parse_reader(reader: impl BufRead) -> Result<MasterVarImport, MasterVarError> {
    let mut out = MasterVarImport {
        reference_build: "GRCh37".to_string(),
        ..Default::default()
    };
    let mut columns: Option<Columns> = None;
    let mut group: Vec<Row> = Vec::new();
    let mut current: Option<u64> = None;

    let flush = |group: &mut Vec<Row>, out: &mut MasterVarImport| {
        if group.is_empty() {
            return;
        }
        out.loci_seen += 1;
        if let Some(call) = locus_call(group) {
            out.snp_loci += 1;
            out.calls.push(call);
        }
        group.clear();
    };

    for line in reader.lines() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        if let Some(meta) = line.strip_prefix('#') {
            // Header key/value pairs are tab-separated: `#SAMPLE\tGS00253-DNA_A01`.
            if let Some((key, val)) = meta.split_once('\t') {
                match key.trim() {
                    "SAMPLE" => out.sample_id = Some(val.trim().to_string()),
                    "GENOME_REFERENCE" => out.reference_build = build_for_reference(val.trim()),
                    _ => {}
                }
            }
            continue;
        }
        if line.starts_with('>') {
            columns = Some(Columns::from_header(&line).ok_or_else(|| {
                MasterVarError::Format("column header is missing required masterVar fields".into())
            })?);
            continue;
        }
        let Some(c) = columns.as_ref() else {
            return Err(MasterVarError::Format(
                "data row before the `>locus …` column header".into(),
            ));
        };
        let fields: Vec<&str> = line.split('\t').collect();
        let Some(row) = Row::parse(&fields, c) else {
            continue; // unparseable row (bad locus/begin) — skip
        };
        match current {
            Some(l) if l == row.locus => group.push(row),
            _ => {
                flush(&mut group, &mut out);
                current = Some(row.locus);
                group.push(row);
            }
        }
    }
    flush(&mut group, &mut out);

    if columns.is_none() {
        return Err(MasterVarError::Format(
            "no `>locus …` column header found (not a masterVar file)".into(),
        ));
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    const HEADER: &str = "#ASSEMBLY_ID\tGS00253-DNA_A01_200_37-ASM\n\
         #GENOME_REFERENCE\tNCBI build 37\n\
         #SAMPLE\tGS00253-DNA_A01\n\
         #GENERATED_BY\tcgatools\n\
         #TYPE\tVAR-ANNOTATION\n\
         >locus\tploidy\tallele\tchromosome\tbegin\tend\tvarType\treference\talleleSeq\tvarScoreVAF\tvarScoreEAF\tvarQuality\thapLink\txRef\n";

    fn parse(body: &str) -> MasterVarImport {
        parse_reader(format!("{HEADER}{body}").as_bytes()).unwrap()
    }

    #[test]
    fn parses_header_metadata() {
        let out = parse("1\t2\tall\tchr1\t0\t10000\tno-ref\t=\t?\t\t\t\t\t\n");
        assert_eq!(out.sample_id.as_deref(), Some("GS00253-DNA_A01"));
        assert_eq!(out.reference_build, "GRCh37");
        assert!(out.calls.is_empty()); // no-ref span yields no SNP
    }

    #[test]
    fn homozygous_snp_two_allele_rows() {
        // Both haplotypes T over ref C → 1/1.
        let out = parse(
            "272\t2\t1\tchr1\t21579\t21580\tsnp\tC\tT\t100\t100\tVQHIGH\t\tdbsnp.83:rs526642\n\
             272\t2\t2\tchr1\t21579\t21580\tsnp\tC\tT\t135\t135\tVQHIGH\t\tdbsnp.83:rs526642\n",
        );
        assert_eq!(out.calls.len(), 1);
        let c = &out.calls[0];
        assert_eq!((c.contig.as_str(), c.position), ("chr1", 21580));
        assert_eq!((c.reference.as_str(), c.alternate.as_str()), ("C", "T"));
        assert_eq!(c.genotype.as_deref(), Some("1/1"));
        assert_eq!(c.rs_id.as_deref(), Some("rs526642"));
    }

    #[test]
    fn heterozygous_snp_ref_plus_snp() {
        // allele 1 ref, allele 2 snp A→G → 0/1.
        let out = parse(
            "896\t2\t1\tchr1\t46669\t46670\tref\tA\tA\t45\t45\tVQHIGH\t\t\n\
             896\t2\t2\tchr1\t46669\t46670\tsnp\tA\tG\t45\t45\tVQHIGH\t\tdbsnp.100:rs2548905\n",
        );
        assert_eq!(out.calls.len(), 1);
        let c = &out.calls[0];
        assert_eq!(c.position, 46670);
        assert_eq!((c.reference.as_str(), c.alternate.as_str()), ("A", "G"));
        assert_eq!(c.genotype.as_deref(), Some("0/1"));
        assert_eq!(c.rs_id.as_deref(), Some("rs2548905"));
    }

    #[test]
    fn heterozygous_snp_with_no_call_partner() {
        // allele 1 snp G→A, allele 2 no-call → 1/. (one alt observed, partner unknown).
        let out = parse(
            "344\t2\t1\tchr1\t23974\t23975\tsnp\tG\tA\t58\t58\tVQHIGH\t\tdbsnp.100:rs2748067\n\
             344\t2\t2\tchr1\t23974\t23975\tno-call\tG\t?\t\t\t\t\t\n",
        );
        assert_eq!(out.calls.len(), 1);
        let c = &out.calls[0];
        assert_eq!((c.reference.as_str(), c.alternate.as_str()), ("G", "A"));
        assert_eq!(c.genotype.as_deref(), Some("1/."));
    }

    #[test]
    fn haploid_y_and_m_snps() {
        // chrY / chrM come through as ploidy 1, allele 1 → hemizygous "1".
        let out = parse(
            "21316594\t1\t1\tchrY\t2661693\t2661694\tsnp\tA\tG\t342\t342\tVQHIGH\t\tdbsnp.100:rs2253109\n\
             21394470\t1\t1\tchrM\t182\t183\tsnp\tA\tG\t5431\t5431\tVQHIGH\t\tdbsnp.132:rs113913230\n",
        );
        assert_eq!(out.calls.len(), 2);
        assert_eq!(out.calls[0].contig, "chrY");
        assert_eq!(out.calls[0].position, 2661694);
        assert_eq!(out.calls[0].genotype.as_deref(), Some("1"));
        assert_eq!(out.calls[1].contig, "chrM");
        assert_eq!(out.calls[1].position, 183);
    }

    #[test]
    fn compound_locus_snp_hidden_behind_same_allele_ref() {
        // A single locus can list several segments for one allele: here allele 1 has a `ref`
        // segment *before* its `snp` segment. The snp must still win (not be masked by the ref),
        // so the locus yields a 0/1 call rather than being dropped.
        let out = parse(
            "500\t2\t1\tchr1\t9000\t9005\tref\tACGTA\tACGTA\t\t\t\t\t\n\
             500\t2\t1\tchr1\t9005\t9006\tsnp\tT\tC\t80\t80\tVQHIGH\t\tdbsnp.1:rs99\n\
             500\t2\t2\tchr1\t9000\t9006\tref\t=\t=\t\t\t\t\t\n",
        );
        assert_eq!(out.calls.len(), 1);
        let c = &out.calls[0];
        assert_eq!((c.position, c.reference.as_str(), c.alternate.as_str()), (9006, "T", "C"));
        assert_eq!(c.genotype.as_deref(), Some("0/1"));
    }

    #[test]
    fn ref_and_indel_spans_are_skipped() {
        let out = parse(
            "3\t2\tall\tchr1\t11099\t11109\tref\t=\t=\t\t\t\t\t\n\
             100\t2\t1\tchr1\t1000\t1002\tdel\tAT\t\t\t\t\t\t\n\
             200\t2\t1\tchr1\t2000\t2001\tins\t\tCG\t\t\t\t\t\n\
             300\t2\t2\tchr1\t3000\t3005\tsub\tACGTA\tGGGGG\t\t\t\t\t\n",
        );
        assert!(out.calls.is_empty());
        assert_eq!(out.loci_seen, 4);
        assert_eq!(out.snp_loci, 0);
    }

    #[test]
    fn rejects_non_mastervar_input() {
        // No `>locus …` header → not a masterVar file.
        let err = parse_reader("rsid\tchrom\tpos\nrs1\t1\t100\n".as_bytes());
        assert!(matches!(err, Err(MasterVarError::Format(_))));
    }

    #[test]
    fn first_rs_id_extracts_first_accession() {
        assert_eq!(first_rs_id("dbsnp.100:rs2748067;dbsnp.131:rs76046194").as_deref(), Some("rs2748067"));
        assert_eq!(first_rs_id("").as_deref(), None);
        assert_eq!(first_rs_id("cosmic:COSM123").as_deref(), None);
    }
}
