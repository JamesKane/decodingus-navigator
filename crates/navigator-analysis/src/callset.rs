//! EIGENSTRAT (Reich-lab / `pileupCaller`) call-set reader — the external autosomal 1240K path.
//!
//! Parses a `.geno`/`.snp`/`.ind` triplet for **one** target individual into diploid allele pairs on
//! the call set's build (the AADR 1240K is GRCh37/hg19). Those pairs feed
//! [`crate::ibd_panel::IbdPanel::resolve_chip`], which re-keys them to canonical CHM13 dosages and
//! **self-orients** against the CHM13 alleles — so the EIGENSTRAT allele labelling need not match the
//! genome reference, and pseudo-haploid calls (pileupCaller emits only `0`/`2`) come through as valid
//! homozygous observations with no het synthesis.
//!
//! EIGENSTRAT format (whitespace-delimited text):
//! - `.ind`: one line per individual — `SampleID  Sex  Population`.
//! - `.snp`: one line per SNP — `SNPName  Chr  GeneticPos  PhysicalPos  RefAllele  VariantAllele`.
//! - `.geno`: one line per SNP (same order as `.snp`), one char per individual —
//!   `0`/`1`/`2` = **count of the first (`.snp` column-5) allele**, `9` = missing.

use std::io::BufRead;
use std::path::Path;

use crate::error::AnalysisError;

/// One target individual's genotypes from an EIGENSTRAT triplet: reference-forward diploid allele
/// pairs on `build`, ready for `IbdPanel::resolve_chip`. No-calls (`9`) are dropped, not emitted.
pub struct CallSet {
    /// Build the `.snp` positions are in. EIGENSTRAT does not encode it; the caller supplies it
    /// (default GRCh37 — the AADR 1240K coordinate system).
    pub build: String,
    /// `(contig, position, allele1, allele2)` per called autosomal site (contig is a bare `"1".."22"`,
    /// matching the panel's GRCh37 loci). Alleles are the actual nucleotides the individual carries.
    pub calls: Vec<(String, i64, char, char)>,
    /// Count of autosomal `.snp` sites that were missing (`9`) for this individual.
    pub missing: usize,
}

/// Map a `.geno` value (count of the first `.snp` allele) to the diploid allele pair, or `None` for
/// missing (`9` / anything else). `2` = two copies of `a1`; `1` = one of each; `0` = two of `a2`.
/// pileupCaller pseudo-haploid emits only `0`/`2`, which land as valid homozygous pairs here.
fn geno_to_pair(g: u8, a1: char, a2: char) -> Option<(char, char)> {
    match g {
        b'2' => Some((a1, a1)),
        b'1' => Some((a1, a2)),
        b'0' => Some((a2, a2)),
        _ => None,
    }
}

/// Bare autosomal contig (`"1".."22"`) for an EIGENSTRAT `Chr` field, or `None` for a sex/mt/unknown
/// contig (EIGENSTRAT uses `23`=X, `24`=Y, `90`/`91`=mt, plus `0`). Accepts an optional `chr` prefix.
fn autosome_contig(chr: &str) -> Option<String> {
    let bare = chr.strip_prefix("chr").or_else(|| chr.strip_prefix("Chr")).unwrap_or(chr);
    match bare.parse::<u8>() {
        Ok(n @ 1..=22) => Some(n.to_string()),
        _ => None,
    }
}

/// Select the target individual's column index in the `.ind` file (0-based, matching the `.geno`
/// character position). `sample` names it; a single-individual file needs no name.
fn select_individual(ind_text: &str, sample: Option<&str>) -> Result<usize, AnalysisError> {
    let ids: Vec<&str> = ind_text
        .lines()
        .filter_map(|l| l.split_whitespace().next())
        .collect();
    if ids.is_empty() {
        return Err(AnalysisError::Message("EIGENSTRAT .ind has no individuals".into()));
    }
    match sample {
        Some(s) => ids
            .iter()
            .position(|id| *id == s)
            .ok_or_else(|| AnalysisError::Message(format!("EIGENSTRAT .ind has no individual {s:?}"))),
        None if ids.len() == 1 => Ok(0),
        None => Err(AnalysisError::Message(format!(
            "EIGENSTRAT .ind has {} individuals — specify which to import (one of: {})",
            ids.len(),
            ids.join(", ")
        ))),
    }
}

/// Streamed core: walk `.snp` and `.geno` in lockstep, emitting the target column's autosomal calls.
/// Split out (over `BufRead`) so it is testable without files.
fn read_eigenstrat_core<S: BufRead, G: BufRead>(
    snp: S,
    geno: G,
    col: usize,
    build: &str,
) -> Result<CallSet, AnalysisError> {
    let mut calls = Vec::new();
    let mut missing = 0usize;
    let mut snp_lines = snp.lines();
    let mut geno_lines = geno.lines();
    loop {
        match (snp_lines.next(), geno_lines.next()) {
            (Some(s), Some(g)) => {
                let s = s.map_err(|e| AnalysisError::Message(format!("reading .snp: {e}")))?;
                let g = g.map_err(|e| AnalysisError::Message(format!("reading .geno: {e}")))?;
                let f: Vec<&str> = s.split_whitespace().collect();
                if f.len() < 6 {
                    continue; // blank / short line
                }
                let Some(contig) = autosome_contig(f[1]) else { continue };
                let Ok(pos) = f[3].parse::<i64>() else { continue };
                let a1 = f[4].chars().next().unwrap_or('N').to_ascii_uppercase();
                let a2 = f[5].chars().next().unwrap_or('N').to_ascii_uppercase();
                let val = g.as_bytes().get(col).copied().unwrap_or(b'9');
                match geno_to_pair(val, a1, a2) {
                    Some((b1, b2)) => calls.push((contig, pos, b1, b2)),
                    None => missing += 1,
                }
            }
            (None, None) => break,
            _ => {
                return Err(AnalysisError::Message(
                    "EIGENSTRAT .snp and .geno have different row counts".into(),
                ))
            }
        }
    }
    Ok(CallSet {
        build: build.to_string(),
        calls,
        missing,
    })
}

/// Read an EIGENSTRAT triplet for one target individual. `sample` selects the individual (required
/// when the `.ind` lists more than one); `build` is the coordinate system of the `.snp` positions
/// (default GRCh37 for the AADR 1240K). Streams `.snp`/`.geno` so a large panel stays cheap.
pub fn read_eigenstrat(
    geno: &Path,
    snp: &Path,
    ind: &Path,
    sample: Option<&str>,
    build: &str,
) -> Result<CallSet, AnalysisError> {
    let ind_text = std::fs::read_to_string(ind).map_err(|e| AnalysisError::io(ind, e))?;
    let col = select_individual(&ind_text, sample)?;
    let snp_rd = std::io::BufReader::new(std::fs::File::open(snp).map_err(|e| AnalysisError::io(snp, e))?);
    let geno_rd = std::io::BufReader::new(std::fs::File::open(geno).map_err(|e| AnalysisError::io(geno, e))?);
    read_eigenstrat_core(snp_rd, geno_rd, col, build)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn geno_value_maps_to_the_right_pair() {
        // geno counts the FIRST (.snp col-5) allele: 2 → hom-a1, 0 → hom-a2, 1 → het.
        assert_eq!(geno_to_pair(b'2', 'A', 'G'), Some(('A', 'A')));
        assert_eq!(geno_to_pair(b'1', 'A', 'G'), Some(('A', 'G')));
        assert_eq!(geno_to_pair(b'0', 'A', 'G'), Some(('G', 'G')));
        assert_eq!(geno_to_pair(b'9', 'A', 'G'), None); // missing
    }

    #[test]
    fn autosomes_only() {
        assert_eq!(autosome_contig("1"), Some("1".to_string()));
        assert_eq!(autosome_contig("chr22"), Some("22".to_string()));
        assert_eq!(autosome_contig("23"), None); // X
        assert_eq!(autosome_contig("90"), None); // mt
        assert_eq!(autosome_contig("0"), None);
    }

    #[test]
    fn selects_the_named_individual() {
        let ind = "SAMPLE_A M PopA\nSAMPLE_B F PopB\n";
        assert_eq!(select_individual(ind, Some("SAMPLE_B")).unwrap(), 1);
        assert!(select_individual(ind, None).is_err()); // ambiguous — must name one
        assert!(select_individual(ind, Some("nope")).is_err());
        assert_eq!(select_individual("ONLY M Pop\n", None).unwrap(), 0); // single → no name needed
    }

    #[test]
    fn reads_the_target_column_and_skips_non_autosomes_and_missing() {
        // Two individuals; import column 1 (SAMPLE_B).
        let snp = "\
rs1 1 0.0 1000 A G
rs2 2 0.0 2000 C T
rsX 23 0.0 3000 A G
rs3 22 0.0 4000 T C
";
        // per row, char 0 = SAMPLE_A, char 1 = SAMPLE_B.
        let geno = "\
20
19
02
21
";
        let cs = read_eigenstrat_core(Cursor::new(snp), Cursor::new(geno), 1, "GRCh37").unwrap();
        // SAMPLE_B: rs1=0→(G,G); rs2=9→missing; rsX skipped (chr23); rs3=1→(T,C).
        assert_eq!(
            cs.calls,
            vec![
                ("1".to_string(), 1000, 'G', 'G'),
                ("22".to_string(), 4000, 'T', 'C'),
            ]
        );
        assert_eq!(cs.missing, 1);
        assert_eq!(cs.build, "GRCh37");
    }

    #[test]
    fn row_count_mismatch_errors() {
        let snp = "rs1 1 0.0 1000 A G\nrs2 1 0.0 2000 C T\n";
        let geno = "2\n"; // one row short
        assert!(read_eigenstrat_core(Cursor::new(snp), Cursor::new(geno), 0, "GRCh37").is_err());
    }
}
