//! The reader of an EIGENSTRAT call set, which comes from the Reich lab and from `pileupCaller`.
//! It is the path for an external autosomal 1240K set.
//!
//! It parses a `.geno`, `.snp` and `.ind` triplet, for **one** target individual, into diploid
//! allele pairs on the build of the call set. The AADR 1240K is on GRCh37, which is hg19.
//!
//! Those pairs go to [`crate::ibd_panel::IbdPanel::resolve_chip`]. That function keys them to
//! canonical CHM13 dosages, and it **orients them itself** against the CHM13 alleles. So the allele
//! labels of the EIGENSTRAT set do not have to match the genome reference. And a pseudo-haploid
//! call, where pileupCaller writes `0` or `2` alone, comes through as a correct homozygous
//! observation, with no het that the code invented.
//!
//! The EIGENSTRAT format is text, with white space between its fields:
//!
//! - `.ind` holds one line for each individual: `SampleID  Sex  Population`.
//! - `.snp` holds one line for each SNP:
//!   `SNPName  Chr  GeneticPos  PhysicalPos  RefAllele  VariantAllele`.
//! - `.geno` holds one line for each SNP, in the same order as `.snp`, with one character for each
//!   individual. A `0`, `1` or `2` is the **count of the first allele**, which is column 5 of
//!   `.snp`. A `9` means that the value is missing.

use std::io::BufRead;
use std::path::Path;

use crate::error::AnalysisError;

/// The genotypes of one target individual, from an EIGENSTRAT triplet. They are diploid allele
/// pairs on `build`, forward on the reference, ready for `IbdPanel::resolve_chip`. The code drops a
/// no-call, which is a `9`, and it emits nothing for it.
pub struct CallSet {
    /// The build that the `.snp` positions sit on. EIGENSTRAT does not record it, so the caller
    /// gives it. The default is GRCh37, which is the coordinate system of the AADR 1240K.
    pub build: String,
    /// A `(contig, position, allele1, allele2)` at each autosomal site with a call. The contig is
    /// bare, from `"1"` to `"22"`, which matches the GRCh37 loci of the panel. The alleles are the
    /// nucleotides that the individual carries.
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

/// The bare autosomal contig, from `"1"` to `"22"`, for a `Chr` field of EIGENSTRAT. It is `None`
/// for a sex contig, for mt, and for a contig that the code does not know. EIGENSTRAT writes `23`
/// for X, `24` for Y, `90` or `91` for mt, and it also uses `0`. This accepts a `chr` prefix, and
/// it does not need one.
fn autosome_contig(chr: &str) -> Option<String> {
    match crate::contig::bare(chr).parse::<u8>() {
        Ok(n @ 1..=22) => Some(n.to_string()),
        _ => None,
    }
}

/// Select the target individual's column index in the `.ind` file (0-based, matching the `.geno`
/// character position). `sample` names it; a single-individual file needs no name.
fn select_individual(ind_text: &str, sample: Option<&str>) -> Result<usize, AnalysisError> {
    let ids: Vec<&str> = ind_text.lines().filter_map(|l| l.split_whitespace().next()).collect();
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

/// The streaming core. It walks `.snp` and `.geno` in step, and it emits the autosomal calls of the
/// target column. It sits apart, over a `BufRead`, so that a test can cover it with no file.
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

/// Read an EIGENSTRAT triplet, for one target individual. `sample` selects that individual, and it
/// is necessary when the `.ind` lists more than one. `build` is the coordinate system of the `.snp`
/// positions, and the default is GRCh37, for the AADR 1240K. It streams `.snp` and `.geno`, so a
/// large panel costs little.
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
        // In each row, character 0 is SAMPLE_A and character 1 is SAMPLE_B.
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
            vec![("1".to_string(), 1000, 'G', 'G'), ("22".to_string(), 4000, 'T', 'C'),]
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
