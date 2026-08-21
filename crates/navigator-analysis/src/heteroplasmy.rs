//! mtDNA heteroplasmy detection (plan §4b reconciliation, phase 6).
//!
//! The haploid consensus caller brings each position down to one base. Heteroplasmy is the
//! opposite: two mitochondrial alleles that live together in one individual.
//!
//! The code finds it with a scan over the A/C/G/T pileup at every chrM position. It flags a site
//! where a second allele sits above a noise floor. That means a minor-allele fraction inside
//! `[min_minor_fraction, 0.5]`, with `min_minor_count` reads behind it or more.
//!
//! This pass looks for candidates. It is not a clinical caller. It reports the observed allele
//! fractions, so that a curator can judge real heteroplasmy against an artifact of the sequencing.
//! Contamination from a NUMT, strand bias and homopolymer noise are such artifacts.
//!
//! chrM is about 16.5 kb, so the code tallies the whole contig in one dense pass, through
//! `tally_region` in the caller.

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::caller::{self, HaploidCallerParams};
use crate::error::AnalysisError;

/// A position that holds two alleles above the noise floor.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HeteroplasmySite {
    /// 1-based position on the contig.
    pub position: i64,
    /// The depth at this position that passes, which is the count of reads that clear the quality
    /// filters.
    pub depth: u32,
    /// The dominant base.
    pub major_base: char,
    /// The count of reads behind the major base.
    pub major_count: u32,
    /// The second-most-common base.
    pub minor_base: char,
    /// The count of reads behind the minor base.
    pub minor_count: u32,
    /// `minor_count / depth`. It is the level of the heteroplasmy.
    pub minor_fraction: f64,
}

/// The thresholds that a site must meet before the code calls it heteroplasmic. The defaults are
/// careful values, for a search over candidates. The parity harness of §4c gates them, as it gates
/// the rest of the caller.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HeteroplasmyParams {
    /// The code skips a site whose depth, over the reads that pass, is below this. Such a site is
    /// too thin to judge.
    pub min_depth: u32,
    /// Minor-allele fraction must reach this to flag (the noise floor).
    pub min_minor_fraction: f64,
    /// The minor allele needs this many reads behind it, or more.
    pub min_minor_count: u32,
    /// The code drops a read below this MAPQ.
    pub min_mapping_quality: u8,
    /// Bases below this quality are not counted.
    pub min_base_quality: u8,
}

impl Default for HeteroplasmyParams {
    fn default() -> Self {
        // A noise floor of 3%, with 3 reads behind the allele or more, is a common default for a
        // search over mtDNA heteroplasmy candidates on short-read data. A depth of 20 keeps the
        // fractions meaningful.
        HeteroplasmyParams {
            min_depth: 20,
            min_minor_fraction: 0.03,
            min_minor_count: 3,
            min_mapping_quality: 20,
            min_base_quality: 20,
        }
    }
}

const BASES: [char; 4] = ['A', 'C', 'G', 'T'];

/// Top two `(base_index, count)` by count; ties keep the earlier base (A<C<G<T).
fn top_two(counts: &[u32; 4]) -> ((usize, u32), (usize, u32)) {
    let mut first = (0usize, counts[0]);
    let mut second = (0usize, 0u32);
    for (i, &c) in counts.iter().enumerate() {
        if c > first.1 {
            second = first;
            first = (i, c);
        } else if i != first.0 && c > second.1 {
            second = (i, c);
        }
    }
    (first, second)
}

/// Scan every position on `contig`, and return the heteroplasmic sites, from the lowest position
/// up. It tallies the whole contig in one pass, which is correct for a contig the size of chrM.
pub fn detect_heteroplasmy(
    bam_path: &Path,
    contig: &str,
    params: &HeteroplasmyParams,
    reference: Option<&Path>,
) -> Result<Vec<HeteroplasmySite>, AnalysisError> {
    let length = caller::read_contig_length(bam_path, contig, reference)?;
    // Use the pileup of the caller again, with the same quality gates. The gates here, on the
    // allele fraction and the minimum depth, belong to heteroplasmy alone. So the parameters of
    // the caller stay open.
    let caller_params = HaploidCallerParams {
        min_depth: 1,
        min_mapping_quality: params.min_mapping_quality,
        min_base_quality: params.min_base_quality,
        min_allele_fraction: 0.0,
        ..HaploidCallerParams::default()
    };
    let (counts, _indel) = caller::tally_region(bam_path, contig, &caller_params, 1, length, reference)?;

    let mut sites = Vec::new();
    for (offset, c) in counts.iter().enumerate() {
        let depth: u32 = c.iter().sum();
        if depth < params.min_depth {
            continue;
        }
        let ((maj_i, maj_n), (min_i, min_n)) = top_two(c);
        if min_n < params.min_minor_count {
            continue;
        }
        let minor_fraction = min_n as f64 / depth as f64;
        if minor_fraction < params.min_minor_fraction {
            continue;
        }
        sites.push(HeteroplasmySite {
            position: (offset + 1) as i64,
            depth,
            major_base: BASES[maj_i],
            major_count: maj_n,
            minor_base: BASES[min_i],
            minor_count: min_n,
            minor_fraction,
        });
    }
    Ok(sites)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn top_two_picks_two_highest() {
        // A=10, C=3, G=1, T=0
        let ((mi, mn), (ni, nn)) = top_two(&[10, 3, 1, 0]);
        assert_eq!((mi, mn), (0, 10));
        assert_eq!((ni, nn), (1, 3));
    }

    #[test]
    fn top_two_ties_keep_earlier_base() {
        // A=5, C=5 → major A, minor C
        let ((mi, _), (ni, _)) = top_two(&[5, 5, 0, 0]);
        assert_eq!(mi, 0);
        assert_eq!(ni, 1);
    }

    #[test]
    fn top_two_single_allele_minor_zero() {
        let (_, (_, nn)) = top_two(&[30, 0, 0, 0]);
        assert_eq!(nn, 0);
    }
}
