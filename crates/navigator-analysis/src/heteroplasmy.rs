//! mtDNA heteroplasmy detection (plan §4b reconciliation, phase 6).
//!
//! Unlike the haploid consensus caller — which collapses each position to a single
//! base — heteroplasmy is the *coexistence* of two mitochondrial alleles in one
//! individual. We detect it by scanning every chrM position's A/C/G/T pileup and
//! flagging sites where a second allele is present above a noise floor: a minor-allele
//! fraction in `[min_minor_fraction, 0.5]` backed by at least `min_minor_count` reads.
//!
//! This is a screening pass, not a clinical caller: it reports observed allele
//! fractions so a curator can judge real heteroplasmy versus sequencing artefacts
//! (NUMT contamination, strand bias, homopolymer noise). chrM is ~16.5 kb, so the
//! whole contig is tallied in a single dense pass via the caller's `tally_region`.

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::caller::{self, HaploidCallerParams};
use crate::error::AnalysisError;

/// A position carrying two alleles above the noise floor.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct HeteroplasmySite {
    /// 1-based position on the contig.
    pub position: i64,
    /// Passing depth (reads clearing the quality filters) at this position.
    pub depth: u32,
    /// The dominant base.
    pub major_base: char,
    /// Reads supporting the major base.
    pub major_count: u32,
    /// The second-most-common base.
    pub minor_base: char,
    /// Reads supporting the minor base.
    pub minor_count: u32,
    /// `minor_count / depth` — the heteroplasmy level.
    pub minor_fraction: f64,
}

/// Thresholds for calling a site heteroplasmic. Defaults are conservative screening
/// values (gated by the §4c parity harness, like the rest of the caller).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HeteroplasmyParams {
    /// Sites below this passing depth are skipped (too shallow to judge).
    pub min_depth: u32,
    /// Minor-allele fraction must reach this to flag (the noise floor).
    pub min_minor_fraction: f64,
    /// Minor allele must be backed by at least this many reads.
    pub min_minor_count: u32,
    /// Reads below this MAPQ are dropped.
    pub min_mapping_quality: u8,
    /// Bases below this quality are not counted.
    pub min_base_quality: u8,
}

impl Default for HeteroplasmyParams {
    fn default() -> Self {
        // 3% noise floor with ≥3 supporting reads is a common screening default for
        // mtDNA heteroplasmy on short-read data; depth 20 keeps fractions meaningful.
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

/// Scan every position on `contig` and return the heteroplasmic sites, ascending by
/// position. Tallies the full contig in one pass (fine for chrM-sized contigs).
pub fn detect_heteroplasmy(
    bam_path: &Path,
    contig: &str,
    params: &HeteroplasmyParams,
    reference: Option<&Path>,
) -> Result<Vec<HeteroplasmySite>, AnalysisError> {
    let length = caller::read_contig_length(bam_path, contig, reference)?;
    // Reuse the caller's pileup with matching quality gates; allele-fraction/min-depth
    // gating here is heteroplasmy-specific, so the caller params stay permissive.
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
