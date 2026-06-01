//! Purpose-built **haploid** variant caller (plan §4b) — the GATK replacement for
//! Y/mtDNA. There is no pure-code Scala caller to port: the legacy app shelled out to
//! GATK `HaplotypeCaller --sample-ploidy 1` (force-call at tree sites) and `Mutect2
//! --mitochondria` / haploid `HaplotypeCaller` (de-novo discovery), then subtracted
//! known tree positions to get private variants. This module reproduces both modes by
//! **pileup-consensus calling**, which is tractable precisely because Y and mtDNA are
//! haploid (ploidy 1) — no diploid local reassembly.
//!
//! Two modes:
//! 1. [`force_call_sites`] — genotype-given-alleles at known tree `Site`s (haplogroup
//!    assignment): pileup, take the consensus base, report whether it is the site's
//!    ref or alt allele.
//! 2. [`call_denovo`] — walk the contig, emit positions whose consensus base differs
//!    from the reference (the candidate private variants). [`subtract_known`] removes
//!    known tree positions to yield the private set.
//!
//! **v1 is SNP-only** (plan §4b): indels/homopolymers are where naive pileup calling
//! diverges from GATK (light local realignment is the planned mitigation), so indel
//! alleles are skipped here and treated as advisory until the §4c parity harness
//! validates them. Defaults are starting points the harness will tune.
//!
//! Memory: de-novo uses a dense per-position tally for the target contig (fine for
//! mtDNA/Y; whole-genome streaming is a later opt). Force-call tallies only the target
//! sites (sparse), so it is cheap regardless of contig size.

use std::collections::HashSet;
use std::path::Path;

use noodles::bam;
use noodles::core::Region;
use noodles::fasta;

use crate::error::AnalysisError;

/// Parameters for haploid calling. Defaults are v1 starting points (gated by §4c).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HaploidCallerParams {
    /// Minimum passing depth (reads clearing the quality filters) to make any call.
    pub min_depth: u32,
    /// Reads below this MAPQ are dropped entirely.
    pub min_mapping_quality: u8,
    /// Bases below this quality are not counted.
    pub min_base_quality: u8,
    /// The consensus base must be at least this fraction of passing depth to call.
    pub min_allele_fraction: f64,
}

impl Default for HaploidCallerParams {
    fn default() -> Self {
        HaploidCallerParams {
            min_depth: 4,
            min_mapping_quality: 20,
            min_base_quality: 20,
            min_allele_fraction: 0.5,
        }
    }
}

/// A known tree/ancestry site to genotype (mirrors the Scala `Locus`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Site {
    pub name: String,
    pub contig: String,
    pub position: i64, // 1-based
    pub reference_allele: String,
    pub alternate_allele: String,
}

/// The allele called at a force-call site (haploid → one allele).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CalledAllele {
    Reference,
    Alternate,
    /// Insufficient depth, below-threshold consensus, or consensus is a third allele.
    NoCall,
}

/// Genotype at a known site.
#[derive(Debug, Clone, PartialEq)]
pub struct GenotypeCall {
    pub name: String,
    pub contig: String,
    pub position: i64,
    pub reference_allele: String,
    pub alternate_allele: String,
    pub called: CalledAllele,
    pub depth: u32,     // passing depth (all bases)
    pub ref_depth: u32,
    pub alt_depth: u32,
    pub allele_fraction: f64, // alt_depth / depth
}

/// A de-novo SNP call (consensus base differs from reference).
#[derive(Debug, Clone, PartialEq)]
pub struct VariantCall {
    pub contig: String,
    pub position: i64, // 1-based
    pub reference_allele: char,
    pub alternate_allele: char,
    pub depth: u32,     // passing depth
    pub alt_depth: u32, // reads supporting the consensus alt
    pub allele_fraction: f64,
}

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

fn base_index(b: u8) -> Option<usize> {
    match b.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// argmax over the four base counts; ties keep the earlier base (A<C<G<T).
fn consensus(counts: &[u32; 4]) -> (usize, u32) {
    let mut bi = 0;
    let mut best = counts[0];
    for i in 1..4 {
        if counts[i] > best {
            best = counts[i];
            bi = i;
        }
    }
    (bi, best)
}

/// Per-position A/C/G/T tallies, dense (whole contig) or sparse (target sites only).
enum AlleleCounts {
    Dense(Vec<[u32; 4]>),
    Sparse(std::collections::HashMap<usize, [u32; 4]>),
}

impl AlleleCounts {
    fn add(&mut self, idx: usize, base: usize) {
        match self {
            AlleleCounts::Dense(v) => v[idx][base] += 1,
            AlleleCounts::Sparse(m) => m.entry(idx).or_insert([0; 4])[base] += 1,
        }
    }

    fn get(&self, idx: usize) -> [u32; 4] {
        match self {
            AlleleCounts::Dense(v) => v[idx],
            AlleleCounts::Sparse(m) => m.get(&idx).copied().unwrap_or([0; 4]),
        }
    }
}

/// Resolve a contig's length from the BAM header.
fn contig_length(header: &noodles::sam::Header, contig: &str) -> Option<usize> {
    header
        .reference_sequences()
        .iter()
        .find(|(name, _)| {
            let n: &[u8] = name.as_ref();
            n == contig.as_bytes()
        })
        .map(|(_, map)| map.length().get())
}

/// One indexed pass over `contig`, tallying passing A/C/G/T observations. When
/// `targets` is `Some`, only those 1-based positions are tallied (sparse).
fn tally_alleles(
    bam_path: &Path,
    contig: &str,
    params: &HaploidCallerParams,
    targets: Option<&HashSet<i64>>,
) -> Result<(usize, AlleleCounts), AnalysisError> {
    let mut reader = bam::io::indexed_reader::Builder::default()
        .build_from_path(bam_path)
        .map_err(|e| AnalysisError::io(bam_path, e))?;
    let header = reader
        .read_header()
        .map_err(|e| AnalysisError::io(bam_path, e))?;

    let length = contig_length(&header, contig)
        .ok_or_else(|| AnalysisError::Message(format!("contig {contig} not in BAM header")))?;

    let mut counts = match targets {
        Some(_) => AlleleCounts::Sparse(std::collections::HashMap::new()),
        None => AlleleCounts::Dense(vec![[0u32; 4]; length]),
    };

    let region: Region = contig
        .parse()
        .map_err(|_| AnalysisError::Message(format!("bad region for contig {contig}")))?;
    let query = reader
        .query(&header, &region)
        .map_err(|e| AnalysisError::io(bam_path, e))?;

    for result in query {
        let record = result.map_err(|e| AnalysisError::io(bam_path, e))?;
        let flags = record.flags();
        if flags.is_secondary() || flags.is_supplementary() || flags.is_duplicate() || flags.is_qc_fail() {
            continue;
        }
        let mapq = record.mapping_quality().map_or(255u8, |m| m.get());
        if mapq < params.min_mapping_quality {
            continue; // read-level MAPQ filter
        }
        let start = match record.alignment_start() {
            Some(p) => p.map_err(|e| AnalysisError::io(bam_path, e))?.get(),
            None => continue,
        };
        let seq = record.sequence();
        let quals = record.quality_scores();
        let quals = quals.as_ref();

        let mut ref_pos = start; // 1-based
        let mut query_off = 0usize;
        for op in record.cigar().iter() {
            let op = op.map_err(|e| AnalysisError::io(bam_path, e))?;
            let kind = op.kind();
            let len = op.len();
            match (kind.consumes_reference(), kind.consumes_read()) {
                (true, true) => {
                    for i in 0..len {
                        let pos = ref_pos + i; // 1-based
                        if pos >= 1 && pos <= length {
                            let base_q = quals.get(query_off + i).copied().unwrap_or(0);
                            if base_q >= params.min_base_quality {
                                if let Some(bi) = seq.get(query_off + i).and_then(base_index) {
                                    let tracked = targets.map_or(true, |t| t.contains(&(pos as i64)));
                                    if tracked {
                                        counts.add(pos - 1, bi);
                                    }
                                }
                            }
                        }
                    }
                    ref_pos += len;
                    query_off += len;
                }
                (true, false) => ref_pos += len,
                (false, true) => query_off += len,
                (false, false) => {}
            }
        }
    }

    Ok((length, counts))
}

/// Force-call (genotype-given-alleles) at known SNP sites on `contig`. Non-SNP sites
/// (multi-base ref/alt) are skipped — v1 is SNP-only.
pub fn force_call_sites(
    bam_path: &Path,
    contig: &str,
    sites: &[Site],
    params: &HaploidCallerParams,
) -> Result<Vec<GenotypeCall>, AnalysisError> {
    let targets: HashSet<i64> = sites
        .iter()
        .filter(|s| s.contig == contig)
        .map(|s| s.position)
        .collect();
    if targets.is_empty() {
        return Ok(Vec::new());
    }
    let (length, counts) = tally_alleles(bam_path, contig, params, Some(&targets))?;

    let mut out = Vec::new();
    for site in sites.iter().filter(|s| s.contig == contig) {
        if site.reference_allele.len() != 1 || site.alternate_allele.len() != 1 {
            continue; // SNP-only
        }
        if site.position < 1 || (site.position as usize) > length {
            continue; // off-contig
        }
        let idx = (site.position - 1) as usize;
        let c = counts.get(idx);
        let depth: u32 = c.iter().sum();
        let ref_bi = base_index(site.reference_allele.as_bytes()[0]);
        let alt_bi = base_index(site.alternate_allele.as_bytes()[0]);
        let ref_depth = ref_bi.map_or(0, |i| c[i]);
        let alt_depth = alt_bi.map_or(0, |i| c[i]);

        let (top_bi, top_count) = consensus(&c);
        let called = if depth < params.min_depth || top_count == 0 {
            CalledAllele::NoCall
        } else if (top_count as f64 / depth as f64) < params.min_allele_fraction {
            CalledAllele::NoCall
        } else if Some(top_bi) == alt_bi {
            CalledAllele::Alternate
        } else if Some(top_bi) == ref_bi {
            CalledAllele::Reference
        } else {
            CalledAllele::NoCall // consensus is a third allele
        };

        out.push(GenotypeCall {
            name: site.name.clone(),
            contig: site.contig.clone(),
            position: site.position,
            reference_allele: site.reference_allele.clone(),
            alternate_allele: site.alternate_allele.clone(),
            called,
            depth,
            ref_depth,
            alt_depth,
            allele_fraction: if depth == 0 { 0.0 } else { alt_depth as f64 / depth as f64 },
        });
    }
    Ok(out)
}

/// De-novo SNP discovery across `contig`: emit positions whose consensus base passes
/// the depth/fraction filters and differs from the reference base.
pub fn call_denovo(
    bam_path: &Path,
    reference_path: &Path,
    contig: &str,
    params: &HaploidCallerParams,
) -> Result<Vec<VariantCall>, AnalysisError> {
    let (length, counts) = tally_alleles(bam_path, contig, params, None)?;

    // Load the contig reference bases for ref-vs-consensus comparison.
    let mut fasta_reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(reference_path)
        .map_err(|e| AnalysisError::io(reference_path, e))?;
    let region: Region = contig
        .parse()
        .map_err(|_| AnalysisError::Message(format!("bad region for contig {contig}")))?;
    let ref_record = fasta_reader
        .query(&region)
        .map_err(|e| AnalysisError::io(reference_path, e))?;
    let ref_bases = ref_record.sequence().as_ref();

    let mut out = Vec::new();
    for idx in 0..length {
        let c = counts.get(idx);
        let depth: u32 = c.iter().sum();
        if depth < params.min_depth {
            continue;
        }
        let (top_bi, top_count) = consensus(&c);
        if top_count == 0 {
            continue;
        }
        let frac = top_count as f64 / depth as f64;
        if frac < params.min_allele_fraction {
            continue;
        }
        let ref_base = ref_bases.get(idx).copied().unwrap_or(b'N');
        // Skip non-ACGT reference (e.g. N) — not a confident SNP site.
        if base_index(ref_base) == Some(top_bi) || base_index(ref_base).is_none() {
            continue; // matches reference, or reference is N/ambiguous
        }
        out.push(VariantCall {
            contig: contig.to_string(),
            position: (idx + 1) as i64,
            reference_allele: ref_base.to_ascii_uppercase() as char,
            alternate_allele: BASES[top_bi] as char,
            depth,
            alt_depth: top_count,
            allele_fraction: frac,
        });
    }
    Ok(out)
}

/// Subtract known tree positions from de-novo calls to yield the private variant set
/// (the role `PrivateSnpProcessor` plays after liftover of the tree loci).
pub fn subtract_known(calls: &[VariantCall], known_positions: &HashSet<i64>) -> Vec<VariantCall> {
    calls
        .iter()
        .filter(|v| !known_positions.contains(&v.position))
        .cloned()
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn consensus_breaks_ties_toward_earlier_base() {
        assert_eq!(consensus(&[3, 3, 0, 0]), (0, 3)); // A wins tie vs C
        assert_eq!(consensus(&[0, 1, 5, 2]), (2, 5)); // G
        assert_eq!(consensus(&[0, 0, 0, 0]), (0, 0)); // empty -> A, count 0
    }

    #[test]
    fn subtract_known_removes_listed_positions() {
        let v = |p| VariantCall {
            contig: "chrM".into(),
            position: p,
            reference_allele: 'C',
            alternate_allele: 'A',
            depth: 4,
            alt_depth: 4,
            allele_fraction: 1.0,
        };
        let calls = vec![v(2), v(3), v(4)];
        let known: HashSet<i64> = [2, 3].into_iter().collect();
        let private = subtract_known(&calls, &known);
        assert_eq!(private.iter().map(|c| c.position).collect::<Vec<_>>(), vec![4]);
    }
}
