//! Build the phased-haplotype reference asset (`ancestry_haps_<build>.bin`) — the substrate the
//! statistical phaser and the parent-split chromosome painter copy from.
//!
//! Unlike the AF panels (which read INFO allele counts or collapse GT to a 0/1/2 dosage), this
//! builder keeps the **phase**: each `0|1` field of the phased 1000G matrix becomes two separate
//! haplotype bits. Every 1000G sample contributes two haplotypes carrying its fine-population label
//! (e.g. GBR, YRI). Only labelled samples enter (the related/unlabelled set is dropped), and only
//! biallelic SNV sites — matching the AIM painting loci the matrix was already sliced to.
//!
//! Input is the **phased, 1000G-only** matrix `$TMP/1kgp.matrix.tsv.gz` (bcftools `%GT` retains the
//! `|` separator) + its sample list, NOT the combined multi-source matrix: only 1000G is phased
//! (AADR is pseudo-haploid, SGDP unphased), so those sources must not enter the copying reference.

use std::collections::{BTreeSet, HashSet};
use std::io::BufRead;
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Parser;
use navigator_analysis::ancestry::{HapSite, HaplotypeReference};

use crate::pca::{first_base, load_fine_map, load_samples, open_maybe_gz, write_bin};

/// The reference build the emitted site coordinates are in (matches the other assets).
const BUILD: &str = "chm13v2.0";

#[derive(Parser)]
pub struct HapPanelArgs {
    /// Phased genotype matrix (`CHROM POS REF ALT [GT...]`, GT cells keep their `|` phase, e.g.
    /// `0|1`). Pass ONLY the phased 1000G matrix, not the combined multi-source one.
    #[arg(long)]
    matrix: PathBuf,
    /// Sample IDs (one per line), positionally aligned to the matrix GT columns.
    #[arg(long)]
    samples: PathBuf,
    /// `sample<TAB>fine-pop` map (e.g. `NA12718  CEU`). Samples absent here are dropped.
    #[arg(long)]
    pops: PathBuf,
    /// Output path for the bincode `HaplotypeReference` (`ancestry_haps_<build>.bin`).
    #[arg(long)]
    out: PathBuf,
    /// Drop a site unless at least this fraction of contributing samples are non-missing (phased
    /// 1000G is fully called, so this only guards against a corrupt slice).
    #[arg(long, default_value_t = 0.99)]
    min_call_rate: f64,
}

/// The two phased alleles of a diploid GT field `a|b` (or `a/b`), each `0` (ref) or `1` (any
/// non-ref), plus whether either allele was missing. A missing allele defaults to ref so the
/// reference stays rectangular (phased 1000G has no true missing calls).
fn parse_gt_phased(gt: &str) -> (u8, u8, bool) {
    let mut it = gt.split(['|', '/']);
    let (a, ma) = allele(it.next());
    let (b, mb) = allele(it.next());
    (a, b, ma || mb)
}

fn allele(s: Option<&str>) -> (u8, bool) {
    match s {
        Some("0") => (0, false),
        None | Some("") | Some(".") => (0, true),
        Some(_) => (1, false),
    }
}

/// A biallelic SNV has single-base ref and alt over `{A,C,G,T}` (skips indels / multiallelic).
fn is_biallelic_snv(ref_s: &str, alt_s: &str) -> bool {
    let single = |s: &str| s.len() == 1 && s.chars().all(|c| matches!(c.to_ascii_uppercase(), 'A' | 'C' | 'G' | 'T'));
    single(ref_s) && single(alt_s)
}

pub fn build_hap_panel(args: HapPanelArgs) -> Result<()> {
    let samples = load_samples(&args.samples)?;
    let fine = load_fine_map(&args.pops)?;

    // Labelled sample columns (index into the matrix GT columns) with their fine population, in
    // matrix order. Unlabelled samples (the related set) are excluded entirely.
    let labelled: Vec<(usize, String)> = samples
        .iter()
        .enumerate()
        .filter_map(|(i, s)| fine.get(s).map(|p| (i, p.clone())))
        .collect();
    anyhow::ensure!(
        !labelled.is_empty(),
        "no matrix samples are in the pop map {} — nothing to build",
        args.pops.display()
    );

    // Population axis (sorted distinct fine pops) and per-haplotype label (two haplotypes/sample).
    let populations: Vec<String> = labelled.iter().map(|(_, p)| p.clone()).collect::<BTreeSet<_>>().into_iter().collect();
    let pop_index = |p: &str| populations.iter().position(|x| x == p).unwrap() as u16;
    let mut hap_pop: Vec<u16> = Vec::with_capacity(labelled.len() * 2);
    for (_, p) in &labelled {
        let idx = pop_index(p);
        hap_pop.push(idx);
        hap_pop.push(idx);
    }
    let n_hap = labelled.len() * 2;
    let n_samples = samples.len();

    // One allele row per haplotype, grown site-by-site; site metadata in parallel.
    let mut rows: Vec<Vec<u8>> = vec![Vec::new(); n_hap];
    let mut sites: Vec<HapSite> = Vec::new();
    let mut seen: HashSet<(String, i64)> = HashSet::new();
    let mut skipped_non_snv = 0usize;
    let mut skipped_low_call = 0usize;

    for line in open_maybe_gz(&args.matrix)?.lines() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let mut f = line.split('\t');
        let contig = f.next().unwrap_or("").to_string();
        let pos: i64 = match f.next().and_then(|s| s.parse().ok()) {
            Some(p) => p,
            None => continue,
        };
        let ref_s = f.next().unwrap_or("N");
        let alt_s = f.next().unwrap_or("N");
        if !is_biallelic_snv(ref_s, alt_s) {
            skipped_non_snv += 1;
            continue;
        }
        if !seen.insert((contig.clone(), pos)) {
            continue; // dedup by position (keep first)
        }
        let gts: Vec<&str> = f.collect();
        anyhow::ensure!(
            gts.len() == n_samples,
            "{}:{} has {} genotype columns, expected {}",
            contig,
            pos,
            gts.len(),
            n_samples
        );

        // Resolve the two phased alleles for each labelled sample.
        let mut per_sample: Vec<(u8, u8)> = Vec::with_capacity(labelled.len());
        let mut missing = 0usize;
        for (ci, _) in &labelled {
            let (a, b, miss) = parse_gt_phased(gts[*ci]);
            if miss {
                missing += 1;
            }
            per_sample.push((a, b));
        }
        if (labelled.len() - missing) < (args.min_call_rate * labelled.len() as f64) as usize {
            skipped_low_call += 1;
            continue;
        }

        for (si, (a, b)) in per_sample.iter().enumerate() {
            rows[si * 2].push(*a);
            rows[si * 2 + 1].push(*b);
        }
        sites.push(HapSite {
            contig,
            position: pos,
            reference_allele: first_base(ref_s),
            alternate_allele: first_base(alt_s),
        });
    }

    anyhow::ensure!(
        !sites.is_empty(),
        "no usable phased biallelic SNV sites in {} (skipped {skipped_non_snv} non-SNV, {skipped_low_call} low-call)",
        args.matrix.display()
    );

    let reference = HaplotypeReference::from_rows(BUILD.to_string(), sites, populations, hap_pop, &rows);
    eprintln!(
        "hap panel: {} haplotypes × {} sites, {} populations (skipped {skipped_non_snv} non-SNV, {skipped_low_call} low-call rows)",
        reference.n_haplotypes,
        reference.n_sites,
        reference.populations.len(),
    );
    let bytes = reference
        .to_bytes()
        .map_err(|e| anyhow::anyhow!("encoding hap reference: {e}"))?;
    write_bin(&args.out, &bytes).with_context(|| format!("writing {}", args.out.display()))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_gt_phased_keeps_phase() {
        assert_eq!(parse_gt_phased("0|1"), (0, 1, false));
        assert_eq!(parse_gt_phased("1|0"), (1, 0, false));
        assert_eq!(parse_gt_phased("1|1"), (1, 1, false));
        assert_eq!(parse_gt_phased("0|0"), (0, 0, false));
        // Unphased separator still parses (defensive), phase just isn't meaningful.
        assert_eq!(parse_gt_phased("0/1"), (0, 1, false));
        // Missing → ref with the missing flag set.
        assert_eq!(parse_gt_phased(".|1"), (0, 1, true));
        assert_eq!(parse_gt_phased(".|."), (0, 0, true));
        // Multiallelic allele index counts as non-ref (matrix should be biallelic, but be safe).
        assert_eq!(parse_gt_phased("0|2"), (0, 1, false));
    }

    #[test]
    fn snv_filter() {
        assert!(is_biallelic_snv("A", "G"));
        assert!(is_biallelic_snv("c", "t"));
        assert!(!is_biallelic_snv("AC", "G")); // indel
        assert!(!is_biallelic_snv("A", "G,T")); // multiallelic
        assert!(!is_biallelic_snv("A", "N"));
    }
}
