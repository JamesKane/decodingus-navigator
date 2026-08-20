//! Build the phased-haplotype reference asset (`ancestry_haps_<build>.bin`). It is the substrate
//! that the statistical phaser, and the parent-split chromosome painter, copy from.
//!
//! An AF panel reads the INFO allele counts, or collapses a GT to a 0/1/2 dosage. This builder is
//! different: it keeps the **phase**. Each `0|1` field of the phased 1000G matrix becomes two
//! separate haplotype bits.
//!
//! Every 1000G sample gives two haplotypes, and both carry its fine-population label, such as GBR
//! or YRI. Only a labelled sample enters, and the builder drops the related and unlabelled set. It
//! also takes only a biallelic SNV site, which matches the AIM painting loci that the matrix
//! already holds.
//!
//! The input is the **phased matrix of 1000G alone**, `$TMP/1kgp.matrix.tsv.gz`, where bcftools
//! `%GT` keeps the `|` separator, together with its sample list. It is NOT the combined matrix of
//! many sources. 1000G alone carries a phase. AADR gives a pseudo-haploid call, and SGDP gives no
//! phase. So those two sources must not enter the copying reference.

use std::collections::{BTreeMap, HashMap};
use std::io::BufRead;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::Parser;
use navigator_analysis::ancestry::{HapSite, HaplotypeReference};

use crate::pca::{first_base, load_fine_map, load_samples, open_maybe_gz, write_bin};

/// The reference build the emitted site coordinates are in (matches the other assets).
const BUILD: &str = "chm13v2.0";

#[derive(Parser)]
pub struct HapPanelArgs {
    /// PHASED genotype matrices, separated by commas, in the form `CHROM POS REF ALT [GT...]`. A GT
    /// keeps its `|` phase, as in `0|1`.
    ///
    /// More than one source, such as 1000G with HGDP, join over the **intersection** of their
    /// biallelic-SNV sites, where the ref and the alt match. Every source must carry a phase. A
    /// pseudo-haplotype from a source with no phase would be a noisier copy template, and it would
    /// bias the copying LAI against that source.
    #[arg(long)]
    matrix: String,
    /// Sample-id files, separated by commas, one for each `--matrix` source, in the same order.
    #[arg(long)]
    samples: String,
    /// A `sample<TAB>fine-pop` map, such as `NA12718  CEU` or `HGDP00511  French`. The builder
    /// drops a sample that this map does not hold. Both haplotypes of a sample carry its population
    /// label.
    #[arg(long)]
    pops: PathBuf,
    /// Output path for the bincode `HaplotypeReference` (`ancestry_haps_<build>.bin`).
    #[arg(long)]
    out: PathBuf,
}

/// `(contig, pos)` → `(ref, alt, per-sample (allele_a, allele_b))` for one phased source.
type PhasedSites = HashMap<(String, i64), (char, char, Vec<(u8, u8)>)>;

/// One phased reference source: its samples, and its phased genotype at each site.
struct Source {
    samples: Vec<String>,
    sites: PhasedSites,
}

fn split_list(s: &str) -> Vec<PathBuf> {
    s.split(',')
        .map(|p| PathBuf::from(p.trim()))
        .filter(|p| !p.as_os_str().is_empty())
        .collect()
}

/// Load a phased matrix + its sample list into a [`Source`] (biallelic SNVs only, dedup by position).
fn load_source(matrix: &Path, samples_path: &Path) -> Result<Source> {
    let samples = load_samples(samples_path)?;
    let n = samples.len();
    let mut sites: PhasedSites = HashMap::new();
    let mut skipped_non_snv = 0usize;
    for line in open_maybe_gz(matrix)?.lines() {
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
        let gts: Vec<&str> = f.collect();
        anyhow::ensure!(
            gts.len() == n,
            "{}:{} has {} genotype columns, expected {}",
            contig,
            pos,
            gts.len(),
            n
        );
        let per: Vec<(u8, u8)> = gts
            .iter()
            .map(|g| {
                let (a, b, _) = parse_gt_phased(g);
                (a, b)
            })
            .collect();
        sites
            .entry((contig, pos))
            .or_insert((first_base(ref_s), first_base(alt_s), per));
    }
    eprintln!(
        "  {} → {} samples, {} SNV sites ({} non-SNV skipped)",
        matrix.display(),
        n,
        sites.len(),
        skipped_non_snv
    );
    Ok(Source { samples, sites })
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
    let single = |s: &str| {
        s.len() == 1
            && s.chars()
                .all(|c| matches!(c.to_ascii_uppercase(), 'A' | 'C' | 'G' | 'T'))
    };
    single(ref_s) && single(alt_s)
}

pub fn build_hap_panel(args: HapPanelArgs) -> Result<()> {
    let matrices = split_list(&args.matrix);
    let sample_files = split_list(&args.samples);
    anyhow::ensure!(
        !matrices.is_empty() && matrices.len() == sample_files.len(),
        "need an equal, non-zero number of --matrix and --samples entries ({} vs {})",
        matrices.len(),
        sample_files.len()
    );
    let fine = load_fine_map(&args.pops)?;

    let sources: Vec<Source> = matrices
        .iter()
        .zip(&sample_files)
        .map(|(m, s)| load_source(m, s))
        .collect::<Result<_>>()?;

    // The site set holds the biallelic-SNV sites that EVERY source has, with a ref and an alt that
    // match. The allele codes then line up across the sources, and the code sorts the set by
    // position. The intersection keeps the packed reference rectangular, and the builder never has
    // to fill in a site that a source lacks.
    let (first, rest) = sources.split_first().expect("non-empty sources");
    let mut site_keys: Vec<(String, i64, char, char)> = first
        .sites
        .iter()
        .filter(|((c, p), (rf, alt, _))| {
            rest.iter().all(|s| {
                s.sites
                    .get(&(c.clone(), *p))
                    .is_some_and(|(r2, a2, _)| r2 == rf && a2 == alt)
            })
        })
        .map(|((c, p), (rf, alt, _))| (c.clone(), *p, *rf, *alt))
        .collect();
    site_keys.sort_by(|a, b| (a.0.as_str(), a.1).cmp(&(b.0.as_str(), b.1)));
    anyhow::ensure!(
        !site_keys.is_empty(),
        "no shared biallelic SNV sites across the sources"
    );
    let sites: Vec<HapSite> = site_keys
        .iter()
        .map(|(c, p, r, a)| HapSite {
            contig: c.clone(),
            position: *p,
            reference_allele: *r,
            alternate_allele: *a,
        })
        .collect();

    // The populations, in the order that they first appear, and the haplotype rows, which are two
    // for each labelled sample. Both join across the sources. The builder drops a sample that the
    // pop map does not hold.
    let mut populations: Vec<String> = Vec::new();
    let mut pop_index: BTreeMap<String, u16> = BTreeMap::new();
    let mut rows: Vec<Vec<u8>> = Vec::new();
    let mut hap_pop: Vec<u16> = Vec::new();
    let mut per_source_labelled = vec![0usize; sources.len()];
    for (src_i, src) in sources.iter().enumerate() {
        // The aligned genotype vector at each site, for this source. That is one map lookup for
        // each site, and not one for each sample.
        let aligned: Vec<&Vec<(u8, u8)>> = site_keys
            .iter()
            .map(|(c, p, _, _)| &src.sites[&(c.clone(), *p)].2)
            .collect();
        for (si, sample) in src.samples.iter().enumerate() {
            let Some(pop) = fine.get(sample) else {
                continue;
            };
            let next = populations.len() as u16;
            let pidx = *pop_index.entry(pop.clone()).or_insert_with(|| {
                populations.push(pop.clone());
                next
            });
            let mut a_row = Vec::with_capacity(sites.len());
            let mut b_row = Vec::with_capacity(sites.len());
            for gts in &aligned {
                let (a, b) = gts[si];
                a_row.push(a);
                b_row.push(b);
            }
            rows.push(a_row);
            rows.push(b_row);
            hap_pop.push(pidx);
            hap_pop.push(pidx);
            per_source_labelled[src_i] += 1;
        }
    }
    anyhow::ensure!(
        !rows.is_empty(),
        "no labelled samples across the sources (check {})",
        args.pops.display()
    );

    let reference = HaplotypeReference::from_rows(BUILD.to_string(), sites, populations, hap_pop, &rows);
    eprintln!(
        "hap panel: {} haplotypes × {} sites, {} populations (labelled samples per source: {:?})",
        reference.n_haplotypes,
        reference.n_sites,
        reference.populations.len(),
        per_source_labelled
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
        // Unphased separator still parses (defensive), phase just is not meaningful.
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
