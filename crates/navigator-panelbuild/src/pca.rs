//! Build the fine-grained (26-population) ancestry assets from a genotype matrix produced by
//! `bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'` over the 1000G genotype VCFs, plus
//! the sample order and sample→population map:
//!
//! * `pca`        — PCA loadings (per-SNP loadings+means, per-population centroids+variances).
//! * `fine-panel` — an [`AncestryPanel`] with per-fine-population alt-allele frequencies.
//!
//! PCA uses the sample-space Gram matrix: with the centred genotype matrix `X` (samples × sites),
//! `X·Xᵀ = U·Σ²·Uᵀ`, so eigendecomposing the small Gram gives `U`/`Σ`; the per-SNP loadings are
//! `V = Xᵀ·U·Σ⁻¹` and reference sample coordinates `R = U·Σ`, from which each population's
//! centroid and per-component variance follow.

use std::collections::{BTreeSet, HashMap};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::Parser;
use flate2::read::MultiGzDecoder;
use nalgebra::{DMatrix, DVector, SymmetricEigen};
use navigator_analysis::ancestry::{AncestryPanel, PanelSite, PcaLoadings};

#[derive(Parser)]
pub struct PcaArgs {
    /// Genotype matrix `CHROM POS REF ALT GT...` per line (bcftools query), optionally .gz.
    #[arg(long)]
    matrix: PathBuf,
    /// Sample IDs, one per line, in the matrix's column order (bcftools query -l).
    #[arg(long)]
    samples: PathBuf,
    /// `sample<TAB>population` (fine 1000G pop, e.g. CEU).
    #[arg(long)]
    pops: PathBuf,
    /// Output PcaLoadings (bincode).
    #[arg(long)]
    out: PathBuf,
    /// Number of principal components to retain.
    #[arg(long, default_value_t = 10)]
    components: usize,
    /// Drop sites whose call rate across samples is below this.
    #[arg(long, default_value_t = 0.9)]
    min_call_rate: f64,
}

#[derive(Parser)]
pub struct FinePanelArgs {
    /// Genotype matrix `CHROM POS REF ALT GT...` per line (bcftools query), optionally .gz.
    #[arg(long)]
    matrix: PathBuf,
    /// Sample IDs, one per line, in the matrix's column order.
    #[arg(long)]
    samples: PathBuf,
    /// `sample<TAB>population` (fine 1000G pop).
    #[arg(long)]
    pops: PathBuf,
    /// Output AncestryPanel (bincode) with per-fine-population allele frequencies.
    #[arg(long)]
    out: PathBuf,
    /// Drop sites whose call rate across samples is below this.
    #[arg(long, default_value_t = 0.5)]
    min_call_rate: f64,
}

/// A genotyped site: coordinates + the biallelic ref/alt the genotypes are relative to.
struct SiteMeta {
    contig: String,
    pos: i64,
    ref_allele: char,
    alt_allele: char,
}

/// Diploid alt-allele dosage from a VCF GT field: 0/1/2, or -1 for a no-call. Counts non-ref
/// alleles (any index > 0), so multiallelic sites collapse to "carries a non-ref allele".
fn parse_gt(gt: &str) -> i8 {
    let mut dosage = 0i8;
    let mut seen = false;
    for a in gt.split(['|', '/']) {
        seen = true;
        match a {
            "." | "" => return -1,
            "0" => {}
            _ => dosage += 1,
        }
    }
    if seen {
        dosage.min(2)
    } else {
        -1
    }
}

fn open_maybe_gz(path: &Path) -> Result<Box<dyn BufRead>> {
    let f = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(f))))
    } else {
        Ok(Box::new(BufReader::new(f)))
    }
}

fn first_base(s: &str) -> char {
    s.chars().next().map(|c| c.to_ascii_uppercase()).unwrap_or('N')
}

fn load_samples(path: &Path) -> Result<Vec<String>> {
    let mut s = String::new();
    open_maybe_gz(path)?.read_to_string(&mut s)?;
    Ok(s.lines().map(|l| l.trim().to_string()).filter(|l| !l.is_empty()).collect())
}

/// `sample → fine population` (e.g. NA12718 → CEU).
fn load_fine_map(path: &Path) -> Result<HashMap<String, String>> {
    let mut s = String::new();
    open_maybe_gz(path)?.read_to_string(&mut s)?;
    Ok(s.lines()
        .filter_map(|l| {
            let mut it = l.split_whitespace();
            Some((it.next()?.to_string(), it.next()?.to_string()))
        })
        .collect())
}

/// The distinct fine populations present among `samples`, sorted for determinism.
fn distinct_fine_pops(samples: &[String], fine: &HashMap<String, String>) -> Vec<String> {
    let set: BTreeSet<String> = samples.iter().filter_map(|s| fine.get(s).cloned()).collect();
    set.into_iter().collect()
}

/// Per-sample index into `pops` (its fine population), or `None` if unmapped.
fn sample_pop_index(samples: &[String], fine: &HashMap<String, String>, pops: &[String]) -> Vec<Option<usize>> {
    samples
        .iter()
        .map(|s| fine.get(s).and_then(|f| pops.iter().position(|p| p == f)))
        .collect()
}

/// Parse the matrix → site metadata + dosage rows (missing = -1), dedup by (contig,pos),
/// keeping sites whose call rate clears `min_call_rate`.
fn load_matrix(path: &Path, n_samples: usize, min_call_rate: f64) -> Result<(Vec<SiteMeta>, Vec<Vec<i8>>)> {
    let mut metas: Vec<SiteMeta> = Vec::new();
    let mut rows: Vec<Vec<i8>> = Vec::new();
    let mut seen: HashMap<(String, i64), ()> = HashMap::new();
    let mut dropped = 0usize;

    for line in open_maybe_gz(path)?.lines() {
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
        let ref_allele = first_base(f.next().unwrap_or("N"));
        let alt_allele = first_base(f.next().unwrap_or("N"));
        if seen.insert((contig.clone(), pos), ()).is_some() {
            continue; // multiallelic split → keep first row
        }
        let row: Vec<i8> = f.map(parse_gt).collect();
        anyhow::ensure!(
            row.len() == n_samples,
            "{}:{} has {} genotype columns, expected {}",
            contig,
            pos,
            row.len(),
            n_samples
        );
        let called = row.iter().filter(|&&d| d >= 0).count();
        if (called as f64) < min_call_rate * n_samples as f64 {
            dropped += 1;
            continue;
        }
        metas.push(SiteMeta { contig, pos, ref_allele, alt_allele });
        rows.push(row);
    }
    eprintln!("matrix: {} sites kept ({dropped} below call rate {min_call_rate})", metas.len());
    Ok((metas, rows))
}

pub fn build_pca(args: PcaArgs) -> Result<()> {
    let samples = load_samples(&args.samples)?;
    let fine = load_fine_map(&args.pops)?;
    let n_samples = samples.len();
    anyhow::ensure!(n_samples > 0, "no samples");
    let pops = distinct_fine_pops(&samples, &fine);
    let sample_pop = sample_pop_index(&samples, &fine, &pops);

    let (metas, rows) = load_matrix(&args.matrix, n_samples, args.min_call_rate)?;
    let n_sites = metas.len();
    anyhow::ensure!(n_sites > 0, "no sites passed the call-rate filter");
    let k = args.components.min(n_samples - 1).min(n_sites);

    // Per-site mean dosage (over called samples); centred matrix X (samples × sites), imputing
    // missing genotypes to the mean (→ centred 0).
    let mut means = vec![0.0f32; n_sites];
    let mut x = DMatrix::<f64>::zeros(n_samples, n_sites);
    for (j, row) in rows.iter().enumerate() {
        let (sum, cnt) = row.iter().filter(|&&d| d >= 0).fold((0.0f64, 0usize), |(s, c), &d| (s + d as f64, c + 1));
        let mean = if cnt > 0 { sum / cnt as f64 } else { 0.0 };
        means[j] = mean as f32;
        for (i, &d) in row.iter().enumerate() {
            x[(i, j)] = if d >= 0 { d as f64 - mean } else { 0.0 };
        }
    }

    eprintln!("computing {n_samples}×{n_samples} Gram + eigendecomposition…");
    let gram = &x * x.transpose();
    let eig = SymmetricEigen::new(gram);
    let mut order: Vec<usize> = (0..eig.eigenvalues.len()).collect();
    order.sort_by(|&a, &b| eig.eigenvalues[b].total_cmp(&eig.eigenvalues[a]));
    order.truncate(k);

    let mut uk = DMatrix::<f64>::zeros(n_samples, k);
    let mut sigma = DVector::<f64>::zeros(k);
    for (c, &idx) in order.iter().enumerate() {
        sigma[c] = eig.eigenvalues[idx].max(0.0).sqrt();
        uk.set_column(c, &eig.eigenvectors.column(idx));
    }

    // Loadings V = Xᵀ·U·Σ⁻¹ (sites × k); reference coords R = U·Σ (samples × k).
    let mut v = x.transpose() * &uk;
    for c in 0..k {
        if sigma[c] > 1e-9 {
            v.column_mut(c).scale_mut(1.0 / sigma[c]);
        }
    }
    let mut r = uk.clone();
    for c in 0..k {
        r.column_mut(c).scale_mut(sigma[c]);
    }

    // Per-population centroid + diagonal variance over reference sample coordinates.
    let n_pops = pops.len();
    let mut centroids = vec![0.0f32; n_pops * k];
    let mut variances = vec![1.0f32; n_pops * k];
    for p in 0..n_pops {
        let members: Vec<usize> = (0..n_samples).filter(|&s| sample_pop[s] == Some(p)).collect();
        if members.is_empty() {
            continue;
        }
        for c in 0..k {
            let vals: Vec<f64> = members.iter().map(|&s| r[(s, c)]).collect();
            let mean = vals.iter().sum::<f64>() / vals.len() as f64;
            let var = if vals.len() > 1 {
                vals.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / (vals.len() as f64 - 1.0)
            } else {
                1.0
            };
            centroids[p * k + c] = mean as f32;
            variances[p * k + c] = (var.max(1e-6)) as f32;
        }
    }

    eprintln!("population centroids (PC1..PC3):");
    for (p, code) in pops.iter().enumerate() {
        let c2 = if k > 1 { centroids[p * k + 1] } else { 0.0 };
        let c3 = if k > 2 { centroids[p * k + 2] } else { 0.0 };
        eprintln!("  {code}: PC1={:8.2} PC2={c2:8.2} PC3={c3:8.2}", centroids[p * k]);
    }

    let loadings: Vec<f32> = (0..n_sites).flat_map(|i| (0..k).map(move |c| (i, c))).map(|(i, c)| v[(i, c)] as f32).collect();
    let pca = PcaLoadings {
        build: "chm13v2.0".to_string(),
        sites: metas.iter().map(|m| (m.contig.clone(), m.pos)).collect(),
        means,
        n_components: k,
        loadings,
        populations: pops,
        centroids,
        variances,
    };
    write_bin(&args.out, &pca.to_bytes().map_err(|e| anyhow::anyhow!("{e}"))?)?;
    eprintln!("wrote {} ({n_sites} sites × {k} components, {n_pops} populations)", args.out.display());
    Ok(())
}

pub fn build_fine_panel(args: FinePanelArgs) -> Result<()> {
    let samples = load_samples(&args.samples)?;
    let fine = load_fine_map(&args.pops)?;
    let n_samples = samples.len();
    anyhow::ensure!(n_samples > 0, "no samples");
    let pops = distinct_fine_pops(&samples, &fine);
    let sample_pop = sample_pop_index(&samples, &fine, &pops);

    let (metas, rows) = load_matrix(&args.matrix, n_samples, args.min_call_rate)?;
    anyhow::ensure!(!metas.is_empty(), "no sites passed the call-rate filter");

    // Per-site, per-population alt-allele frequency = Σ dosage / (2 · called) within the pop.
    let n_pops = pops.len();
    let sites: Vec<PanelSite> = metas
        .iter()
        .zip(&rows)
        .map(|(m, row)| {
            let mut alt = vec![0.0f64; n_pops];
            let mut called = vec![0usize; n_pops];
            for (i, &d) in row.iter().enumerate() {
                if d < 0 {
                    continue;
                }
                if let Some(p) = sample_pop[i] {
                    alt[p] += d as f64;
                    called[p] += 1;
                }
            }
            let freqs = (0..n_pops)
                .map(|p| if called[p] > 0 { (alt[p] / (2.0 * called[p] as f64)) as f32 } else { 0.0 })
                .collect();
            PanelSite {
                contig: m.contig.clone(),
                position: m.pos,
                reference_allele: m.ref_allele,
                alternate_allele: m.alt_allele,
                freqs,
            }
        })
        .collect();

    let panel = AncestryPanel { build: "chm13v2.0".to_string(), populations: pops, sites };
    write_bin(&args.out, &panel.to_bytes().map_err(|e| anyhow::anyhow!("{e}"))?)?;
    eprintln!("wrote {} ({} sites × {n_pops} fine populations)", args.out.display(), panel.len());
    Ok(())
}

fn write_bin(out: &Path, bytes: &[u8]) -> Result<()> {
    if let Some(parent) = out.parent() {
        fs::create_dir_all(parent).ok();
    }
    fs::write(out, bytes).with_context(|| format!("writing {}", out.display()))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gt_parsing() {
        assert_eq!(parse_gt("0|0"), 0);
        assert_eq!(parse_gt("0/0"), 0);
        assert_eq!(parse_gt("0|1"), 1);
        assert_eq!(parse_gt("1/0"), 1);
        assert_eq!(parse_gt("1|1"), 2);
        assert_eq!(parse_gt("1|2"), 2); // multiallelic → capped
        assert_eq!(parse_gt("./."), -1);
        assert_eq!(parse_gt("."), -1);
    }

    #[test]
    fn distinct_pops_are_sorted_and_indexed() {
        let samples = vec!["a".to_string(), "b".to_string(), "c".to_string()];
        let fine: HashMap<String, String> = [("a", "CEU"), ("b", "YRI"), ("c", "CEU")]
            .into_iter()
            .map(|(s, p)| (s.to_string(), p.to_string()))
            .collect();
        let pops = distinct_fine_pops(&samples, &fine);
        assert_eq!(pops, vec!["CEU".to_string(), "YRI".to_string()]);
        let idx = sample_pop_index(&samples, &fine, &pops);
        assert_eq!(idx, vec![Some(0), Some(1), Some(0)]);
    }
}
