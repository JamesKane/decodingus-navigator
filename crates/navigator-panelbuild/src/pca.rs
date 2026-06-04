//! Build PCA loadings from a genotype matrix (bcftools `query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'`
//! over the 1000G genotype VCFs) plus the sample order and sample→population map.
//!
//! PCA via the sample-space Gram matrix: with the centred genotype matrix `X` (samples × sites),
//! `X·Xᵀ = U·Σ²·Uᵀ`, so eigendecomposing the small `samples × samples` Gram gives the sample
//! eigenvectors `U` and `Σ`. The per-SNP loadings (needed to project a *new* sample) are
//! `V = Xᵀ·U·Σ⁻¹`; reference sample coordinates are `R = U·Σ`, from which each population's
//! centroid and per-component variance follow. Output is a [`PcaLoadings`].

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::Parser;
use flate2::read::MultiGzDecoder;
use nalgebra::{DMatrix, DVector, SymmetricEigen};
use navigator_analysis::ancestry::PcaLoadings;

#[derive(Parser)]
pub struct PcaArgs {
    /// Genotype matrix `CHROM POS REF ALT GT...` per line (bcftools query), optionally .gz.
    #[arg(long)]
    matrix: PathBuf,
    /// Sample IDs, one per line, in the matrix's column order (bcftools query -l).
    #[arg(long)]
    samples: PathBuf,
    /// `sample<TAB>population` (fine 1000G pop, e.g. CEU) for per-population centroids.
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

/// Map a fine 1000G population code to its super-population (the panel axis), or `None`.
fn super_pop(fine: &str) -> Option<&'static str> {
    Some(match fine {
        "YRI" | "LWK" | "GWD" | "MSL" | "ESN" | "ASW" | "ACB" => "AFR",
        "MXL" | "PUR" | "CLM" | "PEL" => "AMR",
        "CHB" | "JPT" | "CHS" | "CDX" | "KHV" => "EAS",
        "CEU" | "TSI" | "FIN" | "GBR" | "IBS" => "EUR",
        "GIH" | "PJL" | "BEB" | "STU" | "ITU" => "SAS",
        _ => return None,
    })
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

pub fn build_pca(args: PcaArgs) -> Result<()> {
    // Sample order (matrix columns) + per-sample super-population.
    let samples: Vec<String> = {
        let mut s = String::new();
        open_maybe_gz(&args.samples)?.read_to_string(&mut s)?;
        s.lines().map(|l| l.trim().to_string()).filter(|l| !l.is_empty()).collect()
    };
    let fine: HashMap<String, String> = {
        let mut s = String::new();
        open_maybe_gz(&args.pops)?.read_to_string(&mut s)?;
        s.lines()
            .filter_map(|l| {
                let mut it = l.split_whitespace();
                Some((it.next()?.to_string(), it.next()?.to_string()))
            })
            .collect()
    };
    let n_samples = samples.len();
    anyhow::ensure!(n_samples > 0, "no samples");

    // Per-sample super-population index (into POPS), or None (excluded from centroids).
    let pops = crate::POPS;
    let sample_pop: Vec<Option<usize>> = samples
        .iter()
        .map(|s| {
            fine.get(s)
                .and_then(|f| super_pop(f))
                .and_then(|sp| pops.iter().position(|p| *p == sp))
        })
        .collect();

    // Parse the genotype matrix → site rows of dosages (missing = -1), dedup by (contig,pos).
    let mut sites: Vec<(String, i64)> = Vec::new();
    let mut rows: Vec<Vec<i8>> = Vec::new();
    let mut seen_pos: HashMap<(String, i64), ()> = HashMap::new();
    let mut dropped_callrate = 0usize;

    let reader = open_maybe_gz(&args.matrix)?;
    for line in reader.lines() {
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
        let _ref = f.next();
        let _alt = f.next();
        if seen_pos.insert((contig.clone(), pos), ()).is_some() {
            continue; // multiallelic split → keep first row
        }
        let mut row = Vec::with_capacity(n_samples);
        for gt in f {
            row.push(parse_gt(gt));
        }
        if row.len() != n_samples {
            anyhow::bail!(
                "{}:{} has {} genotype columns, expected {}",
                contig,
                pos,
                row.len(),
                n_samples
            );
        }
        let called = row.iter().filter(|&&d| d >= 0).count();
        if (called as f64) < args.min_call_rate * n_samples as f64 {
            dropped_callrate += 1;
            continue;
        }
        sites.push((contig, pos));
        rows.push(row);
    }
    let n_sites = sites.len();
    eprintln!(
        "matrix: {n_sites} sites kept ({dropped_callrate} dropped below call rate {}), {n_samples} samples",
        args.min_call_rate
    );
    anyhow::ensure!(n_sites > 0, "no sites passed the call-rate filter");
    let k = args.components.min(n_samples - 1).min(n_sites);

    // Per-site mean dosage (over called samples), then centred matrix X (samples × sites),
    // imputing missing genotypes to the site mean (→ centred value 0).
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

    // Sample-space Gram (samples × samples) → eigendecomposition.
    eprintln!("computing {n_samples}×{n_samples} Gram + eigendecomposition…");
    let gram = &x * x.transpose();
    let eig = SymmetricEigen::new(gram);

    // Top-k eigenpairs by eigenvalue, descending.
    let mut order: Vec<usize> = (0..eig.eigenvalues.len()).collect();
    order.sort_by(|&a, &b| eig.eigenvalues[b].total_cmp(&eig.eigenvalues[a]));
    order.truncate(k);

    // U_k (samples × k) and Σ (k).
    let mut uk = DMatrix::<f64>::zeros(n_samples, k);
    let mut sigma = DVector::<f64>::zeros(k);
    for (c, &idx) in order.iter().enumerate() {
        sigma[c] = eig.eigenvalues[idx].max(0.0).sqrt();
        uk.set_column(c, &eig.eigenvectors.column(idx));
    }

    // Loadings V = Xᵀ·U·Σ⁻¹ (sites × k); reference coords R = U·Σ (samples × k).
    let mut v = x.transpose() * &uk; // sites × k
    for c in 0..k {
        let s = sigma[c];
        if s > 1e-9 {
            v.column_mut(c).scale_mut(1.0 / s);
        }
    }
    let mut r = uk.clone();
    for c in 0..k {
        r.column_mut(c).scale_mut(sigma[c]);
    }

    // Per-population centroid + diagonal variance over reference sample coords.
    let n_pops = pops.len();
    let mut centroids = vec![0.0f32; n_pops * k];
    let mut variances = vec![1.0f32; n_pops * k];
    for (p, _code) in pops.iter().enumerate() {
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

    // Sanity check: per-population centroids on the first few PCs should be distinct.
    eprintln!("population centroids (PC1..PC3):");
    for (p, code) in pops.iter().enumerate() {
        let c1 = centroids[p * k];
        let c2 = if k > 1 { centroids[p * k + 1] } else { 0.0 };
        let c3 = if k > 2 { centroids[p * k + 2] } else { 0.0 };
        eprintln!("  {code}: PC1={c1:8.2} PC2={c2:8.2} PC3={c3:8.2}");
    }

    let loadings: Vec<f32> = (0..n_sites).flat_map(|i| (0..k).map(move |c| (i, c))).map(|(i, c)| v[(i, c)] as f32).collect();
    let pca = PcaLoadings {
        build: "chm13v2.0".to_string(),
        sites,
        means,
        n_components: k,
        loadings,
        populations: pops.iter().map(|s| s.to_string()).collect(),
        centroids,
        variances,
    };
    if let Some(parent) = args.out.parent() {
        fs::create_dir_all(parent).ok();
    }
    let bytes = pca.to_bytes().map_err(|e| anyhow::anyhow!("{e}"))?;
    fs::write(&args.out, &bytes).with_context(|| format!("writing {}", args.out.display()))?;
    eprintln!(
        "wrote {} ({} bytes): {} sites × {} components, {} populations",
        args.out.display(),
        bytes.len(),
        n_sites,
        k,
        n_pops
    );
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
    fn super_pop_mapping() {
        assert_eq!(super_pop("CEU"), Some("EUR"));
        assert_eq!(super_pop("YRI"), Some("AFR"));
        assert_eq!(super_pop("CHB"), Some("EAS"));
        assert_eq!(super_pop("ZZZ"), None);
    }
}
