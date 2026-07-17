//! Build the fine-grained (26-population) ancestry assets from a genotype matrix produced by
//! `bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'` over the 1000G genotype VCFs, plus
//! the sample order and sample→population map:
//!
//! * `pca`           — PCA loadings (per-SNP loadings+means, per-population centroids+variances).
//! * `fine-panel`    — an [`AncestryPanel`] with per-fine-population alt-allele frequencies.
//! * `ancient-panel` — an [`AncestryPanel`] over the deep ancestral sources (WHG/ANF/Steppe), with
//!   a per-population call floor so every retained site has genuine frequencies in every source.
//!
//! PCA uses the sample-space Gram matrix: with the centred genotype matrix `X` (samples × sites),
//! `X·Xᵀ = U·Σ²·Uᵀ`, so eigendecomposing the small Gram gives `U`/`Σ`; the per-SNP loadings are
//! `V = Xᵀ·U·Σ⁻¹` and reference sample coordinates `R = U·Σ`, from which each population's
//! centroid and per-component variance follow.

use std::collections::{BTreeSet, HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::Parser;
use flate2::read::MultiGzDecoder;
use nalgebra::{DMatrix, DVector, SymmetricEigen};
use navigator_analysis::ancestry::{AncestryPanel, PanelSite, PcaLoadings};

#[derive(Parser)]
pub struct PcaArgs {
    /// Genotype matrix/matrices `CHROM POS REF ALT GT...` per line (bcftools query), optionally
    /// .gz. Comma-separated to merge several panels by site (e.g. 1000G + SGDP).
    #[arg(long)]
    matrix: String,
    /// Sample-ID files (one per line), comma-separated and parallel to `--matrix`.
    #[arg(long)]
    samples: String,
    /// `sample<TAB>population` for every sample across the matrices.
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
    /// Projection mode: a file of population codes (one per line) whose samples build the PCA
    /// basis. All other labelled samples are *projected* onto that basis rather than shaping it
    /// — use it to keep sparse/biased ancient references (which would distort the axes) out of
    /// the decomposition while still placing them in PC space. Absent → every sample is basis.
    #[arg(long)]
    basis_pops: Option<PathBuf>,
}

#[derive(Parser)]
pub struct AncientPanelArgs {
    /// Genotype matrix/matrices `CHROM POS REF ALT GT...` per line, optionally .gz.
    /// Comma-separated to merge several panels by site.
    #[arg(long)]
    matrix: String,
    /// Sample-ID files (one per line), comma-separated and parallel to `--matrix`.
    #[arg(long)]
    samples: String,
    /// `sample<TAB>population` for every sample across the matrices (the pipeline's pop map;
    /// samples whose population isn't in `--components` are ignored).
    #[arg(long)]
    pops: PathBuf,
    /// The deep source populations, comma-separated and **in panel-axis order**
    /// (e.g. `WHG,ANF,Steppe`). Keep them non-collinear: `Steppe ≈ EHG+CHG`, so listing Steppe
    /// alongside EHG and CHG makes the mixture ill-conditioned.
    #[arg(long, default_value = "WHG,ANF,Steppe")]
    components: String,
    /// Output AncestryPanel (bincode).
    #[arg(long)]
    out: PathBuf,
    /// Keep a site only if **every** component has at least this many called samples there.
    /// This is the point of a separate ancient asset: ancient genomes are sparse, and a site with
    /// no calls in a source has no frequency — it must be dropped, not silently recorded as 0.0.
    #[arg(long, default_value_t = 8)]
    min_called: usize,
    /// **Ascertainment floor (Option A′).** Restrict the panel to the CHM13 `contig<TAB>pos` sites in
    /// this file — a consumer-array manifest. Allele-frequency admixture is only valid when the
    /// sample and the reference share ascertainment; the AADR/1240k universe includes capture sites
    /// consumer chips don't assay, and on those the deep estimate is unstable (a WGS sample reads
    /// ~90% Steppe where its own chip reads ~58%). Intersecting with the sites arrays actually assay
    /// makes the estimate agree across data sources. See `docs/design/ancient-ancestry-rebuild.md` §4.
    /// Optional: omit to build the full (unascertained) panel.
    #[arg(long)]
    ascertain_sites: Option<PathBuf>,
    /// Also write an inspection TSV (contig, pos, ref, alt, per-pop AF and called count).
    #[arg(long)]
    sites_tsv: Option<PathBuf>,
}

#[derive(Parser)]
pub struct FinePanelArgs {
    /// Genotype matrix/matrices `CHROM POS REF ALT GT...` per line, optionally .gz.
    /// Comma-separated to merge several panels by site.
    #[arg(long)]
    matrix: String,
    /// Sample-ID files (one per line), comma-separated and parallel to `--matrix`.
    #[arg(long)]
    samples: String,
    /// `sample<TAB>population` for every sample across the matrices.
    #[arg(long)]
    pops: PathBuf,
    /// Output AncestryPanel (bincode) with per-fine-population allele frequencies.
    #[arg(long)]
    out: PathBuf,
    /// Drop sites whose call rate across samples is below this.
    #[arg(long, default_value_t = 0.5)]
    min_call_rate: f64,
}

/// One matrix indexed by site: `(contig,pos) → (ref, alt, per-sample dosages)`.
type SiteMap = HashMap<(String, i64), (char, char, Vec<i8>)>;
/// Loaded + merged matrices: combined sample IDs, site metadata, and per-site dosage rows.
type LoadedMatrix = (Vec<String>, Vec<SiteMeta>, Vec<Vec<i8>>);

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
    Ok(s.lines()
        .map(|l| l.trim().to_string())
        .filter(|l| !l.is_empty())
        .collect())
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

/// A set of population codes from a file (one per line; `#` comments and blanks skipped).
fn load_pop_set(path: &Path) -> Result<HashSet<String>> {
    let mut s = String::new();
    open_maybe_gz(path)?.read_to_string(&mut s)?;
    Ok(s.lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(String::from)
        .collect())
}

/// Project sample `s` onto the basis loadings `v` (sites × k), centring each site by the basis
/// mean and accumulating `centered · loading`. Missing genotypes are skipped, then the result is
/// un-shrunk by `n_sites / used` — mirroring the runtime `project_pca`, so a sparse ancient
/// reference and a query sample land on the same scale as the basis coordinates.
fn project_sample(rows: &[Vec<i8>], s: usize, basis_means: &[f64], v: &DMatrix<f64>, k: usize) -> Vec<f64> {
    let mut coord = vec![0.0f64; k];
    let mut used = 0usize;
    for (j, row) in rows.iter().enumerate() {
        let d = row[s];
        if d < 0 {
            continue;
        }
        let centered = d as f64 - basis_means[j];
        used += 1;
        for (c, value) in coord.iter_mut().enumerate() {
            *value += centered * v[(j, c)];
        }
    }
    if used > 0 {
        let scale = rows.len() as f64 / used as f64;
        for value in coord.iter_mut() {
            *value *= scale;
        }
    }
    coord
}

/// Per-sample index into `pops` (its fine population), or `None` if unmapped.
fn sample_pop_index(samples: &[String], fine: &HashMap<String, String>, pops: &[String]) -> Vec<Option<usize>> {
    samples
        .iter()
        .map(|s| fine.get(s).and_then(|f| pops.iter().position(|p| p == f)))
        .collect()
}

/// Split a comma-separated path list (`a.tsv,b.tsv`) into paths.
fn split_paths(s: &str) -> Vec<PathBuf> {
    s.split(',')
        .map(|p| PathBuf::from(p.trim()))
        .filter(|p| !p.as_os_str().is_empty())
        .collect()
}

/// Parse one matrix into `(contig,pos) → (ref, alt, dosages)`, dedup by position (keep first).
fn load_one(path: &Path, n_samples: usize) -> Result<SiteMap> {
    let mut map: HashMap<(String, i64), (char, char, Vec<i8>)> = HashMap::new();
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
        let row: Vec<i8> = f.map(parse_gt).collect();
        anyhow::ensure!(
            row.len() == n_samples,
            "{}:{} has {} genotype columns, expected {}",
            contig,
            pos,
            row.len(),
            n_samples
        );
        map.entry((contig, pos)).or_insert((ref_allele, alt_allele, row));
    }
    Ok(map)
}

/// Load and merge one or more matrices by site: combined samples = concatenation of each file's
/// samples (in order); sites = those present in **all** matrices with combined call rate ≥
/// `min_call_rate`; dosages concatenated in the same order. Sorted by (contig, pos).
fn load_combined(matrices: &[PathBuf], sample_files: &[PathBuf], min_call_rate: f64) -> Result<LoadedMatrix> {
    anyhow::ensure!(
        !matrices.is_empty() && matrices.len() == sample_files.len(),
        "need an equal, non-zero number of --matrix and --samples entries"
    );
    let mut all_samples: Vec<String> = Vec::new();
    let mut maps: Vec<SiteMap> = Vec::new();
    for (m, s) in matrices.iter().zip(sample_files) {
        let samples = load_samples(s)?;
        let map = load_one(m, samples.len())?;
        eprintln!("  {} → {} samples, {} sites", m.display(), samples.len(), map.len());
        all_samples.extend(samples);
        maps.push(map);
    }
    let total_n = all_samples.len();

    let mut out: Vec<(SiteMeta, Vec<i8>)> = Vec::new();
    'sites: for (key, (rf, alt, _)) in &maps[0] {
        let mut combined = Vec::with_capacity(total_n);
        for map in &maps {
            match map.get(key) {
                Some((_, _, row)) => combined.extend_from_slice(row),
                None => continue 'sites, // not in every matrix
            }
        }
        let called = combined.iter().filter(|&&d| d >= 0).count();
        if (called as f64) < min_call_rate * total_n as f64 {
            continue;
        }
        out.push((
            SiteMeta {
                contig: key.0.clone(),
                pos: key.1,
                ref_allele: *rf,
                alt_allele: *alt,
            },
            combined,
        ));
    }
    out.sort_by(|a, b| (a.0.contig.as_str(), a.0.pos).cmp(&(b.0.contig.as_str(), b.0.pos)));
    eprintln!(
        "combined: {} samples, {} sites (call rate ≥ {min_call_rate})",
        total_n,
        out.len()
    );
    let (metas, rows): (Vec<_>, Vec<_>) = out.into_iter().unzip();
    Ok((all_samples, metas, rows))
}

pub fn build_pca(args: PcaArgs) -> Result<()> {
    let fine = load_fine_map(&args.pops)?;
    let (samples, metas, rows) = load_combined(
        &split_paths(&args.matrix),
        &split_paths(&args.samples),
        args.min_call_rate,
    )?;
    let n_samples = samples.len();
    anyhow::ensure!(n_samples > 0, "no samples");
    let pops = distinct_fine_pops(&samples, &fine);
    let sample_pop = sample_pop_index(&samples, &fine, &pops);
    let n_sites = metas.len();
    anyhow::ensure!(n_sites > 0, "no sites passed the call-rate filter");

    // Projection mode: only `basis_pops` samples build the PCA basis; all other labelled
    // samples are projected onto it. Absent → every sample is basis (original behaviour).
    let basis_set: Option<HashSet<String>> = match &args.basis_pops {
        Some(p) => Some(load_pop_set(p)?),
        None => None,
    };
    let is_basis = |s: usize| -> bool {
        match (&basis_set, sample_pop[s]) {
            (None, _) => true,
            (Some(set), Some(p)) => set.contains(&pops[p]),
            (Some(_), None) => false, // unlabelled samples can't anchor a basis
        }
    };
    let basis_idx: Vec<usize> = (0..n_samples).filter(|&s| is_basis(s)).collect();
    let n_basis = basis_idx.len();
    anyhow::ensure!(
        n_basis > 1,
        "need >1 basis sample (does --basis-pops match the pop labels?)"
    );
    if basis_set.is_some() {
        eprintln!(
            "projection mode: {n_basis} basis samples, {} projected",
            n_samples - n_basis
        );
    }
    let k = args.components.min(n_basis - 1).min(n_sites);

    // Per-site mean dosage over the BASIS samples only — the centring used both for the basis
    // decomposition and (stored in the asset) for projecting query samples at runtime.
    let mut basis_means = vec![0.0f64; n_sites];
    for (j, row) in rows.iter().enumerate() {
        let (sum, cnt) = basis_idx
            .iter()
            .map(|&s| row[s])
            .filter(|&d| d >= 0)
            .fold((0.0f64, 0usize), |(s, c), d| (s + d as f64, c + 1));
        basis_means[j] = if cnt > 0 { sum / cnt as f64 } else { 0.0 };
    }
    let means: Vec<f32> = basis_means.iter().map(|&m| m as f32).collect();

    // Centred basis matrix X_b (n_basis × sites), missing imputed to the basis mean (→ 0).
    let mut xb = DMatrix::<f64>::zeros(n_basis, n_sites);
    for (bi, &s) in basis_idx.iter().enumerate() {
        for (j, row) in rows.iter().enumerate() {
            let d = row[s];
            xb[(bi, j)] = if d >= 0 { d as f64 - basis_means[j] } else { 0.0 };
        }
    }

    eprintln!("computing {n_basis}×{n_basis} Gram + eigendecomposition…");
    let gram = &xb * xb.transpose();
    let eig = SymmetricEigen::new(gram);
    let mut order: Vec<usize> = (0..eig.eigenvalues.len()).collect();
    order.sort_by(|&a, &b| eig.eigenvalues[b].total_cmp(&eig.eigenvalues[a]));
    order.truncate(k);

    let mut uk = DMatrix::<f64>::zeros(n_basis, k);
    let mut sigma = DVector::<f64>::zeros(k);
    for (c, &idx) in order.iter().enumerate() {
        sigma[c] = eig.eigenvalues[idx].max(0.0).sqrt();
        uk.set_column(c, &eig.eigenvectors.column(idx));
    }

    // Loadings V = X_bᵀ·U·Σ⁻¹ (sites × k); basis coords R_b = U·Σ (n_basis × k).
    let mut v = xb.transpose() * &uk;
    for c in 0..k {
        if sigma[c] > 1e-9 {
            v.column_mut(c).scale_mut(1.0 / sigma[c]);
        }
    }
    let mut rb = uk.clone();
    for c in 0..k {
        rb.column_mut(c).scale_mut(sigma[c]);
    }

    // Unified per-sample coordinates: basis samples take their decomposition rows; every other
    // labelled sample is projected through V (centred by the basis means, with the same
    // missing-data un-shrink as the runtime `project_pca`, so ancient/query coords share a scale).
    let mut coords = DMatrix::<f64>::zeros(n_samples, k);
    for (bi, &s) in basis_idx.iter().enumerate() {
        for c in 0..k {
            coords[(s, c)] = rb[(bi, c)];
        }
    }
    for s in 0..n_samples {
        if is_basis(s) || sample_pop[s].is_none() {
            continue;
        }
        let projected = project_sample(&rows, s, &basis_means, &v, k);
        for (c, &val) in projected.iter().enumerate() {
            coords[(s, c)] = val;
        }
    }

    // Per-population centroid + diagonal variance over the unified coordinates.
    let n_pops = pops.len();
    let mut centroids = vec![0.0f32; n_pops * k];
    let mut variances = vec![1.0f32; n_pops * k];
    for p in 0..n_pops {
        let members: Vec<usize> = (0..n_samples).filter(|&s| sample_pop[s] == Some(p)).collect();
        if members.is_empty() {
            continue;
        }
        for c in 0..k {
            let vals: Vec<f64> = members.iter().map(|&s| coords[(s, c)]).collect();
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

    let loadings: Vec<f32> = (0..n_sites)
        .flat_map(|i| (0..k).map(move |c| (i, c)))
        .map(|(i, c)| v[(i, c)] as f32)
        .collect();
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
    eprintln!(
        "wrote {} ({n_sites} sites × {k} components, {n_pops} populations)",
        args.out.display()
    );
    Ok(())
}

pub fn build_fine_panel(args: FinePanelArgs) -> Result<()> {
    let fine = load_fine_map(&args.pops)?;
    let (samples, metas, rows) = load_combined(
        &split_paths(&args.matrix),
        &split_paths(&args.samples),
        args.min_call_rate,
    )?;
    let n_samples = samples.len();
    anyhow::ensure!(n_samples > 0, "no samples");
    let pops = distinct_fine_pops(&samples, &fine);
    let sample_pop = sample_pop_index(&samples, &fine, &pops);
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
                .map(|p| {
                    if called[p] > 0 {
                        (alt[p] / (2.0 * called[p] as f64)) as f32
                    } else {
                        0.0
                    }
                })
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

    let panel = AncestryPanel {
        build: "chm13v2.0".to_string(),
        populations: pops,
        sites,
    };
    write_bin(&args.out, &panel.to_bytes().map_err(|e| anyhow::anyhow!("{e}"))?)?;
    eprintln!(
        "wrote {} ({} sites × {n_pops} fine populations)",
        args.out.display(),
        panel.len()
    );
    Ok(())
}

/// Build the **ancient** deep-source frequency panel: per-site alt-allele frequency for each of
/// the deep sources (default WHG/ANF/Steppe), from the AADR genotype matrix.
///
/// This is deliberately a *separate asset* from `fine-panel` rather than a column subset of it.
/// `build_fine_panel` writes `0.0` for a population with no called samples at a site, which is
/// indistinguishable from a genuine "alt allele absent". For the 1000G fine populations that is
/// nearly harmless (they are called almost everywhere); for ancient sources it is fatal — they are
/// sparse and pseudo-haploid, so a large fraction of sites would enter the mixture as fake
/// "frequency 0" evidence and the fitted proportions would track *missingness* rather than
/// ancestry. Here a site survives only if **every** source has ≥ `min_called` calls, so every
/// frequency in the emitted panel is backed by real observations.
///
/// Pseudo-haploid genotypes (AADR emits one sampled allele as a homozygous diploid call) still give
/// an unbiased frequency: `E[dosage/2] = f`. Only the *variance* is inflated, which is why the call
/// floor — not the diploid coding — is what matters.
pub fn build_ancient_panel(args: AncientPanelArgs) -> Result<()> {
    let comps: Vec<String> = args
        .components
        .split(',')
        .map(|c| c.trim().to_string())
        .filter(|c| !c.is_empty())
        .collect();
    anyhow::ensure!(comps.len() >= 2, "need at least two source components");

    let pop_of = load_fine_map(&args.pops)?;
    // No global call-rate filter: the AADR matrix is mostly individuals we don't reference, so a
    // matrix-wide call rate says nothing about the sources. The per-component floor below is the
    // filter that matters.
    let (samples, metas, rows) = load_combined(&split_paths(&args.matrix), &split_paths(&args.samples), 0.0)?;
    anyhow::ensure!(!samples.is_empty(), "no samples");
    anyhow::ensure!(!metas.is_empty(), "no sites in the matrix");

    let sample_comp = sample_pop_index(&samples, &pop_of, &comps);
    let k = comps.len();
    let mut n_ref = vec![0usize; k];
    for c in sample_comp.iter().flatten() {
        n_ref[*c] += 1;
    }
    for (i, c) in comps.iter().enumerate() {
        anyhow::ensure!(n_ref[i] > 0, "component {c} has no samples in --pops");
    }

    // Optional ascertainment floor (Option A′): the CHM13 (contig, pos) a consumer array assays.
    let ascertained: Option<std::collections::HashSet<(String, i64)>> = match &args.ascertain_sites {
        Some(p) => {
            let text = std::fs::read_to_string(p).with_context(|| format!("reading {}", p.display()))?;
            let set: std::collections::HashSet<(String, i64)> = text
                .lines()
                .filter(|l| !l.starts_with('#') && !l.is_empty())
                .filter_map(|l| {
                    let mut it = l.split('\t');
                    let contig = it.next()?.trim();
                    let pos: i64 = it.next()?.trim().parse().ok()?;
                    (!contig.eq_ignore_ascii_case("contig")).then(|| (contig.to_string(), pos))
                })
                .collect();
            anyhow::ensure!(!set.is_empty(), "ascertainment file {} had no usable contig<TAB>pos rows", p.display());
            eprintln!("ascertainment floor: {} sites from {}", set.len(), p.display());
            Some(set)
        }
        None => None,
    };

    let mut sites = Vec::new();
    let mut tsv = match &args.sites_tsv {
        Some(p) => {
            let mut w = File::create(p).with_context(|| format!("creating {}", p.display()))?;
            writeln!(
                w,
                "contig\tpos\tref\talt\t{}",
                comps
                    .iter()
                    .map(|c| format!("af_{c}\tn_{c}"))
                    .collect::<Vec<_>>()
                    .join("\t")
            )?;
            Some(w)
        }
        None => None,
    };
    // Per-component running totals, for the build report.
    let mut called_total = vec![0usize; k];

    let mut dropped_unascertained = 0usize;
    for (m, row) in metas.iter().zip(&rows) {
        if let Some(set) = &ascertained {
            if !set.contains(&(m.contig.clone(), m.pos)) {
                dropped_unascertained += 1;
                continue;
            }
        }
        let mut alt = vec![0.0f64; k];
        let mut called = vec![0usize; k];
        for (i, &d) in row.iter().enumerate() {
            if d < 0 {
                continue;
            }
            if let Some(c) = sample_comp[i] {
                alt[c] += d as f64;
                called[c] += 1;
            }
        }
        if called.iter().any(|&n| n < args.min_called) {
            continue;
        }
        let freqs: Vec<f32> = (0..k).map(|c| (alt[c] / (2.0 * called[c] as f64)) as f32).collect();
        if let Some(w) = tsv.as_mut() {
            let cols: Vec<String> = (0..k).map(|c| format!("{:.4}\t{}", freqs[c], called[c])).collect();
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}",
                m.contig,
                m.pos,
                m.ref_allele,
                m.alt_allele,
                cols.join("\t")
            )?;
        }
        for c in 0..k {
            called_total[c] += called[c];
        }
        sites.push(PanelSite {
            contig: m.contig.clone(),
            position: m.pos,
            reference_allele: m.ref_allele,
            alternate_allele: m.alt_allele,
            freqs,
        });
    }
    anyhow::ensure!(
        !sites.is_empty(),
        "no site had ≥{} calls in every component — lower --min-called or widen the source groups",
        args.min_called
    );

    let panel = AncestryPanel {
        build: "chm13v2.0".to_string(),
        populations: comps.clone(),
        sites,
    };
    write_bin(&args.out, &panel.to_bytes().map_err(|e| anyhow::anyhow!("{e}"))?)?;
    let n = panel.len();
    eprintln!(
        "wrote {} ({n} of {} sites survived the ≥{}-call floor in all {k} sources{})",
        args.out.display(),
        metas.len(),
        args.min_called,
        if ascertained.is_some() {
            format!("; {dropped_unascertained} dropped off the ascertainment manifest")
        } else {
            String::new()
        }
    );
    for (i, c) in comps.iter().enumerate() {
        eprintln!(
            "  {c:<8} n={:<4} mean called/site {:.1}",
            n_ref[i],
            called_total[i] as f64 / n as f64
        );
    }
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

    /// Projecting a sample onto a 1-component basis: with loadings all 1.0 and basis mean 1.0,
    /// a fully-genotyped hom-alt sample lands at +n_sites; a half-missing one lands at the same
    /// place after the n_sites/used un-shrink (not pulled toward the origin).
    #[test]
    fn project_sample_centres_and_unshrinks() {
        // rows[site][sample]; one projected sample (index 0), 4 sites.
        let rows: Vec<Vec<i8>> = vec![vec![2], vec![2], vec![2], vec![2]];
        let means = vec![1.0; 4];
        let v = DMatrix::<f64>::from_element(4, 1, 1.0); // sites × k, all loadings 1.0
        let coord = project_sample(&rows, 0, &means, &v, 1);
        assert!((coord[0] - 4.0).abs() < 1e-9, "coord = {}", coord[0]); // (2-1)*1 × 4

        // Two of four sites missing → used=2, raw sum=2, scaled by 4/2 → 4 (same place).
        let sparse: Vec<Vec<i8>> = vec![vec![2], vec![2], vec![-1], vec![-1]];
        let coord = project_sample(&sparse, 0, &means, &v, 1);
        assert!((coord[0] - 4.0).abs() < 1e-9, "coord = {}", coord[0]);
    }

    #[test]
    fn pop_set_skips_comments_and_blanks() {
        let path = std::env::temp_dir().join(format!("panelbuild_pops_{}.txt", std::process::id()));
        fs::write(&path, "# header\nCEU\n\nYRI\n  TSI  \n").unwrap();
        let set = load_pop_set(&path).unwrap();
        let _ = fs::remove_file(&path);
        assert!(set.contains("CEU") && set.contains("YRI") && set.contains("TSI"));
        assert_eq!(set.len(), 3);
    }
}
