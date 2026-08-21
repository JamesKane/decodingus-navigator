//! Build the fine-grained ancestry assets, over 26 populations. The inputs are a genotype matrix,
//! the sample order, and a sample→population map. `bcftools query -f
//! '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n'` makes that matrix from the 1000G genotype VCFs. The assets
//! are:
//!
//! * `pca` gives the PCA loadings: a loading and a mean for each SNP, and a centroid and a
//!   variance for each population.
//! * `fine-panel` gives an [`AncestryPanel`] with the alt-allele frequency of each fine population.
//! * `ancient-panel` gives an [`AncestryPanel`] over the deep ancestral sources, WHG, ANF, and
//!   Steppe. It applies a call floor to each population, so every site that stays has a true
//!   frequency in every source.
//!
//! The PCA works on the Gram matrix in sample space. Take the centred genotype matrix `X`, of
//! samples by sites. Then `X·Xᵀ = U·Σ²·Uᵀ`, so an eigendecomposition of the small Gram gives `U`
//! and `Σ`. The loading of each SNP is `V = Xᵀ·U·Σ⁻¹`, and the reference sample coordinates are
//! `R = U·Σ`. The centroid of each population, and its variance on each component, follow from
//! those.

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
    /// One genotype matrix, or more, with `CHROM POS REF ALT GT...` on each line, from a bcftools
    /// query. Each one can be .gz. Separate them with commas to merge more than one panel by site,
    /// such as 1000G with SGDP.
    #[arg(long)]
    matrix: String,
    /// Sample-ID files, with one id on each line, separated by commas, in the order of
    /// `--matrix`.
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
    /// Projection mode. This is a file of population codes, one on each line, whose samples build
    /// the PCA basis. Every other labelled sample *projects* onto that basis, and does not shape
    /// it.
    ///
    /// Use it to keep a sparse or biased ancient reference out of the decomposition, where it
    /// would bend the axes. Such a reference still gets a place in PC space. With no file, every
    /// sample builds the basis.
    #[arg(long)]
    basis_pops: Option<PathBuf>,
}

#[derive(Parser)]
pub struct AncientPanelArgs {
    /// One genotype matrix, or more, with `CHROM POS REF ALT GT...` on each line, and each one can
    /// be .gz. Separate them with commas to merge more than one panel by site.
    #[arg(long)]
    matrix: String,
    /// Sample-ID files, with one id on each line, separated by commas, in the order of
    /// `--matrix`.
    #[arg(long)]
    samples: String,
    /// A `sample<TAB>population` line for every sample in the matrices. This is the pipeline's pop
    /// map. The builder ignores a sample whose population is not in `--components`.
    #[arg(long)]
    pops: PathBuf,
    /// The deep source (**left**) populations, comma-separated and **in panel-axis order**
    /// (e.g. `WHG,ANF,Steppe`). Keep them non-collinear: `Steppe ≈ EHG+CHG`, so listing Steppe
    /// alongside EHG and CHG makes the mixture ill-conditioned.
    #[arg(long, default_value = "WHG,ANF,Steppe")]
    components: String,
    /// The qpAdm **outgroup**, or right, populations, separated by commas. They go on the panel
    /// axis after the sources, as in `YRI,CHB,GIH,Karitiana,Papuan,Onge`.
    ///
    /// f4 admixture measures how many alleles the target has in common *against* these. They carry
    /// no weight, but each one must relate to the sources by a different amount.
    ///
    /// They have their own call floor, `--outgroup-min-called`, which is lower. A good outgroup is
    /// often one high-quality genome, and Onge with n near 2 is normal.
    ///
    /// With no value, the builder makes a panel of sources alone, which is the old frequency-EM
    /// asset. See documents/design/ancient-ancestry-rebuild.md §7.4.
    #[arg(long, default_value = "")]
    outgroups: String,
    /// Output AncestryPanel (bincode).
    #[arg(long)]
    out: PathBuf,
    /// Keep a site only when **every source** has this many called samples there, or more.
    ///
    /// This is why the ancient asset is separate. An ancient genome is sparse, and a site with no
    /// call in a source has no frequency. The builder must drop such a site. It must not record
    /// 0.0, and say nothing.
    #[arg(long, default_value_t = 8)]
    min_called: usize,
    /// The call floor for an **outgroup**, see `--outgroups`. It is separate from the source floor
    /// `--min-called`, and it is lower.
    ///
    /// A qpAdm outgroup is correctly small: a few present-day genomes for each lineage. The f4
    /// jackknife handles the noise in the frequency. A site stays only when every source passes
    /// `--min-called` **and** every outgroup passes this floor.
    #[arg(long, default_value_t = 2)]
    outgroup_min_called: usize,
    /// **The ascertainment floor (Option A′).** Hold the panel to the CHM13 `contig<TAB>pos` sites
    /// in this file, which is the manifest of a consumer array.
    ///
    /// Allele-frequency admixture holds only when the sample and the reference share their
    /// ascertainment. The AADR and 1240k universe covers capture sites that a consumer chip does
    /// not assay. On those sites the deep estimate is unstable: a WGS sample reads about 90%
    /// Steppe, where that person's own chip reads about 58%.
    ///
    /// An intersection with the sites that an array assays makes the estimate agree across data
    /// sources. See `documents/design/ancient-ancestry-rebuild.md` §4. It is optional: leave it out
    /// to build the full panel, with no ascertainment.
    #[arg(long)]
    ascertain_sites: Option<PathBuf>,
    /// Also write a TSV to read: contig, pos, ref, alt, and the AF and call count of each
    /// population.
    #[arg(long)]
    sites_tsv: Option<PathBuf>,
    /// The CHM13 reference FASTA, with its `.fai` index beside it.
    ///
    /// With this file, the builder orients each site so that its `reference_allele` is the **real
    /// CHM13 base**. Where the input labels are the other way round, it exchanges ref and alt, and
    /// changes each freq to 1−freq.
    ///
    /// Without it, a panel that comes from a *lifted* sites file keeps the allele labels of the
    /// source build. About 30% of those are the other way round from CHM13. Such a panel is
    /// consistent with itself, for its own fit, but nothing can join it to the other assets, which
    /// are CHM13-canonical (docs §7.16).
    #[arg(long)]
    reference: Option<PathBuf>,
}

#[derive(Parser)]
pub struct FinePanelArgs {
    /// One genotype matrix, or more, with `CHROM POS REF ALT GT...` on each line, and each one
    /// can be .gz. Separate them with commas to merge more than one panel by site.
    #[arg(long)]
    matrix: String,
    /// Sample-ID files, with one id on each line, separated by commas, in the order of
    /// `--matrix`.
    #[arg(long)]
    samples: String,
    /// `sample<TAB>population` for every sample across the matrices.
    #[arg(long)]
    pops: PathBuf,
    /// The output AncestryPanel, in bincode, with the allele frequency of each fine population.
    #[arg(long)]
    out: PathBuf,
    /// Drop sites whose call rate across samples is below this.
    #[arg(long, default_value_t = 0.5)]
    min_call_rate: f64,
}

/// One matrix indexed by site: `(contig,pos) → (ref, alt, per-sample dosages)`.
type SiteMap = HashMap<(String, i64), (char, char, Vec<i8>)>;
/// The matrices after a load and a merge: the combined sample IDs, the site metadata, and one
/// dosage row for each site.
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

pub(crate) fn open_maybe_gz(path: &Path) -> Result<Box<dyn BufRead>> {
    let f = File::open(path).with_context(|| format!("opening {}", path.display()))?;
    if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(f))))
    } else {
        Ok(Box::new(BufReader::new(f)))
    }
}

pub(crate) fn first_base(s: &str) -> char {
    s.chars().next().map(|c| c.to_ascii_uppercase()).unwrap_or('N')
}

pub(crate) fn load_samples(path: &Path) -> Result<Vec<String>> {
    let mut s = String::new();
    open_maybe_gz(path)?.read_to_string(&mut s)?;
    Ok(s.lines()
        .map(|l| l.trim().to_string())
        .filter(|l| !l.is_empty())
        .collect())
}

/// `sample → fine population` (e.g. NA12718 → CEU).
pub(crate) fn load_fine_map(path: &Path) -> Result<HashMap<String, String>> {
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

/// A set of population codes from a file, with one code on each line. The parser skips a `#`
/// comment and a blank line.
fn load_pop_set(path: &Path) -> Result<HashSet<String>> {
    let mut s = String::new();
    open_maybe_gz(path)?.read_to_string(&mut s)?;
    Ok(s.lines()
        .map(str::trim)
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .map(String::from)
        .collect())
}

/// Project sample `s` onto the basis loadings `v`, of sites by k, and centre each site by the basis
/// mean.
///
/// It calls the runtime's [`navigator_analysis::ancestry::project_centered`]. So a sparse ancient
/// reference, and a query sample, land on the same scale as the basis coordinates. The un-shrink
/// policy then has exactly one definition.
fn project_sample(rows: &[Vec<i8>], s: usize, basis_means: &[f64], v: &DMatrix<f64>, k: usize) -> Vec<f64> {
    let centered = rows
        .iter()
        .enumerate()
        .filter(|(_, row)| row[s] >= 0)
        .map(|(j, row)| (j, row[s] as f64 - basis_means[j]));
    navigator_analysis::ancestry::project_centered(rows.len(), k, centered, |j, c| v[(j, c)])
}

/// The index of each sample into `pops`, which is its fine population. It is `None` for a sample
/// that the map does not hold.
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

/// Load one matrix, or more, and merge them by site.
///
/// The combined samples are the samples of each file, one file after another, in order. The sites
/// are the ones that **every** matrix holds, and whose combined call rate reaches `min_call_rate`.
/// The dosages follow the same order as the samples. The output sorts by (contig, pos).
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

    // Projection mode. The `basis_pops` samples alone build the PCA basis, and every other
    // labelled sample projects onto it. With no file, every sample builds the basis, which is the
    // first behaviour.
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

    // The mean dosage at each site, over the BASIS samples alone. The basis decomposition uses
    // that centre, and so does the projection of a query sample at run time, which reads it from
    // the asset.
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

    // One set of coordinates for every sample. A basis sample takes its row from the
    // decomposition. Every other labelled sample projects through V, centred by the basis means,
    // with the same un-shrink for missing data that the runtime `project_pca` applies. So the
    // ancient coordinates and the query coordinates share one scale.
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

    // The centroid and the diagonal variance of each population, over the unified coordinates.
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

    // The alt-allele frequency at each site, for each population: Σ dosage / (2 · called) inside
    // that population.
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

/// Build the **ancient** deep-source frequency panel: the alt-allele frequency at each site, for
/// each deep source. The default sources are WHG, ANF, and Steppe, and the input is the AADR
/// genotype matrix.
///
/// This is a *separate asset* from `fine-panel`, on purpose, and not a subset of its columns.
/// `build_fine_panel` writes `0.0` for a population with no called sample at a site, and nothing
/// separates that from a true "no alt allele here".
///
/// For a 1000G fine population that does almost no harm, because such a population has a call
/// almost everywhere. For an ancient source it is fatal. An ancient source is sparse and
/// pseudo-haploid, so a large part of the sites would enter the mixture as false "frequency 0"
/// evidence. The fitted proportions would then follow the *missing data*, and not the ancestry.
///
/// Here a site stays only when **every** source has `min_called` calls or more. So a real
/// observation stands behind every frequency in the panel.
///
/// A pseudo-haploid genotype still gives an unbiased frequency, and AADR writes one sampled allele
/// as a homozygous diploid call: `E[dosage/2] = f`. Only the *variance* grows. That is why the call
/// floor is what matters, and not the diploid form.
pub fn build_ancient_panel(args: AncientPanelArgs) -> Result<()> {
    let parse_list = |s: &str| -> Vec<String> {
        s.split(',')
            .map(|c| c.trim().to_string())
            .filter(|c| !c.is_empty())
            .collect()
    };
    let sources: Vec<String> = parse_list(&args.components);
    let outgroup_comps: Vec<String> = parse_list(&args.outgroups);
    anyhow::ensure!(sources.len() >= 2, "need at least two source components");
    // Panel axis = sources first, then outgroups. The estimator designates roles by index (the
    // committed qpadm_leftpops / rightpops); the asset is just a plain frequency panel over both.
    let comps: Vec<String> = sources.iter().chain(&outgroup_comps).cloned().collect();
    let n_src = sources.len();
    anyhow::ensure!(
        comps.iter().collect::<std::collections::HashSet<_>>().len() == comps.len(),
        "a population appears in both --components and --outgroups"
    );
    // The call floor of each population. A source uses --min-called, and an outgroup uses the
    // lower --outgroup-min-called.
    let floor: Vec<usize> = (0..comps.len())
        .map(|i| {
            if i < n_src {
                args.min_called
            } else {
                args.outgroup_min_called
            }
        })
        .collect();

    let pop_of = load_fine_map(&args.pops)?;
    // There is no call-rate filter over the whole matrix. Most individuals in the AADR matrix are
    // ones that nothing here references, so a call rate across the matrix says nothing about the
    // sources. The floor of each component, below, is the filter that matters.
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
        let role = if i < n_src { "source" } else { "outgroup" };
        anyhow::ensure!(n_ref[i] > 0, "{role} component {c} has no samples in --pops");
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
            anyhow::ensure!(
                !set.is_empty(),
                "ascertainment file {} had no usable contig<TAB>pos rows",
                p.display()
            );
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
    // The cumulative total of each component, for the build report.
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
        if (0..k).any(|c| called[c] < floor[c]) {
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
        "no site cleared every population's call floor — lower --min-called/--outgroup-min-called or widen the groups",
    );

    // Orient every site, so that reference_allele equals the real CHM13 base (docs §7.16). The
    // sites come in matrix order, sorted by contig and then by pos, so the code loads the sequence
    // of each contig once.
    //
    // Where the labels are the other way round, and CHM13 carries the labelled ALT, exchange ref
    // and alt, and change each freq to 1−freq. Where neither allele equals the base, which is a
    // mismatch from the liftover or the alleles, drop the site.
    if let Some(ref_path) = &args.reference {
        let mut cur = String::new();
        let mut seq: Vec<u8> = Vec::new();
        let (mut flipped, mut dropped) = (0usize, 0usize);
        let mut oriented = Vec::with_capacity(sites.len());
        for mut s in sites.into_iter() {
            if s.contig != cur {
                seq = navigator_analysis::reader::read_contig_sequence(ref_path, &s.contig)
                    .map_err(|e| anyhow::anyhow!("reading {} from {}: {e}", s.contig, ref_path.display()))?;
                cur = s.contig.clone();
            }
            let base = seq
                .get((s.position - 1) as usize)
                .map(|b| b.to_ascii_uppercase() as char)
                .unwrap_or('N');
            if base == s.reference_allele {
                // already canonical
            } else if base == s.alternate_allele {
                std::mem::swap(&mut s.reference_allele, &mut s.alternate_allele);
                for f in s.freqs.iter_mut() {
                    *f = 1.0 - *f;
                }
                flipped += 1;
            } else {
                dropped += 1;
                continue;
            }
            oriented.push(s);
        }
        eprintln!(
            "CHM13-oriented {} sites: {flipped} ref/alt-swapped, {dropped} dropped (base matched neither allele)",
            oriented.len()
        );
        sites = oriented;
        anyhow::ensure!(
            !sites.is_empty(),
            "no site survived CHM13 orientation — wrong reference?"
        );
    }

    let panel = AncestryPanel {
        build: "chm13v2.0".to_string(),
        populations: comps.clone(),
        sites,
    };
    write_bin(&args.out, &panel.to_bytes().map_err(|e| anyhow::anyhow!("{e}"))?)?;
    let n = panel.len();
    eprintln!(
        "wrote {} ({n} of {} sites cleared the floor in all {n_src} sources (≥{}) + {} outgroups (≥{}){})",
        args.out.display(),
        metas.len(),
        args.min_called,
        outgroup_comps.len(),
        args.outgroup_min_called,
        if ascertained.is_some() {
            format!("; {dropped_unascertained} dropped off the ascertainment manifest")
        } else {
            String::new()
        }
    );
    for (i, c) in comps.iter().enumerate() {
        let role = if i < n_src { "src" } else { "out" };
        eprintln!(
            "  [{role}] {c:<12} n={:<4} mean called/site {:.1}",
            n_ref[i],
            called_total[i] as f64 / n as f64
        );
    }
    Ok(())
}

pub(crate) fn write_bin(out: &Path, bytes: &[u8]) -> Result<()> {
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

    /// A projection of a sample onto a basis of one component. Every loading is 1.0, and the basis
    /// mean is 1.0. A hom-alt sample with a genotype at every site lands at +n_sites. A sample that
    /// is half missing lands at the same place, after the n_sites/used un-shrink, and nothing pulls
    /// it toward the origin.
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
