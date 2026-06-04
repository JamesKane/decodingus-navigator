//! Build an ancestry reference panel from the 1000 Genomes-on-CHM13 VCFs.
//!
//! The VCFs (`1KGP.CHM13v2.0.chr*.recalibrated.snp_indel.pass.withafinfo.vcf.gz`) carry
//! per-super-population allele counts in INFO (`AC_<POP>_unrel`/`AN_<POP>_unrel` for
//! AFR/AMR/EAS/EUR/SAS), so per-population alt-allele frequency = `AC/AN` straight from INFO —
//! no per-sample genotype parsing. We keep biallelic SNPs where every population has data,
//! score each by Nei's Fst across the five populations, and emit the top-`max-sites` (above
//! `min-fst`) as a [`navigator_analysis::ancestry::AncestryPanel`] (bincode), plus a TSV for
//! inspection.
//!
//! Fetch the VCFs to a local mirror first (multi-GB; see docs/chm13-reference-resources.md):
//!   aws s3 cp --no-sign-request --recursive \
//!     s3://human-pangenomics/T2T/CHM13/assemblies/variants/1000_Genomes_Project/chm13v2.0/unrelated_samples_2504/allele_freq/ <dir>

use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::Parser;
use flate2::read::MultiGzDecoder;
use navigator_analysis::ancestry::{AncestryPanel, PanelSite};

/// The five 1000G super-populations, in panel-axis order.
const POPS: [&str; 5] = ["AFR", "AMR", "EAS", "EUR", "SAS"];

#[derive(Parser)]
#[command(about = "Build an ancestry AIMs panel from the 1000G-on-CHM13 VCFs")]
struct Args {
    /// Directory of `1KGP.CHM13v2.0.chr*.vcf.gz` files (a local mirror).
    #[arg(long)]
    vcf_dir: PathBuf,
    /// Output panel (bincode).
    #[arg(long)]
    out: PathBuf,
    /// Keep at most this many sites (highest Fst first).
    #[arg(long, default_value_t = 20_000)]
    max_sites: usize,
    /// Minimum Nei Fst across the five super-populations to keep a site.
    #[arg(long, default_value_t = 0.10)]
    min_fst: f64,
    /// Restrict to these chromosomes (comma-separated, e.g. `chr1,chr2`); default all.
    #[arg(long)]
    chroms: Option<String>,
    /// Also write an inspection TSV (contig, pos, ref, alt, fst, per-pop AF).
    #[arg(long)]
    sites_tsv: Option<PathBuf>,
}

/// A candidate site held in the bounded top-Fst heap.
#[derive(Debug, Clone)]
struct Candidate {
    contig: String,
    position: i64,
    reference_allele: char,
    alternate_allele: char,
    freqs: Vec<f32>,
    fst: f64,
}

// Min-heap ordering by Fst (so the smallest Fst is popped when the heap overflows).
impl PartialEq for Candidate {
    fn eq(&self, other: &Self) -> bool {
        self.fst == other.fst
    }
}
impl Eq for Candidate {}
impl PartialOrd for Candidate {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for Candidate {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.fst.total_cmp(&other.fst)
    }
}

fn main() -> Result<()> {
    let args = Args::parse();
    let want_chroms: Option<Vec<String>> = args
        .chroms
        .as_ref()
        .map(|s| s.split(',').map(|c| c.trim().to_string()).collect());

    let files = vcf_files(&args.vcf_dir, want_chroms.as_deref())
        .with_context(|| format!("listing VCFs in {}", args.vcf_dir.display()))?;
    if files.is_empty() {
        anyhow::bail!("no .vcf.gz files found in {}", args.vcf_dir.display());
    }

    // Bounded heap of the highest-Fst sites seen so far.
    let mut heap: BinaryHeap<Reverse<Candidate>> = BinaryHeap::new();
    let mut total_seen = 0usize;
    let mut total_kept = 0usize;

    for file in &files {
        eprintln!("scanning {}", file.display());
        let (seen, kept) = scan_file(file, args.min_fst, args.max_sites, &mut heap)?;
        total_seen += seen;
        total_kept += kept;
        eprintln!("  {seen} SNP sites with full population data, {} retained so far", heap.len());
    }

    let mut sites: Vec<Candidate> = heap.into_iter().map(|Reverse(c)| c).collect();
    sites.sort_by(|a, b| (a.contig.as_str(), a.position).cmp(&(b.contig.as_str(), b.position)));

    eprintln!(
        "selected {} sites (from {total_seen} eligible, {total_kept} above min-fst {})",
        sites.len(),
        args.min_fst
    );

    if let Some(tsv) = &args.sites_tsv {
        write_tsv(tsv, &sites).with_context(|| format!("writing {}", tsv.display()))?;
    }

    let panel = AncestryPanel {
        build: "chm13v2.0".to_string(),
        populations: POPS.iter().map(|s| s.to_string()).collect(),
        sites: sites
            .into_iter()
            .map(|c| PanelSite {
                contig: c.contig,
                position: c.position,
                reference_allele: c.reference_allele,
                alternate_allele: c.alternate_allele,
                freqs: c.freqs,
            })
            .collect(),
    };
    if let Some(parent) = args.out.parent() {
        fs::create_dir_all(parent).ok();
    }
    let bytes = panel.to_bytes().map_err(|e| anyhow::anyhow!("{e}"))?;
    fs::write(&args.out, &bytes).with_context(|| format!("writing {}", args.out.display()))?;
    eprintln!("wrote {} ({} bytes) with {} sites", args.out.display(), bytes.len(), panel.len());
    Ok(())
}

/// Per-chromosome VCF files in `dir`, optionally restricted to `chroms` (by filename token).
fn vcf_files(dir: &Path, chroms: Option<&[String]>) -> Result<Vec<PathBuf>> {
    let mut files: Vec<PathBuf> = fs::read_dir(dir)?
        .filter_map(|e| e.ok().map(|e| e.path()))
        .filter(|p| {
            let name = p.file_name().and_then(|n| n.to_str()).unwrap_or("");
            let is_vcf = name.ends_with(".vcf.gz") || name.ends_with(".vcf");
            let chrom_ok = match chroms {
                Some(cs) => cs.iter().any(|c| name.contains(&format!(".{c}."))),
                None => true,
            };
            is_vcf && chrom_ok
        })
        .collect();
    files.sort();
    Ok(files)
}

/// Scan one VCF, pushing eligible high-Fst SNP sites into the bounded heap. Returns
/// `(sites_with_full_population_data, sites_above_min_fst)`.
fn scan_file(
    path: &Path,
    min_fst: f64,
    max_sites: usize,
    heap: &mut BinaryHeap<Reverse<Candidate>>,
) -> Result<(usize, usize)> {
    let file = File::open(path)?;
    let reader: Box<dyn BufRead> = if path.extension().and_then(|e| e.to_str()) == Some("gz") {
        Box::new(BufReader::new(MultiGzDecoder::new(file)))
    } else {
        Box::new(BufReader::new(file))
    };

    let mut seen = 0usize;
    let mut kept = 0usize;
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') || line.is_empty() {
            continue;
        }
        let Some(cand) = parse_record(&line) else { continue };
        seen += 1;
        if cand.fst < min_fst {
            continue;
        }
        kept += 1;
        heap.push(Reverse(cand));
        if heap.len() > max_sites {
            heap.pop(); // drop the lowest-Fst site
        }
    }
    Ok((seen, kept))
}

/// Parse a VCF data line into a [`Candidate`], or `None` if it's not a biallelic SNP with full
/// per-population allele-count data.
fn parse_record(line: &str) -> Option<Candidate> {
    let mut f = line.split('\t');
    let contig = f.next()?;
    let position: i64 = f.next()?.parse().ok()?;
    let _id = f.next()?;
    let reference = f.next()?;
    let alternate = f.next()?;
    let _qual = f.next()?;
    let _filter = f.next()?;
    let info = f.next()?;

    let r = single_base(reference)?;
    let a = single_base(alternate)?;
    let freqs = parse_info_freqs(info, &POPS)?;
    let fst = nei_fst(&freqs);
    Some(Candidate {
        contig: contig.to_string(),
        position,
        reference_allele: r,
        alternate_allele: a,
        freqs,
        fst,
    })
}

/// A single A/C/G/T base (upper-cased), or `None` for multi-base / non-ACGT / multiallelic.
fn single_base(s: &str) -> Option<char> {
    let mut chars = s.chars();
    let c = chars.next()?.to_ascii_uppercase();
    if chars.next().is_some() {
        return None; // length > 1 (indel) or multiallelic "G,T"
    }
    matches!(c, 'A' | 'C' | 'G' | 'T').then_some(c)
}

/// Per-population alt-allele frequency from `AC_<POP>_unrel`/`AN_<POP>_unrel` in INFO. Returns
/// `None` unless every population has `AN > 0` and a parseable `AC`.
fn parse_info_freqs(info: &str, pops: &[&str]) -> Option<Vec<f32>> {
    let map: HashMap<&str, &str> = info
        .split(';')
        .filter_map(|kv| kv.split_once('='))
        .collect();
    let mut freqs = Vec::with_capacity(pops.len());
    for pop in pops {
        let ac: f64 = map.get(format!("AC_{pop}_unrel").as_str())?.parse().ok()?;
        let an: f64 = map.get(format!("AN_{pop}_unrel").as_str())?.parse().ok()?;
        if an <= 0.0 {
            return None;
        }
        freqs.push((ac / an).clamp(0.0, 1.0) as f32);
    }
    Some(freqs)
}

/// Nei's Fst across populations from alt-allele frequencies (equal population weights):
/// `(Ht - Hs) / Ht`, with `Hs` the mean within-population expected heterozygosity and `Ht`
/// the total expected heterozygosity at the mean frequency. 0 when invariant across pops.
fn nei_fst(freqs: &[f32]) -> f64 {
    let k = freqs.len() as f64;
    if k == 0.0 {
        return 0.0;
    }
    let pbar: f64 = freqs.iter().map(|&p| p as f64).sum::<f64>() / k;
    let hs: f64 = freqs.iter().map(|&p| 2.0 * p as f64 * (1.0 - p as f64)).sum::<f64>() / k;
    let ht = 2.0 * pbar * (1.0 - pbar);
    if ht <= 0.0 {
        0.0
    } else {
        ((ht - hs) / ht).clamp(0.0, 1.0)
    }
}

fn write_tsv(path: &Path, sites: &[Candidate]) -> Result<()> {
    let mut w = File::create(path)?;
    write!(w, "contig\tposition\tref\talt\tfst")?;
    for pop in POPS {
        write!(w, "\taf_{pop}")?;
    }
    writeln!(w)?;
    for c in sites {
        write!(
            w,
            "{}\t{}\t{}\t{}\t{:.4}",
            c.contig, c.position, c.reference_allele, c.alternate_allele, c.fst
        )?;
        for f in &c.freqs {
            write!(w, "\t{f:.4}")?;
        }
        writeln!(w)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_per_population_frequencies_from_info() {
        let info = "AC=100;AF=0.5;AN=200;AC_AFR_unrel=180;AN_AFR_unrel=200;\
                    AC_AMR_unrel=10;AN_AMR_unrel=200;AC_EAS_unrel=5;AN_EAS_unrel=200;\
                    AC_EUR_unrel=2;AN_EUR_unrel=200;AC_SAS_unrel=20;AN_SAS_unrel=200";
        let freqs = parse_info_freqs(info, &POPS).unwrap();
        assert_eq!(freqs.len(), 5);
        assert!((freqs[0] - 0.90).abs() < 1e-5); // AFR
        assert!((freqs[3] - 0.01).abs() < 1e-5); // EUR
    }

    #[test]
    fn missing_population_data_is_rejected() {
        // No SAS counts → not eligible.
        let info = "AC_AFR_unrel=1;AN_AFR_unrel=2;AC_AMR_unrel=1;AN_AMR_unrel=2;\
                    AC_EAS_unrel=1;AN_EAS_unrel=2;AC_EUR_unrel=1;AN_EUR_unrel=2";
        assert!(parse_info_freqs(info, &POPS).is_none());
    }

    #[test]
    fn fst_is_high_when_populations_diverge_and_zero_when_uniform() {
        // One population fixed alt, the rest fixed ref → strong differentiation.
        let diverged = nei_fst(&[1.0, 0.0, 0.0, 0.0, 0.0]);
        assert!(diverged > 0.4, "fst = {diverged}");
        // Identical everywhere → no differentiation.
        let uniform = nei_fst(&[0.5, 0.5, 0.5, 0.5, 0.5]);
        assert!(uniform.abs() < 1e-9, "fst = {uniform}");
    }

    #[test]
    fn parse_record_keeps_snps_and_drops_indels() {
        let snp = "chr1\t1000\t.\tA\tG\t.\tPASS\t\
                   AC_AFR_unrel=180;AN_AFR_unrel=200;AC_AMR_unrel=10;AN_AMR_unrel=200;\
                   AC_EAS_unrel=5;AN_EAS_unrel=200;AC_EUR_unrel=2;AN_EUR_unrel=200;\
                   AC_SAS_unrel=20;AN_SAS_unrel=200";
        let c = parse_record(snp).unwrap();
        assert_eq!((c.contig.as_str(), c.position, c.reference_allele, c.alternate_allele), ("chr1", 1000, 'A', 'G'));

        let indel = "chr1\t1000\t.\tA\tAG\t.\tPASS\tAC_AFR_unrel=1;AN_AFR_unrel=2";
        assert!(parse_record(indel).is_none());
    }
}
