//! Build the two **Tier B** assets — the African-outgroup position track and the genome-wide
//! lineage-classification track. Design: `documents/design/ArchaicAncestry_Design.md` §4 (assets 2
//! and 3), §5 (how they are used).
//!
//! Both are delta-varint position streams (see `navigator_analysis::archaic::PositionStream`),
//! because the outgroup track alone is ~67 M positions genome-wide and must neither be held in
//! memory nor shipped as raw integers.
//!
//! As everywhere else in this pipeline, `bcftools` does the VCF decoding in `08_build_archaic.sh`
//! and this crate consumes tab-separated tables, so the logic stays pure and testable.

use std::collections::BTreeMap;
use std::io::BufRead;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::Parser;
use navigator_analysis::archaic::{
    ArchaicCallable, ArchaicClassify, ArchaicOutgroup, CallableContig, ClassifyContig, DiagnosticClass,
    PositionStream,
};

use crate::pca::{open_maybe_gz, write_bin};

const BUILD: &str = "chm13v2.0";

#[derive(Parser)]
pub struct ArchaicOutgroupArgs {
    /// `CHROM<TAB>POS` of every site variable in the African outgroup, on the target build.
    /// Produced by stage 08 from the 1000G-on-CHM13 `AC_AFR_unrel` INFO field.
    #[arg(long)]
    sites: PathBuf,
    /// Minimum outgroup allele count a site needed to be listed (recorded in the asset for
    /// provenance — the filtering itself happens upstream in bcftools).
    #[arg(long, default_value_t = 1)]
    min_allele_count: u32,
    /// Output (bincode `ArchaicOutgroup`).
    #[arg(long)]
    out: PathBuf,
}

/// Read `CHROM<TAB>POS` into per-contig sorted, deduplicated position lists.
fn load_positions(path: &Path) -> Result<BTreeMap<String, Vec<i64>>> {
    let rdr = open_maybe_gz(path).with_context(|| format!("opening {}", path.display()))?;
    let mut by_contig: BTreeMap<String, Vec<i64>> = BTreeMap::new();
    for line in rdr.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let mut it = line.split_whitespace();
        let (Some(c), Some(p)) = (it.next(), it.next()) else { continue };
        let Ok(pos) = p.parse::<i64>() else { continue };
        by_contig.entry(c.to_string()).or_default().push(pos);
    }
    for v in by_contig.values_mut() {
        v.sort_unstable();
        v.dedup();
    }
    Ok(by_contig)
}

pub fn build_archaic_outgroup(args: ArchaicOutgroupArgs) -> Result<()> {
    let by_contig = load_positions(&args.sites)?;
    anyhow::ensure!(!by_contig.is_empty(), "no outgroup sites read from {}", args.sites.display());

    let mut contigs = Vec::with_capacity(by_contig.len());
    let mut total = 0usize;
    for (contig, positions) in &by_contig {
        total += positions.len();
        contigs.push(PositionStream::encode(contig, positions));
    }

    let asset = ArchaicOutgroup {
        build: BUILD.to_string(),
        min_allele_count: args.min_allele_count,
        contigs,
    };
    let encoded = asset.to_bytes().map_err(|e| anyhow::anyhow!("{e}"))?;
    eprintln!(
        "outgroup: {total} positions over {} contigs -> {:.1} MB ({:.2} bytes/site)",
        asset.contigs.len(),
        encoded.len() as f64 / 1_048_576.0,
        encoded.len() as f64 / total.max(1) as f64
    );
    write_bin(&args.out, &encoded)?;
    eprintln!("wrote {}", args.out.display());
    Ok(())
}

#[derive(Parser)]
pub struct ArchaicClassifyArgs {
    /// The candidates table from `archaic-candidates` (GRCh37 payload, carries the per-archaic calls).
    #[arg(long)]
    candidates: PathBuf,
    /// The lifted BED (target build) for those candidates.
    #[arg(long)]
    lifted: PathBuf,
    /// Output (bincode `ArchaicClassify`).
    #[arg(long)]
    out: PathBuf,
    /// Restrict to these contigs (comma-separated, e.g. `chr21,chr22`) — for fast iteration.
    #[arg(long)]
    contigs: Option<String>,
}

pub fn build_archaic_classify(args: ArchaicClassifyArgs) -> Result<()> {
    // The classification track is the genome-wide superset of the marker panel: every polarized
    // candidate with an archaic hom-derived call, BEFORE the frequency filters that select markers.
    // Segment attribution wants maximum diagnostic density, not the ascertained marker subset.
    let want: Option<Vec<String>> = args
        .contigs
        .as_ref()
        .map(|s| s.split(',').map(|c| c.trim().to_string()).collect());

    let rdr = open_maybe_gz(&args.candidates)?;
    let mut payload: BTreeMap<usize, (char, DiagnosticClass)> = BTreeMap::new();
    for line in rdr.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 7 {
            continue;
        }
        let Ok(idx) = f[0].parse::<usize>() else { continue };
        let derived = f[5].chars().next().unwrap_or('N');
        // Re-derive the class from the stored per-genome call tokens: index 3 is Denisova.
        let calls: Vec<char> = f[6].chars().collect();
        let carries = |c: Option<&char>| matches!(c, Some('1') | Some('2'));
        let nea = calls.iter().take(3).any(|c| carries(Some(c)));
        let den = carries(calls.get(3));
        let class = match (nea, den) {
            (true, false) => DiagnosticClass::Neanderthal,
            (false, true) => DiagnosticClass::Denisovan,
            _ => DiagnosticClass::SharedArchaic,
        };
        payload.insert(idx, (derived, class));
    }

    // Join onto the lifted coordinates.
    let rdr = open_maybe_gz(&args.lifted)?;
    let mut by_contig: BTreeMap<String, Vec<(i64, char, DiagnosticClass)>> = BTreeMap::new();
    for line in rdr.lines() {
        let line = line?;
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 4 {
            continue;
        }
        let (Ok(end), Ok(idx)) = (f[2].parse::<i64>(), f[3].trim().parse::<usize>()) else {
            continue;
        };
        if let Some(w) = &want {
            if !w.iter().any(|c| c == f[0]) {
                continue;
            }
        }
        if let Some(&(derived, class)) = payload.get(&idx) {
            by_contig.entry(f[0].to_string()).or_default().push((end, derived, class));
        }
    }
    anyhow::ensure!(!by_contig.is_empty(), "no classification sites survived the join");

    let mut contigs = Vec::with_capacity(by_contig.len());
    let (mut total, mut nea, mut den) = (0usize, 0usize, 0usize);
    for (contig, mut rows) in by_contig {
        rows.sort_by_key(|r| r.0);
        rows.dedup_by_key(|r| r.0);
        let positions: Vec<i64> = rows.iter().map(|r| r.0).collect();
        let derived: Vec<u8> = rows.iter().map(|r| r.1 as u8).collect();
        let classes: Vec<u8> = rows
            .iter()
            .map(|r| match r.2 {
                DiagnosticClass::Neanderthal => 0u8,
                DiagnosticClass::Denisovan => 1,
                DiagnosticClass::SharedArchaic => 2,
            })
            .collect();
        total += rows.len();
        nea += classes.iter().filter(|c| **c == 0).count();
        den += classes.iter().filter(|c| **c == 1).count();
        contigs.push(ClassifyContig {
            positions: PositionStream::encode(&contig, &positions),
            derived,
            classes,
        });
    }

    let asset = ArchaicClassify {
        build: BUILD.to_string(),
        contigs,
    };
    let encoded = asset.to_bytes().map_err(|e| anyhow::anyhow!("{e}"))?;
    eprintln!(
        "classify: {total} sites over {} contigs ({nea} Neanderthal, {den} Denisovan, {} shared) -> {:.1} MB",
        asset.contigs.len(),
        total - nea - den,
        encoded.len() as f64 / 1_048_576.0
    );
    write_bin(&args.out, &encoded)?;
    eprintln!("wrote {}", args.out.display());
    Ok(())
}

#[derive(Parser)]
pub struct ArchaicCallableArgs {
    /// Callable regions as BED on the target build — the intersection of the four archaic genomes'
    /// `FilterBed` masks, lifted. Where all four archaic genomes are callable is where a
    /// private-variant excess can be interpreted at all.
    #[arg(long)]
    bed: PathBuf,
    /// Window width the counts are binned at; must match the segment caller's `window_bp`.
    #[arg(long, default_value_t = 1000)]
    window_bp: i64,
    /// Output (bincode `ArchaicCallable`).
    #[arg(long)]
    out: PathBuf,
}

pub fn build_archaic_callable(args: ArchaicCallableArgs) -> Result<()> {
    anyhow::ensure!(args.window_bp > 0, "--window-bp must be positive");
    let rdr = open_maybe_gz(&args.bed).with_context(|| format!("opening {}", args.bed.display()))?;
    // Accumulate callable bases per window, per contig.
    let mut spans: BTreeMap<String, Vec<(i64, i64)>> = BTreeMap::new();
    for line in rdr.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 3 {
            continue;
        }
        let (Ok(s), Ok(e)) = (f[1].parse::<i64>(), f[2].parse::<i64>()) else {
            continue;
        };
        if e > s {
            spans.entry(f[0].to_string()).or_default().push((s, e));
        }
    }
    anyhow::ensure!(!spans.is_empty(), "no callable intervals read from {}", args.bed.display());

    let mut contigs = Vec::with_capacity(spans.len());
    let mut total_bp = 0f64;
    for (contig, mut v) in spans {
        v.sort_unstable();
        let start = (v[0].0 / args.window_bp) * args.window_bp;
        let last = v.iter().map(|(_, e)| *e).max().unwrap_or(start);
        let n = (((last - start) / args.window_bp) + 1) as usize;
        let mut callable_bp = vec![0u16; n];
        for (s, e) in v {
            // Split the interval across the windows it touches, saturating per window.
            let (mut cur, end) = (s, e);
            while cur < end {
                let idx = ((cur - start) / args.window_bp) as usize;
                let win_end = start + (idx as i64 + 1) * args.window_bp;
                let take = end.min(win_end) - cur;
                if let Some(slot) = callable_bp.get_mut(idx) {
                    *slot = slot.saturating_add(take.clamp(0, u16::MAX as i64) as u16).min(args.window_bp as u16);
                }
                cur = win_end;
            }
        }
        total_bp += callable_bp.iter().map(|&b| b as f64).sum::<f64>();
        contigs.push(CallableContig {
            contig,
            start,
            callable_bp,
        });
    }

    let asset = ArchaicCallable {
        build: BUILD.to_string(),
        window_bp: args.window_bp,
        contigs,
    };
    let encoded = asset.to_bytes().map_err(|e| anyhow::anyhow!("{e}"))?;
    eprintln!(
        "callable: {:.1} Mb over {} contigs, {} bp windows -> {:.1} MB",
        total_bp / 1_000_000.0,
        asset.contigs.len(),
        args.window_bp,
        encoded.len() as f64 / 1_048_576.0
    );
    write_bin(&args.out, &encoded)?;
    eprintln!("wrote {}", args.out.display());
    Ok(())
}
