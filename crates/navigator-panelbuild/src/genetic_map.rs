//! Build the IBD genetic-map asset (`genetic_map_<build>.bin`) from a recombination map.
//!
//! IBD segment lengths in cM (and thus the relationship bands) are only as good as the map, so this
//! replaces the app's flat 1 cM/Mb stand-in with a real sex-averaged map (deCODE 2019 / HapMap II),
//! **already lifted to CHM13**. The lift is coordinate-only (no alleles), so a stage-2 CrossMap BED
//! lift of the map's positions is sufficient — this step just parses the lifted text and serializes
//! it to the bincode [`navigator_analysis::ibd::GeneticMap`] the app loads.
//!
//! Input is whitespace/tab-delimited with columns `chromosome  position(bp)  …  cumulative_cM`
//! (the **last** column is the cumulative genetic position — matching HapMap's
//! `Chromosome Position(bp) Rate(cM/Mb) Map(cM)` and the simple `chrom pos cM` form). A non-numeric
//! first data row is treated as a header and skipped. Positions must be CHM13 coordinates.

use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Parser;
use navigator_analysis::ibd::GeneticMap;

#[derive(Parser)]
pub struct GeneticMapArgs {
    /// Recombination-map text (CHM13 coordinates): `chromosome  position(bp)  …  cumulative_cM`.
    #[arg(long)]
    pub input: PathBuf,
    /// Output asset (bincode), e.g. `~/.decodingus/ancestry/genetic_map_chm13v2.0.bin`.
    #[arg(long)]
    pub out: PathBuf,
}

pub fn build_genetic_map(args: GeneticMapArgs) -> Result<()> {
    let file = File::open(&args.input).with_context(|| format!("open {}", args.input.display()))?;
    let mut by_chrom: BTreeMap<String, Vec<(i32, f64)>> = BTreeMap::new();
    let mut parsed = 0usize;
    for line in BufReader::new(file).lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split_whitespace().collect();
        if f.len() < 3 {
            continue;
        }
        // Position is col 1; cumulative cM is the last column. A header row (non-numeric) is skipped.
        let (Ok(pos), Ok(cm)) = (f[1].parse::<i32>(), f[f.len() - 1].parse::<f64>()) else {
            continue;
        };
        by_chrom.entry(f[0].to_string()).or_default().push((pos, cm));
        parsed += 1;
    }
    if parsed == 0 {
        anyhow::bail!("no (chrom, pos, cM) rows parsed from {}", args.input.display());
    }
    let n_chrom = by_chrom.len();

    let markers = by_chrom.into_iter().map(|(chrom, mut rows)| {
        rows.sort_by_key(|(p, _)| *p);
        rows.dedup_by_key(|(p, _)| *p); // collapse duplicate positions (keep first after sort)
        let positions = rows.iter().map(|(p, _)| *p).collect::<Vec<_>>();
        let cm = rows.iter().map(|(_, c)| *c).collect::<Vec<_>>();
        (chrom, positions, cm)
    });
    let map = GeneticMap::from_markers(markers);
    let bytes = map.to_bytes().map_err(|e| anyhow::anyhow!(e.to_string()))?;
    if let Some(parent) = args.out.parent() {
        fs::create_dir_all(parent).ok();
    }
    File::create(&args.out)?.write_all(&bytes)?;
    eprintln!(
        "wrote {} ({parsed} markers across {n_chrom} chromosomes)",
        args.out.display()
    );
    Ok(())
}
