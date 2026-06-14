//! Build the chip-compatible IBD panel asset (`ibd_panel_<build>.bin`) — a multi-build,
//! palindrome-free SNP set (ancestry-ibd-asset-wiring B2).
//!
//! Input is a tab-separated table with a **named header**; required columns are the rsID and the
//! CHM13 locus; the GRCh37/GRCh38 loci are optional per row (blank ⇒ that build absent for the
//! site). The multi-build coordinates must come from an **allele-aware** liftover (GATK
//! `LiftoverVcf`, which reverse-complements + swaps REF/ALT on inverted chain blocks) so each
//! build's `(REF, ALT)` are the same biological alleles — a CrossMap lift would silently corrupt
//! ~3/4 of sites. This step parses that table, drops strand-ambiguous palindromes, and serializes.
//!
//! ```text
//! rsid  chm13_contig chm13_pos chm13_ref chm13_alt  grch37_contig grch37_pos grch37_ref grch37_alt  grch38_contig grch38_pos grch38_ref grch38_alt
//! ```

use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Parser;
use navigator_analysis::ibd_panel::{IbdPanel, IbdPanelSite, Locus};

#[derive(Parser)]
pub struct IbdPanelArgs {
    /// Multi-build sites table (TSV with a named header; see module docs). CHM13 columns required.
    #[arg(long)]
    pub input: PathBuf,
    /// Output asset (bincode), e.g. `~/.decodingus/ancestry/ibd_panel_chm13v2.0.bin`.
    #[arg(long)]
    pub out: PathBuf,
    /// Canonical build label for the CHM13 loci.
    #[arg(long, default_value = "chm13v2.0")]
    pub build: String,
}

/// One ACGT base, or `None` (multiallelic / indel / blank).
fn base(s: &str) -> Option<char> {
    let c = s.trim().chars().next()?;
    matches!(c.to_ascii_uppercase(), 'A' | 'C' | 'G' | 'T').then_some(c.to_ascii_uppercase()).filter(|_| s.trim().len() == 1)
}

fn locus(cols: &HashMap<String, usize>, f: &[&str], prefix: &str) -> Option<Locus> {
    let get = |name: &str| cols.get(name).and_then(|&i| f.get(i)).map(|s| s.trim()).filter(|s| !s.is_empty());
    let contig = get(&format!("{prefix}_contig"))?;
    let position = get(&format!("{prefix}_pos"))?.parse::<i64>().ok()?;
    let reference = base(get(&format!("{prefix}_ref"))?)?;
    let alternate = base(get(&format!("{prefix}_alt"))?)?;
    Some(Locus { contig: contig.to_string(), position, reference, alternate })
}

pub fn build_ibd_panel(args: IbdPanelArgs) -> Result<()> {
    let file = File::open(&args.input).with_context(|| format!("open {}", args.input.display()))?;
    let mut lines = BufReader::new(file).lines();
    let header = lines
        .next()
        .transpose()?
        .ok_or_else(|| anyhow::anyhow!("empty input {}", args.input.display()))?;
    let cols: HashMap<String, usize> =
        header.split('\t').enumerate().map(|(i, h)| (h.trim().to_ascii_lowercase(), i)).collect();
    for required in ["rsid", "chm13_contig", "chm13_pos", "chm13_ref", "chm13_alt"] {
        anyhow::ensure!(cols.contains_key(required), "missing required column `{required}`");
    }

    let mut sites = Vec::new();
    let mut skipped = 0usize;
    for line in lines {
        let line = line?;
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        let rsid = cols.get("rsid").and_then(|&i| f.get(i)).map(|s| s.trim()).unwrap_or("");
        match locus(&cols, &f, "chm13") {
            Some(chm13) if !rsid.is_empty() => sites.push(IbdPanelSite {
                rsid: rsid.to_string(),
                chm13,
                grch37: locus(&cols, &f, "grch37"),
                grch38: locus(&cols, &f, "grch38"),
            }),
            _ => skipped += 1, // missing/indel/multiallelic CHM13 locus or rsid
        }
    }
    anyhow::ensure!(!sites.is_empty(), "no usable sites parsed from {}", args.input.display());

    let (panel, palindromes) = IbdPanel::from_sites(args.build, sites);
    let bytes = panel.to_bytes().map_err(|e| anyhow::anyhow!(e.to_string()))?;
    if let Some(parent) = args.out.parent() {
        fs::create_dir_all(parent).ok();
    }
    File::create(&args.out)?.write_all(&bytes)?;
    let with37 = panel.sites.iter().filter(|s| s.grch37.is_some()).count();
    let with38 = panel.sites.iter().filter(|s| s.grch38.is_some()).count();
    eprintln!(
        "wrote {} ({} sites; {with37} with GRCh37, {with38} with GRCh38; {palindromes} palindromes dropped, {skipped} rows skipped)",
        args.out.display(),
        panel.sites.len()
    );
    Ok(())
}
