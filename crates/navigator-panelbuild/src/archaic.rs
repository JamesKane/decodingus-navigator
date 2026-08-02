//! Build the archaic-informative marker panel (`archaic_markers_<build>.bin`) — the Tier A asset
//! behind the Neanderthal / Denisovan report. Design: `documents/design/ArchaicAncestry_Design.md`
//! §4 (assets) and §10 (M1).
//!
//! **Why this is two subcommands.** Polarity (which allele is ancestral) must be assigned in
//! **GRCh37**, because that is the build both the EVA archaic VCFs and the Ensembl-75 EPO ancestral
//! sequence are in — doing it there costs no liftover and no allele bookkeeping (design §3a). The
//! panel, however, must ship in CHM13 coordinates. So the work straddles the pipeline's existing
//! lift stage:
//!
//! 1. `archaic-candidates` (GRCh37) — read the per-genome archaic genotype tables, assign
//!    ancestral/derived from the EPO sequence, keep sites with an archaic homozygous-derived call,
//!    and emit a candidates table plus a BED for the lift.
//! 2. *(shell)* `CrossMap bed` GRCh37 → CHM13, as `02_liftover_panel_sites.sh` does.
//! 3. `archaic-panel` (CHM13) — re-join the lifted coordinates, **orient against the CHM13 FASTA**,
//!    apply the African-outgroup filter, drop palindromes, classify, and write the asset.
//!
//! Step 3's orientation is load-bearing: `CrossMap bed` is *not* allele-aware, so a large share of
//! sites arrive with ref/alt reversed relative to CHM13. Shipping them unoriented is precisely the
//! defect that forced an asset rebuild in the ancient-ancestry work (that design's §7.16), and it is
//! silent — f4-style statistics are invariant under a consistent flip, so nothing fails loudly.
//!
//! The heavy lifting (VCF decode, mask intersection, allele-frequency extraction) stays in
//! `bcftools` inside `08_build_archaic.sh`, matching how every other stage feeds this crate: the
//! Rust side consumes tab-separated tables so the selection logic is pure and unit-testable.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use clap::Parser;
use navigator_analysis::archaic::{
    classify_diagnostic, ArchaicCall, ArchaicMarkerPanel, ArchaicPanelThresholds, ArchaicSite, ARCHAIC_GENOMES,
};
use navigator_analysis::ibd_panel::{is_palindromic, Locus};

use crate::pca::{open_maybe_gz, write_bin};

/// The reference build the emitted panel is in (matches the other assets).
const BUILD: &str = "chm13v2.0";

// ─────────────────────────────── stage 1: candidates (GRCh37) ───────────────────────────────

#[derive(Parser)]
pub struct ArchaicCandidatesArgs {
    /// Comma-separated per-genome genotype tables, in [`ARCHAIC_GENOMES`] order (Altai, Vindija,
    /// Chagyrskaya, Denisova). Each is `CHROM POS REF ALT GT`, already restricted to biallelic SNVs
    /// passing that genome's own `FilterBed/` quality mask (bcftools does this in stage 08).
    #[arg(long)]
    archaic: String,
    /// Directory of Ensembl release-75 EPO ancestral FASTAs, one per chromosome
    /// (`homo_sapiens_ancestor_<chr>.fa`). GRCh37 — the same build as the archaic VCFs, which is
    /// why polarity is assigned here and not after the lift.
    #[arg(long)]
    ancestral: PathBuf,
    /// Candidates table (GRCh37) consumed by `archaic-panel` after the lift.
    #[arg(long)]
    out: PathBuf,
    /// BED (GRCh37) to hand to `CrossMap bed`. Its name column carries the candidate's row index, so
    /// the lift can drop sites without disturbing the join.
    #[arg(long)]
    out_bed: PathBuf,
    /// Minimum number of archaic genomes with a call at a site.
    #[arg(long, default_value_t = 1)]
    min_archaic_called: usize,
}

/// A candidate site in GRCh37, before the lift.
#[derive(Debug, Clone, PartialEq)]
struct Candidate {
    contig: String,
    position: i64,
    reference_allele: char,
    alternate_allele: char,
    derived: char,
    calls: [ArchaicCall; 4],
}

/// Split a comma-separated path list (mirrors `hap_panel`'s handling of multi-source inputs).
fn split_list(s: &str) -> Vec<PathBuf> {
    s.split(',')
        .map(|p| PathBuf::from(p.trim()))
        .filter(|p| !p.as_os_str().is_empty())
        .collect()
}

/// The derived allele at a biallelic site, given the inferred ancestral base.
///
/// Returns `None` when the ancestral base is unusable — either low-confidence/absent (the EPO
/// sequence lower-cases low-confidence calls and uses `.`/`-`/`N` for gaps) or matching neither
/// allele, which means the site cannot be polarized and must be dropped rather than guessed.
///
/// Only **upper-case** ancestral bases are accepted: in the EPO alignment lower case marks a
/// low-confidence call, and polarity is the one thing this panel cannot afford to get wrong — an
/// inverted site turns an archaic-derived allele into its opposite.
fn derived_allele(ancestral: char, reference_allele: char, alternate_allele: char) -> Option<char> {
    if !matches!(ancestral, 'A' | 'C' | 'G' | 'T') {
        return None;
    }
    let (r, a) = (
        reference_allele.to_ascii_uppercase(),
        alternate_allele.to_ascii_uppercase(),
    );
    match (ancestral == r, ancestral == a) {
        (true, false) => Some(a),
        (false, true) => Some(r),
        _ => None,
    }
}

/// Parse a VCF `GT` field into its two alleles as bases. Handles `/` and `|`, and treats any
/// missing or non-0/1 allele index as missing (the inputs are pre-restricted to biallelic SNVs).
///
/// **Ploidy is read from the field's arity, not inferred.** A single-allele record (`1`) is a
/// genuinely haploid call and reads as homozygous, but a two-slot record with one slot missing
/// (`0/.`) must stay half-missing. Collapsing the latter to homozygous would let a partially-called
/// archaic genome fabricate the `HomDerived` state that gates site selection.
fn parse_gt(gt: &str, reference_allele: char, alternate_allele: char) -> (Option<char>, Option<char>) {
    let alleles: Vec<Option<char>> = gt
        .split(['/', '|'])
        .map(|a| match a.trim() {
            "0" => Some(reference_allele.to_ascii_uppercase()),
            "1" => Some(alternate_allele.to_ascii_uppercase()),
            _ => None,
        })
        .collect();
    match alleles.as_slice() {
        [a] => (*a, *a),
        [a, b, ..] => (*a, *b),
        [] => (None, None),
    }
}

/// One genome's state relative to `derived`.
fn call_state(alleles: (Option<char>, Option<char>), derived: char) -> ArchaicCall {
    match alleles {
        (Some(a), Some(b)) => match u8::from(a == derived) + u8::from(b == derived) {
            2 => ArchaicCall::HomDerived,
            1 => ArchaicCall::Het,
            _ => ArchaicCall::HomAncestral,
        },
        _ => ArchaicCall::NoCall,
    }
}

/// `(contig, pos)` → `(ref, alt, gt)` for one archaic genome. `alt` is `None` for a
/// **reference-confident** record (`ALT=.`), which the EVA all-sites VCFs emit wherever the genome
/// matches hg19.
type GenomeSites = HashMap<(String, i64), (char, Option<char>, String)>;

/// Read one genome's genotype table.
///
/// The EVA archaic VCFs are **all-sites** files, so a record with `ALT=.` is not junk: it is a
/// positive statement that the genome is homozygous for the hg19 reference base. Keeping those is
/// load-bearing in two ways — it is the only way `HomAncestral` can ever be distinguished from
/// "masked out", and, where the EPO sequence says the *reference* allele is the derived one, a
/// hom-ref archaic genome IS homozygous-derived, i.e. exactly the donor state this panel selects on.
/// Dropping them would silently bias the panel against sites where hg19 carries the archaic allele.
fn load_archaic_table(path: &Path) -> Result<GenomeSites> {
    let rdr = open_maybe_gz(path).with_context(|| format!("opening {}", path.display()))?;
    let mut out = GenomeSites::new();
    for line in rdr.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 5 {
            continue;
        }
        let (Ok(pos), Some(r)) = (f[1].parse::<i64>(), f[2].chars().next()) else {
            continue;
        };
        // A multi-character REF means an indel slipped past the SNV filter.
        if f[2].len() != 1 {
            continue;
        }
        let alt = match f[3] {
            "." => None,
            a if a.len() == 1 => Some(a.chars().next().unwrap_or('N')),
            // Multi-character or multi-allelic ALT: not a biallelic SNV.
            _ => continue,
        };
        out.insert((normalize_contig(f[0]), pos), (r, alt, f[4].to_string()));
    }
    Ok(out)
}

/// Strip a `chr` prefix so GRCh37 sources that disagree about it still join.
fn normalize_contig(c: &str) -> String {
    c.strip_prefix("chr").unwrap_or(c).to_string()
}

/// Read the single record of an Ensembl EPO ancestral FASTA for one chromosome.
fn read_ancestral(dir: &Path, chrom: &str) -> Result<Vec<u8>> {
    let path = dir.join(format!("homo_sapiens_ancestor_{chrom}.fa"));
    let rdr = open_maybe_gz(&path).with_context(|| format!("opening ancestral FASTA {}", path.display()))?;
    let mut seq = Vec::new();
    for line in rdr.lines() {
        let line = line?;
        if line.starts_with('>') {
            continue;
        }
        seq.extend_from_slice(line.trim_end().as_bytes());
    }
    Ok(seq)
}

pub fn build_archaic_candidates(args: ArchaicCandidatesArgs) -> Result<()> {
    let paths = split_list(&args.archaic);
    anyhow::ensure!(
        paths.len() == ARCHAIC_GENOMES.len(),
        "--archaic needs {} comma-separated tables in {:?} order, got {}",
        ARCHAIC_GENOMES.len(),
        ARCHAIC_GENOMES,
        paths.len()
    );

    let genomes: Vec<GenomeSites> = paths
        .iter()
        .zip(ARCHAIC_GENOMES)
        .map(|(p, name)| {
            let sites = load_archaic_table(p)?;
            eprintln!("{name}: {} biallelic SNV sites", sites.len());
            Ok(sites)
        })
        .collect::<Result<_>>()?;

    // Union of every genome's callable sites, keyed by position. A strict four-way intersection
    // would discard most of the panel — the four `FilterBed/` masks differ substantially, and the
    // hom-derived requirement below is what actually defines the introgression donor state.
    let mut keys: Vec<(String, i64)> = genomes.iter().flat_map(|g| g.keys().cloned()).collect();
    keys.sort();
    keys.dedup();
    eprintln!("{} candidate positions across the four genomes", keys.len());

    let mut candidates: Vec<Candidate> = Vec::new();
    let mut cur_chrom = String::new();
    let mut ancestral_seq: Vec<u8> = Vec::new();
    let (mut no_polarity, mut no_hom_derived, mut allele_conflict) = (0usize, 0usize, 0usize);

    for (contig, position) in keys {
        if contig != cur_chrom {
            ancestral_seq = read_ancestral(&args.ancestral, &contig)?;
            cur_chrom = contig.clone();
        }
        let present: Vec<&(char, Option<char>, String)> = genomes
            .iter()
            .filter_map(|g| g.get(&(contig.clone(), position)))
            .collect();
        if present.len() < args.min_archaic_called {
            continue;
        }
        // The site's allele pair comes from the genomes that actually carry a variant; a
        // reference-confident record states only the REF base and cannot define the pair. At least
        // one genome must vary, otherwise the site is invariant across all four and carries no
        // information regardless of polarity.
        let Some((reference_allele, alternate_allele)) = present.iter().find_map(|(r, a, _)| a.map(|alt| (*r, alt)))
        else {
            continue;
        };
        // The varying genomes must agree on that pair, and no genome may contradict the REF base.
        if present
            .iter()
            .any(|(r, a, _)| *r != reference_allele || a.is_some_and(|alt| alt != alternate_allele))
        {
            allele_conflict += 1;
            continue;
        }
        let ancestral = ancestral_seq
            .get((position - 1) as usize)
            .map(|b| *b as char)
            .unwrap_or('N');
        let Some(derived) = derived_allele(ancestral, reference_allele, alternate_allele) else {
            no_polarity += 1;
            continue;
        };

        let mut calls = [ArchaicCall::NoCall; 4];
        for (i, g) in genomes.iter().enumerate() {
            if let Some((r, a, gt)) = g.get(&(contig.clone(), position)) {
                // A reference-confident record (`ALT=.`) means both alleles are the REF base — which
                // is `HomDerived` when the reference itself carries the derived allele.
                calls[i] = call_state(parse_gt(gt, *r, a.unwrap_or(*r)), derived);
            }
        }
        // The donor state: at least one archaic genome homozygous for the derived allele (design §4).
        if !calls.contains(&ArchaicCall::HomDerived) {
            no_hom_derived += 1;
            continue;
        }
        candidates.push(Candidate {
            contig: contig.clone(),
            position,
            reference_allele: reference_allele.to_ascii_uppercase(),
            alternate_allele: alternate_allele.to_ascii_uppercase(),
            derived,
            calls,
        });
    }

    eprintln!(
        "{} candidates kept ({no_polarity} unpolarizable, {no_hom_derived} no hom-derived archaic, \
         {allele_conflict} allele-conflicting)",
        candidates.len()
    );
    anyhow::ensure!(!candidates.is_empty(), "no candidate site survived selection");

    write_candidates(&args.out, &candidates)?;
    write_bed(&args.out_bed, &candidates)?;
    eprintln!(
        "wrote {} and {} — lift the BED to {BUILD}, then run `archaic-panel`",
        args.out.display(),
        args.out_bed.display()
    );
    Ok(())
}

fn call_token(c: ArchaicCall) -> char {
    match c {
        ArchaicCall::HomAncestral => '0',
        ArchaicCall::Het => '1',
        ArchaicCall::HomDerived => '2',
        ArchaicCall::NoCall => '.',
    }
}

fn parse_call_token(c: char) -> ArchaicCall {
    match c {
        '0' => ArchaicCall::HomAncestral,
        '1' => ArchaicCall::Het,
        '2' => ArchaicCall::HomDerived,
        _ => ArchaicCall::NoCall,
    }
}

fn write_candidates(path: &Path, candidates: &[Candidate]) -> Result<()> {
    let mut w = BufWriter::new(File::create(path).with_context(|| format!("creating {}", path.display()))?);
    writeln!(w, "#idx\tcontig\tpos\tref\talt\tderived\tcalls")?;
    for (i, c) in candidates.iter().enumerate() {
        let calls: String = c.calls.iter().map(|c| call_token(*c)).collect();
        writeln!(
            w,
            "{i}\t{}\t{}\t{}\t{}\t{}\t{calls}",
            c.contig, c.position, c.reference_allele, c.alternate_allele, c.derived
        )?;
    }
    w.flush()?;
    Ok(())
}

/// 0-based BED whose name column is the candidate index, so `CrossMap` can drop rows without
/// breaking the join back to the payload.
fn write_bed(path: &Path, candidates: &[Candidate]) -> Result<()> {
    let mut w = BufWriter::new(File::create(path).with_context(|| format!("creating {}", path.display()))?);
    for (i, c) in candidates.iter().enumerate() {
        writeln!(w, "chr{}\t{}\t{}\t{i}", c.contig, c.position - 1, c.position)?;
    }
    w.flush()?;
    Ok(())
}

// ─────────────────────────────── stage 3: panel (CHM13) ───────────────────────────────

#[derive(Parser)]
pub struct ArchaicPanelArgs {
    /// The candidates table written by `archaic-candidates` (GRCh37 payload).
    #[arg(long)]
    candidates: PathBuf,
    /// The lifted BED (CHM13) produced by `CrossMap bed` from `--out-bed`; its name column carries
    /// the candidate index.
    #[arg(long)]
    lifted: PathBuf,
    /// African-outgroup allele frequencies on CHM13: `CHROM POS REF ALT AF_AFR AF_NONAFR`, derived
    /// from the per-super-pop `AC_<POP>_unrel`/`AN_<POP>_unrel` INFO fields the 1000G-on-CHM13 VCFs
    /// already carry (design §9 Q2 — no new data source).
    #[arg(long)]
    outgroup_af: PathBuf,
    /// CHM13 reference FASTA (indexed). Each surviving site is oriented so `reference_allele` is the
    /// actual CHM13 base; sites matching neither allele are dropped. **Not optional** — `CrossMap
    /// bed` is not allele-aware.
    #[arg(long)]
    reference: PathBuf,
    /// Maximum derived-allele frequency in the African outgroup.
    ///
    /// 0.01 confirmed by the checkpoint-A sweep: relaxing it costs precision without buying recall
    /// (0.02 -> 73.8 %, 0.05 -> 63.9 %, against 78.4 % here). This is the criterion doing the actual
    /// archaic-specificity work.
    #[arg(long, default_value_t = 0.01)]
    max_afr_freq: f32,
    /// Minimum derived-allele frequency outside Africa.
    ///
    /// 0.0005, calibrated against the hmmix callset (design §10 checkpoint A): F1 peaks there at
    /// 0.701 (precision 78.4 %, recall 63.4 %). The floor is small but load-bearing — removing it
    /// entirely (0.0) collapses precision to 64.7 % for almost no recall gain. Most confidently
    /// introgressed variants are rare (hmmix median frequency 0.87 %), so a high floor mostly
    /// discards signal: the original 0.05 scored F1 0.473 at 33.8 % recall.
    #[arg(long, default_value_t = 0.0005)]
    min_non_afr_freq: f32,
    /// Output panel (bincode `ArchaicMarkerPanel`).
    #[arg(long)]
    out: PathBuf,
    /// Optional lifted BED (GRCh38) from the same candidates, so the panel carries hg38 coordinates
    /// and a GRCh38 alignment can be genotyped without a runtime liftover.
    #[arg(long)]
    lifted_hg38: Option<PathBuf>,
    /// GRCh38 reference FASTA (indexed). Required with `--lifted-hg38`: the hg38 lift is not
    /// allele-aware either, so its alleles must be oriented the same way the CHM13 ones are.
    #[arg(long)]
    reference_hg38: Option<PathBuf>,
    /// Optional inspection TSV of the final sites.
    #[arg(long)]
    sites_tsv: Option<PathBuf>,
}

fn load_candidates(path: &Path) -> Result<HashMap<usize, Candidate>> {
    let rdr = open_maybe_gz(path).with_context(|| format!("opening {}", path.display()))?;
    let mut out = HashMap::new();
    for line in rdr.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 7 {
            continue;
        }
        let idx: usize = f[0].parse()?;
        let mut calls = [ArchaicCall::NoCall; 4];
        for (i, c) in f[6].chars().take(4).enumerate() {
            calls[i] = parse_call_token(c);
        }
        out.insert(
            idx,
            Candidate {
                contig: f[1].to_string(),
                position: f[2].parse()?,
                reference_allele: f[3].chars().next().unwrap_or('N'),
                alternate_allele: f[4].chars().next().unwrap_or('N'),
                derived: f[5].chars().next().unwrap_or('N'),
                calls,
            },
        );
    }
    Ok(out)
}

/// Lifted coordinates keyed by candidate index: `(contig, 1-based pos)`.
fn load_lifted(path: &Path) -> Result<HashMap<usize, (String, i64)>> {
    let rdr = open_maybe_gz(path).with_context(|| format!("opening {}", path.display()))?;
    let mut out = HashMap::new();
    for line in rdr.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 4 {
            continue;
        }
        let (Ok(end), Ok(idx)) = (f[2].parse::<i64>(), f[3].trim().parse::<usize>()) else {
            continue;
        };
        out.insert(idx, (f[0].to_string(), end));
    }
    Ok(out)
}

/// `(contig, pos)` → `(ref, alt, af_afr, af_non_afr)` on CHM13.
type OutgroupAf = HashMap<(String, i64), (char, char, f32, f32)>;

fn load_outgroup_af(path: &Path) -> Result<OutgroupAf> {
    let rdr = open_maybe_gz(path).with_context(|| format!("opening {}", path.display()))?;
    let mut out = OutgroupAf::new();
    for line in rdr.lines() {
        let line = line?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 6 {
            continue;
        }
        let (Ok(pos), Some(r), Some(a), Ok(afr), Ok(non_afr)) = (
            f[1].parse::<i64>(),
            f[2].chars().next(),
            f[3].chars().next(),
            f[4].parse::<f32>(),
            f[5].parse::<f32>(),
        ) else {
            continue;
        };
        out.insert((f[0].to_string(), pos), (r, a, afr, non_afr));
    }
    Ok(out)
}

/// The derived allele's frequency, given frequencies stated for the outgroup table's own ALT.
///
/// The outgroup table and the candidate need not label ref/alt the same way, so the frequency is
/// re-expressed against the *derived base* rather than against whichever allele happened to be
/// called ALT. Returns `None` when the candidate's derived base is absent from the outgroup's
/// allele pair, which means the two sources disagree about the site and it cannot be filtered.
fn derived_freq(derived: char, og_ref: char, og_alt: char, af_alt: f32) -> Option<f32> {
    let d = derived.to_ascii_uppercase();
    if d == og_alt.to_ascii_uppercase() {
        Some(af_alt)
    } else if d == og_ref.to_ascii_uppercase() {
        Some(1.0 - af_alt)
    } else {
        None
    }
}

/// Lift + orient the candidates onto GRCh38, returning `candidate index -> Locus`.
///
/// Same discipline as the CHM13 pass: `CrossMap bed` is not allele-aware, so each site is oriented
/// against the hg38 reference base (swap ref/alt where reversed, drop where neither matches).
fn build_hg38_loci(
    bed: &Path,
    reference: &Path,
    candidates: &HashMap<usize, Candidate>,
) -> Result<HashMap<usize, Locus>> {
    let lifted = load_lifted(bed)?;
    let mut rows: Vec<(usize, String, i64)> = lifted.into_iter().map(|(i, (c, p))| (i, c, p)).collect();
    rows.sort_by(|a, b| a.1.cmp(&b.1).then(a.2.cmp(&b.2)));
    let (mut out, mut cur, mut seq, mut dropped) = (HashMap::new(), String::new(), Vec::new(), 0usize);
    for (idx, contig, position) in rows {
        let Some(cand) = candidates.get(&idx) else { continue };
        if contig != cur {
            seq = navigator_analysis::reader::read_contig_sequence(reference, &contig)
                .map_err(|e| anyhow::anyhow!("reading {contig} from {}: {e}", reference.display()))?;
            cur = contig.clone();
        }
        let base = seq
            .get((position - 1) as usize)
            .map(|b| b.to_ascii_uppercase() as char)
            .unwrap_or('N');
        let (reference_allele, alternate_allele) = if base == cand.reference_allele {
            (cand.reference_allele, cand.alternate_allele)
        } else if base == cand.alternate_allele {
            (cand.alternate_allele, cand.reference_allele)
        } else {
            dropped += 1;
            continue;
        };
        out.insert(
            idx,
            Locus {
                contig,
                position,
                reference: reference_allele,
                alternate: alternate_allele,
            },
        );
    }
    if dropped > 0 {
        eprintln!("GRCh38: {dropped} sites dropped (reference base matched neither allele)");
    }
    Ok(out)
}

pub fn build_archaic_panel(args: ArchaicPanelArgs) -> Result<()> {
    let candidates = load_candidates(&args.candidates)?;
    let lifted = load_lifted(&args.lifted)?;
    let outgroup = load_outgroup_af(&args.outgroup_af)?;
    eprintln!(
        "{} candidates, {} lifted to {BUILD}, {} outgroup AF rows",
        candidates.len(),
        lifted.len(),
        outgroup.len()
    );

    // GRCh38 loci: lift + orient exactly as the CHM13 pass does. Optional — without them a GRCh38
    // alignment simply falls back to whatever the consensus covers.
    let hg38 = match (&args.lifted_hg38, &args.reference_hg38) {
        (Some(bed), Some(fa)) => build_hg38_loci(bed, fa, &candidates)?,
        (Some(_), None) => anyhow::bail!("--lifted-hg38 requires --reference-hg38 (the lift is not allele-aware)"),
        _ => HashMap::new(),
    };
    if !hg38.is_empty() {
        eprintln!("{} sites carry a GRCh38 locus", hg38.len());
    }

    // Walk in lifted-coordinate order so each contig's reference sequence is loaded once.
    let mut rows: Vec<(usize, String, i64)> = lifted.iter().map(|(i, (c, p))| (*i, c.clone(), *p)).collect();
    rows.sort_by(|a, b| a.1.cmp(&b.1).then(a.2.cmp(&b.2)));

    let (mut no_af, mut afr_common, mut too_rare, mut palindromic, mut flipped, mut unoriented) =
        (0usize, 0usize, 0usize, 0usize, 0usize, 0usize);
    let mut sites: Vec<ArchaicSite> = Vec::new();
    let mut cur = String::new();
    let mut seq: Vec<u8> = Vec::new();

    for (idx, contig, position) in rows {
        let Some(cand) = candidates.get(&idx) else { continue };

        // Outgroup filter — the introgression signature (design §4 step 3).
        let Some((og_ref, og_alt, af_alt, non_afr_alt)) = outgroup.get(&(contig.clone(), position)) else {
            no_af += 1;
            continue;
        };
        let (Some(afr), Some(non_afr)) = (
            derived_freq(cand.derived, *og_ref, *og_alt, *af_alt),
            derived_freq(cand.derived, *og_ref, *og_alt, *non_afr_alt),
        ) else {
            no_af += 1;
            continue;
        };
        if afr > args.max_afr_freq {
            afr_common += 1;
            continue;
        }
        if non_afr < args.min_non_afr_freq {
            too_rare += 1;
            continue;
        }

        // Strand-ambiguous sites cannot be reconciled against a chip's unknown strand.
        if is_palindromic(cand.reference_allele, cand.alternate_allele) {
            palindromic += 1;
            continue;
        }

        // Orient against CHM13 (see the module docs — CrossMap is not allele-aware).
        if contig != cur {
            seq = navigator_analysis::reader::read_contig_sequence(&args.reference, &contig)
                .map_err(|e| anyhow::anyhow!("reading {contig} from {}: {e}", args.reference.display()))?;
            cur = contig.clone();
        }
        let base = seq
            .get((position - 1) as usize)
            .map(|b| b.to_ascii_uppercase() as char)
            .unwrap_or('N');
        let (reference_allele, alternate_allele) = if base == cand.reference_allele {
            (cand.reference_allele, cand.alternate_allele)
        } else if base == cand.alternate_allele {
            flipped += 1;
            (cand.alternate_allele, cand.reference_allele)
        } else {
            unoriented += 1;
            continue;
        };

        sites.push(ArchaicSite {
            contig: contig.clone(),
            position,
            reference_allele,
            alternate_allele,
            // GRCh37 is EXACT — the archaic VCFs' own coordinates and alleles, never lifted, so this
            // build carries no liftover or strand risk at all.
            grch37: Some(Locus {
                contig: cand.contig.clone(),
                position: cand.position,
                reference: cand.reference_allele,
                alternate: cand.alternate_allele,
            }),
            grch38: hg38.get(&idx).cloned(),
            // Unchanged by the swap: the derived allele is stored as a base precisely so orientation
            // cannot invert its meaning.
            archaic_derived_allele: cand.derived,
            calls: cand.calls,
            diagnostic_class: classify_diagnostic(&cand.calls),
            afr_freq: afr,
        });
    }

    eprintln!(
        "{} sites kept ({no_af} no outgroup AF, {afr_common} too common in AFR, {too_rare} too rare \
         outside AFR, {palindromic} palindromic, {unoriented} dropped by CHM13 orientation, \
         {flipped} ref/alt-swapped)",
        sites.len()
    );
    anyhow::ensure!(
        !sites.is_empty(),
        "no site survived — check the thresholds and the reference"
    );

    let n_nea = sites
        .iter()
        .filter(|s| s.diagnostic_class == navigator_analysis::archaic::DiagnosticClass::Neanderthal)
        .count();
    let n_den = sites
        .iter()
        .filter(|s| s.diagnostic_class == navigator_analysis::archaic::DiagnosticClass::Denisovan)
        .count();
    eprintln!(
        "diagnostic split: {n_nea} Neanderthal, {n_den} Denisovan, {} shared",
        sites.len() - n_nea - n_den
    );

    if let Some(tsv) = &args.sites_tsv {
        let mut w = BufWriter::new(File::create(tsv)?);
        writeln!(w, "#contig\tpos\tref\talt\tderived\tclass\tafr_freq")?;
        for s in &sites {
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{:?}\t{:.5}",
                s.contig,
                s.position,
                s.reference_allele,
                s.alternate_allele,
                s.archaic_derived_allele,
                s.diagnostic_class,
                s.afr_freq
            )?;
        }
        w.flush()?;
    }

    let panel = ArchaicMarkerPanel {
        build: BUILD.to_string(),
        thresholds: ArchaicPanelThresholds {
            max_afr_freq: args.max_afr_freq,
            min_non_afr_freq: args.min_non_afr_freq,
        },
        sites,
    };
    let bytes = panel.to_bytes().map_err(|e| anyhow::anyhow!("{e}"))?;
    write_bin(&args.out, &bytes)?;
    eprintln!(
        "wrote {} ({} sites, {} possible copies)",
        args.out.display(),
        panel.len(),
        panel.possible_copies()
    );
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn polarity_needs_a_confident_ancestral_base_matching_one_allele() {
        // Ancestral == REF → ALT is derived, and vice versa.
        assert_eq!(derived_allele('A', 'A', 'G'), Some('G'));
        assert_eq!(derived_allele('G', 'A', 'G'), Some('A'));
        // Matching neither allele cannot be polarized.
        assert_eq!(derived_allele('C', 'A', 'G'), None);
        // Gaps / unknowns.
        assert_eq!(derived_allele('N', 'A', 'G'), None);
        assert_eq!(derived_allele('.', 'A', 'G'), None);
        assert_eq!(derived_allele('-', 'A', 'G'), None);
        // Lower case = low-confidence EPO call, rejected: an inverted polarity would flip the
        // meaning of every archaic call at the site.
        assert_eq!(derived_allele('a', 'A', 'G'), None);
    }

    #[test]
    fn gt_parses_both_separators_and_missingness() {
        assert_eq!(parse_gt("0/1", 'A', 'G'), (Some('A'), Some('G')));
        assert_eq!(parse_gt("1|1", 'A', 'G'), (Some('G'), Some('G')));
        assert_eq!(parse_gt("./.", 'A', 'G'), (None, None));
        assert_eq!(parse_gt("0/.", 'A', 'G'), (Some('A'), None));
        // Haploid record reads as homozygous.
        assert_eq!(parse_gt("1", 'A', 'G'), (Some('G'), Some('G')));
    }

    #[test]
    fn call_state_is_relative_to_the_derived_allele() {
        // Derived is the ALT here.
        assert_eq!(call_state((Some('G'), Some('G')), 'G'), ArchaicCall::HomDerived);
        assert_eq!(call_state((Some('A'), Some('G')), 'G'), ArchaicCall::Het);
        assert_eq!(call_state((Some('A'), Some('A')), 'G'), ArchaicCall::HomAncestral);
        // Derived is the REF — the case a dosage-based reading would inverted.
        assert_eq!(call_state((Some('A'), Some('A')), 'A'), ArchaicCall::HomDerived);
        assert_eq!(call_state((None, Some('G')), 'G'), ArchaicCall::NoCall);
    }

    #[test]
    fn derived_freq_reexpresses_against_the_derived_base() {
        // Outgroup states AF for its own ALT (=G); derived is G → frequency passes through.
        assert_eq!(derived_freq('G', 'A', 'G', 0.02), Some(0.02));
        // Derived is the outgroup's REF → complement.
        assert_eq!(derived_freq('A', 'A', 'G', 0.02), Some(0.98));
        // Disagreeing allele pairs are not silently filtered on.
        assert_eq!(derived_freq('C', 'A', 'G', 0.02), None);
    }

    #[test]
    fn reference_confident_record_is_hom_derived_when_the_reference_carries_the_derived_allele() {
        // The EVA VCFs are all-sites: a genome matching hg19 emits `ALT=.` / `GT=0/0`. When the EPO
        // ancestral base says the REFERENCE is derived, that genome is homozygous-DERIVED — the
        // donor state. Reading such a record as "no call" (or discarding it as non-variant) drops
        // every site where hg19 itself carries the archaic allele, a systematic, one-sided loss.
        let (reference_allele, alternate_allele) = ('A', 'G');
        // Ancestral is G, so the derived allele is the reference base A.
        let derived = derived_allele('G', reference_allele, alternate_allele).expect("polarizable");
        assert_eq!(derived, 'A');

        // A reference-confident record: alt is absent, so both alleles are the REF base.
        let alleles = parse_gt("0/0", reference_allele, alternate_allele);
        assert_eq!(call_state(alleles, derived), ArchaicCall::HomDerived);

        // And the mirror case: ancestral is the REF, so the ALT is derived and a hom-ref genome is
        // homozygous ANCESTRAL — which must be distinguishable from a masked-out NoCall.
        let derived2 = derived_allele('A', reference_allele, alternate_allele).expect("polarizable");
        assert_eq!(derived2, 'G');
        assert_eq!(call_state(alleles, derived2), ArchaicCall::HomAncestral);
    }

    #[test]
    fn call_tokens_round_trip() {
        for c in [
            ArchaicCall::HomAncestral,
            ArchaicCall::Het,
            ArchaicCall::HomDerived,
            ArchaicCall::NoCall,
        ] {
            assert_eq!(parse_call_token(call_token(c)), c);
        }
    }

    #[test]
    fn contig_normalization_joins_chr_prefixed_and_bare_sources() {
        assert_eq!(normalize_contig("chr1"), "1");
        assert_eq!(normalize_contig("1"), "1");
        assert_eq!(normalize_contig("chrX"), "X");
    }
}
