//! Whole-VCF liftover between reference builds — the GATK `LiftoverVcf` replacement (no external
//! tools). Operates on raw VCF lines so INFO/FORMAT/sample columns pass through verbatim; only the
//! parts liftover actually changes are rewritten: CHROM (contig-name normalized to the target),
//! POS (via the UCSC chain), and on a reverse-strand (inverted) lift the REF/ALT alleles are
//! reverse-complemented. REF/ALT-swap recovery reads the target reference base and, when the lifted
//! REF no longer matches it, swaps REF↔ALT (flipping a biallelic single-sample GT). Records whose
//! position doesn't map, whose multi-base REF straddles a chain break, or that can't be safely
//! recovered are dropped and tallied. Output is coordinate-sorted (a lift can reorder/invert).
//!
//! Reuses the `du_bio` chain primitives and mirrors the drop-with-stats shape of
//! [`crate::gateway::ReferenceGateway::lift_hipstr_bed`].

use std::collections::{HashMap, HashSet};
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use du_bio::liftover::Liftover;
use noodles::fasta;

use crate::error::RefgenomeError;

/// Options for [`lift_vcf`].
#[derive(Debug, Clone, Default)]
pub struct VcfLiftOpts {
    /// Drop variants that land in the target chrY pseudoautosomal regions (PAR) — these are
    /// X/Y-ambiguous on the target and usually unwanted in a Y-only lift.
    pub filter_par: bool,
}

/// Lift / drop counts from [`lift_vcf`].
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub struct VcfLiftStats {
    pub total: usize,
    pub lifted: usize,
    /// Position fell in a chain gap / non-syntenic region.
    pub unmapped: usize,
    /// A multi-base REF whose endpoints lifted to different target contigs (straddled a break).
    pub split: usize,
    /// The lifted REF matched neither the target base nor any ALT (couldn't recover).
    pub ref_mismatch: usize,
    /// A REF/ALT swap was needed but couldn't be applied safely (multiallelic or multi-sample).
    pub swap_ambiguous: usize,
    /// Dropped in the target PAR (only when `filter_par`).
    pub par: usize,
    /// An indel/complex allele on a reverse-strand (inverted) lift — not safely representable.
    pub complex_reverse: usize,
}

/// Best-effort source-build inference from a VCF header (`##reference` / `##contig` / assembly
/// tokens). Returns a canonical build label, or `None` when nothing recognizable is found (the
/// caller then requires an explicit source build). Reads only the header lines.
pub fn infer_source_build(path: &Path) -> Option<String> {
    let reader = open_maybe_gz(path).ok()?;
    for line in reader.lines().map_while(Result::ok) {
        if !line.starts_with('#') {
            break; // header ended
        }
        let l = line.to_ascii_lowercase();
        if l.contains("chm13") || l.contains("hs1") || l.contains("t2t") {
            return Some("chm13v2.0".into());
        }
        if l.contains("grch38") || l.contains("hg38") {
            return Some("GRCh38".into());
        }
        if l.contains("grch37") || l.contains("hg19") || l.contains("b37") {
            return Some("GRCh37".into());
        }
    }
    None
}

/// Reverse-complement a DNA allele (ACGTN, case-insensitive → uppercase). Non-ACGTN passes through
/// unchanged (e.g. a `*` or symbolic allele), so callers should gate revcomp to simple alleles.
pub fn revcomp(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c.to_ascii_uppercase() {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'N' => 'N',
            other => other,
        })
        .collect()
}

/// True for a simple single-base ACGT allele (the case we can revcomp + swap-recover precisely).
fn is_snv_allele(a: &str) -> bool {
    a.len() == 1 && a.bytes().all(|b| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T'))
}

/// Candidate names for a lifted (chain-query) contig in the **target** FASTA's naming style, in
/// preference order — covers `chr` prefix presence and `chrM`/`MT`.
fn target_contig_name(q_name: &str, target_names: &HashSet<String>) -> Option<String> {
    let bare = q_name.strip_prefix("chr").unwrap_or(q_name);
    let mut cands = vec![q_name.to_string(), bare.to_string(), format!("chr{bare}")];
    if bare.eq_ignore_ascii_case("M") || bare.eq_ignore_ascii_case("MT") {
        cands.extend(["chrM".into(), "MT".into(), "chrMT".into(), "M".into()]);
    }
    cands.into_iter().find(|c| target_names.contains(c))
}

/// Map a VCF CHROM to the chain's source (`t_name`) naming. T2T/UCSC chains use `chr`-prefixed
/// names, so normalize bare NCBI names (`1`, `MT`) to `chr1` / `chrM`.
fn to_chain_source_name(chrom: &str) -> String {
    if chrom.starts_with("chr") {
        return chrom.to_string();
    }
    match chrom {
        "MT" | "M" => "chrM".to_string(),
        other => format!("chr{other}"),
    }
}

/// Read the 1-based reference base at `contig:pos` from an indexed FASTA, uppercased. `None` if the
/// contig/position can't be queried.
fn ref_base_at(reader: &mut fasta::io::IndexedReader<fasta::io::BufReader<std::fs::File>>, contig: &str, pos: i64) -> Option<u8> {
    let region: noodles::core::Region = format!("{contig}:{pos}-{pos}").parse().ok()?;
    let rec = reader.query(&region).ok()?;
    rec.sequence().as_ref().first().map(|b| b.to_ascii_uppercase())
}

/// A lifted record buffered for the coordinate sort before writing.
struct Lifted {
    contig_rank: usize,
    pos: i64,
    line: String,
}

/// Lift every record in `in_vcf` from the source build to the target build, writing `out_vcf`
/// (gzip when the path ends `.gz`). `lo` is the source→target chain (load it first via the
/// gateway); `target_fa` is the indexed target FASTA (for REF/ALT-swap recovery + the contig set);
/// `target_par` is the target chrY PAR intervals (0-based half-open) used only when
/// `opts.filter_par`. `source_label`/`target_label` are stamped into the header provenance.
#[allow(clippy::too_many_arguments)]
pub fn lift_vcf(
    lo: &Liftover,
    target_fa: &Path,
    target_par: &[(i64, i64)],
    source_label: &str,
    target_label: &str,
    in_vcf: &Path,
    out_vcf: &Path,
    opts: VcfLiftOpts,
) -> Result<VcfLiftStats, RefgenomeError> {
    // Target contig set + order (from the .fai) for naming normalization, ##contig rewrite, sort.
    let fai = read_fai(target_fa)?;
    let target_names: HashSet<String> = fai.iter().map(|(n, _)| n.clone()).collect();
    let contig_rank: HashMap<String, usize> =
        fai.iter().enumerate().map(|(i, (n, _))| (n.clone(), i)).collect();

    let mut fasta_reader = fasta::io::indexed_reader::Builder::default()
        .build_from_path(target_fa)
        .map_err(|e| RefgenomeError::io(target_fa, e))?;

    let in_reader: Box<dyn BufRead> = open_maybe_gz(in_vcf)?;
    let mut header: Vec<String> = Vec::new();
    let mut chrom_line: Option<String> = None;
    let mut lifted: Vec<Lifted> = Vec::new();
    let mut stats = VcfLiftStats::default();

    for line in in_reader.lines() {
        let line = line.map_err(|e| RefgenomeError::io(in_vcf, e))?;
        if line.starts_with("##") {
            if line.starts_with("##contig=") {
                continue; // regenerated from the target .fai below
            }
            header.push(line);
            continue;
        }
        if line.starts_with('#') {
            chrom_line = Some(line);
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        stats.total += 1;
        if let Some(rec) = lift_record(&line, lo, &target_names, &contig_rank, &mut fasta_reader, target_par, &opts, &mut stats) {
            lifted.push(rec);
        }
    }

    // Coordinate-sort (lift can reorder / invert).
    lifted.sort_by(|a, b| a.contig_rank.cmp(&b.contig_rank).then(a.pos.cmp(&b.pos)));

    // Write: passthrough ## header, regenerated ##contig lines + provenance, #CHROM, sorted records.
    let mut out = open_out(out_vcf)?;
    for h in &header {
        writeln!(out, "{h}").map_err(|e| RefgenomeError::io(out_vcf, e))?;
    }
    writeln!(out, "##liftover=navigator; source={source_label}; target={target_label}")
        .map_err(|e| RefgenomeError::io(out_vcf, e))?;
    for (name, len) in &fai {
        writeln!(out, "##contig=<ID={name},length={len}>").map_err(|e| RefgenomeError::io(out_vcf, e))?;
    }
    if let Some(c) = &chrom_line {
        writeln!(out, "{c}").map_err(|e| RefgenomeError::io(out_vcf, e))?;
    }
    for rec in &lifted {
        writeln!(out, "{}", rec.line).map_err(|e| RefgenomeError::io(out_vcf, e))?;
    }
    out.flush().map_err(|e| RefgenomeError::io(out_vcf, e))?;
    Ok(stats)
}

/// Lift one data line; returns the rewritten record or `None` (a drop, with `stats` updated).
#[allow(clippy::too_many_arguments)]
fn lift_record(
    line: &str,
    lo: &Liftover,
    target_names: &HashSet<String>,
    contig_rank: &HashMap<String, usize>,
    fasta_reader: &mut fasta::io::IndexedReader<fasta::io::BufReader<std::fs::File>>,
    target_par: &[(i64, i64)],
    opts: &VcfLiftOpts,
    stats: &mut VcfLiftStats,
) -> Option<Lifted> {
    let mut f: Vec<String> = line.split('\t').map(str::to_string).collect();
    if f.len() < 8 {
        stats.unmapped += 1; // malformed → treat as undroppable-but-skip
        return None;
    }
    let Ok(src_pos) = f[1].parse::<i64>() else {
        stats.unmapped += 1; // unparseable POS
        return None;
    };
    let src_name = to_chain_source_name(&f[0]);

    // Lift the start position; capture the target strand (inverted tracts on the CHM13 Y).
    let lifted = lo
        .chains
        .iter()
        .filter(|c| c.t_name == src_name)
        .find_map(|c| c.lift(src_pos - 1).map(|q| (c.q_name.clone(), q + 1, c.q_strand == '-')));
    let Some((q_name, q_pos, reverse)) = lifted else {
        stats.unmapped += 1; // position fell in a chain gap / non-syntenic region
        return None;
    };
    let Some(target_contig) = target_contig_name(&q_name, target_names) else {
        stats.unmapped += 1;
        return None;
    };

    let ref_allele = f[3].clone();
    let alts: Vec<String> = f[4].split(',').map(str::to_string).collect();
    let simple = is_snv_allele(&ref_allele) && alts.iter().all(|a| is_snv_allele(a));

    // Indel / complex allele on an inverted lift: not safely representable here.
    if reverse && !simple {
        stats.complex_reverse += 1;
        return None;
    }
    // A multi-base REF must lift contiguously (both endpoints on the same target contig).
    if ref_allele.len() > 1 {
        let end = src_pos + ref_allele.len() as i64 - 1;
        match lo.chains.iter().filter(|c| c.t_name == src_name).find_map(|c| c.lift(end - 1)) {
            Some(_) => {}
            None => {
                stats.split += 1;
                return None;
            }
        }
    }

    // Reverse-complement simple alleles on an inverted lift.
    let (mut new_ref, mut new_alts) = if reverse && simple {
        (revcomp(&ref_allele), alts.iter().map(|a| revcomp(a)).collect::<Vec<_>>())
    } else {
        (ref_allele.clone(), alts.clone())
    };

    // REF/ALT-swap recovery against the target reference base (SNVs only).
    let mut swapped = false;
    if simple {
        if let Some(tbase) = ref_base_at(fasta_reader, &target_contig, q_pos) {
            let tb = (tbase as char).to_string();
            if !new_ref.eq_ignore_ascii_case(&tb) {
                // REF doesn't match the target base — try to recover by swapping with a matching ALT.
                if let Some(idx) = new_alts.iter().position(|a| a.eq_ignore_ascii_case(&tb)) {
                    if new_alts.len() != 1 {
                        stats.swap_ambiguous += 1; // multiallelic swap — ambiguous to relabel
                        return None;
                    }
                    std::mem::swap(&mut new_ref, &mut new_alts[idx]);
                    swapped = true;
                } else {
                    stats.ref_mismatch += 1;
                    return None;
                }
            }
        }
    }

    // PAR filter (target chrY).
    let on_y = target_contig.eq_ignore_ascii_case("chrY") || target_contig.eq_ignore_ascii_case("Y");
    let z = q_pos - 1; // 0-based, to match the half-open PAR intervals
    if opts.filter_par && on_y && target_par.iter().any(|&(s, e)| z >= s && z < e) {
        stats.par += 1;
        return None;
    }

    // Rewrite the changed columns.
    f[0] = target_contig.clone();
    f[1] = q_pos.to_string();
    f[3] = new_ref;
    f[4] = new_alts.join(",");
    // Flip a biallelic single-sample GT on a REF/ALT swap (0↔1).
    if swapped {
        flip_biallelic_gt(&mut f);
    }

    stats.lifted += 1;
    Some(Lifted {
        contig_rank: contig_rank.get(&target_contig).copied().unwrap_or(usize::MAX),
        pos: q_pos,
        line: f.join("\t"),
    })
}

/// Flip the allele indices of a biallelic single-sample genotype (0↔1) in the first sample column,
/// after a REF/ALT swap. No-op when there's no FORMAT/sample (sites-only VCF).
fn flip_biallelic_gt(f: &mut [String]) {
    if f.len() < 10 {
        return; // no FORMAT + sample columns
    }
    let gt_idx = f[8].split(':').position(|k| k == "GT");
    let Some(gt_idx) = gt_idx else { return };
    let mut parts: Vec<String> = f[9].split(':').map(str::to_string).collect();
    if let Some(gt) = parts.get_mut(gt_idx) {
        *gt = gt
            .chars()
            .map(|c| match c {
                '0' => '1',
                '1' => '0',
                other => other, // '/', '|', '.'
            })
            .collect();
    }
    f[9] = parts.join(":");
}

/// Parse a `.fai` next to `fa` into `(contig, length)` in file order.
fn read_fai(fa: &Path) -> Result<Vec<(String, i64)>, RefgenomeError> {
    let mut fai = fa.as_os_str().to_os_string();
    fai.push(".fai");
    let fai = std::path::PathBuf::from(fai);
    let text = std::fs::read_to_string(&fai).map_err(|e| RefgenomeError::io(&fai, e))?;
    let mut out = Vec::new();
    for line in text.lines() {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() >= 2 {
            if let Ok(len) = cols[1].parse::<i64>() {
                out.push((cols[0].to_string(), len));
            }
        }
    }
    Ok(out)
}

/// Open a `.vcf` or `.vcf.gz` for buffered line reading.
fn open_maybe_gz(path: &Path) -> Result<Box<dyn BufRead>, RefgenomeError> {
    let file = std::fs::File::open(path).map_err(|e| RefgenomeError::io(path, e))?;
    if path.extension().is_some_and(|e| e == "gz") {
        Ok(Box::new(BufReader::new(flate2::read::MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

/// Open the output VCF (gzip when the path ends `.gz`).
fn open_out(path: &Path) -> Result<Box<dyn Write>, RefgenomeError> {
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent).map_err(|e| RefgenomeError::io(parent, e))?;
    }
    let file = std::fs::File::create(path).map_err(|e| RefgenomeError::io(path, e))?;
    if path.extension().is_some_and(|e| e == "gz") {
        Ok(Box::new(flate2::write::GzEncoder::new(file, flate2::Compression::default())))
    } else {
        Ok(Box::new(std::io::BufWriter::new(file)))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn revcomp_basic() {
        assert_eq!(revcomp("ACGT"), "ACGT");
        assert_eq!(revcomp("A"), "T");
        assert_eq!(revcomp("AC"), "GT");
        assert_eq!(revcomp("acgt"), "ACGT");
    }

    #[test]
    fn target_name_matches_style() {
        let names: HashSet<String> = ["chr1", "chrY", "chrM"].iter().map(|s| s.to_string()).collect();
        assert_eq!(target_contig_name("chrY", &names).as_deref(), Some("chrY"));
        assert_eq!(target_contig_name("Y", &names).as_deref(), Some("chrY"));
        assert_eq!(target_contig_name("MT", &names).as_deref(), Some("chrM"));
        let bare: HashSet<String> = ["1", "Y", "MT"].iter().map(|s| s.to_string()).collect();
        assert_eq!(target_contig_name("chrY", &bare).as_deref(), Some("Y"));
        assert_eq!(target_contig_name("chrM", &bare).as_deref(), Some("MT"));
    }

    #[test]
    fn flip_gt_biallelic() {
        let mut f: Vec<String> =
            "chrY\t100\t.\tA\tG\t.\t.\t.\tGT:DP\t0/0:30".split('\t').map(str::to_string).collect();
        flip_biallelic_gt(&mut f);
        assert_eq!(f[9], "1/1:30");
        let mut het: Vec<String> =
            "chrY\t100\t.\tA\tG\t.\t.\t.\tGT\t0|1".split('\t').map(str::to_string).collect();
        flip_biallelic_gt(&mut het);
        assert_eq!(het[9], "1|0");
    }

    /// End-to-end: a reverse-strand (inverted) chain lift — REF/ALT reverse-complemented, REF/ALT
    /// swapped to match the target reference base, and the single-sample GT flipped accordingly.
    #[test]
    fn lift_vcf_reverse_strand_revcomp_swap_and_gt_flip() {
        let dir = std::env::temp_dir().join(format!("dun-vcflift-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();

        // Target chrY (20 bp): pos1='A'. Build the .fai in place.
        let fa = dir.join("target.fa");
        std::fs::write(&fa, ">chrY\nAAAAACCCCCGGGGGTTTTT\n").unwrap();
        crate::index::decompress_and_index(&fa, &fa).unwrap();

        // Inverted chain: source chrY 0-based t → target 0-based 19-t (q_strand '-').
        let lo = Liftover::parse("chain 1 chrY 20 + 0 20 chrY 20 - 0 20 1\n20\n").unwrap();

        // Source record at chrY:20 (0-based 19 → target 0-based 0 → 1-based 1, reverse).
        // Source REF=C ALT=T → revcomp REF=G ALT=A; target base@1 = A, so swap REF/ALT and flip GT.
        let in_vcf = dir.join("in.vcf");
        std::fs::write(
            &in_vcf,
            "##fileformat=VCFv4.2\n##contig=<ID=chrY,length=20>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\nchrY\t20\trs1\tC\tT\t.\t.\t.\tGT\t0/0\n",
        )
        .unwrap();

        let out_vcf = dir.join("out.vcf");
        let stats = lift_vcf(&lo, &fa, &[], "GRCh38", "chm13v2.0", &in_vcf, &out_vcf, VcfLiftOpts::default()).unwrap();
        assert_eq!(stats.total, 1);
        assert_eq!(stats.lifted, 1);

        let out = std::fs::read_to_string(&out_vcf).unwrap();
        let rec = out.lines().find(|l| !l.starts_with('#')).expect("a record");
        let f: Vec<&str> = rec.split('\t').collect();
        assert_eq!(f[0], "chrY");
        assert_eq!(f[1], "1"); // inverted: source 20 → target 1
        assert_eq!(f[3], "A"); // revcomp(C)=G, then swapped with revcomp(T)=A to match target base A
        assert_eq!(f[4], "G");
        assert_eq!(f[9], "1/1"); // 0/0 flipped on the REF/ALT swap
        // Regenerated provenance + contig header present.
        assert!(out.contains("##liftover=navigator"));
        assert!(out.contains("##contig=<ID=chrY,length=20>"));

        let _ = std::fs::remove_dir_all(&dir);
    }

    /// A position in a chain gap is dropped (unmapped), tallied, and not emitted.
    #[test]
    fn lift_vcf_drops_unmapped() {
        let dir = std::env::temp_dir().join(format!("dun-vcflift-unmapped-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let fa = dir.join("t.fa");
        std::fs::write(&fa, ">chrY\nACGTACGTAC\n").unwrap();
        crate::index::decompress_and_index(&fa, &fa).unwrap();
        // Chain covers only t[0,5); a record at t=8 (1-based 9) falls in the gap.
        let lo = Liftover::parse("chain 1 chrY 10 + 0 5 chrY 10 + 0 5 1\n5\n").unwrap();
        let in_vcf = dir.join("in.vcf");
        std::fs::write(&in_vcf, "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchrY\t9\t.\tA\tG\t.\t.\t.\n").unwrap();
        let out_vcf = dir.join("out.vcf");
        let stats = lift_vcf(&lo, &fa, &[], "GRCh38", "chm13v2.0", &in_vcf, &out_vcf, VcfLiftOpts::default()).unwrap();
        assert_eq!((stats.total, stats.lifted, stats.unmapped), (1, 0, 1));
        let _ = std::fs::remove_dir_all(&dir);
    }
}
