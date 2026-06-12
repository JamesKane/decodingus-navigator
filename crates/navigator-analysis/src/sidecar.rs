//! Parse the `ytree` pipeline's per-sample text sidecars into Navigator's existing result
//! structs, so the fast-path import fills the same caches the CRAM walkers would — without
//! touching the alignment.
//!
//! - `.sex` (`male`/`female`) → [`SexInferenceResult`].
//! - `stats.txt` (`samtools stats`) → [`ReadMetrics`] — **fully** populated: the `SN`
//!   summary gives the scalar counts and the `RL`/`IS` histogram lines give the read-length
//!   and insert-size distributions (median/std/min/max), including the `median_insert_size`
//!   the project report shows.
//! - `coverage.txt` (`samtools coverage`) + `callable.summary.txt` (GATK `CallableLoci`) →
//!   a **lite** [`CoverageResult`]: genome-wide mean depth (length-weighted) + per-contig
//!   stats + callable-base counts. The depth histogram and `pct_Nx` / median need the
//!   per-base walk, so they're left zeroed and the result is flagged `partial` by the caller.
//!
//! Unknown numeric fields are `0.0` (not `NaN`) because the cache round-trips through
//! `serde_json`, which encodes `NaN` as `null` and then fails to read it back.

use std::collections::BTreeMap;

use crate::coverage::{ContigCallableMetrics, ContigCoverageStats, CoverageResult};
use crate::read_metrics::{PairOrientation, ReadMetrics};
use crate::sex::{Confidence, InferredSex, SexInferenceResult};

// ---- .sex --------------------------------------------------------------------

/// Parse a `.sex` sidecar (`male` / `female`, case-insensitive). Anything else → Unknown.
/// Ratio fields are unknown from this file (0.0); the label is what the pipeline computed.
pub fn parse_sex(text: &str) -> SexInferenceResult {
    let (inferred, confidence) = match text.trim().to_ascii_lowercase().as_str() {
        "male" | "m" => (InferredSex::Male, Confidence::High),
        "female" | "f" => (InferredSex::Female, Confidence::High),
        _ => (InferredSex::Unknown, Confidence::Low),
    };
    SexInferenceResult {
        inferred_sex: inferred,
        x_autosome_ratio: 0.0,
        autosome_mean_coverage: 0.0,
        x_coverage: 0.0,
        confidence,
    }
}

// ---- stats.txt (samtools stats) ----------------------------------------------

/// Parse `samtools stats` output into a fully-populated [`ReadMetrics`]. `SN` lines give the
/// scalar counts; `RL`/`IS` lines give the read-length / insert-size histograms (and thus
/// their median/std/min/max). `mean_mapping_quality` isn't emitted by samtools stats → 0.0.
pub fn parse_samtools_stats(text: &str) -> ReadMetrics {
    let mut sn: BTreeMap<&str, f64> = BTreeMap::new();
    let mut rl: BTreeMap<u32, u64> = BTreeMap::new();
    let mut is: BTreeMap<u32, u64> = BTreeMap::new();

    for line in text.lines() {
        let mut col = line.split('\t');
        match col.next() {
            Some("SN") => {
                // SN \t <key>: \t <value> [\t # comment]
                if let (Some(key), Some(val)) = (col.next(), col.next()) {
                    let key = key.trim_end_matches(':').trim();
                    if let Ok(v) = val.trim().parse::<f64>() {
                        sn.insert(key, v);
                    }
                }
            }
            // RL \t <read length> \t <count>
            Some("RL") => add_hist(&mut rl, &mut col),
            // IS \t <insert size> \t <pairs total> \t inward \t outward \t other
            Some("IS") => add_hist(&mut is, &mut col),
            _ => {}
        }
    }

    let get = |k: &str| sn.get(k).copied().unwrap_or(0.0);
    let total_reads = get("raw total sequences") as u64;
    let qc_failed = get("reads QC failed") as u64;
    let pf_reads = total_reads.saturating_sub(qc_failed);
    let pf_reads_aligned = get("reads mapped") as u64;
    let reads_aligned_in_pairs = get("reads mapped and paired") as u64;
    let proper_pairs = get("reads properly paired") as u64;
    let pairs_diff_chrom = get("pairs on different chromosomes") as u64;

    let pct = |num: u64, den: u64| if den == 0 { 0.0 } else { 100.0 * num as f64 / den as f64 };

    let (rl_mean, rl_median, rl_std, rl_min, rl_max) = summarize(&rl);
    let (is_mean, is_median, is_std, is_min, is_max) = summarize(&is);
    // Prefer samtools' own averages where it reports them (it rounds the histogram anyway).
    let mean_insert_size = if sn.contains_key("insert size average") { get("insert size average") } else { is_mean };
    let std_insert_size = if sn.contains_key("insert size standard deviation") { get("insert size standard deviation") } else { is_std };

    ReadMetrics {
        total_reads,
        pf_reads,
        pf_reads_aligned,
        reads_aligned_in_pairs,
        proper_pairs,
        pct_pf_reads_aligned: pct(pf_reads_aligned, total_reads),
        pct_reads_aligned_in_pairs: pct(reads_aligned_in_pairs, pf_reads_aligned),
        pct_proper_pairs: pct(proper_pairs, total_reads),
        median_read_length: rl_median,
        mean_read_length: if rl.is_empty() { get("average length") } else { rl_mean },
        std_read_length: rl_std,
        min_read_length: rl_min,
        max_read_length: rl_max,
        read_length_histogram: rl,
        median_insert_size: is_median,
        mean_insert_size,
        std_insert_size,
        min_insert_size: is_min,
        max_insert_size: is_max,
        insert_size_histogram: is,
        // samtools stats doesn't classify orientation; Illumina paired-end is FR.
        pair_orientation: PairOrientation::Fr,
        // Picard-style chimera rate: read pairs mapping to different chromosomes / total pairs.
        pct_chimeras: pct(pairs_diff_chrom, total_reads / 2),
        mean_mapping_quality: 0.0,
    }
}

/// Read `<bin> \t <count>` from the remaining columns of an `RL`/`IS` line into `hist`.
fn add_hist<'a>(hist: &mut BTreeMap<u32, u64>, col: &mut impl Iterator<Item = &'a str>) {
    if let (Some(bin), Some(count)) = (col.next(), col.next()) {
        if let (Ok(b), Ok(c)) = (bin.trim().parse::<u32>(), count.trim().parse::<u64>()) {
            if c > 0 {
                *hist.entry(b).or_default() += c;
            }
        }
    }
}

/// `(mean, median, std, min, max)` of a value→count histogram. Zeros for an empty histogram.
fn summarize(hist: &BTreeMap<u32, u64>) -> (f64, f64, f64, u32, u32) {
    let total: u64 = hist.values().sum();
    if total == 0 {
        return (0.0, 0.0, 0.0, 0, 0);
    }
    let sum: f64 = hist.iter().map(|(&v, &c)| v as f64 * c as f64).sum();
    let mean = sum / total as f64;
    let var: f64 = hist.iter().map(|(&v, &c)| c as f64 * (v as f64 - mean).powi(2)).sum::<f64>() / total as f64;
    let min = *hist.keys().next().unwrap();
    let max = *hist.keys().next_back().unwrap();
    // Median: first bin whose cumulative count crosses the halfway point.
    let half = total.div_ceil(2);
    let mut cum = 0u64;
    let mut median = min as f64;
    for (&v, &c) in hist {
        cum += c;
        if cum >= half {
            median = v as f64;
            break;
        }
    }
    (mean, median, var.sqrt(), min, max)
}

// ---- coverage.txt + callable.summary.txt -------------------------------------

/// Parse `samtools coverage` TSV → per-contig stats + the length-weighted genome-wide mean
/// depth and total territory (sum of contig lengths). Skips the `#rname …` header.
pub fn parse_samtools_coverage(text: &str) -> (f64, u64, Vec<ContigCoverageStats>) {
    let mut contigs = Vec::new();
    let mut weighted_depth = 0.0f64;
    let mut territory = 0u64;
    for line in text.lines() {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 9 {
            continue;
        }
        let parse_u = |i: usize| f[i].trim().parse::<u64>().unwrap_or(0);
        let parse_f = |i: usize| f[i].trim().parse::<f64>().unwrap_or(0.0);
        let end_pos = parse_u(2);
        let mean_depth = parse_f(6);
        weighted_depth += mean_depth * end_pos as f64;
        territory += end_pos;
        contigs.push(ContigCoverageStats {
            contig: f[0].to_string(),
            start_pos: parse_u(1),
            end_pos,
            num_reads: parse_u(3),
            cov_bases: parse_u(4),
            coverage: parse_f(5),
            mean_depth,
            mean_base_q: parse_f(7),
            mean_map_q: parse_f(8),
            histogram: Vec::new(), // fast-path: samtools coverage has no per-depth histogram
        });
    }
    let mean_coverage = if territory == 0 { 0.0 } else { weighted_depth / territory as f64 };
    (mean_coverage, territory, contigs)
}

/// Parse a GATK `CallableLoci` summary (`state nBases` blocks, one per contig). Contig blocks
/// are introduced by a `--- <contig> ---` line; the **leading headerless block** in a
/// `.chrYM.` summary is the mitochondrion (`chrM`). Returns total CALLABLE bases and the
/// per-contig breakdown.
pub fn parse_callable_summary(text: &str) -> (u64, Vec<ContigCallableMetrics>) {
    let mut out: Vec<ContigCallableMetrics> = Vec::new();
    // The first block precedes any `--- … ---` header → mitochondrion in the chrYM summary.
    let mut current = ContigCallableMetrics {
        contig: "chrM".to_string(),
        ref_n: 0,
        callable: 0,
        no_coverage: 0,
        low_coverage: 0,
        excessive_coverage: 0,
        poor_mapping_quality: 0,
    };
    let mut started = false;
    let mut total_callable = 0u64;

    let flush = |out: &mut Vec<ContigCallableMetrics>, c: &ContigCallableMetrics, started: bool| {
        if started {
            out.push(c.clone());
        }
    };

    for line in text.lines() {
        let t = line.trim();
        if let Some(rest) = t.strip_prefix("---") {
            // New contig section: flush the previous block, start a fresh one.
            flush(&mut out, &current, started);
            let name = rest.trim_end_matches('-').trim().to_string();
            current = ContigCallableMetrics {
                contig: name,
                ref_n: 0,
                callable: 0,
                no_coverage: 0,
                low_coverage: 0,
                excessive_coverage: 0,
                poor_mapping_quality: 0,
            };
            started = true;
            continue;
        }
        // `<STATE> <nBases>` rows (whitespace-padded); the `state nBases` header has a
        // non-numeric second token and is skipped by the parse.
        let mut it = t.split_whitespace();
        let (Some(state), Some(val)) = (it.next(), it.next()) else { continue };
        let Ok(n) = val.parse::<u64>() else { continue };
        started = true;
        match state {
            "REF_N" => current.ref_n = n,
            "CALLABLE" => {
                current.callable = n;
                total_callable += n;
            }
            "NO_COVERAGE" => current.no_coverage = n,
            "LOW_COVERAGE" => current.low_coverage = n,
            "EXCESSIVE_COVERAGE" => current.excessive_coverage = n,
            "POOR_MAPPING_QUALITY" => current.poor_mapping_quality = n,
            _ => {}
        }
    }
    flush(&mut out, &current, started);
    (total_callable, out)
}

/// Assemble a **lite** [`CoverageResult`] from the coverage + callable sidecars. Mean depth
/// and callable counts are real; median/sd/histogram/`pct_Nx` need the per-base walk and are
/// left zeroed — the caller records this artifact as `partial` so the deep pass upgrades it.
pub fn lite_coverage(coverage_txt: &str, callable_summary: Option<&str>) -> CoverageResult {
    let (mean_coverage, genome_territory, contig_coverage_stats) = parse_samtools_coverage(coverage_txt);
    let (callable_bases, contig_callable) =
        callable_summary.map(parse_callable_summary).unwrap_or((0, Vec::new()));
    CoverageResult {
        genome_territory,
        mean_coverage,
        median_coverage: 0.0,
        sd_coverage: 0.0,
        coverage_histogram: Vec::new(),
        pct_1x: 0.0,
        pct_5x: 0.0,
        pct_10x: 0.0,
        pct_15x: 0.0,
        pct_20x: 0.0,
        pct_25x: 0.0,
        pct_30x: 0.0,
        pct_40x: 0.0,
        pct_50x: 0.0,
        callable_bases,
        contig_callable,
        contig_coverage_stats,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sex_parses_label_and_confidence() {
        assert_eq!(parse_sex("male\n").inferred_sex, InferredSex::Male);
        assert_eq!(parse_sex(" Female ").inferred_sex, InferredSex::Female);
        assert_eq!(parse_sex("M").inferred_sex, InferredSex::Male);
        let u = parse_sex("unknown");
        assert_eq!(u.inferred_sex, InferredSex::Unknown);
        assert_eq!(u.confidence, Confidence::Low);
    }

    const STATS: &str = "\
# comment line
SN\traw total sequences:\t1000\t# excluding supplementary
SN\treads QC failed:\t0
SN\treads mapped:\t980
SN\treads mapped and paired:\t970
SN\treads properly paired:\t960
SN\tinsert size average:\t430.6
SN\tinsert size standard deviation:\t97.7
SN\tpairs on different chromosomes:\t10
SN\taverage length:\t150
RL\t150\t1000
IS\t400\t300\t0\t10\t290
IS\t450\t500\t0\t12\t488
IS\t500\t200\t0\t5\t195
";

    #[test]
    fn samtools_stats_scalars_and_histograms() {
        let m = parse_samtools_stats(STATS);
        assert_eq!(m.total_reads, 1000);
        assert_eq!(m.pf_reads_aligned, 980);
        assert_eq!(m.reads_aligned_in_pairs, 970);
        assert_eq!(m.proper_pairs, 960);
        assert!((m.pct_pf_reads_aligned - 98.0).abs() < 1e-9);
        assert_eq!(m.mean_read_length, 150.0);
        // Insert histogram: total 1000 pairs, median bin where cumulative ≥ 500 → 450.
        assert_eq!(m.median_insert_size, 450.0);
        assert_eq!(m.min_insert_size, 400);
        assert_eq!(m.max_insert_size, 500);
        assert!((m.mean_insert_size - 430.6).abs() < 1e-9, "uses samtools' reported average");
        assert_eq!(m.insert_size_histogram.get(&450), Some(&500));
    }

    const COVERAGE: &str = "\
#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq
chr1\t1\t1000\t100\t990\t99.0\t30.0\t29.6\t55.0
chrY\t1\t500\t20\t400\t80.0\t10.0\t29.0\t40.0
";

    #[test]
    fn samtools_coverage_length_weighted_mean() {
        let (mean, territory, contigs) = parse_samtools_coverage(COVERAGE);
        assert_eq!(territory, 1500);
        // (30*1000 + 10*500) / 1500 = 35000/1500.
        assert!((mean - (35000.0 / 1500.0)).abs() < 1e-9);
        assert_eq!(contigs.len(), 2);
        assert_eq!(contigs[1].contig, "chrY");
        assert!((contigs[1].mean_depth - 10.0).abs() < 1e-9);
    }

    const CALLABLE: &str = "\
                         state nBases
                         REF_N 0
                      CALLABLE 16249
                   NO_COVERAGE 0
                  LOW_COVERAGE 0
            EXCESSIVE_COVERAGE 0
          POOR_MAPPING_QUALITY 320
--- chrY ---
                         state nBases
                         REF_N 0
                      CALLABLE 16627537
                   NO_COVERAGE 1127957
                  LOW_COVERAGE 13296366
            EXCESSIVE_COVERAGE 0
          POOR_MAPPING_QUALITY 28612629
";

    #[test]
    fn callable_summary_leading_block_is_mito_then_chr_y() {
        let (total, contigs) = parse_callable_summary(CALLABLE);
        assert_eq!(total, 16249 + 16627537);
        assert_eq!(contigs.len(), 2);
        assert_eq!(contigs[0].contig, "chrM");
        assert_eq!(contigs[0].callable, 16249);
        assert_eq!(contigs[0].poor_mapping_quality, 320);
        assert_eq!(contigs[1].contig, "chrY");
        assert_eq!(contigs[1].callable, 16627537);
        assert_eq!(contigs[1].low_coverage, 13296366);
    }

    #[test]
    fn lite_coverage_combines_sources_and_zeros_deep_fields() {
        let c = lite_coverage(COVERAGE, Some(CALLABLE));
        assert_eq!(c.genome_territory, 1500);
        assert_eq!(c.callable_bases, 16249 + 16627537);
        assert_eq!(c.contig_coverage_stats.len(), 2);
        assert_eq!(c.contig_callable.len(), 2);
        // Deep-walk-only fields stay zeroed (artifact is flagged partial by the caller).
        assert_eq!(c.median_coverage, 0.0);
        assert!(c.coverage_histogram.is_empty());
        assert_eq!(c.pct_10x, 0.0);
    }

    /// Real-data smoke test: parse HG00096's actual pipeline sidecars off the NAS. No-ops
    /// when the share isn't mounted. Run: `cargo test -p navigator-analysis sidecar -- --ignored --nocapture`.
    #[test]
    #[ignore = "reads NAS files; run explicitly"]
    fn real_sidecars_parse() {
        use std::path::Path;
        let dir = Path::new("/Volumes/nas/Genomics/PRJEB31736/HG00096");
        if !dir.exists() {
            eprintln!("skip: {} not mounted", dir.display());
            return;
        }
        let read = |name: &str| std::fs::read_to_string(dir.join(name)).unwrap();
        let sex = parse_sex(&read("HG00096.chm13.sex"));
        assert_eq!(sex.inferred_sex, InferredSex::Male);

        let m = parse_samtools_stats(&read("stats.txt"));
        eprintln!(
            "reads={} %aligned={:.1} mean_rl={:.0} mean_is={:.1} median_is={:.0}",
            m.total_reads, m.pct_pf_reads_aligned, m.mean_read_length, m.mean_insert_size, m.median_insert_size
        );
        assert!(m.total_reads > 0 && m.pf_reads_aligned > 0);
        assert!(m.mean_read_length > 0.0 && m.median_insert_size > 0.0);

        let c = lite_coverage(&read("coverage.txt"), Some(&read("HG00096.chm13.chrYM.callable.summary.txt")));
        eprintln!("mean_cov={:.2} territory={} callable={}", c.mean_coverage, c.genome_territory, c.callable_bases);
        assert!(c.mean_coverage > 0.0 && c.callable_bases > 0);
    }
}
