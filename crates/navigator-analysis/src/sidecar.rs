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

use std::collections::{BTreeMap, HashMap};

use crate::coverage::{ContigCallableMetrics, ContigCoverageStats, CoverageResult};
use crate::read_metrics::{PairOrientation, ReadMetrics};
use crate::sex::{Confidence, InferredSex, SexInferenceResult};

/// Percentage `100·num/den` (the `ReadMetrics` pct convention; 0 when `den == 0`).
fn pct100(num: u64, den: u64) -> f64 {
    if den == 0 {
        0.0
    } else {
        100.0 * num as f64 / den as f64
    }
}

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
    let mean_insert_size = if sn.contains_key("insert size average") {
        get("insert size average")
    } else {
        is_mean
    };
    let std_insert_size = if sn.contains_key("insert size standard deviation") {
        get("insert size standard deviation")
    } else {
        is_std
    };

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
    let var: f64 = hist
        .iter()
        .map(|(&v, &c)| c as f64 * (v as f64 - mean).powi(2))
        .sum::<f64>()
        / total as f64;
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
    let mean_coverage = if territory == 0 {
        0.0
    } else {
        weighted_depth / territory as f64
    };
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
        let (Some(state), Some(val)) = (it.next(), it.next()) else {
            continue;
        };
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
    let (callable_bases, contig_callable) = callable_summary.map(parse_callable_summary).unwrap_or((0, Vec::new()));
    CoverageResult {
        genome_territory,
        mean_coverage,
        median_coverage: 0.0,
        sd_coverage: 0.0,
        mad_coverage: 0.0,
        pct_exc_mapq: 0.0,
        pct_exc_baseq: 0.0,
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

// ---- samtools flagstat -------------------------------------------------------

/// Parse `samtools flagstat` into a [`ReadMetrics`] — **scalar counts only**. flagstat carries no
/// read-length / insert-size distributions or mapping quality, so those stay 0 (an alternative
/// `ReadMetrics` source when `stats.txt` is absent). Lines are `<n> + <qc_failed> <category> [(…)]`;
/// the first number is the QC-passed count, the category is the text before any `(`.
pub fn parse_flagstat(text: &str) -> ReadMetrics {
    let mut cats: Vec<(String, u64)> = Vec::new();
    for line in text.lines() {
        // Split on the first '+' (a later '+' appears inside "(QC-passed + QC-failed)").
        let Some((n_str, rest)) = line.split_once('+') else {
            continue;
        };
        let Ok(n) = n_str.trim().parse::<u64>() else { continue };
        // rest = "<qc_failed> <category> (pct…)" — drop the qc-failed number, strip the "(…)" tail.
        let category = rest
            .trim()
            .split_once(char::is_whitespace)
            .map(|(_, after)| after)
            .unwrap_or("")
            .split('(')
            .next()
            .unwrap_or("")
            .trim()
            .to_ascii_lowercase();
        cats.push((category, n));
    }
    let find = |needle: &str| cats.iter().find(|(c, _)| c == needle).map(|&(_, n)| n).unwrap_or(0);

    let total = find("in total");
    let mapped = find("mapped");
    let proper_pairs = find("properly paired");
    let in_pairs = find("with itself and mate mapped");
    ReadMetrics {
        total_reads: total,
        pf_reads: total, // the QC-passed count is the PF set
        pf_reads_aligned: mapped,
        reads_aligned_in_pairs: in_pairs,
        proper_pairs,
        pct_pf_reads_aligned: pct100(mapped, total),
        pct_reads_aligned_in_pairs: pct100(in_pairs, mapped),
        pct_proper_pairs: pct100(proper_pairs, total),
        ..Default::default()
    }
}

// ---- Picard metrics (CollectWgsMetrics / CollectAlignmentSummaryMetrics) ------

/// Parse a Picard metrics table: skip to the header line beginning with `header_key`, then read the
/// tab-separated data rows until a blank line (Picard appends a histogram section after a blank).
/// Returns `(headers, rows)`. `None` if the header isn't found.
fn parse_picard_rows(text: &str, header_key: &str) -> Option<(Vec<String>, Vec<Vec<String>>)> {
    let mut lines = text.lines();
    let header = lines.by_ref().find(|l| l.trim_start().starts_with(header_key))?;
    let keys: Vec<String> = header.trim().split('\t').map(str::to_string).collect();
    let mut rows = Vec::new();
    for l in lines {
        if l.trim().is_empty() {
            break;
        }
        rows.push(l.trim().split('\t').map(str::to_string).collect());
    }
    Some((keys, rows))
}

/// Zip a header + value row into a name→value lookup.
fn row_map<'a>(keys: &'a [String], row: &'a [String]) -> HashMap<&'a str, &'a str> {
    keys.iter()
        .map(String::as_str)
        .zip(row.iter().map(String::as_str))
        .collect()
}

/// Parse Picard `CollectWgsMetrics` → the genome-wide depth distribution of a [`CoverageResult`]
/// (mean/median/sd/MAD, the MAPQ/baseQ exclusion fractions, and the `pct_Nx` depth thresholds — the
/// fields the lite samtools-coverage path leaves at 0). Per-contig stats / histogram stay empty
/// (Picard is genome-wide); the ingest overlays this onto the lite result's contig breakdown.
/// Picard `PCT_*` are 0–1 fractions, matching `CoverageResult`'s convention. `None` if no table.
pub fn parse_wgs_metrics(text: &str) -> Option<CoverageResult> {
    let (keys, rows) = parse_picard_rows(text, "GENOME_TERRITORY")?;
    let row = rows.first()?;
    let m = row_map(&keys, row);
    let f = |k: &str| m.get(k).and_then(|v| v.parse::<f64>().ok()).unwrap_or(0.0);
    let u = |k: &str| m.get(k).and_then(|v| v.parse::<u64>().ok()).unwrap_or(0);
    Some(CoverageResult {
        genome_territory: u("GENOME_TERRITORY"),
        mean_coverage: f("MEAN_COVERAGE"),
        median_coverage: f("MEDIAN_COVERAGE"),
        sd_coverage: f("SD_COVERAGE"),
        mad_coverage: f("MAD_COVERAGE"),
        pct_exc_mapq: f("PCT_EXC_MAPQ"),
        pct_exc_baseq: f("PCT_EXC_BASEQ"),
        pct_1x: f("PCT_1X"),
        pct_5x: f("PCT_5X"),
        pct_10x: f("PCT_10X"),
        pct_15x: f("PCT_15X"),
        pct_20x: f("PCT_20X"),
        pct_25x: f("PCT_25X"),
        pct_30x: f("PCT_30X"),
        pct_40x: f("PCT_40X"),
        pct_50x: f("PCT_50X"),
        ..Default::default()
    })
}

/// Parse Picard `CollectAlignmentSummaryMetrics` → a [`ReadMetrics`] (the `PAIR` summary row,
/// else `UNPAIRED`, else the first). Counts + alignment percentages + mean read length + chimera
/// rate; read-length / insert-size histograms aren't in this metrics class, so they stay 0. Picard
/// `PCT_*` are 0–1 fractions → scaled to the `ReadMetrics` 0–100 convention. `None` if no table.
pub fn parse_alignment_summary(text: &str) -> Option<ReadMetrics> {
    let (keys, rows) = parse_picard_rows(text, "CATEGORY")?;
    let cat = keys.iter().position(|k| k == "CATEGORY")?;
    let row = rows
        .iter()
        .find(|r| r.get(cat).is_some_and(|c| c == "PAIR"))
        .or_else(|| rows.iter().find(|r| r.get(cat).is_some_and(|c| c == "UNPAIRED")))
        .or_else(|| rows.first())?;
    let m = row_map(&keys, row);
    let f = |k: &str| m.get(k).and_then(|v| v.parse::<f64>().ok());
    let u = |k: &str| m.get(k).and_then(|v| v.parse::<u64>().ok()).unwrap_or(0);
    let aligned = u("PF_READS_ALIGNED");
    let in_pairs = u("READS_ALIGNED_IN_PAIRS");
    Some(ReadMetrics {
        total_reads: u("TOTAL_READS"),
        pf_reads: u("PF_READS"),
        pf_reads_aligned: aligned,
        reads_aligned_in_pairs: in_pairs,
        proper_pairs: 0, // not a direct count in this metrics class
        pct_pf_reads_aligned: f("PCT_PF_READS_ALIGNED").unwrap_or(0.0) * 100.0,
        pct_reads_aligned_in_pairs: pct100(in_pairs, aligned),
        // Picard reports the *improper* fraction; the proper fraction is its complement.
        pct_proper_pairs: f("PCT_PF_READS_IMPROPER_PAIRS")
            .map(|i| (1.0 - i) * 100.0)
            .unwrap_or(0.0),
        mean_read_length: f("MEAN_READ_LENGTH").unwrap_or(0.0),
        pct_chimeras: f("PCT_CHIMERAS").unwrap_or(0.0) * 100.0,
        ..Default::default()
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn flagstat_parses_counts() {
        let f = "49876543 + 0 in total (QC-passed reads + QC-failed reads)\n\
                 123456 + 0 secondary\n\
                 2345678 + 0 duplicates\n\
                 48500000 + 0 mapped (97.24% : N/A)\n\
                 49674186 + 0 paired in sequencing\n\
                 46000000 + 0 properly paired (92.60% : N/A)\n\
                 47000000 + 0 with itself and mate mapped\n\
                 100000 + 0 with mate mapped to a different chr\n";
        let m = parse_flagstat(f);
        assert_eq!(m.total_reads, 49_876_543);
        assert_eq!(m.pf_reads_aligned, 48_500_000);
        assert_eq!(m.proper_pairs, 46_000_000);
        assert_eq!(m.reads_aligned_in_pairs, 47_000_000);
        assert!((m.pct_pf_reads_aligned - 97.24).abs() < 0.01); // computed, ≈ the reported %
        assert!(m.read_length_histogram.is_empty()); // flagstat carries no distributions
    }

    #[test]
    fn wgs_metrics_fills_the_depth_distribution() {
        let f = "## htsjdk.samtools.metrics.StringHeader\n\
                 GENOME_TERRITORY\tMEAN_COVERAGE\tSD_COVERAGE\tMEDIAN_COVERAGE\tMAD_COVERAGE\tPCT_EXC_MAPQ\tPCT_EXC_BASEQ\tPCT_1X\tPCT_10X\tPCT_30X\n\
                 3000000000\t30.5\t8.1\t31\t4\t0.012\t0.034\t0.991\t0.95\t0.6\n\n\
                 ## HISTOGRAM\ncoverage\thigh_quality_coverage_count\n0\t12345\n";
        let c = parse_wgs_metrics(f).unwrap();
        assert_eq!(c.genome_territory, 3_000_000_000);
        assert!((c.mean_coverage - 30.5).abs() < 1e-9);
        assert!((c.mad_coverage - 4.0).abs() < 1e-9);
        assert!((c.pct_exc_mapq - 0.012).abs() < 1e-9); // 0–1 fraction, as CoverageResult wants
        assert!((c.pct_1x - 0.991).abs() < 1e-9);
        assert!((c.pct_30x - 0.6).abs() < 1e-9);
        assert!(c.contig_coverage_stats.is_empty()); // Picard is genome-wide
    }

    #[test]
    fn alignment_summary_prefers_the_pair_row() {
        let f = "## METRICS CLASS\n\
                 CATEGORY\tTOTAL_READS\tPF_READS\tPF_READS_ALIGNED\tPCT_PF_READS_ALIGNED\tREADS_ALIGNED_IN_PAIRS\tMEAN_READ_LENGTH\tPCT_CHIMERAS\tPCT_PF_READS_IMPROPER_PAIRS\n\
                 FIRST_OF_PAIR\t100\t100\t98\t0.98\t96\t150\t0.001\t0.02\n\
                 SECOND_OF_PAIR\t100\t100\t97\t0.97\t96\t150\t0.001\t0.02\n\
                 PAIR\t200\t200\t195\t0.975\t192\t150.5\t0.0012\t0.02\n";
        let m = parse_alignment_summary(f).unwrap();
        assert_eq!(m.total_reads, 200); // the PAIR row, not FIRST/SECOND
        assert_eq!(m.pf_reads_aligned, 195);
        assert!((m.pct_pf_reads_aligned - 97.5).abs() < 1e-6); // 0.975 → 97.5%
        assert!((m.pct_proper_pairs - 98.0).abs() < 1e-6); // 1 - 0.02 → 98%
        assert!((m.mean_read_length - 150.5).abs() < 1e-9);
    }

    #[test]
    fn picard_parsers_return_none_without_a_table() {
        assert!(parse_wgs_metrics("no table here\n").is_none());
        assert!(parse_alignment_summary("nope\n").is_none());
    }

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
        assert!(
            (m.mean_insert_size - 430.6).abs() < 1e-9,
            "uses samtools' reported average"
        );
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

        let c = lite_coverage(
            &read("coverage.txt"),
            Some(&read("HG00096.chm13.chrYM.callable.summary.txt")),
        );
        eprintln!(
            "mean_cov={:.2} territory={} callable={}",
            c.mean_coverage, c.genome_territory, c.callable_bases
        );
        assert!(c.mean_coverage > 0.0 && c.callable_bases > 0);
    }
}
