//! Read called haploid bases at target positions from a GATK ploidy-1 GVCF.
//!
//! The `ytree` pipeline archives a per-sample chrY/chrM GVCF (HaplotypeCaller
//! `--sample-ploidy 1 -ERC GVCF`) next to each CRAM. Those GVCFs already contain exactly
//! what [`crate::caller::call_bases_at`] would recompute by walking the (multi-GB) CRAM at
//! every haplotree position — the *observed haploid base* at each site. Reading the small
//! GVCF instead of the CRAM is the fast path for haplogroup placement.
//!
//! GVCF semantics (ploidy 1, `<NON_REF>` model):
//! - **Variant record** (`ALT` has a real allele besides `<NON_REF>`):
//!   `GT=1` → the sample carries the ALT. SNP → that base is the *derived* observation;
//!   indel → confident but not a usable SNP base (skipped). `GT=0` → confident hom-ref
//!   (an *ancestral* observation at a multiallelic emit site).
//! - **Ref block** (`ALT=<NON_REF>`, `GT=0`, `END=` in INFO): every position in `[POS,END]`
//!   was called hom-ref (ancestral) at the block's confidence.
//!
//! This module decodes the GVCF to two facts per target: the *derived base* where one was
//! called, and whether the site was *callable* at all. [`assemble_calls`] then turns those
//! into the `position → observed base` map [`crate::haplo::score`] consumes — using the
//! tree's ancestral allele for callable-but-not-variant sites (on the native build the
//! reference base == the tree's ancestral allele, so no FASTA lookup is needed).

use std::collections::{HashMap, HashSet};
use std::io::BufRead;
use std::path::Path;

use noodles::bgzf;

use crate::error::AnalysisError;

/// Confidence thresholds for trusting a GVCF call. Ref blocks are gated on `MIN_DP`
/// (falling back to `DP`) and `GQ`; variant records on `DP`/`GQ`. Defaults are permissive
/// enough for low-coverage HiFi (the pipeline's ref blocks carry GQ 70–99) while rejecting
/// genuinely unsupported sites.
#[derive(Debug, Clone, Copy)]
pub struct GvcfReadParams {
    pub min_dp: u32,
    pub min_gq: u32,
}

impl Default for GvcfReadParams {
    fn default() -> Self {
        Self { min_dp: 2, min_gq: 20 }
    }
}

/// The two facts decoded per target position from the GVCF.
#[derive(Debug, Clone, Default)]
pub struct CalledBases {
    /// SNP-derived ALT base (uppercase) at target sites where the sample carries a
    /// single-base ALT (`GT=1`).
    pub variant_bases: HashMap<i64, char>,
    /// Target positions the GVCF confidently called (a passing ref block, a hom-ref
    /// variant emit, or a confident SNP) — i.e. *not* a no-call.
    pub callable: HashSet<i64>,
}

/// Read called bases at `targets` on `contig` from a bgzipped GVCF on disk.
pub fn read_called_bases(
    gvcf: &Path,
    contig: &str,
    targets: &HashSet<i64>,
    params: &GvcfReadParams,
) -> Result<CalledBases, AnalysisError> {
    let file = std::fs::File::open(gvcf).map_err(|e| AnalysisError::io(gvcf, e))?;
    read_called_bases_from(bgzf::Reader::new(file), contig, targets, params)
}

/// Decode core over any `BufRead` (plain-text VCF in tests). Streams the whole file —
/// these GVCFs are small (chrY ~3 MB, chrM ~6 KB) and the targets are a few thousand
/// scattered positions, so a single linear pass beats per-target tabix seeks.
pub fn read_called_bases_from<R: BufRead>(
    mut reader: R,
    contig: &str,
    targets: &HashSet<i64>,
    params: &GvcfReadParams,
) -> Result<CalledBases, AnalysisError> {
    // Sorted targets so a ref block's [POS, END] span resolves by binary search instead of
    // iterating the (potentially thousands-wide) block.
    let mut sorted: Vec<i64> = targets.iter().copied().collect();
    sorted.sort_unstable();

    let mut out = CalledBases::default();
    let mut line = String::new();
    loop {
        line.clear();
        let n = reader
            .read_line(&mut line)
            .map_err(|e| AnalysisError::Message(format!("reading gvcf: {e}")))?;
        if n == 0 {
            break;
        }
        let l = line.trim_end_matches(['\n', '\r']);
        if l.is_empty() || l.starts_with('#') {
            continue;
        }

        let mut col = l.split('\t');
        let chrom = col.next().unwrap_or("");
        if chrom != contig {
            continue;
        }
        let pos: i64 = match col.next().and_then(|s| s.parse().ok()) {
            Some(p) => p,
            None => continue,
        };
        let _id = col.next();
        let refa = col.next().unwrap_or("");
        let alt = col.next().unwrap_or("");
        let _qual = col.next();
        let _filter = col.next();
        let info = col.next().unwrap_or("");
        let format = col.next().unwrap_or("");
        let sample = col.next().unwrap_or("");

        // Haploid GT is the first ':' field; take the first allele (no '/'|'|' for ploidy 1,
        // but be defensive). '.' / missing → skip.
        let gt = format_field(format, sample, "GT").unwrap_or("");
        let allele = gt.split(['/', '|']).next().unwrap_or(gt);
        if allele.is_empty() || allele == "." {
            continue;
        }

        if alt == "<NON_REF>" {
            // Ref block: confident hom-ref over [POS, END].
            if allele != "0" {
                continue;
            }
            let dp = format_field(format, sample, "MIN_DP")
                .or_else(|| format_field(format, sample, "DP"))
                .and_then(|s| s.parse::<u32>().ok())
                .unwrap_or(0);
            let gq = format_field(format, sample, "GQ")
                .and_then(|s| s.parse::<u32>().ok())
                .unwrap_or(0);
            if dp < params.min_dp || gq < params.min_gq {
                continue;
            }
            let end = info_end(info).unwrap_or(pos);
            for &t in targets_in_range(&sorted, pos, end) {
                out.callable.insert(t);
            }
        } else {
            // Variant record. Only target positions matter.
            if !targets.contains(&pos) {
                continue;
            }
            let dp = format_field(format, sample, "DP")
                .and_then(|s| s.parse::<u32>().ok())
                .unwrap_or(0);
            let gq = format_field(format, sample, "GQ")
                .and_then(|s| s.parse::<u32>().ok())
                .unwrap_or(0);
            if dp < params.min_dp || gq < params.min_gq {
                continue;
            }
            if allele == "0" {
                // Confident hom-ref at a multiallelic emit site → ancestral observation.
                out.callable.insert(pos);
                continue;
            }
            // GT carries an ALT. First real ALT (skip the trailing <NON_REF>).
            let alt0 = alt.split(',').find(|a| *a != "<NON_REF>").unwrap_or("");
            if refa.len() == 1 && alt0.len() == 1 {
                let b = alt0.as_bytes()[0].to_ascii_uppercase();
                if matches!(b, b'A' | b'C' | b'G' | b'T') {
                    out.variant_bases.insert(pos, b as char);
                    out.callable.insert(pos);
                }
            }
            // An indel ALT at a (SNP) tree position is left as a no-call rather than
            // asserted ancestral — conservative; avoids a false ancestral refutation.
        }
    }
    Ok(out)
}

/// Assemble the `position → observed base` map [`crate::haplo::score`] consumes from the
/// decoded GVCF facts. A variant (derived) base wins; otherwise a callable hom-ref site takes
/// the **reference genome base** at that position (`ref_base`); otherwise the position is a
/// no-call and is omitted.
///
/// `ref_base` must be the *reference* base, not the tree ancestral — the two differ wherever
/// the reference itself carries a derived allele. CHM13's Y is HG002 (haplogroup J1, deep in
/// the tree), so at every backbone SNP shared by J1 and the sample the GVCF emits a ref block
/// (hom-ref == reference == *derived*), and assuming ancestral there would silently break the
/// descent. This mirrors [`crate::caller::call_bases_at`], which reads the actual base off the
/// reads (== the reference base at a hom-ref site).
pub fn assemble_calls(called: &CalledBases, ref_base: &HashMap<i64, char>) -> HashMap<i64, char> {
    let mut calls: HashMap<i64, char> = HashMap::with_capacity(called.callable.len());
    for &pos in &called.callable {
        if let Some(&r) = ref_base.get(&pos) {
            calls.insert(pos, r);
        }
    }
    // Variant (derived) observations override the reference default.
    for (&pos, &base) in &called.variant_bases {
        calls.insert(pos, base);
    }
    calls
}

/// Value of `key` in a `FORMAT`/sample colon-delimited pair (e.g. `GT:DP:GQ` + `1:18:99`).
fn format_field<'a>(format: &str, sample: &'a str, key: &str) -> Option<&'a str> {
    let idx = format.split(':').position(|k| k == key)?;
    sample.split(':').nth(idx)
}

/// `END=` value from a GVCF ref block's INFO column, if present.
fn info_end(info: &str) -> Option<i64> {
    info.split(';')
        .find_map(|kv| kv.strip_prefix("END="))
        .and_then(|v| v.parse().ok())
}

/// Slice of `sorted` whose values fall in `[lo, hi]` (inclusive).
fn targets_in_range(sorted: &[i64], lo: i64, hi: i64) -> &[i64] {
    let start = sorted.partition_point(|&t| t < lo);
    let end = sorted.partition_point(|&t| t <= hi);
    &sorted[start..end]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn targets(ps: &[i64]) -> HashSet<i64> {
        ps.iter().copied().collect()
    }

    const SAMPLE_GVCF: &str = "\
##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG00096
chrY\t2458321\t.\tG\t<NON_REF>\t.\t.\tEND=2459920\tGT:DP:GQ:MIN_DP:PL\t0:17:99:9:0,288
chrY\t2459921\t.\tG\tA,<NON_REF>\t606\t.\tDP=19\tGT:AD:DP:GQ:PL\t1:0,18,0:18:99:616,0,616
chrY\t2477255\t.\tC\tT,<NON_REF>\t685\t.\tDP=22\tGT:AD:DP:GQ:PL\t1:0,20,0:20:99:695,0,695
chrY\t2481534\t.\tA\tAT,<NON_REF>\t379\t.\tDP=18\tGT:AD:DP:GQ:PL\t1:0,15,0:15:99:389,0,389
chrM\t100\t.\tC\tT,<NON_REF>\t500\t.\tDP=30\tGT:AD:DP:GQ:PL\t1:0,30,0:30:99:510,0,510
";

    #[test]
    fn snp_variant_yields_derived_base_and_callable() {
        let t = targets(&[2459921, 2477255]);
        let c = read_called_bases_from(SAMPLE_GVCF.as_bytes(), "chrY", &t, &GvcfReadParams::default()).unwrap();
        assert_eq!(c.variant_bases.get(&2459921), Some(&'A'));
        assert_eq!(c.variant_bases.get(&2477255), Some(&'T'));
        assert!(c.callable.contains(&2459921));
        assert!(c.callable.contains(&2477255));
    }

    #[test]
    fn ref_block_marks_spanned_targets_callable_only() {
        // 2459000 falls inside the 2458321..2459920 ref block; not a variant → callable, no base.
        let t = targets(&[2459000]);
        let c = read_called_bases_from(SAMPLE_GVCF.as_bytes(), "chrY", &t, &GvcfReadParams::default()).unwrap();
        assert!(c.callable.contains(&2459000));
        assert!(!c.variant_bases.contains_key(&2459000));
    }

    #[test]
    fn indel_variant_is_no_call_not_ancestral() {
        // 2481534 is an insertion (A>AT) → neither a usable SNP base nor a forced ancestral.
        let t = targets(&[2481534]);
        let c = read_called_bases_from(SAMPLE_GVCF.as_bytes(), "chrY", &t, &GvcfReadParams::default()).unwrap();
        assert!(!c.variant_bases.contains_key(&2481534));
        assert!(!c.callable.contains(&2481534));
    }

    #[test]
    fn contig_filter_excludes_other_chromosomes() {
        let t = targets(&[100]);
        let c = read_called_bases_from(SAMPLE_GVCF.as_bytes(), "chrY", &t, &GvcfReadParams::default()).unwrap();
        assert!(c.variant_bases.is_empty(), "chrM:100 must not leak into a chrY read");
        let c = read_called_bases_from(SAMPLE_GVCF.as_bytes(), "chrM", &t, &GvcfReadParams::default()).unwrap();
        assert_eq!(c.variant_bases.get(&100), Some(&'T'));
    }

    #[test]
    fn out_of_region_target_is_no_call() {
        let t = targets(&[99_000_000]);
        let c = read_called_bases_from(SAMPLE_GVCF.as_bytes(), "chrY", &t, &GvcfReadParams::default()).unwrap();
        assert!(c.callable.is_empty() && c.variant_bases.is_empty());
    }

    #[test]
    fn dp_gq_thresholds_reject_weak_calls() {
        let strict = GvcfReadParams { min_dp: 25, min_gq: 20 };
        let t = targets(&[2459921]); // DP=18 < 25
        let c = read_called_bases_from(SAMPLE_GVCF.as_bytes(), "chrY", &t, &strict).unwrap();
        assert!(c.variant_bases.is_empty(), "DP below threshold must be rejected");
    }

    #[test]
    fn assemble_prefers_variant_then_reference() {
        let mut called = CalledBases::default();
        called.variant_bases.insert(2459921, 'A');
        called.callable.insert(2459921);
        called.callable.insert(2459000); // hom-ref-only → takes the reference base
                                         // The reference base at a hom-ref site can be the *derived* allele (CHM13 = J1 Y).
        let ref_base: HashMap<i64, char> = [(2459921, 'G'), (2459000, 'T'), (700, 'C')].into_iter().collect();
        let calls = assemble_calls(&called, &ref_base);
        assert_eq!(calls.get(&2459921), Some(&'A'), "variant (derived) wins over reference");
        assert_eq!(
            calls.get(&2459000),
            Some(&'T'),
            "callable hom-ref takes the reference base"
        );
        assert!(!calls.contains_key(&700), "no-call position is omitted");
    }

    #[test]
    fn targets_in_range_is_inclusive() {
        let s = [10, 20, 30, 40];
        assert_eq!(targets_in_range(&s, 20, 30), &[20, 30]);
        assert_eq!(targets_in_range(&s, 21, 39), &[30]);
        assert_eq!(targets_in_range(&s, 0, 5), &[] as &[i64]);
        assert_eq!(targets_in_range(&s, 40, 100), &[40]);
    }

    /// Real-data smoke test: decode the pipeline's actual bgzipped chrY GVCF for HG00096
    /// over a dense synthetic target grid across the non-PAR span. Validates real bgzf
    /// inflation + record parsing at scale (thousands of records). No-ops when the NAS
    /// file isn't mounted, so it's safe on any machine. Run with:
    ///   cargo test -p navigator-analysis gvcf -- --ignored --nocapture
    #[test]
    #[ignore = "reads a NAS file; run explicitly"]
    fn real_chr_y_gvcf_decodes() {
        let path = Path::new("/Volumes/nas/Genomics/PRJEB31736/HG00096/HG00096.chm13.chrY.g.vcf.gz");
        if !path.exists() {
            eprintln!("skip: {} not mounted", path.display());
            return;
        }
        // Every 50th base across the non-PAR region — a stand-in for tree positions.
        let t: HashSet<i64> = (2_458_321..62_122_809).step_by(50).collect();
        let c = read_called_bases(path, "chrY", &t, &GvcfReadParams::default()).unwrap();
        eprintln!(
            "HG00096 chrY: {} targets, {} variant bases, {} callable",
            t.len(),
            c.variant_bases.len(),
            c.callable.len()
        );
        assert!(
            c.callable.len() > 1000,
            "expected many callable sites on a real chrY GVCF"
        );
        assert!(c.variant_bases.values().all(|b| matches!(b, 'A' | 'C' | 'G' | 'T')));
    }
}
