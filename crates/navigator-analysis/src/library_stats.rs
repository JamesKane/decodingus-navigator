//! Library / instrument inference from an alignment — a Rust port of the Scala
//! `LibraryStatsProcessor`. The header probe ([`crate::probe`]) only reads `@RG PL/PM`, which
//! many vendor BAMs (FGC/YSEQ/Dante…) leave sparse; the instrument serial that identifies the
//! physical sequencer (and, via the crowd-sourced AppView map, the lab) lives in the **read
//! names**. This scans a bounded prefix of records, classifies each read's platform from its
//! qname, extracts the instrument + flowcell, and reports the most-frequent of each plus the
//! `@RG SM/LB/PU` tags.
//!
//! The `instrument_id` is the key datum for resolving the sequencing facility (roadmap D8).

use std::collections::HashMap;
use std::path::Path;

use noodles::sam::header::record::value::map::read_group;

use crate::error::AnalysisError;
use crate::reader::open_seq;

/// Default cap on records scanned — enough to settle the most-frequent instrument without
/// reading a whole multi-GB file.
pub const DEFAULT_MAX_READS: usize = 10_000;

/// What the read-name scan + `@RG` tags inferred. Every field is best-effort (`None`/`Unknown`).
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct LibraryStats {
    /// Primary (non-secondary/supplementary) records scanned.
    pub read_count: usize,
    /// `@RG SM` (sample name).
    pub sample_name: Option<String>,
    /// `@RG LB` (library id).
    pub library_id: Option<String>,
    /// `@RG PU` (platform unit).
    pub platform_unit: Option<String>,
    /// Most-frequent platform inferred from read names: `Illumina`/`PacBio`/`MGI`/`Nanopore`.
    pub platform: Option<String>,
    /// Most-frequent instrument serial from read names (e.g. `A00123`, `m84005`).
    pub instrument_id: Option<String>,
    /// Specific instrument model inferred from the platform + instrument-id prefix.
    pub instrument_model: Option<String>,
    /// Most-frequent flowcell id from read names.
    pub flowcell_id: Option<String>,
    /// `PAIRED` if the scanned reads carry the SAM segmented (0x1) flag, else `SINGLE`; `None`
    /// when no primary reads were scanned.
    pub library_layout: Option<String>,
}

/// Scan up to `max_reads` records of `path` (CRAM needs `reference`) and infer library/instrument
/// metadata from the `@RG` header tags + read-name patterns.
pub fn scan_library_stats(
    path: &Path,
    reference: Option<&Path>,
    max_reads: usize,
) -> Result<LibraryStats, AnalysisError> {
    let (header, reader) = open_seq(path, reference)?;

    // @RG SM/LB/PU from the first informative read group (stable across re-alignments).
    let (mut sample_name, mut library_id, mut platform_unit) = (None, None, None);
    for map in header.read_groups().values() {
        let sm = map.other_fields().get(&read_group::tag::SAMPLE).map(val);
        let lb = map.other_fields().get(&read_group::tag::LIBRARY).map(val);
        let pu = map.other_fields().get(&read_group::tag::PLATFORM_UNIT).map(val).filter(|s| !s.is_empty());
        if sm.is_some() || lb.is_some() || pu.is_some() {
            sample_name = sm;
            library_id = lb;
            platform_unit = pu;
            break;
        }
    }

    let mut read_count = 0usize;
    let mut segmented_count = 0usize;
    let mut platform_counts: HashMap<&'static str, usize> = HashMap::new();
    let mut instruments: HashMap<String, usize> = HashMap::new();
    let mut flowcells: HashMap<String, usize> = HashMap::new();

    let mut reader = reader;
    for rec in reader.records(&header) {
        if read_count >= max_reads {
            break;
        }
        let record = rec?;
        let flags = record.flags();
        if flags.is_secondary() || flags.is_supplementary() {
            continue;
        }
        read_count += 1;
        if flags.is_segmented() {
            segmented_count += 1; // SAM 0x1 = part of a paired/multi-segment template
        }
        let Some(qname) = record.name().map(|n| n.to_string()) else { continue };

        let platform = detect_platform_from_qname(&qname);
        *platform_counts.entry(platform).or_insert(0) += 1;
        let (instrument, flowcell) = parse_instrument_and_flowcell(&qname, platform);
        if let Some(i) = instrument {
            *instruments.entry(i).or_insert(0) += 1;
        }
        if let Some(f) = flowcell {
            *flowcells.entry(f).or_insert(0) += 1;
        }
    }

    let primary_platform = most_frequent(&platform_counts).map(|p| p.to_string());
    let instrument_id = most_frequent(&instruments);
    let flowcell_id = most_frequent(&flowcells);
    let instrument_model = match (primary_platform.as_deref(), instrument_id.as_deref()) {
        (Some(p), Some(i)) => Some(infer_model(p, i)),
        _ => None,
    };

    // Majority vote: PAIRED if most scanned reads are segmented (handles a few stray flags).
    let library_layout = (read_count > 0)
        .then(|| if segmented_count * 2 >= read_count { "PAIRED" } else { "SINGLE" }.to_string());

    Ok(LibraryStats {
        read_count,
        sample_name,
        library_id,
        platform_unit,
        platform: primary_platform,
        instrument_id,
        instrument_model,
        flowcell_id,
        library_layout,
    })
}

/// The key with the highest count (ties → arbitrary stable pick).
fn most_frequent<K: Clone>(counts: &HashMap<K, usize>) -> Option<K> {
    counts.iter().max_by_key(|(_, n)| **n).map(|(k, _)| k.clone())
}

fn val<T: AsRef<[u8]>>(v: &T) -> String {
    String::from_utf8_lossy(v.as_ref()).into_owned()
}

/// Classify a read's platform from its name. Mirrors the Scala heuristics (MGI prefixes /
/// colon-shaped Illumina / UUID Nanopore / `m#####` PacBio); `Unknown` when nothing matches.
fn detect_platform_from_qname(qname: &str) -> &'static str {
    // MGI: a V300/E100/CL100/G400/G99 instrument prefix, or a colon-delimited name whose first
    // field starts V/E/CL/G with a flowcell field starting "L".
    if qname.len() > 15 {
        let prefix = qname.get(0..5).unwrap_or("").to_ascii_uppercase();
        if ["V300", "E100", "CL100", "G400", "G99"].iter().any(|p| prefix.starts_with(p)) {
            return "MGI";
        }
        if qname.matches(':').count() >= 6 {
            let parts: Vec<&str> = qname.split(':').collect();
            let p0 = parts[0];
            if (p0.starts_with('V') || p0.starts_with('E') || p0.starts_with("CL") || p0.starts_with('G'))
                && parts.len() >= 3
                && parts[2].starts_with('L')
            {
                return "MGI";
            }
        }
    }
    if is_illumina_qname(qname) {
        "Illumina"
    } else if is_nanopore_uuid(qname) {
        "Nanopore"
    } else if is_pacbio_qname(qname) {
        "PacBio"
    } else {
        "Unknown"
    }
}

/// Old Casava (`…:N:N:N:N#…`) or modern (`inst:run:flowcell:lane:tile:x:y`) Illumina read names.
fn is_illumina_qname(q: &str) -> bool {
    let num = |s: &str| !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit());
    let alnum = |s: &str| !s.is_empty() && s.bytes().all(|b| b.is_ascii_alphanumeric());
    // Casava 1.7: the four colon fields immediately before '#' are all numeric (lane:tile:x:y).
    if let Some(hash) = q.find('#') {
        let tail: Vec<&str> = q[..hash].rsplit(':').take(4).collect();
        if tail.len() == 4 && tail.iter().all(|f| num(f)) {
            return true;
        }
    }
    // Modern: a window [num, alnum, single-digit, num, num, num] = run:flowcell:lane:tile:x:y.
    let f: Vec<&str> = q.split(':').collect();
    f.windows(6).any(|w| {
        num(w[0]) && alnum(w[1]) && w[2].len() == 1 && num(w[2]) && num(w[3]) && num(w[4]) && num(w[5])
    })
}

/// A leading UUID (`8-4-4-4-12` hex) — Oxford Nanopore.
fn is_nanopore_uuid(q: &str) -> bool {
    let b = q.as_bytes();
    if b.len() < 36 {
        return false;
    }
    let hex = |i: usize| b[i].is_ascii_hexdigit();
    let dash = |i: usize| b[i] == b'-';
    (0..8).all(hex) && dash(8) && (9..13).all(hex) && dash(13) && (14..18).all(hex) && dash(18)
        && (19..23).all(hex) && dash(23) && (24..36).all(hex)
}

/// `m` then ≥5 digits (PacBio movie name, e.g. `m84005_…`).
fn is_pacbio_qname(q: &str) -> bool {
    let b = q.as_bytes();
    b.first() == Some(&b'm') && b[1..].iter().take_while(|c| c.is_ascii_digit()).count() >= 5
}

/// Extract `(instrument, flowcell)` from a read name given its platform.
fn parse_instrument_and_flowcell(qname: &str, platform: &str) -> (Option<String>, Option<String>) {
    match platform {
        // inst:run:flowcell:lane:… → instrument=field0, flowcell=field2.
        "Illumina" => {
            let parts: Vec<&str> = qname.split(':').collect();
            if parts.len() >= 3 {
                (Some(parts[0].to_string()), Some(parts[2].to_string()))
            } else {
                (None, None)
            }
        }
        // movie/zmw → instrument = the part before the first '_' of field 0.
        "PacBio" => {
            let first = qname.split('/').next().unwrap_or("");
            let inst = first.split('_').next().unwrap_or("");
            if inst.is_empty() { (None, None) } else { (Some(inst.to_string()), None) }
        }
        "MGI" => {
            if qname.matches(':').count() >= 3 {
                let parts: Vec<&str> = qname.split(':').collect();
                (Some(parts[0].to_string()), Some(parts[1].to_string()))
            } else if qname.len() > 10 {
                // Concatenated form: <instrument>L<lane>C<col>R<row>…
                if let Some(lpos) = qname.find('L') {
                    if lpos > 0 {
                        let instrument = &qname[..lpos];
                        let rest = &qname[lpos..];
                        if rest.contains('C') {
                            let end = rest.find('R').unwrap_or(rest.len());
                            return (Some(instrument.to_string()), Some(qname[lpos..lpos + end].to_string()));
                        }
                    }
                }
                (None, None)
            } else {
                (None, None)
            }
        }
        _ => (None, None),
    }
}

/// Map a platform + instrument-id prefix to a specific instrument model (Scala `inferPlatform`).
fn infer_model(platform: &str, instrument_id: &str) -> String {
    match platform {
        "Illumina" => match instrument_id.chars().next().map(|c| c.to_ascii_uppercase()) {
            Some('A') => "NovaSeq",
            Some('D') => "HiSeq 2500",
            Some('J') => "HiSeq 3000",
            Some('K') => "HiSeq 4000",
            Some('E') => "HiSeq X",
            Some('N') => "NextSeq",
            Some('M') => "MiSeq",
            Some('V') => "NovaSeq X",
            Some('F') => "iSeq",
            _ => "Unknown Illumina",
        }
        .to_string(),
        "PacBio" => {
            if instrument_id.starts_with("m84") {
                "PacBio Revio"
            } else if instrument_id.starts_with("m64") {
                "PacBio Sequel II/IIe"
            } else if instrument_id.starts_with("m54") {
                "PacBio Sequel"
            } else {
                "PacBio"
            }
        }
        .to_string(),
        "MGI" => {
            if instrument_id.starts_with("V300") {
                "MGI DNBSEQ/MGISEQ-2000"
            } else if instrument_id.starts_with("E100") {
                "MGI MGISEQ-200"
            } else if instrument_id.starts_with("CL100") {
                "MGI MGISEQ-T7"
            } else if instrument_id.starts_with("G400") {
                "MGI DNBSEQ-G400"
            } else if instrument_id.starts_with("G99") {
                "MGI MGISEQ-T1"
            } else {
                "MGI DNBseq"
            }
        }
        .to_string(),
        "Nanopore" => "Oxford Nanopore".to_string(),
        _ => "Unknown".to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classifies_read_names_by_platform() {
        // Modern Illumina: instrument:run:flowcell:lane:tile:x:y
        assert_eq!(detect_platform_from_qname("A00123:45:H7TJ2DSXX:1:1101:1000:1996"), "Illumina");
        // Old Casava with '#index'
        assert_eq!(detect_platform_from_qname("HWUSI-EAS100R:6:73:941:1973#0/1"), "Illumina");
        // PacBio movie
        assert_eq!(detect_platform_from_qname("m84005_230101_000000/1234/ccs"), "PacBio");
        // Nanopore UUID
        assert_eq!(detect_platform_from_qname("abcdef01-2345-6789-abcd-ef0123456789"), "Nanopore");
        // MGI
        assert_eq!(detect_platform_from_qname("V300012345L1C001R0010000123"), "MGI");
        assert_eq!(detect_platform_from_qname("totally random name"), "Unknown");
    }

    #[test]
    fn extracts_instrument_and_flowcell() {
        assert_eq!(
            parse_instrument_and_flowcell("A00123:45:H7TJ2DSXX:1:1101:1000:1996", "Illumina"),
            (Some("A00123".into()), Some("H7TJ2DSXX".into()))
        );
        assert_eq!(
            parse_instrument_and_flowcell("m84005_230101_000000/1234/ccs", "PacBio"),
            (Some("m84005".into()), None)
        );
    }

    #[test]
    fn library_layout_majority_vote() {
        // Helper mirrors the in-scan rule: PAIRED iff segmented*2 >= total.
        let layout = |segmented: usize, total: usize| {
            (total > 0).then(|| if segmented * 2 >= total { "PAIRED" } else { "SINGLE" }.to_string())
        };
        assert_eq!(layout(0, 0), None); // nothing scanned
        assert_eq!(layout(100, 100).as_deref(), Some("PAIRED"));
        assert_eq!(layout(0, 100).as_deref(), Some("SINGLE"));
        assert_eq!(layout(60, 100).as_deref(), Some("PAIRED")); // majority paired
        assert_eq!(layout(40, 100).as_deref(), Some("SINGLE"));
    }

    #[test]
    fn infers_specific_models_from_instrument_prefix() {
        assert_eq!(infer_model("Illumina", "A00123"), "NovaSeq");
        assert_eq!(infer_model("Illumina", "M01234"), "MiSeq");
        assert_eq!(infer_model("PacBio", "m84005"), "PacBio Revio");
        assert_eq!(infer_model("PacBio", "m64012"), "PacBio Sequel II/IIe");
        assert_eq!(infer_model("MGI", "V300012345"), "MGI DNBSEQ/MGISEQ-2000");
    }
}
