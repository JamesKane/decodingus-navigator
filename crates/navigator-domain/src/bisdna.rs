//! BISDNA chromo2 Y-chromosome raw-data parsing. The export is a tab-delimited table —
//! `SNPID`, an Illumina TOP-strand `genotype`, and a `result` verdict (positive/negative/
//! no_call/back-mutated) — preceded by a multi-line prose preamble. Crucially it carries
//! **no positions or alleles**: a SNP name plus a derived/ancestral verdict. Turning those
//! into placeable variant calls needs an external name→locus dictionary (see the design
//! `docs/design/bisdna-import.md`); this module is only the faithful, IO-free file parse.
//!
//! Strand note: the genotype is on the Illumina TOP strand, which need not match the
//! reference + strand, so it is kept verbatim and is *not* the source of truth for
//! derived/ancestral — the `result` column is (see [`Verdict`]).

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::variants::{self, VariantCall};
use crate::ysnp_dict::YsnpDictionary;

/// The `result`-column verdict for one marker.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Verdict {
    /// Derived allele carried (includes the parenthetical `(positive)` form BISDNA uses for
    /// a positive call on a back-mutation-prone marker).
    Positive,
    /// Ancestral allele carried.
    Negative,
    /// Undetermined — the genotype is `00` and BISDNA could not call the marker.
    NoCall,
    /// The lineage is derived but the base reads ancestral (a documented back-mutation, e.g.
    /// S163). The placement layer flags and excludes these — a position→base call can't
    /// represent "derived lineage showing the ancestral base".
    BackMutated,
}

/// One parsed BISDNA row: the SNP name, its raw (TOP-strand) genotype, and the verdict.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct BisdnaCall {
    /// SNP name as exported (e.g. `CTS10003`, `Apt`). Resolved to a locus later.
    pub name: String,
    /// Raw Illumina TOP-strand genotype, doubled for haploid Y calls (`GG`, `AG`, `00`).
    pub genotype: String,
    pub verdict: Verdict,
}

/// Trim whitespace and one layer of surrounding double-quotes from a cell.
fn clean_cell(s: &str) -> &str {
    s.trim().trim_matches('"').trim()
}

/// Map a `result`-column token to a [`Verdict`], case- and whitespace-insensitive. Returns
/// `None` for an unrecognized token (a malformed row, skipped by [`parse`]).
fn parse_verdict(token: &str) -> Option<Verdict> {
    let t = token.trim().trim_matches(|c| c == '"').to_ascii_lowercase();
    // Strip the parentheses BISDNA wraps a back-mutation-prone positive in: `(positive)`.
    let core = t.trim_matches(|c| c == '(' || c == ')').trim();
    match core {
        "positive" => Some(Verdict::Positive),
        "negative" => Some(Verdict::Negative),
        "no_call" | "nocall" | "no call" => Some(Verdict::NoCall),
        "back-mutated" | "back_mutated" | "backmutated" | "back mutated" => Some(Verdict::BackMutated),
        _ => None,
    }
}

/// Is this tab-split row the `SNPID / genotype / result` header? (Case-insensitive; the prose
/// preamble lines are single-column and never match.)
fn is_header(cols: &[&str]) -> bool {
    cols.len() >= 3
        && clean_cell(cols[0]).eq_ignore_ascii_case("snpid")
        && clean_cell(cols[1]).eq_ignore_ascii_case("genotype")
        && clean_cell(cols[2]).eq_ignore_ascii_case("result")
}

/// Parse a BISDNA chromo2 export into calls. Skips the prose preamble by seeking the
/// `SNPID<TAB>genotype<TAB>result` header, then reads each tab-delimited data row
/// (`name`, `genotype`, `result`). Blank lines and rows with an unrecognized verdict are
/// skipped; every recognized row is kept verbatim (including `NoCall`/`BackMutated` — the
/// importer, not the parser, decides what to drop). Errors only if the header is missing or
/// no data rows follow it.
pub fn parse(text: &str) -> Result<Vec<BisdnaCall>, String> {
    let mut lines = text.lines();

    // Seek the header, skipping the multi-line prose preamble.
    let header_found = lines.by_ref().any(|line| {
        let cols: Vec<&str> = line.split('\t').collect();
        is_header(&cols)
    });
    if !header_found {
        return Err("not a BISDNA export (missing `SNPID\\tgenotype\\tresult` header)".into());
    }

    let mut calls = Vec::new();
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 3 {
            continue;
        }
        let name = clean_cell(cols[0]);
        let genotype = clean_cell(cols[1]);
        let Some(verdict) = parse_verdict(cols[2]) else {
            continue;
        };
        if name.is_empty() {
            continue;
        }
        calls.push(BisdnaCall { name: name.to_string(), genotype: genotype.to_string(), verdict });
    }

    if calls.is_empty() {
        return Err("BISDNA header found but no marker rows followed".into());
    }
    Ok(calls)
}

/// Watson–Crick complement (non-ACGT passes through). Local copy so the resolver stays in
/// `navigator-domain`; the genotype QC needs it because BISDNA calls the Illumina TOP strand.
fn complement(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other,
    }
}

/// Does `genotype` carry `allele` (or its complement)? QC only — a miss on both strands flags
/// a likely dictionary/name mismatch, but the verdict (not the genotype) decides the call.
fn genotype_supports(genotype: &str, allele: &str) -> bool {
    let Some(want) = allele.bytes().next().map(|b| b.to_ascii_uppercase()) else {
        return true;
    };
    let comp = complement(want);
    genotype.bytes().map(|b| b.to_ascii_uppercase()).any(|b| b == want || b == comp)
}

/// The result of resolving BISDNA calls against the Y-SNP dictionary on a given build: the
/// emitted variant calls (positives only) plus a per-category tally.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct ResolveOutcome {
    /// Positive (derived) calls resolved to a locus, as carried `VariantCall`s.
    pub calls: Vec<VariantCall>,
    /// Negative (ancestral) markers — not variants, so not emitted.
    pub ancestral: usize,
    /// `no_call` markers.
    pub no_call: usize,
    /// Back-mutated markers — flagged, excluded from placement.
    pub back_mutated: usize,
    /// Positive markers whose name the dictionary couldn't place on this build.
    pub unresolved: usize,
    /// A capped sample of unresolved names (for diagnostics).
    pub unresolved_names: Vec<String>,
    /// Positive calls whose genotype disagreed with the dictionary alleles on either strand.
    pub strand_mismatches: usize,
}

/// Resolve parsed BISDNA `calls` to carried Y-SNP variant calls on `build`, using `dict` for
/// name→locus. Only **positive** (derived) markers are emitted (`reference` = ancestral,
/// `alternate` = derived, genotype `"1"`); a negative is not a variant, and the variant-level
/// reconciler weights every stored call as a carried allele. Negative/no_call/back-mutated and
/// dictionary-unresolved markers are tallied, not emitted. `unresolved_cap` bounds the sample
/// of unresolved names kept. Pure — no IO.
pub fn resolve_calls(
    calls: &[BisdnaCall],
    dict: &YsnpDictionary,
    build: &str,
    unresolved_cap: usize,
) -> ResolveOutcome {
    let mut out = ResolveOutcome::default();
    for c in calls {
        match c.verdict {
            Verdict::Negative => out.ancestral += 1,
            Verdict::NoCall => out.no_call += 1,
            Verdict::BackMutated => out.back_mutated += 1,
            Verdict::Positive => {
                let Some(resolved) = dict.resolve(&c.name, build) else {
                    out.unresolved += 1;
                    if out.unresolved_names.len() < unresolved_cap {
                        out.unresolved_names.push(c.name.clone());
                    }
                    continue;
                };
                let coord = resolved.coord;
                if !genotype_supports(&c.genotype, &coord.derived) {
                    out.strand_mismatches += 1;
                }
                if let Some(call) = variants::snp_call(
                    &coord.chrom,
                    coord.position,
                    &coord.ancestral,
                    &coord.derived,
                    None,
                    Some("1".into()),
                ) {
                    out.calls.push(call);
                }
            }
        }
    }
    out
}

/// Build the position→base map for **haplogroup placement** from BISDNA calls resolved on
/// `build`. Unlike [`resolve_calls`] (which emits only carried variants, for storage and the
/// allele-weighted reconciler), this includes **negatives** too: a negative is genuine
/// ancestral evidence that prunes over-deep branches in the Kulczynski scorer. Positive →
/// derived base, negative → ancestral base; `no_call`/back-mutated/dictionary-unresolved
/// markers are omitted (no confident base). Bases are uppercased; on duplicate positions the
/// last call wins. The result feeds `haplo::score` directly (`HashMap<position, base>`).
pub fn placement_calls(calls: &[BisdnaCall], dict: &YsnpDictionary, build: &str) -> HashMap<i64, char> {
    let mut map = HashMap::new();
    for c in calls {
        let want_derived = match c.verdict {
            Verdict::Positive => true,
            Verdict::Negative => false,
            Verdict::NoCall | Verdict::BackMutated => continue,
        };
        let Some(resolved) = dict.resolve(&c.name, build) else {
            continue;
        };
        let allele = if want_derived { &resolved.coord.derived } else { &resolved.coord.ancestral };
        if let Some(base) = allele.chars().next() {
            map.insert(resolved.coord.position, base.to_ascii_uppercase());
        }
    }
    map
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ysnp_dict::YsnpDictionary;

    /// The real file's shape: a long prose preamble, then the header, then rows.
    const SAMPLE: &str = "This is your Y chromosome raw data for the chromo2 chip. We only advise expert genetic genealogists to use this file. Alleles are called to the Illumina TOP strand. One such example is the marker S163.\n\
        SNPID\tgenotype\tresult\n\
        Apt\tGG\tnegative\n\
        CTS10003\tCC\tnegative\n\
        CTS10149\tGG\tpositive\n\
        CTS12633\tAT\tpositive\n\
        CTS3281\t00\tno_call\n\
        S163\tAA\t(positive)\n";

    #[test]
    fn parses_real_shape_skipping_preamble() {
        let calls = parse(SAMPLE).unwrap();
        assert_eq!(calls.len(), 6);
        assert_eq!(calls[0], BisdnaCall { name: "Apt".into(), genotype: "GG".into(), verdict: Verdict::Negative });
        assert_eq!(calls[2].verdict, Verdict::Positive);
    }

    #[test]
    fn het_positive_genotype_is_kept_verbatim() {
        let calls = parse(SAMPLE).unwrap();
        let cts = calls.iter().find(|c| c.name == "CTS12633").unwrap();
        assert_eq!(cts.genotype, "AT"); // apparent het, still a positive call
        assert_eq!(cts.verdict, Verdict::Positive);
    }

    #[test]
    fn no_call_and_parenthetical_positive() {
        let calls = parse(SAMPLE).unwrap();
        let nc = calls.iter().find(|c| c.name == "CTS3281").unwrap();
        assert_eq!(nc.genotype, "00");
        assert_eq!(nc.verdict, Verdict::NoCall);
        // `(positive)` (back-mutation-prone marker called derived) parses as Positive.
        assert_eq!(calls.iter().find(|c| c.name == "S163").unwrap().verdict, Verdict::Positive);
    }

    #[test]
    fn back_mutated_label_is_distinct() {
        let f = "SNPID\tgenotype\tresult\nS163\tAA\tback-mutated\n";
        assert_eq!(parse(f).unwrap()[0].verdict, Verdict::BackMutated);
    }

    #[test]
    fn quoted_and_padded_cells_are_cleaned() {
        let f = "SNPID\tgenotype\tresult\n\" CTS10003 \"\t\" CC \"\t\" negative \"\n";
        assert_eq!(parse(f).unwrap()[0], BisdnaCall { name: "CTS10003".into(), genotype: "CC".into(), verdict: Verdict::Negative });
    }

    #[test]
    fn unrecognized_verdict_row_is_skipped() {
        let f = "SNPID\tgenotype\tresult\nGood\tAA\tnegative\nBad\tAA\tmystery\n";
        let calls = parse(f).unwrap();
        assert_eq!(calls.len(), 1);
        assert_eq!(calls[0].name, "Good");
    }

    #[test]
    fn missing_header_errors() {
        assert!(parse("SNP,value\nCTS1,AA\n").is_err());
        assert!(parse("just prose, no table here\n").is_err());
    }

    #[test]
    fn header_without_data_errors() {
        assert!(parse("SNPID\tgenotype\tresult\n\n").is_err());
    }

    // ---- resolve_calls ------------------------------------------------------

    const DICT_TSV: &str = "\
name\tbuild\tchrom\tposition\tstrand\tancestral\tderived
CTS10149\ths1\tchrY\t14800000\t+\tA\tG
CTS12633\ths1\tchrY\t14900000\t-\tA\tT
S163\ths1\tchrY\t15000000\t+\tA\tC
";

    fn dict() -> YsnpDictionary {
        YsnpDictionary::from_text(DICT_TSV, "").unwrap()
    }

    #[test]
    fn resolves_only_positives_into_carried_calls() {
        let calls = parse(SAMPLE).unwrap(); // Apt-, CTS10003-, CTS10149+, CTS12633+, CTS3281 no_call, S163 (positive)
        let out = resolve_calls(&calls, &dict(), "hs1", 10);

        // Apt + CTS10003 are negative → counted, not emitted; CTS10003 also isn't in the dict.
        assert_eq!(out.ancestral, 2);
        assert_eq!(out.no_call, 1); // CTS3281
        // Three positives (CTS10149, CTS12633, S163) are all in the dict → three calls.
        assert_eq!(out.calls.len(), 3);
        assert_eq!(out.unresolved, 0);

        let c = out.calls.iter().find(|v| v.position == 14800000).unwrap();
        assert_eq!((c.reference.as_str(), c.alternate.as_str()), ("A", "G")); // ancestral→ref, derived→alt
        assert_eq!(c.genotype.as_deref(), Some("1")); // carried (derived)
        assert_eq!(c.contig, "chrY");
    }

    #[test]
    fn unresolved_positive_is_tallied_not_emitted() {
        // A positive whose name isn't in the dictionary.
        let f = "SNPID\tgenotype\tresult\nUNKNOWNSNP\tGG\tpositive\nCTS10149\tGG\tpositive\n";
        let out = resolve_calls(&parse(f).unwrap(), &dict(), "hs1", 10);
        assert_eq!(out.calls.len(), 1);
        assert_eq!(out.unresolved, 1);
        assert_eq!(out.unresolved_names, vec!["UNKNOWNSNP".to_string()]);
    }

    #[test]
    fn missing_build_makes_positives_unresolved() {
        // Dict only has hs1 coords; asking for GRCh38 resolves nothing.
        let out = resolve_calls(&parse(SAMPLE).unwrap(), &dict(), "GRCh38", 10);
        assert!(out.calls.is_empty());
        assert_eq!(out.unresolved, 3); // the three positives
    }

    #[test]
    fn genotype_strand_mismatch_is_flagged_but_still_emitted() {
        // CTS10149 derived=G; a positive genotype of "AA" carries neither G nor its complement C.
        let f = "SNPID\tgenotype\tresult\nCTS10149\tAA\tpositive\n";
        let out = resolve_calls(&parse(f).unwrap(), &dict(), "hs1", 10);
        assert_eq!(out.calls.len(), 1); // still emitted — we trust the verdict
        assert_eq!(out.strand_mismatches, 1);
    }

    #[test]
    fn placement_calls_include_ancestral_for_negatives() {
        // CTS10003 is negative & not in the dict → omitted. The three positives map to derived;
        // a synthetic negative on a dict marker maps to ancestral.
        let f = "SNPID\tgenotype\tresult\n\
                 CTS10149\tGG\tpositive\n\
                 CTS12633\tAT\tpositive\n\
                 S163\tAA\tnegative\n\
                 CTS3281\t00\tno_call\n";
        let calls = parse(f).unwrap();
        let map = placement_calls(&calls, &dict(), "hs1");
        assert_eq!(map.get(&14800000), Some(&'G')); // CTS10149 positive → derived G
        assert_eq!(map.get(&14900000), Some(&'T')); // CTS12633 positive → derived T
        assert_eq!(map.get(&15000000), Some(&'A')); // S163 negative → ancestral A
        assert_eq!(map.len(), 3); // no_call omitted
    }

    #[test]
    fn top_strand_complement_genotype_is_not_a_mismatch() {
        // CTS12633 derived=T; a "AA"… no. Use a marker derived=G with genotype "CC" (C=comp(G)).
        let f = "SNPID\tgenotype\tresult\nCTS10149\tCC\tpositive\n"; // derived G, genotype C = comp(G)
        let out = resolve_calls(&parse(f).unwrap(), &dict(), "hs1", 10);
        assert_eq!(out.strand_mismatches, 0);
    }
}
