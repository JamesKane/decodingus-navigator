//! Vendor mtDNA FASTA sequences (Scala's `DataType.MtdnaFasta`) — a full mitochondrial
//! sequence (~16,569 bp, aligned to rCRS) imported from a `.fa`/`.fasta` export. Unlike a
//! chip, an mtDNA sequence is tiny, so we keep the sequence itself; calling variants vs
//! rCRS for haplogroup analysis is a later step. [`parse_fasta`] is a pure validator.

use du_domain::ids::SampleGuid;
use serde::{Deserialize, Serialize};

/// Plausible mtDNA length window (rCRS is 16,569 bp); guards against importing the wrong file.
const MIN_LEN: usize = 16_000;
const MAX_LEN: usize = 17_000;

/// A subject's imported mtDNA sequence.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct MtdnaSequence {
    pub id: i64,
    pub biosample_guid: SampleGuid,
    /// The FASTA header line (without the leading `>`), if any.
    pub defline: Option<String>,
    /// The full sequence, uppercased (A/C/G/T/N).
    pub sequence: String,
    /// Number of `N` (ambiguous) bases.
    pub n_count: i64,
    pub source_file_name: Option<String>,
}

impl MtdnaSequence {
    pub fn length(&self) -> usize {
        self.sequence.len()
    }
}

/// Fields for creating an mtDNA sequence (the store assigns the id).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NewMtdnaSequence {
    pub biosample_guid: SampleGuid,
    pub defline: Option<String>,
    pub sequence: String,
    pub n_count: i64,
    pub source_file_name: Option<String>,
}

/// A validated mtDNA FASTA: its header (sans `>`) and concatenated uppercase sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParsedMtdna {
    pub defline: Option<String>,
    pub sequence: String,
    pub n_count: i64,
}

/// Parse and validate a single-record mtDNA FASTA: must start with a `>` header; the
/// concatenated sequence must be ~16,569 bp (16,000–17,000) and contain only A/C/G/T/N.
/// Only the first record is read. Returns the sequence + `N` count.
pub fn parse_fasta(text: &str) -> Result<ParsedMtdna, String> {
    let mut lines = text.lines().map(str::trim).filter(|l| !l.is_empty());

    let Some(first) = lines.next() else {
        return Err("empty file".into());
    };
    if !first.starts_with('>') {
        return Err("file does not start with a FASTA header (>)".into());
    }
    let defline = {
        let d = first[1..].trim();
        (!d.is_empty()).then(|| d.to_string())
    };

    let mut sequence = String::new();
    for line in lines {
        if line.starts_with('>') {
            break; // stop at the next record — single-sequence import
        }
        sequence.push_str(&line.to_ascii_uppercase());
    }

    if sequence.len() < MIN_LEN {
        return Err(format!(
            "sequence is too short ({} bp); expected ~16,569 bp for mtDNA",
            sequence.len()
        ));
    }
    if sequence.len() > MAX_LEN {
        return Err(format!(
            "sequence is too long ({} bp); expected ~16,569 bp for mtDNA",
            sequence.len()
        ));
    }
    let mut n_count = 0i64;
    for b in sequence.bytes() {
        match b {
            b'A' | b'C' | b'G' | b'T' => {}
            b'N' => n_count += 1,
            other => return Err(format!("sequence contains an invalid base: {:?}", other as char)),
        }
    }
    Ok(ParsedMtdna {
        defline,
        sequence,
        n_count,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fasta(seq_body: &str) -> String {
        format!(">rCRS test\n{seq_body}\n")
    }

    #[test]
    fn parses_and_counts_ns() {
        // 16,569 bp: mostly A, with 3 Ns and a wrapped line.
        let mut body = "A".repeat(16_566);
        body.insert_str(8000, "NNN");
        // wrap into 70-col lines to exercise multi-line concatenation
        let wrapped: String = body
            .as_bytes()
            .chunks(70)
            .map(|c| format!("{}\n", std::str::from_utf8(c).unwrap()))
            .collect();
        let p = parse_fasta(&fasta(&wrapped)).unwrap();
        assert_eq!(p.defline.as_deref(), Some("rCRS test"));
        assert_eq!(p.sequence.len(), 16_569);
        assert_eq!(p.n_count, 3);
    }

    #[test]
    fn rejects_missing_header() {
        assert!(parse_fasta(&"A".repeat(16_569)).is_err());
    }

    #[test]
    fn rejects_wrong_length() {
        assert!(parse_fasta(&fasta(&"ACGT".repeat(100))).is_err()); // far too short
    }

    #[test]
    fn rejects_invalid_base() {
        let body = format!("{}Z", "A".repeat(16_569));
        assert!(parse_fasta(&fasta(&body)).is_err());
    }
}
