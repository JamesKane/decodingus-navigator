//! Derive mtDNA variants by comparing a sample sequence to a reference (rCRS, the revised
//! Cambridge Reference Sequence, NC_012920.1, 16,569 bp). Vendor mtDNA FASTAs are already
//! aligned to rCRS, so this is a position-wise substitution diff over the common length —
//! the same approach as the Scala `MtDnaFastaProcessor.compareToRcrs`. Positions where
//! either base is `N` (or non-A/C/G/T) are skipped. Indels are not derived in v1.
//!
//! Pure: callers provide both sequences; the rCRS reference itself is supplied externally.

/// A single mtDNA substitution relative to the reference (1-based rCRS position).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MtVariant {
    pub position: i64,
    pub reference: char,
    pub alternate: char,
}

impl MtVariant {
    /// Compact mtDNA notation, e.g. `263A>G`.
    pub fn notation(&self) -> String {
        format!("{}{}>{}", self.position, self.reference, self.alternate)
    }
}

fn is_base(b: u8) -> bool {
    matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

/// Compare `sample` to `reference` position-by-position and return the substitutions. Both
/// are treated as rCRS-aligned (same coordinate origin); only the common-length prefix is
/// compared. Positions where either base is `N`/ambiguous are skipped (no confident call).
pub fn derive(reference: &str, sample: &str) -> Vec<MtVariant> {
    let r = reference.as_bytes();
    let s = sample.as_bytes();
    let n = r.len().min(s.len());
    let mut variants = Vec::new();
    for i in 0..n {
        let rb = r[i].to_ascii_uppercase();
        let sb = s[i].to_ascii_uppercase();
        if rb != sb && is_base(rb) && is_base(sb) {
            variants.push(MtVariant { position: (i + 1) as i64, reference: rb as char, alternate: sb as char });
        }
    }
    variants
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn finds_substitutions_and_skips_ns() {
        //          1234567
        let refseq = "ACGTACG";
        let sample = "AGGTNCG"; // pos2 C>G; pos5 ref A vs sample N -> skipped
        let v = derive(refseq, sample);
        assert_eq!(v.len(), 1);
        assert_eq!(v[0], MtVariant { position: 2, reference: 'C', alternate: 'G' });
        assert_eq!(v[0].notation(), "2C>G");
    }

    #[test]
    fn identical_sequences_have_no_variants() {
        assert!(derive("ACGTACGT", "acgtacgt").is_empty()); // case-insensitive
    }

    #[test]
    fn compares_only_the_common_length() {
        // trailing extra bases in the sample are ignored (indels not derived in v1)
        let v = derive("ACGT", "ACGTAAAA");
        assert!(v.is_empty());
    }
}
