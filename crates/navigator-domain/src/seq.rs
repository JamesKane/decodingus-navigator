//! Tiny shared sequence helpers used across the desktop crates (navigator-domain, -analysis, -app
//! all depend on navigator-domain, so this is their common home — no extra dependencies). Keep it to
//! pure, allocation-free base math.

/// Watson–Crick complement of a single base (`char`); non-ACGT (incl. lowercase non-matches and `N`)
/// passes through unchanged. Used for strand reconciliation against a reference/tree.
pub fn complement_base(b: char) -> char {
    match b.to_ascii_uppercase() {
        'A' => 'T',
        'T' => 'A',
        'C' => 'G',
        'G' => 'C',
        other => other,
    }
}

/// Byte (`u8`) variant of [`complement_base`] for callers working on raw allele bytes.
pub fn complement_base_u8(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        other => other,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complements_acgt_and_passes_through_others() {
        assert_eq!(complement_base('A'), 'T');
        assert_eq!(complement_base('c'), 'G'); // case-insensitive, uppercased result
        assert_eq!(complement_base('N'), 'N');
        assert_eq!(complement_base('-'), '-');
        assert_eq!(complement_base_u8(b'G'), b'C');
        assert_eq!(complement_base_u8(b'N'), b'N');
    }
}
