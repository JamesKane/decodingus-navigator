//! Contig-name normalization — the one definition of "strip the `chr` prefix".
//!
//! Contig naming is build-determined: GRCh37 uses bare names (`22`, `X`, `MT`), GRCh38/CHM13 use a
//! `chr` prefix (`chr22`, `chrX`, `chrM`). Anything that matches loci across builds — panels,
//! liftover, callsets, chip/vendor imports, charts — has to normalize first. This lives in
//! `navigator-domain` because every other crate depends on it, so there is exactly one
//! implementation instead of a per-call-site closure.

/// `name` without a leading `chr` prefix, in any case (`chr7` / `Chr7` / `CHR7` → `7`). Names with
/// no prefix (`7`, `MT`, `HLA-A`) come back unchanged.
pub fn bare(name: &str) -> &str {
    match name.get(..3) {
        Some(p) if p.eq_ignore_ascii_case("chr") => &name[3..],
        _ => name,
    }
}

/// [`bare`] uppercased — the canonical key for matching a contig across builds, so a source's
/// `chr1` lines up with a panel locus stored as `1`.
pub fn bare_upper(name: &str) -> String {
    bare(name).to_ascii_uppercase()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strips_the_prefix_in_any_case() {
        for (input, want) in [
            ("chr7", "7"),
            ("Chr7", "7"),
            ("CHR7", "7"),
            ("chrX", "X"),
            ("chrM", "M"),
            ("7", "7"),
            ("MT", "MT"),
            ("HLA-A", "HLA-A"),
            ("chr1_KI270706v1_random", "1_KI270706v1_random"),
        ] {
            assert_eq!(bare(input), want, "bare({input})");
        }
    }

    #[test]
    fn leaves_short_and_lookalike_names_alone() {
        // Shorter than the prefix, or merely starting with some of its letters.
        for name in ["", "1", "ch", "chX"] {
            assert_eq!(bare(name), name, "bare({name})");
        }
    }

    #[test]
    fn is_utf8_safe() {
        // `get(..3)` returns None on a non-char-boundary rather than panicking.
        assert_eq!(bare("é1"), "é1");
    }

    #[test]
    fn bare_upper_uppercases() {
        assert_eq!(bare_upper("chrx"), "X");
        assert_eq!(bare_upper("mt"), "MT");
    }
}
