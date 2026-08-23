//! Contig-name normalization: the one definition of "strip the `chr` prefix".
//!
//! The build determines the contig names. GRCh37 uses bare names (`22`, `X`, `MT`), and GRCh38 and
//! CHM13 use a `chr` prefix (`chr22`, `chrX`, `chrM`). Anything that matches loci across builds
//! must normalize first: panels, liftover, callsets, chip and vendor imports, and charts. This
//! lives in `navigator-domain` because every other crate depends on it, so there is exactly one
//! implementation, and not a closure at each call site.

/// `name` without a `chr` prefix at the start, in any case (`chr7` / `Chr7` / `CHR7` → `7`). A
/// name with no prefix (`7`, `MT`, `HLA-A`) comes back unchanged.
pub fn bare(name: &str) -> &str {
    match name.get(..3) {
        Some(p) if p.eq_ignore_ascii_case("chr") => &name[3..],
        _ => name,
    }
}

/// [`bare`] in upper case. This is the canonical key to match a contig across builds. The `chr1`
/// of a source then lines up with a panel locus that the store holds as `1`.
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
        // Shorter than the prefix, or it only starts with some of its letters.
        for name in ["", "1", "ch", "chX"] {
            assert_eq!(bare(name), name, "bare({name})");
        }
    }

    #[test]
    fn is_utf8_safe() {
        // `get(..3)` returns None on a non-char-boundary, and does not panic.
        assert_eq!(bare("é1"), "é1");
    }

    #[test]
    fn bare_upper_uppercases() {
        assert_eq!(bare_upper("chrx"), "X");
        assert_eq!(bare_upper("mt"), "MT");
    }
}
