//! The class of a contig name, which every walker shares. It follows the regular expressions of
//! the Scala code: an autosome matches `^(chr)?([1-9]|1[0-9]|2[0-2])$`, and X, Y and M or MT are
//! the others.
//!
//! The removal of the prefix lives in [`navigator_domain::contig`], because every crate needs it,
//! and that includes crates below this one. This module exports it again, so that a caller has one
//! import for all of its work on a contig.
//!
//! The class does not depend on the case, of the prefix or of the name. `chrx`, `Chr7` and `mt` all
//! work.

pub use navigator_domain::contig::{bare, bare_upper};

/// An autosome, from 1 to 22. The number carries no zero in front.
pub fn is_autosome(name: &str) -> bool {
    let c = bare(name);
    c.parse::<u32>()
        .map(|n| (1..=22).contains(&n) && c == n.to_string())
        .unwrap_or(false)
}

pub fn is_chr_x(name: &str) -> bool {
    bare(name).eq_ignore_ascii_case("X")
}

pub fn is_chr_y(name: &str) -> bool {
    bare(name).eq_ignore_ascii_case("Y")
}

pub fn is_chr_m(name: &str) -> bool {
    let c = bare(name);
    c.eq_ignore_ascii_case("M") || c.eq_ignore_ascii_case("MT")
}

/// Main assembly: autosomes + X/Y/M(T). Excludes alts, decoys, HLA, etc.
pub fn is_main_assembly(name: &str) -> bool {
    is_autosome(name) || is_chr_x(name) || is_chr_y(name) || is_chr_m(name)
}

/// The **haploid** contigs. chrY, and chrM or MT, each carry one allele. So the diploid model, with
/// its het `0/1` and hom-alt `1/1`, does not apply to them. The haploid caller owns them, and
/// so does the placement of a Y or mt haplogroup.
///
/// chrX is haploid in a male alone. The refinement that knows the sex decides that, and this
/// function does not.
pub fn is_haploid(name: &str) -> bool {
    is_chr_y(name) || is_chr_m(name)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classifies_contigs_like_the_scala_regexes() {
        for ok in [
            "1", "22", "X", "Y", "M", "MT", "chr1", "chr22", "chrX", "chrM", "chrx", "Chr7", "mt",
        ] {
            assert!(is_main_assembly(ok), "{ok} should be main assembly");
        }
        for no in [
            "0",
            "23",
            "01",
            "chr0",
            "chrUn",
            "chr1_KI270706v1_random",
            "HLA-A",
            "M1",
            "chrUn_KI270302v1",
            "GL000220.1",
            "",
        ] {
            assert!(!is_main_assembly(no), "{no} should not be main assembly");
        }
        assert!(is_autosome("chr21") && !is_autosome("chrX"));
        assert!(is_chr_x("X") && is_chr_y("chrY") && is_chr_m("MT"));
        assert!(is_haploid("chrY") && is_haploid("MT") && !is_haploid("chrX") && !is_haploid("1"));
    }
}
