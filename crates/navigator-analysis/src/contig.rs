//! Contig-name classification shared across walkers. Mirrors the Scala regexes:
//! autosomes `^(chr)?([1-9]|1[0-9]|2[0-2])$`, plus X / Y / M|MT.
//!
//! Prefix stripping itself lives in [`navigator_domain::contig`] (every crate needs it, including
//! ones below this one) and is re-exported here so callers have a single import for contig work.
//! Classification is case-insensitive on both the prefix and the name (`chrx`, `Chr7`, `mt`).

pub use navigator_domain::contig::{bare, bare_upper};

/// Autosome 1-22 (no leading zeros).
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

/// **Haploid** contigs: chrY and chrM/MT carry a single allele, so the diploid (het `0/1` +
/// hom-alt `1/1`) model doesn't apply — the haploid caller and Y/mt haplogroup placement own them.
/// (chrX is haploid only in a male; that's left to the sex-aware refinement, not decided here.)
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
