//! Contig-name classification shared across walkers. Mirrors the Scala regexes:
//! autosomes `^(chr)?([1-9]|1[0-9]|2[0-2])$`, plus X / Y / M|MT.

fn core(name: &str) -> &str {
    name.strip_prefix("chr").unwrap_or(name)
}

/// Autosome 1-22 (no leading zeros).
pub fn is_autosome(name: &str) -> bool {
    let c = core(name);
    c.parse::<u32>()
        .map(|n| (1..=22).contains(&n) && c == n.to_string())
        .unwrap_or(false)
}

pub fn is_chr_x(name: &str) -> bool {
    core(name) == "X"
}

pub fn is_chr_y(name: &str) -> bool {
    core(name) == "Y"
}

pub fn is_chr_m(name: &str) -> bool {
    matches!(core(name), "M" | "MT")
}

/// Main assembly: autosomes + X/Y/M(T). Excludes alts, decoys, HLA, etc.
pub fn is_main_assembly(name: &str) -> bool {
    is_autosome(name) || is_chr_x(name) || is_chr_y(name) || is_chr_m(name)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classifies_contigs_like_the_scala_regexes() {
        for ok in ["1", "22", "X", "Y", "M", "MT", "chr1", "chr22", "chrX", "chrM"] {
            assert!(is_main_assembly(ok), "{ok} should be main assembly");
        }
        for no in ["0", "23", "01", "chr0", "chrUn", "chr1_KI270706v1_random", "HLA-A", "M1"] {
            assert!(!is_main_assembly(no), "{no} should not be main assembly");
        }
        assert!(is_autosome("chr21") && !is_autosome("chrX"));
        assert!(is_chr_x("X") && is_chr_y("chrY") && is_chr_m("MT"));
    }
}
