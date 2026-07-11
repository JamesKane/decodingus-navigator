//! File-type detection for the unified "Add data" flow (Scala's `FileTypeDetector`).
//! Binary/structured formats are detected by extension; ambiguous text tables (STR vs
//! chip) are scored by content fingerprint. Pure: callers pass the name + a head sample.

/// What a dropped/picked file looks like.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DetectedData {
    /// BAM/CRAM aligned reads (attaches to a sequencing test, not directly to a subject).
    Alignment,
    /// VCF variant calls.
    Variants,
    /// CompleteGenomics `masterVar` whole-genome variant table (`var-*-ASM.tsv[.bz2]`).
    CompleteGenomicsVar,
    /// FTDNA Big Y CSV variant report (Named/Private Variants) — chrY derived calls, the
    /// "lesser access" substitute for the BAM/CRAM/VCF.
    FtdnaCsvVariants,
    /// Y-STR profile table.
    StrProfile,
    /// Named Y-SNP panel (e.g. BISDNA chromo2) — name + genotype + positive/negative verdict.
    YSnpPanel,
    /// Genotyping-array (chip) export.
    ChipData,
    /// mtDNA FASTA sequence.
    MtdnaFasta,
    /// Unrecognized.
    Unknown,
}

impl DetectedData {
    pub fn description(self) -> &'static str {
        match self {
            DetectedData::Alignment => "Alignment file",
            DetectedData::Variants => "VCF variants",
            DetectedData::CompleteGenomicsVar => "CompleteGenomics masterVar",
            DetectedData::FtdnaCsvVariants => "FTDNA Big Y variant CSV",
            DetectedData::StrProfile => "Y-STR profile",
            DetectedData::YSnpPanel => "Y-SNP panel",
            DetectedData::ChipData => "Chip / array data",
            DetectedData::MtdnaFasta => "mtDNA FASTA",
            DetectedData::Unknown => "Unknown format",
        }
    }
}

/// Detect a file's type from its `file_name` and a `head` sample of its text content
/// (ignored for binary formats). Mirrors the Scala detector: extension first, then a
/// STR-vs-chip content score with the same thresholds.
pub fn detect(file_name: &str, head: &str) -> DetectedData {
    let name = file_name.to_ascii_lowercase();
    let ends = |s: &str| name.ends_with(s);

    if ends(".bam") || ends(".cram") {
        return DetectedData::Alignment;
    }
    if ends(".vcf") || ends(".vcf.gz") || ends(".vcf.bgz") {
        return DetectedData::Variants;
    }
    if ends(".fasta")
        || ends(".fa")
        || ends(".fna")
        || ends(".fas")
        || ends(".fasta.gz")
        || ends(".fa.gz")
        || ends(".fna.gz")
    {
        return DetectedData::MtdnaFasta;
    }

    // CompleteGenomics masterVar — a whole-genome variant TSV (`.tsv[.bz2]`) with an unambiguous
    // `>locus ploidy allele chromosome …` column header and a `cgatools`/`VAR-ANNOTATION` preamble.
    // Checked here (before the STR/chip scorer) on the head, which the caller has decompressed.
    if looks_like_cg_master_var(head) {
        return DetectedData::CompleteGenomicsVar;
    }

    // Text content: score STR vs chip.
    let lines: Vec<&str> = head.lines().take(50).collect();
    let data_lines: Vec<&str> = lines
        .iter()
        .filter(|l| !l.trim_start().starts_with('#') && !l.trim().is_empty())
        .copied()
        .collect();
    if data_lines.is_empty() {
        return DetectedData::Unknown;
    }

    // FTDNA Big Y Named/Private Variants CSV — an exact header signature, checked before the
    // STR/chip scorer (which would otherwise mis-score the named report as chip).
    if crate::ftdna_csv::looks_like_ftdna_variant_csv(head) {
        return DetectedData::FtdnaCsvVariants;
    }

    // A named Y-SNP panel (BISDNA chromo2) is unambiguous — check it before the STR/chip
    // scorer, which would otherwise mis-score it as chip.
    if looks_like_ysnp_panel(&lines) {
        return DetectedData::YSnpPanel;
    }

    let str_score = str_score(&data_lines, &name);
    let chip_score = chip_score(&lines, &data_lines, &name);
    if str_score > chip_score && str_score >= 3 {
        DetectedData::StrProfile
    } else if chip_score > str_score && chip_score >= 3 {
        DetectedData::ChipData
    } else if str_score >= 2 {
        DetectedData::StrProfile
    } else if chip_score >= 2 {
        DetectedData::ChipData
    } else {
        DetectedData::Unknown
    }
}

/// Recognize a CompleteGenomics masterVar table from its head text. The `>locus … chromosome …
/// varType …` column header is the unambiguous signature; the `cgatools` / `VAR-ANNOTATION`
/// comment preamble corroborates it. Tolerant of the file being uncompressed here (the caller
/// decompresses `.bz2` / `.gz` before sniffing).
fn looks_like_cg_master_var(head: &str) -> bool {
    let mut has_column_header = false;
    let mut has_preamble = false;
    for line in head.lines().take(64) {
        if let Some(cols) = line.strip_prefix('>') {
            let l = cols.to_ascii_lowercase();
            if l.starts_with("locus\t") && l.contains("chromosome") && l.contains("vartype") && l.contains("alleleseq")
            {
                has_column_header = true;
            }
        } else if line.starts_with('#') {
            let l = line.to_ascii_lowercase();
            if l.contains("cgatools") || l.contains("var-annotation") {
                has_preamble = true;
            }
        }
    }
    // The column header alone is definitive; the preamble alone (without it) is not enough.
    has_column_header || (has_preamble && head.contains("masterVar"))
}

fn split_line(line: &str) -> Vec<&str> {
    let sep = if line.contains('\t') { '\t' } else { ',' };
    line.split(sep).map(|s| s.trim().trim_matches('"')).collect()
}

/// Count occurrences of `prefix` immediately followed by `min_digits`–`max_digits` digits.
fn count_token(haystack: &str, prefix: &str, min_digits: usize, max_digits: usize) -> usize {
    let bytes = haystack.as_bytes();
    let p = prefix.as_bytes();
    let mut count = 0;
    let mut i = 0;
    while i + p.len() <= bytes.len() {
        if bytes[i..].starts_with(p) {
            let mut d = 0;
            while i + p.len() + d < bytes.len() && bytes[i + p.len() + d].is_ascii_digit() {
                d += 1;
            }
            if d >= min_digits && d <= max_digits {
                count += 1;
            }
            i += p.len() + d.max(1);
        } else {
            i += 1;
        }
    }
    count
}

/// Recognize a named Y-SNP panel (BISDNA chromo2): either the exact
/// `SNPID<TAB>genotype<TAB>result` header, or — lacking it — several tab rows whose third
/// column is a positive/negative/no_call/back-mutated verdict. Tolerant of the multi-line
/// prose preamble BISDNA prepends (those lines aren't tab-delimited and never match).
fn looks_like_ysnp_panel(lines: &[&str]) -> bool {
    let is_verdict = |s: &str| {
        let v = s.trim().trim_matches(|c| c == '"').to_ascii_lowercase();
        let core = v.trim_matches(|c| c == '(' || c == ')');
        matches!(core, "positive" | "negative" | "no_call" | "back-mutated")
    };

    let mut verdict_rows = 0;
    for line in lines {
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 3 {
            continue;
        }
        let (a, b, c) = (cols[0].trim(), cols[1].trim(), cols[2].trim());
        if a.eq_ignore_ascii_case("snpid") && b.eq_ignore_ascii_case("genotype") && c.eq_ignore_ascii_case("result") {
            return true; // exact header — unambiguous
        }
        if is_verdict(c) {
            verdict_rows += 1;
        }
    }
    verdict_rows >= 3
}

fn str_score(data_lines: &[&str], file_name: &str) -> i32 {
    let mut score = 0;
    if file_name.contains("str") || file_name.contains("ystr") {
        score += 2;
    }
    if file_name.contains("ftdna") || file_name.contains("yseq") {
        score += 1;
    }

    let content = data_lines.join("\n").to_ascii_uppercase();
    let dys = count_token(&content, "DYS", 2, 3);
    score += match dys {
        d if d >= 10 => 4,
        d if d >= 5 => 3,
        d if d >= 2 => 2,
        d if d >= 1 => 1,
        _ => 0,
    };
    if ["DYF", "GATA", "YCAII", "CDY"].iter().any(|m| content.contains(m)) {
        score += 2;
    }

    // Two-column `marker,value` with a small-integer (or "a-b") value is STR-shaped.
    if let Some(first) = data_lines.first() {
        let cols = split_line(first);
        if cols.len() == 2 {
            let v = cols[1];
            let small_int = !v.is_empty() && v.len() <= 2 && v.bytes().all(|b| b.is_ascii_digit());
            let range = v.split_once('-').is_some_and(|(a, b)| {
                !a.is_empty() && !b.is_empty() && a.bytes().chain(b.bytes()).all(|x| x.is_ascii_digit())
            });
            if small_int || range {
                score += 1;
            }
        }
    }
    score
}

fn chip_score(lines: &[&str], data_lines: &[&str], file_name: &str) -> i32 {
    let mut score = 0;
    for (token, pts) in [("23andme", 3), ("ancestry", 3), ("myheritage", 3), ("livingdna", 3)] {
        if file_name.contains(token) {
            score += pts;
        }
    }
    if file_name.contains("ftdna") && file_name.contains("raw") {
        score += 2;
    }
    if file_name.contains("snp") || file_name.contains("chip") || file_name.contains("array") {
        score += 1;
    }

    let comments: String = lines
        .iter()
        .filter(|l| l.trim_start().starts_with('#') || l.trim_start().starts_with("\"#"))
        .map(|l| l.to_ascii_lowercase())
        .collect::<Vec<_>>()
        .join("\n");
    for token in [
        "23andme",
        "ancestrydna",
        "ancestry dna",
        "myheritage",
        "living dna",
        "livingdna",
    ] {
        if comments.contains(token) {
            score += 3;
            break;
        }
    }

    let content_lower = data_lines
        .iter()
        .take(20)
        .copied()
        .collect::<Vec<_>>()
        .join("\n")
        .to_ascii_lowercase();
    let rsid = count_token(&content_lower, "rs", 4, 12);
    score += match rsid {
        d if d >= 10 => 4,
        d if d >= 5 => 3,
        d if d >= 2 => 2,
        d if d >= 1 => 1,
        _ => 0,
    };

    let header = data_lines.first().map(|l| l.to_ascii_uppercase()).unwrap_or_default();
    let indicators = [
        "CHROMOSOME",
        "CHROM",
        "CHR",
        "POSITION",
        "POS",
        "GENOTYPE",
        "ALLELE1",
        "ALLELE2",
        "RSID",
    ];
    let header_hits = indicators.iter().filter(|h| header.contains(*h)).count();
    score += match header_hits {
        h if h >= 3 => 3,
        2 => 2,
        1 => 1,
        _ => 0,
    };

    // A 4–6 column row whose last field is a genotype-ish token.
    if let Some(first) = data_lines.first() {
        let cols = split_line(first);
        if (4..=6).contains(&cols.len()) {
            let last = cols.last().unwrap().to_ascii_uppercase();
            let is_gt = (last.len() == 2 && last.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T')))
                || last == "--"
                || last == "NC";
            if is_gt {
                score += 2;
            }
        }
    }
    score
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extensions_win_first() {
        assert_eq!(detect("HG002.bam", ""), DetectedData::Alignment);
        assert_eq!(detect("x.cram", ""), DetectedData::Alignment);
        assert_eq!(detect("calls.vcf", ""), DetectedData::Variants);
        assert_eq!(detect("calls.vcf.gz", ""), DetectedData::Variants);
        assert_eq!(detect("calls.vcf.bgz", ""), DetectedData::Variants);
        assert_eq!(detect("seq.fasta", ""), DetectedData::MtdnaFasta);
        assert_eq!(detect("seq.fa", ""), DetectedData::MtdnaFasta);
    }

    #[test]
    fn detects_23andme_chip_by_content() {
        let head = "# This data file generated by 23andMe\nrsid\tchromosome\tposition\tgenotype\n\
                    rs4477212\t1\t82154\tAA\nrs3094315\t1\t752566\tAG\n";
        assert_eq!(detect("genome_data.txt", head), DetectedData::ChipData);
    }

    #[test]
    fn detects_completegenomics_master_var_by_content() {
        let head = "#ASSEMBLY_ID\tGS00253-DNA_A01_200_37-ASM\n\
                    #GENOME_REFERENCE\tNCBI build 37\n\
                    #GENERATED_BY\tcgatools\n\
                    #TYPE\tVAR-ANNOTATION\n\
                    >locus\tploidy\tallele\tchromosome\tbegin\tend\tvarType\treference\talleleSeq\tvarScoreVAF\tvarScoreEAF\tvarQuality\thapLink\txRef\n\
                    1\t2\tall\tchr1\t0\t10000\tno-ref\t=\t?\t\t\t\t\t\n";
        // Both the raw name and a `.tsv.bz2` (extension isn't consulted for this format) detect.
        assert_eq!(detect("var-GS00253-DNA_A01_200_37-ASM.tsv", head), DetectedData::CompleteGenomicsVar);
        assert_eq!(detect("var-GS00253-DNA_A01_200_37-ASM.tsv.bz2", head), DetectedData::CompleteGenomicsVar);
    }

    #[test]
    fn detects_str_profile_by_content() {
        let head = "Marker,Value\nDYS393,13\nDYS390,24\nDYS19,14\nDYS391,11\nDYS385,11-14\nDYS426,12\n";
        assert_eq!(detect("markers.csv", head), DetectedData::StrProfile);
    }

    #[test]
    fn detects_bisdna_by_header_after_preamble() {
        let head = "This is your Y chromosome raw data for the chromo2 chip. Expert use only.\n\
                    SNPID\tgenotype\tresult\nApt\tGG\tnegative\nCTS10149\tGG\tpositive\nCTS3281\t00\tno_call\n";
        assert_eq!(detect("results.txt", head), DetectedData::YSnpPanel);
    }

    #[test]
    fn detects_bisdna_by_verdict_column_without_header() {
        let head = "M269\tCC\tpositive\nCTS10003\tGG\tnegative\nL21\tAA\t(positive)\nDF27\tTT\tnegative\n";
        assert_eq!(detect("snps.txt", head), DetectedData::YSnpPanel);
    }

    #[test]
    fn ysnp_panel_not_confused_with_chip_or_str() {
        // A real 23andMe chip and an STR table must NOT read as a Y-SNP panel.
        let chip = "# 23andMe\nrsid\tchromosome\tposition\tgenotype\nrs1\t1\t100\tAA\nrs2\t1\t200\tAG\n";
        assert_eq!(detect("genome.txt", chip), DetectedData::ChipData);
        let str_tbl = "DYS393,13\nDYS390,24\nDYS19,14\nDYS391,11\n";
        assert_eq!(detect("markers.csv", str_tbl), DetectedData::StrProfile);
    }

    #[test]
    fn junk_is_unknown() {
        assert_eq!(
            detect("notes.txt", "hello world\nthis is not genetic data\n"),
            DetectedData::Unknown
        );
    }
}
