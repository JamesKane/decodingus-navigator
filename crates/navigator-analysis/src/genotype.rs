//! Genotype-likelihood model for biallelic genotyping at a known site (the foundation
//! the population / ancestry / IBD paths need — they consume dosage 0/1/2 + a quality).
//!
//! Standard GATK/bcftools model. For a site with reference allele R and alternate A,
//! and per-read base observations with phred base qualities, the per-base error is
//! `e = 10^(-Q/10)` and `P(base | allele) = 1-e` on a match, `e/3` on a mismatch. For a
//! genotype carrying `g` alt copies out of `ploidy` P, the allele pool gives
//! `P(base | g) = [ (P-g)·P(base|R) + g·P(base|A) ] / P`. The genotype log-likelihoods
//! are summed over reads; the call is the argmax, with phred-scaled likelihoods (PL,
//! best = 0) and genotype quality `GQ` = the second-smallest PL.
//!
//! `ploidy` is supplied by the caller (sex → `sex::ploidy_for_contig`): 2 for autosomes
//! / female chrX, 1 for chrY / chrM / male chrX. Biallelic (ref + one alt) for v1.

/// Result of genotyping one site.
#[derive(Debug, Clone, PartialEq)]
pub struct GenotypeResult {
    /// Alt-allele count of the called genotype (0..=ploidy), or -1 for a no-call.
    pub dosage: i32,
    /// Phred-scaled likelihoods indexed by alt count (length `ploidy + 1`); best is 0.
    pub pls: Vec<u8>,
    /// Genotype quality (phred), capped at 99.
    pub gq: u8,
    /// Passing observations (ACGT bases clearing the quality filters).
    pub depth: u32,
    pub ref_depth: u32,
    pub alt_depth: u32,
}

const MAX_PL: f64 = 255.0;
const MAX_GQ: u8 = 99;

fn no_call(ploidy: u8, depth: u32, ref_depth: u32, alt_depth: u32) -> GenotypeResult {
    GenotypeResult { dosage: -1, pls: vec![0; ploidy as usize + 1], gq: 0, depth, ref_depth, alt_depth }
}

/// Call a biallelic genotype from passing `(base, phred_qual)` observations.
pub fn call_genotype(
    observations: &[(u8, u8)],
    reference_allele: u8,
    alternate_allele: u8,
    ploidy: u8,
    min_depth: u32,
) -> GenotypeResult {
    let r = reference_allele.to_ascii_uppercase();
    let a = alternate_allele.to_ascii_uppercase();
    let depth = observations.len() as u32;
    let ref_depth = observations.iter().filter(|(b, _)| b.to_ascii_uppercase() == r).count() as u32;
    let alt_depth = observations.iter().filter(|(b, _)| b.to_ascii_uppercase() == a).count() as u32;

    if ploidy == 0 || depth < min_depth {
        return no_call(ploidy.max(1), depth, ref_depth, alt_depth);
    }

    let p = ploidy as f64;
    let mut logl = vec![0.0f64; ploidy as usize + 1];
    for &(base, qual) in observations {
        let e = 10f64.powf(-(qual as f64) / 10.0);
        let b = base.to_ascii_uppercase();
        let p_ref = if b == r { 1.0 - e } else { e / 3.0 };
        let p_alt = if b == a { 1.0 - e } else { e / 3.0 };
        for (g, slot) in logl.iter_mut().enumerate() {
            let gf = g as f64;
            let p_bg = ((p - gf) * p_ref + gf * p_alt) / p;
            *slot += p_bg.max(1e-300).ln();
        }
    }

    // Phred-scale relative to the best genotype.
    let max_logl = logl.iter().cloned().fold(f64::MIN, f64::max);
    let ln10 = std::f64::consts::LN_10;
    let pls: Vec<u8> = logl
        .iter()
        .map(|&l| ((-10.0 * (l - max_logl) / ln10).round()).clamp(0.0, MAX_PL) as u8)
        .collect();

    let dosage = pls.iter().position(|&pl| pl == 0).unwrap_or(0) as i32;
    // GQ = second-smallest PL (confidence of the best call over the next-best).
    let mut sorted = pls.clone();
    sorted.sort_unstable();
    let gq = (*sorted.get(1).unwrap_or(&0)).min(MAX_GQ);

    GenotypeResult { dosage, pls, gq, depth, ref_depth, alt_depth }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn obs(base: u8, n: usize) -> Vec<(u8, u8)> {
        vec![(base, 40); n] // Phred 40
    }

    #[test]
    fn diploid_homozygous_reference() {
        let r = call_genotype(&obs(b'A', 20), b'A', b'G', 2, 4);
        assert_eq!(r.dosage, 0);
        assert_eq!(r.ref_depth, 20);
        assert_eq!(r.alt_depth, 0);
        assert!(r.gq > 50, "gq {}", r.gq);
        assert_eq!(r.pls[0], 0); // best is hom-ref
    }

    #[test]
    fn diploid_heterozygous() {
        let mut o = obs(b'C', 10); // ref
        o.extend(obs(b'G', 10)); // alt
        let r = call_genotype(&o, b'C', b'G', 2, 4);
        assert_eq!(r.dosage, 1);
        assert_eq!(r.ref_depth, 10);
        assert_eq!(r.alt_depth, 10);
        assert!(r.gq > 50);
        assert_eq!(r.pls[1], 0);
    }

    #[test]
    fn diploid_homozygous_alternate() {
        let r = call_genotype(&obs(b'A', 20), b'T', b'A', 2, 4);
        assert_eq!(r.dosage, 2);
        assert_eq!(r.alt_depth, 20);
        assert_eq!(r.pls[2], 0);
    }

    #[test]
    fn haploid_calls_zero_or_one() {
        assert_eq!(call_genotype(&obs(b'A', 10), b'A', b'G', 1, 4).dosage, 0);
        assert_eq!(call_genotype(&obs(b'G', 10), b'A', b'G', 1, 4).dosage, 1);
    }

    #[test]
    fn low_depth_is_a_no_call() {
        let r = call_genotype(&obs(b'A', 2), b'A', b'G', 2, 4);
        assert_eq!(r.dosage, -1);
        assert_eq!(r.gq, 0);
    }

    #[test]
    fn low_base_quality_lowers_confidence() {
        // a 10/10 split of low-quality bases is much less confident than high-quality.
        let mut hi = vec![(b'C', 40u8); 10];
        hi.extend(vec![(b'G', 40u8); 10]);
        let mut lo = vec![(b'C', 3u8); 10];
        lo.extend(vec![(b'G', 3u8); 10]);
        let g_hi = call_genotype(&hi, b'C', b'G', 2, 4);
        let g_lo = call_genotype(&lo, b'C', b'G', 2, 4);
        assert_eq!(g_hi.dosage, 1);
        assert!(g_lo.gq < g_hi.gq, "low-qual gq {} should be < high-qual gq {}", g_lo.gq, g_hi.gq);
    }
}
