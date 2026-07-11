//! Multi-build, chip-compatible IBD reference panel (ancestry-ibd-asset-wiring B2/B2c).
//!
//! IBD matching needs a neutral, dense SNP set that's also **assayed by consumer arrays** — chip
//! kits outnumber WGS by orders of magnitude, so the panel must be where chip and WGS overlap.
//! Each site carries its `(contig, pos, REF, ALT)` on **CHM13, GRCh37, and GRCh38** (built once via
//! allele-aware GATK liftover, offline), so a chip genotype on *any* build resolves to the canonical
//! CHM13 site + orientation with **no runtime liftover** — the panel pre-computes it.
//!
//! Two correctness rules:
//! - The per-build loci carry the **same biological alleles** (GATK reverse-complements / swaps
//!   REF↔ALT on inverted chain blocks), so "count of that build's ALT" == "count of the CHM13 ALT".
//!   The dosage is therefore build-agnostic.
//! - **Strand-ambiguous palindromes (A/T, C/G) are excluded** ([`is_palindromic`]) — `rc(A)=T` is
//!   also a valid allele, so a chip's strand can't be disambiguated by allele comparison.

use std::collections::HashMap;

use serde::{Deserialize, Serialize};

use crate::ancestry::dosage_from_alleles;
use crate::caller::SiteGenotype;
use crate::error::AnalysisError;

/// A site's locus on one reference build: coordinates + the `(REF, ALT)` on that build's strand.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Locus {
    pub contig: String,
    pub position: i64,
    pub reference: char,
    pub alternate: char,
}

/// One IBD panel site (a chip-assayed biallelic SNP). The CHM13 locus is canonical; GRCh37/38 are
/// present when the site lifts cleanly to those builds.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct IbdPanelSite {
    pub rsid: String,
    pub chm13: Locus,
    #[serde(default)]
    pub grch37: Option<Locus>,
    #[serde(default)]
    pub grch38: Option<Locus>,
}

impl IbdPanelSite {
    /// The locus for a build name (`GRCh37`/`hg19`/`b37`, `GRCh38`/`hg38`, `chm13`/`hs1`/`t2t`).
    pub fn locus(&self, build: &str) -> Option<&Locus> {
        let b = build.to_ascii_lowercase();
        if b.contains("38") || b == "hg38" {
            self.grch38.as_ref()
        } else if b.contains("37") || b == "hg19" || b == "b37" {
            self.grch37.as_ref()
        } else if b.contains("chm13") || b == "hs1" || b == "t2t" {
            Some(&self.chm13)
        } else {
            None
        }
    }
}

/// A multi-build IBD reference panel. `build` is the canonical build of the `chm13` loci.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct IbdPanel {
    pub build: String,
    pub sites: Vec<IbdPanelSite>,
}

impl IbdPanel {
    /// Deserialize the built asset (bincode).
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, AnalysisError> {
        bincode::deserialize(bytes).map_err(|e| AnalysisError::Message(format!("ibd panel decode: {e}")))
    }

    /// Serialize to the binary asset form (bincode).
    pub fn to_bytes(&self) -> Result<Vec<u8>, AnalysisError> {
        bincode::serialize(self).map_err(|e| AnalysisError::Message(format!("ibd panel encode: {e}")))
    }

    /// Build from sites, **retaining** strand-ambiguous palindromes (A/T, C/G). The panel is a probe
    /// superset: WGS + ancestry genotype palindromic sites fine (a read gives the reference-strand
    /// base), and only the CHIP path can't orient them — so [`resolve_chip`] skips palindromes at
    /// resolve time rather than excluding them from the panel here. Returns `(panel, n_palindromic)`
    /// (the retained palindrome count, for the build log).
    pub fn from_sites(build: impl Into<String>, sites: Vec<IbdPanelSite>) -> (Self, usize) {
        let palindromic = sites
            .iter()
            .filter(|s| is_palindromic(s.chm13.reference, s.chm13.alternate))
            .count();
        (
            IbdPanel {
                build: build.into(),
                sites,
            },
            palindromic,
        )
    }

    /// Resolve chip calls (on `build`, as `(contig, pos, a1, a2)`) to canonical CHM13 dosages.
    /// Indexes the panel by the build's `(contig, position)`, then counts copies of the **canonical
    /// CHM13 ALT** directly from each observed pair — direct or reverse-complement
    /// ([`dosage_from_alleles`]) — and emits it as a [`SiteGenotype`] at the CHM13 locus. No runtime
    /// liftover; no alignment. Unmatched / no-call / non-reconciling calls are dropped.
    ///
    /// We deliberately score against the **CHM13** `(REF, ALT)`, not the build locus's: chip allele
    /// letters are absolute, and the asset's per-build `(REF, ALT)` labels are *not* reliably oriented
    /// to the CHM13 ALT (a large fraction are ref/alt-swapped relative to CHM13 — the GRCh37 reference
    /// allele is often the CHM13 ALT), so scoring against the build ALT flips the dosage 0↔2 at those
    /// sites. Comparing the chip alleles to the CHM13 alleles (with the rc retry for strand) is
    /// orientation-bug-proof; the build locus is used only to look the site up by position.
    pub fn resolve_chip(&self, build: &str, calls: &[(String, i64, char, char)]) -> Vec<SiteGenotype> {
        let mut index: HashMap<(&str, i64), &IbdPanelSite> = HashMap::new();
        for s in &self.sites {
            if let Some(l) = s.locus(build) {
                index.insert((l.contig.as_str(), l.position), s);
            }
        }
        let mut out = Vec::new();
        for (contig, pos, a1, a2) in calls {
            let Some(site) = index.get(&(contig.as_str(), *pos)) else {
                continue;
            };
            // Strand-ambiguous palindromes (A/T, C/G) can't be oriented from a chip's reported
            // alleles — skip them for the chip path (WGS/ancestry still use them via direct base
            // calls). The probe panel retains them; this is where the chip-only exclusion lives.
            if is_palindromic(site.chm13.reference, site.chm13.alternate) {
                continue;
            }
            let Some(dosage) = dosage_from_alleles(*a1, *a2, site.chm13.reference, site.chm13.alternate) else {
                continue;
            };
            out.push(panel_site_genotype(site, dosage));
        }
        out
    }

    /// Resolve a **whole-genome, variant-only** source (a WGS VCF or CompleteGenomics masterVar) to
    /// canonical CHM13 dosages over the panel. Unlike a chip — which reports a genotype at every
    /// array site — such a source lists *only* the non-reference sites, so every panel site the
    /// source could have called but didn't is taken as **homozygous reference** (dosage 0). That
    /// assumption is valid **only** for a source that genotyped the whole genome (absent ⇒ hom-ref,
    /// not no-call); never pass a targeted panel (Big Y / Sanger) here.
    ///
    /// `variant_calls` are the source's variant sites on `build` as `(contig, pos, a1, a2)`
    /// reference-forward allele pairs. Contigs match `chr`-insensitively (a source's `chr1` lines up
    /// with a panel `grch37` locus stored as `1`). A variant whose alleles don't reconcile to the
    /// site (multiallelic mismatch) is dropped, not mis-called hom-ref. Palindromic (A/T, C/G) sites
    /// are skipped — strand-ambiguous across builds, exactly as [`resolve_chip`].
    pub fn resolve_whole_genome(&self, build: &str, variant_calls: &[(String, i64, char, char)]) -> Vec<SiteGenotype> {
        let norm = |c: &str| c.strip_prefix("chr").unwrap_or(c).to_ascii_uppercase();
        let variants: HashMap<(String, i64), (char, char)> = variant_calls
            .iter()
            .map(|(c, p, a1, a2)| ((norm(c), *p), (*a1, *a2)))
            .collect();

        let mut out = Vec::new();
        for site in &self.sites {
            let Some(l) = site.locus(build) else {
                continue; // no coordinate on this build → the source can't have called it
            };
            if is_palindromic(site.chm13.reference, site.chm13.alternate) {
                continue;
            }
            let dosage = match variants.get(&(norm(&l.contig), l.position)) {
                Some((a1, a2)) => match dosage_from_alleles(*a1, *a2, site.chm13.reference, site.chm13.alternate) {
                    Some(d) => d,
                    None => continue, // observed a variant but it didn't reconcile — drop it
                },
                None => 0, // whole-genome source didn't list this site ⇒ homozygous reference
            };
            out.push(panel_site_genotype(site, dosage));
        }
        out
    }

    /// The canonical CHM13 `(contig, position)` sites — the targets a WGS caller genotypes so its
    /// dosages line up with the chip path.
    pub fn chm13_sites(&self) -> Vec<(&str, i64)> {
        self.sites
            .iter()
            .map(|s| (s.chm13.contig.as_str(), s.chm13.position))
            .collect()
    }
}

/// Build a diploid [`SiteGenotype`] at a panel site's canonical CHM13 locus with the given alt
/// dosage — the shared emit for the chip and whole-genome resolvers (no depth/quality; dosage only).
fn panel_site_genotype(site: &IbdPanelSite, dosage: i32) -> SiteGenotype {
    SiteGenotype {
        name: site.rsid.clone(),
        contig: site.chm13.contig.clone(),
        position: site.chm13.position,
        reference_allele: site.chm13.reference.to_string(),
        alternate_allele: site.chm13.alternate.to_string(),
        ploidy: 2,
        dosage,
        gq: 0,
        depth: 0,
        ref_depth: 0,
        alt_depth: 0,
        pls: Vec::new(),
        gt: None,
        allele_depths: None,
    }
}

/// Whether `(a, b)` is a strand-ambiguous palindrome (A/T or C/G) — excluded from a chip-compatible
/// panel because reverse-complement can't disambiguate the array's strand.
pub fn is_palindromic(a: char, b: char) -> bool {
    matches!(
        (a.to_ascii_uppercase(), b.to_ascii_uppercase()),
        ('A', 'T') | ('T', 'A') | ('C', 'G') | ('G', 'C')
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    fn site(rsid: &str, chm13: (i64, char, char), g37: Option<(i64, char, char)>) -> IbdPanelSite {
        IbdPanelSite {
            rsid: rsid.into(),
            chm13: Locus {
                contig: "chr1".into(),
                position: chm13.0,
                reference: chm13.1,
                alternate: chm13.2,
            },
            grch37: g37.map(|(p, r, a)| Locus {
                contig: "1".into(),
                position: p,
                reference: r,
                alternate: a,
            }),
            grch38: None,
        }
    }

    #[test]
    fn palindromes_retained_in_panel_skipped_for_chip() {
        assert!(is_palindromic('A', 'T') && is_palindromic('C', 'G') && is_palindromic('g', 'c'));
        assert!(!is_palindromic('A', 'G') && !is_palindromic('C', 'T'));
        let sites = vec![
            site("rs1", (100, 'A', 'G'), Some((500, 'A', 'G'))), // non-palindromic
            site("rs2", (200, 'A', 'T'), Some((600, 'A', 'T'))), // palindrome
            site("rs3", (300, 'C', 'G'), Some((700, 'C', 'G'))), // palindrome
        ];
        // The probe panel RETAINS palindromes (count reported); WGS/ancestry use them.
        let (panel, palindromic) = IbdPanel::from_sites("chm13v2.0", sites);
        assert_eq!(palindromic, 2);
        assert_eq!(panel.sites.len(), 3);
        // The chip path skips palindromes (can't orient strand) but resolves the non-palindromic one.
        let g = panel.resolve_chip(
            "GRCh37",
            &[
                ("1".into(), 500, 'A', 'G'),
                ("1".into(), 600, 'A', 'T'),
                ("1".into(), 700, 'C', 'G'),
            ],
        );
        assert_eq!(g.len(), 1);
        assert_eq!(g[0].name, "rs1");
    }

    #[test]
    fn resolve_chip_same_and_opposite_strand() {
        // rs1: GRCh37 1:500 A/G (same alleles as CHM13 chr1:100 A/G).
        // rs2: GRCh37 1:600 T/C — CHM13 chr1:200 A/G (a strand flip: GRCh37 alleles are rc).
        let (panel, _) = IbdPanel::from_sites(
            "chm13v2.0",
            vec![
                site("rs1", (100, 'A', 'G'), Some((500, 'A', 'G'))),
                site("rs2", (200, 'A', 'G'), Some((600, 'T', 'C'))),
            ],
        );
        // Chip on GRCh37: rs1 het AG → dosage 1; rs2 het TC → reconciles via rc → dosage 1.
        let calls = vec![
            ("1".to_string(), 500, 'A', 'G'),
            ("1".to_string(), 600, 'T', 'C'),
            ("1".to_string(), 999, 'A', 'G'), // no panel site → dropped
        ];
        let g = panel.resolve_chip("GRCh37", &calls);
        assert_eq!(g.len(), 2);
        // Output is at the canonical CHM13 loci with build-agnostic ALT dosage.
        let by_pos: std::collections::HashMap<i64, i32> = g.iter().map(|s| (s.position, s.dosage)).collect();
        assert_eq!(by_pos.get(&100), Some(&1)); // AG → one ALT(G)
        assert_eq!(by_pos.get(&200), Some(&1)); // TC == rc(AG) → one ALT
        assert!(g.iter().all(|s| s.contig == "chr1")); // canonical CHM13 contig
    }

    #[test]
    fn resolve_chip_ref_alt_swapped_against_chm13() {
        // The asset's GRCh37 locus is ref/alt-SWAPPED vs CHM13: chm13 chr1:100 G/T (ALT=T) but
        // grch37 1:500 T/G (ALT=G). A chip hom for G is hom-CHM13-REF → dosage 0. Scoring against the
        // build ALT (G) would wrongly give 2; scoring against the CHM13 ALT (T) gives the correct 0.
        let (panel, _) = IbdPanel::from_sites(
            "chm13v2.0",
            vec![site("rs_swap", (100, 'G', 'T'), Some((500, 'T', 'G')))],
        );
        let g = panel.resolve_chip("GRCh37", &[("1".to_string(), 500, 'G', 'G')]);
        assert_eq!(g.len(), 1);
        assert_eq!(g[0].dosage, 0); // hom CHM13-ref, not 2
        assert_eq!(
            (g[0].reference_allele.as_str(), g[0].alternate_allele.as_str()),
            ("G", "T")
        ); // canonical CHM13 alleles
           // The other homozygote (T/T) is hom CHM13-ALT → dosage 2.
        assert_eq!(
            panel.resolve_chip("GRCh37", &[("1".to_string(), 500, 'T', 'T')])[0].dosage,
            2
        );
    }

    #[test]
    fn resolve_chip_hom_alt_and_unknown_build() {
        let (panel, _) = IbdPanel::from_sites("chm13v2.0", vec![site("rs1", (100, 'A', 'G'), Some((500, 'A', 'G')))]);
        let calls = vec![("1".to_string(), 500, 'G', 'G')]; // hom-alt → dosage 2
        assert_eq!(panel.resolve_chip("GRCh37", &calls)[0].dosage, 2);
        // A build with no loci in the panel resolves nothing.
        assert!(panel.resolve_chip("GRCh38", &calls).is_empty());
    }

    #[test]
    fn resolve_whole_genome_densifies_hom_ref_and_skips_palindromes() {
        // Panel: rs1 (non-palindromic, GRCh37 strand-flipped), rs2 (non-palindromic), rs3 (palindrome).
        let (panel, _) = IbdPanel::from_sites(
            "chm13v2.0",
            vec![
                site("rs1", (100, 'A', 'G'), Some((500, 'T', 'C'))), // GRCh37 alleles are rc of CHM13
                site("rs2", (200, 'A', 'G'), Some((600, 'A', 'G'))),
                site("rs3", (300, 'A', 'T'), Some((700, 'A', 'T'))), // palindrome — always skipped
            ],
        );
        // A whole-genome source that lists ONLY rs1 as a het (GRCh37 forward alleles T/C, contig
        // "chr1" to prove chr-insensitive matching). rs2 is unlisted ⇒ hom-ref; rs3 skipped.
        let calls = vec![("chr1".to_string(), 500, 'T', 'C')];
        let g = panel.resolve_whole_genome("GRCh37", &calls);
        let by_pos: std::collections::HashMap<i64, i32> = g.iter().map(|s| (s.position, s.dosage)).collect();
        assert_eq!(g.len(), 2, "rs1 + rs2 densified; the palindrome rs3 is excluded");
        assert_eq!(by_pos.get(&100), Some(&1)); // rs1 het via rc reconciliation
        assert_eq!(by_pos.get(&200), Some(&0)); // rs2 unlisted ⇒ homozygous reference
        assert!(!by_pos.contains_key(&300)); // palindrome never emitted
    }

    #[test]
    fn resolve_whole_genome_hom_alt_and_non_reconciling() {
        let (panel, _) = IbdPanel::from_sites(
            "chm13v2.0",
            vec![
                site("rs1", (100, 'A', 'G'), Some((500, 'A', 'G'))),
                site("rs2", (200, 'C', 'T'), Some((600, 'C', 'T'))),
            ],
        );
        // rs1 hom-alt (G/G → dosage 2); rs2 listed but the alleles are internally inconsistent with
        // the biallelic site: C matches ref directly, A only matches alt(T) under rc — neither a
        // pure direct nor a pure rc pair, so it doesn't reconcile → dropped, NOT called hom-ref.
        let calls = vec![("1".to_string(), 500, 'G', 'G'), ("1".to_string(), 600, 'C', 'A')];
        let g = panel.resolve_whole_genome("GRCh37", &calls);
        let by_pos: std::collections::HashMap<i64, i32> = g.iter().map(|s| (s.position, s.dosage)).collect();
        assert_eq!(by_pos.get(&100), Some(&2));
        assert_eq!(by_pos.get(&200), None); // non-reconciling variant dropped, not hom-ref
    }

    #[test]
    fn round_trips_through_bincode() {
        let (panel, _) = IbdPanel::from_sites("chm13v2.0", vec![site("rs1", (100, 'A', 'G'), Some((500, 'A', 'G')))]);
        assert_eq!(IbdPanel::from_bytes(&panel.to_bytes().unwrap()).unwrap(), panel);
    }
}
