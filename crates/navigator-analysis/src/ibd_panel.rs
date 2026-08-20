//! The IBD reference panel. It covers more than one build, and a chip can use it. See the
//! ancestry-ibd asset design, B2 and B2c.
//!
//! IBD matching needs a neutral, dense set of SNPs that a **consumer array also assays**. Chip kits
//! outnumber WGS runs by orders of magnitude, so the panel must sit where a chip and a WGS run
//! overlap.
//!
//! Each site carries its `(contig, pos, REF, ALT)` on **CHM13, GRCh37 and GRCh38**. An offline
//! GATK liftover that knows about alleles builds those once. So a chip genotype on *any* build
//! resolves to the canonical CHM13 site, with its orientation. The app needs **no liftover at run
//! time**, because the panel already holds the answer.
//!
//! Two rules keep this correct:
//!
//! - The locus of each build carries the **same biological alleles**. On an inverted chain block,
//!   GATK takes the reverse complement, and it exchanges REF with ALT. The count of the ALT of
//!   that build then equals the count of the CHM13 ALT. The dosage does not depend on the
//!   build.
//! - The panel **leaves out a palindrome that is ambiguous about the strand**, which is A/T or
//!   C/G. See [`is_palindromic`]. There `rc(A)=T` is also a correct allele, so a comparison of the
//!   alleles can not tell you the strand of a chip.

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

    /// Build from a set of sites, and **keep** a palindrome that is ambiguous about the strand,
    /// which is A/T or C/G.
    ///
    /// The panel is a superset of probes. A WGS run and the ancestry path genotype a palindromic
    /// site without trouble, because a read gives the base on the reference strand. The CHIP path
    /// alone can not orient one. So [`resolve_chip`] skips a palindrome when it resolves, and this
    /// function does not remove one from the panel.
    ///
    /// Returns `(panel, n_palindromic)`, where the second value is the count of palindromes that
    /// stayed, for the build log.
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

    /// Resolve the chip calls, which are on `build` and come as `(contig, pos, a1, a2)`, to
    /// canonical CHM13 dosages.
    ///
    /// It indexes the panel by the `(contig, position)` of that build. It then counts the copies of
    /// the **canonical CHM13 ALT** straight from each observed pair, either directly or through the
    /// reverse complement. See [`dosage_from_alleles`]. It emits the result as a [`SiteGenotype`]
    /// at the CHM13 locus. There is no liftover at run time, and no alignment. It drops a call that
    /// does not match a site, a no-call, and a call that does not reconcile.
    ///
    /// The score goes against the **CHM13** `(REF, ALT)`, and not against the locus of the build.
    /// That is deliberate. The allele letters of a chip are absolute. But the `(REF, ALT)` labels
    /// of each build, in the asset, do not reliably point to the CHM13 ALT. A large share of them
    /// have REF and ALT the other way round, and the reference allele of GRCh37 is often the CHM13
    /// ALT. At those sites, a score against the ALT of the build turns the dosage from 0 to 2, and
    /// from 2 to 0.
    ///
    /// A comparison of the chip alleles against the CHM13 alleles, with the reverse-complement
    /// retry for the strand, can not go wrong that way. The locus of the build has one use: to look
    /// the site up by position.
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
            // Take a palindrome that is ambiguous about the strand, which is A/T or C/G. The
            // reported alleles of a chip give it no orientation. So skip it on the chip path. A
            // WGS run and the ancestry path still use it, from a direct base call. The probe panel
            // keeps it, and this line is where the chip path alone leaves it out.
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

    /// Resolve a source that covers the **whole genome and lists variants alone** to canonical
    /// CHM13 dosages over the panel. A WGS VCF and a CompleteGenomics masterVar are such
    /// sources.
    ///
    /// A chip reports a genotype at every array site. Such a source instead lists *only* the sites
    /// that are not reference. So every panel site that the source could have called, and did not,
    /// counts as **homozygous reference**, at dosage 0.
    ///
    /// That assumption holds **only** for a source that genotyped the whole genome, where an
    /// absent site means hom-ref and not a no-call. Never give this function a targeted panel,
    /// such as a Big Y or a Sanger run.
    ///
    /// `variant_calls` holds the variant sites of the source on `build`, as `(contig, pos, a1, a2)`
    /// allele pairs, forward on the reference. A contig matches whatever `chr` prefix it carries,
    /// so a `chr1` from the source lines up with a `grch37` panel locus stored as `1`.
    ///
    /// A variant whose alleles do not reconcile to the site, which is a mismatch at a site with
    /// more than two alleles, goes away. The code does not call it hom-ref by mistake. It also
    /// skips a palindromic site, A/T or C/G, because that site is ambiguous about the strand across
    /// builds, exactly as in [`resolve_chip`].
    pub fn resolve_whole_genome(&self, build: &str, variant_calls: &[(String, i64, char, char)]) -> Vec<SiteGenotype> {
        let norm = crate::contig::bare_upper;
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

    /// Take genotypes that the caller made at the loci of **this build**, and key them back to
    /// canonical CHM13 dosages.
    ///
    /// The caller genotypes the BAM of an alignment that is not on CHM13. It works at the
    /// `locus(build)` of each site, which holds the contig, the position, the REF and the ALT of
    /// that build. The `dosage` that comes back counts the ALT of the *build*.
    ///
    /// The REF and ALT of a build are often the other way round from CHM13. See [`resolve_chip`].
    /// So this function builds the observed alleles again from that dosage. It scores them against
    /// the **CHM13** REF and ALT, either directly or through the reverse complement. It emits each
    /// one at its CHM13 locus, with the depth and the GQ of the alignment.
    ///
    /// It skips a palindrome, A/T or C/G. That site is ambiguous about the strand across builds,
    /// exactly as in the chip resolver and the whole-genome one.
    ///
    /// This is the counterpart of [`resolve_chip`] for an alignment. It lets a GRCh37 or GRCh38 WGS
    /// run reach the ancestry panel, whose coordinates are CHM13, with no liftover at run time. The
    /// panel already carries the coordinates of every build.
    pub fn resolve_alignment(&self, build: &str, genotypes: &[SiteGenotype]) -> Vec<SiteGenotype> {
        let norm = crate::contig::bare_upper;
        let mut index: HashMap<(String, i64), &IbdPanelSite> = HashMap::new();
        for s in &self.sites {
            if let Some(l) = s.locus(build) {
                index.insert((norm(&l.contig), l.position), s);
            }
        }
        let mut out = Vec::new();
        for g in genotypes {
            if g.dosage < 0 {
                continue; // no-call
            }
            let Some(site) = index.get(&(norm(&g.contig), g.position)) else {
                continue;
            };
            if is_palindromic(site.chm13.reference, site.chm13.alternate) {
                continue;
            }
            // The caller made this genotype against the REF and ALT of the build, which are
            // g.reference_allele and g.alternate_allele. Build the observed diploid alleles again
            // from the dosage, and then score them against the CHM13 alleles.
            let br = g.reference_allele.chars().next().unwrap_or('N');
            let ba = g.alternate_allele.chars().next().unwrap_or('N');
            let (a1, a2) = match g.dosage {
                0 => (br, br),
                1 => (br, ba),
                _ => (ba, ba),
            };
            let Some(dosage) = dosage_from_alleles(a1, a2, site.chm13.reference, site.chm13.alternate) else {
                continue;
            };
            out.push(SiteGenotype {
                name: site.rsid.clone(),
                contig: site.chm13.contig.clone(),
                position: site.chm13.position,
                reference_allele: site.chm13.reference.to_string(),
                alternate_allele: site.chm13.alternate.to_string(),
                ploidy: 2,
                dosage,
                gq: g.gq,
                depth: g.depth,
                ref_depth: 0,
                alt_depth: 0,
                pls: Vec::new(),
                gt: None,
                allele_depths: None,
            });
        }
        out
    }

    /// The canonical CHM13 `(contig, position)` sites. A WGS caller genotypes these targets, and
    /// its dosages then line up with those of the chip path.
    pub fn chm13_sites(&self) -> Vec<(&str, i64)> {
        self.sites
            .iter()
            .map(|s| (s.chm13.contig.as_str(), s.chm13.position))
            .collect()
    }
}

/// Build a diploid [`SiteGenotype`] at the canonical CHM13 locus of a panel site, with the given
/// alt dosage. The chip resolver and the whole-genome one share this. It carries no depth and no
/// quality, and the dosage alone.
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

/// True when `(a, b)` is a palindrome that is ambiguous about the strand, which is A/T or C/G. A
/// panel that a chip can use leaves those out, because the reverse complement can not tell you the
/// strand of the array.
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
        // rs1: GRCh37 1:500 A/G, which holds the same alleles as CHM13 chr1:100 A/G.
        // rs2: GRCh37 1:600 T/C, against CHM13 chr1:200 A/G. That is a strand flip, and the GRCh37
        // alleles are the reverse complement.
        let (panel, _) = IbdPanel::from_sites(
            "chm13v2.0",
            vec![
                site("rs1", (100, 'A', 'G'), Some((500, 'A', 'G'))),
                site("rs2", (200, 'A', 'G'), Some((600, 'T', 'C'))),
            ],
        );
        // A chip on GRCh37. rs1 is het AG, so the dosage is 1. rs2 is het TC. That reconciles
        // through the reverse complement, so its dosage is 1 too.
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
        // The GRCh37 locus of the asset has REF and ALT the OTHER WAY ROUND from CHM13. CHM13
        // chr1:100 is G/T, with ALT=T. GRCh37 1:500 is T/G, with ALT=G. A chip that is hom for G is
        // then hom for the CHM13 REF, and the dosage is 0. A score against the ALT of the build,
        // which is G, would wrongly give 2. A score against the CHM13 ALT, which is T, gives the
        // correct 0.
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
        // A whole-genome source that lists ONLY rs1, as a het. Its GRCh37 forward alleles are
        // T/C. Its contig is "chr1", which shows that the match ignores the `chr` prefix. rs2 is
        // not in the list, so it is hom-ref. The code skips rs3.
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
        // rs1 is hom-alt, G/G, so its dosage is 2.
        //
        // The rs2 record is in the list, but its alleles do not agree with the site, which has two
        // alleles. C matches the ref directly. A matches the
        // alt, T, only under the reverse complement. The pair is then neither a pure direct match
        // nor a pure reverse-complement one. It does not reconcile. The code drops it, and it does
        // NOT call it hom-ref.
        let calls = vec![("1".to_string(), 500, 'G', 'G'), ("1".to_string(), 600, 'C', 'A')];
        let g = panel.resolve_whole_genome("GRCh37", &calls);
        let by_pos: std::collections::HashMap<i64, i32> = g.iter().map(|s| (s.position, s.dosage)).collect();
        assert_eq!(by_pos.get(&100), Some(&2));
        assert_eq!(by_pos.get(&200), None); // non-reconciling variant dropped, not hom-ref
    }

    // A panel site with an explicit build-locus contig (e.g. "chr1" for the b38 column vs bare "1").
    fn site_b(rsid: &str, chm13: (i64, char, char), build: (&str, i64, char, char), which: &str) -> IbdPanelSite {
        let locus = Locus {
            contig: build.0.into(),
            position: build.1,
            reference: build.2,
            alternate: build.3,
        };
        IbdPanelSite {
            rsid: rsid.into(),
            chm13: Locus {
                contig: "chr1".into(),
                position: chm13.0,
                reference: chm13.1,
                alternate: chm13.2,
            },
            grch37: (which == "37").then(|| locus.clone()),
            grch38: (which == "38").then_some(locus),
        }
    }

    // A SiteGenotype as `caller::genotype_sites_all_contigs` would produce it at a build locus.
    fn geno(name: &str, contig: &str, pos: i64, r: &str, a: &str, dosage: i32) -> SiteGenotype {
        SiteGenotype {
            name: name.into(),
            contig: contig.into(),
            position: pos,
            reference_allele: r.into(),
            alternate_allele: a.into(),
            ploidy: 2,
            dosage,
            gq: 30,
            depth: 20,
            ref_depth: 0,
            alt_depth: 0,
            pls: Vec::new(),
            gt: None,
            allele_depths: None,
        }
    }

    #[test]
    fn resolve_alignment_rekeys_swap_strand_and_skips_palindrome() {
        // rs1: same orientation; rs_swap: grch38 REF/ALT swapped vs CHM13; rs_pal: palindrome.
        let sites = vec![
            site_b("rs1", (100, 'A', 'G'), ("chr1", 500, 'A', 'G'), "38"),
            site_b("rs_swap", (200, 'G', 'T'), ("chr1", 600, 'T', 'G'), "38"),
            site_b("rs_pal", (300, 'A', 'T'), ("chr1", 700, 'A', 'T'), "38"),
        ];
        let (panel, _) = IbdPanel::from_sites("chm13v2.0", sites);
        // The genotypes at the grch38 loci, whose dosage points at the build.
        //
        // rs1 is het, so the CHM13 dosage is 1. rs_swap is hom for the grch38 ALT, which is G/G. G
        // is the CHM13 REF, so the CHM13 dosage is 0.
        let raw = vec![
            geno("rs1", "chr1", 500, "A", "G", 1),
            geno("rs_swap", "chr1", 600, "T", "G", 2),
            geno("rs_pal", "chr1", 700, "A", "T", 1),
        ];
        let out = panel.resolve_alignment("GRCh38", &raw);
        let by_pos: HashMap<i64, i32> = out.iter().map(|s| (s.position, s.dosage)).collect();
        assert_eq!(by_pos.get(&100), Some(&1));
        assert_eq!(by_pos.get(&200), Some(&0), "grch38 G/G == CHM13 hom-ref, not hom-alt");
        assert!(!by_pos.contains_key(&300), "palindrome skipped");
        // Emitted at CHM13 loci with CHM13 alleles; depth preserved.
        let swap = out.iter().find(|s| s.position == 200).unwrap();
        assert_eq!(
            (swap.reference_allele.as_str(), swap.alternate_allele.as_str()),
            ("G", "T")
        );
        assert_eq!(swap.depth, 20);
        assert!(out.iter().all(|s| s.contig == "chr1"));
    }

    #[test]
    fn resolve_alignment_matches_contig_chr_insensitively() {
        // Panel grch37 locus stored as "chr1"; a b37 BAM genotyped with bare "1" still resolves.
        let sites = vec![site_b("rs1", (100, 'A', 'G'), ("chr1", 500, 'A', 'G'), "37")];
        let (panel, _) = IbdPanel::from_sites("chm13v2.0", sites);
        let out = panel.resolve_alignment("GRCh37", &[geno("rs1", "1", 500, "A", "G", 2)]);
        assert_eq!(out.len(), 1);
        assert_eq!(
            (out[0].position, out[0].dosage, out[0].contig.as_str()),
            (100, 2, "chr1")
        );
        // The code drops a no-call, which is a dosage below 0.
        assert!(panel
            .resolve_alignment("GRCh37", &[geno("rs1", "1", 500, "A", "G", -1)])
            .is_empty());
    }

    #[test]
    fn round_trips_through_bincode() {
        let (panel, _) = IbdPanel::from_sites("chm13v2.0", vec![site("rs1", (100, 'A', 'G'), Some((500, 'A', 'G')))]);
        assert_eq!(IbdPanel::from_bytes(&panel.to_bytes().unwrap()).unwrap(), panel);
    }
}
