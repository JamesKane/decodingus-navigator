//! The Y-SNP name → locus dictionary that gives a BISDNA export, or any name-only Y panel, its
//! missing coordinates. A SNP name like `CTS10003` resolves to a position, plus its ancestral and
//! derived alleles, **for each reference build**. The codebase then stays build-agnostic. The
//! build label is the key of `coordinates` (`"GRCh38"`, `"GRCh37"`, `"hs1"`, …), which is exactly
//! the convention the DecodingUs Y-tree uses. The caller gives the importer the build it places
//! against, and the importer reads that coordinate. Nothing here is specific to CHM13.
//!
//! The bulk data is a generated asset, which `scripts/ysnp-dictionary/` builds from YBrowse and a
//! liftover. A small checked-in chromo2 panel manifest uses the same format. This module is pure
//! over text that is already in memory ([`YsnpDictionary::from_text`]), with a thin
//! [`YsnpDictionary::load`] IO boundary that reads the asset files. See
//! `documents/design/bisdna-import.md`.

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

/// The locus of one SNP on a specific reference build. The alleles are on the + strand of that
/// build, so a liftover that flips the strand stores its own complemented alleles. The alleles
/// belong to a coordinate, and not to a SNP.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Coord {
    pub chrom: String,
    pub position: i64,
    /// `+`/`-` orientation of the SNP on this build (for the importer's strand cross-check).
    pub strand: char,
    pub ancestral: String,
    pub derived: String,
}

/// A SNP and its coordinates across builds (canonical name + a build-keyed coordinate map).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SnpEntry {
    /// Canonical name, original case (e.g. `CTS10003`).
    pub name: String,
    /// build label → coordinate on that build.
    pub coordinates: HashMap<String, Coord>,
}

/// A successful resolution: the canonical name plus the coordinate on the requested build.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResolvedSnp<'a> {
    pub canonical: &'a str,
    pub coord: &'a Coord,
}

/// The loaded dictionary: canonical entries, plus an alias → canonical index. A lookup ignores
/// the case of the SNP name, and matches a build key verbatim.
#[derive(Debug, Clone, Default)]
pub struct YsnpDictionary {
    /// lowercased canonical name → entry.
    by_name: HashMap<String, SnpEntry>,
    /// lowercased alias → lowercased canonical name.
    alias_to_canonical: HashMap<String, String>,
}

/// Dictionary asset dir: `$NAVIGATOR_YSNP_DIR`, else `$NAVIGATOR_REFGENOME_DIR/ysnp`, else
/// `~/.decodingus/ysnp` (mirrors the cache convention in `navigator-refgenome`).
pub fn asset_dir() -> PathBuf {
    if let Some(dir) = std::env::var_os("NAVIGATOR_YSNP_DIR") {
        return PathBuf::from(dir);
    }
    if let Some(base) = std::env::var_os("NAVIGATOR_REFGENOME_DIR") {
        return PathBuf::from(base).join("ysnp");
    }
    crate::paths::decodingus_dir().join("ysnp")
}

/// Split a TSV line into trimmed cells, and drop an empty last cell that a final tab makes.
fn cells(line: &str) -> Vec<&str> {
    line.split('\t').map(str::trim).collect()
}

/// True for a `#`-comment or blank line.
fn is_skippable(line: &str) -> bool {
    let t = line.trim();
    t.is_empty() || t.starts_with('#')
}

impl YsnpDictionary {
    /// Build from the two asset texts (no IO).
    ///
    /// A `dictionary` row is
    /// `name<TAB>build<TAB>chrom<TAB>position<TAB>strand<TAB>ancestral<TAB>derived`. An `aliases`
    /// row is `alias<TAB>canonical`, and that file is optional and can be empty.
    ///
    /// This drops a first header row whose first cell is `name` or `alias`. It also drops a `#`
    /// comment, a blank line, and a row whose position does not parse. For one (name, build), the
    /// first coordinate wins, and later duplicates go. That is deterministic over a sorted asset.
    /// Errors only if no usable entry comes out.
    pub fn from_text(dictionary: &str, aliases: &str) -> Result<Self, String> {
        let mut by_name: HashMap<String, SnpEntry> = HashMap::new();
        for line in dictionary.lines() {
            if is_skippable(line) {
                continue;
            }
            let c = cells(line);
            if c.len() < 7 {
                continue;
            }
            if c[0].eq_ignore_ascii_case("name") {
                continue; // header row
            }
            let (name, build) = (c[0], c[1]);
            let Ok(position) = c[3].parse::<i64>() else { continue };
            let strand = c[4].chars().next().unwrap_or('+');
            let coord = Coord {
                chrom: c[2].to_string(),
                position,
                strand,
                ancestral: c[5].to_ascii_uppercase(),
                derived: c[6].to_ascii_uppercase(),
            };
            let entry = by_name.entry(name.to_ascii_lowercase()).or_insert_with(|| SnpEntry {
                name: name.to_string(),
                coordinates: HashMap::new(),
            });
            entry.coordinates.entry(build.to_string()).or_insert(coord);
        }

        let mut alias_to_canonical = HashMap::new();
        for line in aliases.lines() {
            if is_skippable(line) {
                continue;
            }
            let c = cells(line);
            if c.len() < 2 || c[0].is_empty() || c[1].is_empty() {
                continue;
            }
            if c[0].eq_ignore_ascii_case("alias") {
                continue; // header row
            }
            // Only index aliases that point at a known canonical entry, and never shadow a
            // real canonical name.
            let (alias, canonical) = (c[0].to_ascii_lowercase(), c[1].to_ascii_lowercase());
            if by_name.contains_key(&canonical) && !by_name.contains_key(&alias) {
                alias_to_canonical.entry(alias).or_insert(canonical);
            }
        }

        if by_name.is_empty() {
            return Err(
                "Y-SNP dictionary is empty (no `name<TAB>build<TAB>chrom<TAB>pos<TAB>strand<TAB>anc<TAB>der` rows)"
                    .into(),
            );
        }
        Ok(YsnpDictionary {
            by_name,
            alias_to_canonical,
        })
    }

    /// Candidate dictionary filenames, in the preference order of `load`. The full ~200 MB catalog
    /// of ~2M names comes first, then the small panel of one chip as a fallback. The chromo2 chip
    /// panel is a stale subset of ~14k names, and it would hide current names that the full
    /// catalog holds. So the catalog wins whenever the machine has it, and it is the one the app
    /// downloads on first use.
    pub const ASSET_FILENAMES: &'static [&'static str] = &["dictionary.tsv", "chromo2-panel.tsv"];

    /// Read the asset from `dir`: the first of [`Self::ASSET_FILENAMES`] that exists, plus an
    /// optional sibling `aliases.tsv`. It prefers the full catalog, for the widest and most
    /// current name coverage. It uses the chromo2 panel only when the catalog is absent.
    pub fn load(dir: &Path) -> Result<Self, String> {
        let dict_path = Self::ASSET_FILENAMES
            .iter()
            .map(|f| dir.join(f))
            .find(|p| p.is_file())
            .ok_or_else(|| {
                format!(
                    "no Y-SNP dictionary in {} (looked for {})",
                    dir.display(),
                    Self::ASSET_FILENAMES.join(", ")
                )
            })?;
        let dictionary =
            std::fs::read_to_string(&dict_path).map_err(|e| format!("reading {}: {e}", dict_path.display()))?;
        let aliases = std::fs::read_to_string(dir.join("aliases.tsv")).unwrap_or_default();
        Self::from_text(&dictionary, &aliases)
    }

    /// Resolve a SNP `name` to its coordinate on `build`. Case-insensitive on the name;
    /// follows one alias hop to the canonical entry. `None` if the name is unknown or the
    /// entry has no coordinate on `build`.
    pub fn resolve(&self, name: &str, build: &str) -> Option<ResolvedSnp<'_>> {
        let key = name.trim().to_ascii_lowercase();
        let canonical_key = self.alias_to_canonical.get(&key).map(String::as_str).unwrap_or(&key);
        let entry = self.by_name.get(canonical_key)?;
        let coord = entry.coordinates.get(build)?;
        Some(ResolvedSnp {
            canonical: &entry.name,
            coord,
        })
    }

    /// Build a reverse index `position → canonical name` for one reference `build`. This is the
    /// inverse of [`resolve`](Self::resolve). A position-only call is a novel or private Y
    /// variant. This index lets a caller add the catalogued Y-SNP name at that site, when such a
    /// name exists.
    ///
    /// At one position the first name wins, which is deterministic over a sorted asset. This
    /// leaves out a position that `build` does not have. In practice every entry is chrY, so the
    /// key is the position alone. A caller that resolves the correct build avoids the cross-build
    /// integer collision, which is in any case very improbable.
    pub fn position_index(&self, build: &str) -> HashMap<i64, &str> {
        let mut idx = HashMap::new();
        for entry in self.by_name.values() {
            if let Some(c) = entry.coordinates.get(build) {
                idx.entry(c.position).or_insert(entry.name.as_str());
            }
        }
        idx
    }

    /// Number of canonical SNP entries.
    pub fn len(&self) -> usize {
        self.by_name.len()
    }

    pub fn is_empty(&self) -> bool {
        self.by_name.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const DICT: &str = "\
name\tbuild\tchrom\tposition\tstrand\tancestral\tderived
CTS10003\tGRCh38\tchrY\t15000000\t+\tc\tt
CTS10003\ths1\tchrY\t14800000\t+\tC\tT
M269\tGRCh38\tchrY\t22739367\t+\tT\tC
M269\ths1\tchrY\t21200000\t-\tA\tG
OnlyHg38\tGRCh38\tchrY\t9999999\t+\tA\tG
";

    const ALIASES: &str = "\
alias\tcanonical
PF6517\tM269
S163\tNoSuchSnp
M269\tCTS10003
";

    fn dict() -> YsnpDictionary {
        YsnpDictionary::from_text(DICT, ALIASES).unwrap()
    }

    #[test]
    fn resolves_per_build_with_distinct_coords() {
        let d = dict();
        let g38 = d.resolve("CTS10003", "GRCh38").unwrap();
        assert_eq!(g38.coord.position, 15000000);
        assert_eq!((g38.coord.ancestral.as_str(), g38.coord.derived.as_str()), ("C", "T")); // upcased
        let hs1 = d.resolve("CTS10003", "hs1").unwrap();
        assert_eq!(hs1.coord.position, 14800000); // different build, different position
    }

    #[test]
    fn name_lookup_is_case_insensitive_canonical_kept() {
        let d = dict();
        let r = d.resolve("cts10003", "hs1").unwrap();
        assert_eq!(r.canonical, "CTS10003"); // original case preserved
    }

    #[test]
    fn load_prefers_full_dictionary_over_chromo2_panel() {
        // Both present → the full `dictionary.tsv` wins; the stale ~14k-name chromo2 chip panel must
        // not shadow current names in the full catalog. With only the panel, it is the fallback.
        let dir = std::env::temp_dir().join(format!("dun-ysnp-pref-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        let hdr = "name\tbuild\tchrom\tposition\tstrand\tancestral\tderived\n";
        std::fs::write(
            dir.join("dictionary.tsv"),
            format!("{hdr}FullOnlySnp\ths1\tchrY\t123\t+\tA\tG\n"),
        )
        .unwrap();
        std::fs::write(
            dir.join("chromo2-panel.tsv"),
            format!("{hdr}PanelOnlySnp\ths1\tchrY\t456\t+\tA\tG\n"),
        )
        .unwrap();

        let d = YsnpDictionary::load(&dir).unwrap();
        assert!(d.resolve("FullOnlySnp", "hs1").is_some(), "loaded the full catalog");
        assert!(
            d.resolve("PanelOnlySnp", "hs1").is_none(),
            "did not load the chromo2 panel"
        );

        std::fs::remove_file(dir.join("dictionary.tsv")).unwrap();
        let d2 = YsnpDictionary::load(&dir).unwrap();
        assert!(
            d2.resolve("PanelOnlySnp", "hs1").is_some(),
            "fell back to the chromo2 panel"
        );

        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn alias_resolves_to_canonical() {
        let d = dict();
        // PF6517 is an alias of M269, so resolve through the alias.
        let r = d.resolve("PF6517", "GRCh38").unwrap();
        assert_eq!(r.canonical, "M269");
        assert_eq!(r.coord.position, 22739367);
    }

    #[test]
    fn alias_to_unknown_canonical_is_ignored() {
        let d = dict();
        // S163 -> NoSuchSnp (not a real entry). The code drops the alias, so S163 does not
        // resolve.
        assert!(d.resolve("S163", "GRCh38").is_none());
    }

    #[test]
    fn alias_never_shadows_a_real_canonical_name() {
        let d = dict();
        // The aliases file maps M269 -> CTS10003, but M269 is itself canonical: the real
        // entry must win, not the alias.
        let r = d.resolve("M269", "GRCh38").unwrap();
        assert_eq!(r.canonical, "M269");
        assert_eq!(r.coord.position, 22739367);
    }

    #[test]
    fn missing_build_yields_none() {
        let d = dict();
        assert!(d.resolve("OnlyHg38", "hs1").is_none()); // no hs1 coordinate (didn't lift)
        assert!(d.resolve("OnlyHg38", "GRCh38").is_some());
    }

    #[test]
    fn strand_is_captured() {
        let d = dict();
        assert_eq!(d.resolve("M269", "hs1").unwrap().coord.strand, '-');
    }

    #[test]
    fn position_index_reverse_lookup_per_build() {
        let d = dict();
        let hs1 = d.position_index("hs1");
        assert_eq!(hs1.get(&14800000).copied(), Some("CTS10003")); // canonical case preserved
        assert_eq!(hs1.get(&21200000).copied(), Some("M269"));
        assert_eq!(hs1.get(&15000000), None); // that's the GRCh38 position, not hs1
        let g38 = d.position_index("GRCh38");
        assert_eq!(g38.get(&15000000).copied(), Some("CTS10003"));
        assert_eq!(g38.get(&9999999).copied(), Some("OnlyHg38"));
    }

    #[test]
    fn empty_dictionary_errors() {
        assert!(YsnpDictionary::from_text("# only comments\n", "").is_err());
    }

    #[test]
    fn tolerates_missing_aliases_and_short_rows() {
        let d = YsnpDictionary::from_text("CTS1\ths1\tchrY\t100\t+\tA\tG\nbadrow\tonly\ttwo\n", "").unwrap();
        assert_eq!(d.len(), 1);
        assert!(d.resolve("CTS1", "hs1").is_some());
    }
}
