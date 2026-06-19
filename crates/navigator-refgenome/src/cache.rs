//! On-disk cache layout. Mirrors the `~/.decodingus` convention (and the `NAVIGATOR_TREE_DIR`
//! override pattern) used elsewhere: references live under `<base>/references/`, liftover
//! chains under `<base>/liftover/`.

use std::path::{Path, PathBuf};

use crate::registry::Build;

/// Cache root: `$NAVIGATOR_REFGENOME_DIR`, else `~/.decodingus`, else the current dir.
pub fn base_dir() -> PathBuf {
    if let Some(dir) = std::env::var_os("NAVIGATOR_REFGENOME_DIR") {
        return PathBuf::from(dir);
    }
    let home = std::env::var("HOME").unwrap_or_else(|_| ".".to_string());
    PathBuf::from(home).join(".decodingus")
}

/// The decompressed, indexed reference FASTA path for a build.
pub fn reference_path(base: &Path, build: Build) -> PathBuf {
    base.join("references").join(format!("{}.fa", build.as_str()))
}

/// The companion `.fai` index path.
pub fn reference_fai(base: &Path, build: Build) -> PathBuf {
    base.join("references").join(format!("{}.fa.fai", build.as_str()))
}

/// The cached liftover chain path for a build pair.
pub fn chain_path(base: &Path, from: Build, to: Build) -> PathBuf {
    base.join("liftover").join(format!("{}-to-{}.chain", from.as_str(), to.as_str()))
}

/// The cached annotation-mask BED path for a named mask (e.g. the curated CHM13 Y palindrome /
/// amplicon BEDs). Stored under `<base>/masks/<name>.bed`.
pub fn mask_path(base: &Path, name: &str) -> PathBuf {
    base.join("masks").join(format!("{name}.bed"))
}

/// The parsed genome-regions JSON for a build, under `<base>/regions/<build>.json`.
pub fn regions_path(base: &Path, build: Build) -> PathBuf {
    base.join("regions").join(format!("{}.json", build.as_str()))
}

/// Age of a cached file in days (for TTL checks); `None` if it doesn't exist or its mtime is
/// unreadable / in the future.
pub fn age_days(path: &Path) -> Option<f64> {
    let modified = std::fs::metadata(path).ok()?.modified().ok()?;
    let elapsed = std::time::SystemTime::now().duration_since(modified).ok()?;
    Some(elapsed.as_secs_f64() / 86_400.0)
}

/// Whether `path` exists and is non-empty (the cache-hit predicate).
pub fn is_present(path: &Path) -> bool {
    std::fs::metadata(path).map(|m| m.len() > 0).unwrap_or(false)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn paths_are_under_the_base() {
        let base = Path::new("/tmp/dun");
        assert_eq!(reference_path(base, Build::Chm13v2), Path::new("/tmp/dun/references/chm13v2.0.fa"));
        assert_eq!(reference_fai(base, Build::Chm13v2), Path::new("/tmp/dun/references/chm13v2.0.fa.fai"));
        assert_eq!(chain_path(base, Build::Grch38, Build::Chm13v2), Path::new("/tmp/dun/liftover/GRCh38-to-chm13v2.0.chain"));
    }

    #[test]
    fn base_dir_honors_env_override() {
        // Safe: this test only reads the var via base_dir; set+remove around the assertion.
        std::env::set_var("NAVIGATOR_REFGENOME_DIR", "/tmp/refcache");
        assert_eq!(base_dir(), Path::new("/tmp/refcache"));
        std::env::remove_var("NAVIGATOR_REFGENOME_DIR");
    }
}
