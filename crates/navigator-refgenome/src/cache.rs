//! On-disk cache layout. Mirrors the `~/.decodingus` convention (and the `NAVIGATOR_TREE_DIR`
//! override pattern) used elsewhere: references live under `<base>/references/`, liftover
//! chains under `<base>/liftover/`.

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use crate::registry::Build;

/// Write `bytes` to `path` **atomically**: write a uniquely-named sibling temp file, flush it to
/// disk, then `rename` it over the target. `rename` is atomic on a POSIX filesystem, so a reader —
/// or a concurrent/ crashing writer — never sees a torn, half-written, or head-of-new + tail-of-old
/// file. This is the safe replacement for `std::fs::write` on **config** files (`reference_sources.json`,
/// `settings.json`), which are written from spawned tasks that can race: a plain `fs::write` truncates
/// then streams, so two racing writes interleave into a corrupt mix (the classic short-head/long-tail
/// garbage). The parent dir is created if absent. A stale temp from a crashed write is harmless (never
/// read; overwritten by name reuse or left as an obvious `*.tmp.*`).
pub fn atomic_write(path: &Path, bytes: &[u8]) -> std::io::Result<()> {
    use std::io::Write;
    static COUNTER: AtomicU64 = AtomicU64::new(0);
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    // Unique temp name **in the same directory** (rename must stay on one filesystem). pid + a
    // process-wide counter keeps two concurrent writers from sharing — and clobbering — a temp file.
    let uniq = format!("{}.{}", std::process::id(), COUNTER.fetch_add(1, Ordering::Relaxed));
    let tmp = path.with_extension(format!("tmp.{uniq}"));
    let res = (|| {
        let mut f = std::fs::File::create(&tmp)?;
        f.write_all(bytes)?;
        f.sync_all()?; // durability: don't expose an empty/partial file after a crash
        drop(f); // close before renaming — Windows will not move a handle it still owns
        retry_while_replacing(|| std::fs::rename(&tmp, path))
    })();
    if res.is_err() {
        let _ = std::fs::remove_file(&tmp); // best-effort cleanup of our own temp on failure
    }
    res
}

/// Read a file that [`atomic_write`] may be replacing concurrently.
///
/// POSIX `rename` swaps a directory entry, so a reader opens either the old inode or the new one and
/// **always succeeds**. Windows offers no such guarantee: while `MoveFileEx` replaces the target the
/// name is briefly *delete-pending*, and opening a delete-pending file fails with
/// `ERROR_ACCESS_DENIED`. A plain `fs::read` therefore fails spuriously whenever a writer happens to
/// be mid-replace — and every caller here treats an unreadable config as "no config", so the user's
/// reference overrides or settings would silently revert to defaults. Retry through that window; it
/// lasts microseconds.
///
/// A genuinely **missing** file still returns `NotFound` immediately: that is the ordinary
/// no-config-yet case and must stay fast.
pub fn read_atomic(path: &Path) -> std::io::Result<Vec<u8>> {
    retry_while_replacing(|| std::fs::read(path))
}

/// Whether an error is Windows' transient "this path is being replaced right now" family:
/// `ERROR_ACCESS_DENIED` (5, delete-pending) and the sharing violations (32/33). Unix produces none
/// of these for renaming over, or opening, a file being replaced — so the retries are inert there.
fn is_replace_race(e: &std::io::Error) -> bool {
    matches!(e.raw_os_error(), Some(5) | Some(32) | Some(33)) || e.kind() == std::io::ErrorKind::PermissionDenied
}

/// Run `op`, retrying briefly while it fails with a [`is_replace_race`] error. Capped at ~0.4 s
/// total: long enough to outlast any replace window, short enough that a genuine permission problem
/// still surfaces promptly.
fn retry_while_replacing<T>(mut op: impl FnMut() -> std::io::Result<T>) -> std::io::Result<T> {
    const ATTEMPTS: u32 = 25;
    let mut delay = std::time::Duration::from_micros(200);
    for attempt in 1..=ATTEMPTS {
        match op() {
            Ok(v) => return Ok(v),
            Err(e) if attempt < ATTEMPTS && is_replace_race(&e) => {
                std::thread::sleep(delay);
                delay = (delay * 2).min(std::time::Duration::from_millis(20));
            }
            Err(e) => return Err(e),
        }
    }
    unreachable!("the final attempt returns")
}

/// Cache root: `$NAVIGATOR_REFGENOME_DIR`, else `~/.decodingus`, else the current dir.
pub fn base_dir() -> PathBuf {
    if let Some(dir) = std::env::var_os("NAVIGATOR_REFGENOME_DIR") {
        return PathBuf::from(dir);
    }
    navigator_domain::paths::decodingus_dir()
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
    base.join("liftover")
        .join(format!("{}-to-{}.chain", from.as_str(), to.as_str()))
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
    fn atomic_write_replaces_wholly_and_leaves_no_torn_file() {
        let dir = std::env::temp_dir().join(format!("atomicw_{}_{}", std::process::id(), "seq"));
        let _ = std::fs::remove_dir_all(&dir);
        let path = dir.join("cfg.json");
        // A long write then a SHORT write: the short must fully replace — no leftover tail (the exact
        // corruption a non-atomic write leaves when the new content is shorter than the old).
        atomic_write(&path, b"{\"references\":{\"a\":1,\"b\":2,\"cccccccccc\":3}}").unwrap();
        atomic_write(&path, b"{\"x\":1}").unwrap();
        assert_eq!(std::fs::read_to_string(&path).unwrap(), "{\"x\":1}");
        // The rename consumes the temp, so only the target remains — no stray `*.tmp.*` files.
        let entries: Vec<_> = std::fs::read_dir(&dir).unwrap().map(|e| e.unwrap().file_name()).collect();
        assert_eq!(entries, vec![std::ffi::OsString::from("cfg.json")]);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn atomic_write_reads_are_never_torn_under_concurrency() {
        // Many threads hammer the same path with different-length payloads. Every concurrent read
        // must observe *exactly one* whole payload — never head-of-one + tail-of-another, which is
        // what racing `fs::write` calls produce (issue #26).
        use std::sync::Arc;
        let dir = std::env::temp_dir().join(format!("atomicw_{}_conc", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let path = Arc::new(dir.join("cfg.json"));
        let payloads: Arc<Vec<String>> =
            Arc::new((0..8).map(|i| format!("[{i}{}]", ",0".repeat(i * 400))).collect());
        atomic_write(&path, payloads[0].as_bytes()).unwrap();
        let mut handles = Vec::new();
        for _ in 0..24 {
            let (path, payloads) = (path.clone(), payloads.clone());
            handles.push(std::thread::spawn(move || {
                for p in payloads.iter() {
                    atomic_write(&path, p.as_bytes()).unwrap();
                    // `read_atomic`, not `fs::read`: on Windows the target is briefly delete-pending
                    // mid-replace and a plain open fails with "Access is denied" (this test caught
                    // exactly that on windows-msvc). The invariant under test is the *content* one —
                    // whatever we read is one whole payload, never a splice of two.
                    let got = String::from_utf8(read_atomic(&path).unwrap()).unwrap();
                    assert!(payloads.contains(&got), "torn read: {got:?}");
                }
            }));
        }
        for h in handles {
            h.join().unwrap();
        }
        let _ = std::fs::remove_dir_all(&dir);
    }

    /// The retry must not turn the ordinary "no config yet" case into a stall: a missing file
    /// reports `NotFound` on the first attempt, and a present one reads straight through.
    #[test]
    fn read_atomic_is_immediate_for_missing_and_present_files() {
        let dir = std::env::temp_dir().join(format!("atomicr_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let path = dir.join("cfg.json");
        let started = std::time::Instant::now();
        let err = read_atomic(&path).expect_err("missing file must error");
        assert_eq!(err.kind(), std::io::ErrorKind::NotFound);
        assert!(started.elapsed() < std::time::Duration::from_millis(100), "missing file was retried");
        atomic_write(&path, b"{\"x\":1}").unwrap();
        assert_eq!(read_atomic(&path).unwrap(), b"{\"x\":1}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn paths_are_under_the_base() {
        let base = Path::new("/tmp/dun");
        assert_eq!(
            reference_path(base, Build::Chm13v2),
            Path::new("/tmp/dun/references/chm13v2.0.fa")
        );
        assert_eq!(
            reference_fai(base, Build::Chm13v2),
            Path::new("/tmp/dun/references/chm13v2.0.fa.fai")
        );
        assert_eq!(
            chain_path(base, Build::Grch38, Build::Chm13v2),
            Path::new("/tmp/dun/liftover/GRCh38-to-chm13v2.0.chain")
        );
    }

    #[test]
    fn base_dir_honors_env_override() {
        // Safe: this test only reads the var via base_dir; set+remove around the assertion.
        std::env::set_var("NAVIGATOR_REFGENOME_DIR", "/tmp/refcache");
        assert_eq!(base_dir(), Path::new("/tmp/refcache"));
        std::env::remove_var("NAVIGATOR_REFGENOME_DIR");
    }
}
