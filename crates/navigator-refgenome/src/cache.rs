//! On-disk cache layout. Mirrors the `~/.decodingus` convention (and the `NAVIGATOR_TREE_DIR`
//! override pattern) used elsewhere: references live under `<base>/references/`, liftover
//! chains under `<base>/liftover/`.

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use crate::registry::Build;

/// Write `bytes` to `path` **atomically**. It writes a temp file with a unique name, in the same
/// directory, flushes it to disk, and then renames it over the target.
///
/// `rename` is atomic on a POSIX filesystem. So a reader never sees a file that is torn, half
/// written, or made of a new head and an old tail. That holds while another writer runs at the same
/// time, and it holds if a writer stops in the middle.
///
/// This is the safe replacement for `std::fs::write` on a **config** file, such as
/// `reference_sources.json` or `settings.json`. Spawned tasks write those, and two tasks can write
/// at the same time. A plain `fs::write` truncates the file and then streams into it. Two such
/// writes then mix into a corrupt file: a short head from one, and a long tail from the other.
///
/// It creates the parent directory when that is absent. A temp file left by a write that stopped in
/// the middle does no harm. Nothing reads it. A later write with the same name replaces it, or it
/// stays on disk as a `*.tmp.*` that anyone can see.
pub fn atomic_write(path: &Path, bytes: &[u8]) -> std::io::Result<()> {
    use std::io::Write;
    static COUNTER: AtomicU64 = AtomicU64::new(0);
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)?;
    }
    // A unique temp name **in the same directory**, because a rename must stay on one file system.
    // The pid, plus a counter for the whole process, keeps two writers that run at the same time
    // apart. Neither one can use the other's temp file.
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

/// Read a file while [`atomic_write`] can replace it at the same time.
///
/// A POSIX `rename` swaps a directory entry. So a reader opens the old inode or the new one, and
/// **always succeeds**.
///
/// Windows gives no such guarantee. While `MoveFileEx` replaces the target, the name is
/// *delete-pending* for a moment, and an open of a delete-pending file fails with
/// `ERROR_ACCESS_DENIED`. So a plain `fs::read` fails for no good reason whenever a writer is in the
/// middle of a replace.
///
/// Every caller here reads an unreadable config as "no config". So the user's reference overrides,
/// or their settings, would go back to the defaults with no warning. Try again through that window.
/// It lasts microseconds.
///
/// A genuinely **missing** file still returns `NotFound` immediately: that is the ordinary
/// no-config-yet case and must stay fast.
pub fn read_atomic(path: &Path) -> std::io::Result<Vec<u8>> {
    retry_while_replacing(|| std::fs::read(path))
}

/// Whether an error belongs to the transient Windows family that means "another process replaces
/// this path now". That family is `ERROR_ACCESS_DENIED` (5, delete-pending) and the two
/// share-violation codes, 32 and 33. Unix gives none of these when a process renames over a file,
/// or opens a file that another process replaces. So the second try does nothing there.
fn is_replace_race(e: &std::io::Error) -> bool {
    matches!(e.raw_os_error(), Some(5) | Some(32) | Some(33)) || e.kind() == std::io::ErrorKind::PermissionDenied
}

/// Run `op`, and try again for a short time while it fails with an [`is_replace_race`] error. The
/// cap is about 0.4 s in total. That is long enough to outlast any replace window, and short enough
/// that a real permission problem still shows quickly.
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

/// Age of a cached file in days (for TTL checks); `None` if it does not exist or its mtime is
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
        // A long write, then a SHORT write. The short one must replace the file completely, and
        // leave no tail. That tail is exactly what a write which is not atomic leaves behind, when
        // the new content is shorter than the old.
        atomic_write(&path, b"{\"references\":{\"a\":1,\"b\":2,\"cccccccccc\":3}}").unwrap();
        atomic_write(&path, b"{\"x\":1}").unwrap();
        assert_eq!(std::fs::read_to_string(&path).unwrap(), "{\"x\":1}");
        // The rename consumes the temp file, so only the target stays. It leaves no `*.tmp.*`.
        let entries: Vec<_> = std::fs::read_dir(&dir)
            .unwrap()
            .map(|e| e.unwrap().file_name())
            .collect();
        assert_eq!(entries, vec![std::ffi::OsString::from("cfg.json")]);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn atomic_write_reads_are_never_torn_under_concurrency() {
        // Many threads write the same path, fast, with payloads of different lengths. Every read
        // that runs at the same time must see *exactly one* whole payload. It must never see the
        // head of one and the tail of another, which is what `fs::write` calls give when they run
        // together (issue #26).
        use std::sync::Arc;
        let dir = std::env::temp_dir().join(format!("atomicw_{}_conc", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let path = Arc::new(dir.join("cfg.json"));
        let payloads: Arc<Vec<String>> = Arc::new((0..8).map(|i| format!("[{i}{}]", ",0".repeat(i * 400))).collect());
        atomic_write(&path, payloads[0].as_bytes()).unwrap();
        let mut handles = Vec::new();
        for _ in 0..24 {
            let (path, payloads) = (path.clone(), payloads.clone());
            handles.push(std::thread::spawn(move || {
                for p in payloads.iter() {
                    atomic_write(&path, p.as_bytes()).unwrap();
                    // Use `read_atomic`, and not `fs::read`. On Windows, a replace leaves the
                    // target in a delete-pending state for a moment, and a plain open then fails
                    // with `Access is denied`. This test caught exactly that on windows-msvc.
                    //
                    // The rule under test is about the *content*: what we read is one whole
                    // payload, and never two payloads joined.
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

    /// The second try must not make the ordinary "no config yet" case slow. A file that is absent
    /// gives `NotFound` on the first try, and a file that is present reads straight through.
    #[test]
    fn read_atomic_is_immediate_for_missing_and_present_files() {
        let dir = std::env::temp_dir().join(format!("atomicr_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        let path = dir.join("cfg.json");
        let started = std::time::Instant::now();
        let err = read_atomic(&path).expect_err("missing file must error");
        assert_eq!(err.kind(), std::io::ErrorKind::NotFound);
        assert!(
            started.elapsed() < std::time::Duration::from_millis(100),
            "missing file was retried"
        );
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
        // This is safe. The test reads the var through base_dir alone, and it sets and removes the
        // var around the assertion.
        std::env::set_var("NAVIGATOR_REFGENOME_DIR", "/tmp/refcache");
        assert_eq!(base_dir(), Path::new("/tmp/refcache"));
        std::env::remove_var("NAVIGATOR_REFGENOME_DIR");
    }
}
