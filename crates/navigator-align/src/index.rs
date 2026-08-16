//! The minimap2 index (`.mmi`) cache, and building one part by part.
//!
//! An index is specific to both the reference build *and* the preset — `sr`, `map-hifi`, and
//! `map-ont` disagree on k-mer and window size, so an index built for one is wrong for another.
//! The cache key is therefore `(build, preset)`, laid out alongside the reference cache the
//! refgenome crate already owns:
//!
//! ```text
//! <base>/minimap2_index/<build>/<preset>.mmi
//! ```
//!
//! ## Part by part
//!
//! [`build_index`] streams: read one index part from the FASTA, write it, drop it, repeat. That is
//! what keeps peak memory at the [`BatchSize`] rather than at the whole genome — building a single
//! resident index for CHM13 costs ~19 GiB, and building it in 1 Gbase parts costs 11.7 GiB for the
//! same output in the same wall time.
//!
//! The `.mmi` is written to a temporary path and renamed on success, so an interrupted build can
//! never leave a half-written index that a later run would load as if it were whole. The file is
//! 8.93 GB for CHM13 — large enough that "just rebuild it if it looks wrong" is not a strategy.

use std::path::{Path, PathBuf};

use minimap2::flags::IdxFlags;
use minimap2::index::reader::IdxReader;
use minimap2::index::MmIdx;
use minimap2::options::IdxOpt;

use crate::batch::BatchSize;
use crate::error::AlignError;
use crate::preset::Preset;

/// Progress during a long index build: `(parts_done, bases_done)`. Parts arrive as they are
/// written, so a caller can report "part 2 of ~4" against [`BatchSize::part_estimate`].
pub type ProgressFn<'a> = &'a mut dyn FnMut(usize, u64);

/// The cache root the aligner index lives under: `$NAVIGATOR_REFGENOME_DIR`, else `~/.decodingus`.
///
/// Deliberately the same answer `navigator-refgenome::cache::base_dir` gives, reached the same way
/// — through `navigator_domain::paths::decodingus_dir`, the one definition of the cache root — so
/// `minimap2_index/` lands beside `references/` and `liftover/` rather than in a second location
/// that only this crate knows about. This crate is a leaf and cannot depend on `navigator-refgenome`
/// (that would invert the layering), which is why the resolution is repeated rather than imported;
/// the shared *definition* is what stops the two drifting.
pub fn cache_root() -> PathBuf {
    if let Some(dir) = std::env::var_os("NAVIGATOR_REFGENOME_DIR") {
        return PathBuf::from(dir);
    }
    navigator_domain::paths::decodingus_dir()
}

/// Where the cached index for `(build, preset)` lives under `base`.
///
/// `base` is the refgenome cache root — [`cache_root`] resolves it, or a caller that already has
/// one (the app, which resolves it once) passes it in. Tests point it anywhere.
pub fn index_path(base: &Path, build: &str, preset: Preset) -> PathBuf {
    base.join("minimap2_index")
        .join(build)
        .join(format!("{}.mmi", preset.as_str()))
}

/// Build the index for `reference` into the cache, unless it is already there.
///
/// Returns the cached path. Idempotent: an existing index is returned untouched, which is what
/// makes this safe to call at the top of every realignment job.
pub fn ensure_index(
    base: &Path,
    build: &str,
    reference: &Path,
    preset: Preset,
    batch: BatchSize,
    progress: ProgressFn<'_>,
) -> Result<PathBuf, AlignError> {
    let path = index_path(base, build, preset);
    if path.is_file() {
        return Ok(path);
    }
    build_index(reference, &path, preset, batch, progress)?;
    Ok(path)
}

/// [`ensure_index`] against the real cache root, sizing the index for this machine.
///
/// This is the call a job should make: it resolves where the cache lives, picks a batch size from
/// the machine's RAM, and returns a ready index — none of which the caller should have to know how
/// to do. See [`BatchSize::for_this_machine`] for why the sizing is detected rather than asked.
pub fn ensure_cached_index(
    build: &str,
    reference: &Path,
    preset: Preset,
    progress: ProgressFn<'_>,
) -> Result<PathBuf, AlignError> {
    ensure_index(
        &cache_root(),
        build,
        reference,
        preset,
        BatchSize::for_this_machine(),
        progress,
    )
}

/// Build a `.mmi` for `reference` at `out`, one part at a time.
///
/// Exposed separately from [`ensure_index`] so a caller can build to an arbitrary location — the
/// tests do, and so would a "rebuild this index" maintenance action.
pub fn build_index(
    reference: &Path,
    out: &Path,
    preset: Preset,
    batch: BatchSize,
    progress: ProgressFn<'_>,
) -> Result<PathBuf, AlignError> {
    let parent = out.parent().unwrap_or(Path::new("."));
    std::fs::create_dir_all(parent).map_err(|e| AlignError::io(parent, e))?;

    let opts = idx_opts(preset)?;
    let reference_str = path_str(reference)?;

    // `is_idx = false`: the input is a FASTA to sketch, not an existing `.mmi` to load.
    let mut reader = IdxReader::open(
        &reference_str,
        false,
        opts.w as i32,
        opts.k as i32,
        opts.bucket_bits,
        IdxFlags::empty(),
        opts.mini_batch_size,
        batch.bases(),
    )
    .map_err(|e| AlignError::io(reference, e))?;

    // Write to a sibling temp path and rename at the end: a torn 8.93 GB index that looks complete
    // is a far worse outcome than a build that has to be repeated.
    let tmp = out.with_extension("mmi.partial");
    let file = std::fs::File::create(&tmp).map_err(|e| AlignError::io(&tmp, e))?;
    // Paced, like every other multi-GB write in the pipeline: an index build is a one-off, but it
    // is nine gigabytes in one uninterrupted push, and it happens on the machine of a user who is
    // still using it. It also puts those bytes in the counter the resource watch reports, so the
    // stage stops looking idle in the log.
    let mut writer = std::io::BufWriter::with_capacity(1 << 20, navigator_resource::PacedFile::new(file));

    let mut parts = 0usize;
    let mut bases = 0u64;
    loop {
        let Some(part) = reader.read_next().map_err(|e| AlignError::io(reference, e))? else {
            break;
        };
        bases += part_bases(&part);
        parts += 1;
        write_part(&mut writer, &part, &tmp)?;
        // `part` is dropped here — this is the line that bounds peak memory to one part.
        progress(parts, bases);
    }

    use std::io::Write as _;
    writer.flush().map_err(|e| AlignError::io(&tmp, e))?;
    // Sync before the rename. The rename is what publishes this as a complete index, and a cache
    // entry whose contents are still only a page-cache promise is the torn-index case the temp path
    // exists to prevent.
    writer.get_ref().sync().map_err(|e| AlignError::io(&tmp, e))?;
    drop(writer);

    if parts == 0 {
        let _ = std::fs::remove_file(&tmp);
        return Err(AlignError::Message(format!(
            "{} produced no index parts — is it a FASTA?",
            reference.display()
        )));
    }

    std::fs::rename(&tmp, out).map_err(|e| AlignError::io(out, e))?;
    Ok(out.to_path_buf())
}

/// Total reference bases in one index part.
fn part_bases(part: &MmIdx) -> u64 {
    part.seqs.iter().map(|s| s.len as u64).sum()
}

fn write_part<W: std::io::Write>(writer: &mut W, part: &MmIdx, path: &Path) -> Result<(), AlignError> {
    minimap2::index::io::idx_dump(writer, part).map_err(|e| AlignError::io(path, e))
}

/// The indexing options a preset implies (k, w, bucket bits, mini-batch).
fn idx_opts(preset: Preset) -> Result<IdxOpt, AlignError> {
    let (io, _mo) = minimap2::prelude::preset(preset.as_str())
        .map_err(|e| AlignError::Message(format!("preset {}: {e}", preset.as_str())))?;
    Ok(io)
}

/// minimap2's API takes paths as `&str`. A non-UTF-8 path is a real possibility on both Unix and
/// Windows, so it is refused with a clear message rather than lossily converted into a path that
/// does not exist.
pub(crate) fn path_str(path: &Path) -> Result<String, AlignError> {
    path.to_str()
        .map(str::to_string)
        .ok_or_else(|| AlignError::Message(format!("path is not valid UTF-8: {}", path.display())))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn scratch(tag: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("dun-align-{}-{tag}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).unwrap();
        dir
    }

    /// A small FASTA with `contigs` sequences of `len` bases each, deterministic so a test can
    /// assert on part counts.
    fn write_fasta(dir: &Path, contigs: usize, len: usize) -> PathBuf {
        let path = dir.join("ref.fa");
        let mut text = String::new();
        for c in 0..contigs {
            text.push_str(&format!(">contig{c}\n"));
            // Non-repetitive enough to produce minimizers rather than one degenerate bucket.
            let seq: String = (0..len)
                .map(|i| match (i * 7 + c * 13) % 4 {
                    0 => 'A',
                    1 => 'C',
                    2 => 'G',
                    _ => 'T',
                })
                .collect();
            for chunk in seq.as_bytes().chunks(60) {
                text.push_str(std::str::from_utf8(chunk).unwrap());
                text.push('\n');
            }
        }
        std::fs::write(&path, text).unwrap();
        path
    }

    #[test]
    fn the_cache_path_keys_on_both_build_and_preset() {
        let base = Path::new("/cache");
        let sr = index_path(base, "chm13v2.0", Preset::ShortRead);
        let ont = index_path(base, "chm13v2.0", Preset::MapOnt);
        let other_build = index_path(base, "GRCh38", Preset::ShortRead);

        assert_eq!(sr, Path::new("/cache/minimap2_index/chm13v2.0/sr.mmi"));
        assert_ne!(sr, ont, "presets disagree on k/w, so they need separate indexes");
        assert_ne!(sr, other_build);
    }

    #[test]
    fn building_a_small_reference_produces_a_loadable_index() {
        let dir = scratch("build");
        let fasta = write_fasta(&dir, 2, 5_000);
        let out = dir.join("out.mmi");

        let mut parts_seen = 0;
        build_index(
            &fasta,
            &out,
            Preset::ShortRead,
            BatchSize::default(),
            &mut |parts, _bases| parts_seen = parts,
        )
        .unwrap();

        assert!(out.is_file(), "index written");
        assert!(std::fs::metadata(&out).unwrap().len() > 0);
        assert_eq!(parts_seen, 1, "10 kb fits in one part");
        assert!(
            !out.with_extension("mmi.partial").exists(),
            "the temp file is renamed away, never left behind"
        );
    }

    /// The memory bound in action: a reference larger than the batch must yield several parts, and
    /// the resulting `.mmi` must still be one loadable file with every base accounted for.
    ///
    /// Sizing this fixture is fussier than it looks, and the shape is worth recording. A part
    /// accumulates whole sequences until the total *exceeds* the batch, so parts overshoot by up
    /// to one sequence and a reference only a little larger than the batch still comes out as one
    /// part. The fixture must also clear [`BatchSize`]'s 1 Mbase floor, which exists because
    /// smaller parts are pathological in production. 3 Mbase against a 1 Mbase batch clears both.
    #[test]
    fn a_reference_larger_than_the_batch_splits_into_several_parts() {
        let dir = scratch("split");
        let fasta = write_fasta(&dir, 6, 500_000);
        let out = dir.join("split.mmi");

        let mut parts_seen = 0;
        let mut bases_seen = 0;
        build_index(
            &fasta,
            &out,
            Preset::ShortRead,
            BatchSize::new(1_000_000),
            &mut |parts, bases| {
                parts_seen = parts;
                bases_seen = bases;
            },
        )
        .unwrap();

        assert!(parts_seen > 1, "expected a split index, got {parts_seen} part(s)");
        assert_eq!(bases_seen, 3_000_000, "every base accounted for across parts");
        assert!(out.is_file());
    }

    /// `ensure_index` is called at the top of every job, so a second call must not rebuild an
    /// 8.93 GB artifact.
    #[test]
    fn ensure_index_is_idempotent() {
        let dir = scratch("ensure");
        let fasta = write_fasta(&dir, 1, 4_000);
        let base = dir.join("cache");

        let first = ensure_index(
            &base,
            "testbuild",
            &fasta,
            Preset::ShortRead,
            BatchSize::default(),
            &mut |_, _| {},
        )
        .unwrap();
        let stamp = std::fs::metadata(&first).unwrap().modified().unwrap();

        let mut rebuilt = false;
        let second = ensure_index(
            &base,
            "testbuild",
            &fasta,
            Preset::ShortRead,
            BatchSize::default(),
            &mut |_, _| rebuilt = true,
        )
        .unwrap();

        assert_eq!(first, second);
        assert!(!rebuilt, "a cached index must not be rebuilt");
        assert_eq!(stamp, std::fs::metadata(&second).unwrap().modified().unwrap());
    }

    /// The aligner index has to land beside the reference cache, not in a second place only this
    /// crate knows about. `navigator-refgenome` resolves its root the same way — env override
    /// first, then the shared `decodingus_dir` — and this pins that agreement.
    #[test]
    fn the_cache_root_follows_the_refgenome_override() {
        let _guard = ENV_LOCK.lock().unwrap_or_else(|e| e.into_inner());
        // SAFETY: guarded by ENV_LOCK, so no other test reads the environment concurrently.
        unsafe { std::env::set_var("NAVIGATOR_REFGENOME_DIR", "/tmp/dun-cache-probe") };
        let overridden = cache_root();
        unsafe { std::env::remove_var("NAVIGATOR_REFGENOME_DIR") };
        let default = cache_root();

        assert_eq!(overridden, Path::new("/tmp/dun-cache-probe"));
        assert_eq!(
            default,
            navigator_domain::paths::decodingus_dir(),
            "without the override it must be the shared cache root, not an invention"
        );
        assert!(
            index_path(&default, "chm13v2.0", Preset::ShortRead).starts_with(&default),
            "the index lives under the cache root"
        );
    }

    /// `set_var` mutates process-global state, so the tests that touch it must not overlap.
    static ENV_LOCK: std::sync::Mutex<()> = std::sync::Mutex::new(());

    /// Feeding something that is not a reference must fail loudly rather than caching an empty
    /// index that every later job would load and quietly map nothing against.
    #[test]
    fn a_reference_with_no_sequences_is_an_error() {
        let dir = scratch("empty");
        let fasta = dir.join("empty.fa");
        std::fs::write(&fasta, "").unwrap();
        let out = dir.join("empty.mmi");

        let err = build_index(&fasta, &out, Preset::ShortRead, BatchSize::default(), &mut |_, _| {});
        assert!(err.is_err());
        assert!(!out.exists(), "nothing cached from a failed build");
    }
}
