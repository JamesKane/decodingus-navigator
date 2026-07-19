//! Answering *why* an alignment can't be read, instead of guessing.
//!
//! The reader helpers all funnel their failures through [`AnalysisError::io`], which formats the
//! path **the caller passed**, not the file that actually failed. That is fine for a hot loop but
//! actively misleading at the edge: `open_indexed` hands it the CRAM, yet the open it performs
//! also autoloads the sibling `.crai`, resolves the reference FASTA and that FASTA's `.fai`. So an
//! unreadable index reports `io error on sample.cram`, and whoever reads that message goes looking
//! at the CRAM. Worse, [`crate::reader::has_region_index`] is built on `Path::exists`, which
//! answers `false` for *both* "no index here" and "the OS refused to tell me" — the two cases with
//! completely different fixes.
//!
//! This module takes the opposite approach: probe each participating file **separately**, name it
//! explicitly, and keep the raw `errno` rather than collapsing everything to a `bool`. The errno is
//! the whole diagnosis on macOS, where the three failures look identical in a status bar but mean
//! unrelated things:
//!
//! | errno | name | meaning | fix |
//! |---|---|---|---|
//! | 2 | `ENOENT` | the file is not there | create/fetch it |
//! | 13 | `EACCES` | Unix mode bits deny it | `chmod` / `chown` |
//! | 1 | `EPERM` | **macOS privacy (TCC) denied it** | grant Full Disk Access, or move the file |
//!
//! `EPERM` is the one that motivated this module. It is not a Unix permission failure — those are
//! `EACCES` — it is macOS refusing the process regardless of mode bits, which is why a `chmod 777`
//! file in `~/Desktop` still fails. Reading a directory listing is enough to distinguish the cases,
//! so [`diagnose`] does that too: an index that `stat` denies but that shows up in the parent
//! directory is a privacy denial, full stop.
//!
//! Nothing here mutates, downloads, or decodes more than a header and one region query, so it is
//! always safe to run — including on a file that is already failing.

use std::fmt;
use std::fs::File;
use std::path::{Path, PathBuf};

use noodles::core::Region;

use crate::reader::{self, detect_format, Format};

/// How a single check came out. `Warn` is for a condition that degrades behaviour but has a
/// working fallback (a missing index still reads sequentially); `Fail` is for one that stops the
/// operation outright.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Status {
    Ok,
    Warn,
    Fail,
}

impl Status {
    fn marker(self) -> &'static str {
        match self {
            Status::Ok => "ok  ",
            Status::Warn => "WARN",
            Status::Fail => "FAIL",
        }
    }
}

/// One named check against one named file. `path` is the file *this* check actually touched — the
/// point of the whole module — so a failure is never attributed to a file that was merely nearby.
#[derive(Debug, Clone, serde::Serialize)]
pub struct Check {
    pub name: String,
    pub path: Option<PathBuf>,
    pub status: Status,
    pub detail: String,
    /// The raw OS error number, when the check failed on a syscall. Kept unmapped because the
    /// interpretation is platform-specific and the number is what makes a bug report actionable.
    pub errno: Option<i32>,
}

impl Check {
    fn ok(name: impl Into<String>, path: Option<PathBuf>, detail: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            path,
            status: Status::Ok,
            detail: detail.into(),
            errno: None,
        }
    }

    fn new(name: impl Into<String>, path: Option<PathBuf>, status: Status, detail: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            path,
            status,
            detail: detail.into(),
            errno: None,
        }
    }
}

/// The outcome of diagnosing one alignment.
#[derive(Debug, Clone, serde::Serialize)]
pub struct Report {
    pub alignment: PathBuf,
    pub reference: Option<PathBuf>,
    pub checks: Vec<Check>,
}

impl Report {
    /// Whether any check failed outright (warnings don't count — they have fallbacks).
    pub fn failed(&self) -> bool {
        self.checks.iter().any(|c| c.status == Status::Fail)
    }

    /// The first failing check — the one whose fix unblocks the rest, since later checks depend on
    /// earlier ones succeeding.
    pub fn first_failure(&self) -> Option<&Check> {
        self.checks.iter().find(|c| c.status == Status::Fail)
    }

    fn push(&mut self, c: Check) {
        self.checks.push(c);
    }
}

impl fmt::Display for Report {
    /// A pasteable plain-text report. This is the format a user drops into a bug report, so it
    /// leads with the failing check rather than making the reader scan for it.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "alignment: {}", self.alignment.display())?;
        match &self.reference {
            Some(r) => writeln!(f, "reference: {}", r.display())?,
            None => writeln!(f, "reference: (none supplied)")?,
        }
        writeln!(f)?;
        for c in &self.checks {
            write!(f, "  [{}] {}", c.status.marker(), c.name)?;
            if let Some(p) = &c.path {
                write!(f, " — {}", p.display())?;
            }
            writeln!(f)?;
            if !c.detail.is_empty() {
                writeln!(f, "         {}", c.detail)?;
            }
        }
        if let Some(first) = self.first_failure() {
            writeln!(f)?;
            writeln!(f, "diagnosis: {}", first.name)?;
            if let Some(p) = &first.path {
                writeln!(f, "  file: {}", p.display())?;
            }
            writeln!(f, "  {}", first.detail)?;
        }
        Ok(())
    }
}

/// Explain an I/O error in terms of what the user has to *do*, keyed on the raw errno. The
/// distinction that matters is `EPERM` vs `EACCES`: they render almost identically in a status bar
/// ("Operation not permitted" vs "Permission denied") and have nothing to do with each other.
fn explain(path: &Path, e: &std::io::Error) -> (Status, String, Option<i32>) {
    let errno = e.raw_os_error();
    let detail = match errno {
        Some(2) => format!("not found: {}", path.display()),
        Some(13) => format!(
            "denied by Unix permissions ({e}). Check the mode bits and owner on this file and every \
             directory above it."
        ),
        Some(1) if cfg!(target_os = "macos") => format!(
            "macOS denied access to this file ({e}). This is the privacy layer (TCC), not file \
             permissions — mode bits are irrelevant and chmod will not help. Either grant the app \
             Full Disk Access in System Settings › Privacy & Security › Full Disk Access, or move \
             the file somewhere unprotected (not Desktop/Documents/Downloads, not iCloud Drive, and \
             not an external or network volume). Note that a grant is tied to the app's code \
             signature, so replacing or rebuilding the binary revokes it."
        ),
        _ => format!("{e}"),
    };
    (Status::Fail, detail, errno)
}

/// What a single file looks like to this process: does it exist, and can we actually open it?
///
/// Both halves are necessary. `metadata` alone answers a different question than `open` on macOS —
/// a privacy denial can let `stat` through and refuse the `open`, or refuse both — so the check
/// that matters is the one the reader will actually perform, which is opening it.
fn probe_file(name: &str, path: &Path) -> Check {
    match File::open(path) {
        Ok(_) => {
            let size = std::fs::metadata(path).map(|m| m.len()).unwrap_or(0);
            Check::ok(name, Some(path.to_path_buf()), format!("readable, {size} bytes"))
        }
        Err(e) => {
            let (status, mut detail, errno) = explain(path, &e);
            // A file the OS won't open but that is visible in its own directory listing is being
            // withheld, not absent — worth saying, because "not found" would send the user looking
            // for a file that is sitting right there.
            if errno == Some(2) && directory_lists(path) {
                detail = format!(
                    "{detail}\n         (the parent directory lists this name, so it exists but \
                     cannot be opened)"
                );
            }
            Check {
                name: name.to_string(),
                path: Some(path.to_path_buf()),
                status,
                detail,
                errno,
            }
        }
    }
}

/// Whether `path`'s own parent directory lists it. Distinguishes "absent" from "withheld".
fn directory_lists(path: &Path) -> bool {
    let (Some(dir), Some(file)) = (path.parent(), path.file_name()) else {
        return false;
    };
    let Ok(entries) = std::fs::read_dir(dir) else {
        return false;
    };
    entries.flatten().any(|e| e.file_name() == file)
}

/// Every path that could serve as the coordinate index for `path`, in the order the readers accept
/// them: the dotted `foo.cram.crai` spelling `samtools` writes, then the replaced `foo.crai`.
pub fn index_candidates(path: &Path) -> Vec<PathBuf> {
    match detect_format(path) {
        Format::Bam => vec![path.with_extension("bam.bai"), path.with_extension("bai")],
        Format::Cram => vec![path.with_extension("cram.crai"), path.with_extension("crai")],
    }
}

/// Diagnose an alignment, in dependency order: the file itself, then its index, then the reference
/// and the reference's own index, then the operations built on all of them (header read, indexed
/// open, one region query).
///
/// The order is the point — each check presupposes the previous one, so [`Report::first_failure`]
/// names the thing to fix rather than the last thing to fall over. Reads at most one header and one
/// region's records; never writes, never downloads.
pub fn diagnose(alignment: &Path, reference: Option<&Path>) -> Report {
    let mut report = Report {
        alignment: alignment.to_path_buf(),
        reference: reference.map(Path::to_path_buf),
        checks: Vec::new(),
    };
    let format = detect_format(alignment);
    report.push(Check::ok(
        "format",
        None,
        match format {
            Format::Bam => "BAM (detected from the extension)",
            Format::Cram => "CRAM (detected from the extension) — a reference FASTA is required",
        },
    ));

    let file = probe_file("alignment file", alignment);
    let alignment_ok = file.status == Status::Ok;
    report.push(file);
    if !alignment_ok {
        return report;
    }

    // The index. Its absence is a warning, not a failure: sequential walks (read metrics, coverage,
    // sex) fall back and succeed, which is exactly why an alignment can look healthy in the UI
    // right up until something needs a region query. An index that exists but won't open is a
    // failure, and is the case `has_region_index` silently reports as "no index".
    let candidates = index_candidates(alignment);
    let found = candidates.iter().find(|p| directory_lists(p));
    let has_index = found.is_some();
    match found {
        None => report.push(Check::new(
            "coordinate index",
            None,
            Status::Warn,
            format!(
                "no index found. Looked for: {}. Sequential passes (read metrics, coverage, sex) \
                 still work; anything needing a region query — Y haplogroup, mtDNA, SV, callable \
                 intervals — does not. Build one with `samtools index {}`.",
                candidates
                    .iter()
                    .map(|p| p.display().to_string())
                    .collect::<Vec<_>>()
                    .join(", "),
                alignment.display()
            ),
        )),
        Some(idx) => {
            let probe = probe_file("coordinate index", idx);
            let index_ok = probe.status == Status::Ok;
            report.push(probe);
            if !index_ok {
                // No point attempting the indexed open — it would fail and, worse, blame the CRAM.
                return report;
            }
        }
    }

    // The reference. Required to decode CRAM at all; for BAM it is optional here.
    let reference = match (format, reference) {
        (Format::Cram, None) => {
            report.push(Check::new(
                "reference FASTA",
                None,
                Status::Fail,
                "a CRAM cannot be decoded without its reference FASTA, and none was supplied or \
                 resolved for this alignment's build.",
            ));
            return report;
        }
        (_, r) => r,
    };
    if let Some(r) = reference {
        let probe = probe_file("reference FASTA", r);
        let reference_ok = probe.status == Status::Ok;
        report.push(probe);
        if !reference_ok {
            return report;
        }
        // CRAM decode goes through an *indexed* FASTA reader, so a missing `.fai` fails the open
        // just as hard as a missing FASTA — and reports the FASTA's path when it does.
        let fai = PathBuf::from(format!("{}.fai", r.display()));
        let probe = probe_file("reference index (.fai)", &fai);
        let fai_ok = probe.status == Status::Ok;
        report.push(probe);
        if !fai_ok {
            return report;
        }
    }

    // Now the composite operations, in the order the analysis paths perform them.
    match reader::read_header(alignment, reference) {
        Ok(h) => report.push(Check::ok(
            "read header",
            Some(alignment.to_path_buf()),
            format!("{} reference sequences", h.reference_sequences().len()),
        )),
        Err(e) => {
            report.push(Check::new(
                "read header",
                Some(alignment.to_path_buf()),
                Status::Fail,
                e.to_string(),
            ));
            return report;
        }
    }

    let (header, mut idx) = match reader::open_indexed(alignment, reference) {
        Ok(v) => {
            report.push(Check::ok(
                "open indexed",
                Some(alignment.to_path_buf()),
                "the index loaded and the file is ready for region queries",
            ));
            v
        }
        Err(e) => {
            // The whole point of this module: don't repeat the upstream message's mistake of
            // blaming the alignment. If we already established there is no index, *that* is the
            // finding — an `ENOENT` naming the CRAM here means the reader could not autoload a
            // sibling index, not that the CRAM went missing between two reads of it.
            let check = if has_index {
                Check::new(
                    "open indexed",
                    Some(alignment.to_path_buf()),
                    Status::Fail,
                    format!(
                        "{e}\n         (this message names the alignment, but the file itself \
                         opened fine above — the failure is in its index or the reference)"
                    ),
                )
            } else {
                Check::new(
                    "open indexed",
                    Some(alignment.to_path_buf()),
                    Status::Fail,
                    format!(
                        "there is no coordinate index, so region queries cannot run — this is the \
                         missing `.crai`/`.bai` reported above, not a problem with the alignment \
                         itself. Build one with `samtools index {}`.\n         (underlying: {e})",
                        alignment.display()
                    ),
                )
            };
            report.push(check);
            return report;
        }
    };

    // One real region query. Everything above can pass on a file whose index is stale or truncated;
    // this is the check that actually exercises a seek, which is what the Y/mtDNA/SV paths do.
    let Some(contig) = header
        .reference_sequences()
        .keys()
        .next()
        .map(|k| String::from_utf8_lossy(k.as_ref()).into_owned())
    else {
        report.push(Check::new(
            "region query",
            None,
            Status::Fail,
            "the header declares no reference sequences",
        ));
        return report;
    };
    let region = Region::new(contig.as_bytes().to_vec(), ..);
    let probe = match idx.query(&header, &region) {
        Ok(mut records) => match records.next() {
            Some(Err(e)) => Check::new(
                "region query",
                Some(alignment.to_path_buf()),
                Status::Fail,
                format!("decoding the first record of {contig} failed: {e}"),
            ),
            Some(Ok(_)) => Check::ok(
                "region query",
                None,
                format!("seeked to {contig} and decoded a record"),
            ),
            None => Check::new(
                "region query",
                None,
                Status::Warn,
                format!("seeked to {contig} but it holds no records"),
            ),
        },
        Err(e) => Check::new(
            "region query",
            Some(alignment.to_path_buf()),
            Status::Fail,
            format!("querying {contig} failed: {e}"),
        ),
    };
    report.push(probe);
    report
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn index_candidates_cover_both_spellings() {
        let cram = index_candidates(Path::new("/d/s.hg38.sorted.cram"));
        assert_eq!(cram[0], PathBuf::from("/d/s.hg38.sorted.cram.crai"));
        assert_eq!(cram[1], PathBuf::from("/d/s.hg38.sorted.crai"));

        let bam = index_candidates(Path::new("/d/s.bam"));
        assert_eq!(bam[0], PathBuf::from("/d/s.bam.bai"));
        assert_eq!(bam[1], PathBuf::from("/d/s.bai"));
    }

    /// EPERM and EACCES are the two that a status bar makes look alike and that have unrelated
    /// fixes, so the explanation must separate them — chmod advice on a TCC denial sends the user
    /// down a dead end.
    #[test]
    fn explains_tcc_denial_separately_from_unix_permissions() {
        let p = Path::new("/d/s.cram");

        let eacces = std::io::Error::from_raw_os_error(13);
        let (_, detail, errno) = explain(p, &eacces);
        assert_eq!(errno, Some(13));
        assert!(detail.contains("Unix permissions"), "{detail}");

        let eperm = std::io::Error::from_raw_os_error(1);
        let (_, detail, errno) = explain(p, &eperm);
        assert_eq!(errno, Some(1));
        if cfg!(target_os = "macos") {
            assert!(detail.contains("Full Disk Access"), "{detail}");
            assert!(detail.contains("chmod will not help"), "{detail}");
        }
    }

    /// A readable file with no index must warn about the *index* and name both accepted spellings,
    /// and must not let that warning masquerade as a problem with the alignment — the confusion
    /// that made the original bug report unreadable.
    #[test]
    fn missing_index_is_reported_against_the_index_not_the_alignment() {
        let dir = std::env::temp_dir().join("navigator-preflight-noindex");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("sample.bam");
        std::fs::write(&bam, b"not really a bam").unwrap();

        let report = diagnose(&bam, None);
        let index = report
            .checks
            .iter()
            .find(|c| c.name == "coordinate index")
            .expect("index is always checked");
        assert_eq!(index.status, Status::Warn, "a missing index has a sequential fallback");
        assert!(index.detail.contains("sample.bam.bai"), "{}", index.detail);
        assert!(index.detail.contains("sample.bai"), "{}", index.detail);

        // The alignment itself opened fine, so nothing may blame it for the missing index.
        let file = report
            .checks
            .iter()
            .find(|c| c.name == "alignment file")
            .expect("file is always checked");
        assert_eq!(file.status, Status::Ok, "{report}");

        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn missing_alignment_fails_fast_and_names_itself() {
        let report = diagnose(Path::new("/nonexistent/sample.cram"), None);
        let first = report.first_failure().expect("missing file must fail");
        assert_eq!(first.name, "alignment file");
        assert_eq!(first.errno, Some(2));
        // The reference check must not run — the report stops at the first real blocker.
        assert!(!report.checks.iter().any(|c| c.name == "reference FASTA"), "{report}");
    }
}
