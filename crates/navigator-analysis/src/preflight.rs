//! This module says *why* the code can not read an alignment. It does not guess.
//!
//! Every reader helper sends its failure through [`AnalysisError::io`], which prints the path
//! **that the caller gave**, and not the file that failed. That is correct in a hot loop, but at
//! the edge it points the user the wrong way. `open_indexed` gives it the CRAM. But the open that
//! it does also loads the `.crai` beside that CRAM, resolves the reference FASTA, and resolves
//! the `.fai` of that FASTA. So an index that the code can not read reports
//! `io error on sample.cram`, and the person who reads that message looks at the CRAM.
//!
//! [`crate::reader::has_region_index`] is worse. It stands on `Path::exists`, which answers
//! `false` for *both* "there is no index here" and "the OS refused to tell me". Those two cases
//! need completely different fixes.
//!
//! This module does the opposite. It probes each file that takes part **on its own**, it names
//! that file, and it keeps the raw `errno`. It does not collapse everything into a `bool`. On
//! macOS the errno is the whole diagnosis. The three failures look the same in a status bar, and
//! they have nothing to do with each other:
//!
//! | errno | name | what it says | the fix |
//! |---|---|---|---|
//! | 2 | `ENOENT` | the file is not there | create it, or fetch it |
//! | 13 | `EACCES` | the Unix mode bits deny it | `chmod` or `chown` |
//! | 1 | `EPERM` | **macOS privacy (TCC) denied it** | grant Full Disk Access |
//!
//! `EPERM` is the reason that this module exists. It is not a Unix permission failure, because
//! those give `EACCES`. It is macOS that refuses the process, whatever the mode bits say. That is
//! why a file at `chmod 777` in `~/Desktop` still fails.
//!
//! A read of a directory listing is enough to separate the cases, and [`diagnose`] does that too.
//! Take an index that `stat` refuses, and that the listing of the parent directory shows. That is
//! a privacy denial, and nothing else.
//!
//! Nothing here changes a file, downloads a file, or decodes more than a header and one region
//! query. It is always safe to run, and that includes a file that already fails.

use std::fmt;
use std::fs::File;
use std::path::{Path, PathBuf};

use noodles::core::Region;

use crate::reader::{self, detect_format, Format};

/// The result of one check. `Warn` marks a condition that makes the behaviour worse, but that has
/// a fallback which works. A missing index is one: the code still reads the file sequentially.
/// `Fail` marks a condition that stops the operation completely.
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

/// Which check this is. A caller branches on the identity, and not on the string that the UI
/// shows. A batch that decides whether to skip a sample must not depend on prose. Somebody can
/// change that prose, or translate it.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum CheckId {
    Format,
    AlignmentFile,
    CoordinateIndex,
    ReferenceFasta,
    ReferenceIndex,
    ReadHeader,
    OpenIndexed,
    RegionQuery,
}

impl CheckId {
    /// The human label. It is the single source of truth, so the name of a check and its
    /// identity always agree.
    pub fn label(self) -> &'static str {
        match self {
            CheckId::Format => "format",
            CheckId::AlignmentFile => "alignment file",
            CheckId::CoordinateIndex => "coordinate index",
            CheckId::ReferenceFasta => "reference FASTA",
            CheckId::ReferenceIndex => "reference index (.fai)",
            CheckId::ReadHeader => "read header",
            CheckId::OpenIndexed => "open indexed",
            CheckId::RegionQuery => "region query",
        }
    }

    /// True when a failure of this check makes the file unreadable *completely*, and that includes
    /// a sequential pass.
    ///
    /// This is what decides whether a caller may skip a sample. A broken index, and anything that
    /// stands on it, blocks a region query and nothing else. Read metrics, coverage and sex fall
    /// back to a sequential walk, and they still succeed. To read that as "nobody can analyze this
    /// sample" would throw away results that do work.
    ///
    /// The list is narrow, and that is deliberate. Only the open of the file, and the read of its
    /// header, count. Those are the minimum that every sequential path does.
    ///
    /// A problem with the reference is *not* on the list, although some steps need one. How much
    /// it matters depends on the format and on the step. And a read of the header of a CRAM
    /// already needs the reference, so a reference that is truly unusable fails
    /// [`CheckId::ReadHeader`] in any case.
    ///
    /// To skip a sample throws away work that could have succeeded. So it must follow only from a
    /// failure that leaves nothing to try.
    pub fn blocks_sequential_reads(self) -> bool {
        matches!(self, CheckId::AlignmentFile | CheckId::ReadHeader)
    }
}

/// One named check against one named file. `path` is the file that *this* check touched, which is
/// the point of the whole module. A failure never goes to a file that only sat nearby.
#[derive(Debug, Clone, serde::Serialize)]
pub struct Check {
    pub id: CheckId,
    pub name: String,
    pub path: Option<PathBuf>,
    pub status: Status,
    pub detail: String,
    /// The raw OS error number, when the check failed on a syscall. Kept unmapped because the
    /// interpretation is platform-specific and the number is what makes a bug report actionable.
    pub errno: Option<i32>,
}

impl Check {
    fn ok(id: CheckId, path: Option<PathBuf>, detail: impl Into<String>) -> Self {
        Self::new(id, path, Status::Ok, detail)
    }

    fn new(id: CheckId, path: Option<PathBuf>, status: Status, detail: impl Into<String>) -> Self {
        Self {
            id,
            name: id.label().to_string(),
            path,
            status,
            detail: detail.into(),
            errno: None,
        }
    }
}

/// The result of a diagnosis of one alignment.
#[derive(Debug, Clone, serde::Serialize)]
pub struct Report {
    pub alignment: PathBuf,
    pub reference: Option<PathBuf>,
    pub checks: Vec<Check>,
}

impl Report {
    /// True when a check failed completely. A warning does not count, because a warning has a
    /// fallback.
    pub fn failed(&self) -> bool {
        self.checks.iter().any(|c| c.status == Status::Fail)
    }

    /// The first check that failed. A fix to that one clears the rest, because a later check needs
    /// the earlier checks to pass.
    pub fn first_failure(&self) -> Option<&Check> {
        self.checks.iter().find(|c| c.status == Status::Fail)
    }

    /// True when nothing can read the file, and that includes a sequential pass.
    ///
    /// A batch must answer this question before it skips a sample. Take a failure that blocks a
    /// region query alone, such as an index that is missing, or one that the code can not read.
    /// That failure must not skip the sample. Read metrics, coverage and sex still finish through
    /// the sequential fallback. To throw those away because the Y step can not run would lose
    /// results that the user would otherwise get.
    pub fn blocks_sequential_reads(&self) -> bool {
        self.checks
            .iter()
            .any(|c| c.status == Status::Fail && c.id.blocks_sequential_reads())
    }

    fn push(&mut self, c: Check) {
        self.checks.push(c);
    }
}

impl fmt::Display for Report {
    /// A plain-text report that a user can paste. This is the form that goes into a bug report,
    /// so the check that failed comes first. The reader does not have to look for it.
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

/// Explain an I/O error by what the user must *do*. The key is the raw errno. The difference that
/// matters is `EPERM` against `EACCES`. A status bar shows the two almost the same way,
/// "Operation not permitted" against "Permission denied", and they have nothing to do with each
/// other.
fn explain(path: &Path, e: &std::io::Error) -> (Status, String, Option<i32>) {
    let errno = e.raw_os_error();
    // The key of "not found" is the portable `ErrorKind`, and not a raw errno. Unix returns
    // ENOENT (2). Windows returns ERROR_FILE_NOT_FOUND (2) *or* ERROR_PATH_NOT_FOUND (3), and
    // which one depends on the component of the path that is absent. Both of those map to
    // `NotFound`. The other branches keep the errno as their key, because they draw a difference,
    // EPERM against EACCES, that `ErrorKind` collapses.
    if e.kind() == std::io::ErrorKind::NotFound {
        return (Status::Fail, format!("not found: {}", path.display()), errno);
    }
    let detail = match errno {
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

/// What one file looks like to this process. Does it exist, and can the process open it?
///
/// Both halves are necessary. On macOS, `metadata` answers a different question from `open`. A
/// privacy denial can let `stat` through and refuse the `open`, or it can refuse both. So the
/// check that matters is the one that the reader itself does, which is the open.
fn probe_file(id: CheckId, path: &Path) -> Check {
    match File::open(path) {
        Ok(_) => {
            let size = std::fs::metadata(path).map(|m| m.len()).unwrap_or(0);
            Check::ok(id, Some(path.to_path_buf()), format!("readable, {size} bytes"))
        }
        Err(e) => {
            let (status, mut detail, errno) = explain(path, &e);
            // Take a file that the OS will not open, and that its own directory listing shows.
            // The OS holds that file back, and the file is not absent. Say so, because "not
            // found" would send the user to look for a file that is right there.
            if e.kind() == std::io::ErrorKind::NotFound && directory_lists(path) {
                detail = format!(
                    "{detail}\n         (the parent directory lists this name, so it exists but \
                     cannot be opened)"
                );
            }
            Check {
                id,
                name: id.label().to_string(),
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

/// Every path that can be the coordinate index of `path`, in the order that the readers accept
/// them. First comes the `foo.cram.crai` form that `samtools` writes, and then `foo.crai`, where
/// the extension replaces the old one.
pub fn index_candidates(path: &Path) -> Vec<PathBuf> {
    match detect_format(path) {
        Format::Bam => vec![path.with_extension("bam.bai"), path.with_extension("bai")],
        Format::Cram => vec![path.with_extension("cram.crai"), path.with_extension("crai")],
    }
}

/// Diagnose an alignment, in the order of the dependencies. First the file itself, then its index,
/// then the reference and the index of that reference. Last come the operations that stand on all
/// of those: the read of the header, the indexed open, and one region query.
///
/// The order is the point. Each check needs the check before it, so [`Report::first_failure`]
/// names the thing to fix, and not the last thing to fall over. This function reads one header at
/// most, and the records of one region. It never writes, and it never downloads.
pub fn diagnose(alignment: &Path, reference: Option<&Path>) -> Report {
    let mut report = Report {
        alignment: alignment.to_path_buf(),
        reference: reference.map(Path::to_path_buf),
        checks: Vec::new(),
    };
    let format = detect_format(alignment);
    report.push(Check::ok(
        CheckId::Format,
        None,
        match format {
            Format::Bam => "BAM (detected from the extension)",
            Format::Cram => "CRAM (detected from the extension) — a reference FASTA is required",
        },
    ));

    let file = probe_file(CheckId::AlignmentFile, alignment);
    let alignment_ok = file.status == Status::Ok;
    report.push(file);
    if !alignment_ok {
        return report;
    }

    // The index. Its absence is a warning and not a failure. The three sequential walks, which
    // are the read metrics, the coverage and the sex, fall back and succeed. That is exactly why an alignment
    // can look healthy in the UI until something needs a region query.
    //
    // An index that exists but that will not open is a failure. That is the case where
    // `has_region_index` reports "no index", and nobody sees the real cause.
    let candidates = index_candidates(alignment);
    let found = candidates.iter().find(|p| directory_lists(p));
    let has_index = found.is_some();
    match found {
        None => report.push(Check::new(
            CheckId::CoordinateIndex,
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
            let probe = probe_file(CheckId::CoordinateIndex, idx);
            let index_ok = probe.status == Status::Ok;
            report.push(probe);
            if !index_ok {
                // Do not try the indexed open. It would fail, and it would put the blame on the
                // CRAM.
                return report;
            }
        }
    }

    // The reference. Required to decode CRAM at all; for BAM it is optional here.
    let reference = match (format, reference) {
        (Format::Cram, None) => {
            report.push(Check::new(
                CheckId::ReferenceFasta,
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
        let probe = probe_file(CheckId::ReferenceFasta, r);
        let reference_ok = probe.status == Status::Ok;
        report.push(probe);
        if !reference_ok {
            return report;
        }
        // A CRAM decode goes through an FASTA reader that uses an *index*. So a missing `.fai`
        // fails the open as hard as a missing FASTA does, and it reports the path of the FASTA
        // when it fails.
        let fai = PathBuf::from(format!("{}.fai", r.display()));
        let probe = probe_file(CheckId::ReferenceIndex, &fai);
        let fai_ok = probe.status == Status::Ok;
        report.push(probe);
        if !fai_ok {
            return report;
        }
    }

    // Now the operations that put those parts together, in the order that the analysis paths do
    // them.
    match reader::read_header(alignment, reference) {
        Ok(h) => report.push(Check::ok(
            CheckId::ReadHeader,
            Some(alignment.to_path_buf()),
            format!("{} reference sequences", h.reference_sequences().len()),
        )),
        Err(e) => {
            report.push(Check::new(
                CheckId::ReadHeader,
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
                CheckId::OpenIndexed,
                Some(alignment.to_path_buf()),
                "the index loaded and the file is ready for region queries",
            ));
            v
        }
        Err(e) => {
            // This is the whole point of the module. Do not repeat the mistake of the message
            // above, which puts the blame on the alignment. If the code already showed that there
            // is no index, *that* is the answer. An `ENOENT` that names the CRAM here means that
            // the reader could not load an index beside it. It does not mean that the CRAM went
            // away between two reads of it.
            let check = if has_index {
                Check::new(
                    CheckId::OpenIndexed,
                    Some(alignment.to_path_buf()),
                    Status::Fail,
                    format!(
                        "{e}\n         (this message names the alignment, but the file itself \
                         opened fine above — the failure is in its index or the reference)"
                    ),
                )
            } else {
                Check::new(
                    CheckId::OpenIndexed,
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

    // One real region query. Everything above can pass on a file whose index is stale or cut
    // short. This is the check that does a seek, and a seek is what the Y, mtDNA and SV paths
    // do.
    let Some(contig) = header
        .reference_sequences()
        .keys()
        .next()
        .map(|k| String::from_utf8_lossy(k.as_ref()).into_owned())
    else {
        report.push(Check::new(
            CheckId::RegionQuery,
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
                CheckId::RegionQuery,
                Some(alignment.to_path_buf()),
                Status::Fail,
                format!("decoding the first record of {contig} failed: {e}"),
            ),
            Some(Ok(_)) => Check::ok(
                CheckId::RegionQuery,
                None,
                format!("seeked to {contig} and decoded a record"),
            ),
            None => Check::new(
                CheckId::RegionQuery,
                None,
                Status::Warn,
                format!("seeked to {contig} but it holds no records"),
            ),
        },
        Err(e) => Check::new(
            CheckId::RegionQuery,
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

    /// A status bar makes EPERM and EACCES look alike, and their fixes have nothing to do with
    /// each other. So the explanation must separate them. Advice to run chmod, on a TCC denial,
    /// sends the user nowhere.
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

    /// Take a file that the code can read, and that has no index. The report must warn about the
    /// *index*, and it must name both of the forms that the readers accept. It must not let that
    /// warning look like a problem with the alignment. That confusion is what made the original
    /// bug report impossible to read.
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
            .find(|c| c.id == CheckId::CoordinateIndex)
            .expect("index is always checked");
        assert_eq!(index.status, Status::Warn, "a missing index has a sequential fallback");
        assert!(index.detail.contains("sample.bam.bai"), "{}", index.detail);
        assert!(index.detail.contains("sample.bai"), "{}", index.detail);

        // The alignment itself opened fine, so nothing may blame it for the missing index.
        let file = report
            .checks
            .iter()
            .find(|c| c.id == CheckId::AlignmentFile)
            .expect("file is always checked");
        assert_eq!(file.status, Status::Ok, "{report}");

        std::fs::remove_dir_all(&dir).ok();
    }

    /// The invariant that the batch depends on. A file with no index still analyzes correctly in
    /// a sequential pass, for coverage, read metrics and sex. The missing index, and the
    /// `open indexed` failure that follows it, must not read as "nobody can analyze this sample".
    /// The wrong answer here would drop the results of every CRAM without an index in a project,
    /// and nobody would see it happen.
    #[test]
    fn a_missing_index_does_not_block_sequential_reads() {
        let dir = std::env::temp_dir().join("navigator-preflight-blocking");
        std::fs::create_dir_all(&dir).unwrap();
        let bam = dir.join("sample.bam");
        std::fs::write(&bam, b"not really a bam").unwrap();

        let report = diagnose(&bam, None);
        assert!(report.failed(), "a garbage BAM fails somewhere: {report}");
        assert!(
            !report
                .checks
                .iter()
                .any(|c| c.status == Status::Fail && c.id == CheckId::CoordinateIndex),
            "a missing index is a warning, never a failure: {report}"
        );

        std::fs::remove_dir_all(&dir).ok();

        // An index-only failure is not a sequential blocker; an unreadable file is.
        assert!(!CheckId::CoordinateIndex.blocks_sequential_reads());
        assert!(!CheckId::OpenIndexed.blocks_sequential_reads());
        assert!(!CheckId::RegionQuery.blocks_sequential_reads());
        // A problem with the reference does not skip the sample on its own. A BAM reads without a
        // reference. And a CRAM that truly can not use one fails the header read, and that check
        // does skip the sample.
        assert!(!CheckId::ReferenceFasta.blocks_sequential_reads());
        assert!(!CheckId::ReferenceIndex.blocks_sequential_reads());
        assert!(CheckId::AlignmentFile.blocks_sequential_reads());
        assert!(CheckId::ReadHeader.blocks_sequential_reads());
    }

    /// An unreadable alignment blocks everything, so a batch may skip the sample outright.
    #[test]
    fn an_unreadable_alignment_blocks_sequential_reads() {
        let report = diagnose(Path::new("/nonexistent/sample.cram"), None);
        assert!(report.blocks_sequential_reads(), "{report}");
    }

    #[test]
    fn missing_alignment_fails_fast_and_names_itself() {
        let report = diagnose(Path::new("/nonexistent/sample.cram"), None);
        let first = report.first_failure().expect("missing file must fail");
        assert_eq!(first.id, CheckId::AlignmentFile);
        // Assert "not found" in a portable way. Unix reports ENOENT (2). Windows reports 2 or 3,
        // and which one depends on the component of the path that is absent. So key on the
        // message, and not on a Unix errno.
        assert!(first.detail.starts_with("not found"), "{}", first.detail);
        // The check on the reference must not run. The report stops at the first real block.
        assert!(
            !report.checks.iter().any(|c| c.id == CheckId::ReferenceFasta),
            "{report}"
        );
    }
}
