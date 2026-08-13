//! The last step of stage C: put the marked alignment where it belongs and index it.
//!
//! ## Why this emits BAM rather than CRAM
//!
//! CRAM is the better container on paper — reference-based compression makes a 30x WGS roughly
//! 17 GB against 60–80 GB of BAM — and stage C was built to produce it. Two defects in
//! `noodles-cram` 0.94 make that unreachable today, and both surfaced only at real scale:
//!
//! - **Writing** panics on a secondary alignment with `SEQ: *`, which is legal SAM and what
//!   minimap2 emits, because CRAM encodes a read as its difference from the reference and there
//!   are no bases to difference. [`super::cram`] handles it by dropping those records.
//! - **Indexing** then panics on any *multi-reference* slice: `cram::fs::index` decodes records
//!   with `fasta::Repository::default()` — an empty one, still carrying its `// TODO` upstream —
//!   and the reader `expect`s a name that cannot be there. With 25 contigs, every slice that
//!   straddles a contig boundary is multi-reference, so a whole-genome CRAM cannot be indexed at
//!   all. Only the coordinates are actually needed, which is why the decode is a bug rather than
//!   a requirement.
//!
//! Emitting BAM avoids both, and costs less than it looks: duplicate marking has already written
//! exactly the bytes that belong in the output, so finalising is a rename and an index instead of
//! a whole re-compression pass over the alignment. The trade is disk, and for the handful of
//! whole genomes a desktop user holds that is the cheaper side.
//!
//! The CRAM path is kept and still tested — the defects are upstream and fixable, and the choice
//! should be revisitable without rebuilding the stage.

use std::path::{Path, PathBuf};

use noodles::bam;

use crate::error::AnalysisError;

/// What finalising produced.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FinalizedAlignment {
    pub bam: PathBuf,
    pub index: PathBuf,
}

/// Move `input` to `output` and write its `.bai`.
///
/// A rename when both sit on one filesystem, which is the normal case — the scratch directory
/// lives beside the output precisely so this costs nothing. Falls back to a copy across devices.
pub fn finalize_bam(input: &Path, output: &Path) -> Result<FinalizedAlignment, AnalysisError> {
    if let Some(parent) = output.parent() {
        std::fs::create_dir_all(parent).map_err(|e| AnalysisError::io(parent, e))?;
    }

    if std::fs::rename(input, output).is_err() {
        std::fs::copy(input, output).map_err(|e| AnalysisError::io(output, e))?;
        std::fs::remove_file(input).map_err(|e| AnalysisError::io(input, e))?;
    }

    let index = index_bam(output)?;
    Ok(FinalizedAlignment {
        bam: output.to_path_buf(),
        index,
    })
}

/// Build the `.bai` for a finished BAM, returning its path.
///
/// Separate from [`finalize_bam`] because it reads the completed file back, so it is also what a
/// "reindex this alignment" repair would call.
pub fn index_bam(bam_path: &Path) -> Result<PathBuf, AnalysisError> {
    let index = bam::fs::index(bam_path).map_err(|e| AnalysisError::io(bam_path, e))?;
    let index_path = bai_path(bam_path);
    let mut writer = std::fs::File::create(&index_path)
        .map(bam::bai::io::Writer::new)
        .map_err(|e| AnalysisError::io(&index_path, e))?;
    writer
        .write_index(&index)
        .map_err(|e| AnalysisError::io(&index_path, e))?;
    Ok(index_path)
}

/// `foo.bam` → `foo.bam.bai`, the name every reader looks for.
pub fn bai_path(bam_path: &Path) -> PathBuf {
    let mut name = bam_path.as_os_str().to_os_string();
    name.push(".bai");
    PathBuf::from(name)
}
