//! The last step of stage C: put the marked alignment where it belongs and index it.
//!
//! ## Why this gives a BAM, and not a CRAM
//!
//! On paper a CRAM is the better container. Its compression against a reference makes a 30x WGS
//! about 17 GB, against 60 to 80 GB for a BAM. Stage C exists to make one. Two defects in
//! `noodles-cram` 0.94 put that out of reach today, and both appeared only at real scale:
//!
//! - The **write** panics on a secondary alignment with a `SEQ` of `*`. That is legal SAM, and it
//!   is what minimap2 gives. A CRAM encodes a read as its difference from the reference, and there
//!   is no base to take a difference of. [`super::cram`] handles it, and drops those records.
//! - The **index** then panics on any slice that holds *more than one reference*.
//!   `cram::fs::index` decodes the records with a `fasta::Repository::default()`, which is empty
//!   and still carries its `// TODO` upstream. The reader then calls `expect` on a name that can
//!   not be there. With 25 contigs, every slice that crosses a contig boundary holds more than one
//!   reference, so nothing can index a whole-genome CRAM at all. The index needs the coordinates
//!   alone, and that is why the decode is a bug, and not a requirement.
//!
//! A BAM avoids both, and it costs less than it looks. The duplicate mark already wrote exactly
//! the bytes that belong in the output. So this last step is a rename and an index, and not a
//! second compression pass over the whole alignment. The cost is disk space. For the few whole
//! genomes that a desktop user holds, that is the cheaper side.
//!
//! The CRAM path stays, and the tests still cover it. The defects are upstream, and somebody can
//! fix them. The choice must then be open again, and nobody should have to build the stage
//! anew.

use std::path::{Path, PathBuf};

use noodles::bam;

use crate::error::AnalysisError;

/// What this last step made.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FinalizedAlignment {
    pub bam: PathBuf,
    pub index: PathBuf,
}

/// Move `input` to `output` and write its `.bai`.
///
/// This is a rename when both paths sit on one filesystem, and that is the usual case. The scratch
/// directory lives beside the output for exactly that reason, so the move costs nothing. Across
/// two devices it falls back to a copy.
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

/// Build the `.bai` of a finished BAM. Returns its path.
///
/// It sits apart from [`finalize_bam`] because it reads the completed file back. So it is also the
/// function that a repair which makes the index again would call.
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
