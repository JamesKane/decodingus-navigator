//! Project-directory scanner (port of the Scala `ProjectDirectoryScanner`). NAS layout:
//! `{projectRoot}/{sampleId}/files…` — each immediate subdirectory is one sample, and the
//! files within are classified by role. The app turns the result into Project → Biosample
//! → SequenceRun → Alignment rows.
//!
//! Pure filesystem classification: no DB, no noodles. Only alignment/index/variant files
//! drive import this slice; `coverage.txt`/`stats.txt`/`*.dragstr.model` are recognized
//! but not consumed (coverage is recomputed from the alignment).

use std::fs;
use std::path::{Path, PathBuf};

use crate::error::AnalysisError;

/// A file discovered in a sample directory, classified by role.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DiscoveredFileType {
    /// `.bam` / `.cram`.
    Alignment,
    /// `.bai` / `.crai` / `.tbi` / `.csi`.
    Index,
    /// `.vcf` / `.vcf.gz` / `.g.vcf.gz` / `.gvcf.gz`.
    Variant,
    /// A `coverage.txt` (precomputed; ignored — coverage is recomputed).
    Coverage,
    /// A `stats.txt` (precomputed; ignored).
    Stats,
    /// A `*.dragstr.model` (recorded; ignored).
    DragstrModel,
    /// Anything else.
    Other,
}

/// Classify a file by its (case-insensitive) name. Multi-part extensions are checked
/// first so `.g.vcf.gz` is a Variant, not matched by a bare `.gz`.
pub fn classify(name: &str) -> DiscoveredFileType {
    let lower = name.to_ascii_lowercase();
    const VARIANT: [&str; 4] = [".g.vcf.gz", ".gvcf.gz", ".vcf.gz", ".vcf"];
    if VARIANT.iter().any(|p| lower.ends_with(p)) {
        return DiscoveredFileType::Variant;
    }
    if lower == "coverage.txt" {
        return DiscoveredFileType::Coverage;
    }
    if lower == "stats.txt" {
        return DiscoveredFileType::Stats;
    }
    if lower.ends_with(".dragstr.model") {
        return DiscoveredFileType::DragstrModel;
    }
    match lower.rsplit('.').next().unwrap_or("") {
        "bam" | "cram" => DiscoveredFileType::Alignment,
        "bai" | "crai" | "tbi" | "csi" => DiscoveredFileType::Index,
        _ => DiscoveredFileType::Other,
    }
}

/// A discovered file with its classified role.
#[derive(Debug, Clone)]
pub struct DiscoveredFile {
    pub path: PathBuf,
    pub kind: DiscoveredFileType,
}

/// The `ytree` pipeline's per-sample sidecars, matched by name suffix. Present only when the
/// sample was processed by that workflow; absent for a plain alignment-only directory. The
/// app's fast-path ingest reads these instead of walking the CRAM.
#[derive(Debug, Clone, Default)]
pub struct SampleSidecars {
    /// `*.chrY.g.vcf.gz` — ploidy-1 chrY GVCF (males).
    pub chr_y_gvcf: Option<PathBuf>,
    /// `*.chrM.g.vcf.gz` — ploidy-1 chrM GVCF.
    pub chr_m_gvcf: Option<PathBuf>,
    /// `*.callable.bed` — CallableLoci track.
    pub callable_bed: Option<PathBuf>,
    /// `*.callable.summary.txt` — per-state base counts.
    pub callable_summary: Option<PathBuf>,
    /// `*.sex` — `male` / `female`.
    pub sex: Option<PathBuf>,
    /// `coverage.txt` — samtools coverage.
    pub coverage: Option<PathBuf>,
    /// `stats.txt` — samtools stats.
    pub stats: Option<PathBuf>,
    /// `*.flagstat[.txt]` — samtools flagstat (an alternative read-metrics source).
    pub flagstat: Option<PathBuf>,
    /// Picard `CollectWgsMetrics` output (`*wgs*metric*`) — the genome-wide depth distribution.
    pub wgs_metrics: Option<PathBuf>,
    /// Picard `CollectAlignmentSummaryMetrics` (`*alignment_summary*`).
    pub alignment_summary: Option<PathBuf>,
    /// Build token parsed from the GVCF name (e.g. `chm13`), for confirming the GVCF and the
    /// alignment share a build before the liftover-free fast path is taken.
    pub build_hint: Option<String>,
}

impl SampleSidecars {
    /// True when the haplogroup fast path is available (at least one GVCF present).
    pub fn has_haplogroup_gvcf(&self) -> bool {
        self.chr_y_gvcf.is_some() || self.chr_m_gvcf.is_some()
    }
}

/// Detect pipeline sidecars among a sample's files by (case-insensitive) name. Specific
/// multi-part suffixes are matched against the full file name, so `*.chrY.g.vcf.gz.tbi`
/// (an index) does not match `*.chrY.g.vcf.gz`.
fn detect_sidecars(files: &[DiscoveredFile]) -> SampleSidecars {
    let by_suffix = |suffix: &str| {
        files
            .iter()
            .find(|f| {
                f.path
                    .file_name()
                    .and_then(|n| n.to_str())
                    .is_some_and(|n| n.to_ascii_lowercase().ends_with(suffix))
            })
            .map(|f| f.path.clone())
    };
    let by_name = |name: &str| {
        files
            .iter()
            .find(|f| {
                f.path
                    .file_name()
                    .and_then(|n| n.to_str())
                    .is_some_and(|n| n.eq_ignore_ascii_case(name))
            })
            .map(|f| f.path.clone())
    };
    // Picard/flagstat outputs have no fixed name — match a substring of the (lower-cased) filename
    // (mirrors the Scala scanner's loose patterns).
    let by_pred = |pred: &dyn Fn(&str) -> bool| {
        files
            .iter()
            .find(|f| {
                f.path
                    .file_name()
                    .and_then(|n| n.to_str())
                    .is_some_and(|n| pred(&n.to_ascii_lowercase()))
            })
            .map(|f| f.path.clone())
    };

    let chr_y_gvcf = by_suffix(".chry.g.vcf.gz");
    let chr_m_gvcf = by_suffix(".chrm.g.vcf.gz");
    let build_hint = chr_y_gvcf.as_ref().or(chr_m_gvcf.as_ref()).and_then(|p| build_token(p));

    SampleSidecars {
        chr_y_gvcf,
        chr_m_gvcf,
        callable_bed: by_suffix(".callable.bed"),
        callable_summary: by_suffix(".callable.summary.txt"),
        sex: by_suffix(".sex"),
        coverage: by_name("coverage.txt"),
        stats: by_name("stats.txt"),
        flagstat: by_pred(&|n| n.contains("flagstat")),
        wgs_metrics: by_pred(&|n| n.contains("wgs") && n.contains("metric")),
        alignment_summary: by_pred(&|n| n.contains("alignment_summary")),
        build_hint,
    }
}

/// The build segment of a GVCF name, e.g. `HG00096.chm13.chrY.g.vcf.gz` → `chm13`.
fn build_token(gvcf: &Path) -> Option<String> {
    let name = gvcf.file_name()?.to_str()?.to_ascii_lowercase();
    let stem = name
        .strip_suffix(".chry.g.vcf.gz")
        .or_else(|| name.strip_suffix(".chrm.g.vcf.gz"))?;
    stem.rsplit('.').next().filter(|s| !s.is_empty()).map(|s| s.to_string())
}

/// A sample subdirectory holding at least one alignment or variant file.
#[derive(Debug, Clone)]
pub struct DiscoveredSample {
    /// Subdirectory name (typically a sample alias, e.g. `HG00096`).
    pub sample_id: String,
    pub directory: PathBuf,
    pub alignment_files: Vec<PathBuf>,
    pub index_files: Vec<PathBuf>,
    pub variant_files: Vec<PathBuf>,
    pub all_files: Vec<DiscoveredFile>,
    /// Pipeline sidecars for this sample, if present (drives the fast-path ingest).
    pub sidecars: SampleSidecars,
}

/// A project directory and its discovered samples.
#[derive(Debug, Clone)]
pub struct DiscoveredProject {
    /// ENA accession or directory name (e.g. `PRJEB31736`).
    pub project_id: String,
    pub directory: PathBuf,
    pub samples: Vec<DiscoveredSample>,
}

fn is_hidden(path: &Path) -> bool {
    path.file_name()
        .and_then(|n| n.to_str())
        .is_some_and(|n| n.starts_with('.'))
}

/// Recursively collect files under `dir` up to `max_depth`, skipping hidden directories.
fn list_files_recursive(dir: &Path, max_depth: usize, depth: usize, out: &mut Vec<PathBuf>) {
    let entries = match fs::read_dir(dir) {
        Ok(e) => e,
        Err(_) => return,
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_file() {
            out.push(path);
        } else if path.is_dir() && depth < max_depth && !is_hidden(&path) {
            list_files_recursive(&path, max_depth, depth + 1, out);
        }
    }
}

/// Scan one sample directory into a [`DiscoveredSample`] (alignment/index/variant files + pipeline
/// sidecars). Always returns a sample; the caller decides whether it holds usable data ([`scan`]
/// drops samples with neither an alignment nor a variant file). Used directly by the app to ingest
/// a single staged sample directory onto an existing subject.
pub fn scan_sample(dir: &Path) -> DiscoveredSample {
    let mut files = Vec::new();
    list_files_recursive(dir, 2, 0, &mut files);
    files.sort();

    let all_files: Vec<DiscoveredFile> = files
        .into_iter()
        .map(|path| {
            let kind = path
                .file_name()
                .and_then(|n| n.to_str())
                .map_or(DiscoveredFileType::Other, classify);
            DiscoveredFile { path, kind }
        })
        .collect();

    let collect = |k: DiscoveredFileType| {
        all_files
            .iter()
            .filter(|f| f.kind == k)
            .map(|f| f.path.clone())
            .collect::<Vec<_>>()
    };

    let sidecars = detect_sidecars(&all_files);
    DiscoveredSample {
        sample_id: dir.file_name().and_then(|n| n.to_str()).unwrap_or_default().to_string(),
        directory: dir.to_path_buf(),
        alignment_files: collect(DiscoveredFileType::Alignment),
        index_files: collect(DiscoveredFileType::Index),
        variant_files: collect(DiscoveredFileType::Variant),
        all_files,
        sidecars,
    }
}

/// Scan a project directory: each immediate (non-hidden) subdirectory is a sample. Samples
/// with neither an alignment nor a variant file are dropped. Errors if the path is missing,
/// not a directory, has no subdirectories, or yields no samples with data.
pub fn scan(project_dir: &Path) -> Result<DiscoveredProject, AnalysisError> {
    if !project_dir.exists() {
        return Err(AnalysisError::Message(format!(
            "directory does not exist: {}",
            project_dir.display()
        )));
    }
    if !project_dir.is_dir() {
        return Err(AnalysisError::Message(format!(
            "not a directory: {}",
            project_dir.display()
        )));
    }
    let project_id = project_dir
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("project")
        .to_string();

    let mut subdirs: Vec<PathBuf> = fs::read_dir(project_dir)
        .map_err(|e| AnalysisError::io(project_dir, e))?
        .flatten()
        .map(|e| e.path())
        .filter(|p| p.is_dir() && !is_hidden(p))
        .collect();
    subdirs.sort();
    if subdirs.is_empty() {
        return Err(AnalysisError::Message(format!(
            "no sample subdirectories in {}",
            project_dir.display()
        )));
    }

    let samples: Vec<DiscoveredSample> = subdirs
        .iter()
        .map(|d| scan_sample(d))
        .filter(|s| !s.alignment_files.is_empty() || !s.variant_files.is_empty())
        .collect();
    if samples.is_empty() {
        return Err(AnalysisError::Message(format!(
            "no samples with data files in {}",
            project_dir.display()
        )));
    }

    Ok(DiscoveredProject {
        project_id,
        directory: project_dir.to_path_buf(),
        samples,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classify_rules() {
        use DiscoveredFileType::*;
        assert_eq!(classify("HG00096.chm13.g.vcf.gz"), Variant); // multi-part before .gz
        assert_eq!(classify("HG00096.chm13.mito.vcf.gz"), Variant);
        assert_eq!(classify("x.vcf"), Variant);
        assert_eq!(classify("HG00096.chm13.cram"), Alignment);
        assert_eq!(classify("HG00096.chm13.CRAM"), Alignment); // case-insensitive
        assert_eq!(classify("HG00096.chm13.cram.crai"), Index);
        assert_eq!(classify("HG00096.chm13.g.vcf.gz.tbi"), Index);
        assert_eq!(classify("coverage.txt"), Coverage);
        assert_eq!(classify("stats.txt"), Stats);
        assert_eq!(classify("HG00096.dragstr.model"), DragstrModel);
        assert_eq!(classify("HG00096.chm13.mito.vcf.gz.stats"), Other);
        assert_eq!(classify("notes.md"), Other);
    }

    /// Unique scratch dir under the system temp dir (no tempfile dep).
    fn scratch(tag: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("dun-scan-{}-{}", std::process::id(), tag));
        let _ = fs::remove_dir_all(&dir);
        dir
    }

    fn touch(path: PathBuf) {
        fs::create_dir_all(path.parent().unwrap()).unwrap();
        fs::write(path, b"x").unwrap();
    }

    #[test]
    fn scans_project_tree_grouping_and_skipping() {
        let root = scratch("prj");
        // Two real samples + an empty dir + a hidden dir.
        for s in ["HG00096", "HG00097"] {
            touch(root.join(s).join(format!("{s}.chm13.cram")));
            touch(root.join(s).join(format!("{s}.chm13.cram.crai")));
            touch(root.join(s).join(format!("{s}.chm13.mito.vcf.gz")));
            touch(root.join(s).join("coverage.txt"));
            touch(root.join(s).join("stats.txt"));
        }
        fs::create_dir_all(root.join("EMPTY")).unwrap();
        touch(root.join(".hidden").join("HGXXXX.cram")); // hidden dir → skipped

        let project = scan(&root).unwrap();
        assert_eq!(project.project_id, root.file_name().unwrap().to_str().unwrap());
        assert_eq!(project.samples.len(), 2, "empty + hidden dirs excluded");

        let s = &project.samples[0];
        assert_eq!(s.sample_id, "HG00096");
        assert_eq!(s.alignment_files.len(), 1);
        assert_eq!(s.index_files.len(), 1);
        assert_eq!(s.variant_files.len(), 1);
        // coverage.txt / stats.txt are classified but not in the alignment/variant lists.
        assert!(s.all_files.iter().any(|f| f.kind == DiscoveredFileType::Coverage));
        assert!(s.all_files.iter().any(|f| f.kind == DiscoveredFileType::Stats));

        let _ = fs::remove_dir_all(&root);
    }

    #[test]
    fn detects_pipeline_sidecars() {
        let root = scratch("sidecars");
        let s = "HG00096";
        let d = root.join(s);
        for f in [
            "HG00096.chm13.cram",
            "HG00096.chm13.cram.crai",
            "HG00096.chm13.chrY.g.vcf.gz",
            "HG00096.chm13.chrY.g.vcf.gz.tbi",
            "HG00096.chm13.chrM.g.vcf.gz",
            "HG00096.chm13.chrM.g.vcf.gz.tbi",
            "HG00096.chm13.chrYM.callable.bed",
            "HG00096.chm13.chrYM.callable.summary.txt",
            "HG00096.chm13.sex",
            "coverage.txt",
            "stats.txt",
        ] {
            touch(d.join(f));
        }

        let project = scan(&root).unwrap();
        let sc = &project.samples[0].sidecars;
        assert!(sc.has_haplogroup_gvcf());
        // The GVCF is matched, not its .tbi index.
        assert!(sc
            .chr_y_gvcf
            .as_ref()
            .unwrap()
            .to_str()
            .unwrap()
            .ends_with("chrY.g.vcf.gz"));
        assert!(sc
            .chr_m_gvcf
            .as_ref()
            .unwrap()
            .to_str()
            .unwrap()
            .ends_with("chrM.g.vcf.gz"));
        assert!(sc.callable_bed.is_some());
        assert!(sc.callable_summary.is_some());
        assert!(sc.sex.is_some());
        assert!(sc.coverage.is_some());
        assert!(sc.stats.is_some());
        assert_eq!(sc.build_hint.as_deref(), Some("chm13"));

        let _ = fs::remove_dir_all(&root);
    }

    #[test]
    fn no_sidecars_for_plain_alignment_dir() {
        let root = scratch("plain");
        touch(root.join("S1").join("S1.cram"));
        touch(root.join("S1").join("S1.cram.crai"));
        let project = scan(&root).unwrap();
        assert!(!project.samples[0].sidecars.has_haplogroup_gvcf());
        assert!(project.samples[0].sidecars.build_hint.is_none());
        let _ = fs::remove_dir_all(&root);
    }

    #[test]
    fn errors_when_no_samples_have_data() {
        let root = scratch("nodata");
        touch(root.join("README").join("notes.md")); // a subdir with no alignment/variant
        assert!(scan(&root).is_err());
        let _ = fs::remove_dir_all(&root);
    }
}
