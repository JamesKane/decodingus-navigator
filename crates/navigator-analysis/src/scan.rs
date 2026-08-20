//! The scanner over a project directory. It is the port of the Scala `ProjectDirectoryScanner`.
//!
//! The layout on the NAS is `{projectRoot}/{sampleId}/files…`. Each subdirectory of the root is
//! one sample, and the code puts the files inside it into classes by role. The app turns the
//! result into rows: a Project, a Biosample, a SequenceRun and an Alignment.
//!
//! This module classifies files and nothing more. It uses no database and no noodles. In this
//! slice, only an alignment file, an index file and a variant file drive an import. The code
//! recognizes `coverage.txt`, `stats.txt` and `*.dragstr.model`, and it reads none of them,
//! because it computes the coverage again from the alignment.

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
    /// A `coverage.txt`, which somebody computed before. The code ignores it, because it computes
    /// the coverage again.
    Coverage,
    /// A `stats.txt` (precomputed; ignored).
    Stats,
    /// A `*.dragstr.model` (recorded; ignored).
    DragstrModel,
    /// Anything else.
    Other,
}

/// Put a file into a class, by its name. The case does not matter. The code checks an extension of
/// more than one part first, so `.g.vcf.gz` gives a Variant, and a bare `.gz` does not match it
/// first.
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

/// The sidecars that the `ytree` pipeline writes for each sample. The code matches them by the end
/// of the name. They are there only when that workflow ran on the sample, and absent from a
/// directory that holds an alignment alone. The fast-path ingest of the app reads these, and it
/// does not walk the CRAM.
///
/// This type serializes, so the app can record which files an alignment came from. The code finds
/// the sidecars in a directory scan, and it runs that scan once, at the import. Without a record
/// of the result, a later run of the fast path would have nothing to run against. Such a run
/// places a haplogroup again, against a newer tree.
#[derive(Debug, Clone, Default, serde::Serialize, serde::Deserialize)]
pub struct SampleSidecars {
    /// `*.chrY.g.vcf.gz`: the chrY GVCF at ploidy 1, for a male.
    pub chr_y_gvcf: Option<PathBuf>,
    /// `*.chrM.g.vcf.gz`: the chrM GVCF at ploidy 1.
    pub chr_m_gvcf: Option<PathBuf>,
    /// `*.callable.bed`: the CallableLoci track.
    pub callable_bed: Option<PathBuf>,
    /// `*.callable.summary.txt`: the count of bases in each state.
    pub callable_summary: Option<PathBuf>,
    /// `*.sex`: it holds `male` or `female`.
    pub sex: Option<PathBuf>,
    /// `coverage.txt`: the output of samtools coverage.
    pub coverage: Option<PathBuf>,
    /// `stats.txt`: the output of samtools stats.
    pub stats: Option<PathBuf>,
    /// `*.flagstat[.txt]`: the output of samtools flagstat. It is another source of read
    /// metrics.
    pub flagstat: Option<PathBuf>,
    /// The output of Picard `CollectWgsMetrics`, at `*wgs*metric*`. It holds the depth
    /// distribution over the whole genome.
    pub wgs_metrics: Option<PathBuf>,
    /// Picard `CollectAlignmentSummaryMetrics` (`*alignment_summary*`).
    pub alignment_summary: Option<PathBuf>,
    /// The build token that the code read out of the GVCF name, such as `chm13`. Use it to check
    /// that the GVCF and the alignment sit on the same build. That check comes before the fast
    /// path, which needs no liftover.
    pub build_hint: Option<String>,
}

impl SampleSidecars {
    /// True when the haplogroup fast path is available (at least one GVCF present).
    pub fn has_haplogroup_gvcf(&self) -> bool {
        self.chr_y_gvcf.is_some() || self.chr_m_gvcf.is_some()
    }
}

/// Find the pipeline sidecars among the files of a sample, by name. The case does not matter. The
/// code matches an exact suffix of more than one part against the whole file name. So
/// `*.chrY.g.vcf.gz.tbi`, which is an index, does not match `*.chrY.g.vcf.gz`.
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
    // The output of Picard, and that of flagstat, have no fixed name. So the code matches a part
    // of the file name, in lower case. The Scala scanner used the same open patterns.
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

    // Match the GVCF in both forms. One carries a sample name in front, as
    // `HG00096.chm13.chrY.g.vcf.gz`, which is the flat name that ytree writes to its work area.
    // The other is a bare file of one analysis, as `chrY.g.vcf.gz`, which is the GATK layout
    // inside the repo.
    //
    // The pattern has no dot in front, so both match. It still leaves out `chrY.vcf.gz`, which
    // holds called variants and is not a GVCF, and it leaves out the `.tbi` index.
    let chr_y_gvcf = by_suffix("chry.g.vcf.gz");
    let chr_m_gvcf = by_suffix("chrm.g.vcf.gz");
    let build_hint = chr_y_gvcf.as_ref().or(chr_m_gvcf.as_ref()).and_then(|p| build_token(p));

    SampleSidecars {
        chr_y_gvcf,
        chr_m_gvcf,
        // The output of CallableLoci carries the name `*.callable.bed` in the work area, and
        // `callable_status.bed` in the GATK repo layout. So match any `.bed` whose name holds
        // "callable".
        callable_bed: by_pred(&|n| n.ends_with(".bed") && n.contains("callable")),
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

/// A subdirectory of one sample, which holds one alignment file or variant file, or more.
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

/// Collect the files under `dir`, down to `max_depth`. It walks into a subdirectory, and it skips
/// a directory that is hidden.
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

/// Scan a project directory. Each subdirectory of it that is not hidden holds one sample. The
/// code drops a sample that has no alignment file and no variant file. It returns an error when
/// the path is absent,
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
        // The code classifies coverage.txt and stats.txt, and it puts neither into the list of
        // alignments or the list of variants.
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
        // The code matches the GVCF, and not its .tbi index.
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
    fn detects_bare_named_gatk_repo_sidecars() {
        // The real layout of the D2C repo. A sample directory holds a `CP086569.2/` analysis
        // subtree, and the bare-named files sit two levels down inside it: `chrY.g.vcf.gz` and
        // `callable_status.bed`. Those are not the `<sample>.chrY.g.vcf.gz` names of the work
        // area. scan_sample must still find them.
        let dir = scratch("bare-repo").join("1aceb711");
        for f in [
            "CP086569.2/chrYM.cram",
            "CP086569.2/chrYM.cram.crai",
            "CP086569.2/coverage.txt",
            "CP086569.2/stats.txt",
            "CP086569.2/gatk3/callable_status.bed",
            "CP086569.2/gatk4/chrY.g.vcf.gz",
            "CP086569.2/gatk4/chrY.g.vcf.gz.tbi",
            "CP086569.2/gatk4/chrY.vcf.gz", // called variants (not a GVCF) — must NOT be the sidecar
        ] {
            touch(dir.join(f));
        }

        let sample = scan_sample(&dir);
        let sc = &sample.sidecars;
        assert!(
            sc.has_haplogroup_gvcf(),
            "bare chrY.g.vcf.gz must be detected as the Y GVCF"
        );
        assert!(sc.chr_y_gvcf.as_ref().unwrap().ends_with("chrY.g.vcf.gz"));
        assert!(sc
            .callable_bed
            .as_ref()
            .is_some_and(|p| p.ends_with("callable_status.bed")));
        assert!(sc.coverage.is_some() && sc.stats.is_some());
        assert_eq!(sample.alignment_files.len(), 1, "the chrYM.cram");

        let _ = fs::remove_dir_all(dir.parent().unwrap());
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
