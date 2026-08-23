//! The desktop workspace aggregate: Project → Biosample → SequenceRun → Alignment, plus analysis
//! artifacts. It is a graph of links. The legacy Scala model used string `atUri` and ref fields;
//! here a link is a typed foreign key.
//!
//! Read metrics live as flat fields, and not as a 22-tuple JSONB blob, as plan §3 says. Each entity
//! has a `New*` form for inserts, without the id that the DB assigns.

use chrono::{DateTime, Utc};
use du_domain::ids::SampleGuid;
use serde::{Deserialize, Serialize};

/// A research project that groups samples.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Project {
    pub id: i64,
    pub name: String,
    pub description: Option<String>,
    pub administrator: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct NewProject {
    pub name: String,
    pub description: Option<String>,
    pub administrator: String,
}

/// A biosample (donor sample). Identity is a stable cross-system `SampleGuid`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Biosample {
    pub guid: SampleGuid,
    pub sample_accession: Option<String>,
    pub donor_identifier: String,
    pub description: Option<String>,
    pub center_name: Option<String>,
    pub sex: Option<String>,
    pub project_id: Option<i64>,
}

impl Biosample {
    /// A biosample with only its identity set, and every descriptive field empty.
    ///
    /// This and its siblings ([`SequenceRun::new`], [`NewSequenceRun::new`], [`Alignment::new`],
    /// [`NewAlignment::new`]) take exactly the fields that have no sensible empty value. Some
    /// callers truly know only the identity: fixtures, and imports that fill the rest in later.
    /// These constructors stop those callers from a column of `None` values. A caller that does
    /// know more must say so with functional-update syntax:
    /// `Biosample { sex: Some("M".into()), ..Biosample::new(g, id) }`.
    pub fn new(guid: SampleGuid, donor_identifier: impl Into<String>) -> Self {
        Biosample {
            guid,
            sample_accession: None,
            donor_identifier: donor_identifier.into(),
            description: None,
            center_name: None,
            sex: None,
            project_id: None,
        }
    }
}

/// A sequencing run for a biosample, with summary read metrics as flat fields.
///
/// The lab and instrument identity block is `instrument_id`, `sample_name`, `library_id`,
/// `platform_unit` and `flowcell_id`. The import infers it from the alignment, with a read-name
/// scan and the `@RG` tags. It is the crowd-source input that resolves the sequencing facility.
/// `sequencing_facility` is the lab (FGC/FTDNA/YSEQ/Dante/Nebula…). A person sets it by hand for
/// now, and it will come from `instrument_id` when the AppView lookup endpoint ships (roadmap D8).
/// All of these are `None` until something fills them.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SequenceRun {
    pub id: i64,
    pub biosample_guid: SampleGuid,
    pub platform_name: String,
    pub instrument_model: Option<String>,
    pub test_type: String,
    pub library_layout: Option<String>,
    pub total_reads: Option<i64>,
    pub pf_reads_aligned: Option<i64>,
    pub mean_read_length: Option<f64>,
    pub mean_insert_size: Option<f64>,
    /// Exact total sequenced yield in base pairs (Σ read_length_histogram). This is the "Gbases"
    /// figure of the standardized test label. It is `None` until a read-metrics pass runs.
    pub total_bases: Option<i64>,
    /// Read chemistry or mode that the import inferred
    /// (`SHORT`/`HIFI`/`CLR`/`ONT_SIMPLEX`/`ONT_DUPLEX`). This is the long-read arm of the
    /// standardized test label. `None` until a library-stats scan runs.
    pub read_type: Option<String>,
    /// The sequencing laboratory (a [`crate::labs`] display name), e.g. "YSEQ", "Dante Labs".
    pub sequencing_facility: Option<String>,
    /// Most-frequent instrument serial from the read names / `@RG` (e.g. `A00123`, `m84…`).
    pub instrument_id: Option<String>,
    /// `@RG SM`: the sample name in the alignment tags (it can differ from the biosample).
    pub sample_name: Option<String>,
    /// `@RG LB`: library id (stable across realignments).
    pub library_id: Option<String>,
    /// `@RG PU`: platform unit (flowcell.lane.barcode).
    pub platform_unit: Option<String>,
    /// Most-frequent flowcell id from the read names.
    pub flowcell_id: Option<String>,
}

impl SequenceRun {
    /// A run with only the fields the database needs. The whole metrics block and lab-identity
    /// block stay `None`, which is exactly their state until an analysis pass fills them. See
    /// [`Biosample::new`] for why these constructors exist.
    pub fn new(
        id: i64,
        biosample_guid: SampleGuid,
        platform_name: impl Into<String>,
        test_type: impl Into<String>,
    ) -> Self {
        SequenceRun {
            id,
            biosample_guid,
            platform_name: platform_name.into(),
            instrument_model: None,
            test_type: test_type.into(),
            library_layout: None,
            total_reads: None,
            pf_reads_aligned: None,
            mean_read_length: None,
            mean_insert_size: None,
            total_bases: None,
            read_type: None,
            sequencing_facility: None,
            instrument_id: None,
            sample_name: None,
            library_id: None,
            platform_unit: None,
            flowcell_id: None,
        }
    }

    /// The standardized, vendor-neutral test label (`WGS150 45Gbases`, `HiFi 90Gbases`,
    /// `BigY-700`). `None` when this is not a yield or product test we standardize (chips, panels),
    /// and the caller then falls back to the raw `test_type`. See [`du_domain::testprofile`].
    pub fn standardized_label(&self) -> Option<String> {
        du_domain::testprofile::standardized_label(&du_domain::testprofile::RunProfile {
            test_type: Some(self.test_type.as_str()),
            platform: Some(self.platform_name.as_str()),
            read_type: self.read_type.as_deref(),
            mean_read_len: self.mean_read_length,
            total_bases: self.total_bases,
        })
    }
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NewSequenceRun {
    pub biosample_guid: SampleGuid,
    pub platform_name: String,
    pub instrument_model: Option<String>,
    pub test_type: String,
    pub library_layout: Option<String>,
    pub total_reads: Option<i64>,
    pub pf_reads_aligned: Option<i64>,
    pub mean_read_length: Option<f64>,
    pub mean_insert_size: Option<f64>,
}

impl NewSequenceRun {
    /// A run to insert, with only the required fields set. See [`Biosample::new`].
    pub fn new(biosample_guid: SampleGuid, platform_name: impl Into<String>, test_type: impl Into<String>) -> Self {
        NewSequenceRun {
            biosample_guid,
            platform_name: platform_name.into(),
            instrument_model: None,
            test_type: test_type.into(),
            library_layout: None,
            total_reads: None,
            pf_reads_aligned: None,
            mean_read_length: None,
            mean_insert_size: None,
        }
    }
}

/// An alignment of a sequence run to a reference build. `bam_path`/`reference_path`
/// locate the files so analysis can be run directly from the record.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Alignment {
    pub id: i64,
    pub sequence_run_id: i64,
    pub reference_build: String,
    pub aligner: String,
    pub variant_caller: Option<String>,
    pub bam_path: Option<String>,
    pub reference_path: Option<String>,
    /// SHA-256 of the content of the alignment file (hex). The import computes it, or the first
    /// analysis does for a file that came from a batch import. It is the content identity of the
    /// file, and it invalidates a cached analysis only when the file itself changes. `None` until
    /// something computes it.
    pub content_sha256: Option<String>,
    /// The alignment this one came from, for a row Navigator derived and did not import. `None`
    /// means this is an original: the alignment of a vendor, or anything a direct import made.
    ///
    /// Realignment sets it. Realignment maps the reads of a vendor alignment to another reference,
    /// and registers the result under the same `sequence_run_id`. That is the same physical
    /// library, mapped a different way. Without this field, a subject that has both builds has two
    /// alignments, and no way to tell which came from which.
    pub derived_from_alignment_id: Option<i64>,
    /// How it was derived, as `realign:<backend>-<preset>` (e.g. `realign:minimap2-sr`). `None`
    /// alongside a `None` parent. [`Alignment::aligner`] still carries the mapper alone.
    pub derivation: Option<String>,
}

impl Alignment {
    /// An alignment with only the fields the database needs: no file paths, no caller, and no
    /// derivation. That is an original, and not something Navigator made. See [`Biosample::new`].
    pub fn new(id: i64, sequence_run_id: i64, reference_build: impl Into<String>, aligner: impl Into<String>) -> Self {
        Alignment {
            id,
            sequence_run_id,
            reference_build: reference_build.into(),
            aligner: aligner.into(),
            variant_caller: None,
            bam_path: None,
            reference_path: None,
            content_sha256: None,
            derived_from_alignment_id: None,
            derivation: None,
        }
    }

    /// True when Navigator made this alignment from another one, and did not import it.
    ///
    /// The user sees this distinction. A derived alignment can go, and the code can build it again
    /// from its source. The UI must also say where it came from, and must not show it as something
    /// the vendor supplied.
    pub fn is_derived(&self) -> bool {
        self.derived_from_alignment_id.is_some()
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct NewAlignment {
    pub sequence_run_id: i64,
    pub reference_build: String,
    pub aligner: String,
    pub variant_caller: Option<String>,
    pub bam_path: Option<String>,
    pub reference_path: Option<String>,
    /// Content SHA-256 if already known at creation (else `None`; filled in lazily).
    pub content_sha256: Option<String>,
    /// The alignment this one came from; see [`Alignment::derived_from_alignment_id`].
    pub derived_from_alignment_id: Option<i64>,
    /// How it was derived; see [`Alignment::derivation`].
    pub derivation: Option<String>,
}

impl NewAlignment {
    /// An alignment to insert, with only the required fields set. See [`Biosample::new`].
    pub fn new(sequence_run_id: i64, reference_build: impl Into<String>, aligner: impl Into<String>) -> Self {
        NewAlignment {
            sequence_run_id,
            reference_build: reference_build.into(),
            aligner: aligner.into(),
            variant_caller: None,
            bam_path: None,
            reference_path: None,
            content_sha256: None,
            derived_from_alignment_id: None,
            derivation: None,
        }
    }
}

/// A persisted analysis result, with the key `(alignment, kind, algorithm_version)`. The version
/// is part of the key, so a change of the algorithm invalidates a cache entry (plan §6, the cache
/// version fix). `payload` is JSON of the result type.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AnalysisArtifact {
    pub id: i64,
    pub alignment_id: i64,
    pub kind: String,
    pub algorithm_version: String,
    pub created_at: DateTime<Utc>,
    pub payload: String,
    /// What made this result: `navigator-walk` (CRAM walk) or `pipeline-sidecar` (fast-path
    /// ingest). `None` for a row from before this field → the code reads it as `navigator-walk`.
    pub source: Option<String>,
    /// `full` or `partial` (e.g. lite coverage from sidecars, upgradeable by the deep pass).
    /// `None` → treated as `full`.
    pub completeness: Option<String>,
    /// The signature (`mtime:size`) of the source file when the code computed this artifact, for
    /// a staleness check. A BAM or CRAM that changed invalidates it. `None` for a row from before
    /// this field, and for a source that is not a file. The code reads those as fresh.
    pub source_sig: Option<String>,
}
