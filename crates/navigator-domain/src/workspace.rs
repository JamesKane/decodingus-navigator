//! The desktop workspace aggregate: Project → Biosample → SequenceRun → Alignment,
//! plus analysis artifacts. A reference-linked graph (the legacy Scala model used
//! string `atUri`/ref fields); here links are typed foreign keys.
//!
//! Read metrics live as flat fields (not a 22-tuple JSONB blob), per plan §3. Each
//! entity has a `New*` form without the DB-assigned id for inserts.

use chrono::{DateTime, Utc};
use du_domain::ids::SampleGuid;
use serde::{Deserialize, Serialize};

/// A research project grouping samples.
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

/// A sequencing run for a biosample, with summary read metrics as flat fields.
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

/// An alignment of a sequence run to a reference build.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Alignment {
    pub id: i64,
    pub sequence_run_id: i64,
    pub reference_build: String,
    pub aligner: String,
    pub variant_caller: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct NewAlignment {
    pub sequence_run_id: i64,
    pub reference_build: String,
    pub aligner: String,
    pub variant_caller: Option<String>,
}

/// A persisted analysis result, keyed by `(alignment, kind, algorithm_version)`. The
/// version is part of the key so a cache entry is invalidated when the algorithm
/// changes (plan §6 cache-versioning fix). `payload` is JSON of the result type.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AnalysisArtifact {
    pub id: i64,
    pub alignment_id: i64,
    pub kind: String,
    pub algorithm_version: String,
    pub created_at: DateTime<Utc>,
    pub payload: String,
}
