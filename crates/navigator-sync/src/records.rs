//! AT Proto record contracts Navigator publishes.
//!
//! **No floats:** atproto records are DAG-CBOR, which has no float type — the PDS
//! rejects them. So every f64 metric (mean depth, % at depth, …) is encoded as a
//! string (lossless shortest round-trip) and parsed back by the consumer; only genuine
//! integers stay numeric. See documents/atmosphere/13-Local-PDS-Testing.md.

use serde::{Deserialize, Serialize};

// NOTE: the per-sample coverage summary record now lives in the shared
// `du_domain::fed::AlignmentRecord` (collection `com.decodingus.atmosphere.alignment`),
// so the publisher and the AppView's Jetstream consumer share one contract. The old
// `com.decodingus.navigator.coverageSummary` shape was a drift the AppView never ingested.

/// Collection NSID for a sample's de-novo / private variant calls on one contig.
pub const PRIVATE_VARIANTS_COLLECTION: &str = "com.decodingus.navigator.privateVariants";

/// One variant call. `allele_fraction` is a string (no floats); positions/depths integers.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct VariantCallEntry {
    pub position: i64,
    pub reference: String,
    pub alternate: String,
    pub depth: i64,
    pub alt_depth: i64,
    pub allele_fraction: String,
}

impl VariantCallEntry {
    pub fn new(
        position: i64,
        reference: char,
        alternate: char,
        depth: u32,
        alt_depth: u32,
        allele_fraction: f64,
    ) -> Self {
        VariantCallEntry {
            position,
            reference: reference.to_string(),
            alternate: alternate.to_string(),
            depth: depth as i64,
            alt_depth: alt_depth as i64,
            allele_fraction: allele_fraction.to_string(),
        }
    }
}

/// A sample's de-novo / private variant calls for one contig (the writer-side payload;
/// branch-creation *proposals* to the AppView curation API are a separate, later path).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct PrivateVariantsRecord {
    #[serde(rename = "$type")]
    pub record_type: String,
    pub contig: String,
    pub caller_version: String,
    pub created_at: String,
    pub variants: Vec<VariantCallEntry>,
}

impl PrivateVariantsRecord {
    pub fn new(
        contig: impl Into<String>,
        caller_version: impl Into<String>,
        created_at: impl Into<String>,
        variants: Vec<VariantCallEntry>,
    ) -> Self {
        PrivateVariantsRecord {
            record_type: PRIVATE_VARIANTS_COLLECTION.to_string(),
            contig: contig.into(),
            caller_version: caller_version.into(),
            created_at: created_at.into(),
            variants,
        }
    }
}

// ---- haplogroup reconciliation -------------------------------------------

/// Collection NSID for the multi-run haplogroup reconciliation record (donor-level
/// consensus across biosamples/runs). Consumed by the AppView (see
/// documents/atmosphere/09-Reconciliation-Records.md).
pub const HAPLOGROUP_RECONCILIATION_COLLECTION: &str = "com.decodingus.atmosphere.haplogroupReconciliation";

/// `com.decodingus.atmosphere.defs#recordMeta`: version + timestamps for sync.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct RecordMeta {
    pub version: i64,
    pub created_at: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub updated_at: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub last_modified_field: Option<String>,
}

/// `#reconciliationStatus`: the summary consensus across runs. The lexicon types
/// `confidence`/`branchCompatibilityScore`/`snpConcordance` as floats, but DAG-CBOR has
/// no float type (module docs), so they ride as strings like every other f64 metric.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ReconciliationStatusRecord {
    pub compatibility_level: String,
    pub consensus_haplogroup: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub confidence: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub divergence_point: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub branch_compatibility_score: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub snp_concordance: Option<String>,
    pub run_count: i64,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub warnings: Vec<String>,
}

/// `#runHaplogroupCall`: one run's call with quality metrics (floats → strings).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct RunHaplogroupCallRecord {
    pub source_ref: String,
    pub haplogroup: String,
    pub confidence: String,
    pub call_method: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub score: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub supporting_snps: Option<i64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub conflicting_snps: Option<i64>,
}

/// `#heteroplasmyObservation`: a mixed mtDNA position (`majorAlleleFrequency` → string).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct HeteroplasmyObservationRecord {
    pub position: i64,
    pub major_allele: String,
    pub minor_allele: String,
    pub major_allele_frequency: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub depth: Option<i64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub is_defining_snp: Option<bool>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub affected_haplogroup: Option<String>,
}

/// `#identityVerification`: same-individual metrics (kinship/concordance → strings).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct IdentityVerificationRecord {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub kinship_coefficient: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub fingerprint_snp_concordance: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub y_str_distance: Option<i64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub verification_status: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub verification_method: Option<String>,
}

/// The manual-override sub-object of the reconciliation record.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct ManualOverrideRecord {
    pub overridden_haplogroup: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub reason: Option<String>,
    pub overridden_at: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub overridden_by: Option<String>,
}

/// `#auditEntry`: one reconciliation-history entry.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct AuditEntryRecord {
    pub timestamp: String,
    pub action: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub previous_consensus: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub new_consensus: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub run_ref: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub notes: Option<String>,
}

/// The donor-level multi-run haplogroup reconciliation record. The app maps its domain
/// types into these primitives; `at_uri` is assigned by the PDS on create, so it is
/// written empty and the consumer reads the record's own URI.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct HaplogroupReconciliationRecord {
    #[serde(rename = "$type")]
    pub record_type: String,
    pub at_uri: String,
    pub meta: RecordMeta,
    pub specimen_donor_ref: String,
    pub dna_type: String,
    pub status: ReconciliationStatusRecord,
    pub run_calls: Vec<RunHaplogroupCallRecord>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub heteroplasmy_observations: Vec<HeteroplasmyObservationRecord>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub identity_verification: Option<IdentityVerificationRecord>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub last_reconciliation_at: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub manual_override: Option<ManualOverrideRecord>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub audit_log: Vec<AuditEntryRecord>,
}

impl HaplogroupReconciliationRecord {
    /// Assemble a reconciliation record. `created_at`/`last_reconciliation_at` are RFC3339.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        specimen_donor_ref: impl Into<String>,
        dna_type: impl Into<String>,
        created_at: impl Into<String>,
        status: ReconciliationStatusRecord,
        run_calls: Vec<RunHaplogroupCallRecord>,
        heteroplasmy_observations: Vec<HeteroplasmyObservationRecord>,
        identity_verification: Option<IdentityVerificationRecord>,
        manual_override: Option<ManualOverrideRecord>,
        audit_log: Vec<AuditEntryRecord>,
    ) -> Self {
        let created_at = created_at.into();
        HaplogroupReconciliationRecord {
            record_type: HAPLOGROUP_RECONCILIATION_COLLECTION.to_string(),
            at_uri: String::new(),
            meta: RecordMeta {
                version: 1,
                created_at: created_at.clone(),
                updated_at: None,
                last_modified_field: None,
            },
            specimen_donor_ref: specimen_donor_ref.into(),
            dna_type: dna_type.into(),
            status,
            run_calls,
            heteroplasmy_observations,
            identity_verification,
            last_reconciliation_at: Some(created_at),
            manual_override,
            audit_log,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// No JSON value anywhere in the record is a float (the atproto constraint).
    fn assert_no_floats(v: &serde_json::Value) {
        match v {
            serde_json::Value::Number(n) => assert!(!n.is_f64(), "float in record: {n}"),
            serde_json::Value::Array(a) => a.iter().for_each(assert_no_floats),
            serde_json::Value::Object(o) => o.values().for_each(assert_no_floats),
            _ => {}
        }
    }

    #[test]
    fn private_variants_encode_allele_fraction_as_string() {
        let rec = PrivateVariantsRecord::new(
            "chrM",
            "haploid-denovo-2",
            "2026-06-02T00:00:00Z",
            vec![
                VariantCallEntry::new(2, 'C', 'A', 4, 4, 1.0),
                VariantCallEntry::new(16302, 'T', 'C', 66, 35, 0.5303030303030303),
            ],
        );
        let v = serde_json::to_value(&rec).unwrap();
        assert_eq!(v["$type"], PRIVATE_VARIANTS_COLLECTION);
        assert_eq!(v["contig"], "chrM");
        assert_eq!(v["variants"][0]["position"], 2); // integer
        assert_eq!(v["variants"][0]["alleleFraction"], "1"); // string
        assert_eq!(v["variants"][1]["alleleFraction"], "0.5303030303030303");
        assert_no_floats(&v); // including inside the array

        let back: PrivateVariantsRecord = serde_json::from_value(v).unwrap();
        assert_eq!(back, rec);
    }

    #[test]
    fn reconciliation_record_shape_and_no_floats() {
        let status = ReconciliationStatusRecord {
            compatibility_level: "COMPATIBLE".into(),
            consensus_haplogroup: "R-FGC29067".into(),
            confidence: Some(0.94_f64.to_string()),
            divergence_point: None,
            branch_compatibility_score: Some(1.0_f64.to_string()),
            snp_concordance: Some(0.997_f64.to_string()),
            run_count: 2,
            warnings: vec![],
        };
        let calls = vec![RunHaplogroupCallRecord {
            source_ref: "aln:7".into(),
            haplogroup: "R-FGC29067".into(),
            confidence: 0.94_f64.to_string(),
            call_method: "SNP_PHYLOGENETIC".into(),
            score: Some(0.94_f64.to_string()),
            supporting_snps: Some(118),
            conflicting_snps: Some(2),
        }];
        let het = vec![HeteroplasmyObservationRecord {
            position: 16093,
            major_allele: "T".into(),
            minor_allele: "C".into(),
            major_allele_frequency: 0.82_f64.to_string(),
            depth: Some(120),
            is_defining_snp: Some(false),
            affected_haplogroup: None,
        }];
        let rec = HaplogroupReconciliationRecord::new(
            "biosample:abc",
            "Y_DNA",
            "2026-06-03T00:00:00Z",
            status,
            calls,
            het,
            None,
            Some(ManualOverrideRecord {
                overridden_haplogroup: "R-FGC29067".into(),
                reason: Some("Sanger-confirmed".into()),
                overridden_at: "2026-06-03T00:00:00Z".into(),
                overridden_by: Some("did:plc:xyz".into()),
            }),
            vec![AuditEntryRecord {
                timestamp: "2026-06-03T00:00:00Z".into(),
                action: "MANUAL_OVERRIDE".into(),
                previous_consensus: Some("R-FGC29071".into()),
                new_consensus: Some("R-FGC29067".into()),
                run_ref: None,
                notes: None,
            }],
        );
        let v = serde_json::to_value(&rec).unwrap();
        assert_eq!(v["$type"], HAPLOGROUP_RECONCILIATION_COLLECTION);
        assert_eq!(v["dnaType"], "Y_DNA");
        assert_eq!(v["meta"]["version"], 1);
        assert_eq!(v["status"]["consensusHaplogroup"], "R-FGC29067");
        assert_eq!(v["status"]["confidence"], "0.94"); // string, not float
        assert_eq!(v["runCalls"][0]["supportingSnps"], 118); // integer stays numeric
        assert_eq!(v["heteroplasmyObservations"][0]["majorAlleleFrequency"], "0.82");
        assert_eq!(v["manualOverride"]["reason"], "Sanger-confirmed");
        assert_eq!(v["auditLog"][0]["action"], "MANUAL_OVERRIDE");
        assert_no_floats(&v);

        let back: HaplogroupReconciliationRecord = serde_json::from_value(v).unwrap();
        assert_eq!(back, rec);
    }
}
