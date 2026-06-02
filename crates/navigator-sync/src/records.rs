//! AT Proto record contracts Navigator publishes.
//!
//! **No floats:** atproto records are DAG-CBOR, which has no float type — the PDS
//! rejects them. So every f64 metric (mean depth, % at depth, …) is encoded as a
//! string (lossless shortest round-trip) and parsed back by the consumer; only genuine
//! integers stay numeric. See documents/atmosphere/13-Local-PDS-Testing.md.

use serde::{Deserialize, Serialize};

/// Collection NSID for per-sample coverage summaries (under the navigator namespace).
pub const COVERAGE_SUMMARY_COLLECTION: &str = "com.decodingus.navigator.coverageSummary";

/// A public per-sample coverage summary record. Float metrics are strings (see module
/// docs); `genome_territory`/`callable_bases` are integers.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub struct CoverageSummaryRecord {
    #[serde(rename = "$type")]
    pub record_type: String,
    pub reference_build: String,
    pub mean_coverage: String,
    pub median_coverage: String,
    pub sd_coverage: String,
    pub pct_10x: String,
    pub pct_20x: String,
    pub pct_30x: String,
    pub genome_territory: i64,
    pub callable_bases: i64,
    pub created_at: String,
}

impl CoverageSummaryRecord {
    /// Build from coverage metrics, encoding floats as strings. `created_at` is RFC3339.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        reference_build: impl Into<String>,
        mean_coverage: f64,
        median_coverage: f64,
        sd_coverage: f64,
        pct_10x: f64,
        pct_20x: f64,
        pct_30x: f64,
        genome_territory: u64,
        callable_bases: u64,
        created_at: impl Into<String>,
    ) -> Self {
        CoverageSummaryRecord {
            record_type: COVERAGE_SUMMARY_COLLECTION.to_string(),
            reference_build: reference_build.into(),
            mean_coverage: mean_coverage.to_string(),
            median_coverage: median_coverage.to_string(),
            sd_coverage: sd_coverage.to_string(),
            pct_10x: pct_10x.to_string(),
            pct_20x: pct_20x.to_string(),
            pct_30x: pct_30x.to_string(),
            genome_territory: genome_territory as i64,
            callable_bases: callable_bases as i64,
            created_at: created_at.into(),
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
    fn coverage_summary_encodes_floats_as_strings() {
        let rec = CoverageSummaryRecord::new(
            "chm13v2.0", 178.81308467620255, 182.0, 28.9, 1.0, 1.0, 1.0, 16569, 16292, "2026-06-02T00:00:00Z",
        );
        let v = serde_json::to_value(&rec).unwrap();
        assert_eq!(v["$type"], COVERAGE_SUMMARY_COLLECTION);
        assert_eq!(v["referenceBuild"], "chm13v2.0");
        assert_eq!(v["meanCoverage"], "178.81308467620255"); // string, lossless
        assert_eq!(v["genomeTerritory"], 16569); // integer stays numeric
        assert_eq!(v["callableBases"], 16292);
        assert_no_floats(&v);

        // round-trips back to the same record.
        let back: CoverageSummaryRecord = serde_json::from_value(v).unwrap();
        assert_eq!(back, rec);
    }
}
