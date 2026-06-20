//! Asset integrity manifest (ancestry-ibd-asset-wiring, cross-cutting). `navigator-panelbuild`
//! writes `ancestry_manifest_<build>.json` listing each built `.bin`'s SHA-256; the app verifies a
//! loaded asset against it and refuses a mismatch — a cheap integrity guard for CDN-delivered assets.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::error::AnalysisError;

/// Lowercase-hex SHA-256 of `bytes`. Re-exported from the shared `du-bio` helper so existing
/// `manifest::sha256_hex` callers (e.g. `navigator-panelbuild`) keep working.
pub use du_bio::hash::sha256_hex;

/// One asset's integrity record.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AssetEntry {
    pub sha256: String,
    pub bytes: u64,
}

/// Per-build asset manifest: asset filename (on-disk name) → integrity record.
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct AssetManifest {
    pub build: String,
    #[serde(default)]
    pub generated_at: String,
    #[serde(default)]
    pub assets: BTreeMap<String, AssetEntry>,
}

impl AssetManifest {
    pub fn from_json(s: &str) -> Result<Self, AnalysisError> {
        serde_json::from_str(s).map_err(|e| AnalysisError::Message(format!("manifest decode: {e}")))
    }

    pub fn to_json(&self) -> Result<String, AnalysisError> {
        serde_json::to_string_pretty(self).map_err(|e| AnalysisError::Message(format!("manifest encode: {e}")))
    }

    /// Record `bytes` for `filename`.
    pub fn insert(&mut self, filename: impl Into<String>, bytes: &[u8]) {
        self.assets.insert(
            filename.into(),
            AssetEntry {
                sha256: sha256_hex(bytes),
                bytes: bytes.len() as u64,
            },
        );
    }

    /// Verify `bytes` for `filename`. `Ok` when the manifest has no entry for the file (advisory —
    /// unlisted assets aren't gated) or the digest matches; `Err(expected, got)` on a mismatch.
    pub fn verify(&self, filename: &str, bytes: &[u8]) -> Result<(), (String, String)> {
        if let Some(e) = self.assets.get(filename) {
            let got = sha256_hex(bytes);
            if got != e.sha256 {
                return Err((e.sha256.clone(), got));
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn insert_verify_and_json_round_trip() {
        let mut m = AssetManifest {
            build: "chm13v2.0".into(),
            generated_at: String::new(),
            assets: BTreeMap::new(),
        };
        m.insert("ancestry_panel_chm13v2.0.bin", b"hello");
        assert_eq!(m.assets["ancestry_panel_chm13v2.0.bin"].bytes, 5);
        // Matching bytes verify; tampered bytes are rejected; unlisted files pass (advisory).
        assert!(m.verify("ancestry_panel_chm13v2.0.bin", b"hello").is_ok());
        assert!(m.verify("ancestry_panel_chm13v2.0.bin", b"hELLo").is_err());
        assert!(m.verify("not_listed.bin", b"anything").is_ok());
        let back = AssetManifest::from_json(&m.to_json().unwrap()).unwrap();
        assert_eq!(back, m);
    }
}
