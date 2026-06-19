//! Asset integrity manifest (ancestry-ibd-asset-wiring, cross-cutting). `navigator-panelbuild`
//! writes `ancestry_manifest_<build>.json` listing each built `.bin`'s SHA-256; the app verifies a
//! loaded asset against it and refuses a mismatch — a cheap integrity guard for CDN-delivered assets.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

use crate::error::AnalysisError;

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
        self.assets.insert(filename.into(), AssetEntry { sha256: sha256_hex(bytes), bytes: bytes.len() as u64 });
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

/// Lowercase hex SHA-256 of `bytes`.
pub fn sha256_hex(bytes: &[u8]) -> String {
    Sha256::digest(bytes).iter().map(|b| format!("{b:02x}")).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn insert_verify_and_json_round_trip() {
        let mut m = AssetManifest { build: "chm13v2.0".into(), generated_at: String::new(), assets: BTreeMap::new() };
        m.insert("ancestry_panel_chm13v2.0.bin", b"hello");
        assert_eq!(m.assets["ancestry_panel_chm13v2.0.bin"].bytes, 5);
        // Matching bytes verify; tampered bytes are rejected; unlisted files pass (advisory).
        assert!(m.verify("ancestry_panel_chm13v2.0.bin", b"hello").is_ok());
        assert!(m.verify("ancestry_panel_chm13v2.0.bin", b"hELLo").is_err());
        assert!(m.verify("not_listed.bin", b"anything").is_ok());
        let back = AssetManifest::from_json(&m.to_json().unwrap()).unwrap();
        assert_eq!(back, m);
    }

    #[test]
    fn sha256_is_stable_hex() {
        // Known SHA-256 of the empty input.
        assert_eq!(sha256_hex(b""), "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
    }
}
