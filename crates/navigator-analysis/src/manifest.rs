//! The manifest that holds the integrity of each asset. See the ancestry-ibd asset design, which
//! this crosses.
//!
//! `navigator-panelbuild` writes `ancestry_manifest_<build>.json`. That file lists the SHA-256 of
//! each `.bin` that the build made. The app checks a loaded asset against it, and it refuses one
//! that does not match. That guard costs little, and it covers an asset that came over a CDN.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::error::AnalysisError;

/// The SHA-256 of `bytes`, in lower-case hex. It comes from the shared helper in `du-bio`, and
/// this module exports it again, so that a caller of `manifest::sha256_hex` still compiles.
/// `navigator-panelbuild` is such a caller.
pub use du_bio::hash::sha256_hex;

/// One asset's integrity record.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AssetEntry {
    pub sha256: String,
    pub bytes: u64,
}

/// The asset manifest of one build. It maps the file name of an asset, as it sits on disk, to its
/// integrity record.
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

    /// Check `bytes` for `filename`. It gives `Ok` in two cases. The manifest holds no entry for
    /// that file: this check is advisory, and it gates nothing that the manifest does not list. Or
    /// the digest matches. It gives `Err(expected, got)` when the two differ.
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
        // Bytes that match pass the check. Bytes that somebody changed do not. A file that the
        // manifest does not list passes, because the check is advisory.
        assert!(m.verify("ancestry_panel_chm13v2.0.bin", b"hello").is_ok());
        assert!(m.verify("ancestry_panel_chm13v2.0.bin", b"hELLo").is_err());
        assert!(m.verify("not_listed.bin", b"anything").is_ok());
        let back = AssetManifest::from_json(&m.to_json().unwrap()).unwrap();
        assert_eq!(back, m);
    }
}
