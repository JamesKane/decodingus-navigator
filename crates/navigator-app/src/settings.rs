//! Persisted application settings at `~/.decodingus/config/settings.json`.
//!
//! These are consulted by the resolvers in [`crate`] **below** any environment variable (env wins →
//! settings → built-in default), so the Settings UI can change app behavior — AppView URL, Y-tree
//! provider, tree-cache TTL, theme — without env vars or a relaunch. The file is small; resolvers
//! re-read it per call (they run per-analysis, not in a hot loop), so edits apply immediately.

use serde::{Deserialize, Serialize};
use std::path::PathBuf;

#[derive(Debug, Clone, Default, PartialEq, Deserialize, Serialize)]
pub struct AppSettings {
    /// Y-tree provider: `"decodingus"` or `"ftdna"`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub y_tree_provider: Option<String>,
    /// AppView base URL (tree API + sequencer-lab lookup).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub appview_url: Option<String>,
    /// Haplotree cache TTL in days (`0` = always refetch).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tree_ttl_days: Option<u64>,
    /// UI theme: `"dark"` or `"light"`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub theme: Option<String>,
    /// Ask before downloading large reference files.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub prompt_before_download: Option<bool>,
    /// UI scale (egui zoom factor) — raise it on a native-4K / HiDPI display where the OS reports a
    /// 1.0 scale factor and the default text is tiny. `None` = 1.0.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ui_scale: Option<f32>,
}

impl AppSettings {
    /// `~/.decodingus/config/settings.json` (honoring `NAVIGATOR_REFGENOME_DIR`, same base as the
    /// reference-source overrides).
    pub fn path() -> PathBuf {
        navigator_refgenome::cache::base_dir().join("config").join("settings.json")
    }

    /// Load settings; a missing or unreadable/invalid file yields the empty default.
    pub fn load() -> Self {
        std::fs::read_to_string(Self::path())
            .ok()
            .and_then(|s| serde_json::from_str(&s).ok())
            .unwrap_or_default()
    }

    /// The `~/.decodingus` base directory (honoring `NAVIGATOR_REFGENOME_DIR`).
    pub fn cache_base_dir() -> PathBuf {
        navigator_refgenome::cache::base_dir()
    }

    /// Persist to disk (creating the `config/` dir), pretty-printed.
    pub fn save(&self) -> std::io::Result<()> {
        let path = Self::path();
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent)?;
        }
        let json = serde_json::to_string_pretty(self).map_err(std::io::Error::other)?;
        std::fs::write(path, json)
    }
}
