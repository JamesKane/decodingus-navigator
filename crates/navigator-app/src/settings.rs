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
    /// Prefer a trusted external caller (GATK4 GVCF / 1240K call set, imported via the sidecar fast
    /// path) over Navigator's own genotyping. When on (the built-in default), an external Y/mt/
    /// autosomal call wins reconciliation and Navigator's internal caller does not re-walk that
    /// alignment. `None` = the default (on); set `Some(false)` to always run the internal caller.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub prefer_external_calls: Option<bool>,
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
    /// Interface mode: `"simple"` (casual single-person briefs) or `"advanced"` (full power-user UI).
    /// `None` = the user has never pinned a mode, so the UI applies its first-run heuristic.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ui_mode: Option<String>,
    /// Enable local-LLM assisted narration / chat. Off until the user opts in.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub llm_enabled: Option<bool>,
    /// OpenAI-compatible base URL of the *local* model server, e.g. "http://localhost:1234/v1".
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub llm_base_url: Option<String>,
    /// Model id to request (as reported by `GET /models`), e.g. "llama-3.1-8b-instruct".
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub llm_model: Option<String>,
    /// Max response (completion) tokens to request. Reasoning models spend most of this on their
    /// chain-of-thought, so it must be large enough for the thinking *and* the answer. `None` = default.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub llm_max_tokens: Option<u32>,
    /// Check GitHub Releases for a newer installer at startup and notify. `None` = the built-in
    /// default (enabled); set `Some(false)` to opt out. Never auto-installs.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub check_for_updates: Option<bool>,
    /// A version the user asked not to be reminded about (the exact `latest_version` string). A
    /// *newer* release than this still notifies.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub skip_update_version: Option<String>,
}

impl AppSettings {
    /// `~/.decodingus/config/settings.json` (honoring `NAVIGATOR_REFGENOME_DIR`, same base as the
    /// reference-source overrides).
    pub fn path() -> PathBuf {
        navigator_refgenome::cache::base_dir()
            .join("config")
            .join("settings.json")
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
