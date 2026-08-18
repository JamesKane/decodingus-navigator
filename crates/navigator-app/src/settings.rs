//! Persisted application settings at `~/.decodingus/config/settings.json`.
//!
//! The resolvers in [`crate`] read these settings. An environment variable has priority over a
//! setting, and a setting has priority over the built-in default.
//!
//! So the Settings UI can change the behaviour of the app with no environment variable and no
//! restart. The user can change the AppView URL, the Y-tree provider, the TTL of the tree cache,
//! and the theme.
//!
//! The file is small, and a resolver reads it again at each call. A resolver runs once for each
//! analysis, not in a loop, so a change applies immediately.

use serde::{Deserialize, Serialize};
use std::path::PathBuf;

#[derive(Debug, Clone, Default, PartialEq, Deserialize, Serialize)]
pub struct AppSettings {
    /// Y-tree provider: `"decodingus"` or `"ftdna"`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub y_tree_provider: Option<String>,
    /// Use a trusted external caller before the internal caller of Navigator. The external caller
    /// is a GATK4 GVCF or a 1240K call set, and the sidecar fast path imports it.
    ///
    /// When this setting is on, which is the built-in default, an external Y call, mt call, or
    /// autosomal call wins the reconciliation. The internal caller does not walk that alignment
    /// again. `None` means the default, which is on. `Some(false)` makes the internal caller always
    /// run.
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
    /// Ask the user before a download of a large reference file.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub prompt_before_download: Option<bool>,
    /// The scale of the UI, which is the zoom factor of egui. Increase it on a 4K display or a
    /// HiDPI display when the operating system reports a scale factor of 1.0 and the text is very
    /// small. `None` means 1.0.
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
    /// The maximum count of response tokens to request. A model that reasons uses most of these
    /// tokens for its internal steps. So the value must be large enough for those steps and for the
    /// answer. `None` means the default.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub llm_max_tokens: Option<u32>,
    /// Check GitHub Releases for a newer installer at startup and notify. `None` = the built-in
    /// default (enabled); set `Some(false)` to opt out. Never auto-installs.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub check_for_updates: Option<bool>,
    /// A version that the user does not want a reminder about. The value is the exact
    /// `latest_version` string. The app still gives a notification for a newer release.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub skip_update_version: Option<String>,
    /// The last inner size of the window as `[width, height]` in egui points. The app keeps this
    /// value between launches. It is `None` until the first run writes it.
    ///
    /// At the next start, the app fits the size to the current monitor. A stored size can be too
    /// large, for example from a larger display, and the app then makes it smaller.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub window_size: Option<[f32; 2]>,
    /// The navigation view that the user selected last. The value is `"dashboard"`, `"subjects"`,
    /// `"projects"`, or `"community"`. The app opens this view at the next start.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub last_nav: Option<String>,
    /// Last focused subject (biosample GUID), restored once the subject list has loaded.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub last_subject: Option<String>,
    /// Last selected subject detail tab (`"overview"` / `"ydna"` / `"mtdna"` / `"autosomal"` /
    /// `"ancestry"` / `"sources"` / `"ibd"`).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub last_detail_tab: Option<String>,

    // ── Calibration values for the chromosome painter (copy-LAI). `None` = the built-in default ──
    /// The switch intensity of a reference haplotype for each cM. This value models recombination
    /// in the copy step. A lower value gives longer copied tracts. Longer tracts give a cleaner
    /// population call, and they attract a drifted isolate less. The default is 0.1.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lai_recomb_per_cm: Option<f64>,
    /// The maximum count of reference haplotypes for each population. This limit balances the
    /// panel, so a large 1000G sample does not win only by its count. The default is 50.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lai_max_ref_haps: Option<u32>,
    /// Global-composition gate: drop super-populations below this genome-wide fraction. Default 0.05.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lai_min_ancestry: Option<f64>,
    /// The switch intensity of an ancestry segment for each cM, which the Viterbi step uses to
    /// make the track smooth. A lower value gives longer segments. The default is 0.05.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lai_switch_per_cm: Option<f64>,
    /// The minimum length of a segment in centiMorgans. The code joins a shorter run to its
    /// neighbour. The default is 4.0.
    ///
    /// The unit is genetic distance, not a count of sites. So the value stays correct when the
    /// marker density of the panel changes.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lai_min_segment_cm: Option<f64>,
    /// The exponent that corrects for the size of each population. A value of 0 turns the
    /// correction off. A value of 1 gives the full average for each haplotype. The default is
    /// 0.5.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lai_size_normalize: Option<f64>,
    /// Copy mismatch/mutation rate μ. Default 0.02.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lai_mismatch: Option<f64>,
}

impl AppSettings {
    /// `~/.decodingus/config/settings.json`. The path obeys `NAVIGATOR_REFGENOME_DIR` and uses the
    /// same base as the
    /// reference-source overrides).
    pub fn path() -> PathBuf {
        navigator_refgenome::cache::base_dir()
            .join("config")
            .join("settings.json")
    }

    /// Read the settings. If the file is absent, unreadable, or invalid, the function gives the
    /// empty default.
    ///
    /// The read goes through [`navigator_refgenome::cache::read_atomic`]. So a [`Self::save`] call
    /// at the same time can not make the read fail. On Windows, no process can open the file for a
    /// short time during a replace.
    ///
    /// Without this protection, the function returns the default and the app removes the settings
    /// of the user without a message. The painter calibration is one of those settings, because
    /// `copying_lai_params()` reads the settings at paint time.
    pub fn load() -> Self {
        navigator_refgenome::cache::read_atomic(&Self::path())
            .ok()
            .and_then(|b| String::from_utf8(b).ok())
            .and_then(|s| serde_json::from_str(&s).ok())
            .unwrap_or_default()
    }

    /// The `~/.decodingus` base directory. The path obeys `NAVIGATOR_REFGENOME_DIR`.
    pub fn cache_base_dir() -> PathBuf {
        navigator_refgenome::cache::base_dir()
    }

    /// Write the settings to disk in a readable layout, and make the `config/` directory if it
    /// does not exist.
    ///
    /// The write is **atomic**. The function writes a temporary file and then renames it. Some UI
    /// paths save the settings at the same time. A write that is not atomic can give a torn file,
    /// as it did with `reference_sources.json`.
    pub fn save(&self) -> std::io::Result<()> {
        let json = serde_json::to_string_pretty(self).map_err(std::io::Error::other)?;
        navigator_refgenome::cache::atomic_write(&Self::path(), json.as_bytes())
    }
}
