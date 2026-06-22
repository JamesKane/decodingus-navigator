//! Local-LLM client (OpenAI-compatible) — **M0**: configuration, resolvers, and health/model
//! discovery only. No generation yet (that is M1).
//!
//! The entire feature is **local-only** by design (see `docs/design/local-llm-integration.md`):
//! Navigator is a *client* of a model server the user runs (LM Studio / Ollama / llama.cpp). There is
//! no hosted-provider path and no API key. The transport is the OpenAI Chat Completions wire format,
//! the common denominator across local runtimes, spoken over the app's existing `reqwest` client.

use crate::{App, AppError, AppSettings};
use serde::Deserialize;

/// LM Studio's default OpenAI-compatible base URL — the happy-path local server.
pub const DEFAULT_LLM_BASE_URL: &str = "http://localhost:1234/v1";

/// Resolved local-LLM configuration (env → settings → default, like the other resolvers).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LlmConfig {
    pub enabled: bool,
    /// Base URL including the OpenAI-compatible path prefix (e.g. `.../v1`), no trailing slash.
    pub base_url: String,
    /// Model id to request, or `None` to let the server use its single loaded model.
    pub model: Option<String>,
}

fn resolve_enabled(env: Option<&str>, settings: Option<bool>) -> bool {
    match env.map(str::trim) {
        Some(v) => matches!(v.to_ascii_lowercase().as_str(), "1" | "true" | "yes" | "on"),
        None => settings.unwrap_or(false),
    }
}

fn resolve_base_url(env: Option<String>, settings: Option<String>) -> String {
    env.or(settings)
        .map(|s| s.trim().trim_end_matches('/').to_string())
        .filter(|s| !s.is_empty())
        .unwrap_or_else(|| DEFAULT_LLM_BASE_URL.to_string())
}

fn resolve_model(env: Option<String>, settings: Option<String>) -> Option<String> {
    env.or(settings).map(|s| s.trim().to_string()).filter(|s| !s.is_empty())
}

/// The configured local-LLM settings, honoring `NAVIGATOR_LLM_*` over the persisted values.
pub fn llm_config() -> LlmConfig {
    let s = AppSettings::load();
    LlmConfig {
        enabled: resolve_enabled(std::env::var("NAVIGATOR_LLM_ENABLED").ok().as_deref(), s.llm_enabled),
        base_url: resolve_base_url(std::env::var("NAVIGATOR_LLM_BASE_URL").ok(), s.llm_base_url),
        model: resolve_model(std::env::var("NAVIGATOR_LLM_MODEL").ok(), s.llm_model),
    }
}

/// Is `base_url`'s host a loopback address? Drives the Settings warning when a user points the client
/// at a non-local server (results would leave the machine). Conservative: anything we can't confirm
/// is loopback is treated as remote.
pub fn is_loopback_url(base_url: &str) -> bool {
    let after_scheme = base_url.split_once("://").map(|(_, r)| r).unwrap_or(base_url);
    let authority = after_scheme.split(['/', '?', '#']).next().unwrap_or("");
    // Strip an IPv6 bracket or a trailing :port to get the bare host.
    let host = if let Some(rest) = authority.strip_prefix('[') {
        rest.split(']').next().unwrap_or("")
    } else {
        authority.rsplit_once(':').map(|(h, _)| h).unwrap_or(authority)
    };
    host.eq_ignore_ascii_case("localhost") || host == "127.0.0.1" || host == "::1"
}

#[derive(Deserialize)]
struct ModelsResponse {
    #[serde(default)]
    data: Vec<ModelEntry>,
}

#[derive(Deserialize)]
struct ModelEntry {
    id: String,
}

impl App {
    /// Health check + model discovery against the configured local server.
    pub async fn llm_models(&self) -> Result<Vec<String>, AppError> {
        let cfg = llm_config();
        self.llm_models_at(&cfg.base_url).await
    }

    /// Health check + model discovery against an explicit base URL — used by the Settings
    /// "Test connection" button so the user can verify a URL *before* saving it. `GET {base}/models`
    /// (the OpenAI-compatible discovery endpoint). Errors are plain-language for the UI.
    pub async fn llm_models_at(&self, base_url: &str) -> Result<Vec<String>, AppError> {
        let base = base_url.trim().trim_end_matches('/');
        if base.is_empty() {
            return Err(AppError::Llm("No server URL set.".into()));
        }
        let url = format!("{base}/models");
        let resp = self
            .auth
            .http
            .get(&url)
            .timeout(std::time::Duration::from_secs(5))
            .send()
            .await
            .map_err(|_| {
                AppError::Llm(format!(
                    "No local model server found at {base}. Start LM Studio's server (or Ollama) and try again."
                ))
            })?;
        let resp = resp
            .error_for_status()
            .map_err(|e| AppError::Llm(format!("The model server returned an error: {e}")))?;
        let body = resp.text().await.map_err(|e| AppError::Llm(e.to_string()))?;
        let parsed: ModelsResponse = serde_json::from_str(&body)
            .map_err(|e| AppError::Llm(format!("Unexpected response from the server: {e}")))?;
        Ok(parsed.data.into_iter().map(|m| m.id).collect())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn enabled_resolution_env_over_settings() {
        assert!(resolve_enabled(Some("true"), Some(false)));
        assert!(resolve_enabled(Some("1"), None));
        assert!(!resolve_enabled(Some("off"), Some(true)));
        assert!(resolve_enabled(None, Some(true)));
        assert!(!resolve_enabled(None, None));
    }

    #[test]
    fn base_url_normalization() {
        assert_eq!(
            resolve_base_url(Some("http://host:1234/v1/".into()), None),
            "http://host:1234/v1"
        );
        assert_eq!(resolve_base_url(None, Some("  ".into())), DEFAULT_LLM_BASE_URL);
        assert_eq!(resolve_base_url(None, None), DEFAULT_LLM_BASE_URL);
        // env wins over settings
        assert_eq!(
            resolve_base_url(Some("http://a/v1".into()), Some("http://b/v1".into())),
            "http://a/v1"
        );
    }

    #[test]
    fn loopback_detection() {
        assert!(is_loopback_url("http://localhost:1234/v1"));
        assert!(is_loopback_url("http://127.0.0.1:11434/v1"));
        assert!(is_loopback_url("http://[::1]:8080/v1"));
        assert!(is_loopback_url("http://LocalHost:1234"));
        assert!(!is_loopback_url("http://192.168.1.50:1234/v1"));
        assert!(!is_loopback_url("https://api.example.com/v1"));
    }
}
