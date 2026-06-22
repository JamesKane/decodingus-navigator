//! Local-LLM client (OpenAI-compatible): configuration + resolvers, health/model discovery (M0),
//! and brief narration via chat completions (M1).
//!
//! The entire feature is **local-only** by design (see `docs/design/local-llm-integration.md`):
//! Navigator is a *client* of a model server the user runs (LM Studio / Ollama / llama.cpp). There is
//! no hosted-provider path and no API key. The transport is the OpenAI Chat Completions wire format,
//! the common denominator across local runtimes, spoken over the app's existing `reqwest` client.

use crate::{App, AppError, AppSettings};
use navigator_domain::brief::SubjectBrief;
use navigator_domain::du_domain::ids::SampleGuid;
use navigator_domain::llm_prompt;
use navigator_refgenome::cache as refgenome_cache;
use serde::{Deserialize, Serialize};

/// LM Studio's default OpenAI-compatible base URL — the happy-path local server.
pub const DEFAULT_LLM_BASE_URL: &str = "http://localhost:1234/v1";

/// Default max response tokens. Generous so a reasoning model has room for its full chain-of-thought
/// plus the answer (a small cap is consumed entirely by reasoning and `content` comes back empty).
/// It is a ceiling, not a target — non-reasoning models stop well before it.
pub const DEFAULT_LLM_MAX_TOKENS: u32 = 8192;

/// Resolved local-LLM configuration (env → settings → default, like the other resolvers).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct LlmConfig {
    pub enabled: bool,
    /// Base URL including the OpenAI-compatible path prefix (e.g. `.../v1`), no trailing slash.
    pub base_url: String,
    /// Model id to request, or `None` to let the server use its single loaded model.
    pub model: Option<String>,
    /// Max response (completion) tokens.
    pub max_tokens: u32,
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

fn resolve_max_tokens(env: Option<String>, settings: Option<u32>) -> u32 {
    env.and_then(|v| v.trim().parse::<u32>().ok())
        .or(settings)
        .filter(|n| *n > 0)
        .unwrap_or(DEFAULT_LLM_MAX_TOKENS)
}

/// The configured local-LLM settings, honoring `NAVIGATOR_LLM_*` over the persisted values.
pub fn llm_config() -> LlmConfig {
    let s = AppSettings::load();
    LlmConfig {
        enabled: resolve_enabled(std::env::var("NAVIGATOR_LLM_ENABLED").ok().as_deref(), s.llm_enabled),
        base_url: resolve_base_url(std::env::var("NAVIGATOR_LLM_BASE_URL").ok(), s.llm_base_url),
        model: resolve_model(std::env::var("NAVIGATOR_LLM_MODEL").ok(), s.llm_model),
        max_tokens: resolve_max_tokens(std::env::var("NAVIGATOR_LLM_MAX_TOKENS").ok(), s.llm_max_tokens),
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

// ---- chat completion (M1) ----------------------------------------------------

#[derive(Serialize)]
struct ChatMessage {
    role: String,
    content: String,
}

/// One prior turn of an "ask my results" conversation, carried by the UI and replayed as context.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ChatTurn {
    pub from_user: bool,
    pub text: String,
}

#[derive(Serialize)]
struct ChatRequest {
    model: String,
    messages: Vec<ChatMessage>,
    temperature: f32,
    max_tokens: u32,
    stream: bool,
}

#[derive(Deserialize)]
struct ChatResponse {
    #[serde(default)]
    model: Option<String>,
    #[serde(default)]
    choices: Vec<ChatChoice>,
}

#[derive(Deserialize)]
struct ChatChoice {
    message: ChatChoiceMessage,
    /// `"stop"` | `"length"` | … — `"length"` on a reasoning model means it never reached the answer.
    #[serde(default)]
    finish_reason: Option<String>,
}

#[derive(Deserialize)]
struct ChatChoiceMessage {
    #[serde(default)]
    content: String,
}

/// Strip a reasoning model's chain-of-thought from the visible answer: drop any `<think>…</think>`
/// blocks (some servers inline reasoning in `content` this way) and anything before a closing
/// `</think>`. Returns the trimmed final answer.
fn strip_reasoning(content: &str) -> String {
    // Reasoning models emit a leading `<think>…</think>` block then the answer; keep only what
    // follows the last close tag.
    let after = match content.rfind("</think>") {
        Some(pos) => &content[pos + "</think>".len()..],
        None => content,
    };
    // A remaining opening tag means reasoning was truncated with no answer (it ran out of budget) →
    // drop it so the result is empty and the caller reports the reasoning-ran-out error.
    let answer = match after.find("<think>") {
        Some(open) => &after[..open],
        None => after,
    };
    answer.trim().to_string()
}

/// An AI-assisted narration of a [`SubjectBrief`], with the model that produced it (for labelling).
/// Always rendered *alongside* the structured cards, never instead of them.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NarratedBrief {
    pub prose: String,
    pub model: String,
}

fn narration_cache_path(key: &str) -> std::path::PathBuf {
    refgenome_cache::base_dir()
        .join("briefs")
        .join("narration")
        .join(format!("{key}.json"))
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

    /// Build the subject's brief and narrate it via the local model (M1 entry point used by the UI).
    pub async fn narrate_subject(&self, guid: SampleGuid) -> Result<NarratedBrief, AppError> {
        let brief = self.subject_brief(guid).await?;
        self.narrate_brief(&brief).await
    }

    /// Ask the local model to rewrite the brief's facts as casual-reader prose. Grounded by
    /// [`llm_prompt`] (facts-only, no health, preserve uncertainty); cached on disk keyed by the
    /// fact sheet + model (so changing inputs or model regenerates, and re-opening is free). Returns
    /// `Err(AppError::Llm)` on disabled / unreachable / bad / unsafe output — the UI then keeps the
    /// deterministic brief unchanged. Never the *only* output the user sees.
    pub async fn narrate_brief(&self, brief: &SubjectBrief) -> Result<NarratedBrief, AppError> {
        let cfg = llm_config();
        if !cfg.enabled {
            return Err(AppError::Llm("The AI assistant is turned off.".into()));
        }

        let model = self.resolve_model_id(&cfg).await?;

        let system = llm_prompt::narrate_system_prompt();
        let facts = llm_prompt::narrate_fact_sheet(brief);
        // Key on model + system prompt + facts, so changing the prompt (or facts/model) regenerates
        // rather than serving a stale narration written under the old instructions.
        let key = crate::sha256_str(&format!("{model}\n{system}\n{facts}"));
        let cache_path = narration_cache_path(&key);
        if let Some(cached) = std::fs::read_to_string(&cache_path)
            .ok()
            .and_then(|s| serde_json::from_str::<NarratedBrief>(&s).ok())
        {
            return Ok(cached);
        }

        let messages = vec![
            ChatMessage {
                role: "system".into(),
                content: system,
            },
            ChatMessage {
                role: "user".into(),
                content: facts,
            },
        ];
        let (prose, used_model) = self.chat_complete(&cfg, &model, messages).await?;

        // Post-generation guardrail: reject anything that strays into health/clinical language.
        if llm_prompt::mentions_health(&prose) {
            return Err(AppError::Llm(
                "The AI response was withheld (it strayed outside ancestry).".into(),
            ));
        }

        let result = NarratedBrief {
            prose,
            model: used_model,
        };
        if let Some(parent) = cache_path.parent() {
            let _ = std::fs::create_dir_all(parent);
        }
        if let Ok(json) = serde_json::to_string(&result) {
            let _ = std::fs::write(&cache_path, json);
        }
        Ok(result)
    }

    /// Answer one "ask my results" question, grounded in the subject's brief (M2). The brief's fact
    /// sheet is the only source of facts; chat `history` is replayed for continuity. Off / unreachable
    /// / bad output → `Err`. A clearly health/medical question (in or out) is met with a fixed
    /// ancestry-only deflection instead of a model answer.
    pub async fn answer_question(
        &self,
        guid: SampleGuid,
        history: Vec<ChatTurn>,
        question: String,
    ) -> Result<String, AppError> {
        let cfg = llm_config();
        if !cfg.enabled {
            return Err(AppError::Llm("The AI assistant is turned off.".into()));
        }
        // Incoming scope guard: don't even ask the model a medical question.
        if llm_prompt::mentions_health(&question) {
            return Ok(llm_prompt::health_deflection().to_string());
        }

        let model = self.resolve_model_id(&cfg).await?;
        let brief = self.subject_brief(guid).await?;
        let facts = llm_prompt::narrate_fact_sheet(&brief);
        // One system message carries the guardrails + the results (the only allowed facts).
        let system = format!(
            "{}\n\nRESULTS (your only source of facts):\n{}",
            llm_prompt::answer_system_prompt(),
            facts
        );

        let mut messages = vec![ChatMessage {
            role: "system".into(),
            content: system,
        }];
        // Replay the recent conversation (bounded) for continuity.
        for turn in history.iter().rev().take(12).collect::<Vec<_>>().into_iter().rev() {
            messages.push(ChatMessage {
                role: if turn.from_user { "user" } else { "assistant" }.into(),
                content: turn.text.clone(),
            });
        }
        messages.push(ChatMessage {
            role: "user".into(),
            content: question,
        });

        let (answer, _) = self.chat_complete(&cfg, &model, messages).await?;
        // Outgoing scope guard: a strayed answer is replaced by the deflection, not shown.
        if llm_prompt::mentions_health(&answer) {
            return Ok(llm_prompt::health_deflection().to_string());
        }
        Ok(answer)
    }

    /// Resolve a concrete model id for a request — `cfg.model` if set, else the server's single
    /// loaded model (servers like Ollama require a name).
    async fn resolve_model_id(&self, cfg: &LlmConfig) -> Result<String, AppError> {
        match cfg.model.clone() {
            Some(m) => Ok(m),
            None => self
                .llm_models_at(&cfg.base_url)
                .await?
                .into_iter()
                .next()
                .ok_or_else(|| AppError::Llm("No model is loaded on the server.".into())),
        }
    }

    /// POST a chat completion and return `(answer, model)` with reasoning stripped. Shared by
    /// narration and Q&A. Errors on unreachable server, bad response, or empty/reasoning-truncated
    /// output (callers apply their own health guard).
    async fn chat_complete(
        &self,
        cfg: &LlmConfig,
        model: &str,
        messages: Vec<ChatMessage>,
    ) -> Result<(String, String), AppError> {
        let req = ChatRequest {
            model: model.to_string(),
            messages,
            temperature: 0.4,
            max_tokens: cfg.max_tokens,
            stream: false,
        };
        let url = format!("{}/chat/completions", cfg.base_url.trim_end_matches('/'));
        let resp = self
            .auth
            .http
            .post(&url)
            .json(&req)
            // Generous — a reasoning model on a large local model can take minutes to think + answer.
            .timeout(std::time::Duration::from_secs(300))
            .send()
            .await
            .map_err(|e| AppError::Llm(format!("Could not reach the local model server: {e}")))?;
        let resp = resp
            .error_for_status()
            .map_err(|e| AppError::Llm(format!("The model server returned an error: {e}")))?;
        let parsed: ChatResponse = resp
            .json()
            .await
            .map_err(|e| AppError::Llm(format!("Unexpected response from the server: {e}")))?;

        let choice = parsed.choices.first();
        let answer = choice.map(|c| strip_reasoning(&c.message.content)).unwrap_or_default();
        if answer.is_empty() {
            let truncated = choice.and_then(|c| c.finish_reason.as_deref()) == Some("length");
            return Err(AppError::Llm(if truncated {
                "The model ran out of room before answering — it spent its budget reasoning. Raise \
                 'Max response tokens' in Settings, or use a smaller / non-reasoning model."
                    .into()
            } else {
                "The model returned an empty response.".into()
            }));
        }
        Ok((answer, parsed.model.unwrap_or_else(|| model.to_string())))
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
    fn strip_reasoning_keeps_final_answer() {
        assert_eq!(strip_reasoning("<think>pondering…</think>The answer."), "The answer.");
        assert_eq!(strip_reasoning("  plain answer  "), "plain answer");
        // Reasoning closed, then a stray re-open with no answer.
        assert_eq!(strip_reasoning("<think>a</think>answer <think>b"), "answer");
        // Unterminated reasoning (hit the cap) → nothing usable.
        assert_eq!(strip_reasoning("<think>still thinking"), "");
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
