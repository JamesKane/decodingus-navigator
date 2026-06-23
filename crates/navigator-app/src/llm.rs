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
use navigator_domain::results_context::{
    IbdFact, IbdMatchFact, MtMutationsFact, PrivateYFact, ResultsContext, SexFact, YStrPanelFact,
};
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
    /// llama.cpp/LM Studio Jinja-template kwargs. We pass `{"enable_thinking": false}` so reasoning
    /// models (Gemma 4, Qwen, DeepSeek-R1) never emit a thinking channel — saving the tokens/latency
    /// `strip_reasoning` would otherwise discard. Skipped from the body when unset.
    #[serde(skip_serializing_if = "Option::is_none")]
    chat_template_kwargs: Option<serde_json::Value>,
}

/// One parsed SSE `data:` chunk from a streamed chat completion.
struct StreamDelta {
    content: Option<String>,
    /// `"stop"` | `"length"` | … — `"length"` on a reasoning model means it never reached the answer.
    finish_reason: Option<String>,
    model: Option<String>,
}

/// Parse one SSE payload (already stripped of the `data:` prefix). `None` for `[DONE]` / unparseable.
fn parse_stream_event(data: &str) -> Option<StreamDelta> {
    if data == "[DONE]" {
        return None;
    }
    let v: serde_json::Value = serde_json::from_str(data).ok()?;
    let choice = v.get("choices").and_then(|c| c.get(0));
    Some(StreamDelta {
        content: choice
            .and_then(|c| c.get("delta"))
            .and_then(|d| d.get("content"))
            .and_then(|x| x.as_str())
            .map(String::from),
        finish_reason: choice
            .and_then(|c| c.get("finish_reason"))
            .and_then(|x| x.as_str())
            .map(String::from),
        model: v.get("model").and_then(|x| x.as_str()).map(String::from),
    })
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

fn narration_cache_path_in(subdir: &str, key: &str) -> std::path::PathBuf {
    refgenome_cache::base_dir()
        .join("briefs")
        .join(subdir)
        .join(format!("{key}.json"))
}

fn narration_cache_path(key: &str) -> std::path::PathBuf {
    narration_cache_path_in("narration", key)
}

impl App {
    /// Health check + model discovery against the configured local server.
    pub async fn llm_models(&self) -> Result<Vec<String>, AppError> {
        let cfg = llm_config();
        self.llm_models_at(&cfg.base_url).await
    }

    /// The on-disk cached narration for a brief, if one exists for the currently-configured model —
    /// a **no-network** lookup (only when a model is explicitly set) used to fold the AI story into
    /// the exported "DNA Story" without triggering generation.
    pub fn cached_narration(&self, brief: &SubjectBrief) -> Option<NarratedBrief> {
        let cfg = llm_config();
        let model = cfg.model?; // explicit only — avoid resolving the loaded model (a network call)
        let system = llm_prompt::narrate_system_prompt();
        let facts = llm_prompt::narrate_fact_sheet(brief);
        let key = crate::sha256_str(&format!("{model}\n{system}\n{facts}"));
        std::fs::read_to_string(narration_cache_path(&key))
            .ok()
            .and_then(|s| serde_json::from_str(&s).ok())
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
        self.narrate_brief_streaming(brief, |_| {}).await
    }

    /// Build the subject's brief and narrate it, streaming visible answer text to `on_chunk` as it
    /// arrives (M3). The final return value is authoritative; chunks are a live preview.
    pub async fn narrate_subject_streaming(
        &self,
        guid: SampleGuid,
        on_chunk: impl FnMut(&str),
    ) -> Result<NarratedBrief, AppError> {
        let brief = self.subject_brief(guid).await?;
        self.narrate_brief_streaming(&brief, on_chunk).await
    }

    /// Streaming core for narration: cache-first (a hit replays the cached prose through `on_chunk`),
    /// else stream a completion, health-guard, and cache.
    pub async fn narrate_brief_streaming(
        &self,
        brief: &SubjectBrief,
        on_chunk: impl FnMut(&str),
    ) -> Result<NarratedBrief, AppError> {
        let cfg = llm_config();
        if !cfg.enabled {
            return Err(AppError::Llm("The AI assistant is turned off.".into()));
        }

        let model = self.resolve_model_id(&cfg).await?;
        let system = llm_prompt::narrate_system_prompt();
        let facts = llm_prompt::narrate_fact_sheet(brief);
        self.run_cached_narration(&cfg, &model, system, facts, "narration", on_chunk)
            .await
    }

    /// Explain a single result signal (M5 per-tab "Explain this") in plain language, streaming the
    /// prose to `on_chunk`. Grounded in only that signal's curated section (see
    /// [`navigator_domain::results_context::signal_section`]); cached and health-guarded like brief
    /// narration. `Err` when the assistant is off / unreachable / the subject has nothing for that
    /// signal — the UI then just doesn't show an explanation.
    pub async fn narrate_signal_streaming(
        &self,
        guid: SampleGuid,
        kind: navigator_domain::results_context::SignalKind,
        on_chunk: impl FnMut(&str),
    ) -> Result<NarratedBrief, AppError> {
        use navigator_domain::results_context as rc;
        let cfg = llm_config();
        if !cfg.enabled {
            return Err(AppError::Llm("The AI assistant is turned off.".into()));
        }
        let model = self.resolve_model_id(&cfg).await?;
        let ctx = self.results_context(guid).await?;
        let section = rc::signal_section(&ctx, kind)
            .ok_or_else(|| AppError::Llm("There's nothing to explain here yet.".into()))?;
        let system = llm_prompt::narrate_signal_system_prompt(kind.label());
        self.run_cached_narration(&cfg, &model, system, section, "signals", on_chunk)
            .await
    }

    /// Shared cached-streaming narration core for [`narrate_brief_streaming`] and
    /// [`narrate_signal_streaming`]: cache-first (a hit replays the cached prose through `on_chunk`),
    /// else stream a completion, reject health-straying output, and cache. The cache key is
    /// `model + system + facts`, so changing the prompt, the facts, or the model regenerates rather
    /// than serving a stale narration written under the old instructions. `subdir` separates brief
    /// vs per-signal caches under `briefs/`.
    async fn run_cached_narration(
        &self,
        cfg: &LlmConfig,
        model: &str,
        system: String,
        facts: String,
        subdir: &str,
        mut on_chunk: impl FnMut(&str),
    ) -> Result<NarratedBrief, AppError> {
        let key = crate::sha256_str(&format!("{model}\n{system}\n{facts}"));
        let cache_path = narration_cache_path_in(subdir, &key);
        if let Some(cached) = std::fs::read_to_string(&cache_path)
            .ok()
            .and_then(|s| serde_json::from_str::<NarratedBrief>(&s).ok())
        {
            on_chunk(&cached.prose);
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
        let (prose, used_model) = self
            .chat_complete_streaming(cfg, model, messages, &mut on_chunk)
            .await?;

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
        self.answer_question_streaming(guid, history, question, |_| {}).await
    }

    /// Streaming variant of [`answer_question`]: visible answer text is sent to `on_chunk` as it
    /// arrives (the final return value is authoritative). The scope-guard deflections do not stream.
    pub async fn answer_question_streaming(
        &self,
        guid: SampleGuid,
        history: Vec<ChatTurn>,
        question: String,
        mut on_chunk: impl FnMut(&str),
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
        let ctx = self.results_context(guid).await?;
        let facts = navigator_domain::results_context::results_fact_sheet(&ctx);
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

        let (answer, _) = self
            .chat_complete_streaming(&cfg, &model, messages, &mut on_chunk)
            .await?;
        // Outgoing scope guard: a strayed answer is replaced by the deflection, not shown.
        if llm_prompt::mentions_health(&answer) {
            return Ok(llm_prompt::health_deflection().to_string());
        }
        Ok(answer)
    }

    /// Assemble the broader grounding context for the M4 chat: the subject brief plus curated,
    /// summary-level facts for the other signals (genetic sex, Y-STR panels, private-Y variants,
    /// mtDNA mutations, IBD matches). Every signal is best-effort — a missing or un-run one is simply
    /// omitted, so the chat grounds in whatever the subject actually has without erroring out.
    pub async fn results_context(&self, guid: SampleGuid) -> Result<ResultsContext, AppError> {
        use navigator_analysis::mtvariants::MtRegion;
        use navigator_analysis::sex::{Confidence, InferredSex};

        let brief = self.subject_brief(guid).await?;

        // Genetic sex (needs an alignment; only a definite call is grounded).
        let sex = match self.default_alignment_for_subject(guid).await? {
            Some((_, aln_id)) => self.cached_sex(aln_id).await?.and_then(|r| {
                let label = match r.inferred_sex {
                    InferredSex::Male => "Male (XY)",
                    InferredSex::Female => "Female (XX)",
                    InferredSex::Unknown => return None,
                };
                let confidence = match r.confidence {
                    Confidence::High => "high",
                    Confidence::Medium => "medium",
                    Confidence::Low => "low",
                };
                Some(SexFact {
                    label: label.into(),
                    confidence: confidence.into(),
                })
            }),
            None => None,
        };

        // Y-STR panels — name + marker count only (never the raw values: token cost, no answerable
        // gain, and they're lineage patterns not facts to recite).
        let ystr: Vec<YStrPanelFact> = self
            .list_str_profiles(guid)
            .await
            .unwrap_or_default()
            .into_iter()
            .map(|p| YStrPanelFact {
                panel: p.panel_name,
                markers: p.markers.len(),
            })
            .collect();

        // Private Y variants — the PrivateBucket confidence split.
        let private_y = self.donor_private_y(guid).await?.map(|b| PrivateYFact {
            novel_unique: b.novel_in_unique_sequence(),
            off_path: b.off_path(),
            structural: b.in_structural_region(),
        });

        // mtDNA mutations relative to rCRS, bucketed by region with a few example notations.
        let mt_mutations = match self.list_mtdna_sequences(guid).await?.first() {
            Some(seq) => {
                let vars = self.mtdna_variants(seq.id).await?;
                if vars.is_empty() {
                    None
                } else {
                    let (mut hvr1, mut hvr2, mut coding) = (0usize, 0usize, 0usize);
                    for v in &vars {
                        match v.region() {
                            MtRegion::Hvr1 => hvr1 += 1,
                            MtRegion::Hvr2 => hvr2 += 1,
                            MtRegion::Coding => coding += 1,
                        }
                    }
                    let examples = vars.iter().take(8).map(|v| v.notation()).collect();
                    Some(MtMutationsFact {
                        total: vars.len(),
                        hvr1,
                        hvr2,
                        coding,
                        examples,
                    })
                }
            }
            None => None,
        };

        // IBD matches (completed exchanges): count + closest by shared cM, partner identity withheld.
        let exchanges = self.list_ibd_exchanges_for_subject(guid).await.unwrap_or_default();
        let ibd = if exchanges.is_empty() {
            None
        } else {
            let closest = exchanges
                .iter()
                .max_by(|a, b| {
                    a.total_shared_cm
                        .partial_cmp(&b.total_shared_cm)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
                .map(|e| IbdMatchFact {
                    relationship: e.relationship.clone(),
                    total_shared_cm: e.total_shared_cm,
                    segment_count: e.segment_count,
                });
            Some(IbdFact {
                match_count: exchanges.len(),
                closest,
            })
        };

        Ok(ResultsContext {
            brief,
            sex,
            ystr,
            private_y,
            mt_mutations,
            ibd,
        })
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

    /// Stream a chat completion (`stream: true`), forwarding visible answer text to `on_chunk` as it
    /// arrives and returning the final `(answer, model)` with reasoning stripped. Reasoning is
    /// suppressed live: only the post-`</think>` answer streams (so a thinking model shows nothing
    /// until it starts answering). Uses `Response::chunk()` (no extra dependency). Shared by
    /// narration and Q&A; callers apply their own health guard.
    async fn chat_complete_streaming(
        &self,
        cfg: &LlmConfig,
        model: &str,
        messages: Vec<ChatMessage>,
        mut on_chunk: impl FnMut(&str),
    ) -> Result<(String, String), AppError> {
        let req = ChatRequest {
            model: model.to_string(),
            messages,
            temperature: 0.4,
            max_tokens: cfg.max_tokens,
            stream: true,
            // Grounded "explain my results" never needs chain-of-thought — disable it at the server
            // so Gemma 4 et al. don't waste tokens/latency on a reasoning channel we'd discard.
            chat_template_kwargs: Some(serde_json::json!({ "enable_thinking": false })),
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
        let mut resp = resp
            .error_for_status()
            .map_err(|e| AppError::Llm(format!("The model server returned an error: {e}")))?;

        let mut buf: Vec<u8> = Vec::new();
        let mut full = String::new();
        let mut prev_visible = String::new();
        let mut finish_reason: Option<String> = None;
        let mut used_model: Option<String> = None;

        while let Some(bytes) = resp
            .chunk()
            .await
            .map_err(|e| AppError::Llm(format!("The connection to the model was interrupted: {e}")))?
        {
            buf.extend_from_slice(&bytes);
            // Process every complete `\n`-terminated SSE line (decoded only once whole, so a
            // multibyte char split across network chunks stays intact).
            while let Some(pos) = buf.iter().position(|&b| b == b'\n') {
                let line_bytes: Vec<u8> = buf.drain(..=pos).collect();
                let line = String::from_utf8_lossy(&line_bytes);
                let Some(data) = line.trim().strip_prefix("data:") else {
                    continue;
                };
                let Some(ev) = parse_stream_event(data.trim()) else {
                    continue;
                };
                if used_model.is_none() {
                    used_model = ev.model;
                }
                if ev.finish_reason.is_some() {
                    finish_reason = ev.finish_reason;
                }
                if let Some(c) = ev.content {
                    full.push_str(&c);
                    let visible = strip_reasoning(&full);
                    if visible.len() > prev_visible.len() && visible.starts_with(&prev_visible) {
                        on_chunk(&visible[prev_visible.len()..]);
                        prev_visible = visible;
                    }
                }
            }
        }

        let answer = strip_reasoning(&full);
        if answer.is_empty() {
            let truncated = finish_reason.as_deref() == Some("length");
            return Err(AppError::Llm(if truncated {
                "The model ran out of room before answering — it spent its budget reasoning. Raise \
                 'Max response tokens' in Settings, or use a smaller / non-reasoning model."
                    .into()
            } else {
                "The model returned an empty response.".into()
            }));
        }
        Ok((answer, used_model.unwrap_or_else(|| model.to_string())))
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
    fn stream_event_parsing() {
        let e = parse_stream_event(r#"{"model":"m","choices":[{"delta":{"content":"Hi"}}]}"#).unwrap();
        assert_eq!(e.content.as_deref(), Some("Hi"));
        assert_eq!(e.model.as_deref(), Some("m"));
        let e = parse_stream_event(r#"{"choices":[{"delta":{},"finish_reason":"length"}]}"#).unwrap();
        assert_eq!(e.content, None);
        assert_eq!(e.finish_reason.as_deref(), Some("length"));
        assert!(parse_stream_event("[DONE]").is_none());
        assert!(parse_stream_event("not json").is_none());
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
