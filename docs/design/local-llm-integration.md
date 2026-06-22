# Local LLM integration (LM Studio / Ollama) — design

Status: **design** (2026-06-22). Implementation not started. Pilot.

## Goal

Let DUNavigator optionally use a **locally-running large language model** to make results easier to
understand — first by narrating the [Subject Brief](subject-brief-simple-mode.md) in warmer,
casual-reader prose, later as an "ask my results" chat. The model runs on the user's own machine via
a host they already control (**LM Studio**, **Ollama**, `llama.cpp`'s server, or any
OpenAI-compatible endpoint). Navigator talks to it over `http://localhost`.

This is scoped as a **pilot**: prove the plumbing and the guardrails on one tightly-grounded use case
(brief narration, M1) before opening up the more open-ended chat (M2).

### Why local-only — and why that is the whole point

DUNavigator processes **genetic data**, the most sensitive personal data there is. The entire app is
built around *local* analysis (no JVM, no external tools, results in a local SQLite workspace). An LLM
feature has to honor that posture:

- **No cloud LLM fallback, ever.** If the configured local endpoint is unreachable, the feature
  degrades to the existing deterministic output — it does **not** silently send a prompt containing a
  person's haplogroups/ancestry to a hosted API. There is no API-key field and no hosted provider in
  this design. (A future opt-in hosted mode would be a separate, explicit, loudly-consented design.)
- **The user owns the model.** Navigator is a *client* of a model server the user started. We don't
  bundle, download, or manage model weights — that is LM Studio / Ollama's job.
- **Off by default.** The feature is disabled until the user enables it in Settings and points it at a
  running endpoint. Nothing changes for users who never touch it.

This makes "local LLM" a genuine product differentiator, not just a technical choice: *AI-assisted
interpretation of your genome that never leaves your computer.*

---

## What we talk to: the OpenAI-compatible surface

We standardize on the **OpenAI Chat Completions** wire format as the common denominator. It is what
LM Studio serves, what Ollama serves at `/v1`, and what `llama.cpp`'s server and most local runtimes
expose. Targeting one schema means one client and no per-vendor branches.

| Host        | Default base URL              | Models endpoint        | Notes                                |
|-------------|-------------------------------|------------------------|--------------------------------------|
| LM Studio   | `http://localhost:1234/v1`    | `GET /models`          | Local server tab; OpenAI-compatible. |
| Ollama      | `http://localhost:11434/v1`   | `GET /models`          | Native API also at `:11434/api/*`.   |
| llama.cpp   | `http://localhost:8080/v1`    | `GET /models`          | `server` binary.                     |

Endpoints used:

- `GET  {base}/models` — health check + model discovery (populate the model dropdown, confirm the
  server is up).
- `POST {base}/chat/completions` — the actual generation. Request body: `model`, `messages`
  (`system`/`user`), `temperature`, `max_tokens`, and `stream: true` for token streaming. Response is
  SSE chunks (`data: {json}` lines, terminated by `data: [DONE]`).

We do **not** depend on any vendor SDK — the bodies are small, hand-rolled serde structs over the
existing `reqwest` client. No new heavyweight dependency for the pilot.

---

## Architecture & layering

Crate rule (`ui → app → {analysis, store, sync, refgenome} → {domain, du-*}`) is respected. The LLM
is network I/O that composes app-level signals, so it lives in **`navigator-app`** alongside the brief
composition that already does exactly this kind of "pull signals + call out + assemble" work.

```
   navigator-domain::llm_prompt   (pure: SubjectBrief + facts -> prompt strings; grounding rules)
            ▲
            │  builds messages
   navigator-app::llm             (the client)
       ├─ LlmConfig (resolved from settings/env)
       ├─ LlmClient: health()/models()  + complete()/complete_streaming()
       └─ narrate_brief(&SubjectBrief)  ← M1 entry point
       └─ answer_question(ctx, q)       ← M2 entry point
            ▲
            │  Command/Event (worker thread, like LoadSubjectBrief)
   navigator-ui                  (Settings panel + "polish with AI" toggle / chat panel)
```

- **Pure prompt construction + grounding lives in `navigator-domain`** (`llm_prompt.rs`): given a
  `SubjectBrief` (and, in M2, a question + context bundle), produce the `system`/`user` message text.
  No I/O, unit-testable, so the *exact* instructions we send (especially the guardrails) are reviewable
  and stable — the same discipline as the deterministic `brief` templating.
- **The HTTP client + config lives in `navigator-app::llm`**: health check, model list, blocking and
  streaming completion, and the two task functions (`narrate_brief`, `answer_question`). It uses the
  app's existing `reqwest::Client`.
- **The UI** only renders and dispatches: a Settings section, an M1 "polish with AI" affordance on the
  brief, and an M2 chat panel. No HTTP, no prompt text in widgets.

### Why a module, not a `navigator-llm` crate (for now)

A dedicated crate would be cleaner long-term, but for a pilot the surface is small and only
`navigator-app` consumes it. Start as `navigator-app::llm`; promote to its own crate if/when analysis
or sync also need it. Keeping the pure prompt logic in `navigator-domain` already isolates the part
most worth testing.

---

## Settings

Add to `AppSettings` (`navigator-app/src/settings.rs`), all optional so existing files stay valid:

```rust
/// Enable local-LLM assisted narration / chat. Off until the user opts in.
#[serde(default, skip_serializing_if = "Option::is_none")]
pub llm_enabled: Option<bool>,
/// OpenAI-compatible base URL of the *local* model server, e.g. "http://localhost:1234/v1".
#[serde(default, skip_serializing_if = "Option::is_none")]
pub llm_base_url: Option<String>,
/// Model id to request (as reported by GET /models), e.g. "llama-3.1-8b-instruct".
#[serde(default, skip_serializing_if = "Option::is_none")]
pub llm_model: Option<String>,
```

Resolvers follow the established **env → settings → default** precedence (like `appview_url`):
`NAVIGATOR_LLM_ENABLED`, `NAVIGATOR_LLM_BASE_URL`, `NAVIGATOR_LLM_MODEL`. Defaults: disabled;
base URL `http://localhost:1234/v1` (LM Studio); model empty (forces the user to pick from
`/models`, or pass the server's single loaded model through). A **non-localhost** base URL is allowed
but surfaces a clear warning in the UI ("this will send your results to a non-local address") — the
guardrail against accidental exfiltration is explicit and visible, not hidden.

### Settings UI (`settings-ui` modal, new "AI assistant" section)

- A toggle: **Use a local AI model** (off by default).
- Base URL field, prefilled with the LM Studio default, with quick-pick chips for LM Studio / Ollama /
  llama.cpp ports.
- **Test connection** button → `GET /models`: shows ✓ with the server's loaded model(s) or a plain-
  language error ("No local model server found at …. Start LM Studio's server and try again.").
- Model dropdown populated from `/models`.
- A short privacy line: *"Prompts and your results are sent only to this local address and never leave
  your computer."* Turns to a warning if the URL isn't a loopback host.

---

## M1 — Narrate the Subject Brief (the pilot)

The brief is already a fully-resolved `SubjectBrief` (deterministic, computed off-thread, cached). M1
adds an **optional rephrasing pass**: feed the brief's *facts* to the local model and ask it to write
the casual-reader prose, replacing/augmenting the template strings the brief already produced.

### Grounding — the core safety mechanism

The model is a **rewriter, not a source of facts**. The prompt gives it the structured signals and
forbids inventing any others. Concretely, `navigator-domain::llm_prompt::narrate_system_prompt()`
fixes the rules; the user message is a compact JSON-ish fact sheet built from the `SubjectBrief`:

System prompt (sketch — the real text lives in code, unit-tested):

> You rewrite genetic-genealogy results into warm, plain-language prose for a non-expert reader.
> **Use only the facts provided.** Do not add haplogroup ages, place names, migration claims, dates,
> or numbers that are not in the facts. **Never** make health, medical, disease, trait, or clinical
> statements — this is genealogy and ancestry only. Preserve the stated confidence: if a placement is
> "tentative", say so. Keep it to N short paragraphs. Do not address the reader as a patient. If a fact
> is missing, omit that point — do not guess.

User message = the brief's facts, e.g. the headline, each `LineageBrief` (haplogroup, age_phrase,
origin_phrase, story, confidence_phrase, matched_ancestor), and the `TestBrief` (test_name,
what_it_tells, limitations, quality_phrase). Crucially these are **already-curated, already-rounded**
strings from the deterministic pipeline (e.g. "formed roughly 4,200 years ago", "tentative placement,
from a single test"), so even a faithful rewrite can only restate vetted content.

The reference-pack `story` text remains the factual backbone; the LLM's job is cohesion and tone, not
new genealogy. This keeps M1 firmly inside the existing [brief guardrails](subject-brief-simple-mode.md#tone--guardrails).

### Output handling & fallback

`narrate_brief` returns a `NarratedBrief { prose: Vec<Section>, model: String }` or — on **any** of
{feature disabled, server unreachable, timeout, HTTP error, empty/garbage response} — `None`. The UI
then shows the existing deterministic brief unchanged. Narration is **always** strictly additive: the
structured cards remain; the AI prose is shown as a clearly-labelled "DNA Story (AI-assisted)" block
above or beside them, with a small note naming the local model and a one-click **"Show the facts"** /
**regenerate** control. A user can never end up with *only* model output and no underlying data.

### Trigger & caching

- Lazy / explicit for the pilot: a **"Polish with AI"** button on the brief (Simple Mode overview),
  not automatic, so the first run is a deliberate user action and the latency is expected. Later it can
  default-on when enabled.
- `Command::NarrateBrief(guid)` → worker runs `app.narrate_brief(&brief)` → streamed
  `Event::BriefNarrationChunk { guid, text }` … `Event::BriefNarrationDone { guid, model }` (or
  `Event::BriefNarrationUnavailable { guid, reason }`). Streaming gives the immediate-mode UI live
  token output (the brief already has a worker/loading pattern to mirror).
- Cache the result on disk keyed by **(guid, pack version, analysis provenance, model id)** — same
  invalidation inputs as the brief cache plus the model — so re-opening a subject doesn't re-spend
  tokens, and changing model/inputs regenerates.

### Determinism note

LLM output is non-deterministic, so it is **not** unit-tested for exact wording (unlike the
deterministic templates). What *is* tested is the pure prompt builder (`llm_prompt`) — given a
`SubjectBrief`, assert the fact sheet contains the expected fields and the guardrail instructions, and
contains **no** data the brief didn't have. The client gets a thin integration test behind `#[ignore]`
(needs a running server), consistent with the repo's live-BAM/network test convention.

---

## M2 — Ask-my-results Q&A chat (same plumbing)

Once the client + grounding + settings exist, M2 adds a conversational panel without new transport.

- **Context bundle**: the subject's `SubjectBrief` + relevant reference-pack entries (the same facts
  M1 uses) become the grounding context. RAG-lite: we already have the curated, structured content;
  the "retrieval" is selecting the brief sections + pack lookups for the subject in view. No vector DB
  for the pilot.
- **Prompt**: `llm_prompt::answer_system_prompt()` reuses the M1 guardrails (facts-only, no health,
  preserve uncertainty) plus "if the answer isn't in the provided context, say you don't know rather
  than guessing." The user's question + context + short chat history go as messages.
- **Scope guard**: questions that fish for medical/clinical interpretation get a fixed deflection
  ("Navigator covers ancestry and lineage, not health"). This is enforced in the system prompt and
  reinforced by a lightweight keyword check on the way out.
- **UI**: a chat panel (Advanced first; Simple Mode later) scoped to the selected subject, with the
  model name and a persistent "answers are AI-generated from your results — verify against the data
  tabs" banner. Streamed via the same `Event` chunk pattern.

M2 inherits the local-only posture, the off-by-default toggle, and the streaming/cache plumbing from
M1. It is explicitly *not* an "agent" — no tool-calling into the analysis engine in this design; it
answers from the already-computed results. (Tool-use over the `App` query API is a tempting future
direction but out of pilot scope and a much bigger safety surface.)

---

## Plumbing summary

| Layer            | Addition                                                                          |
|------------------|----------------------------------------------------------------------------------|
| `navigator-domain` | `llm_prompt.rs`: pure prompt/fact-sheet builders + guardrail text (unit-tested) |
| `navigator-app`    | `llm.rs`: `LlmConfig`, `LlmClient` (`health`/`models`/`complete`/streaming), `narrate_brief`, `answer_question`; settings fields + resolvers; on-disk narration cache |
| `navigator-ui`     | Settings "AI assistant" section (toggle/URL/test/model); `Command::{TestLlmConnection, NarrateBrief, AskQuestion}` + matching streamed `Event`s; "Polish with AI" button (M1); chat panel (M2); i18n keys (en/es parity) |

All new user-facing strings route through the locale system (`self.tr(...)`), matching the i18n rule.

---

## Guardrails & tone (carried from the brief, made explicit for generation)

- **No health/clinical/trait output** — enforced in the system prompt *and* a post-generation keyword
  check; violations fall back to the deterministic brief.
- **Facts-only / no hallucinated genealogy** — the model only restates the vetted fact sheet; missing
  facts are omitted, not invented.
- **Honest uncertainty preserved** — confidence phrasing ("tentative", "from a single test") is part
  of the facts and must survive the rewrite.
- **Clear labelling** — AI prose is always badged as AI-assisted, names the local model, and sits
  *alongside* (never instead of) the structured facts, with a one-click reveal.
- **Local-only, visible** — loopback URL is the default and the happy path; a non-local URL triggers a
  warning; there is no hosted-provider path in this design.
- **Attribution intact** — reference-pack `sources` and any AppView-enriched markers stay on the
  structured cards regardless of narration.

---

## Phasing

1. **M0 — Client + settings + health.** `LlmConfig`/resolvers, `LlmClient::health`/`models`, the
   Settings "AI assistant" section with Test-connection. No generation yet; de-risks discovery and the
   local-only posture. Verifiable against a real LM Studio/Ollama server.
2. **M1 — Brief narration (the pilot).** `llm_prompt` builders + guardrails, `narrate_brief`, streamed
   `Command/Event`, the "Polish with AI" button, additive labelled rendering, cache + fallback.
3. **M2 — Q&A chat.** Context-bundle grounding, chat panel, scope guard, history. Same transport.
4. **M3 — Polish.** Optional default-on narration when enabled; Simple-Mode chat entry; expanded
   prompt tuning per model size; export of the AI story alongside the brief export.

## Open questions

- **Default base URL**: ship LM Studio's `:1234/v1` as the default, or detect a live server by probing
  the common ports on first enable?
- **Model size guidance**: small local models (3–8B) can drift from grounding instructions. Do we
  ship a recommended-model note, and/or tighten the fact sheet (fewer free-text fields) for robustness?
- **Streaming vs blocking** in the egui immediate-mode loop: streamed chunks are nicer but add Event
  volume; is a single blocking completion with a spinner enough for the pilot?
- **M2 history scope**: per-subject ephemeral chat (cleared on switch) vs persisted conversation in the
  workspace DB.
- **Telemetry**: none proposed (consistent with local-first), but do we want a local-only count of
  narration successes/fallbacks to tune prompts?
```
