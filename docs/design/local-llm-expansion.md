# Local-LLM expansion — design

Status: **Proposed** (2026-06-22). Follows on from the shipped pilot in
[local-llm-integration.md](local-llm-integration.md) (M0–M3, on `main`). No code yet — this is the
plan for where the integration goes next.

## Where the pilot left us

The pilot gave us a clean, reusable spine:

- One streaming transport — `navigator-app::llm::chat_complete_streaming` (SSE over the existing
  `reqwest` client; reasoning suppressed server-side via `chat_template_kwargs`, stripped as fallback).
- Pure, unit-tested grounding in `navigator-domain::llm_prompt` — `narrate_system_prompt`,
  `answer_system_prompt`, `narrate_fact_sheet(&SubjectBrief)`, `mentions_health`, `health_deflection`.
- Two surfaces: **M1 "Polish with AI"** (narrate the brief) and **M2 "Ask my results"** (chat).

The constraint that makes it safe is also what limits it: **both surfaces are grounded in exactly one
fact sheet — `narrate_fact_sheet(&SubjectBrief)`** (`llm.rs:362`). The `SubjectBrief` carries only
four sections (paternal line, maternal line, ancestry, test & quality), so that is all the AI can
narrate *or* answer about. Everything else the app computes — Y-STR panels, IBD/network matches,
private-Y variants, mtDNA mutations, the genetic sex call — is invisible to the model. Ask the chat
"how many novel Y variants do I have?" or "who am I related to in here?" and it correctly says it
doesn't know, because those facts were never in its context.

## Goal

Widen what the AI can speak to, **without** widening the safety surface. Every expansion stays inside
the pilot's guardrails: facts-only (the model restates already-computed, already-curated values, never
sources new ones), no health/clinical language, honest uncertainty preserved, local-only, additive,
off-by-default. No tool-calling / agent access to the analysis engine — still explicitly out of scope.

The unlock is a single idea: **a richer grounded fact sheet, assembled from curated per-signal
summaries.** Build that once and three things fall out of it — a smarter chat (M4), cheap per-tab
"Explain this" buttons (M5), and a project-level summary (M6).

---

## M4 — Broaden the chat grounding (recommended first)

**One change, amplifies everything.** Today the chat's system message is `answer_system_prompt()` +
`narrate_fact_sheet(&brief)` (`llm.rs:364`). M4 replaces the fact sheet with a broader
**results context** that appends curated summaries of the other signals, so the existing "Ask my
results" panel can answer about all of them with no new UI plumbing.

### Keep narration lean, make chat broad

Narration (M1) deliberately produces a tight 2–4 paragraph *story* — it should **not** grow to cover
STR marker counts and IBD match lists. So we split the two consumers rather than fattening
`SubjectBrief`:

- **M1 narration** keeps using `narrate_fact_sheet(&SubjectBrief)` — unchanged.
- **M4 chat** uses a new `results_fact_sheet(&ResultsContext)` that is *`narrate_fact_sheet` plus the
  extra signal sections*.

This keeps the pilot's narration grounding byte-for-byte stable (no regression risk to the shipped
feature) while the chat gets the wider context.

### New pure type + builder (navigator-domain)

A new `ResultsContext` aggregate in `navigator-domain` (sibling to `SubjectBrief`), holding the brief
plus optional, **already-curated** per-signal summaries — strings, not raw analysis structs, so the
domain layer stays pure and the exact text we send is unit-testable:

```rust
pub struct ResultsContext {
    pub brief: SubjectBrief,
    pub sex: Option<String>,            // "Genetic sex: XY (high confidence)"
    pub ystr: Option<YStrSummary>,      // panels + marker counts + conflict note
    pub private_y: Option<String>,      // "12 novel variants (new-branch candidates), 3 off-tree"
    pub mt_mutations: Option<String>,   // "41 differences from rCRS across HVR1/HVR2/coding"
    pub ibd: Option<IbdSummary>,        // match count + closest-relationship band
}
```

`results_fact_sheet(&ResultsContext) -> String` starts from `narrate_fact_sheet(&ctx.brief)` and
appends a labelled block per present signal (absent signals are simply omitted, exactly as the brief
sections are today). Each section is **summary-level**, not a data dump — see the per-signal notes
below for why.

### Assembling the context (navigator-app)

A new `App::results_context(guid)` builds the `ResultsContext` by calling the existing accessors and
reducing each to its curated summary string:

| Signal | Accessor (verified) | Reduce to |
|--------|---------------------|-----------|
| Sex | `analysis.rs:169 cached_sex(...)` | call + confidence phrase |
| Y-STR | `import_profiles.rs:84 list_str_profiles(guid)` | per-panel name + marker count + any conflict (reuse `strpanel` classification); **not** the 37–700 raw marker values |
| Private-Y | `queries.rs:123 donor_private_y(guid) -> PrivateBucket` | counts: novel (new-branch candidates) vs off-tree, with the paralog caveat |
| mtDNA mutations | `import_profiles.rs:489 mtdna_variants(mtdna_id)` | count + region spread; notable named mutations only |
| IBD | `sync.rs:338 ibd_suggestions()` + `ibd_exchange.rs:548 list_ibd_exchanges_for_subject(guid)` | match count + closest relationship band; PII handled carefully (see below) |

The chat path (`answer_question_streaming`, `llm.rs:344`) changes one line: build `results_context`
instead of `subject_brief`, and call `results_fact_sheet` instead of `narrate_fact_sheet`.

### Per-signal grounding notes

- **Y-STR** — summarize, don't dump. A Big-Y profile has hundreds of markers; sending them all blows
  the token budget and adds no answerable value. Send "FTDNA Y-111 panel (111 markers); YSEQ Alpha;
  2 markers in conflict between panels." Frame STR results as lineage-pattern, never trait/health.
- **Private-Y** — counts only, with the existing distinction novel (candidate new branch) vs off-tree
  (unknown branch depth) from `PrivateBucket`. Carry the paralog/structural caveat so the model can't
  overstate ("candidates for", not "you have a new haplogroup").
- **mtDNA mutations** — count + region spread relative to rCRS; surface only notable named mutations.
  Heteroplasmy framed as a lineage-depth marker, never clinical. `mentions_health` already guards the
  output.
- **IBD** — the PII-sensitive one. These are *the user's own workspace relatives*, so it's their data,
  but the chat should answer about **relationship structure** ("3 close matches; the closest shares
  ~210 cM, consistent with a 2nd–3rd cousin") rather than reciting names unprompted. A short scope
  line in the section keeps it framed as DNA-inferred, not genealogically verified.
- **Sex** — trivial: the call + confidence. Also a candidate to fold into the brief's test section
  later, but cheapest to add here.

### Cost / caching

`subject_brief` is already cached; the new accessors add per-question DB work. Build
`results_context` **once when the chat panel opens** (or cache the assembled fact sheet keyed by the
same provenance inputs as the brief) and reuse it across the conversation, rather than re-assembling on
every turn. The fact sheet is plain text, so this is a cheap memo.

### UI

Almost nothing. The "Ask about your results" panel already exists. Two small touches:

- Update the suggested-question chips / placeholder to hint the new coverage ("How many novel Y
  variants do I have?", "Who am I most related to here?").
- The existing "answers are AI-generated — verify against the data tabs" banner already covers the
  broader scope.

### Testing

Same discipline as the pilot: unit-test the pure builder. Given a `ResultsContext`, assert each
present section appears, absent signals don't, confidence/caveats survive, and `mentions_health` is
clean on the assembled sheet. The accessors get the repo's usual `#[ignore]` live-data treatment.

---

## M5 — Per-tab "Explain this" narration (follow-on, cheap after M4)

Once M4 exists, each signal already has a curated summary builder. M5 reuses those to add an
**"Explain this"** affordance on individual tabs (Y-DNA, IBD, mtDNA), each narrating *just* that
signal in plain language — the highest-ROI Tier-1 targets from the recon:

1. **Y-STR panel** (Y-DNA tab) — "what these markers mean for your line."
2. **IBD / network matches** (IBD tab) — "what these shared segments imply."
3. **Private-Y variants** (Y-DNA → Private-Y) — novel vs off-tree explained.

The pattern per target is mechanical and mirrors M1: a focused system prompt + the signal's fact
section (already built in M4) → `chat_complete_streaming`, with a new `Command::Narrate{Signal}` /
`Event::{Signal}Narration` pair and a button in the tab. Cache keyed by
`sha256(model_id + system_prompt + section_text)`.

M5 is deferred behind M4 precisely because **M4 builds the per-signal summaries M5 needs** — doing the
grounding first makes the buttons nearly free.

---

## M6 — Project-level AI summary (separate audience)

Narrate a whole **project / batch** (the coverage-report / batch-import workflow) rather than a single
subject — e.g. "this project has 24 samples; coverage ranges 4–38×; haplogroups cluster in R-M269."
Grounded in the existing project report + per-subject briefs. Distinct from M4/M5 (group-level, not
per-subject) and lower priority, but a natural home for the batch users.

---

## Sequencing & rationale

1. **M4 — broaden chat grounding.** One surface, reuses all existing plumbing, makes every signal
   answerable, and produces the per-signal summary builders the rest depends on. Highest leverage.
2. **M5 — per-tab narrate buttons.** Cheap once M4's summaries exist; most visible per-tab payoff.
3. **M6 — project summary.** New audience; do last.

## Scope boundaries (unchanged from the pilot)

Facts-only restatement; no health/clinical/trait output (`mentions_health` on every new section and
output); local-only, no cloud fallback; additive (AI never replaces the structured tabs); grounded
context only (no open-ended retrieval); off by default; **no tool-calling / agent access to the
analysis engine.**

## Open questions

- **Token budget.** All signals present on a Big-Y + IBD-rich subject could make a large system
  message. Do we cap sections, or tier the context by what the question seems to be about?
- **IBD PII in prompts.** Names vs relationship-bands-only — confirm the default. Leaning bands-only
  unless the user names a match.
- **Context freshness.** Cache the assembled `results_context` per chat session, or rebuild when any
  underlying analysis is re-run mid-session?
- **Sex placement.** Add via M4's context, or fold into the brief's test section so narration covers
  it too?
