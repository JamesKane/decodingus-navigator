# Subject Brief & Simple Mode — design

Status: **design** (2026-06-22). Implementation not started.

## Goal

Most DUNavigator users are **casual genetic genealogists**, not bioinformaticians. They tested one
person (themselves or a relative), uploaded the raw data, and want to know — in plain language —
*what their DNA says about their ancestry and lineages*. The current UI is built for power users:
multi-biosample projects, per-alignment analysis tabs, reconciliation audits, source provenance.
That is overkill and intimidating for the majority.

This subsystem adds:

1. **Simple Mode** — an app-level interface toggle that reduces Navigator to a single-person
   experience (hides Projects/Community and the advanced per-DNA-type tabs), defaulting on for new
   users. Power users flip to **Advanced** for today's full UI. The setting persists.
2. **Subject Brief** — a narrative, plain-language summary of one subject built from the raw
   analysis results, shown as the Overview in Simple Mode (and available as an exportable report in
   Advanced).

Decisions locked with the product owner (2026-06-22):

- **Structure**: app-level Simple/Advanced mode toggle (not just an Overview redesign).
- **Narrative content**: a **bundled offline reference pack** that makes briefs work with no
  network, **enriched from AppView/DecodingUs** (ages, provenance, descriptions) when online.
- **First-version sections**: Paternal & maternal lines · Ancestry composition · Genetic relatives
  · Your test & quality. (All four.)

Out of scope, explicitly: any **health / clinical / trait** interpretation. Briefs are genealogical
and ancestral only. Uncertainty is always surfaced, never hidden (see *Tone & guardrails*).

---

## Architecture overview

```
                       ┌─────────────────────────────────────────────┐
   reference pack ───► │ navigator-app::brief                         │
   (bundled JSON)      │   BriefBuilder::subject_brief(guid)          │ ──► SubjectBrief
   AppView enrich ───► │   composes existing queries + ref content    │     (pure, render-ready)
                       └─────────────────────────────────────────────┘
                                          │  Command/Event
                                          ▼
   navigator-ui  ──  Simple Mode shell  ──  renders SubjectBrief (read-only cards)
                 └─  Advanced Mode  ──  today's tabs + "Brief" export
```

The brief is **computed once off the UI thread** (worker), like the project Y-STR chart, and cached
per subject. It is a precomputed, fully-resolved model — the renderer does zero interpretation.

Crate layering is respected: `ui → app → {analysis, store, refgenome, sync} → {domain, du-*}`.
- **Reference-pack data model + plain-language templating**: `navigator-domain::brief` (pure, no I/O,
  unit-testable). This is where narrative strings are assembled from structured inputs.
- **Composition (pull the signals, load/enrich the pack, persist cache)**: `navigator-app::brief`.
- **Mode state + rendering**: `navigator-ui`.

---

## Part 1 — Simple Mode

### State & persistence

Add to `AppSettings` (`navigator-app/src/settings.rs`):

```rust
/// UI mode: "simple" (casual, single-person briefs) or "advanced" (full power-user UI).
#[serde(default, skip_serializing_if = "Option::is_none")]
pub ui_mode: Option<String>,
```

Resolver `ui_mode()` with precedence env (`NAVIGATOR_UI_MODE`) → settings → default. **Default for a
fresh workspace = `simple`**; once a user has created a Project or flipped to Advanced, persist their
choice. (First-run heuristic: if the workspace has 0 projects and ≤1 biosample, default Simple.)

UI state in `NavigatorApp`: `ui_mode: UiMode` (`enum UiMode { Simple, Advanced }`), loaded at
startup, toggled from the app bar (`[Simple ▾]` next to Settings), written through to settings via the
existing settings Command/Event path (see `settings-ui`).

### Nav gating (`chrome.rs`)

- **Simple**: nav bar shows only **My DNA** (the single-subject experience) and **Dashboard**
  (welcome / import). Projects and Community buttons hidden. The Subjects side-panel becomes a simple
  "who are you looking at" selector (most users have one subject — auto-select it and hide the list).
- **Advanced**: unchanged — Dashboard / Subjects / Projects / Community, all detail tabs.

The detail tabs (`DetailTab::{Overview, YDna, MtDna, Autosomal, Ancestry, Sources, IbdMatches}`) are
**hidden in Simple Mode**; the subject view *is* the brief (with a discreet "See the data →" link that
flips to Advanced on that subject for the curious). No analysis code changes — Simple Mode is a
presentation layer over the same `App`.

### Import flow

Simple Mode keeps the existing single-file Add Data path (already auto-detects type, creates the
run/alignment, resolves the reference, runs analysis). The casual user drops a file and waits; when
analysis completes the brief appears. No project required. The deep-analyze and batch/FTDNA project
import paths live only in Advanced.

---

## Part 2 — The Subject Brief

### Model (`navigator-domain::brief`)

A `SubjectBrief` is a render-ready tree of sections; each section is a small struct with already-
formatted strings plus the structured numbers the UI needs for any visual (donut, map pin, timeline).
Every narrative field is `Option` — a section degrades gracefully when its data or reference content
is missing (e.g. Y-only test → no ancestry section; haplogroup not in the pack → show the lineage
without the origin story).

```rust
pub struct SubjectBrief {
    pub headline: Headline,          // name, test summary chip, one-line "who you are"
    pub paternal: Option<LineageBrief>,   // Y; None for female sample or no Y data
    pub maternal: Option<LineageBrief>,   // mtDNA
    pub ancestry: Option<AncestryBrief>,
    pub relatives: Option<RelativesBrief>,
    pub test: TestBrief,             // always present
    pub caveats: Vec<String>,        // global uncertainty notes
    pub enriched: bool,              // true if AppView content was folded in
}

pub struct LineageBrief {
    pub haplogroup: String,          // terminal, e.g. "R-FGC29071"
    pub lineage_path: Vec<String>,   // root→tip (for an optional expandable trail)
    pub formed_ybp: Option<i32>,     // from du-domain Haplogroup (AppView) or pack
    pub age_phrase: Option<String>,  // "formed roughly 4,000 years ago"
    pub origin_phrase: Option<String>,// "associated with the spread of ... across NW Europe"
    pub story: Option<String>,       // 2–4 sentence curated narrative
    pub confidence_phrase: String,   // plain-language confidence ("strong / tentative placement")
    pub sources: Vec<String>,        // attribution for the narrative content
}

pub struct AncestryBrief {
    pub summary_phrase: String,      // "Predominantly Northwest European"
    pub super_pops: Vec<(String, f64)>,   // label, fraction (for the donut)
    pub fine_pops: Vec<(String, f64)>,
    pub interpretation: Option<String>,   // plain-language note on the mix
    pub method_note: String,         // "estimated from N genome-wide markers"
}

pub struct RelativesBrief {
    pub count: usize,
    pub items: Vec<RelativeItem>,    // who, why (plain-language signal), strength
    pub note: String,                // how to act on a match / privacy framing
}

pub struct TestBrief {
    pub test_name: String,           // "Whole-genome sequence", "Big Y-700", "23andMe chip"…
    pub what_it_tells: String,       // plain-language capability summary
    pub limitations: Option<String>, // "covers only the Y chromosome — no ancestry"
    pub quality_phrase: String,      // "high-quality (30× average depth)"
    pub quality_ok: bool,            // drives a ✓ / ⚠ chip
}
```

`navigator-domain::brief` owns the **templating**: pure functions like
`age_phrase(formed_ybp) -> String`, `quality_phrase(&CoverageResult, &TestType) -> (String, bool)`,
`ancestry_summary(&[SuperPopulationSummary]) -> String`. These are deterministic and unit-tested
(snapshot-style: given inputs → expected sentence), so the wording is reviewable and stable.

### Composition (`navigator-app::brief`)

`App::subject_brief(guid) -> Result<SubjectBrief, AppError>` pulls the already-available signals and
joins them to reference content:

| Section   | Signals (existing queries)                                              | Reference content needed |
|-----------|------------------------------------------------------------------------|--------------------------|
| Headline  | biosample, `default_alignment_for_subject`, test type                  | —                        |
| Paternal  | `haplogroup_consensus(Y)` → terminal + lineage + confidence            | Y haplogroup origin/age/story |
| Maternal  | `haplogroup_consensus(Mt)`; `cached_mt_profile` mutation count         | mt haplogroup origin/age/story |
| Ancestry  | `donor_ancestry`, `consensus_ancestry("FINE_ADMIXTURE")`               | population plain names + interp rules |
| Relatives | `ibd_suggestions`, `list_ibd_exchanges_for_subject`                    | signal→phrase map        |
| Test      | `SequenceRun.test_type`, `cached_coverage`, `cached_read_metrics`, `cached_sex` | test-type descriptions; quality thresholds |

All signal queries already exist (see the exploration map). No new analysis. The builder is mostly
glue + reference lookups + `navigator-domain::brief` templating.

### Reference pack (`brief-pack`)

A bundled, versioned, offline JSON pack — the missing narrative content. Shipped in the binary
(via `include_str!` or an on-disk `~/.decodingus/briefs/` cache seeded on first run) and **manifest-
verified** like the ancestry/IBD assets (`asset-manifest-verification`).

Schema (sketch):

```json
{
  "version": "2026.06",
  "y_haplogroups": {
    "R-M269":   { "formed_ybp": 6400, "origin": "the Pontic-Caspian steppe / early Europe",
                  "story": "R-M269 is the most common paternal lineage in Western Europe…",
                  "sources": ["YFull", "..."] },
    "...": {}
  },
  "mt_haplogroups": { "U5a": { ... } },
  "populations": { "EUR": { "name": "European", "blurb": "..." },
                   "GBR": { "name": "British", "blurb": "..." } },
  "test_types": { "WGS": { "what": "...", "limits": null },
                  "BIG_Y_700": { "what": "...", "limits": "Y chromosome only — no ancestry" } },
  "ancestry_rules": [ /* mix → interpretation phrases */ ]
}
```

Lookup is by haplogroup name with **ancestor fallback**: if the terminal (`R-FGC29071`) isn't in the
pack, walk up `lineage_path` to the nearest covered ancestor (`R-M269`) and narrate that ("On your
paternal line you descend from R-M269 … your specific branch is R-FGC29071"). This guarantees a
useful story for almost everyone even with a compact pack.

**AppView enrichment** (when online): `du-domain::Haplogroup` already carries `formed_ybp`, `tmrca_ybp`,
and a `provenance` blob. A `GET /api/v1/haplogroup/{name}` fetch (cached, TTL like the tree cache)
overrides/fills `formed_ybp`/age and can supply richer provenance. `SubjectBrief.enriched` records
whether this happened so the UI can show a subtle "live data" marker. Offline → pack values stand.

Pack authoring is content work, not code — start with the top ~100 Y and ~50 mt haplogroups (covers
the long tail via ancestor fallback) and the test-type catalog we already enumerate in
`navigator-domain::testtype`.

### Plumbing

- `Command::LoadSubjectBrief(guid)` → worker runs `app.subject_brief(guid)` → `Event::SubjectBrief
  { guid, brief }`. UI stores `Option<SubjectBrief>` + a loading flag; spinner while building.
- Rebuild triggers (mirror the project chart): on subject select, after analysis completes
  (`YHaplogroup`, `MtProfile`, ancestry estimate, `ReconciliationChanged`, coverage), and on
  Add-Data import for that subject. Cheap — it's cache reads + pack lookups.
- Optional persistence: cache the built brief as an on-disk artifact keyed by (guid, pack version,
  analysis provenance) so re-opening a subject is instant; invalidate when any input changes.

### Rendering (Simple Mode)

Card stack, one card per section, generous whitespace, minimal jargon, every technical term hover-
explained:

- **Headline**: name + test chip + a single friendly sentence.
- **Paternal / Maternal**: lineage name (big), age phrase, origin phrase, the 2–4 sentence story,
  a confidence chip, an expandable "lineage trail" (the root→tip path) for the curious.
- **Ancestry**: the existing donut (reuse the Ancestry tab's plot) + plain-language summary + fine
  breakdown rows.
- **Relatives**: match count + cards (who / why / strength) + a privacy-minded "what to do" note.
- **Your test**: test name, what it tells you, limitations, a ✓/⚠ quality chip.

Each card has an **Export** affordance (reuses `navigator-app::export` patterns → HTML/PDF) so a
casual user can share/print a "DNA Story" report. In Advanced Mode the whole brief is reachable as a
**"Brief" export** from the Overview without changing the existing tabs.

---

## Tone & guardrails

- **No health/clinical/trait claims.** Genealogy and ancestry only.
- **Honest uncertainty.** Confidence is always phrased ("strong placement" vs "tentative — based on a
  shallow test"). Low-confidence ancestry or haplogroup placements are labelled, not dropped. Coverage
  below useful thresholds yields a ⚠ quality chip and a caveat.
- **Attribution.** Narrative content carries `sources`; AppView-enriched values are marked.
- **Geography is broad, not deterministic.** Origin phrases describe lineage *associations and spread*,
  not "you are from X". Avoid nationalist/ethnic over-claiming.
- **i18n.** All templated phrases route through the locale system (en/es parity-tested). Reference-pack
  prose is English-first; translatable packs are a later concern.

---

## Phasing

1. **M1 — Simple Mode shell**: `ui_mode` setting + toggle + nav/tab gating + single-subject selection.
   No briefs yet; Overview shows today's `overview_dashboard`. De-risks the mode plumbing alone.
2. **M2 — Brief model + Test & Paternal/Maternal sections** off a minimal bundled pack (ages + stories
   for top haplogroups; test-type descriptions). Worker plumbing + caching + Simple-Mode cards.
3. **M3 — Ancestry + Relatives sections**; reuse the donut; signal→phrase mapping for matches.
4. **M4 — AppView enrichment** (haplogroup age/provenance fetch + cache) and the **Export** report.
5. **M5 — Pack breadth + polish**: expand reference coverage, hover-glossary, first-run defaulting,
   "See the data →" bridge to Advanced.

---

## Open questions

- First-run mode default heuristic — purely "0 projects & ≤1 subject", or an explicit welcome choice?
- Reference-pack distribution: in-binary `include_str!` vs a CDN-downloaded, manifest-verified asset
  (lets content update without an app release, at the cost of a first-run fetch).
- Relatives in Simple Mode: show only when the user is signed in to the network, or always with a
  "connect to find matches" prompt?
- Does Simple Mode need a maternal/paternal **migration map** (PCA/`population_lonlat` already exists)
  in M3, or is that M5 polish?
