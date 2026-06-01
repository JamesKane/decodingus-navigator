# AppView Lifecycle Management

> **Scope reduction (2026-06).** The AppView is **no longer a network mirror.** It does
> not ingest `biosample` / `sequencerun` / `alignment` / `genotype` / etc. via a
> full-CRUD firehose subscription. Its network-facing role collapses to **two narrow
> flows**: (1) maintain the **known-variant catalog** via direct proposal submission, and
> (2) produce **coverage aggregates on demand** from Navigator-published summary records.
> Everything else moves to the researcher's Navigator workspace. The previous
> per-collection ingestion design is preserved in git history; what it's replaced with is
> below.

---

## 1. AppView responsibilities (reduced)

| Responsibility | Mechanism | What the AppView stores |
|:---|:---|:---|
| **Known-variant catalog** | Direct proposal submission from Navigator (curation API) | The variant/tree catalog (curator-owned: `core.variant`, `tree.haplogroup`) |
| **Coverage aggregation** | On-demand aggregation over Navigator-published public summary records | Nothing persistent per-sample — aggregates are computed at query time |

Explicitly **not** an AppView job anymore: mirroring biosamples, sequence runs, alignments,
genotypes, STR profiles, match records; computing/storing per-sample coverage centrally;
firehose-driven CRUD with orphan handling.

---

## 2. Known-variant catalog — direct submission

Variants reach the catalog by **Navigator submitting proposals directly** to an AppView
curation API — decoupled from PDS records, no firehose harvest.

1. **Submit.** Navigator posts a variant/branch proposal (candidate position, ref/alt,
   supporting sample evidence, proposed placement) to the AppView submission endpoint,
   authenticated as the Navigator user.
2. **Pool & consensus.** The AppView aggregates proposals for the same candidate across
   submitters; a confidence threshold gates promotion.
3. **Curator review.** Above threshold → flagged for curator review and naming.
4. **Promote.** Accepted proposals update the named catalog (`naming_status: Named`) with a
   new catalog version.

Versioning/conflict resolution (§5) applies **here** — this is the one place the AppView
owns authoritative, mutable network state.

> Privacy: proposals carry only the variant evidence needed for placement, not raw reads.
> Consistent with the core principle — only computed summaries leave the device.

---

## 3. Coverage aggregation — on demand over published summaries

The AppView does not ingest or store per-sample coverage. Instead:

1. **Navigator publishes** per-sample coverage **summaries** as public PDS records
   (biosample/alignment-level metrics — mean/median coverage, callable fractions per
   contig — not raw data).
2. **AppView aggregates on demand.** When a population/network-level coverage view is
   requested, the AppView reads the relevant published summary records and aggregates at
   query time. It keeps **no normalized mirror** and runs no CRUD ingestion handlers.

> **Open mechanism:** on-demand aggregation needs a way to *discover* which public summary
> records exist for a cohort. Options: a lightweight firehose-derived **URI index** (record
> pointers only, not a full mirror), or a query against a relay/AppView-of-record. To
> decide with the decodingus team — this is the only residual use of the firehose, and it
> is discovery-only, not state synchronization.

---

## 4. What was removed (the firehose mirror)

The following are **removed from the decodingus codebase** under this re-scope:

- The `com.atproto.sync.subscribeRepos` full-CRUD subscription across all collections.
- Per-collection event handlers and normalized tables for `biosample`, `sequencerun`,
  `alignment`, `genotype`, `imputation`, `populationBreakdown`, `instrumentObservation`,
  `matchConsent`, `matchList`, `matchRequest`, `strProfile`.
- Orphan handling and `at_*` sync-tracking fields on those mirrored tables.
- The branch-discovery harvest that read `privateVariants` from ingested biosample records
  — superseded by the §2 submission flow.

**Backfeed (AppView → user PDS), unchanged in principle:** records "written BY DecodingUs"
(`matchList`, `haplogroupAncestralStr`) still require an explicit user-granted scope; prefer
AppView-owned records the client reads. See
[11-Auth-and-Permissions.md](./11-Auth-and-Permissions.md) §5.

---

## 5. Versioning & conflict resolution (variant catalog only)

These concepts now apply solely to the variant/tree catalog (§2), not to a record mirror:

| Field | Purpose |
|:---|:---|
| catalog `version` | Optimistic locking on catalog updates |
| proposal `cid` | De-duplicate resubmitted proposals |
| `naming_status` | `Unnamed` → `PendingReview` → `Named` lifecycle |

On a conflicting catalog update: higher version wins; on tie, later timestamp; log for
curator review.
