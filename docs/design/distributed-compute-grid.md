# DecodingUs Grid — distributed community realignment & analysis

Status: **design / specification only** (no code). Cross-repo: **Navigator** (edge worker) +
**AppView** (`decodingus`, coordinator) + **shared** (`decodingus-shared`, wire records).

A Seti@Home / Folding@Home–style layer. The AppView publishes a list of **work units** — public
ENA samples. Navigator instances volunteer to reserve a unit for a bounded lease, fetch the data
from ENA, (re)align it to **CHM13v2 / hs1**, run the full analysis stack, submit signed results,
and release the lease. Validated contributions earn **compute credit** on a public **leaderboard**
and a capped, positive bump to the contributor's **reputation**.

The payoff: a growing, uniformly-hs1-aligned, community-computed corpus of Y/mt haplogroups,
ancestry, coverage, and callability over public ENA data — derived once, verified, and shared —
without any central compute cost.

---

## 1. Decisions locked (with rationale)

These four forks were decided before drafting; the doc is built on them.

| # | Decision | Choice | Why |
|---|----------|--------|-----|
| D1 | **Aligner integration** | **minimap2 via `minimap2-rs` FFI** (`static` + `simde`) | Single-binary spirit preserved; independently the chosen backend in [`realignment-module.md`](realignment-module.md). The Grid *consumes* that module's engine — it does not re-decide it. |
| D2 | **Trust model** | **Adaptive replication** | Untrusted nodes run in shadow/quorum; reputation graduates them to trusted single-run + random spot-recheck. BOINC-proven; K× cost only where trust is unearned. |
| D3 | **First cut** | **Staged — CRAM-passthrough first** | Phase 1 claims ENA samples that already have a CRAM/BAM, skips realignment, and just runs the stack. Proves the whole lease→submit→validate→leaderboard loop with zero aligner risk. Phase 2 adds FASTQ→minimap2. |
| D4 | **Result home** | **Contributor PDS + AppView canonical** | Contributor publishes fed records into their *own* repo, tagged with the ENA accession as subject + a `computedBy`/provenance block; AppView ingests, dedups by `(accession, method)`, promotes a canonical copy. Keeps federation; requires the new subject≠contributor split. |

---

## 2. What already exists (reuse) vs. what's greenfield

Reconnaissance across the three repos. **Reuse aggressively; the coordination substrate is
mostly already in the AppView DB.**

### Reuse — AppView (`decodingus`)
- **Dormant lease/node/submission scaffold** — `migrations/0008_fed.sql`, built for almost exactly
  this and never wired:
  - `fed.pds_node` — node registry: `capabilities JSONB`, `status`, `last_heartbeat`, `software_version`.
  - `fed.pds_registration` — a **lease**: `leased_by_instance_id`, `lease_expires_at`, `processing_status`, index on expiry.
  - `fed.pds_heartbeat_log` — `load_metrics`, `processing_queue_size`.
  - `fed.pds_submission` — a work queue with a status lifecycle (`PENDING/ACCEPTED/REJECTED/SUPERSEDED`).
- **Edge auth, solved** — `du-web/src/sig.rs::verify_signed(pool, did, message, sig)` + `fed.device_key`
  (Ed25519 `did:key`), with a ±300 s replay guard (`ensure_fresh_ts`). Every `/exchange/*` handler
  already uses this pattern; the Grid reuses it verbatim.
- **ENA study catalog** — `pubs.genomic_study` (`accession`, `source ENA|NCBI_*`) + `du-external`'s
  `EnaClient::study()` (ENA Portal API, no creds). Study-*metadata* only today.
- **Reputation subsystem** — `social.reputation_event` append-only ledger + cached
  `social.user_reputation_score`, `record_event`/`record_once`, seeded event types. Per-**user**.
- **Jetstream ingest** — `du-jobs/src/jetstream.rs` already dispatches `com.decodingus.*` records by
  NSID into `fed.*` upserts. Grid result records ride this same pipe.
- **`du-jobs` scheduler** — in-process interval jobs (`scheduler.rs`); the natural home for the
  lease-reaper and the validation/canonicalization job.

### Reuse — Navigator (`DUNavigator`)
- **Full analysis stack, per-alignment** — `navigator-app/src/analysis.rs`: `run_unified_metrics`
  (coverage+read_metrics+sex), `run_sv`, `run_denovo_caller`, + `haplogroup.rs`
  (`assign_y_haplogroup`, `assign_mtdna_haplogroup_from_alignment`, `place_{y,mt}_consensus`,
  `estimate_ancestry_from_consensus`). Each computes + persists a versioned artifact.
- **The realignment engine** — [`realignment-module.md`](realignment-module.md) fully specifies
  revert → minimap2-rs align → sort/markdup/CRAM → register, incl. the `.mmi` index cache and
  per-tech presets. **The Grid is that module's first heavy consumer.**
- **Reference fetch/cache** — `navigator-refgenome::Gateway::resolve_reference("chm13v2", …)`
  (streaming download, SHA-pinned, on-disk cache).
- **Durable publish outbox** — `publish_*` (`ibd_exchange.rs`) → `enqueue_publish` → `sync_outbox`
  → `drain_outbox`, idempotent via `sync_state`. Grid results publish through this.
- **Device-key signing + signed AppView client** — `navigator-sync::DeviceKey` (Ed25519, keychain);
  the `exchange_get_poll` pattern (`ibd_exchange.rs:594`) signs `did/ts/sig` params. Copy it.
- **CLI** — clap subcommands in `navigator-ui/src/cli.rs`; add `contribute`.

### Reuse — shared (`decodingus-shared`)
- `du-domain::fed` wire records — `BiosampleRecord`, `SequenceRunRecord`, `AlignmentRecord`
  (coverage), `PopulationBreakdownRecord`; `RecordMeta` + `$type` + `WireF64` envelope;
  `ExternalId { namespace, value }` (already lists **ENA** as a namespace); `at://`-ref linking.
- `du-atproto::signature::verify_did_key` + `did.rs` (`did:key` ↔ Ed25519).

### Greenfield (net-new)
1. **ENA sequence fetch** — Navigator has *zero* download code; AppView pulls study metadata only.
   Need run-level file resolution (FASTQ/CRAM URLs, md5, bytes) + a resilient downloader.
2. **Work-unit coordination** — the lease-*acquisition* SQL (`… WHERE lease_expires_at < now() …
   FOR UPDATE SKIP LOCKED RETURNING`), the `grid.work_unit` catalog, and the claim/submit/validate
   endpoints. The tables partly exist; the logic does not.
3. **Subject ≠ contributor** — today a fed record is authored-by *and about* the same repo DID.
   "DID A computed a result about ownerless ENA sample X" needs a `subject` + `computedBy` +
   structured `provenance` (software/version/reference/aligner) on the records.
4. **Compute credit + leaderboard** — a `WORK_UNIT_COMPLETED` reputation event + a dedicated
   compute-credit tally and a public leaderboard endpoint/view.
5. **Adaptive-replication validator** — digest comparison, quorum, trust tiers, spot-recheck,
   divergence penalties.

---

## 3. Architecture & lifecycle

```
                         ┌──────────────────────── AppView (decodingus) ────────────────────────┐
   ENA Portal API  ──►   │ curate: pubs.genomic_study + filereport → grid.work_unit (AVAILABLE)  │
                         │ coordinate: claim (lease) · heartbeat · submit · release              │
                         │ validate: adaptive replication over result DIGESTS → CANONICAL        │
                         │ credit: reputation event + compute-credit tally → leaderboard         │
                         └───────▲───────────────────────────┬──────────────────────────────────┘
              signed (device key)│                           │ signed claim / lease
                                 │                           ▼
   ┌───────────────────────── Navigator node (edge worker) ─────────────────────────┐
   │ register(capabilities) → claim(N) → for each unit:                             │
   │   fetch from ENA (FASTQ | CRAM)  ──►  P1: CRAM passthrough (skip align)         │
   │                                       P2: FASTQ → minimap2-rs → CHM13 (realign) │
   │   → run full analysis stack (coverage/sex/SV/Y/mt/ancestry)                     │
   │   → build signed result DIGEST + fed records (subject=ENA acc, computedBy=self) │
   │   → publish records (sync_outbox) + POST /grid/submit → release lease           │
   │   → clean scratch                                                              │
   └────────────────────────────────────────────────────────────────────────────────┘
```

**Work-unit state machine (AppView):**

```
AVAILABLE ──claim──► LEASED ──submit──► SUBMITTED ──validate──►┬─(quorum agree)─► CANONICAL
    ▲                   │                                       └─(diverge)──────► CONTESTED ─► (re-quorum)
    │            lease expiry / release                                                   │
    └───────────────────────────────────────────────────────────────────────────────────┘
CANONICAL ──(later contradicting result)──► CONTESTED    ·    any state ──curator──► RETIRED
```

A unit needs `required_replicas` (default 2, computed from the trust of submitters — §6). It
reaches **CANONICAL** when a quorum of *agreeing* digests exists; a trusted node can satisfy the
quorum alone, with ~5 % of such units randomly re-queued for a shadow replica.

---

## 4. AppView — data model & coordination

New Postgres schema **`grid`** (migration `0059_grid.sql`, next in sequence). Reuse `fed.pds_node`
+ `fed.device_key`; everything work-specific is new so we don't overload `fed.pds_submission`'s
existing semantics.

### 4.1 Tables (sketch)

```sql
CREATE SCHEMA IF NOT EXISTS grid;

CREATE TYPE grid.unit_state AS ENUM
  ('AVAILABLE','LEASED','SUBMITTED','CANONICAL','CONTESTED','RETIRED');
CREATE TYPE grid.data_kind  AS ENUM ('CRAM','BAM','FASTQ');

-- One ENA sample = one work unit (covers its runs; analysis is per-biosample).
CREATE TABLE grid.work_unit (
  id                BIGINT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
  ena_sample_acc    TEXT NOT NULL UNIQUE,      -- SAMEA… / ERS…
  ena_study_acc     TEXT,                       -- PRJEB… / ERP…
  data_kind         grid.data_kind NOT NULL,    -- P1 curates CRAM/BAM; P2 opens FASTQ
  run_manifest      JSONB NOT NULL,             -- [{run_acc, urls[], md5[], bytes, layout, platform, read_type}]
  est_bases         BIGINT,                     -- ENA base_count → credit weight & size preflight
  est_download_bytes BIGINT,
  reference_build   TEXT NOT NULL DEFAULT 'chm13v2.0',
  stack_floor       TEXT,                       -- min analysis stack semver a result must meet
  required_replicas SMALLINT NOT NULL DEFAULT 2,
  state             grid.unit_state NOT NULL DEFAULT 'AVAILABLE',
  canonical_digest  JSONB,                       -- set on CANONICAL
  priority          INT NOT NULL DEFAULT 0,
  created_at        TIMESTAMPTZ NOT NULL DEFAULT now()
);
CREATE INDEX ON grid.work_unit (state, priority DESC) WHERE state = 'AVAILABLE';

-- Lease. Purpose-built (cleaner than repurposing fed.pds_registration's PDS-cursor semantics).
CREATE TABLE grid.lease (
  id            BIGINT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
  work_unit_id  BIGINT NOT NULL REFERENCES grid.work_unit(id),
  node_did      TEXT NOT NULL,                  -- resolves to a user via fed.device_key
  instance_id   TEXT NOT NULL,                  -- device/install id from the node
  claimed_at    TIMESTAMPTZ NOT NULL DEFAULT now(),
  expires_at    TIMESTAMPTZ NOT NULL,           -- claimed_at + requested lease (bounded, §4.3)
  last_heartbeat TIMESTAMPTZ NOT NULL DEFAULT now(),
  state         TEXT NOT NULL DEFAULT 'ACTIVE'  -- ACTIVE | COMPLETED | RELEASED | EXPIRED
);
CREATE INDEX ON grid.lease (state, expires_at);
-- At most one ACTIVE lease per (unit, node); many nodes may hold replicas of one unit.
CREATE UNIQUE INDEX ON grid.lease (work_unit_id, node_did) WHERE state = 'ACTIVE';

-- A submitted result = the signed digest + pointers to the contributor's PDS fed records.
CREATE TABLE grid.submission (
  id             BIGINT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
  work_unit_id   BIGINT NOT NULL REFERENCES grid.work_unit(id),
  node_did       TEXT NOT NULL,
  contributor_user_id UUID REFERENCES ident.users(id),  -- resolved from node_did
  digest         JSONB NOT NULL,                -- canonical discrete calls (§5.2)
  digest_sig     TEXT NOT NULL,                 -- device-key signature over the canonical digest bytes
  stack_version  TEXT NOT NULL,
  aligner        TEXT,                          -- "minimap2-rs <ver>/<preset>" or NULL (passthrough)
  record_refs    JSONB NOT NULL,               -- {biosample: at://…, coverage: at://…, ancestry: …, …}
  verdict        TEXT NOT NULL DEFAULT 'PENDING', -- PENDING | AGREED | DIVERGENT | SUPERSEDED
  submitted_at   TIMESTAMPTZ NOT NULL DEFAULT now()
);
CREATE INDEX ON grid.submission (work_unit_id, verdict);

-- Compute credit ledger (leaderboard), distinct from social reputation so it can't be
-- farmed to dominate social gates; a capped slice ALSO feeds social.reputation_event.
CREATE TABLE grid.credit (
  id            BIGINT GENERATED ALWAYS AS IDENTITY PRIMARY KEY,
  user_id       UUID NOT NULL REFERENCES ident.users(id),
  work_unit_id  BIGINT NOT NULL REFERENCES grid.work_unit(id),
  cobblestones  BIGINT NOT NULL,               -- credit magnitude (§6.3)
  reason        TEXT NOT NULL,                 -- CANONICAL_FIRST | QUORUM_AGREE | SPOTCHECK_PASS
  awarded_at    TIMESTAMPTZ NOT NULL DEFAULT now(),
  UNIQUE (user_id, work_unit_id)               -- one credit per user per unit
);
CREATE INDEX ON grid.credit (user_id);
```

Node liveness reuses **`fed.pds_node`** (register capabilities/status/heartbeat there) — no new
node table. The migrations test (`du-db/tests/migrations.rs`) gets the new tables added to its list.

### 4.2 The claim (the one piece of genuinely new concurrency)

Atomic multi-worker claim — pick available units the node can handle, lease them, all in one
statement:

```sql
WITH picked AS (
  SELECT id FROM grid.work_unit
  WHERE state = 'AVAILABLE'
    AND data_kind = ANY($caps_kinds)          -- node capability filter (CRAM only? FASTQ too?)
    AND est_download_bytes <= $max_bytes
  ORDER BY priority DESC, id
  LIMIT $n
  FOR UPDATE SKIP LOCKED                        -- the BOINC-scale multi-worker claim primitive
)
UPDATE grid.work_unit w SET state = 'LEASED'
FROM picked WHERE w.id = picked.id
RETURNING w.id, w.ena_sample_acc, w.run_manifest, w.reference_build;
-- then INSERT grid.lease rows (unit_id, node_did, instance_id, expires_at = now()+lease)
```

For **replication**, a unit may be re-offered to *additional* nodes while `SUBMITTED` but below
quorum: the reaper/validator flips such units back to `AVAILABLE` with a bumped `required_replicas`
so a second independent node picks them up (never the same `node_did` — enforced by the lease
unique index + a "not already a submitter" filter).

### 4.3 Lease honesty & reclamation
- **Bounded TTL.** Node requests a lease of *X* days; AppView clamps to `[min, max]` (e.g. 1–14 d)
  sized against `est_download_bytes`/`est_bases` so a node can't hoard the pool.
- **Heartbeat renewal.** `POST /grid/heartbeat` (signed) updates `last_heartbeat` and may extend
  `expires_at` while progress continues (carries `stage` + `pct` for UI/telemetry).
- **Reaper** (`du-jobs` interval): `UPDATE grid.lease SET state='EXPIRED' WHERE state='ACTIVE' AND
  expires_at < now()`; the freed unit returns to `AVAILABLE`. Straggler mitigation: a unit one
  replica short of quorum with a stale lease is re-offered early.
- **Voluntary release** on shutdown/cancel so units recycle fast.

### 4.4 Signed edge endpoints (`/api/v1/grid/*`)

All mutations verify via `sig::verify_signed(did, canonical_message, sig)` + `ensure_fresh_ts`,
exactly like `/exchange/*`. Canonical messages get byte-for-byte twins in a shared
`messages::grid` module (mirroring `exchange::messages`) so Navigator and AppView agree.

| Method | Path | Auth | Purpose |
|--------|------|------|---------|
| POST | `/grid/node/register` | signed | Upsert `fed.pds_node` capabilities (data kinds, threads, disk budget, OS/arch, stack version). |
| POST | `/grid/node/heartbeat` | signed | Node liveness + load (not per-lease). |
| POST | `/grid/claim` | signed | Lease up to N units matching capabilities → manifests + `expires_at`. |
| POST | `/grid/heartbeat` | signed | Per-lease progress + optional TTL extension. |
| POST | `/grid/submit` | signed | Digest + `digest_sig` + `record_refs`; marks lease `COMPLETED`, unit `SUBMITTED`. |
| POST | `/grid/release` | signed | Abandon a lease; unit → `AVAILABLE`. |
| GET  | `/grid/leaderboard` | public | Ranked contributors (see §6.4). |
| GET  | `/grid/work/{ena_acc}` | public | Canonical community result for a sample. |
| GET  | `/grid/stats` | public | Grid throughput / units-remaining / active nodes. |

### 4.5 Curating the work list (`du-jobs`)
A new interval job (or `run-once` backfill) turns ENA metadata into `grid.work_unit` rows:
- Enumerate candidate samples from `pubs.genomic_study` (source `ENA`) and/or a curated study
  allow-list.
- For each run, hit the **ENA Portal `filereport`** endpoint:
  `.../filereport?accession=<run>&result=read_run&fields=run_accession,sample_accession,
  fastq_ftp,fastq_bytes,fastq_md5,submitted_ftp,submitted_md5,library_layout,instrument_platform,
  instrument_model,read_count,base_count&format=tsv`.
- **P1 curation:** only emit units where `submitted_ftp` (or an ENA analysis object) exposes a
  **CRAM/BAM** → `data_kind = CRAM`, skip-align path. **P2:** open FASTQ-only samples
  (`data_kind = FASTQ`) once realignment is wired.
- Store the run manifest (URLs + md5 + bytes + layout + inferred `read_type`) on the unit; set
  `est_bases`/`est_download_bytes` for weighting and preflight.

---

## 5. Results — records, provenance, and the digest

### 5.1 Subject ≠ contributor (new)
Grid results are *about* an ownerless public sample but *computed by* the contributor. The fed
records gain (in `du-domain::fed`, additive, back-compatible via `#[serde(default)]`):

```rust
// New shared block, embedded on grid-produced records.
pub struct Provenance {
    pub computed_by: String,       // contributor did:plc / did:key
    pub software: String,          // "navigator"
    pub stack_version: String,     // analysis stack semver (matches grid.submission.stack_version)
    pub reference_build: String,   // "chm13v2.0"
    pub aligner: Option<String>,   // "minimap2-rs <ver>/<preset>" | None (passthrough)
    pub source: String,            // "ena:read_run" — how the input was obtained
}
```

- **Subject** = the ENA accession, carried as an `ExternalId { namespace: "ENA", value: "SAMEA…" }`
  on `BiosampleRecord` (that field already exists). This is what makes the record *about* the
  public sample rather than about the contributor's own genome.
- **Contributor** = the publishing repo DID (as today) **plus** the explicit `provenance.computed_by`
  so the AppView can attribute credit even after canonicalization pools across contributors.
- Records published: `BiosampleRecord` (anchor, with the ENA `ExternalId`), `SequenceRunRecord`(s),
  `AlignmentRecord` (coverage), `PopulationBreakdownRecord` (ancestry), and haplogroup calls — the
  **same builders Navigator already has** (`publish.rs`), extended with `Provenance`.
- **Community tier.** These records are marked community/public (never personal-genome tier). The
  contributor's PDS holds them, but the AppView's canonical copy is the citable community asset.

### 5.2 The result digest (what validation compares)
Realignment is **not byte-deterministic** across thread counts/builds, so we never hash the BAM. We
compare a small canonical digest of **discrete calls** plus **bucketed** continuous metrics:

```jsonc
{
  "unit": "SAMEA0000000",
  "reference_build": "chm13v2.0",
  "stack_version": "1.7.0",
  "aligner": "minimap2-rs 2.28/sr",       // null for CRAM passthrough
  "calls": {
    "sex": "XY",
    "y_terminal": "R-FGC29071",           // exact match required
    "mt_terminal": "U5a1b1g",             // exact match required
    "ancestry_superpop_argmax": "EUR",    // exact match required
    "coverage_mean_bucket": 30,           // bucketed (e.g. round to nearest 2×) — float drift tolerant
    "callable_fraction_bucket": 0.94      // bucketed to 2 decimals
  }
}
```

- **Digest is signed** with the device key (`grid.submission.digest_sig`) — the same
  `verify_did_key` path proves *this node* produced *this digest*.
- **Comparison rule:** two digests **agree** iff all discrete calls match exactly and every bucketed
  metric matches its bucket. Only digests with **compatible** `(reference_build, stack_version-major)`
  are compared; a stack-major bump can re-open units (define a compatibility window per metric).
- Continuous fields stay *out* of the agreement test but are still published in the full records for
  downstream use; only the digest gates canonicalization.

---

## 6. Validation, credit, reputation, leaderboard

### 6.1 Adaptive replication (the validator job, `du-jobs`)
Trust tiers derived from the contributor's grid history (not social score alone):

| Tier | Entry condition | Replication policy |
|------|-----------------|--------------------|
| **Untrusted** | new node / < N agreed units | Submissions only *contribute to* quorum; never canonical alone. Unit needs ≥2 agreeing digests from distinct DIDs. |
| **Provisional** | ≥ N agreed, 0 recent divergence | Quorum = 2, but its agreement can pair with one Untrusted to canonicalize. |
| **Trusted** | ≥ M agreed, sustained agreement | A single submission canonicalizes; ~5 % of units randomly flagged for a shadow replica (spot-recheck). |

Validator loop, per `SUBMITTED` unit:
1. Gather `grid.submission` digests for the unit.
2. Cluster by agreement (§5.2). If a cluster meets the unit's `required_replicas` **and** trust
   policy → set unit `CANONICAL`, store `canonical_digest`, mark those submissions `AGREED`, credit
   their contributors (§6.3).
3. If clusters conflict (divergence) → unit `CONTESTED`, bump `required_replicas`, re-offer to a
   fresh node; mark minority submissions `DIVERGENT` and apply the divergence penalty (§6.2).
4. Trusted-node single-run canonicalized units: with 5 % probability, still re-offer once for a
   shadow check; a contradicting shadow flips to `CONTESTED`.

### 6.2 Anti-abuse
- **Sybil resistance.** A node's DID resolves to a `ident.users` account via `fed.device_key`;
  Untrusted submissions can't self-canonicalize, so a lone attacker can't inject canonical results.
  Rate-limit `claim` per user; cap concurrent leases per user.
- **Poisoning.** Quorum over signed digests + divergence penalties. A `DIVERGENT` submission costs
  reputation (`SPAM_REPORT_VALIDATED`-style negative event) and demotes the node's tier; repeated
  divergence → cooldown.
- **Free-riding / duplicate submit.** One credit per `(user, unit)` (unique index); resubmitting a
  known canonical digest without independent compute earns nothing (can't beat the first-submitter
  timestamp, and shadow re-checks are AppView-chosen, not self-selected).
- **Replay.** `ensure_fresh_ts` ±300 s on every signed call.
- **ENA fair-use.** Respect ENA endpoint etiquette (md5-verify, resumable, bounded concurrency);
  the work list is curated centrally so nodes don't hammer ENA discovering files.

### 6.3 Credit formula (cobblestones)
Credit ∝ work magnitude, awarded **only** on `AGREED`/canonical:
```
cobblestones = base
             + realign_factor * (est_bases / 1e9)     // P2 realign: paid per Gbp mapped
             + analysis_factor                          // fixed for running the full stack
first-to-canonical  → CANONICAL_FIRST bonus
quorum agreement    → QUORUM_AGREE (full)
shadow spot-check   → SPOTCHECK_PASS (small)
divergent           → 0 (+ reputation penalty)
```
P1 (passthrough, no realign) pays `base + analysis_factor` — lighter, reflecting the smaller compute.

### 6.4 Reputation & leaderboard
- **Compute-credit leaderboard** — the primary artifact of the ask. `grid.credit` summed per user:
  ```sql
  SELECT u.handle, SUM(c.cobblestones) AS score, COUNT(*) AS units
  FROM grid.credit c JOIN ident.users u ON u.id = c.user_id
  GROUP BY u.id, u.handle ORDER BY score DESC LIMIT 100;
  ```
  Served at `GET /api/v1/grid/leaderboard` (public), with all-time / rolling-30d windows.
- **Reputation bump** — each canonical contribution also fires a **capped** `WORK_UNIT_COMPLETED`
  `social.reputation_event` (seeded like `0042_reputation_seed.sql`) so grid work has the requested
  "positive impact on reputation" — but capped/diminishing so it can't be farmed to blow past social
  gates (DM/GROUP/RECRUIT thresholds). Compute standing lives mainly in the `grid.credit` board;
  reputation gets a bounded, honest boost.

---

## 7. Navigator — the edge worker

### 7.1 Placement (respects `ui → app → {analysis, store, sync, refgenome} → {domain, du-*}`)
- **`navigator-analysis`** — realignment (`align`/revert modules) per
  [`realignment-module.md`](realignment-module.md). Pure compute; no new upward deps.
- **`navigator-app`** — new modules:
  - `ena.rs` — ENA fetch client (Portal `filereport` already resolved server-side; the node just
    downloads the manifest's URLs). Resilient/resumable/md5-verified, mirroring `refgenome::download`.
  - `grid.rs` — the coordination client (register/claim/heartbeat/submit/release) using the
    device-key-signed pattern from `ibd_exchange.rs`, and the **per-unit driver**.
  - Lift `run_full_analysis_streaming` (currently in `navigator-ui/src/worker.rs`) into an
    app-level `App::run_full_analysis(alignment_id, progress)` so the same sequence runs headless.
    This is a prerequisite refactor (the worker keeps calling the lifted function).
- **`navigator-ui`** — a "Contribute / Grid" panel (claimed units, per-stage progress, credits,
  leaderboard rank, pause/resume, resource budget) + `cli.rs` `contribute` subcommand.

Note the crate-graph reality: `grid.rs` needs analysis + refgenome + sync + store — all already
below `app`, so it lives *in* `app`. (A dedicated `navigator-grid` crate is an option if the surface
grows, but starting in-app matches how `import_unified`/`ibd_exchange`/`sync` already live there.)

### 7.2 Per-unit driver (headless, cancellable)
```
claim unit → preflight (disk budget for download + scratch + output; refuse early)
  → download manifest files from ENA (resumable, md5-verify)      [ena.rs]
  → P1  CRAM/BAM: index if needed; register Alignment on its stored build
     P2  FASTQ:    minimap2-rs → CHM13 → sort/markdup/CRAM → register Alignment
                    (realignment-module.md Stages B–D; no revert — inputs are already unaligned)
  → App::run_full_analysis(alignment_id)  → coverage/sex/read_metrics/SV/Y/mt
  → estimate_ancestry_from_consensus (biosample level)
  → build result DIGEST; sign with DeviceKey
  → build fed records (Provenance{computed_by=self, …}, ExternalId ENA=<acc>) → publish_* (sync_outbox)
  → POST /grid/submit (digest + digest_sig + record at:// refs) → lease COMPLETED
  → clean scratch dir; loop to next unit
```
- Reuses the existing streaming/cancellable spawn-loop discipline (progress events, honor
  `CancelAnalysis`, `await` between stages).
- **Resource governance:** `--max-units`, `--lease-days`, `--max-disk`, data-kind filter, threads
  (`NAVIGATOR_ANALYSIS_THREADS` / `NAVIGATOR_REALIGN_THREADS`), pause on AC/thermal (nice-to-have).
- **CLI:** `navigator contribute --data-kind cram --max-units 4 --lease-days 3 --max-disk 200G`.

### 7.3 Platform reality
P1 (passthrough) runs everywhere Navigator does. P2 (realign) inherits the realignment module's
**macOS + Linux (incl. Apple Silicon via `simde`)** target; Windows nodes can still contribute in P1
(passthrough) and get FASTQ realignment once the Windows FFI spike lands. `/grid/claim` filters by
the node's advertised capabilities so Windows nodes are simply never offered FASTQ units.

---

## 8. Phasing / milestones

| Phase | Deliverable | Proves |
|-------|-------------|--------|
| **P0** | Shared: `Provenance` block + subject/`computedBy` on records; `messages::grid` canonical strings; AppView `grid` schema (`0059`) + reaper; `du-jobs` ENA curation (CRAM-only). | Wire contracts + coordination substrate. |
| **P1** | **CRAM passthrough, end-to-end.** Navigator `ena.rs` + `grid.rs` + lifted `run_full_analysis`; `contribute` CLI; register/claim/heartbeat/submit/release; validator (adaptive replication) + `grid.credit` + `/grid/leaderboard`. No aligner. | The **whole distributed loop** (lease→compute→submit→validate→canonical→credit→board) with zero aligner risk. |
| **P2** | **FASTQ → minimap2 realign.** Wire the realignment engine ([`realignment-module.md`](realignment-module.md)) into the driver for `data_kind=FASTQ`; open FASTQ curation; per-Gbp credit. | The real vision — uniform hs1 realignment of arbitrary ENA reads. |
| **P3** | GUI Grid panel (progress, credits, rank, budget); rolling leaderboards; public `/grid/work/{acc}` result pages; grid-wide stats. | Community-facing polish + the visible leaderboard. |
| **P4** | Hardening: trust-tier tuning, divergence-penalty calibration, spot-check rate tuning, ENA fair-use throttles, Windows FASTQ (realignment P5 spike). | Robustness at scale. |

---

## 9. Open questions

- **Work-unit granularity** — per ENA *sample* (merge runs; matches per-biosample analysis) vs per
  *run* (finer leases, simpler downloads, but multiple runs per sample need re-merging for
  consensus haplogroups). Leaning sample-level; runs listed in the manifest.
- **Target reference** — `Chm13v2` vs the analysis-tuned `Chm13v2MaskedRcrs` (PAR-masked + rCRS).
  Must match whatever the ancestry/IBD panels are built against; realignment-module.md flags the
  same question. The digest's `reference_build` must pin the exact choice.
- **Stack-version compatibility window** — when does a stack bump *invalidate* an existing canonical
  digest vs. remain comparable? Per-metric policy (e.g. haplogroup tree version matters; a coverage
  refactor may not). Needs a compatibility matrix.
- **Contributor PDS storage cost** — publishing community records into a volunteer's own repo grows
  their PDS. Acceptable, or should the contributor publish only a lightweight *attestation* while
  the AppView holds the full canonical records? (D4 says both; revisit if repos bloat.)
- **ENA data governance** — public consented research data, but confirm per-study data-use notes;
  mark provenance so downstream consumers can honor any study-specific terms.
- **Consensus across contributors** — when two contributors produce the *records* (not just digests)
  for one canonical unit, which record set becomes canonical? Propose: first-to-canonical's records,
  with the others retained as `AGREED` corroboration (and credited).
- **Credit calibration** — cobblestone constants (`realign_factor`, bonuses) and the reputation cap
  need real throughput data before they're fair; ship P1 with conservative placeholders.

---

## 10. Cross-references
- [`realignment-module.md`](realignment-module.md) — the minimap2-rs realignment engine this Grid
  consumes for `data_kind=FASTQ` (Stages B–D; revert is skipped since ENA FASTQ is already unaligned).
- [`academic-ena-import.md`](academic-ena-import.md) — the (design-only) single-sample ENA import
  path; the Grid generalizes its ENA-fetch idea to a coordinated fleet.
- AppView `migrations/0008_fed.sql` (`fed.pds_*`) — the dormant lease/node/submission scaffold reused
  here; `du-web/src/sig.rs` + `fed.device_key` — the reused edge-auth primitive;
  `social.reputation_event` — the reused credit ledger.
- `du-domain::fed` — the wire records extended with `Provenance` + subject/`computedBy`.
