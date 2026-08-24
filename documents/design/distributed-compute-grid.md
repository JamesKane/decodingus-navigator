# DecodingUs Grid — distributed community realignment & analysis

Status: **design / specification only** (no code — re-verified by grep 2026-08-24). Cross-repo:
**Navigator** (edge worker) + **AppView** (`decodingus`, coordinator) + **shared**
(`decodingus-shared`, wire records).

> **Read [§11](#11-reconnaissance-refresh-2026-08-24) before costing any of this, and
> [§12](#12-what-was-built-and-what-the-build-changed) for what has since been built and the three
> places where building it changed the design.** The doc was
> drafted while realignment was still a plan. Realignment has since shipped, and it shipped on a
> **different aligner backend than D1 locks in** — which is why D1 and §7.3 below are struck and
> corrected. Three of the four §2 "greenfield" items also moved. The design's shape survives the
> refresh intact; its estimates do not.

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
| D1 | **Aligner integration** | ~~minimap2 via `minimap2-rs` FFI (`static` + `simde`)~~ → **`minimap2-pure-rs` (pure Rust)** — corrected 2026-08-24, [§11](#11-reconnaissance-refresh-2026-08-24) | The rationale held and the answer changed under it. This row said the Grid *consumes* the realignment module's engine rather than re-deciding it — correct, and that module shipped on a **pure-Rust translation of minimap2 v2.31**, not an FFI binding. No C toolchain, so every Rust target builds, **Windows included**. Measured 99.74 % byte-identical to the C implementation, with zero disagreements at MAPQ > 0. |
| D2 | **Trust model** | **Adaptive replication** | Untrusted nodes run in shadow/quorum; reputation graduates them to trusted single-run + random spot-recheck. BOINC-proven; K× cost only where trust is unearned. |
| D3 | **First cut** | ~~Staged — CRAM-passthrough first~~ → **both data kinds in P1** (amended 2026-08-24, [§12](#12-what-was-built-and-what-the-build-changed)) | The staging existed to retire aligner risk before the coordination loop. That risk evaporated when realignment shipped on a pure-Rust mapper (D1), so the reason for the stage went with it. P1 now claims both CRAM (passthrough) and FASTQ (realign), and `data_kind` is on the work unit from the first migration rather than bolted on later. |
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
- **`du-jobs` job runner** — ~~in-process interval jobs (`scheduler.rs`)~~ **that scheduler is
  retired** (found 2026-08-24): it "fired every job on startup and let same-period jobs overlap,
  spiking DB + memory". Jobs now run as `run-once <job>` under systemd timers, serialized by one
  Postgres advisory lock (`DU_JOBS_LOCK`). The lease-reaper and the validator are therefore
  `run-once` jobs, not loops — a better fit, since both are naturally batch.

### Reuse — Navigator (`DUNavigator`)
- **Full analysis stack, per-alignment** — `navigator-app/src/analysis.rs`: `run_unified_metrics`
  (coverage+read_metrics+sex), `run_sv`, `run_denovo_caller`, + `haplogroup.rs`
  (`assign_y_haplogroup`, `assign_mtdna_haplogroup_from_alignment`, `place_{y,mt}_consensus`,
  `estimate_ancestry_from_consensus`). Each computes + persists a versioned artifact.
- **The realignment engine — SHIPPED, no longer merely specified.**
  [`realignment-module.md`](realignment-module.md) merged as `bf576ab` and released in
  `v0.1.0-alpha.17`: revert → map → sort / duplicate-mark / CRAM → register, with the aligner-index
  cache and the per-technology presets, in the `navigator-align` crate. **The Grid is that module's
  first heavy consumer** — and it is no longer waiting on it. This is what unblocked the Grid.
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
1. **ENA sequence fetch — the AppView half is largely built already** (corrected 2026-08-24, §11).
   `EnaClient::run_files` resolves the run-level `filereport` and returns exactly the fields a work
   unit needs; `du-jobs/crawl_project.rs` groups those runs by sample and materializes the files. The
   genuinely missing piece is the **Navigator-side downloader** — resumable and **md5**-verified,
   which is *not* what `refgenome::download` gives you (SHA-256, no `Range`, one blind retry).
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

> **The unit state machine above is superseded** — see
> [§12.1](#121-claimability-is-derived-the-state-machine-above-cannot-express-replication).
> `LEASED` and `SUBMITTED` cannot coexist with `required_replicas > 1`, which the sentence directly
> above this note requires. What shipped keeps only the exclusive milestones as states and derives
> claimability. The diagram is kept because the *lifecycle* it draws is still right.

---

## 4. AppView — data model & coordination

New Postgres schema **`grid`** (migration `0075_grid.sql` — next in sequence as of 2026-08-24; the
doc originally guessed `0059`). Reuse `fed.pds_node`
+ `fed.device_key`; everything work-specific is new so we don't overload `fed.pds_submission`'s
existing semantics.

### 4.1 Tables (sketch)

> **Superseded by `rust/migrations/0075_grid.sql` in the `decodingus` repo**, which is the truth:
> applied migrations are checksummed by `sqlx::migrate!`, so the SQL that ran is the SQL that is.
> The sketch below is kept as the intent. Two things came out different — the state model
> ([§12.1](#121-claimability-is-derived-the-state-machine-above-cannot-express-replication)) and
> the credit units ([§12.2](#122-the-credit-ledger-is-an-integer)).

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

> **Built as `du_db::grid::claim`.** The sketch below is close but wrong in one way that matters:
> it flips the unit to `LEASED`, which makes the unit unclaimable by the *second* replica the same
> design requires. What shipped never mutates the unit on claim —
> [§12.1](#121-claimability-is-derived-the-state-machine-above-cannot-express-replication).

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
- **P1 curation (amended 2026-08-24):** emit **both** kinds. `submitted_ftp` exposing a CRAM/BAM →
  `data_kind = CRAM` (skip-align); FASTQ-only samples → `data_kind = FASTQ` (realign). The original
  text deferred FASTQ to P2; D3's staging is gone, so the split is now only a per-unit label, and
  the node's advertised kinds decide what it is offered.
- Much of this is **already built AppView-side** — `EnaClient::run_files` and
  `du-jobs/crawl_project.rs`. See [§11.5](#115-the-appview-already-resolves-ena-at-run-level-and-already-curates).
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
    // "mt_terminal" was here. REMOVED 2026-08-24 — see the note below and §12.3.
    "ancestry_superpop_argmax": "EUR",    // exact match required
    "coverage_mean_bucket": 30,           // bucketed (e.g. round to nearest 2×) — float drift tolerant
    "callable_fraction_bucket": 0.94      // bucketed to 2 decimals
  }
}
```

> **`mt_terminal` was a problem, not a field — and it is now DECIDED (2026-08-24).** It demanded an
> exact match, but the analysis path deliberately does not produce it: `App::analyze_biosample`
> states in its own doc comment that it "does not assign mtDNA, by design. That value is not final
> on CHM13." Since the whole Grid realigns *to* CHM13, that was not a wiring gap to route around.
> **Resolution: mt is out of the agreement test**, exactly as §5.2 already treats continuous
> fields — still published in the full records, simply not gating canonicalization. `y_terminal`,
> sex and the ancestry argmax carry the discrete signal. See [§12.3](#123-mt_terminal-is-out-of-the-agreement-test).

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

> The validator is a **`run-once` job under a systemd timer**, not an in-process loop — the
> interval scheduler this section assumed is retired (§2). Batch suits it: the loop below is
> already written as "per SUBMITTED unit", which is a pass, not a daemon.
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

> **Ledger units (built 2026-08-24):** `grid.credit.cobblestones_milli` is an exact integer in
> **thousandths**, and `du_db::grid::COBBLESTONE = 1_000`. `NUMERIC` would need a decimal feature
> this workspace's sqlx is not built with, and `f64` is the wrong shape for a column that gets
> SUMmed over every contribution ever made. Three decimals was the intended precision anyway, so
> the integer *is* the value. Divide at the point of rendering, never before —
> [§12.2](#122-the-credit-ledger-is-an-integer).

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
A **passthrough** unit pays `base + analysis_factor` — lighter, reflecting the smaller compute. (This
read "P1 … pays" when P1 was passthrough-only; since D3 was amended it is a property of the unit's
`data_kind`, not of the phase.)

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
  - ~~Lift `run_full_analysis_streaming` into an app-level `App::run_full_analysis`~~ — **this
    prerequisite mostly landed already** (2026-08-24, §11). `App::analyze_biosample`
    (`queries.rs:975`) is the headless, cancellable, resumable, preflighted **unit of work for one
    sample**, built for batch analysis in PR #47. The driver calls it, plus
    `estimate_ancestry_from_consensus` for the biosample-level ancestry. It covers coverage, Y, sex
    and read metrics — **not** mtDNA and **not** SV, both by deliberate choice, which §11 unpacks.
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
  → App::analyze_biosample(biosample)     → coverage/sex/read_metrics/Y  (NOT SV — see §11)
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

### 7.3 Platform reality (rewritten 2026-08-24 — the constraint this section was built on is gone)

**Both phases run everywhere Navigator does, Windows included.** The original text below the fold
assumed the aligner reached the node through a C FFI binding, which is why it split the fleet by
operating system. `navigator-align` ships a **pure-Rust** mapper (D1), so there is no C toolchain,
no `simde` target list, and no Windows spike to wait for. A Windows node can take FASTQ units on the
day P2 opens.

`/grid/claim` should still filter on advertised capabilities — but on **RAM, free disk and thread
count**, which is what actually decides whether a node can finish a 30× WGS unit. Not on OS.

> ~~*Superseded:* "P2 (realign) inherits the realignment module's macOS + Linux (incl. Apple Silicon
> via `simde`) target; Windows nodes can still contribute in P1 (passthrough) and get FASTQ
> realignment once the Windows FFI spike lands."~~

---

## 8. Phasing / milestones

| Phase | Deliverable | Proves |
|-------|-------------|--------|
| **P0** | Shared: `Provenance` block + subject/`computedBy` on records; `messages::grid` canonical strings; AppView `grid` schema (`0075`) + reaper; `du-jobs` ENA curation (both kinds). **Partly built — see [§12](#12-what-was-built-and-what-the-build-changed).** | Wire contracts + coordination substrate. |
| **P1** | **The whole loop, both data kinds** (amended 2026-08-24 — P2 folded in, D3). Navigator `ena.rs` + `grid.rs` + the driver over `App::analyze_biosample`; `contribute` CLI; register/claim/heartbeat/submit/release; validator (adaptive replication) + `grid.credit` + `/grid/leaderboard`; the realignment engine wired for `data_kind = FASTQ`; per-Gbp credit live. | lease→compute→submit→validate→canonical→credit→board, **and** the uniform-hs1 payoff. |
| ~~**P2**~~ | *Folded into P1.* Kept as a heading so existing references resolve. The aligner risk that justified staging it separately no longer exists ([§11.1](#111-the-aligner-is-pure-rust-so-the-fleet-is-not-split-by-os)). | — |
| **P3** | GUI Grid panel (progress, credits, rank, budget); rolling leaderboards; public `/grid/work/{acc}` result pages; grid-wide stats. | Community-facing polish + the visible leaderboard. |
| **P4** | Hardening: trust-tier tuning, divergence-penalty calibration, spot-check rate tuning, ENA fair-use throttles. (~~Windows FASTQ~~ — no longer a milestone; the pure-Rust mapper made it free. See §11.) | Robustness at scale. |

---

## 9. Open questions

- ~~**Work-unit granularity**~~ — **RESOLVED 2026-08-24: per ENA *sample***, with the runs listed in
  the unit's manifest. It is the grain Navigator already analyses at (`App::analyze_biosample`),
  the grain consensus haplogroups need, and the grain `du-jobs/crawl_project.rs` already groups ENA
  runs into. A multi-run sample is one lease; merging its runs is work the node must do anyway
  before consensus means anything. `grid.work_unit.sample_accession` is the unit's identity.
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

---

## 11. Reconnaissance refresh (2026-08-24)

§2 was written against the three repos as they stood when the Grid was still blocked on
realignment. Realignment shipped on 2026-08-14, and this section re-walks §2's claims against the
tree. **Every claim below was checked by grep, not by reading a status header.**

The headline: the design's *shape* survives intact — staged phases, adaptive replication, digest
comparison, the lease state machine. What moved is the **cost**, and it moved down. One locked
decision was overtaken, one prerequisite refactor turned out to be mostly built, and one greenfield
item turned out to be half-built in the AppView.

### 11.1 The aligner is pure Rust, so the fleet is not split by OS

D1 locked "minimap2 via `minimap2-rs` FFI (`static` + `simde`)". `navigator-align` instead uses
**`minimap2-pure-rs`**, a pure-Rust translation of minimap2 v2.31 — the crate's own module doc says
it "does not link the C library through FFI. It needs no C toolchain, so Windows and every other
Rust target build unchanged," with a parity measurement of **99.74 % byte-identical output and zero
disagreements at MAPQ > 0**.

This is the single largest correction in this refresh, because §7.3 built a two-tier fleet on top of
the assumption it contradicts:

| Assumed | Actual |
|---|---|
| P2 realign is macOS + Linux only | P2 realign runs on every Rust target |
| Windows nodes are P1-only until an FFI spike lands | Windows nodes take FASTQ units the day P2 opens |
| `/grid/claim` filters by OS capability | `/grid/claim` should filter on **RAM / disk / threads** |
| P4 carries "Windows FASTQ (realignment P5 spike)" | That milestone does not exist |

§7.3 is rewritten and the P4 row is corrected in place.

### 11.2 The §7.1 prerequisite refactor is mostly done — PR #47 did it for other reasons

§7.1 called for lifting `run_full_analysis_streaming` out of `navigator-ui/src/worker.rs` into the
app so the sequence could run headless. That function is still in `worker.rs`, but the Grid no
longer needs it lifted: **`App::analyze_biosample` (`queries.rs:975`) is already that unit of work**,
built in `331e8cb` (PR #47, "Make batch analysis tractable"). Its own doc comment calls it "the unit
of work for one sample. The project pass and the deep-analyze job both call it." It is headless,
cancellable mid-sample, resumable (skips what the store already holds), and runs a preflight before
any I/O so a batch does not discover a file problem the slow way — all properties the per-unit driver
in §7.2 would otherwise have had to grow itself.

**But read what it deliberately omits, because §5.2 depends on one of them.**

### 11.3 `mt_terminal` in the digest contradicts the analysis path

`analyze_biosample` states: *"The method does not assign mtDNA, by design. That value is not final on
CHM13. See the notes on the reconciliation and the liftover."* §5.2's digest lists `mt_terminal` as a
field requiring an **exact match** for two submissions to agree.

The Grid realigns everything **to** CHM13. So this is not a wiring gap that the driver can route
around by calling one more method — the value the digest wants to compare is the value the analysis
path declines to state on this reference. Three ways out, in the order I would consider them:

1. **Drop `mt_terminal` from the agreement test** and publish it in the full records only, exactly as
   §5.2 already does for continuous fields. Cheapest, and loses little: `y_terminal`, sex, and the
   ancestry argmax already carry the discrete signal.
2. **Resolve CHM13 mt placement first**, and make it a P1 prerequisite rather than a Grid concern.
   Correct, but it puts a research question on the critical path of an infrastructure milestone.
3. **Pin the digest to `Chm13v2MaskedRcrs`** (the rCRS-masked analysis reference), which is the same
   fork §9 already flags as an open question under "Target reference". If mt is to stay in the
   digest, this is the coherent way — and it makes §9's question load-bearing rather than academic.

**Recommendation: (1) for P1, and let §9's reference question settle on its own timeline.** Shipping
the coordination loop should not wait on an mtDNA placement decision. — **Taken, 2026-08-24. See
[§12.3](#123-mt_terminal-is-out-of-the-agreement-test); the field is gone from §5.2's digest.**

### 11.4 SV does not belong in the per-unit driver

§7.2's driver line listed SV. `analyze_biosample` excludes it, and `CLAUDE.md` is explicit that SV is
opt-in and never automatic because it walks every read for its own sake: **2–5 h per WGS sample,
against ~1 h for the entire rest of the stack**. Putting it in the driver roughly triples what a
volunteer pays per unit — to compute a signal the digest never compares. Removed from §7.2.

### 11.5 The AppView already resolves ENA at run level, and already curates

§2 greenfield #1 said "AppView pulls study metadata only". That has not been true for some time:

- **`EnaClient::run_files(accession)`** (`du-external/src/ena.rs`) enumerates every run for a study
  or BioProject through the ENA `filereport` endpoint, requesting exactly the columns a work unit
  needs: `submitted_ftp`, `submitted_md5`, `submitted_bytes`, `submitted_format`, `fastq_ftp`,
  `fastq_md5`, `fastq_bytes`, `instrument_model`, `library_layout`, `read_count`, `first_public`.
  Parsing is header-driven, so it tolerates ENA reordering fields.
- **`du-jobs/src/crawl_project.rs`** already groups those runs by sample, upserts each biosample with
  `source = EXTERNAL`, links it to the study's publications, and materializes the files through
  `du_db::sequence::ingest_libraries` — idempotent at sample granularity. It already separates
  `ALIGNED = [BAM, CRAM]` from `INDEX = [CRAI, BAI]`, **which is precisely D3's CRAM-passthrough
  selector**.

So §4.5's "curating the work list" is closer to a **query over `genomics.sequence_file` plus a
`grid.work_unit` projection** than to a new ENA integration. The politeness discipline §6.2 asks for
(bounded batches, a request gap) is already the crawl job's practice.

### 11.6 The Navigator downloader is a real build, not a copy

§7.1 said `ena.rs` should mirror `refgenome::download`. It cannot, quite. That function streams to a
`.part` file, computes **SHA-256** inline so the caller needs no second read, and retries exactly
once on a transient error — but it sends **no `Range` header, so it cannot resume**, and ENA
publishes **md5**. On multi-GB run files over a volunteer's connection, resume is the feature that
decides whether a unit ever completes. Budget `ena.rs` as new code: resumable, md5-verified,
bounded-concurrency. The `.part` + rename discipline is worth copying; the hash and the retry are not.

### 11.7 Everything else in §2 still holds

Re-verified unchanged:

- **The `fed.pds_*` lease/node/submission scaffold is still dormant** — the only references anywhere
  in the AppView are a migrations test and the dedup repointer. Still available, still unwired.
- **`Provenance` does not exist** in `du-domain` (§5.1 remains fully greenfield). `ExternalId
  { namespace, value }` does exist and is ready to carry the ENA accession.
- **No `grid` schema exists.** The AppView is at migration **`0074`**, so the Grid schema is `0075`,
  not the `0059` §8 guessed.
- **Navigator's signing and publishing spine is intact** — `navigator-sync::DeviceKey`, the
  `exchange_get_poll` signed-GET pattern (`ibd_exchange.rs:613`), `enqueue_publish` → `sync_outbox` →
  `drain_outbox`.
- **Edge auth got stronger, not weaker.** The AppView now has `verify_signed_fresh` alongside
  `verify_signed`/`ensure_fresh_ts`, plus a `signed_request_seen` table (migration `0067`) — a
  replay-**seen** guard, where §6.2 assumed only the ±300 s freshness window.

### Where this leaves the plan

**P0/P1 as written, with the corrections above folded in.** D3 staged CRAM-passthrough first to
retire aligner risk before the coordination loop — and that risk has since evaporated on its own, so
what P1 now buys is purely the value of proving lease → compute → submit → validate → canonical →
credit end to end on the cheapest possible unit. That is still worth doing first.

**The critical path is the AppView half**, not Navigator's. `grid.work_unit` + the
`FOR UPDATE SKIP LOCKED` claim + the validator job are the parts with no existing analogue. On the
Navigator side, P1 is now roughly: a signed HTTP client (copy `exchange_get_poll`), a resumable md5
downloader (§11.6), and a driver that calls `analyze_biosample` and
`estimate_ancestry_from_consensus` and signs a digest.

~~**One decision to settle before writing P1 code:** `mt_terminal` in the digest (§11.3).~~
**Settled 2026-08-24**, along with two others that the build then surfaced. What that decision was,
and what has since been built on it, is [§12](#12-what-was-built-and-what-the-build-changed) — read
that for the current state; this section is the reconnaissance that preceded it.

---

## 12. What was built, and what the build changed (2026-08-24)

The first increment is the AppView coordination substrate — the half §11 identified as the critical
path, because it is the half with no existing analogue anywhere in the three repos.

**Landed:**

| Artefact | What it is |
|---|---|
| `rust/migrations/0075_grid.sql` | `grid.work_unit` / `lease` / `submission` / `credit`. Reuses `fed.pds_node` (registry) and `fed.device_key` (auth); does **not** reuse `fed.pds_submission`, whose status lifecycle means curator review of a proposed call — a different thing from digest quorum, and overloading it would make both unreadable. |
| `du-db/src/grid.rs` | `claim` · `heartbeat` · `release` · `reap_expired` · `submit` · `award_credit` · `leaderboard` · `register_node`, plus `messages` — the canonical signed strings. |
| `du-db/tests/grid.rs` | Nine live-Postgres tests: replica bounds, self-replication, the data-kind filter, reclamation, submit/resubmit, credit idempotence, and the three curation cases. |
| `du-db::grid::curation_candidates` + `du-jobs/src/grid_curate.rs` | The `run-once grid-curate` job (§4.5), projecting crawled samples into the work list. |

Building it settled three things the design had left ambiguous or wrong. Each is recorded in the
migration header as well, because `sqlx::migrate!` checksums applied migrations — the SQL that ran
can never be edited, so its comments are the one explanation that cannot drift from it.

### 12.1 Claimability is derived; the state machine above cannot express replication

§3 gives the work unit `AVAILABLE → LEASED → SUBMITTED → CANONICAL`, and §4.2's sketch flips the
unit to `LEASED` on claim. Two paragraphs after that diagram, the same document says
`required_replicas` defaults to **2** — so a unit routinely needs a second independent result while
a first node still holds a lease. Under the sketch, that second node can never get it: the unit is
no longer `AVAILABLE`.

`LEASED` and `SUBMITTED` would each have to mean "…and also still claimable", which is not a state.
So what shipped keeps on `work_unit.state` only the milestones that are genuinely exclusive —
`AVAILABLE` / `CANONICAL` / `CONTESTED` / `RETIRED` — and **derives** claimability:

```
claimable  ⇔  state IN ('AVAILABLE','CONTESTED')
              AND (active leases + non-divergent submissions) < required_replicas
              AND the calling DID holds no lease and no submission on the unit
```

One `SELECT … FOR UPDATE SKIP LOCKED` answers that, which is what §4.2 wanted in the first place.
The lease and submission tables are the source of truth for how many replicas are in flight, so no
counter is maintained and no counter can drift from the rows it summarises. A claim never mutates
the work unit at all.

The third clause is a **correctness** rule, not a rate limit: quorum means *independent* results, so
a contributor must never be handed a unit it already holds or has already answered.

### 12.2 The credit ledger is an integer

`grid.credit.cobblestones_milli BIGINT`, with `du_db::grid::COBBLESTONE = 1_000`, rather than the
`NUMERIC` §4.1 sketched. This workspace builds sqlx without any decimal feature and nothing else in
the repo uses one, so `NUMERIC` has no Rust mapping — and adding a workspace-wide dependency for a
single column is out of proportion to what it buys. `f64` was the other option and is the wrong
shape for a column summed over every contribution ever made. Since `NUMERIC(12,3)` had already
chosen three decimal places, the integer **is** the value at the intended precision, and it sums
exactly. Divide by 1000 when rendering, never before.

### 12.3 `mt_terminal` is out of the agreement test

Decided as §11.3 recommended: option (1). The digest keeps `sex`, `y_terminal` and
`ancestry_superpop_argmax` as exact-match fields; mtDNA is published in the full records and does
not gate canonicalization. This unblocks P1 without putting a research question — CHM13 mt
placement — on an infrastructure milestone's critical path, and it leaves §9's `Chm13v2` vs
`Chm13v2MaskedRcrs` question free to settle on its own timeline.

### 12.4 Curation is a projection, and it exposed a hole in the credit formula

`run-once grid-curate` publishes the work list, and it makes **no network calls at all**. §11.5 was
right that the AppView already resolves ENA at run level: `crawl_project` has already grouped runs
by sample and written every file URL, md5 and size into `genomics.sequence_file`. Curation is one
query over tables we have. That also *is* the ENA fair-use control §6.2 asks for — a node receives a
finished manifest and never goes discovering files for itself.

`data_kind` is decided per sample and **the manifest is then filtered to match it**: a sample with
any CRAM/BAM is a passthrough unit carrying only aligned files, otherwise a FASTQ unit carrying only
reads. `build_libraries` already prefers aligned over FASTQ per sample, so the two normally agree —
deciding it again here means a manifest can never list a file the data kind says the node will not
open, and the download budget cannot be inflated by files nobody fetches.

**The hole: `est_bases` is `NULL` for essentially every unit.** It is `reads × read_length`, and the
crawl sets `read_length` to `None` — ENA's `filereport` does expose `base_count`, but
`du-external`'s `RUN_FIELDS` does not request it and `genomics.sequence_library` has nowhere to put
it. So the **per-Gbp term of the credit formula (§6.3) has nothing to weigh a FASTQ unit by**.

A byte-derived estimate was the tempting fix and is the wrong one: a fabricated number in a ledger
that pays people is worse than an honest null. So curation publishes the null and the job *warns*
with a count, rather than leaving a silent hole. **The fix, before credit goes live:** add
`base_count` to `RUN_FIELDS`, carry it through `crawl_project`, and store it — either as a column on
`sequence_library` or in the `atproto` provenance slot that already holds `run_accession`.

### What is NOT yet built

The `grid-reap` and `grid-validate` `run-once` jobs, the `/api/v1/grid/*` signed endpoints,
`Provenance` in `du-domain` (§5.1), and the whole Navigator edge (`ena.rs`, `grid.rs`, the driver,
the `contribute` CLI).

**One honest caveat on what did land:** the nine integration tests **compile but have not been run**.
The development host has no reachable Postgres — no local server, no Docker, and the Apple
`container` runtime's published port accepts TCP but resets on protocol traffic, the same
limitation that blocked live-PDS validation earlier. The claim SQL is therefore **unverified
against a real database**, which is exactly the part most worth verifying: `FOR UPDATE SKIP
LOCKED`, a partial-index `ON CONFLICT`, and replica arithmetic across two other tables only mean
anything inside a real transaction. The curation query is in the same position, and is if anything
more exposed: it reads JSONB paths (`http_locations->0->>'file_url'`,
`checksums->0->>'checksum'`) whose shape is defined only by `sequence::ingest_libraries` — which is
why the curation tests seed through that function rather than writing their own rows, so a test
cannot agree with the query while both disagree with the crawl. Two type ambiguities were removed
pre-emptively for the same reason (`LIMIT` takes a bigint; `bigint * interval` has no operator, so
`make_interval` is used instead).
**Run them before building anything on top:**

```
DATABASE_URL='postgres://…@localhost:5432/postgres?sslmode=disable' \
  cargo test -p du-db --test grid -- --nocapture
```
