# DUNavigator Backlog

Last reviewed: 2026-03-07

Source documents: `Edge_Client_Implementation_Status.md`, `UI_Redesign_Proposal.md`,
`HaplogroupReconciliation_Implementation.md`, `GenomeRegionsAPI_Specification.md`,
`UI_i18n_Guidelines.md`, `documents/atmosphere/` (multi-file Atmosphere Lexicon),
`Atmosphere_Lexicon_Alignment_Plan.md`, `YProfile_and_RegionAnnotations.md`

---

## Tier 1 — Partially Built (code exists, needs completion)

### 1.1 Ancestry Estimation
- **Source:** Edge_Client, UI_Redesign
- **Status:** Reference data pipeline stubbed; `AncestryEstimator`, `AncestryProcessor`, `AncestryReferenceCache` exist but need reference data deployed
- **Blocked on:** Genome Regions API / reference data hosting
- **Scope:** Wire reference download, validate estimator against known samples, build results UI

### 1.2 Chip Data to Haplogroup ✓
- **Source:** Edge_Client
- **Status:** Complete. `ChipHaplogroupAdapter` implemented and validated (Dec 2025).
- **Scope:** SNP-to-tree mapping for common chip panels (23andMe, AncestryDNA, FTDNA, LivingDNA)

### 1.3 Multi-Run Reconciliation — Phase 3
- **Source:** HaplogroupReconciliation_Implementation
- **Status:** Core algorithm implemented (branch compatibility + SNP concordance + cascade cleanup). Manual override UI not started.
- **Scope:**
  - [x] Branch compatibility scoring (LCA-based) — ancestor/descendant detection, divergence point, compatibility levels
  - [x] SNP-level conflict detection — concordance from supporting/conflicting counts, warnings for low concordance
  - [x] Identity verification warnings — flags incompatible haplogroups, suggests sample verification
  - [x] Cascade cleanup on run/profile deletion — SequenceDataManager + WorkspaceOperations wired
  - [ ] Manual override UI — accept/reject specific calls, override consensus

### 1.4 Genome Regions API (Server Deployment)
- **Source:** GenomeRegionsAPI_Specification
- **Status:** Client ready and feature-flagged (`GenomeRegionService` falls back to bundled `grch38.json`)
- **Scope:** Deploy server-side API; enable client to fetch region updates

### 1.5 Granular Record Sync
- **Source:** Edge_Client, Atmosphere_Lexicon
- **Status:** Sync queue tables exist (V003 migration), `SyncService`/`AsyncSyncService` partially wired
- **Scope:** Per-record PDS CRUD, conflict resolution, end-to-end sync flow

### 1.6 Atmosphere Lexicon Alignment
- **Source:** `documents/atmosphere/`, `Atmosphere_Lexicon_Alignment_Plan.md`
- **Status:** Complete (all phases A–D)
- **Scope:**
  - [x] Phase A: ChipProfile alignment (`vendor`→`provider`, new fields, V009 migration)
  - [x] Phase B: AncestryResult → PopulationBreakdown promotion (type renames, CI restructure, PopulationBreakdown model + entity + repository, V010 migration, SyncEntityType)
  - [x] Phase C: HaplogroupReconciliation enrichment (HeteroplasmyObservation, IdentityVerification, ManualOverride, AuditEntry types + audit logging in withRunCall/removeRunCall/recalculate, V011 migration, PdsClient codecs)
  - [x] Phase D: Sync-time validation & test type mapping (PdsSyncValidation module, SequenceRun.toSyncTestType, 14 tests)

---

## Tier 2 — Designed but Not Started

### 2.1 IBD Matching System
- **Source:** Edge_Client, UI_Redesign
- **Lexicon records:** `matchConsent`, `matchList`, `matchRequest`, `confirmedMatch`
- **Scope:** Consent management, match discovery, chromosome browser, relationship estimation
- **Depends on:** Ancestry estimation, granular sync
- **Implementation plan:** `documents/IBD_Matching_Implementation_Plan.md` (6 phases)
- **AppView backlog:** `decodingus/documents/planning/ibd-matching-appview-backlog.md` (5 items: IBD-AV-1 through IBD-AV-5)
- **Library decisions:** JDK 17 built-in crypto (X25519, AES-GCM, Ed25519), STTP WebSocket (already in project), custom IBS/IBD detector using HTSJDK (already in project) — no new dependencies
- **Existing scaffolding:** IBD tab placeholder in SubjectDetailView (lines 744-877), 19 i18n keys, `supportsAutosomalIbd` in TestTypeDefinition

### 2.2 Phase 4 UI — Ancestry & IBD Views
- **Source:** UI_Redesign
- **Scope:**
  - [ ] Ancestry composition visualization (pie/bar charts, admixture plot)
  - [ ] IBD matches table
  - [ ] Chromosome browser for shared segments
  - [ ] Relationship estimation display

### 2.3 OAuth Client Auth (replaces app passwords + REST/Kafka relay)
- **Source:** Edge_Client, `documents/atmosphere/11-Auth-and-Permissions.md`
- **Scope:** Replace `AuthenticationService.loginAtProto` (app-password `createSession`)
  with the OAuth authorization-code flow; add DPoP-bound tokens + refresh rotation
  (today `refreshJwt` is dropped); publish the `com.decodingus.atmosphere.navigatorCore`
  permission set; point `AsyncSyncService.pushCreate/Update/Delete` (stubbed) at direct
  PDS writes under the OAuth session.
- **Removes:** Kafka relay plan (cut); REST relay demoted to legacy bootstrap path
- **Note:** OAuth covers *writes* only — the AppView firehose ingest path is unaffected
  (reads/subscriptions are out of the permission spec)

### 2.4 Firehose Ingest (AppView, read path)
- **Source:** Edge_Client
- **Scope:** AppView subscribes to a **standard relay / Jetstream** for
  `com.decodingus.atmosphere.*` collections — no custom relay infrastructure
- **Depends on:** Granular sync completion (write side)

### 2.5 AppView Backfeed
- **Source:** Edge_Client
- **Scope:** Generate AppView records for 7 record types (public profile data)
- **Depends on:** Granular sync completion
- **Rework:** Under OAuth, writing into a user's repo needs an explicit granted scope.
  Prefer AppView-owned records the client reads; reserve user-repo writes for cases that
  require them (see `11-Auth-and-Permissions.md` §5)

---

## Tier 3 — Future / Exploratory

### 3.1 Internationalization — Additional Locales
- **Source:** UI_i18n_Guidelines
- **Status:** Architecture defined, only English (`messages.properties`) implemented
- **Scope:** German, Spanish, French translations; RTL layout support

### 3.2 Imputation
- **Source:** Edge_Client
- **Status:** No code
- **Scope:** Genotype imputation from chip/low-coverage data

### 3.3 Instrument Observation Records
- **Source:** Edge_Client
- **Status:** Lexicon defined, no implementation

### 3.4 Ancestral STR Reconstruction
- **Source:** Edge_Client
- **Status:** Lexicon defined, no implementation

### 3.5 Interactive Haplogroup Tree Visualization
- **Source:** UI_Redesign
- **Scope:** Zoomable/searchable tree with subject placement highlighted

### 3.6 Cross-Subject IBD Network View
- **Source:** UI_Redesign
- **Scope:** Graph visualization of IBD relationships across all subjects

---

## UI Polish (can be interleaved with other work)

- [x] SNP export from variant tables
- [x] Heteroplasmy indicators in mtDNA view
- [x] FASTA export for mtDNA sequences
- [x] Comparison report export (PDF/CSV)
- [x] Dashboard quick actions (recent files, pinned subjects)
- [x] Batch import workflow

---

## Code Quality (from codebase analysis, 2026-03-06)

- [x] ~~Consolidate duplicate `HaplogroupResult` types~~ — Analysis-layer renamed to `ScoredHaplogroup`, import aliases eliminated
- [x] ~~Resolve metrics model fragmentation~~ — `ContigSummary` eliminated in favor of `ContigMetrics`; manual conversion code removed
- [x] ~~Decide on SV pipeline integration~~ — SV pipeline is feature-toggled and integrated into `AnalysisCoordinator`; not orphaned
- [x] ~~Migrate `HaplogroupReconciliationRepository` to use `SyncableRepositoryBase`~~ — Extends base class, ~40 lines of duplicate sync logic removed
- [x] ~~Standardize error handling~~ — All GATK processors now return `Either[String, T]`; redundant `.left.map(_.getMessage)` conversions removed
- [x] ~~Standardize progress callback pattern~~ — All processor public APIs use `(String, Double, Double)` standard signature
