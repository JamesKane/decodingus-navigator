# DUNavigator Backlog

Last reviewed: 2026-08-05 (§3.3 block tree built; §1.1 archaic status corrected; four
pre-existing bugs found during that work recorded under Cross-cutting). Rewritten against the Rust
tree 2026-07-26.

> **This file was rewritten.** The previous version was the **Scala-era** inventory (last reviewed
> 2026-03-07) and described classes that no longer exist — `AncestryEstimator`, `SyncService`,
> `messages.properties`, V0xx Flyway migrations, HTSJDK/JDK-17 library decisions. The ScalaFX app was
> removed at cutover (commit `0dee32c`, 2026-06-19); that inventory is in git history if needed.
>
> **How to use this file.** It is the *cross-cutting* list — one entry per open thread, pointing at
> the design doc that carries the detail. The design docs in `documents/design/` are authoritative
> for scope and status; agent memory (`~/.claude/projects/.../memory/`, indexed by `MEMORY.md`) is
> the most current running record. Orientation lives in
> [`documents/design/HANDOFF.md`](design/HANDOFF.md).

---

## Tier 1 — Active / next up

Code exists or the design is settled; these are the near-term threads.

### 1.1 Archaic ancestry (Neanderthal / Denisovan)
- **Design:** [`design/ArchaicAncestry_Design.md`](design/ArchaicAncestry_Design.md)
- **Status (corrected 2026-08-02):** **Tier A shipped** (`230353b`, `#34`) and reports a *count*,
  never a % Neanderthal. **Tier B is built but gated OFF** (`#35`, `#40`) — the diagnosis is that it
  measured the wrong observable, not that the HMM is broken; read `#41`/`#42` before reopening. The
  "design draft, no code" status below was already stale when this file was written.
- **Scope:** Phase 1 = compute our own marker panel (EVA archaic VCFs + Ensembl-75 ancestral alleles
  + 1kGP AFR outgroup) + Tier A `count_archaic_markers` + domain/store/UI card — the 23andMe
  equivalent, for chip *and* WGS, reusing the ancestry-panel machinery. Phase 2 = Tier B segment HMM
  for WGS (% genome, chromosome browser, Nea/Den attribution). Phase 3 = trait associations, export
  / PDS records.
- **Blocked on (three open questions, §9):** informative-site overlap with 23andMe/AncestryDNA v5
  chip content after lift to CHM13; confirmation of Tier B inputs (diploid-caller VCF as the private
  variant source, bundled AFR outgroup + callability mask); choice of percentile cohort.
- **Note:** no fabricated Denisovan estimate for Europeans — report "none reliably detected" below
  the confidence threshold.

### 1.2 Deep-ancestry plumbing
- **Design:** [`design/ancient-ancestry-rebuild.md`](design/ancient-ancestry-rebuild.md) §7.18
- **Status:** The feature itself **ships and is enabled** (`ANCIENT_ANCESTRY_ENABLED = true`); this
  is the leftover plumbing.
- **Scope:** a GUI trigger for panel genotyping (`genotype_panel_for_subject` is CLI-only today —
  `navigator genotype-panel`), and folding panel batch mode into the project-wide analyze /
  deep-analyze streaming flow with progress.

### 1.3 Chromosome painter v2 — tail
- **Memory:** `chromosome-painter-v2.md`
- **Status:** Shipped in v0.1.0-alpha.13 (M1–M6 + the dense haplotype asset, published on demand).
- **Scope:** **PBWT acceleration** — `ReferencePhaser` stands in an IBS-match-length heuristic for
  full PBWT, which needs to land before the real ~5008-haplotype asset; `TrioPhaser` (Mendelian
  phasing) is still unwritten; trio validation on PRJEB36890 never run; the asset
  replace/quarantine/restore I/O path has no automated coverage (hand-verification recipe in memory).

### 1.4 Reassembly caller — calibration
- **Design:** [`design/haploid-reassembly-caller.md`](design/haploid-reassembly-caller.md)
- **Status:** Implemented and **default-on** (`reassembly.rs`, driven from `caller.rs`). What remains
  is tuning, not construction.
- **Scope:** active-region window size + merge tolerance (POC used ±40 bp; multi-SNV windows want
  ~150–250 bp); genotype threshold τ and its Phred-scaling into a GQ the publish gate trusts;
  fragment dedup on disagreement (consensus vs drop); whether active-region reassembly should
  subsume the homopolymer indel realignment in `realign.rs`; mtDNA heteroplasmy (fractional
  genotyping vs haploid argmax).

### 1.5 External-caller precedence — deferred follow-ups
- **Design:** [`design/external-caller-precedence.md`](design/external-caller-precedence.md) §13–14
- **Status:** All five phases landed; these were explicitly deferred.
- **Scope:** autosomal provenance-gating (§5.3 — needs a `CallProvenance` in `reconcile_diploid` so
  CRAM-genotyped dosages are skipped for a subject with an external autosomal source); PLINK
  `.bed/.bim/.fam`; rsID-based join (position-based today); project-scan sidecar auto-discovery for
  both call-set types; the GUI "Compare callers" button.

### 1.6 IBD network matching — consent & discovery UX
- **Design:** [`design/ancestry-ibd-asset-wiring.md`](design/ancestry-ibd-asset-wiring.md),
  [`IBD_Matching_Implementation_Plan.md`](IBD_Matching_Implementation_Plan.md)
- **Status:** Detection, identity math, the encrypted exchange channel (X3DH/AES-GCM), signed
  attestations, and the pairwise consensus-IBD chromosome browser are **all built and live-validated**.
  The consent/discovery surface landed 2026-08-02 (branch `feat/ibd-matching-ux`): a durable
  `ibd_request` ledger (mig 0041) + `App::refresh_matching`, a top-level **Matching** tab
  (Suggestions / Requests / Results) replacing the per-subject discovery cards, an informed-consent
  modal, and the two previously unwired AppView endpoints (`/ibd/dismiss`, `/ibd/attest`).
  Attest needed a companion AppView change — `/ibd/suggestions` now returns the caller's own
  `target_sample_guid`, without which `owns_sample` could never be satisfied from the edge.
- **Scope remaining:** background polling + an unread badge for inbound consent requests (the
  Community 🔔 pattern); a Settings discoverability opt-in; a UI path for the *direct*
  `exchange_request(partner_did, …)` initiator (still test-only); the segment ideogram for persisted
  exchange results; and **live two-peer validation** of the whole flow against a running AppView.
- **Note:** both subject pickers in this area (Matching, and `consensus_ibd_section` in `ui/ibd.rs`)
  are now filter + virtualized `show_rows` rather than `ComboBox` — a workspace can hold 10k
  subjects, and a `ComboBox` builds a widget per entry per frame. Follow that pattern for any new
  roster-wide picker.

### 1.7 Packaging & release — open items
- **Design:** [`design/packaging-and-release.md`](design/packaging-and-release.md)
- **Status:** Shipping (all four installers build on a `v*` tag; assets fetched on demand from the
  GitHub asset release).
- **Scope:** code signing + notarization (Apple Developer ID $99/yr; a Windows cert) — deferred for
  alpha with a documented Gatekeeper work-around; the Linux glibc-2.28 container CI is authored but
  **has never been run**; `default_reference_sha` is still `None` for all four builds
  (`navigator-refgenome/src/registry.rs:172`), awaiting confirmed publisher checksums.

---

## Tier 2 — Designed, not started

Verified 2026-07-26 to have no implementation in the tree.

### 2.1 Realignment module — **in progress** (phase 1 landed 2026-08-08)
- **Design:** [`design/realignment-module.md`](design/realignment-module.md) — revised 2026-08-08
  after a phase 0 spike that **retracted the module's motivating premise** (ancestry is *not*
  build-locked; off-build samples already estimate ancestry through the multi-build IBD panel) and
  reversed the backend decision to pure-Rust `minimap2-pure-rs`. Read the correction blocks before
  planning further work — whether the remaining payoff justifies the module is an open product
  question.
- **Scope:** revert + realign GRCh37/38 vendor WGS to CHM13v2 / hs1; aligner-index cache in
  `navigator-refgenome`, job orchestration + provenance, opt-in background job with warnings.
- **Done:** phase 0 spikes; **phase 1** — `navigator-analysis/src/revert/` (stage A): primaries-only
  revert with orientation restore and `OQ` preference, a disk-backed external merge sort that
  collates by read name, and synchronized paired-FASTQ output. 19 tests including BAM/CRAM parity.
- **Next:** phase 2 — a `navigator-align` crate wrapping the mapper, with part-by-part index
  build/map (`-I` sized from RAM; see Decision 4 — a monolithic index costs ~19 GB and is the
  failure mode to avoid).
- **Do not confuse** with `navigator-analysis/src/realign.rs`, which is *indel local realignment*
  (plan §4b) and is a different thing entirely.

### 2.2 Distributed compute grid
- **Design:** [`design/distributed-compute-grid.md`](design/distributed-compute-grid.md)
- **Scope:** a Seti@Home-style layer — the AppView publishes public-ENA work units, Navigator
  instances reserve a lease, fetch, realign to CHM13, run the analysis stack, submit signed results,
  and earn capped compute credit. Cross-repo (Navigator worker + AppView coordinator + shared wire
  records). Depends on 2.1.

### 2.3 Academic / public-dataset (ENA) import
- **Design:** [`design/academic-ena-import.md`](design/academic-ena-import.md)
- **Scope:** the IRB- and publishing-facing import profile. Note the Scala-era **`EnaClient`**
  (ENA/1KG metadata resolution) was never re-ported to Rust, so metadata resolution is part of this
  work rather than a separate item.

### 2.4 Local-LLM M6 — project-level summary
- **Design:** [`design/local-llm-expansion.md`](design/local-llm-expansion.md) §M6
- **Status:** M0–M5 shipped (narration, ask-my-results chat, per-tab "Explain this"). No
  `project_summary` exists.
- **Scope:** a project-level AI summary — a different audience from the per-subject brief; do last.

### 2.5 Pangenome / GAM as a data source
- **Design:** [`design/pangenome-gam-data-sources.md`](design/pangenome-gam-data-sources.md),
  [`design/PangenomeExpansion.md`](design/PangenomeExpansion.md)
- **Scope:** turn graph alignments (GAM/GAF) into the records Navigator and AppView already model.
  Explicitly a post-launch horizon item — the AppView *storage* side is modelled, the producer side
  is not.

---

## Tier 3 — Future / exploratory

### 3.1 Imputation
- **Status:** No implementation.
- **Scope:** genotype imputation from chip / low-coverage data. Lexicon record `imputation` is
  defined but unused.

### 3.2 Ancestral STR reconstruction
- **Status:** No implementation. Lexicon defined only.

### 3.3 Interactive haplogroup tree visualization
- **Status:** **Mostly built.** Per-subject: `ui/descent.rs` draws the root→terminal path
  (YFull-YReport style, Simple and Advanced densities) and `ui/branch.rs` gives a per-marker branch
  report with TSV export. Cohort: the project **block tree** below.
- **Built, cohort-scoped** — [`design/project-block-tree.md`](design/project-block-tree.md), branch
  `feat/project-block-tree` (13 commits, unpushed). A project Y **block tree**: induced subtree over
  the members' terminals, equivalent-SNP blocks, and **candidate branches** inferred from private
  variants two or more members share — the thing a published tree cannot show. A `ProjectTab::Tree`
  canvas draws it; clicking a candidate opens the per-carrier read evidence behind it.
- **Scope remaining:** a zoomable/searchable *whole-tree* view (this is cohort-scoped by design), and
  the open items in that doc's §11 — chiefly **no GUI trigger for the private-Y batch** (CLI only,
  and candidates cannot fire without it) and 162 `:`-suffixed terminals that fall to `unplaced`.
- **Note:** phase 3 grew well past its original scope because private-Y existed for exactly one
  subject workspace-wide. That pulled in a `private-y --project` batch, a **VCF-backed private-Y
  engine** for the majority of members who have no alignment, and an artefact-filter stack. Four
  pre-existing bugs surfaced on the way — see the design doc and the entries below.

### 3.4 Cross-subject IBD network view
- **Status:** **Partial.** Cross-subject Y ranking (`ymatch`), federated
  `network_suggestions_section`, per-pair consensus IBD, and a relatives section all exist.
- **Scope remaining:** the graph visualization of IBD relationships across all workspace subjects.

### 3.5 Additional locales
- **Status:** en + es at key parity (`crates/navigator-domain/locales/{en,es}.txt`, parity-tested) —
  the old "English only" entry was stale.
- **Scope:** German, French; RTL layout support.

---

## Cross-cutting / smaller

- **i18n tail** — `self.status` transient strings and `format!` dynamics are still English; the
  key-based UI is at en/es parity.
- **IBD sliding-window rolling update** — `find_candidate_segments`
  (`navigator-analysis/src/ibd.rs:326`) recomputes the whole window at every position, i.e. O(n·w).
  A rolling update makes it O(n) and, per
  [`design/simd-optimization-targets.md`](design/simd-optimization-targets.md), dwarfs any
  vectorization gain. The cheap `target-cpu` win from that doc is already done.
- **Unified-walker perf plateau** — stuck near 5×: a serial unmapped-tail sweep plus the
  single-largest-contig floor set wall time. Split big contigs / parallelize the sweep to push further.
- **STR loose ends** — the CHM13 lift dropped 33 named chrY markers (incl. DYS19/391/426; the table
  retains their GRCh38 values for the BAM path); multi-copy marker aggregation; reference auto-download.
- **Compare multi-select** — the Compare view still takes a limited selection.
- **Tree cache staleness** — `fetch_tree` is cache-first with no freshness check, so a stale cached
  haplotree silently under-places. Open.
- **Fixed while building the block tree** (branch `feat/project-block-tree`), all pre-existing and
  costing continuously:
  - a **chrY-only VCF was routed to the autosomal 1240K path** (`looks_like_genotyped_callset_vcf`
    asked only "does it emit `0/0` rows?" and never looked at contigs), so a Big Y import produced no
    Y variant set, no placement and no private-Y source;
  - the **VCF import discarded QUAL/FILTER/DP/GQ/AD**, leaving nothing downstream able to tell a 40×
    hom-alt call from a 2-read artefact (migration 0042);
  - **GVCF sidecar discovery looked in one directory for one filename spelling**, missing
    `gatk4/chrY.g.vcf.gz` entirely — every affected subject decoded a whole chromosome a sidecar
    could have answered (50–90 s → 1.3–3.3 s);
  - **`localize` leaked its alignment copies** — cleanup lived in one caller while three created
    them; reached 687 files / **145 GB** and filled the volume mid-run. Now RAII with a refcount.
- **`project_report` reads every artifact payload** — it still uses the unfiltered
  `artifact::list_for_alignments`, which pulled gigabytes of `tree-genotype` JSON to read a few small
  kinds. `artifact::list_for_alignments_of_kind` exists; convert the caller.
- **mtDNA FASTA export** — the Scala app had it; the Rust `export.rs` covers coverage / read-metrics
  / ancestry / mtDNA-variants / IBD-segments / branch / descent / callable-BED / subject-brief
  (TSV + HTML), but **not** FASTA. Carry it over if still wanted. No PDF export exists either.
- **Live-PDS validation** — §4 attestation publish and §5 idempotent putRecord + PULL are confirmed
  at the protocol level by curl, but the over-the-wire run of the App code is environment-blocked
  (the agent host can't reach the Apple-container subnet). Not code-blocked.

---

## Known gaps left by deliberate removals

- **Autosomal VCF cross-source concordance.** Variant-level reconciliation (`reconcile_variants` /
  `VariantStatus` / `ReconciledVariant` + the "Cross-source concordance" card) was removed at merge
  `c471ede` because it predated the DNA-type-agnostic consensus engine. Y (name-keyed, tree-aware),
  mt (vs rCRS) and autosomal (IBD-panel diploid genotyping) consensus all have cards — but
  position-level concordance of ≥2 imported **autosomal VCF** variant sets at their own positions has
  no replacement. The intended fix is an autosomal consensus *adapter* feeding the consensus engine,
  **not** reviving the old path. See [`design/MultiRunReconciliation.md`](design/MultiRunReconciliation.md).

---

## Social layer — deferred slices

Tiers 3a (peer DMs), 3b (publish `feed.post`) and 3c (recruitment, respond-only) are **all built and
merged**. Deferred by design:

- Navigator-native campaign **creation** — blocked on the groupProject-PDS project bridge (campaign
  creation is gated to an AppView group_project admin, which Navigator cannot act as yet).
- Feed voting / report / block actions.
- Threaded federated replies — `FeedPostRecord.reply` is modelled, but the Feed UI only posts
  top-level; needs parent/root at-uri tracking.
- DM follow-ups — truly-async handshake (persisted ephemeral → derive-on-arrival, so "Connect"
  needn't be simultaneous), background auto-poll of conversations, typing indicators / read receipts,
  in-DM block and delete.

---

## AppView-side (tracked in the `decodingus` repo, listed for completeness)

- **`fitDistance` ingest** — needs a shared-crate rev bump.
- **IBD attestation indexing** — the Navigator side publishes; the AppView must index.
- **Firehose ingest** — subscribe a standard relay / Jetstream for
  `com.decodingus.atmosphere.*`; no custom relay infrastructure (the Kafka relay plan was cut).
- **AppView backfeed** — prefer AppView-owned records the client reads; reserve user-repo writes for
  cases that require them, since under OAuth that needs an explicit granted scope
  (see [`atmosphere/11-Auth-and-Permissions.md`](atmosphere/11-Auth-and-Permissions.md) §5).
- **Genome Regions API** — the client is ready and falls back to a bundled region set; the server
  side is an AppView deployment task
  ([`GenomeRegionsAPI_Specification.md`](GenomeRegionsAPI_Specification.md)).

---

## Closed since the Scala-era backlog

Recorded so these are not re-scoped as gaps. Per-feature detail is in agent memory and the design docs.

| Old item | Outcome |
|---|---|
| 1.1 Ancestry estimation | Shipped — super-pop + fine admixture, PCA, painting, and qpAdm deep ancestry; panels rebuilt at 200k depth |
| 1.2 Chip → haplogroup | Shipped |
| 1.3 Multi-run reconciliation Phase 3 | Shipped — including the **manual override UI** (`ui/sources.rs`: override controls, reason, audit log) that was the last open checkbox |
| 1.4 Genome Regions API (client) | Shipped; server deployment is AppView-side (above) |
| 1.5 Granular record sync | Shipped — `sync_outbox`, idempotent putRecord at TID, PULL reconcile |
| 1.6 Atmosphere lexicon alignment | Shipped (phases A–D) |
| 2.1 IBD matching system | Built end-to-end; only the consent/discovery UX remains (1.6) |
| 2.2 Phase 4 UI — ancestry & IBD views | Shipped — pie/bar breakdowns, PCA scatter, match browser, chromosome browser |
| 2.3 OAuth client auth | Shipped — PKCE + DPoP-bound tokens via `du-atproto` / `navigator-sync`; app passwords gone |
| 3.3 Instrument observation records | Shipped as sequencing-lab + instrument/platform inference with AppView lookup |
| UI polish (SNP export, heteroplasmy indicators, comparison export, dashboard actions, batch import) | Shipped in Rust, **except** mtDNA FASTA export and PDF output (see Cross-cutting) |
| Code quality (2026-03-06 list) | Obsolete — those were Scala types. The Rust equivalents are the rustfmt adoption, the God-file split, and the dedup/simplification round |
