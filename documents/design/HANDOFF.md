# Rust rewrite — handoff / resume notes

Last updated: 2026-07-26 (status re-check; the build log below is still as of 2026-06-20). The
rewrite is now **trunk on `main`** (the legacy ScalaFX app was removed at cutover and lives in git
history only). Pick up here next session.

> **Landed since the 06-20 body was written** (not reflected in the sections below, which are kept
> as a build log): ROH / F_ROH endogamy detection · external-caller precedence + autosomal external
> ingest (EIGENSTRAT 1240K + GATK4 gVCF) · **deep/ancient ancestry shipped and ENABLED** (qpAdm) ·
> fine/PCA panels rebuilt at 200k depth with continental-European populations · chromosome painter
> v2 (parent-split phasing + haplotype-copying LAI) · progressive autosomal consensus · window/session
> state persistence · packaging through GitHub-release on-demand assets, shipping as
> `v0.1.0-alpha.13`.

> The detailed running record is **agent memory** (`~/.claude/projects/.../memory/`, indexed by
> `MEMORY.md`) — it is more current and granular than this file. This doc is the orientation +
> active-work pointer; the per-topic memory files carry the specifics. The pre-06-07 theme list that
> used to live here is superseded; see memory for ancestry, UQMW, tree-provider, etc.

## Three-repo topology

- **DUNavigator** (this repo) — the desktop edge app (egui). Pulls the shared crates from the
  sibling repo by **path**, so working-tree changes build immediately (no rev bump needed).
- **decodingus-shared** (`/Development/decodingus-shared`) — `du-domain` / `du-atproto` / `du-bio`.
  Federated record contracts live in `du-domain::fed`; shared SHA-256 helpers in `du-bio::hash`
  (added + merged to its `main` 2026-06-20). NB its working tree often carries WIP on a feature
  branch — stage only your own files when committing here.
- **decodingus** (`/Development/decodingus`, the AppView) — the PostgreSQL hub + web (`rust/`).
  Pulls the shared crates by pinned **git rev**; **bump the rev (or add a local `[patch]`) to pick
  up newer shared changes**. The social layer (all three roadmap tiers) is built here on branch
  `feat/social-layer-orchestration` — signed `/api/v1/social/*` Edge endpoints, `du_db::social` +
  `du_db::notification`, web inbox/feed.

## State

Workspace builds clean (`cargo clippy --all-targets -- -D warnings`); `cargo test --workspace`
green except a **known parallel-isolation flake** (`y_profile_build_persists_and_reloads`,
`import_23andme_*`) that passes in isolation. `cargo fmt` clean is a per-commit gate.

### Social layer (Navigator/Edge side) — DONE, on `main` (not pushed)

All social Edge tiers shipped and **merged to `main`** (the `feat/social-community-tab` branch was
fast-forwarded into `main` and deleted; nothing pushed yet). Roadmap:
`decodingus/documents/planning/social-layer-roadmap.md`. Per-feature memory:
[[social-feedpost-publish]] (3b), [[social-peer-dms]] (3a), [[social-recruitment-3c]] (3c).
Commits: `96b8577` signed Edge client · `5a0cb5e` Community tab · `06a754c` 3b feed.post ·
`a142081` 3a peer DMs · `3e943e7` 3c recruitment. (AppView `decodingus` `main`: `c4cc15c` recruitment
Edge — plus the user's own follow-on extending it to full create+respond.)

A later **maintenance feature** also landed on `main`: **Clear subject data** — a "Clear data" button
(subject header, next to Delete) + confirm modal that resets a subject's analysis (runs/alignments/
artifacts, Y/mt haplogroups + consensus + reconciliation audit/overrides, ancestry, IBD results,
chip/STR/variant/mtDNA profiles) while keeping the subject (name/sex/IDs, project memberships, MDKA).
`biosample::clear_data` (one transaction). Also closed the orphan root-cause: `purge_alignment_derived`
now clears the reconciliation audit + per-alignment ancestry, and `delete_biosample` sweeps via
`clear_data` so a delete can't leave dangling rows (fixed subject 103589's stale Y consensus).

Below is the original communication-core build log (kept for reference):

- `96b8577` — **signed Edge client**: `navigator-sync::social::messages` (canonical signing strings
  mirroring `du_db::social::messages` byte-for-byte) + `navigator-app::social` (device-key-signed
  POST/GET helpers like the IBD `exchange` client, + response DTOs). Methods: `support_threads` /
  `support_thread` / `open_support_thread` / `reply_support_thread` / `community_feed` /
  `post_community` / `notifications` / `mark_notification_read`. Unit-tested (canonical strings +
  DTO wire-shape round-trip).
- `5a0cb5e` — **Community tab UI**: top-level `Nav::Community` (Support / Feed / Notifications
  sub-tabs) + app-bar unread **🔔** bell; `ui/community.rs`; worker Commands/Events; en/es i18n.
  Sign-in gated.

- **3b — publish `feed.post`** (DONE, uncommitted on this branch): Navigator now publishes
  `com.decodingus.atmosphere.feed.post` to the signed-in PDS via the durable sync outbox, completing
  the federated feed loop the AppView already ingests. `FeedPostRecord` lives in shared
  `du-domain::fed` (top-level `createdAt`, optional `topic` + `reply.{root,parent}.uri`; PII-free,
  no `WireF64`); `App::publish_feed_post` enqueues under `NS_FEED_POST` with a fresh per-post
  `entity_ref` (append-only — never coalesced, deliberately **not** in `PUBLISHED_COLLECTIONS` so a
  PULL can't resurrect a deleted post; errors for a `did:key` identity). Wired as an **opt-in
  checkbox** ("Publish publicly to my PDS") on the Feed composer, gated to PDS accounts; the native
  signed-Edge post still happens, the federated copy mirrors back badged "via Atmosphere". en/es
  i18n + unit tests (shared wire-shape round-trip, app builder). Not yet live-tested against a PDS.

- **3a — peer DMs** (DONE, uncommitted on this branch): end-to-end-encrypted peer DMs over the D1
  relay, reusing the generic exchange crypto/session driver unchanged. New `navigator-store` mig
  `0031_dm` (`dm_conversation` persists the derived **session key** + outgoing seq counter so a
  conversation is async + restart-safe; `dm_message` with `UNIQUE(session_id,from_did,seq)` dedupe)
  + `store/dm.rs`. New `navigator-app/src/dm.rs`: `dm_initiate`/`dm_incoming`/`dm_ready`/
  `dm_consent`/`dm_connect` (one-time handshake → persist key)/`dm_conversations`/`dm_messages`/
  `dm_send`/`dm_sync`, purpose `GENEALOGY_PII` (AppView already titles it — no AppView change). UI:
  **Messages** sub-tab in Community (start-by-DID + incoming inbox + ready-to-connect + conversation
  list + transcript/composer) **and** a "Message" button on IBD exchange results that opens a DM +
  jumps to Messages. en/es i18n + store/app unit tests. **Crux:** the session key is persisted so
  only the *first* connect needs both peers online; messaging is then fully async. Not yet live-tested
  (needs AppView + 2 peers, like the IBD exchange).

- **3c — recruitment invitations** (DONE, respond-only; **cross-repo**, uncommitted on both branches):
  testers can view + accept/decline recruitment invitations from Navigator. The AppView recruitment
  engine + web UI existed but had **no signed-Edge API** — so this ADDED endpoints on the AppView
  (`decodingus`, on `main`): `du_db::recruitment::messages` (poll/respond) + `routes/recruitment_edge.rs`
  (`GET /api/v1/recruitment/invitations`, `POST /api/v1/recruitment/respond` — notifies the researcher
  on accept) + router registration; DB-gated test **passes against the real DB**. Navigator side:
  `navigator-sync::recruitment::messages` mirror, `navigator-app/src/recruitment.rs`
  (`recruitment_invitations`/`recruitment_respond`, reusing `social_get`/`social_post` now `pub(crate)`),
  worker Commands/Events, and a "Recruitment invitations" section atop the Community → **Notifications**
  sub-tab (Accept/Decline). en/es i18n. **Scope = respond-only** (user-confirmed): campaign *creation*
  is gated to an AppView group_project admin, which Navigator can't act as until the groupProject-PDS
  project bridge exists. NB: `decodingus` `main` has pre-existing clippy debt in `du-db/src/ystr.rs`
  (a `0..=10` index loop) unrelated to 3c — my files are clippy-clean.

**Deferred** (later slices, all KEPT): Navigator-native campaign **creation** (needs a "list my
recruitable group_projects" Edge endpoint + the groupProject-PDS project bridge); feed voting/report/block
actions; threaded federated replies (`FeedPostRecord.reply` is modelled but the Feed UI only posts
top-level — needs parent/root at-uri tracking); **DM follow-ups** — truly-async handshake (persisted
ephemeral → derive-on-arrival, so "Connect" needn't be simultaneous), background auto-poll of
conversations (MVP syncs on open + manual refresh; request/consent already notify via the bell),
typing indicators / read receipts, in-DM block/delete.

### Other work landed this session (on `main`)

- **FTDNA project import** (PR #6, merged) — roster/ancestry/Y-STR CSV import, match/dedup,
  review→commit, Y-STR autoclustering. See `memory/ftdna-import-platform.md`.
- **Project report membership fix** (`fb0f186`) — `project_report`/members/count now read the M:N
  `biosample_project` table ∪ legacy home column, so an FTDNA-merged subject shows in the report.
- **sha256 dedup** (`383d6d5` + shared) — consolidated scattered SHA-256 impls onto `du_bio::hash`.
- **Run-delete derived purge** (`0c252cd`) + **source_file FK unlink** (`9f974bf`).

Untracked: `CLAUDE.md`, `GEMINI.md` (leave). A stray `crates/.claude/settings.local.json` is
recreated by the environment — handled by `exclude = ["crates/.claude"]` in the root `Cargo.toml`.

## Build / validate

```bash
cargo build && cargo test --workspace
cargo run -p navigator-ui            # desktop app
```

### Live (`#[ignore]`) tests — real data, run explicitly
Test sample: `/Users/jkane/Genomics/GFX0457637/GFX0457637.pbmm2.chm13v2.bam` (CHM13 HiFi, ~4×,
male, Y=R-FGC29071, mtDNA=U5a1b1g, European). Reference: `/Users/jkane/Genomics/chm13v2.0/chm13v2.0.fa`.

```bash
GFX_CHM13_BAM=/Users/jkane/Genomics/GFX0457637/GFX0457637.pbmm2.chm13v2.bam \
GFX_CHM13_REF=/Users/jkane/Genomics/chm13v2.0/chm13v2.0.fa \
NAVIGATOR_ANCESTRY_PANEL=/Users/jkane/.decodingus/ancestry/ancestry_panel_chm13v2.0.bin \
NAVIGATOR_ANCESTRY_PCA=/Users/jkane/.decodingus/ancestry/ancestry_pca_chm13v2.0.bin \
  cargo test -p navigator-app --release \
  validate_gfx_chm13_ancestry local_ancestry_paints_gfx gfx_sex_is_male -- --ignored --nocapture
```
Expected: European ~98% (admixture), DNA painting EUR-dominant, sex=Male. Other ignored live tests:
`validate_gfx_chm13_haplogroups` (Y/mt), parity_real.rs (HG002 env), PDS publish (PDS_TEST_URL).
`NAVIGATOR_ANCESTRY_PCA_ANCIENT` points the PCA-GMM at an ancient-component asset when present.

## Ancestry assets (regenerable; not committed)

Installed at `~/.decodingus/ancestry/`:
- `ancestry_panel_chm13v2.0.bin` — AF panel (genotyping + admixture; the default the app loads)
- `ancestry_pca_chm13v2.0.bin` — PCA loadings + per-pop centroids (drives PCA-GMM + nMonte)

Today's assets come from the archived genotype matrix `~/Genomics/archive/1kgp_chm13_pca_build/`
(`gt_all.tsv.gz` 1000G + `sgdp_gt.tsv.gz` + `combined_pops.txt`):
```bash
A=~/Genomics/archive/1kgp_chm13_pca_build; O=~/.decodingus/ancestry
navigator-panelbuild fine-panel --matrix $A/gt_all.tsv.gz,$A/sgdp_gt.tsv.gz \
  --samples $A/samples.txt,$A/sgdp_subset_samples.txt --pops $A/combined_pops.txt --out $O/ancestry_panel_chm13v2.0.bin
navigator-panelbuild pca        --matrix $A/gt_all.tsv.gz,$A/sgdp_gt.tsv.gz \
  --samples $A/samples.txt,$A/sgdp_subset_samples.txt --pops $A/combined_pops.txt [--basis-pops modern.txt] --out $O/ancestry_pca_chm13v2.0.bin
```
The **next-gen** asset path is the global pipeline in `scripts/ancestry-panel/` (modern + ancient
deep components over a 1240k-restricted panel, projection-mode PCA, CDN publish) — needs the
datasets fetched (verify `# VERIFY` URLs; slice panel sites to avoid the multi-TB pull).

## Architecture / design pointers

- `documents/design/SubjectCentricModel.md` — donor-centric tab model (P1–P3 implemented).
- `documents/design/AncestryAnalysis.md` — the 3 estimators + ancient-asset build + nMonte/G25.
- `documents/atmosphere/` — the lexicon spec; `du-domain::fed` is the implemented write subset.
- `decodingus/documents/planning/social-layer-roadmap.md` — the social-layer build plan (AppView
  side built; Navigator/Edge side is the active work, communication core done).
- `documents/design/` — the design backlog (FTDNA import, BISDNA, realignment, packaging, SIMD,
  pangenome-GAM, scala-rust-gap-analysis, …).
- `documents/BACKLOG.md` — **Scala-era** feature inventory (March 2026, pre-rewrite); use as the
  master feature list, not current status.
- Agent memory (`~/.claude/projects/.../memory/`) is the most current running record.

## Remaining gaps

The 06-07 audit is superseded — most of it shipped over 06-10 → 06-13 (UQMW + parallel walker,
DecodingUs Y-tree provider, BISDNA + chip-haplogroup import, vendor/mtDNA import, diploid SNV+indel
caller, settings UI, Y-STR reporting, report exports, genome-region/ideogram, federated IBD phases
1–2 + the encrypted exchange channel, sync durability, FTDNA project import). Per-feature status
lives in agent memory (`MEMORY.md` index) — treat that as authoritative, not this file.

Still open, broadly (corrected 2026-07-26 against the tree):
- **Social layer (Edge)** — 3a/3b/3c are all **DONE and merged to `main`** (see the section above;
  the earlier "still to build" line here was stale). What remains are the *deferred slices*:
  Navigator-native campaign **creation** (blocked on the groupProject-PDS bridge), feed
  voting/report/block, threaded federated replies, and the DM follow-ups (truly-async handshake,
  background auto-poll, typing/read receipts, in-DM block).
- **IBD network matching** — detection, identity math, and the encrypted exchange channel are built;
  the consent/discovery/chromosome-browser UX is the remaining surface.
- **Design backlog** in `documents/design/` — still genuinely design-only: **realignment module**
  (revert + realign to CHM13 — note `analysis/realign.rs` is *indel* local realignment, a different
  thing), **distributed-compute-grid**, **academic-ENA import**, **pangenome-GAM**, and local-LLM
  **M6** (project-level summary). Packaging and SIMD have both since been partly acted on — see
  their headers. **Archaic ancestry** (`documents/design/ArchaicAncestry_Design.md`) is the newest
  design-only item and the shortest path to a visible feature.
- **Smaller** — i18n `self.status`/`format!` tails; Compare multi-select; the unified-walker perf
  plateau (~5×: serial unmapped-tail sweep + the single largest contig floor wall time — split big
  contigs / parallelize the sweep to push further); the IBD sliding-window rolling update
  (`ibd.rs:326`, still O(n·w)). *Ancestry genotype-level pooling is now **done*** — `reconcile_diploid`
  + the shared `DiploidProfile` consensus.
- **AppView side** — `fitDistance` ingest (needs a shared rev bump); IBD-matching backlog; IBD
  attestation indexing.
- **Live-PDS validation** — §4/§5 are protocol-confirmed by curl only; the over-the-wire run of the
  App code is environment-blocked, not code-blocked.

## Recommended next steps (pick one)

*(Rewritten 2026-07-26 — the previous list was completed.)*

1. **Archaic ancestry (Neanderthal/Denisovan)** — `documents/design/ArchaicAncestry_Design.md`.
   The only design-only item that is a *new user-visible feature*, and Phase 1 reuses the ancestry
   panel machinery almost verbatim. Three open questions gate it (chip overlap after lift to CHM13,
   Tier B inputs, percentile cohort).
2. **Finish the deep-ancestry plumbing** — a GUI trigger for panel genotyping (CLI-only today) and
   folding panel batch mode into the project-wide analyze flow with progress (§7.18).
3. **IBD network matching** — the consent/discovery/chromosome-browser UX is the remaining
   user-facing surface.
4. **Reassembly-caller calibration** — the caller ships default-on; its Open Questions (window size,
   τ/GQ, fragment dedup) are unresolved tuning on the full truth set.
5. **i18n tail** — `self.status` transient strings + `format!` dynamics are still English (the
   key-based UI is at en/es parity).

For the broader unported-from-Scala inventory and per-feature status, the authoritative source is
agent memory (`MEMORY.md` index) — not this file.
