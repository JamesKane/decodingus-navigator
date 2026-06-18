# Scala → Rust gap analysis (functionality still missing)

A systematic pass over the legacy Scala/ScalaFX source (`src/main/scala/com/decodingus/`,
238 files / 66k LOC, still on this branch) against the Rust rewrite (`crates/`). This catalogs
**capabilities present in Scala that are absent or incomplete in Rust**, grouped by subsystem and
prioritized.

> **Revised 2026-06-17.** Much of the original (2026-06-12) backlog has since landed. This revision
> reconciles every section against the current tree + project memory, marks completed work with its
> commit, and re-derives the priority summary from what genuinely remains. Net: §1 (STR), §2, §3,
> §5(P1), §8a fully done; §4, §6, §7, §8 substantially advanced. The freshest remaining slices are
> **the live federated arc (§4 IBD exchange / segment attestation)** and **sync conflict/PULL (§5
> P2)** — both to be scoped against the running AppView — plus small §8-misc UI items.

## Method & exclusions

Six parallel sweeps (analysis, haplogroup/Y/STR, ancestry/IBD, workspace/sync, UI, refgenome),
then reality-checked against the tree. **Deliberately excluded** (not gaps — by design):

- **External-tool wrappers.** The Rust engine is GATK/samtools/bcftools-free by design
  (`CLAUDE.md`). So `GatkRunner`, `GatkToolProcessor`, GATK `HaplotypeCaller`/`Mutect2`/`LiftoverVcf`
  orchestration, `samtools faidx`/`flagstat` shelling, etc. are **not** counted as missing — the
  Rust caller/walkers/noodles replace them. Where a wrapper produced a real *output capability* the
  Rust side lacked (diploid/indel calling, VCF liftover), that capability is tracked below.
- **Already ported.** coverage+callable (+per-contig histograms, MAD + per-base exclusion
  fractions), read metrics, sex, SV, the unified quality-metrics walker (parallel, CRAM zero-copy,
  smooth bp-progress), the haploid caller, **diploid SNV+indel calling + whole-genome diploid VCF**,
  **subject-level consensus joint diploid genotype** (force-call merge), heteroplasmy, header probe,
  lab/instrument/library-stats inference, pipeline-sidecar fast-path import, **flagstat/Picard
  metrics importers**, Y+mt haplogroup placement (parsimony guard, FTDNA+DecodingUs providers) +
  **genome-level consensus placement**, chip-raw-data haplogroups (23andMe/AncestryDNA), BISDNA
  Y-SNP import, **vendor VCF/mtFull-FASTA → variant import**, **chip → autosomal ancestry**, ancestry
  (AF-likelihood / PCA-GMM / nMonte / ADMIXTURE / fine-admixture / local-ancestry painting — now
  **consensus-driven** + diploid pair-state painting; federated **from the consensus**), refgenome
  retrieval+cache, Y/mt rotation-aware liftover, masked-rCRS build, the **local** IBD math + chip
  panel + genetic map + manifest verification, **genome-region service + ideogram tab**, **Y-region
  quality modifiers**, **settings UI**, **report/TSV/HTML/BED exports**, **persistent sync outbox +
  history**, **federated IBD device-key signing + encrypted exchange channel (X25519/AES-GCM)**.

Status legend: **MISSING** = no Rust equivalent · **PARTIAL** = some behavior, notable holes · **DONE**.

---

## 1. STR calling & reference

### §1a STR reporting + panel classification — **DONE** (committed e5977cd)
`strpanel.rs` ports `str-panels.conf` (FTDNA Y-12..Y-111 + YSEQ tiers, classify/badges, multi-copy
order-independent match); `compare_profiles` conflict detection; `ystr_report_section` Y-DNA tab
(By-Panel transposed grid + filterable All-Markers + Consensus). Validated on real FTDNA(→Y-111)/
YSEQ(→Alpha) exports. *Remaining:* Big-Y Y-500/700 enumeration; Y-STR genetic distance / TMRCA /
match ranking (→ §2 cross-subject).

### §1b STR calling from BAM/CRAM — **PARTIAL** (caller foundation DONE; was MISSING)

Native STR calling now exists via the **enclosing-read** model (HipSTR/GangSTR), committed 986e00b:

| Capability | Scala | Status |
|---|---|---|
| STR repeat-count **calling from BAM/CRAM** (CIGAR-aware span, stutter model, HIGH/MED/LOW tiers) | `analysis/StrCaller.scala` | **DONE** — `strcaller.rs` enclosing-read genotyper (per-read count off the CIGAR, geometric stutter ML, haploid chrY / diploid elsewhere). Validated on GRCh38 chrY (4581 genotypes, ref loci measure ref_copies exactly) |
| STR reference parse — HipSTR BED → tract loci (period, ref_copies, motif) | `refgenome/StrAnnotator.scala` | **DONE** — `strref.rs` (end-inclusive tracts; per-contig + min-period filter) |
| STR reference **gateway** — download + per-build liftover + cache | `refgenome/StrReferenceGateway.scala`, `StrReferenceCache.scala` | **PARTIAL** — cross-build **liftover DONE** (b9a7eed: `ReferenceGateway::lift_hipstr_bed`; GRCh38→CHM13 chrY validated on a real CHM13 CRAM, offsets build-independent); reads via `NAVIGATOR_STR_REFERENCE` / `~/.decodingus/str/{build}.hipstr_reference.bed.gz`. *Remaining:* auto-download of the GRCh38 HipSTR BED |
| STR marker comparison (Simple/MultiCopy/Complex value matching) | `str/StrMarkerComparator.scala` | PARTIAL (panel compare done in §1a) |

**Vendor bridge (§1b-vendor) — BUILT + VALIDATED + WIRED.** The DYS→coordinate mapping was **free**
(the HipSTR BED names ~206 chrY DYS markers; `strref` parses them). The FTDNA *convention* is
calibrated in `strmarker.rs` (5fc6641) against a **14-kit Big Y corpus**: 70 reliable (offset 0) + 22
convention-offset (±1–3, constant across kits) + 29 excluded (tract mismatch / multi-copy). `examples/
str_calibrate.rs` is the harness. `str_concordance` joins caller→convention→imported vendor profile;
**34/34 calibrated markers agree** on James's GRCh38 chrY. CHM13 support via the lifted reference
(b9a7eed) — validated on his CHM13 CRAM, offsets build-independent. **UI wired** (c142ae7): Y-DNA tab
"Y-STR from sequence (HipSTR)" card. *Remaining:* widen the offset table with the ~300-kit CHM13 corpus
(+ harness QC: build/naming, BAM↔CSV swap detection, coverage) when it lands on the NAS; multi-copy/
nested aggregation (DYS385/DYS464/CDY/DYS389II — mostly uncallable anyway, tracts exceed read length);
auto-download the GRCh38 HipSTR BED. (Full status in project memory `str-caller`.)

> **History:** a 2026-06-12 length÷period port over *feature regions* (a bundled 24-locus catalog)
> was reverted — feature coords aren't tight tracts, so it was systematically offset. **Resolved
> 2026-06-16** by using the HipSTR reference's tight, end-inclusive tracts + the enclosing-read CIGAR
> model (measure each read against the known ref allele). HiFi (~4×) still keeps most loci LOW; the
> value is highest on 20–30× short-read WGS (validated there).

## 2. Y-chromosome profile management (variant-level, multi-source) — **DONE** (was MISSING)

The multi-source Y variant profile (combine WGS/chip/STR/private observations of the same position,
quality-weighted consensus, provenance, conflict detection, persistence) **is now built**:

| Capability | Scala | Status |
|---|---|---|
| Y-variant concordance / quality-weighted consensus (callable-state + region modifiers) | `yprofile/concordance/YVariantConcordance.scala` | **DONE** — `navigator-domain/consensus.rs` (`reconcile`/`obs_weight`) + `yprofile.rs` adapter; weighted by SourceType × depth × mapq × callable × region modifier |
| Y-profile persistence (profile, sources, variant calls, novel) | `HaplogroupProcessor.populateYProfile` | **DONE** — `consensus_profile` table (mig 0022), `build_y_profile`/`cached_y_profile`; private-Y union persisted |
| Genome-level consensus **placement** (pool all sources, place once) | — | **DONE** (fd599d9) — `place_y_consensus`, one tree/coord, authoritative terminal |
| Y-profile source-type weighting (method tiers × SNP/STR weights) | `yprofile/YProfileSourceType` | **DONE** — `SourceType::snp_weight` |
| Y-SNP profile **comparison / FTDNA-Big-Y-style match list** (cross-subject) | `yprofile/YProfileService.scala` | **DONE** (e0e44bf) — `navigator-domain/ymatch.rs` shared-derived/novel ranking; app `y_matches` one-vs-all over cached profiles; Y-DNA-tab match card |
| Y-STR genetic distance / TMRCA / match ranking | `yprofile/YProfileService.scala` | **DONE** (e0e44bf) — STR-GD via `values_match` (multi-copy aware); SNP TMRCA (~83 yr/SNP) + STR TMRCA (stepwise), shown with caveats |

**DONE (e0e44bf):** cross-subject Y matching — `navigator-domain/ymatch.rs` (pure: shared-derived/
novel set ranking, divergence = longest-common-prefix of lineages, STR-GD via `values_match`, SNP +
STR TMRCA) + app `y_matches(query, project_id)` (one-vs-all over **cached** profiles, no re-genotype;
tree-fetch degrades gracefully) + worker/UI "Y matches (workspace)" card. Local-only v1; the
federated match surface stays under §4 (the AppView IBD hub).

## 3. Vendor & mtDNA data import — **DONE** (committed; §3 complete per memory)

| Capability | Scala | Status |
|---|---|---|
| Vendor VCF import (FTDNA Big Y, YSEQ) + metadata/source tagging | `analysis/VcfCache.scala` | **DONE** (`##source=aengine` tagging; generic VCF → variant import) |
| mtDNA FASTA → variants vs rCRS → haplogroup | `analysis/MtDnaFastaProcessor.scala` | **DONE** (`mtvariants::derive`; FASTA import on add) |
| Vendor mtDNA FASTA import (FTDNA mtFull, YSEQ) | `analysis/MtDnaFastaProcessor.scala` | **DONE** |
| Chip genotypes → autosomal ancestry | `analysis/ChipAncestryAdapter.scala` | **DONE** (3eb9a07; liftover GRCh37→CHM13 AIMs; ~99% EUR on real 23andMe) |
| Pre-computed metrics importers — flagstat, Picard CollectWgsMetrics/AlignmentSummaryMetrics | `analysis/MetricsFileLoader.scala` | **DONE** (a92d29e) |
| Vendor test-type ID from BAI coverage shape (Big Y / Y Elite enrichment) | — | **DONE** |

**Remaining:** only the *import-UX dialogs* (multi-file / drag-drop / vendor-VCF/FASTA pickers) — the
backend ingests all of these, but the GUI is single-file-add only → tracked under §8.

## 4. Federated IBD matching — **DONE** (transport + payload + UI, validated live; AppView-side indexing remains)

Re-derived against the **AppView-mediated** model (not the Scala P2P relay). Device-key identity +
the encrypted edge-to-edge channel are in place and **validated end-to-end against a live AppView**
(localhost:9000, 2026-06-17): two `did:key` peers completed publish → request → symmetric-blind
discovery → dual consent → X3DH-lite handshake → AES-GCM round-trip. The whole transport works.

| Capability | Scala equiv | Status |
|---|---|---|
| Ed25519 device-key signing (published, verified via did:key) | — | **DONE** (a289f7a) `navigator-sync/device_key.rs`; **did:key self-cert validated live** |
| **did:key local identity (desktop bootstrap — no PDS/OAuth)** | — | **DONE** (1e43f12) `App::use_local_identity` / `set_active_account`; `ensure_device_key` skips PDS publish for did:key |
| Encrypted channel (X25519 IK/EK + X3DH-lite → HKDF → AES-256-GCM Envelope) | `ibd/IbdCryptoService.scala` | **DONE** (df988c0) `navigator-sync/exchange.rs` + app `/exchange/*` driver; **live round-trip confirmed** |
| Direct `exchange_request(partner_did)` initiator | — | **DONE** (1e43f12) — the counterpart to suggestion-based `ibd_introduce`; live request→consent→session confirmed |
| Signed AppView `ibd_suggestions` / `ibd_introduce` + UI card | `IbdMatchingCoordinator.scala` | **DONE** (a289f7a); `ibd_suggestions` signed poll validated live |
| Local refinements: chip IBD panel, real genetic map, compare over panel | `PairwiseIbdDetector.scala` | **DONE** (533b6ef/02cdfb7); consensus-driven `compare_ibd_consensus` |
| IBD **segment detection over the exchange** + attestation sign/verify | `IbdAttestation.scala` | **DONE** (816fcea) — `ibd_attest.rs` (IbdSite/IbdAttestation/IbdExchangeMsg) + app `exchange_ibd`; dosages over the channel → symmetric detect → signed dual-hash attestations. **Validated live** (two did:key peers: identical 75 cM / ThirdCousin, both verified + agreed) |
| Consent / request / result storage + lifecycle | `Match{Consent,Request,Result}Repository.scala` | **DONE** (02efee5) — `ibd_exchange_result` (mig 0025) persists results + both attestations; `exchange_ibd_for_subject` (real panel dosages, gzip + decimation to fit the 1 MiB envelope) + best-effort PDS publish; **Encrypted-exchange UI** (consent → run → results) + a did:key local-identity button. Validated live on James's 1.23M real sites |
| ROH detection, HalfSibling category | `RelationshipEstimator.scala` | PARTIAL |

**Done (1e43f12, 816fcea, 02efee5), all validated live:** transport (did:key auth + X3DH + AES-GCM),
segment exchange + dual-hash attestation, real-data resolver (`ibd_dosages_for_subject`, panel-
restricted, gzip + decimation to fit the 1 MiB envelope), persistence (`ibd_exchange_result`,
mig 0025) + best-effort PDS publish, and the **Encrypted-exchange UI** + did:key bootstrap button.
**Remaining (Navigator-side, small):** the AppView consuming the published attestations (Jetstream
index → discovery), and surfacing the initiator side in the UI (today the UI consents/runs; requests
are opened via the existing network-suggestions `ibd_introduce`). The PII-posture decision still
governs what content rides the channel.

## 5. Sync durability — **DONE** (Phase 1 + Phase 2; live-PDS validation pending)

| Capability | Scala | Status |
|---|---|---|
| Persistent outbox (survive restart/offline, batched, backoff-capped) + drain | `SyncQueueRepository.scala`, `AsyncSyncService.scala` | **DONE** (7213269) `sync_outbox` (mig 0021) + background drain |
| Sync history / audit trail | `SyncHistoryRepository.scala` | **DONE** (`sync_history`; `record_dir` adds PULL/CONFLICT) |
| **Idempotent publish** (no duplicate PDS records on re-publish) | — | **DONE** (e38f0b3) — `sync_state` (mig 0026) keeps the PDS-assigned TID; drain `putRecord`s at it |
| Conflict detection + resolution (local↔remote divergence, LWW) | `SyncConflictRepository.scala`, `PdsSyncValidation.scala` | **DONE** (e38f0b3) — pure `sync_reconcile::plan` (CID compare, remote-wins LWW, conflict logged) |
| PULL (ingest own PDS records back; reconcile) | — | **DONE** (e38f0b3) — `list_records`/`pull_list` + `pull_sync`; apply remote→local for biosample (the PII-free summaries limit apply to the fields they carry); derived summaries tracked |
| Source-file tracking by checksum (stable identity if path moves) | `SourceFileRepository.scala` | **DONE** (e38f0b3) — `source_file` (mig 0027) + registration on content-hash + `verify_source_files` |
| Per-entity at-uri/at-cid columns | `Repository.scala` | **DONE** — `sync_state` carries rkey/at_uri/at_cid/payload_hash per entity |

**Done (e38f0b3):** idempotent publish, PULL reconcile, conflict detection, `source_file`, per-entity
at-uri. **NB live-PDS validation deferred** — putRecord/listRecords need a real did:plc repo (did:key
has none; dev PDS down); the reconcile logic stands on unit tests. The fed records are PII-free
*summaries* (no local guid), so remote→local apply is inherently limited to the fields they carry — a
real per-entity remote→local mirror would need richer records (an AppView-contract decision).

## 6. Analysis caching/resume & report completeness — **PARTIAL**

| Capability | Scala | Status |
|---|---|---|
| Diploid / indel variant calling, whole-genome diploid VCF | `WholeGenomeVariantCaller.scala` | **DONE** (598226a/50ecaa7/06fae27); + consensus joint genotype (98adfe5) |
| Report/exports — TSV/HTML (ancestry, metrics), callable BED | `*ReportWriter` | **DONE** (ec3c3e1) |
| WGS-metrics completeness — MAD coverage, per-base exclusion fractions | `WgsMetricsProcessor.scala` | PARTIAL→mostly (83ea6d0: MAD + pct_exc_mapq/baseq; dup/unpaired/overlap/capped deferred) |
| **Multi-step checkpoint/resume** (skip completed steps; BAM-mtime invalidation) | `AnalysisCheckpoint.scala`, `AnalysisCache.scala` | **DONE** (192a939) — artifacts carry a `source_sig` (BAM/CRAM `mtime:size`, mig 0024); `load_analysis` invalidates on source change so every `cached_*` is stale-aware; `run_sv`/`run_denovo` now cache-first → full-analysis resumes fresh steps. *Remaining:* content-hash (vs mtime) option; a "resumed N/5" progress readout |
| Multiallelic indel calling, left-normalization edge cases | — | PARTIAL (multiallelic SNV done; multiallelic indel deferred) |
| Callable-loci **SVG track** + haplogroup-report CSV | `BioVisualizationUtil.scala` | MISSING (BED + tables done; no SVG) |

**Remaining (small-medium):** cross-step checkpoint/resume + BAM-mtime invalidation (high value now
that consensus/diploid passes are heavy — avoids recompute); SVG track; multiallelic indels.

## 7. refgenome breadth — **PARTIAL**

| Capability | Scala | Status |
|---|---|---|
| Genome-region API + 2-layer cache (centromere/telomere/Y regions, offline fallback) | `refgenome/GenomeRegionService.scala` | **DONE** (99351de) + ideogram tab (5682976) |
| Full Y-region annotation (PAR/ampliconic/palindrome + quality modifiers) | `refgenome/YRegion*.scala` | **PARTIAL→mostly** (4de7ff8: PAR/palindrome/amplicon/heterochromatin + modifier ladder; XTR/STR/centromere data still thin) |
| Asset integrity (sha256 manifest verify) for ancestry/IBD assets | — | **DONE** (4ec09be) |
| Genotype liftover (single-position SNP/STR, strand-flip + rev-comp) | `liftover/GenotypeLiftover.scala` | PARTIAL (haplo placement lifts via du-bio; no general batch API) |
| **VCF liftover orchestration** (contig UCSC↔NCBI norm, PAR filtering, REF/ALT swap recovery) | `liftover/LiftoverProcessor.scala` | **MISSING** |
| **Reference download checksum/integrity verification** | `refgenome/ReferenceGateway.scala` | **MISSING** (asset manifests verified; raw reference/chain downloads not) |

**Remaining (low-medium):** VCF liftover orchestration + reference-download checksums (both mostly
STR/VCF-workflow enablers).

## 8. UI — **PARTIAL** (settings + ideogram + painting + consensus tabs landed)

| Capability | Scala | Status |
|---|---|---|
| Settings/preferences dialog (tree-provider, reference config, cache) | `SettingsDialog.scala` | **DONE** (d221b34) |
| Ideogram / cytoband visualization | `YChromosomeIdeogramPanel.scala` | **DONE** (5682976) karyotype ideogram tab |
| Chromosome painting track (local ancestry) | — | **DONE** diploid two-copy painting (consensus-driven) |
| Consensus-driven Y/mt/Autosomal/Ancestry/IBD tabs | — | **DONE** (this arc) |
| **Batch / project-bundle / vendor-VCF / vendor-FASTA import dialogs** (multi-file, drag-drop, auto-detect) | `{BatchImport,ProjectImport,ImportVendorVcf,ImportVendorFasta}Dialog.scala` | **DONE** (59b5696) — per-subject Add Data is now a multi-file + folder picker; drag-drop routes files+folders through one batch (`add_data_batch`, auto-detect each) → import-summary modal. Project-bundle import (`import_project_dir`) already had a folder picker. *Remaining:* explicit vendor presets (FTDNA Big Y / mtFull labels) are cosmetic — auto-detect already routes those formats |
| Y-profile management/detail + source-reconciliation dialogs | `YProfile*Dialog.scala`, `SourceReconciliationPanel.scala` | PARTIAL (Y-profile card + consensus block exist; no dedicated management/audit dialog) |
| IBD match-detail browser — chromosome ideogram with segment painting; segment CSV export | `MatchDetailDialog.scala`, `ChromosomeBrowserPanel.scala` | MISSING (downstream of §4) |
| **PCA scatter** (PC1×PC2 projection plot) | `ui/…` | **MISSING** (loadings/projection computed; donut + composition + map exist; no scatter widget) |
| Haplogroup report dialog (scored candidates / lineage / SNPs / private) | `HaplogroupReportDialog.scala` | PARTIAL (Y-DNA tab shows terminal + branches; no full scored-candidate dialog) |
| Fingerprint-match / merge-sequence-runs dialogs | `{FingerprintMatch,MergeSequenceRuns}Dialog.scala` | MISSING |

**Remaining:** the **import dialogs** are the standout — every backend (batch `import_project_dir`,
vendor VCF, mtFull FASTA, sidecar/Picard) already exists, so this is pure UX surfacing of shipped
capability. Plus PCA scatter (small), Y-profile management dialog, IBD match browser (→§4),
fingerprint/merge dialogs.

---

## Priority summary (2026-06-17)

| # | Subsystem | Impact | Size | Notes |
|---|---|---|---|---|
| ~~8-import~~ | ~~Import UX dialogs~~ | — | — | **DONE 59b5696** — multi-file + folder Add Data, drag-drop, auto-detect, summary modal |
| ~~6-resume~~ | ~~Analysis checkpoint/resume~~ | — | — | **DONE 192a939** — source-sig invalidation + cache-first SV/denovo |
| 1b-caller | ~~STR calling from sequence~~ | — | — | **DONE 986e00b** — enclosing-read genotyper + HipSTR-reference parse, validated on GRCh38 chrY |
| 1b-vendor | STR vendor bridge — convention layer + concordance + CHM13 ref | High | — | **DONE** (5fc6641 convention layer; b9a7eed cross-build lift; c142ae7 UI; **b631d79 216-kit CHM13 recalibration**). Offset table rebuilt on 216 CHM13 Big Y kits (swap-QC + per-kit panic isolation in the harness); 6 build-dependent markers handled via `StrBuild`+`GRCH38_DELTA`. Validated held-out: CHM13 1001615 44/44+14/14, GRCh38 27520 55/55+15/15, zero mismatches. *Remaining (low):* the CHM13 lift dropped 33 named chrY markers (incl DYS19/391/426 — table retains their GRCh38 values for the BAM path); multi-copy aggregation; auto-download |
| ~~2-match~~ | ~~Cross-subject Y matching~~ | — | — | **DONE e0e44bf** — `ymatch.rs` (shared-SNP/novel ranking, divergence LCA, STR-GD, SNP+STR TMRCA) + app `y_matches` one-vs-all + Y-DNA-tab match card. Local v1; federated surface under §4 |
| ~~5-p2~~ | ~~Sync conflict detection + PULL + `source_file`~~ | — | — | **DONE e38f0b3** — idempotent publish (sync_state keeps the PDS TID → putRecord), pure reconcile planner + pull_sync, source_file (mig 0027). Live-PDS validation pending (needs did:plc repo) |
| ~~4-live~~ | ~~IBD live exchange~~ | — | — | **DONE** (1e43f12, 816fcea, 02efee5) — transport + segment payload + attestation + real-data resolver + persistence + Encrypted-exchange UI + did:key bootstrap, validated live (James's 1.23M sites → ParentChild, verified+agreed). Remaining: AppView attestation indexing |
| 7 | VCF liftover orchestration + reference-download checksums | Low-Med | Medium | STR/VCF-workflow enablers |
| 8-misc | PCA scatter, Y-profile management dialog, IBD match browser, fingerprint/merge dialogs | Low-Med | Mixed | Several small; IBD browser downstream of §4 |

**Recently shipped:** import UX (59b5696), checkpoint/resume (192a939), **STR caller foundation**
(986e00b — the hard, twice-attempted part), **STR vendor bridge fully landed** (b631d79 — 216-kit
CHM13 recalibration + build-aware offsets), **§2 cross-subject Y matching** (e0e44bf), and **§4
federated IBD end-to-end** (1e43f12 / 816fcea / 02efee5 — transport, attestation, real-data exchange,
persistence, UI; validated live), and **§5 sync durability Phase 2** (e38f0b3 — idempotent publish +
PULL reconcile + source_file). **Best next steps:** the **§8-misc** small UI items (PCA scatter, IBD
match browser), §7 VCF-liftover orchestration + reference checksums, or §3-style polish. The big
functional gaps (STR, cross-subject Y, federated IBD, sync) are now all closed.
**Live-PDS-gated:** §5-p2's putRecord/listRecords + §4's attestation publish both need a running dev
PDS + a did:plc account to validate end-to-end; until then they rest on unit tests + the did:key
AppView validation.
