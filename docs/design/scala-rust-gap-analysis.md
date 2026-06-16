# Scala → Rust gap analysis (functionality still missing)

A systematic pass over the legacy Scala/ScalaFX source (`src/main/scala/com/decodingus/`,
238 files / 66k LOC, still on this branch) against the Rust rewrite (`crates/`). This catalogs
**capabilities present in Scala that are absent or incomplete in Rust**, grouped by subsystem and
prioritized.

> **Revised 2026-06-16.** Much of the original (2026-06-12) backlog has since landed. This revision
> reconciles every section against the current tree + project memory, marks completed work with its
> commit, and re-derives the priority summary from what genuinely remains. Net: §1a, §3, §5(P1),
> §8a fully done; §2, §4, §6, §7, §8 substantially advanced. The standing big gap is **§1b STR
> calling from sequence**; the freshest remaining slices are **import UX dialogs (§8)**, **analysis
> checkpoint/resume (§6)**, and **sync conflict/PULL (§5 P2)**.

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

**Remaining (the vendor-grade bridge) — researched + validated, TABLED pending corpus on NAS:** the
DYS→coordinate mapping turned out to be **free** — the HipSTR BED already names ~206 chrY DYS markers
(`strref` parses the name). Validated against a real Big Y kit (FTDNA 27520, hg38): the caller
recovers ~100 single-copy markers, **59/102 exact** vs FTDNA; the rest are small fixed convention
offsets (DYS19 −1, DYS438 +2…) or large tract-mismatches to exclude (Y-GATA-H4, DYS484…). What's
left: a **corpus-calibrated offset table** (offsets are conventions only if constant across samples —
needs ~50 Big Y kits, not n=1), then the convention layer (name-normalize + single-copy offset +
multi-copy aggregation + DYS389 nesting) + §1a By-Panel report join + UI. Multi-copy/palindromic
(DYS385/DYS464/CDY) mostly not callable by the enclosing-read method (tracts exceed read length).
Plus a download/cache gateway + CHM13/GRCh37 ref or liftover. **Blocked on:** moving the few-hundred
Big Y corpus to the NAS. (Resume plan in project memory `str-caller`.)

> **History:** a 2026-06-12 length÷period port over *feature regions* (a bundled 24-locus catalog)
> was reverted — feature coords aren't tight tracts, so it was systematically offset. **Resolved
> 2026-06-16** by using the HipSTR reference's tight, end-inclusive tracts + the enclosing-read CIGAR
> model (measure each read against the known ref allele). HiFi (~4×) still keeps most loci LOW; the
> value is highest on 20–30× short-read WGS (validated there).

## 2. Y-chromosome profile management (variant-level, multi-source) — **PARTIAL** (was MISSING)

The multi-source Y variant profile (combine WGS/chip/STR/private observations of the same position,
quality-weighted consensus, provenance, conflict detection, persistence) **is now built**:

| Capability | Scala | Status |
|---|---|---|
| Y-variant concordance / quality-weighted consensus (callable-state + region modifiers) | `yprofile/concordance/YVariantConcordance.scala` | **DONE** — `navigator-domain/consensus.rs` (`reconcile`/`obs_weight`) + `yprofile.rs` adapter; weighted by SourceType × depth × mapq × callable × region modifier |
| Y-profile persistence (profile, sources, variant calls, novel) | `HaplogroupProcessor.populateYProfile` | **DONE** — `consensus_profile` table (mig 0022), `build_y_profile`/`cached_y_profile`; private-Y union persisted |
| Genome-level consensus **placement** (pool all sources, place once) | — | **DONE** (fd599d9) — `place_y_consensus`, one tree/coord, authoritative terminal |
| Y-profile source-type weighting (method tiers × SNP/STR weights) | `yprofile/YProfileSourceType` | **DONE** — `SourceType::snp_weight` |
| Y-SNP profile **comparison / FTDNA-Big-Y-style match list** (cross-subject) | `yprofile/YProfileService.scala` | **MISSING** — no shared-SNP / shared-novel match ranking between subjects |
| Y-STR genetic distance / TMRCA / match ranking | `yprofile/YProfileService.scala` | **MISSING** |

**Remaining (medium):** cross-subject Y matching — shared-derived/shared-novel ranking, genetic
distance / TMRCA. The single-subject profile backbone is complete; this is the *between-subjects*
layer (relates to §4 federated matching + the AppView IBD hub).

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

## 4. Federated IBD matching — **PARTIAL** (Phase 1+2 built; was MISSING)

Re-derived against the **AppView-mediated** model (not the Scala P2P relay). Device-key identity +
the encrypted edge-to-edge channel are in place.

| Capability | Scala equiv | Status |
|---|---|---|
| Ed25519 device-key signing (published, verified via did:key) | — | **DONE** (a289f7a) `navigator-sync/device_key.rs` |
| Encrypted channel (X25519 IK/EK + X3DH-lite → HKDF → AES-256-GCM Envelope) | `ibd/IbdCryptoService.scala` | **DONE** (df988c0) `navigator-sync/exchange.rs` + app `/exchange/*` driver |
| Signed AppView `ibd_suggestions` / `ibd_introduce` + UI card | `IbdMatchingCoordinator.scala` | **DONE** (a289f7a) |
| Local refinements: chip IBD panel, real genetic map, compare over panel | `PairwiseIbdDetector.scala` | **DONE** (533b6ef/02cdfb7); consensus-driven `compare_ibd_consensus` |
| IBD **segment detection over the exchange** + attestation sign/verify | `IbdAttestation.scala` | **MISSING/deferred** (channel exists; the match-segment exchange + attest not wired) |
| Consent / request / result storage + lifecycle | `Match{Consent,Request,Result}Repository.scala` | **MISSING** (no tables; AppView-side decision pending) |
| ROH detection, HalfSibling category | `RelationshipEstimator.scala` | PARTIAL |

**Remaining (medium-large):** the live two-peer exchange flow (segment exchange + attestation) and
consent/request/result lifecycle — both gated on a running AppView + the PII-posture decision.

## 5. Sync durability — **PARTIAL** (Phase 1 done; was MISSING)

| Capability | Scala | Status |
|---|---|---|
| Persistent outbox (survive restart/offline, batched, backoff-capped) + drain | `SyncQueueRepository.scala`, `AsyncSyncService.scala` | **DONE** (7213269) `sync_outbox` (mig 0021) + background drain |
| Sync history / audit trail | `SyncHistoryRepository.scala` | **DONE** (`sync_history`) |
| Conflict detection + resolution state (local↔remote divergence, strategy, snapshots) | `SyncConflictRepository.scala`, `PdsSyncValidation.scala` | **MISSING** |
| PULL (ingest own PDS records back; reconcile) | — | **MISSING** (publishes are write-only today) |
| Source-file tracking by checksum (stable identity if path moves) | `SourceFileRepository.scala` | PARTIAL (deferred content-hash; no `source_file` table) |
| Per-entity at-uri/at-cid columns | `Repository.scala` | PARTIAL (outbox carries at_uri; not per-entity) |

**Remaining (medium):** conflict detection + PULL + the `source_file` table. Ties to the AppView design.

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

## Priority summary (2026-06-16)

| # | Subsystem | Impact | Size | Notes |
|---|---|---|---|---|
| ~~8-import~~ | ~~Import UX dialogs~~ | — | — | **DONE 59b5696** — multi-file + folder Add Data, drag-drop, auto-detect, summary modal |
| ~~6-resume~~ | ~~Analysis checkpoint/resume~~ | — | — | **DONE 192a939** — source-sig invalidation + cache-first SV/denovo |
| 1b-caller | ~~STR calling from sequence~~ | — | — | **DONE 986e00b** — enclosing-read genotyper + HipSTR-reference parse, validated on GRCh38 chrY |
| 1b-vendor | STR vendor bridge — convention layer + concordance + CHM13 ref | High | — | **CORE DONE** (5fc6641 convention layer; b9a7eed cross-build lift). 34/34 calibrated markers agree on James's GRCh38 chrY; CHM13-lifted ref validated on his CHM13 CRAM (49 markers, offsets build-independent). **UI DONE** (c142ae7: Y-DNA tab "Y-STR from sequence" card). *Remaining:* multi-copy aggregation, widen offset table with the ~300-kit CHM13 corpus + harness QC, auto-download |
| 2-match | Cross-subject Y matching (shared-SNP/novel ranking, genetic distance/TMRCA) | Med | Medium | Between-subjects layer atop the done single-subject Y profile; relates to §4 |
| 5-p2 | Sync conflict detection + PULL + `source_file` table | Med | Medium | Ties to AppView design |
| 4-live | IBD live exchange (segment exchange + attestation) + consent/request/result | Med | Large | Gated on running AppView + PII-posture decision; don't port P2P verbatim |
| 7 | VCF liftover orchestration + reference-download checksums | Low-Med | Medium | STR/VCF-workflow enablers |
| 8-misc | PCA scatter, Y-profile management dialog, IBD match browser, fingerprint/merge dialogs | Low-Med | Mixed | Several small; IBD browser downstream of §4 |

**Recently shipped:** import UX (59b5696), checkpoint/resume (192a939), **STR caller foundation**
(986e00b — the hard, twice-attempted part). **Best next steps:** finish STR as a user feature
(**§1b-vendor**: DYS-name mapping + §1a panel integration + UI), or **§2-match cross-subject Y
matching** (medium, builds on the done Y profile). Both are medium and high-value.
**Biggest coherent feature:** **§1b STR calling** — still the standout functional gap, still hard.
**Verify-before-building:** §4 (IBD live) and §5-p2 (sync conflict/PULL) must be scoped against the
current AppView-mediated architecture, not ported verbatim.
