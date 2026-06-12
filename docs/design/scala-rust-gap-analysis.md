# Scala → Rust gap analysis (functionality still missing)

A systematic pass over the legacy Scala/ScalaFX source (`src/main/scala/com/decodingus/`,
238 files / 66k LOC, still on this branch) against the Rust rewrite (`crates/`). This catalogs
**capabilities present in Scala that are absent or incomplete in Rust**, grouped by subsystem and
prioritized.

## Method & exclusions

Six parallel sweeps (analysis, haplogroup/Y/STR, ancestry/IBD, workspace/sync, UI, refgenome),
then reality-checked against the tree. **Deliberately excluded** (not gaps — by design):

- **External-tool wrappers.** The Rust engine is GATK/samtools/bcftools-free by design
  (`CLAUDE.md`). So `GatkRunner`, `GatkToolProcessor`, GATK `HaplotypeCaller`/`Mutect2`/`LiftoverVcf`
  orchestration, `samtools faidx`/`flagstat` shelling, etc. are **not** counted as missing — the
  Rust caller/walkers/noodles replace them. Where a wrapper produced a real *output capability* the
  Rust side lacks (e.g. diploid/indel calling, VCF liftover), that capability is listed below.
- **Already ported** (per project memory): coverage+callable (+per-contig histograms), read
  metrics, sex, SV, unified quality-metrics walker, haploid caller (force-call + de-novo SNP),
  heteroplasmy, header probe, lab/instrument/library-stats inference, pipeline-sidecar fast-path
  import, Y+mt haplogroup placement (parsimony guard, FTDNA+DecodingUs providers), **chip-raw-data
  haplogroups** (23andMe/AncestryDNA), BISDNA Y-SNP import, ancestry (AF-likelihood / PCA-GMM /
  nMonte / ADMIXTURE / local-ancestry painting — ~95%), refgenome retrieval+cache, Y/mt
  rotation-aware liftover, masked-rCRS build. The **local** IBD math (segment detection, genetic
  map, relationship estimation) is also ported.

Status legend: **MISSING** = no Rust equivalent · **PARTIAL** = some behavior, notable holes.

---

## 1. STR calling & reference — MISSING (largest clean gap)

The whole Y-STR pipeline *from sequence data* is absent. Only vendor-CSV parsing has landed
(`strprofile`, the wide FTDNA/YSEQ layout) — i.e. we can read someone else's STR results but can't
produce our own.

| Capability | Scala | Status |
|---|---|---|
| STR repeat-count **calling from BAM/CRAM** (CIGAR-aware span, stutter filtering, HIGH/MED/LOW tiers) | `analysis/StrCaller.scala` | MISSING |
| STR reference gateway — HipSTR BED download + per-build liftover + cache | `refgenome/StrReferenceGateway.scala`, `StrReferenceCache.scala` | MISSING |
| STR annotator — position → locus name/period (indexed BED, binary search) | `refgenome/StrAnnotator.scala` | MISSING |
| STR panel config (FTDNA/YSEQ tiers, marketing-vs-actual counts, multi-copy markers, cumulative/exclusive) | `str/StrPanelConfig.scala` + `resources/str-panels.conf` | MISSING |
| STR panel classification + cumulative-threshold matching (FTDNA Y-12⊂Y-25⊂… vs YSEQ non-cumulative) | `str/StrPanelService.scala` | PARTIAL (CSV parse only) |
| STR marker comparison (Simple/MultiCopy/Complex value matching, conflict detect) | `str/StrMarkerComparator.scala` | MISSING |
| **Y-STR reporting tables (FTDNA/YSEQ style)** — summary card (provider toggle, tier badges, conflict badge) + detail dialog (All-Markers searchable table + By-Panel grouped Y-12/25/37/67/111 view, conflict highlighting) | `ui/components/YStrSummaryPanel.scala`, `YStrDetailDialog.scala` | MISSING |
| Y-STR genetic distance / TMRCA / match ranking | `yprofile/YProfileService.scala` | MISSING |

**Impact:** high — blocks native Y-STR profiling and STR-based matching for the WGS cohort. There is
a `feature/str-calling` branch from the Scala era; this is a coherent, self-contained subsystem to
port. The *caller* depends on the STR reference gateway (HipSTR BEDs, §7) landing first.

**Decoupled near-term win — STR reporting.** The FTDNA/YSEQ-style tables (summary card + the
By-Panel grouped marker view) and the underlying panel config (`str-panels.conf`: tier definitions,
marketing-vs-actual counts for multi-value markers like DYS385×2 / DYS464×4 / YCAII×2 / CDY×2,
cumulative + exclusive-marker rules) render from **vendor-CSV data we already parse** (`strprofile`).
They do **not** need the BAM caller or the STR reference gateway. So STR reporting + panel
classification can land independently and immediately, giving users a real Y-STR view from their
FTDNA/YSEQ exports while native calling is still pending.

## 2. Y-chromosome profile management (variant-level, multi-source) — MISSING

Distinct from haplogroup *placement* (done) and run-level haplogroup *reconciliation* (the
`0010_reconciliation` migration, present). This is the per-variant, multi-source Y profile: combine
Sanger/WGS/chip/STR observations of the same Y position, vote a consensus, track provenance/audit,
detect concordance conflicts, persist as a first-class profile.

| Capability | Scala | Status |
|---|---|---|
| Y-variant concordance / quality-weighted consensus (SNP + STR tiers, callable-state modifiers) | `yprofile/concordance/YVariantConcordance.scala` | MISSING |
| Y-profile persistence (profile, sources, variant calls, novel variants) + write API | `yprofile/…`, `HaplogroupProcessor.populateYProfile` | MISSING (no store tables) |
| Y-SNP profile comparison / FTDNA-Big-Y-style match list | `yprofile/YProfileService.scala` | MISSING |
| Y-profile source-type weighting (method tiers × SNP/STR weights) | `yprofile/YProfileSourceType` | PARTIAL (enum sketched in `navigator-domain/variants.rs`) |

**Impact:** high — this is the backbone for serious Y research (multi-test integration, private
variants, match attribution) and underpins much of the missing Y UI (§8). Large: ~13 Scala files +
dialogs.

## 3. Vendor & mtDNA data import — PARTIAL/MISSING

| Capability | Scala | Status |
|---|---|---|
| Vendor VCF import (FTDNA Big Y, YSEQ) with BED target regions + metadata tracking | `analysis/VcfCache.scala` (vendor ops) | MISSING |
| mtDNA FASTA → variant extraction vs rCRS → haplogroup/VCF | `analysis/MtDnaFastaProcessor.scala` | PARTIAL (`mtdna.rs` validates FASTA shape only; no variant calling) |
| Vendor mtDNA FASTA import (FTDNA mtFull, YSEQ) | `analysis/MtDnaFastaProcessor.scala` | MISSING |
| Chip genotypes → ancestry (autosomal) | `analysis/ChipAncestryAdapter.scala` | MISSING (chip → *haplogroup* is done; chip → *ancestry* is not) |
| Pre-computed metrics importers — samtools `flagstat`, Picard `CollectWgsMetrics`/`AlignmentSummaryMetrics` | `analysis/MetricsFileLoader.scala`, `FlagstatParser.scala` | MISSING (sidecar path covers samtools `coverage` + GATK `CallableLoci` only) |

**Impact:** medium-high — many users arrive with vendor exports (Big Y VCF, mtFull FASTA) rather
than BAMs; today those can't be ingested as variants.

## 4. Federated IBD matching — MISSING (local math present)

Only the *pure math* is ported (segment detection, genetic map, relationship estimation). The entire
**federated matching protocol** is absent: ECDH (X25519) key exchange, AES-256-GCM payload
encryption, Ed25519 attestations, the WebSocket relay client, the protocol state machine, and the
consent/request/result storage + service.

| Capability | Scala | Status |
|---|---|---|
| IBD crypto (X25519 / AES-256-GCM / Ed25519) | `ibd/IbdCryptoService.scala` | MISSING |
| Match protocol coordinator + state machine + relay client | `ibd/IbdMatchingCoordinator.scala`, `MatchingProtocolState.scala`, `IbdRelayClient.scala` | MISSING |
| Attestation sign/verify + deterministic match-hash | `ibd/IbdAttestation.scala`, `MatchSummary.scala` | MISSING |
| Consent / request / result storage + service | `repository/Match{Consent,Request,Result}Repository.scala`, `service/IbdMatchService.scala` | MISSING (no tables) |
| Local refinements: ROH detection, HalfSibling category, HapMap genetic-map loader, compact binary encode | `ibd/PairwiseIbdDetector.scala`, `RelationshipEstimator.scala`, `GeneticMap.scala` | PARTIAL |

**Architectural caveat:** the federated design was deliberately reworked toward an **AppView-mediated**
model (see project memory: writes-only spec, relay cut, AppView is the IBD-suggestion hub). The Scala
peer-to-peer crypto/relay stack should **not** be ported verbatim — re-derive the consent/request/
result + attestation pieces against the current AppView design. The *local* refinements (ROH,
HalfSibling, HapMap loader) are straight ports.

## 5. Sync durability — MISSING

`navigator-sync` has the publish engine + retry, but only in-memory. The Scala durability layer is
absent.

| Capability | Scala | Status |
|---|---|---|
| Persistent outbox/queue (survive restart/offline→online, batched, backoff-capped) | `repository/SyncQueueRepository.scala`, `sync/AsyncSyncService.scala` | MISSING |
| Conflict detection + resolution state (local↔remote divergence, suggested strategy, snapshots) | `repository/SyncConflictRepository.scala`, `pds/PdsSyncValidation.scala` | MISSING |
| Sync history / audit trail (direction, result, version/CID before-after) | `repository/SyncHistoryRepository.scala` | MISSING |
| Source-file tracking by checksum (stable identity if path moves; analyzed-state) | `repository/SourceFileRepository.scala` | PARTIAL (deferred content-hash exists, no `source_file` table) |
| Per-entity sync metadata (atUri/atCid, SyncStatus) columns | `repository/Repository.scala` | PARTIAL |

**Impact:** medium — edit-offline-then-reconnect currently loses pending syncs; no "why did sync
fail" visibility. Ties to the AppView design decision.

## 6. Analysis caching/resume & report completeness — PARTIAL

| Capability | Scala | Status |
|---|---|---|
| Multi-step checkpoint/resume (skip completed steps; BAM-mtime invalidation) | `analysis/AnalysisCheckpoint.scala`, `AnalysisCache.scala` | PARTIAL (artifact store caches per `algorithm_version`; no cross-step resume) |
| Diploid / indel variant calling, whole-genome diploid VCF | `analysis/WholeGenomeVariantCaller.scala` | MISSING (caller is haploid SNP-only by design; indels/diploid are a real capability gap) |
| WGS-metrics completeness — MAD coverage, exclusion fractions (`pctExcMapq/Dupe/Unpaired/Baseq/Overlap/Capped`) | `analysis/WgsMetricsProcessor.scala` | PARTIAL |
| Report/exports — text/HTML/TSV (ancestry, metrics), callable-loci BED + SVG track | `analysis/*ReportWriter`, `BioVisualizationUtil.scala` | PARTIAL (JSON only) |

## 7. refgenome breadth — PARTIAL

Core retrieval/cache/liftover is done; the breadth pieces below are missing.

| Capability | Scala | Status |
|---|---|---|
| Genotype liftover (single-position SNP/STR lift, strand-flip + rev-comp) | `liftover/GenotypeLiftover.scala` | PARTIAL (haplo placement lifts via du-bio; no general batch API) |
| VCF liftover orchestration (contig UCSC↔NCBI normalization, PAR filtering, REF/ALT swap recovery) | `liftover/LiftoverProcessor.scala` | MISSING |
| Genome-region API + 2-layer cache (centromere/telomere/Y regions, versioned, offline fallback) | `refgenome/GenomeRegionService.scala` | MISSING |
| Full Y-region annotation (PAR/XTR/STR/ampliconic/cytoband + quality modifiers + overlap consolidation) | `refgenome/YRegion{Gateway,Annotator}.scala` | PARTIAL (curated CHM13 masks: amplicon/palindrome/AZF only) |
| Download checksum/integrity verification | `refgenome/ReferenceGateway.scala` | MISSING |

## 8. UI — MISSING (much of it downstream of §1–4)

No charting/visualization beyond the new coverage histogram, and **no settings UI at all**.

| Capability | Scala | Status |
|---|---|---|
| **Settings/preferences dialog** — tree-provider select, reference config, cache management | `ui/components/SettingsDialog.scala`, `ReferenceConfigDialog.scala` | MISSING (config is startup/env/CLI only) |
| Batch import + project-bundle import + vendor VCF/FASTA import dialogs (multi-file, drag-drop, auto-detect) | `ui/components/{BatchImport,ProjectImport,ImportVendorVcf,ImportVendorFasta}Dialog.scala` | MISSING (single-file add only) |
| Y-profile management/detail + Y-STR detail dialogs; Y-ideogram; source-reconciliation panel | `ui/components/YProfile*Dialog.scala`, `YChromosomeIdeogramPanel.scala`, `SourceReconciliationPanel.scala` | MISSING (downstream of §2) |
| IBD match-detail browser — chromosome ideogram with segment painting; segment CSV export | `ui/components/MatchDetailDialog.scala`, `ChromosomeBrowserPanel.scala` | MISSING (downstream of §4) |
| **Y-STR reporting tables (FTDNA/YSEQ style)** — see §1; summary card + By-Panel grouped table. Renders from already-parsed vendor CSV, no backend dependency | `ui/components/YStrSummaryPanel.scala`, `YStrDetailDialog.scala` | MISSING |
| Haplogroup report dialog (scored candidates / lineage / SNPs / private), mtDNA-variants export, callable-loci result dialog | `ui/components/{HaplogroupReport,MtdnaVariants,CallableLociResult}Dialog.scala` | MISSING/PARTIAL |
| Fingerprint-match / merge-sequence-runs dialogs | `ui/components/{FingerprintMatch,MergeSequenceRuns}Dialog.scala` | MISSING |
| PCA scatter / chromosome-painting / ancestry-pie visualizations | `ui/…` | PARTIAL (donut + composition bar + map exist; no PCA scatter / painting track) |

---

## Priority summary

| # | Subsystem | Impact | Size | Notes |
|---|---|---|---|---|
| 1a | **STR reporting + panel classification** | Med-High | **Small** | **Near-term win** — FTDNA/YSEQ tables render from already-parsed vendor CSV (`strprofile`); no caller/reference dependency |
| 1b | STR calling & reference (from BAM) | High | Large | Self-contained; needs HipSTR BED gateway first. `feature/str-calling` branch exists |
| 2 | Y-profile (multi-source variant) management | High | Large | Backbone for Y research; unlocks much Y UI |
| 3 | Vendor/mtDNA import (Big Y VCF, mtFull FASTA, flagstat/Picard) | Med-High | Medium | Many users arrive with vendor exports |
| 8a | Settings/config UI | Med | Small | Cheapest high-visibility win — no settings dialog today |
| 4 | Federated IBD matching | Med | Large | **Re-derive vs AppView design, don't port P2P crypto/relay** |
| 5 | Sync durability (outbox/conflict/history) | Med | Medium | Ties to AppView design |
| 6 | Diploid/indel calling, metrics & report completeness | Med | Medium | Indels/diploid are the real caller gap |
| 7 | refgenome breadth (VCF liftover, region API, checksums) | Low-Med | Medium | Mostly STR/VCF-workflow enablers |

**Cheapest meaningful wins:** §8a settings/config UI and §1a STR reporting tables (both small, high
visibility, no backend dependency — STR reporting rides the vendor CSV we already parse).
**Biggest coherent feature:** §1b STR calling (+ §7 STR reference) — a clean, self-contained subsystem.
**Verify-before-building:** §4 IBD and §5 sync must be re-scoped against the current AppView-mediated
architecture, not ported verbatim from the Scala peer-to-peer design.
