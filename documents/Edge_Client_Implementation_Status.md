# Edge Client Implementation Status

Navigator Desktop implementation status against the Atmosphere Lexicon specification.

**Overall Completion: ~55%**

Last updated: 2025-12-08

---

## Implementation Summary

### Core Records

| Record Type | Status | Notes |
|:------------|:-------|:------|
| `workspace` | ✅ Complete | Full CRUD, local JSON persistence, PDS sync via AT Protocol |
| `biosample` | ✅ Complete | All fields supported including haplogroups, refs to child records |
| `sequenceRun` | ✅ Complete | BAM/CRAM/FASTQ support, multi-platform detection, test type taxonomy |
| `alignment` | ✅ Complete | Full metrics, contig stats, reference build detection |
| `genotype` | 🚧 In Development | Chip parsing (23andMe, AncestryDNA, FTDNA, MyHeritage, LivingDNA), Y/mtDNA marker counts |
| `imputation` | ⬜ Planned | Not started |
| `project` | ✅ Complete | Sample grouping, metadata |

### Ancestry & Population Records

| Record Type | Status | Notes |
|:------------|:-------|:------|
| `populationBreakdown` | 🚧 In Development | PCA projection algorithm defined, 33 reference populations, awaiting reference data |

### Discovery & Inference Records

| Record Type | Status | Notes |
|:------------|:-------|:------|
| `instrumentObservation` | ⬜ Planned | Future crowdsourced lab discovery |

### STR Records

| Record Type | Status | Notes |
|:------------|:-------|:------|
| `strProfile` | ✅ Complete | Multi-panel support, complex STR values, FTDNA/YSEQ/WGS-derived sources |
| `haplogroupAncestralStr` | ⬜ Planned | Future ancestral STR reconstruction |

### IBD Matching Records

| Record Type | Status | Notes |
|:------------|:-------|:------|
| `matchConsent` | ⬜ Planned | Future IBD matching opt-in |
| `matchList` | ⬜ Planned | Future IBD match results |
| `matchRequest` | ⬜ Planned | Future match contact requests |

### Multi-Run Reconciliation Records

| Record Type | Status | Notes |
|:------------|:-------|:------|
| `haplogroupReconciliation` | ⬜ Planned | Multi-run conflict resolution |

### Backfeed Records (AppView → PDS)

| Record Type | Status | Notes |
|:------------|:-------|:------|
| `haplogroupUpdate` | ⬜ Planned | AppView haplogroup refinement notifications |
| `branchDiscovery` | ⬜ Planned | Novel branch discovery notifications |
| `treeVersionUpdate` | ⬜ Planned | Tree version change notifications |
| `potentialMatchList` | ⬜ Planned | Potential match candidates |
| `confirmedMatch` | ⬜ Planned | Confirmed match stamps |
| `syncStatus` | ⬜ Planned | PDS-AppView sync tracking |
| `updateDigest` | ⬜ Planned | Daily update digests |

### Common Definitions

| Definition | Status | Notes |
|:-----------|:-------|:------|
| `recordMeta` | ✅ Complete | Version tracking, timestamps |
| `fileInfo` | ✅ Complete | Metadata only (no file content transmission) |
| `haplogroupResult` | ✅ Complete | Full scoring with private variants |
| `privateVariants` | ✅ Complete | Novel variant tracking for branch discovery |
| `alignmentMetrics` | ✅ Complete | WGS metrics, callable loci |
| `strValue` (simple/complex) | ✅ Complete | Simple, multi-copy, and complex multi-allelic STRs |
| `populationComponent` | 🚧 In Development | 33 populations, super-population grouping |
| `pcaCoordinates` | 🚧 In Development | PCA projection for ancestry visualization |

### Analysis Capabilities

| Feature | Status | Notes |
|:--------|:-------|:------|
| Y-DNA Haplogroup Analysis | ✅ Complete | FTDNA + DecodingUs trees, two-pass calling, private variants |
| MT-DNA Haplogroup Analysis | ✅ Complete | FTDNA tree support |
| Private Variant Discovery | ✅ Complete | Integrated with haplogroup pipeline |
| WGS Metrics | ✅ Complete | Coverage, callable loci, contig stats |
| Library Stats | ✅ Complete | Platform/instrument detection, reference build |
| Chip Data Parsing | 🚧 In Development | 5 vendors supported, Y/mtDNA marker extraction |
| Ancestry Estimation | 🚧 In Development | PCA projection algorithm, awaiting reference panel |
| STR Import | ✅ Complete | Multi-vendor, multi-panel support |

### PDS Integration

| Feature | Status | Notes |
|:--------|:-------|:------|
| AT Protocol Authentication | 🚧 In Development | Feature-flagged, not yet enabled by default |
| Workspace Sync | 🚧 In Development | putRecord/getRecord implemented |
| Granular Record Sync | ⬜ Planned | Individual record CRUD not yet implemented |
| Firehose Publishing | ⬜ Planned | Event publishing for AppView consumption |

### Legend

- ✅ Complete - Feature implemented and working
- 🚧 In Development - Partially implemented or actively being developed
- ⬜ Planned - Designed in lexicon but not yet implemented in Edge client

---

## Features Remaining to Develop

### Priority 1: Complete In-Development Features

#### Ancestry Estimation (populationBreakdown)
**What it does:** Estimates population percentages using PCA projection against 1000 Genomes + HGDP/SGDP reference panels. Supports quick AIMs panel (~5k SNPs, 2-5 min) and detailed genome-wide analysis (~500k SNPs, 20-30 min).

**Remaining work:**
- Prepare reference data: Pre-compute PCA loadings and population centroids from reference panel
- Package reference data for download (~15MB AIMs, ~250MB genome-wide)
- Implement `AncestryReferenceGateway` for reference data download/caching
- Wire up UI in WorkbenchView for triggering analysis
- Generate HTML/JSON reports with pie charts and confidence intervals

**Files:** `src/main/scala/com/decodingus/ancestry/`

#### Chip Data Processing (genotype)
**What it does:** Parses raw genotype exports from consumer testing companies (23andMe, AncestryDNA, FTDNA, MyHeritage, LivingDNA). Extracts Y-DNA and mtDNA markers for haplogroup estimation.

**Remaining work:**
- Complete haplogroup estimation from chip Y-DNA/mtDNA markers
- Wire up ancestry analysis integration (chip → populationBreakdown)
- UI for chip file import and results display
- Workspace persistence for ChipProfile records

**Files:** `src/main/scala/com/decodingus/genotype/`

---

### Priority 2: Network Integration

#### Granular Record Sync
**What it does:** Syncs individual records (biosample, sequenceRun, alignment, etc.) to PDS instead of monolithic workspace sync. Enables fine-grained updates and reduced data transfer.

**Remaining work:**
- Implement per-record createRecord/putRecord/deleteRecord calls
- Track sync status per record (local-only, synced, pending)
- Handle conflict resolution when local and remote diverge
- Add SyncStatus model to workspace

#### Firehose Publishing
**What it does:** Publishes record changes to AT Protocol firehose for AppView consumption. Enables real-time network-wide features like branch discovery consensus.

**Remaining work:**
- Implement firehose event generation on record changes
- Define event payload schemas
- Add opt-in consent flow for network participation

---

### Priority 3: IBD Matching System

#### Match Consent (matchConsent)
**What it does:** Opt-in record for IBD matching participation. Controls what data is shared and matching thresholds.

**Remaining work:**
- Define consent UI with granular privacy controls
- Implement matchConsent record creation/management
- Integrate with ancestry analysis for population context

#### Match Discovery (matchList)
**What it does:** Stores potential IBD matches discovered by AppView. User explores matches locally before confirmation.

**Remaining work:**
- Implement matchList record consumption from AppView
- Build match explorer UI with segment visualization
- Calculate relationship estimates from shared cM

#### Match Confirmation (matchRequest, confirmedMatch)
**What it does:** Two-way handshake for confirming matches. Both parties must agree before contact info exchange.

**Remaining work:**
- Implement matchRequest creation and response flow
- Build confirmedMatch stamping when both agree
- Notification system for incoming match requests

---

### Priority 4: Multi-Run Reconciliation

#### Haplogroup Reconciliation (haplogroupReconciliation)
**What it does:** Reconciles haplogroup calls across multiple sequence runs (e.g., Big Y + WGS + chip). Identifies conflicts, resolves discrepancies, verifies sample identity.

**Remaining work:**
- Implement reconciliation algorithm comparing calls across runs
- Detect SNP conflicts and heteroplasmy observations
- Identity verification via autosomal concordance
- Build reconciliation report UI
- Integrate with existing haplogroup analysis pipeline

---

### Priority 5: AppView Backfeed

> **OAuth note:** "AppView → PDS" writes are no longer implicit. Writing into a user's
> repo now requires an explicit user-granted scope (`rpc`/service auth). Prefer modeling
> these as **AppView-owned records the client reads/subscribes to**, rather than the
> AppView mutating the user's repo. See `documents/atmosphere/11-Auth-and-Permissions.md` §5.

#### Haplogroup Update Notifications (haplogroupUpdate)
**What it does:** AppView notifies user when haplogroup assignment is refined based on new tree version or network consensus.

**Remaining work:**
- Subscribe to haplogroupUpdate records from AppView
- Display notifications in UI
- Optionally trigger re-analysis with updated tree

#### Branch Discovery (branchDiscovery)
**What it does:** AppView notifies when user's private variants contributed to new branch discovery.

**Remaining work:**
- Subscribe to branchDiscovery records
- Display discovery notifications with attribution
- Link to updated tree showing new branch

#### Tree Version Updates (treeVersionUpdate)
**What it does:** AppView notifies when haplogroup tree is updated.

**Remaining work:**
- Subscribe to treeVersionUpdate records
- Prompt user to re-run haplogroup analysis
- Show changelog of tree updates

---

### Priority 6: Future Capabilities

#### Imputation (imputation)
**What it does:** Imputes missing genotypes from chip data using reference panels. Enables ancestry analysis on sparse chip data.

**Scope:** Significant infrastructure - may require external imputation service or local Minimac4 integration.

#### Instrument Observation (instrumentObservation)
**What it does:** Crowdsourced sequencer/lab discovery. Users report instrument IDs and lab names to build community database.

**Scope:** Requires AppView aggregation logic and curator workflow.

#### Ancestral STR Reconstruction (haplogroupAncestralStr)
**What it does:** Reconstructs ancestral STR states at haplogroup branch points using network-wide STR profiles.

**Scope:** Requires significant network data aggregation. AppView feature primarily.

---

## Completion Metrics

| Category | Complete | In Development | Planned | Total |
|:---------|:--------:|:--------------:|:-------:|:-----:|
| Core Records | 5 | 1 | 1 | 7 |
| Ancestry Records | 0 | 1 | 0 | 1 |
| Discovery Records | 0 | 0 | 1 | 1 |
| STR Records | 1 | 0 | 1 | 2 |
| IBD Records | 0 | 0 | 3 | 3 |
| Reconciliation Records | 0 | 0 | 1 | 1 |
| Backfeed Records | 0 | 0 | 7 | 7 |
| Common Definitions | 5 | 2 | 0 | 7 |
| Analysis Capabilities | 6 | 2 | 0 | 8 |
| PDS Integration | 0 | 2 | 2 | 4 |
| **Total** | **17** | **8** | **16** | **41** |

**Percentage Complete:** 17/41 = 41% fully complete, 25/41 = 61% at least started

---

## Changelog

| Date | Changes |
|:-----|:--------|
| 2025-12-08 | Initial status document created from Atmosphere_Lexicon.md |
