# Atmosphere Lexicon - Executive Summary

Current status and milestones for the DecodingUs Atmosphere Lexicon implementation across all teams.

**Last Updated:** 2025-12-09

---

## Overview

The Atmosphere Lexicon defines decentralized, user-owned genomic records for the AT Protocol (Bluesky) ecosystem. This enables citizens to own their genomic data in Personal Data Stores (PDS) while DecodingUs operates as an AppView for network-wide aggregation and analysis.

**Core Principle:** Raw genomic data (BAM, CRAM, VCF, FASTQ, genotype files) **never** leaves the user's device. All analysis is performed locally in Navigator Workbench. Only computed summaries and metadata flow through the PDS to DecodingUs.

---

## MVP Status

### **AppView MVP: SHIPPABLE**

**Status Date:** 2025-12-09

> **⚠️ Superseded by the 2026-06 scope reduction.** The AppView is **no longer a network
> mirror**. The full-CRUD firehose ingestion and per-collection event handlers described
> below are being **removed** from the decodingus codebase. The AppView's role narrows to
> (1) the **known-variant catalog** via direct proposal submission and (2) **on-demand
> coverage aggregation** from Navigator-published summary records. Per-sample data lives in
> the researcher's Navigator workspace. See [08-AppView-Lifecycle.md](./08-AppView-Lifecycle.md).
> The completed-handler status below is retained as historical record.

The DecodingUs AppView backend has reached MVP status for Phases 1 and 2. All core record types are fully implemented with event handlers ready to process inbound firehose events.

| Component | Status | Notes |
|:----------|:-------|:------|
| Lexicon Schema Definitions | ✅ Complete | All record types defined in v1.8 |
| Database Migrations | ✅ Complete | Migrations 37-40 applied |
| Domain Models (Scala) | ✅ Complete | JSONB consolidation for Slick 22-tuple limit |
| DAL Tables (Slick) | ✅ Complete | All tables with nested projections |
| Repositories | ✅ Complete | Full CRUD for all MVP entities |
| Event Handlers | ✅ Complete | All core + extended record types |
| Firehose Controller | ✅ Complete | JSON discriminator-based routing |

### MVP Scope

**Included:**
- Inbound event processing (Create/Update/Delete) for all core record types
- Full data persistence with AT URI/CID tracking
- Optimistic locking and conflict detection
- Soft deletes with orphan handling

**Post-MVP (Backlog):**
- REST API query endpoints (back-flow channels)
- Integration test suite expansion
- OAuth client auth + standard-relay firehose subscription

---

## Record Implementation Status

| Record Type | Schema | DAL | Repository | Handler | Status |
|:------------|:-------|:----|:-----------|:--------|:-------|
| `biosample` | ✅ | ✅ | ✅ | ✅ | **MVP** |
| `sequencerun` | ✅ | ✅ | ✅ | ✅ | **MVP** |
| `alignment` | ✅ | ✅ | ✅ | ✅ | **MVP** |
| `project` | ✅ | ✅ | ✅ | ✅ | **MVP** |
| `genotype` | ✅ | ✅ | ✅ | ✅ | **MVP** |
| `populationBreakdown` | ✅ | ✅ | ✅ | ✅ | **MVP** |
| `haplogroupReconciliation` | ✅ | ✅ | ✅ | ✅ | **MVP** |
| `strProfile` | ✅ | 📋 | 📋 | 📋 | Future |
| `matchConsent` | ✅ | 📋 | 📋 | 📋 | Future |
| `matchList` | ✅ | 📋 | 📋 | 📋 | Future |
| `instrumentObservation` | ✅ | 📋 | 📋 | 📋 | Future |
| `imputation` | ✅ | 📋 | 📋 | 📋 | Future |

---

## Team Milestones

### DecodingUs (AppView Backend)

| Milestone | Status | Description |
|:----------|:-------|:------------|
| Core Record Schema | ✅ Complete | `biosample`, `sequencerun`, `alignment`, `project`, `workspace` |
| Firehose Event Handlers | ✅ Complete | Full CRUD for all core + new record types |
| Haplogroup Reconciliation | ✅ Complete | Multi-run consensus, conflict resolution, audit trail |
| Genotype Record Schema | ✅ Complete | Multi-test-type support with taxonomy codes |
| Population Breakdown Schema | ✅ Complete | 33 populations, 9 super-populations, PCA coordinates |
| Database Tables | ✅ Complete | `genotype_data`, `population_breakdown`, `haplogroup_reconciliation` |
| Atmosphere Records (Scala) | ✅ Complete | All record types in `AtmosphereRecords.scala` |
| Repositories | ✅ Complete | Full CRUD for all MVP entities |
| Event Handler Routing | ✅ Complete | `AtmosphereEventHandler` routes all events |
| **MVP Release** | ✅ **SHIPPABLE** | Ready for Phase 1/2 integration |

**Next Focus:** Integration testing, then OAuth client auth + standard-relay firehose subscription (see [11-Auth-and-Permissions.md](./11-Auth-and-Permissions.md)).

---

### Navigator Workbench (Edge App)

| Milestone | Status | Description |
|:----------|:-------|:------------|
| Chip File Parsing | 🚧 In Progress | 23andMe, AncestryDNA, FTDNA, MyHeritage, LivingDNA |
| Haplogroup Calling (Chip) | 🚧 In Progress | Y-DNA and mtDNA from ~3-4K chip markers |
| Ancestry Analysis | 🚧 In Progress | PCA projection + GMM onto 1000G + HGDP reference |
| PDS Sync (Genotype) | 📋 Planned | Sync genotype metadata to user's PDS |
| PDS Sync (Ancestry) | 📋 Planned | Sync population breakdown to user's PDS |
| Multi-Run Reconciliation | 📋 Planned | Local reconciliation UI and logic |

**Current Focus:** Multi-test-type genotype parsing and ancestry analysis pipeline.

---

### Nexus (BGS Node)

| Milestone | Status | Description |
|:----------|:-------|:------------|
| WGS Pipeline | ✅ Complete | FASTQ → BAM/CRAM → VCF pipeline |
| Haplogroup Calling (WGS) | ✅ Complete | Full Y-DNA/mtDNA SNP-based calling |
| Biosample Sync | ✅ Complete | Push biosample metadata to DecodingUs |
| Sequence Run Sync | ✅ Complete | Push sequencing metadata to DecodingUs |
| Alignment Metrics Sync | ✅ Complete | Push coverage/quality metrics |
| AT Protocol Integration | 📋 Planned | OAuth-scoped direct PDS writes |

**Current Focus:** Production stability and OAuth-scoped direct-to-PDS writes (Kafka integration dropped — superseded by OAuth permission sets).

---

## Integration Phases

> **Revised (OAuth/permissions landed).** AT Protocol now supports granular,
> per-collection write authorization via OAuth permission sets. This removes the reason
> the REST/Kafka relay existed — a backend holding the user's full-access app password.
> Clients can now write **directly to the user's PDS** under a narrow scope. The phased
> relay plan below is collapsed accordingly. See
> [11-Auth-and-Permissions.md](./11-Auth-and-Permissions.md).

### Phase 1: MVP (Legacy bootstrap) - READY

- BGS Node → REST API → DecodingUs
- Navigator → REST API → DecodingUs
- Retained only as a bootstrap/import path; not the target architecture
- No PDS integration

### ~~Phase 2: Hybrid (Kafka)~~ - CUT

- **Removed.** Kafka was a scaled version of the credential-holding relay. With OAuth
  direct-to-PDS writes, the intermediary is unnecessary and is dropped from the plan.

### Phase 2 (was Phase 3): Full Atmosphere (AppView) - Target

- Clients authenticate via **OAuth** and request the `com.decodingus.atmosphere.navigatorCore`
  permission set (narrow per-collection write scopes) — app passwords deprecated
- All clients write directly to the user's PDS under that scope
- DecodingUs ingests via a **standard relay / Jetstream** subscription (no custom relay
  infrastructure to build)
- Full record compliance with this Lexicon

> **Firehose note:** the permission spec covers *writes* only — reads and subscriptions
> are explicitly out of scope, so the AppView's firehose ingest path remains required.
> "Custom firehose" referred to the REST/Kafka relay, which is what we are removing.

---

## Key Schema Versions

| Version | Date | Changes |
|:--------|:-----|:--------|
| 1.5 | 2025-12-08 | Multi-run reconciliation (`haplogroupReconciliation`) |
| 1.6 | 2025-12-08 | Enhanced ancestry: 33 populations, 9 super-populations |
| 1.7 | 2025-12-08 | Multi-test-type: `testTypeCode` taxonomy |
| 1.8 | 2025-12-09 | AppView implementation complete |
| 1.9 | 2025-12-09 | MVP marked shippable |
| **2.0** | **2026-06-01** | **OAuth/permissions: direct-to-PDS writes, Kafka relay cut, app passwords deprecated** |

---

## Reference Documents

| Document | Location | Purpose |
|:---------|:---------|:--------|
| Atmosphere Lexicon | `documents/atmosphere/` | Full schema specification |
| Multi-Test-Type Roadmap | `documents/multi-test-type-roadmap.md` | Genotype support planning |
| Ancestry Analysis | `documents/AncestryAnalysis.md` | PCA/GMM algorithm details |
| Multi-Run Reconciliation | `documents/MultiRunReconciliation.md` | Haplogroup consensus planning |
| IBD Matching System | `documents/ibd-matching-system.md` | Match system planning |
| Edge Client Status | `documents/Edge_Client_Implementation_Status.md` | Navigator implementation tracking |
| Auth & Permissions | `documents/atmosphere/11-Auth-and-Permissions.md` | OAuth migration, permission set, relay removal |

---

## Contact

- **DecodingUs Backend:** [Backend Team]
- **Navigator Workbench:** [Navigator Team]
- **Nexus BGS Node:** [Nexus Team]
