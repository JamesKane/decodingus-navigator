# Atmosphere Lexicon - Edge App Reference

Navigator Desktop implementation guide for the Atmosphere Lexicon. For complete schema definitions, see
the [global Atmosphere Lexicon](../../decodingus/documents/atmosphere/).

## Overview

Navigator Desktop is an **Edge App** in the Atmosphere ecosystem - it performs local analysis and syncs metadata
summaries to the user's Personal Data Store (PDS).

### Edge Computing Model

| Local Only              | Syncs to PDS                              |
|:------------------------|:------------------------------------------|
| BAM/CRAM/VCF files      | Haplogroup assignments + private variants |
| Genotype chip files     | Coverage/quality metrics                  |
| Full analysis artifacts | Ancestry percentages                      |
| IBD segment data        | STR marker values                         |
|                         | File metadata (not content)               |

### Key Concepts

- **Biosample ≈ Subject**: In Navigator, each Subject maps 1:1 to a Biosample record
- **Specimen Donor**: AppView-level entity linking multiple biosamples from the same individual. For Edge App users,
  `donorIdentifier` field serves as the link
- **Reconciliation**: Lives at Specimen Donor level (not biosample) to handle multi-kit scenarios

---

## Record Types Used by Navigator

### Core Records (Implemented)

| Record       | NSID                                    | Navigator Model | Status     |
|:-------------|:----------------------------------------|:----------------|:-----------|
| Workspace    | `com.decodingus.atmosphere.workspace`   | `Workspace`     | ✅ Complete |
| Biosample    | `com.decodingus.atmosphere.biosample`   | `Biosample`     | ✅ Complete |
| Sequence Run | `com.decodingus.atmosphere.sequencerun` | `SequenceRun`   | ✅ Complete |
| Alignment    | `com.decodingus.atmosphere.alignment`   | `Alignment`     | ✅ Complete |
| Project      | `com.decodingus.atmosphere.project`     | `Project`       | ✅ Complete |

### Analysis Records (In Development)

| Record               | NSID                                            | Navigator Model  | Status     |
|:---------------------|:------------------------------------------------|:-----------------|:-----------|
| Genotype             | `com.decodingus.atmosphere.genotype`            | `ChipProfile`    | 🚧 In Dev  |
| Population Breakdown | `com.decodingus.atmosphere.populationBreakdown` | `AncestryResult` | 🚧 In Dev  |
| STR Profile          | `com.decodingus.atmosphere.strProfile`          | `StrProfile`     | ✅ Complete |

### Future Records (Not Yet Implemented)

| Record                    | NSID                                                 | Purpose             |
|:--------------------------|:-----------------------------------------------------|:--------------------|
| Haplogroup Reconciliation | `com.decodingus.atmosphere.haplogroupReconciliation` | Multi-run consensus |
| Match Consent             | `com.decodingus.atmosphere.matchConsent`             | IBD matching opt-in |
| Match List                | `com.decodingus.atmosphere.matchList`                | IBD match results   |

---

## Navigator Model Mappings

### Workspace → `com.decodingus.atmosphere.workspace`

```scala
case class Workspace(
                      meta: RecordMeta,
                      main: WorkspaceContent // Contains sampleRefs, projectRefs
                    )
```

### Biosample → `com.decodingus.atmosphere.biosample`

```scala
case class Biosample(
                      atUri: Option[String],
                      meta: RecordMeta,
                      sampleAccession: String, // Unique ID
                      donorIdentifier: String, // Links to Specimen Donor
                      description: Option[String],
                      centerName: String,
                      sex: Option[String],
                      haplogroups: Option[HaplogroupAssignments],
                      sequenceRunRefs: List[String],
                      genotypeRefs: List[String],
                      populationBreakdownRef: Option[String],
                      strProfileRef: Option[String]
                    )
```

### SequenceRun → `com.decodingus.atmosphere.sequencerun`

```scala
case class SequenceRun(
                        atUri: Option[String],
                        meta: RecordMeta,
                        biosampleRef: String,
                        platformName: String, // ILLUMINA, PACBIO, NANOPORE, etc.
                        instrumentModel: Option[String],
                        testType: String, // WGS, WES, BIG_Y_700, etc.
                        libraryLayout: Option[String],
                        readLength: Option[Int],
                        files: List[FileInfo],
                        alignmentRefs: List[String]
                      )
```

### Alignment → `com.decodingus.atmosphere.alignment`

```scala
case class Alignment(
                      atUri: Option[String],
                      meta: RecordMeta,
                      sequenceRunRef: String,
                      referenceBuild: String, // GRCh38, GRCh37, T2T-CHM13
                      aligner: Option[String],
                      files: List[FileInfo],
                      metrics: Option[AlignmentMetrics]
                    )
```

### AlignmentMetrics (embedded in Alignment)

```scala
case class AlignmentMetrics(
                             genomeTerritory: Option[Long],
                             meanCoverage: Option[Double],
                             medianCoverage: Option[Double],
                             pct10x: Option[Double],
                             pct20x: Option[Double],
                             pct30x: Option[Double],
                             callableBases: Option[Long],
                             contigs: List[ContigMetrics]
                           )
```

### HaplogroupResult (embedded in Biosample.haplogroups)

```scala
case class HaplogroupResult(
                             haplogroupName: String,
                             score: Double,
                             matchingSnps: Int,
                             ancestralMatches: Int,
                             treeDepth: Int,
                             lineagePath: Option[List[String]],
                             privateVariants: Option[List[PrivateVariant]]
                           )
```

---

## Test Type Taxonomy

Navigator supports these test type codes for sequence runs and genotypes:

### Sequencing Test Types

| Code        | Display Name            | Y-DNA | mtDNA |
|:------------|:------------------------|:-----:|:-----:|
| `WGS`       | Whole Genome Sequencing |   ✓   |   ✓   |
| `WGS_LP`    | Low-Pass WGS            |   ✓   |   ✓   |
| `WES`       | Whole Exome Sequencing  |   ✗   |   ✗   |
| `HIFI_WGS`  | PacBio HiFi WGS         |   ✓   |   ✓   |
| `ONT_WGS`   | Nanopore WGS            |   ✓   |   ✓   |
| `BIG_Y_700` | FTDNA Big Y-700         |   ✓   |   ✗   |
| `BIG_Y_500` | FTDNA Big Y-500         |   ✓   |   ✗   |
| `Y_ELITE`   | YSEQ Y Elite            |   ✓   |   ✗   |
| `MT_FULL`   | mtDNA Full Sequence     |   ✗   |   ✓   |

### Genotype (Chip) Test Types

| Code                | Display Name   | Vendor        |
|:--------------------|:---------------|:--------------|
| `ARRAY_23ANDME_V5`  | 23andMe v5     | 23andMe       |
| `ARRAY_ANCESTRY_V2` | AncestryDNA v2 | Ancestry      |
| `ARRAY_FTDNA_V3`    | FTDNA v3       | FamilyTreeDNA |
| `ARRAY_MYHERITAGE`  | MyHeritage     | MyHeritage    |
| `ARRAY_LIVINGDNA`   | LivingDNA      | LivingDNA     |

---

## Local Storage Structure

Navigator stores analysis artifacts locally:

```
~/.decodingus/
├── config/
│   ├── workspace.json              # Workspace state
│   └── reference_config.json       # Reference genome paths
└── cache/
    ├── references/                 # Downloaded reference genomes
    ├── trees/                      # Haplogroup tree data
    │   ├── ftdna-ytree.json
    │   └── *-GRCh38-sites.vcf
    ├── {sha256}.json               # Analysis cache by file hash
    └── subjects/{sampleAccession}/
        └── runs/{runId}/
            └── alignments/{alignmentId}/
                ├── wgs_metrics.txt
                ├── callable_loci/
                │   └── chr*.callable.svg
                └── haplogroup/
                    ├── ydna_tree_sites.vcf
                    ├── ydna_private_variants.vcf
                    └── ydna_report.txt
```

---

## PDS Sync Strategy

### Current Implementation (Feature-Flagged)

1. **Workspace-level sync**: Entire workspace JSON synced as single record
2. **On-demand**: User explicitly triggers sync
3. **Read-only from AppView**: AppView indexes but doesn't modify PDS records

### Future Implementation

1. **Granular record sync**: Individual record CRUD operations
2. **Automatic sync**: Background sync on record changes
3. **Conflict resolution**: Version-based merge strategy

---

## Implementation Status

See [Edge_Client_Implementation_Status.md](Edge_Client_Implementation_Status.md) for detailed tracking.

**Overall: ~55% complete**

- ✅ Core records (workspace, biosample, sequenceRun, alignment, project)
- ✅ Haplogroup analysis (Y-DNA, mtDNA, private variants)
- ✅ STR profile support
- 🚧 Genotype (chip) parsing
- 🚧 Ancestry estimation
- ⬜ IBD matching
- ⬜ Multi-run reconciliation

---

## Changelog

| Version | Date       | Changes                                                          |
|:--------|:-----------|:-----------------------------------------------------------------|
| 1.0     | 2025-12-05 | Initial design                                                   |
| 1.5     | 2025-12-08 | Added reconciliation, ancestry, STR support                      |
| 1.7     | 2025-12-08 | Multi-test-type taxonomy                                         |
| 2.0     | 2025-12-09 | Streamlined to Edge App reference; full schema in global lexicon |
