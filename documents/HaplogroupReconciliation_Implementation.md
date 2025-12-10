# Haplogroup Reconciliation Implementation

## Overview

This document describes the multi-run haplogroup reconciliation system implemented to prevent lower-quality results (e.g., chip data) from overwriting higher-quality results (e.g., WGS data), and to properly track results across different tree providers.

## Problem Statement

Previously, running haplogroup analysis would unconditionally overwrite existing results:
- Chip haplogroup analysis would overwrite more accurate WGS results
- Switching tree providers (FTDNA ↔ DecodingUs) would lose previous results
- No tracking of which run/source produced which result

## Implemented Solution

### Schema Alignment

The implementation aligns with the global Atmosphere Lexicon schema at `/decodingus/documents/atmosphere/`:

| Component | Location | Purpose |
|-----------|----------|---------|
| `HaplogroupAssignments` | `Biosample.haplogroups` | Consensus/best result (syncs to PDS) |
| `HaplogroupReconciliation` | `WorkspaceContent.haplogroupReconciliations` | Multi-run tracking (local) |

### Model Changes

#### HaplogroupAssignments (Simplified)

```scala
case class HaplogroupAssignments(
  yDna: Option[HaplogroupResult] = None,
  mtDna: Option[HaplogroupResult] = None
)
```

Stores only the consensus/best result for each DNA type. This is what syncs to PDS.

#### HaplogroupResult (Extended)

```scala
case class HaplogroupResult(
  haplogroupName: String,
  score: Double,
  matchingSnps: Option[Int],
  mismatchingSnps: Option[Int],
  ancestralMatches: Option[Int],
  treeDepth: Option[Int],
  lineagePath: Option[List[String]],
  privateVariants: Option[PrivateVariantData],
  // Provenance tracking
  source: Option[String],       // "wgs", "bigy", "chip"
  sourceRef: Option[String],    // AT URI of SequenceRun or ChipProfile
  treeProvider: Option[String], // "ftdna", "decodingus"
  treeVersion: Option[String],
  analyzedAt: Option[Instant]
)
```

Added provenance fields to track where the result came from.

#### HaplogroupReconciliation (New)

```scala
case class HaplogroupReconciliation(
  atUri: Option[String],
  meta: RecordMeta,
  biosampleRef: String,
  dnaType: DnaType,              // Y_DNA or MT_DNA
  status: ReconciliationStatus,
  runCalls: List[RunHaplogroupCall],
  snpConflicts: List[SnpConflict],
  lastReconciliationAt: Option[Instant]
)
```

Tracks all haplogroup calls from different analyses and maintains reconciliation status.

### Supporting Types

| Type | Values | Purpose |
|------|--------|---------|
| `DnaType` | Y_DNA, MT_DNA | DNA type for reconciliation |
| `CompatibilityLevel` | COMPATIBLE, MINOR_DIVERGENCE, MAJOR_DIVERGENCE, INCOMPATIBLE | Cross-run agreement |
| `HaplogroupTechnology` | WGS, WES, BIG_Y, SNP_ARRAY, AMPLICON, STR_PANEL | Source technology |
| `CallMethod` | SNP_PHYLOGENETIC, STR_PREDICTION, VENDOR_REPORTED | How haplogroup was determined |
| `ConflictResolution` | ACCEPT_MAJORITY, ACCEPT_HIGHER_QUALITY, ACCEPT_HIGHER_COVERAGE, UNRESOLVED, HETEROPLASMY | Conflict resolution method |

### Reconciliation Logic

Quality tier ranking (higher = better):
1. **WGS** (tier 3) - Highest quality, full coverage
2. **Big Y / Y Elite** (tier 2) - Targeted Y-DNA sequencing
3. **Chip / SNP Array** (tier 1) - Limited SNP coverage

When selecting the consensus result:
1. Higher quality tier wins
2. At same tier, deeper tree depth wins
3. At same depth, more recent analysis wins

## Files Modified

| File | Changes |
|------|---------|
| `workspace/model/HaplogroupAssignments.scala` | Simplified to single yDna/mtDna fields |
| `workspace/model/HaplogroupResult.scala` | Added provenance fields and qualityTier method |
| `workspace/model/HaplogroupReconciliation.scala` | **NEW** - Full reconciliation model |
| `workspace/model/Workspace.scala` | Added `haplogroupReconciliations` to WorkspaceContent |
| `workspace/WorkspaceService.scala` | Added Circe codecs for new types |
| `workspace/services/AnalysisCoordinator.scala` | Updated to use new model |
| `workspace/WorkbenchViewModel.scala` | Updated WGS and chip haplogroup paths |
| `pds/PdsClient.scala` | Added codecs for PDS sync |
| `documents/Atmosphere_Lexicon.md` | Updated documentation |

## Current State

### What Works Now
- Haplogroup analysis stores results with provenance tracking
- Results include source type, source reference, tree provider, and timestamp
- Model structure matches global Atmosphere schema

### What's NOT Yet Implemented

#### 1. Reconciliation Record Creation
The analysis paths update `Biosample.haplogroups` but don't yet create/update `HaplogroupReconciliation` records.

**Required changes:**
- `AnalysisCoordinator.runHaplogroupInternal()` - Create/update Y-DNA or mtDNA reconciliation
- `WorkbenchViewModel.runChipHaplogroupAnalysis()` - Create/update reconciliation for chip results
- `WorkbenchViewModel.runHaplogroupAnalysis()` - Create/update reconciliation for WGS results

#### 2. Consensus Selection from Reconciliation
Currently, analysis directly updates `Biosample.haplogroups`. Should instead:
1. Add run call to `HaplogroupReconciliation`
2. Recalculate consensus using `recalculate()` method
3. Update `Biosample.haplogroups` with the consensus result

#### 3. UI for Reconciliation Status
- Show reconciliation status indicator on subject card
- Traffic light: Green (compatible), Yellow (minor divergence), Red (major/incompatible)
- Detailed view showing all run calls and any conflicts

#### 4. Cleanup on Run/Profile Deletion
When a SequenceRun or ChipProfile is deleted:
- Remove its call from `HaplogroupReconciliation.runCalls`
- Recalculate consensus
- Update `Biosample.haplogroups`

## Next Steps

### Phase 1: Wire Up Reconciliation Records
1. Update `AnalysisCoordinator` to create/update reconciliation records
2. Update `WorkbenchViewModel` chip path similarly
3. Move consensus selection logic to use reconciliation

### Phase 2: UI Integration
1. Add reconciliation status to subject card
2. Create reconciliation detail view
3. Add ability to view all run calls

### Phase 3: Advanced Features
1. SNP-level conflict detection
2. Branch compatibility scoring (LCA analysis)
3. Identity verification metrics
4. Manual override capability

## Related Documents

- [MultiRunReconciliation.md](design/MultiRunReconciliation.md) - Original design document
- [Atmosphere_Lexicon.md](Atmosphere_Lexicon.md) - Edge App schema reference
- Global schema: `/decodingus/documents/atmosphere/09-Reconciliation-Records.md`
