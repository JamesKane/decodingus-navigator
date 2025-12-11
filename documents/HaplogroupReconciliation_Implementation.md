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
| `ui/components/WorkbenchView.scala` | Added reconciliation status indicator to subject detail |
| `ui/components/ReconciliationDetailDialog.scala` | **NEW** - Dialog showing all run calls and reconciliation status |
| `documents/Atmosphere_Lexicon.md` | Updated documentation |

## Current State

### What Works Now (Phase 1 & 2 Complete)

#### Phase 1: Reconciliation Records
- `WorkspaceOperations` has helper methods for reconciliation management:
  - `getOrCreateReconciliation()` - Creates/retrieves Y-DNA or mtDNA reconciliation
  - `addHaplogroupCall()` - Adds a run call and auto-selects consensus
  - `removeHaplogroupCall()` - Removes a call when a run is deleted
- `AnalysisCoordinator.runHaplogroupInternal()` creates `RunHaplogroupCall` and uses reconciliation
- `WorkbenchViewModel` chip haplogroup path uses reconciliation
- Consensus selection uses quality tier ranking (WGS > Big Y > Chip)

#### Phase 2: UI Integration
- Subject detail view shows haplogroup results with reconciliation status indicator
- Traffic light colors: Green (compatible), Orange (minor divergence), Red (major divergence), Purple (incompatible)
- Clickable button shows run count and opens detail dialog
- `ReconciliationDetailDialog` displays:
  - Y-DNA and mtDNA reconciliation panels
  - Consensus haplogroup with confidence
  - Status indicator with tooltip
  - Table of all individual run calls with technology, haplogroup, SNPs, confidence, tree provider

### What's NOT Yet Implemented

#### 1. Cleanup on Run/Profile Deletion
When a SequenceRun or ChipProfile is deleted:
- Remove its call from `HaplogroupReconciliation.runCalls`
- Recalculate consensus
- Update `Biosample.haplogroups`

#### 2. Advanced Reconciliation Logic
- Branch compatibility scoring (LCA analysis between different haplogroups)
- Automatic detection of COMPATIBLE vs MINOR_DIVERGENCE vs MAJOR_DIVERGENCE
- Currently always sets COMPATIBLE; needs tree structure comparison

## Next Steps

### Phase 3: Advanced Features
1. Cleanup on run/profile deletion
2. SNP-level conflict detection
3. Branch compatibility scoring (LCA analysis)
4. Identity verification metrics
5. Manual override capability

## Related Documents

- [MultiRunReconciliation.md](design/MultiRunReconciliation.md) - Original design document
- [Atmosphere_Lexicon.md](Atmosphere_Lexicon.md) - Edge App schema reference
- Global schema: `/decodingus/documents/atmosphere/09-Reconciliation-Records.md`
