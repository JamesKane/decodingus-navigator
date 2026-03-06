# Sequence Run Fingerprinting - Design Document

## Status: Implemented (Phase 3)

## Problem Statement

A single sequencing run (same reads from the same library prep) can be aligned to multiple reference genomes:
- GRCh38 (current standard)
- CHM13v2.0 / T2T (telomere-to-telomere)
- GRCh37/hg19 (legacy)

When a user loads both the GRCh38 and CHM13 alignments of the same sample, we need to:
1. Recognize they represent the same underlying sequencing run
2. Group them under a single `SequenceRun` with multiple `Alignment` records
3. Avoid creating duplicate `SequenceRun` entries

### Current Limitation

The current duplicate detection uses SHA-256 checksum of the BAM/CRAM file. This fails for multi-reference scenarios because:
- Different reference = different aligned positions = different file content = different checksum
- The underlying reads are identical, but the aligned coordinates differ

## Fingerprinting Strategy

### Tier 1: Read Group (@RG) Header Matching

BAM/CRAM files contain `@RG` (Read Group) headers with metadata about the sequencing run.

**GATK Required Fields** (reliable in practice):

| Field | Description | GATK Required | Stability |
|-------|-------------|---------------|-----------|
| `ID` | Read group identifier | Yes | May change during re-alignment |
| `SM` | Sample name | Yes | **Stable** - set at sequencing |
| `LB` | Library identifier | Yes | **Stable** - identifies library prep |
| `PL` | Platform (ILLUMINA, PACBIO, etc.) | Yes | **Stable** |
| `PU` | Platform unit (flowcell.lane.barcode) | No* | **Stable** - most unique |

*PU is not strictly required by GATK but is used by BQSR when present.

**Other fields** (CN, DT, PM, etc.) are inconsistent in the wild and cannot be relied upon.

**Fingerprint Priority:**
1. `PU` alone (if present) - most unique identifier
2. `LB + SM + PL` combination - GATK-required fields
3. Read statistics fallback - when headers are incomplete

**Confidence Levels:**
- HIGH: PU present, or exact fingerprint match
- MEDIUM: LB + SM match (no PU)
- LOW: Stats-only match (needs user confirmation)

### Tier 2: Read Statistics Fingerprinting

When @RG headers are incomplete or missing, use read-level statistics:

| Statistic | Tolerance | Notes |
|-----------|-----------|-------|
| Total read count | Exact | Must match exactly |
| Read length distribution | Exact | Histogram of read lengths |
| Mean insert size | ±1% | Slight variation from sampling |
| Paired read ratio | ±0.1% | Should be identical |

**Composite Fingerprint:**
```scala
case class ReadStatsFingerprint(
  totalReads: Long,
  readLengthHistogramHash: String,  // Hash of sorted histogram
  meanInsertSize: Double,           // Rounded to nearest 0.1
  pairedReadRatio: Double           // Rounded to nearest 0.001
)
```

**Confidence Level:** HIGH if all match, requires user confirmation otherwise

### Tier 3: Read Name Sampling

For highest confidence, sample a subset of read names (QNAMEs):

```scala
// Sample first N read names from first 1000 reads
val sampleReadNames: Set[String] = reads.take(1000).map(_.getReadName).toSet.take(100)
```

If 95%+ of sampled read names match between two files, they're from the same run.

**Confidence Level:** DEFINITIVE

## Data Model Changes

### SequenceRun - Add Fingerprint Fields

```scala
case class SequenceRun(
  // ... existing fields ...

  // Fingerprint fields for matching across references
  libraryId: Option[String] = None,        // @RG LB - library prep identifier
  platformUnit: Option[String] = None,     // @RG PU - flowcell.lane.barcode
  runFingerprint: Option[String] = None,   // Computed composite fingerprint hash

  // Already have these:
  // sampleName: Option[String]  // @RG SM
  // instrumentId: Option[String] // from QNAME parsing
  // flowcellId: Option[String]
  // totalReads: Option[Long]
)
```

### RunFingerprint Model

```scala
/**
 * Fingerprint data for identifying identical sequencing runs
 * across different reference alignments.
 */
case class RunFingerprint(
  // Tier 1: @RG header fields
  sampleName: Option[String],
  libraryId: Option[String],
  platformUnit: Option[String],
  platform: Option[String],
  sequencingCenter: Option[String],
  runDate: Option[String],

  // Tier 2: Read statistics
  totalReads: Long,
  readLengthHistogramHash: String,
  pairedReadRatio: Double,

  // Tier 3: Read name sample (optional, for verification)
  sampledReadNamesHash: Option[String]
) {
  /**
   * Compute fingerprint hash for comparison.
   * Uses available fields in priority order.
   */
  def computeHash: String = {
    // Primary: PU is most unique
    platformUnit.map(pu => sha256(s"PU:$pu")).getOrElse {
      // Secondary: LB + SM combination
      (libraryId, sampleName) match {
        case (Some(lb), Some(sm)) => sha256(s"LB:$lb:SM:$sm")
        // Tertiary: Read stats
        case _ => sha256(s"READS:$totalReads:$readLengthHistogramHash:$pairedReadRatio")
      }
    }
  }

  def confidenceLevel: FingerprintConfidence = {
    if (platformUnit.isDefined) FingerprintConfidence.High
    else if (libraryId.isDefined && sampleName.isDefined) FingerprintConfidence.Medium
    else FingerprintConfidence.Low
  }
}

enum FingerprintConfidence:
  case High    // PU or PU+LB+SM - very reliable
  case Medium  // LB+SM only - usually reliable
  case Low     // Stats only - may need user confirmation
```

## Matching Algorithm

```scala
def findMatchingSequenceRun(
  newFingerprint: RunFingerprint,
  existingRuns: List[SequenceRun],
  biosampleRef: String
): Option[(SequenceRun, MatchConfidence)] = {

  // Only consider runs for same biosample
  val candidateRuns = existingRuns.filter(_.biosampleRef == biosampleRef)

  // Tier 1: Exact PU match
  newFingerprint.platformUnit.flatMap { pu =>
    candidateRuns.find(_.platformUnit.contains(pu))
      .map(run => (run, MatchConfidence.Definitive))
  }.orElse {
    // Tier 2: LB + SM match
    (newFingerprint.libraryId, newFingerprint.sampleName) match {
      case (Some(lb), Some(sm)) =>
        candidateRuns.find(r => r.libraryId.contains(lb) && r.sampleName.contains(sm))
          .map(run => (run, MatchConfidence.High))
      case _ => None
    }
  }.orElse {
    // Tier 3: Read stats match
    candidateRuns.find { run =>
      run.totalReads.contains(newFingerprint.totalReads) &&
      // Could add more stats comparison here
      true
    }.map(run => (run, MatchConfidence.NeedsConfirmation))
  }
}
```

## User Experience

### Automatic Grouping (High Confidence)

When fingerprint match confidence is HIGH or DEFINITIVE:
1. Automatically add new Alignment to existing SequenceRun
2. Show notification: "Detected CHM13 alignment of existing GRCh38 sample - grouped automatically"

### User Confirmation (Low Confidence)

When confidence is LOW or stats-only:
1. Show dialog: "This file appears to be a different alignment of an existing sample. Group them together?"
2. Options: "Group Together" / "Keep Separate" / "Compare Details"
3. "Compare Details" shows side-by-side comparison of metadata

### Manual Grouping

User can always:
1. Right-click on an Alignment
2. Select "Link to different Sequence Run..."
3. Choose from list of runs for same biosample

## Implementation Phases

### Phase 1: Extract Fingerprint Data ✅
- [x] Add @RG parsing to LibraryStatsProcessor (LB, PU)
- [x] Compute read length histogram hash
- [x] Add fingerprint fields to SequenceRun model
- [x] Add fingerprint fields to LibraryStats

### Phase 2: Fingerprint Matching ✅
- [x] Implement findMatchingSequenceRun algorithm
- [x] Update addFileAndAnalyze to use fingerprint matching
- [x] Auto-group HIGH/MEDIUM confidence matches

### Phase 3: User Experience ✅
- [x] Add confirmation dialog for low-confidence matches (FingerprintMatchDialog)
- [ ] Add "Link to different Sequence Run" context menu (deferred)
- [x] Update SequenceDataTable to show multiple alignments (Ref column shows all references)
- [x] Update SequenceDataTable to show multiple files (File column shows count)

### Phase 4: Verification
- [ ] Add read name sampling for definitive verification
- [ ] Add "Verify Match" action for uncertain groupings

## Edge Cases

1. **Split FASTQ processing**: Same sample sequenced on multiple lanes, processed separately
   - Different PU values but same LB
   - Solution: Group by LB when PUs differ but other metadata matches

2. **Re-sequencing**: Same sample sequenced multiple times
   - Same SM but different PU, DT, maybe different LB
   - Solution: Keep separate - these are genuinely different runs

3. **Merged BAMs**: Multiple runs merged into single file
   - Multiple @RG entries
   - Solution: Extract all RGs, create/link to multiple SequenceRuns

4. **Missing headers**: BAM with no @RG or minimal headers
   - Solution: Fall back to read stats, require user confirmation

## References

- [SAM Format Specification - Read Groups](https://samtools.github.io/hts-specs/SAMv1.pdf)
- [GATK Read Groups](https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups)
