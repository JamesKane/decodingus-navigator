# Unified Quality Metrics Walker

## Status: Active Development

## Motivation

**CollectMultipleMetrics has a hard dependency on R** for generating insert size histograms. There is no way to disable chart generation, making it unsuitable for a self-contained desktop application.

This elevates the unified walker from a performance optimization to a **required replacement**.

## Problem Statement

Currently, gathering comprehensive quality metrics for a BAM/CRAM file requires **three separate passes** through the file:

1. **CollectWgsMetrics** - Coverage depth, genome coverage percentages
2. **CollectMultipleMetrics** - Alignment summary (read counts, alignment rates) + Insert size distribution (**requires R**)
3. **CallableLoci** - Per-position callable base analysis with contig-level summaries

Each of these tools reads the entire BAM/CRAM file sequentially. For a typical 30x WGS file (~100GB), this means reading ~300GB of data when a single pass could collect everything.

### Additional Constraint

CollectMultipleMetrics/CollectInsertSizeMetrics requires R to be installed for histogram generation. This is unacceptable for a standalone desktop application - we cannot require users to install R just to get insert size metrics.

### Current Tool Matrix

| Metric Category | Current Tool | Pass Required |
|-----------------|--------------|---------------|
| Total reads, PF reads | CollectAlignmentSummaryMetrics | Pass 1 |
| Aligned reads, paired reads | CollectAlignmentSummaryMetrics | Pass 1 |
| Mean read length | CollectAlignmentSummaryMetrics | Pass 1 |
| Insert size distribution | CollectInsertSizeMetrics | Pass 1 |
| Mean/median coverage | CollectWgsMetrics | Pass 2 |
| Coverage histogram | CollectWgsMetrics | Pass 2 |
| PCT_1X, PCT_10X, PCT_20X, etc. | CollectWgsMetrics | Pass 2 |
| Callable bases per contig | CallableLoci | Pass 3 |
| Callable regions BED | CallableLoci | Pass 3 |

## Proposed Solution

Implement a custom HTSJDK-based Walker that collects all metrics in a single pass through the BAM/CRAM file.

### Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                  UnifiedMetricsWalker                       │
├─────────────────────────────────────────────────────────────┤
│  Input: BAM/CRAM + Reference                                │
│                                                             │
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐         │
│  │ ReadCounter │  │ Coverage    │  │ Callable    │         │
│  │ Collector   │  │ Accumulator │  │ Loci        │         │
│  │             │  │             │  │ Tracker     │         │
│  │ - total     │  │ - depth[]   │  │ - state per │         │
│  │ - aligned   │  │ - histogram │  │   position  │         │
│  │ - paired    │  │ - by contig │  │ - BED out   │         │
│  │ - insert sz │  │             │  │             │         │
│  └─────────────┘  └─────────────┘  └─────────────┘         │
│         │                │                │                 │
│         └────────────────┼────────────────┘                 │
│                          ▼                                  │
│               UnifiedMetricsResult                          │
└─────────────────────────────────────────────────────────────┘
```

### Implementation Options

#### Option A: Pure HTSJDK Walker

Use `SamReader` directly with a custom traversal loop:

```scala
class UnifiedMetricsWalker {
  def process(bamPath: String, refPath: String): UnifiedMetricsResult = {
    val reader = SamReaderFactory.makeDefault()
      .referenceSequence(new File(refPath))
      .open(new File(bamPath))

    val iterator = reader.iterator()
    while (iterator.hasNext) {
      val read = iterator.next()
      readCounter.accept(read)
      coverageAccumulator.accept(read)
      callableLociTracker.accept(read)
    }

    buildResult()
  }
}
```

**Pros:**
- No GATK dependency for this tool
- Full control over memory management
- Can stream results for very large files

**Cons:**
- Must implement pileup logic ourselves for coverage
- More code to maintain

#### Option B: GATK Walker Extension

Extend GATK's `ReadWalker` or `LocusWalker`:

```scala
class UnifiedMetricsWalker extends ReadWalker[UnifiedMetricsResult] {
  override def apply(read: GATKRead, ref: ReferenceContext, fc: FeatureContext): Unit = {
    // Collect read-level metrics
  }

  override def onTraversalDone(result: UnifiedMetricsResult): UnifiedMetricsResult = {
    // Finalize and return
  }
}
```

**Pros:**
- Reuses GATK infrastructure (indexing, parallelization, pileup)
- Can use GATK's Spark support for distributed execution
- Well-tested read traversal logic

**Cons:**
- Tied to GATK version
- Heavier dependency
- GATK Walker API has restrictions

#### Option C: Hybrid - HTSJDK with Pileup Library

Use HTSJDK for reading but leverage a lightweight pileup implementation:

```scala
// Use HTSJDK's SamLocusIterator for pileup
val locusIterator = new SamLocusIterator(reader)
locusIterator.setEmitUncoveredLoci(false)

for (pileup <- locusIterator.asScala) {
  val depth = pileup.getRecordAndOffsets.size()
  coverageAccumulator.addDepth(pileup.getSequenceIndex, pileup.getPosition, depth)
  callableLociTracker.update(pileup)
}
```

**Pros:**
- Simpler than full GATK dependency
- `SamLocusIterator` handles pileup complexity
- Good balance of control and reuse

**Cons:**
- `SamLocusIterator` may not expose all needed info
- Still single-threaded without extra work

### Recommended Approach: Option C (Hybrid)

The hybrid approach provides the best balance:
1. Use HTSJDK's `SamLocusIterator` for position-by-position pileup
2. Track read-level stats on first encounter of each read
3. Accumulate coverage and callable state per position
4. Output results compatible with existing GATK tool formats

## Data Structures

### Memory-Efficient Coverage Tracking

For 3.1Gb genome at single-base resolution:
- `byte[]` array: 3.1GB (impractical)
- Sparse map: Variable, but high overhead
- **Windowed/binned**: Fixed reasonable size

```scala
case class CoverageAccumulator(
  windowSize: Int = 1000,  // 1kb windows
  histogram: Array[Long] = new Array[Long](256),  // depth 0-255
  contigWindows: Map[String, Array[Short]]  // per-contig windowed coverage
)
```

For typical 30x WGS:
- ~3.1M windows at 1kb resolution
- 2 bytes per window = ~6MB per contig
- ~150MB total for main assembly

### Callable Loci State Machine

Based on GATK's `CallableLoci` which extends `LocusWalker`:

```scala
/**
 * Callable state enum matching GATK CallableLoci states.
 * Evaluation order matters - checked hierarchically.
 */
enum CallableState:
  case RefN              // Reference base is N (non-callable by definition)
  case NoCoverage        // Zero reads at locus
  case PoorMappingQuality // High fraction of reads with low MAPQ
  case LowCoverage       // Insufficient QC+ reads
  case ExcessiveCoverage // Exceeds max depth threshold
  case Callable          // Meets all criteria

/**
 * Parameters aligned with GATK CallableLoci defaults.
 */
class CallableLociTracker(
  minDepth: Int = 4,                    // GATK default: 4
  maxDepth: Option[Int] = None,         // GATK default: unlimited
  minMappingQuality: Int = 10,          // GATK default: 10
  minBaseQuality: Int = 20,             // GATK default: 20
  maxLowMapQ: Int = 1,                  // MAPQ threshold for "low" reads
  maxFractionLowMapQ: Double = 0.1      // GATK default: 0.1 (10%)
) {
  def stateAt(refBase: Byte, pileup: SamLocusIterator.LocusInfo): CallableState = {
    // 1. Check reference
    if (refBase == 'N' || refBase == 'n') return CallableState.RefN

    val allReads = pileup.getRecordAndOffsets
    val rawDepth = allReads.size()

    // 2. No coverage
    if (rawDepth == 0) return CallableState.NoCoverage

    // 3. Count QC-passing and low-MAPQ reads
    var qcPassCount = 0
    var lowMapQCount = 0
    val iter = allReads.iterator()
    while (iter.hasNext) {
      val rec = iter.next()
      val read = rec.getRead
      val mapQ = read.getMappingQuality
      val baseQ = rec.getBaseQuality

      if (mapQ >= minMappingQuality && baseQ >= minBaseQuality) {
        qcPassCount += 1
      }
      if (mapQ <= maxLowMapQ) {
        lowMapQCount += 1
      }
    }

    // 4. Poor mapping quality (too many low-MAPQ reads)
    val lowMapQFraction = lowMapQCount.toDouble / rawDepth
    if (lowMapQFraction > maxFractionLowMapQ) return CallableState.PoorMappingQuality

    // 5. Low coverage (after QC filtering)
    if (qcPassCount < minDepth) return CallableState.LowCoverage

    // 6. Excessive coverage
    if (maxDepth.exists(qcPassCount > _)) return CallableState.ExcessiveCoverage

    // 7. Callable
    CallableState.Callable
  }
}
```

**Key insight from GATK source**: The state evaluation is hierarchical - a locus is classified by the *first* failing condition, not a combination. This matches how GATK reports summary counts.

## Output Format

### UnifiedMetricsResult

```scala
case class UnifiedMetricsResult(
  // Read-level (from CollectAlignmentSummaryMetrics)
  totalReads: Long,
  pfReads: Long,
  pfReadsAligned: Long,
  pctPfReadsAligned: Double,
  readsPaired: Long,
  pctReadsPaired: Double,
  pctProperPairs: Double,
  meanReadLength: Double,

  // Insert size (from CollectInsertSizeMetrics)
  medianInsertSize: Double,
  meanInsertSize: Double,
  stdInsertSize: Double,
  pairOrientation: String,

  // Coverage (from CollectWgsMetrics) - genome-wide
  meanCoverage: Double,
  medianCoverage: Double,
  sdCoverage: Double,
  pct1x: Double,
  pct5x: Double,
  pct10x: Double,
  pct20x: Double,
  pct30x: Double,
  coverageHistogram: Array[Long],

  // Coverage - per chromosome (for coverage histograms/visualizations)
  contigCoverage: Map[String, ContigCoverageMetrics],

  // Callable bases (from CallableLoci) - genome-wide
  callableBases: Long,
  callableRegionsBed: Option[Path],

  // Callable bases - per chromosome (for callable loci visualizations)
  contigCallableLoci: Map[String, ContigCallableMetrics]
)

/**
 * Per-chromosome coverage metrics for coverage histogram visualizations.
 */
case class ContigCoverageMetrics(
  contig: String,
  length: Long,
  meanCoverage: Double,
  medianCoverage: Double,
  pct1x: Double,
  pct10x: Double,
  pct20x: Double,
  pct30x: Double,
  coverageHistogram: Array[Long]  // depth 0-255+ binned counts
)

/**
 * Per-chromosome callable loci metrics matching current CallableLociResult.
 */
case class ContigCallableMetrics(
  contig: String,
  length: Long,
  callableBases: Long,
  noCoverageBases: Long,
  lowCoverageBases: Long,
  excessiveCoverageBases: Long,
  poorMappingQualityBases: Long,
  refNBases: Long,
  pctCallable: Double
)
```

## Performance Expectations

| Scenario | Current (3 passes) | Unified (1 pass) | Speedup |
|----------|-------------------|------------------|---------|
| 30x WGS (~100GB) | ~90 min | ~30 min | 3x |
| 60x WGS (~200GB) | ~180 min | ~60 min | 3x |
| Exome (~10GB) | ~10 min | ~4 min | 2.5x |

*Note: I/O bound, so speedup is roughly linear with pass reduction*

## Implementation Plan

### Phase 1: Core Walker (Priority: Immediate)
- [ ] Implement `UnifiedMetricsWalker` with HTSJDK `SamReader`
- [ ] Read counter collector (total, aligned, paired, proper pairs)
- [ ] Insert size accumulator with histogram (replaces R-dependent CollectInsertSizeMetrics)
- [ ] Mean read length calculation
- [ ] Basic coverage depth tracking

**This phase eliminates the R dependency by replacing CollectMultipleMetrics.**

### Phase 2: Coverage & Callable
- [ ] Windowed coverage accumulator
- [ ] Callable loci state machine (per GATK CallableLoci logic)
- [ ] BED file output for callable regions
- [ ] Coverage percentile calculations (PCT_1X, PCT_10X, etc.)

### Phase 3: Integration
- [ ] Create `UnifiedMetricsProcessor` following existing processor pattern
- [ ] Replace `MultipleMetricsProcessor` (remove R dependency)
- [ ] Add to WorkbenchViewModel as replacement
- [ ] Update SequenceRun/Alignment models if needed

### Phase 4: Optimization
- [ ] Multi-threaded contig processing
- [ ] Memory profiling and tuning
- [ ] Progress reporting with accurate ETA
- [ ] Optional Spark-based distributed mode

## Risks and Mitigations

| Risk | Mitigation |
|------|------------|
| Memory pressure on large files | Windowed accumulation, streaming output |
| Pileup edge cases | Extensive testing against GATK outputs |
| HTSJDK API changes | Pin to specific version, integration tests |
| Different results than GATK | Validation suite comparing outputs |

## Success Criteria

1. Single-pass collection of all metrics currently requiring 3 passes
2. Results within 1% of GATK tool outputs
3. 2.5x+ speedup on typical WGS files
4. Memory usage under 2GB for 30x WGS
5. Progress reporting comparable to current tools
6. **Chromosome-level metrics retained** for both coverage and callable loci:
   - Per-contig coverage histograms for visualization
   - Per-contig callable/uncallable base counts by state
   - Compatible with existing `ContigSummary` visualization in UI

## References

- [GATK CallableLoci Source](https://github.com/broadinstitute/gatk/blob/master/src/main/java/org/broadinstitute/hellbender/tools/walkers/coverage/CallableLoci.java) - Reference implementation for callable state logic
- [HTSJDK SamLocusIterator](https://samtools.github.io/htsjdk/javadoc/htsjdk/htsjdk/samtools/util/SamLocusIterator.html)
- [GATK CollectWgsMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360036856051-CollectWgsMetrics-Picard)
- [GATK CollectMultipleMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360036480312-CollectMultipleMetrics-Picard)
- [Picard Metrics Source](https://github.com/broadinstitute/picard/tree/master/src/main/java/picard/analysis)
