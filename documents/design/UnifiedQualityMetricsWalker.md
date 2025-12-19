# Unified Quality Metrics Walker

## Status: Phase 1 & 2 Complete, Phase 3 Ready

**Last Updated**: 2025-12-18

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

### Phase 1: Core Walker ✅ COMPLETE
- [x] Implement `UnifiedMetricsWalker` with HTSJDK `SamReader`
- [x] Read counter collector (total, aligned, paired, proper pairs)
- [x] Insert size accumulator with histogram (replaces R-dependent CollectInsertSizeMetrics)
- [x] Mean read length calculation
- [x] `UnifiedMetricsProcessor` with artifact output

**Implemented in:**
- `src/main/scala/com/decodingus/analysis/UnifiedMetricsWalker.scala`
- `src/main/scala/com/decodingus/analysis/UnifiedMetricsProcessor.scala`

**This phase eliminates the R dependency by replacing CollectMultipleMetrics.**

### Phase 2: Coverage & Callable Loci ✅ COMPLETE
- [x] Windowed coverage accumulator
- [x] Callable loci state machine (per GATK CallableLoci logic)
- [x] BED file output for callable regions (per-contig for chrY analysis)
- [x] Coverage percentile calculations (PCT_1X, PCT_10X, etc.)

**Implemented in:**
- `src/main/scala/com/decodingus/analysis/CoverageCallableWalker.scala`
- `src/main/scala/com/decodingus/analysis/CoverageCallableProcessor.scala`

#### Phase 2 Detailed Design

The key insight is that Phase 2 requires **position-level iteration** (pileup), not just read-level iteration as in Phase 1. This is fundamentally different because:

1. **Read metrics** (Phase 1): One iteration per read, ~800M reads for 30x WGS
2. **Coverage/Callable** (Phase 2): One iteration per position, ~3.1B positions for GRCh38

The existing `CallableLociProcessor` runs GATK's `CallableLoci` per-contig (25 separate GATK invocations). A unified walker must handle this more efficiently.

##### Two-Pass Strategy (Recommended)

Given the fundamentally different iteration patterns, the optimal approach is **two passes with shared I/O**:

```
Pass 1: Read-level (UnifiedMetricsWalker - COMPLETE)
        ├── Read counts, alignment stats
        ├── Insert size histogram
        └── Read length distribution

Pass 2: Position-level (CoverageCallableWalker - NEW)
        ├── Coverage histogram per position
        ├── Callable state per position
        ├── Per-contig summaries
        └── PCT_1X, PCT_10X, etc.
```

This still achieves the goal of eliminating the **3-pass problem** because:
- Original: `CollectMultipleMetrics (Pass 1)` + `CollectWgsMetrics (Pass 2)` + `CallableLoci (Pass 3)`
- Unified: `ReadMetrics (Pass 1)` + `CoverageCallable (Pass 2)`

But more importantly, Pass 2 replaces **both** `CollectWgsMetrics` AND `CallableLoci` in a single traversal.

##### CoverageCallableWalker Design

```scala
/**
 * Single-pass position-level walker that collects both coverage metrics
 * and callable loci state. Replaces both CollectWgsMetrics and CallableLoci.
 */
class CoverageCallableWalker {

  /**
   * Process a BAM/CRAM file using SamLocusIterator for pileup.
   * Collects coverage histogram and callable state simultaneously.
   */
  def collectCoverageAndCallable(
    bamPath: String,
    referencePath: String,
    intervals: Option[List[String]] = None,  // Optional: restrict to contigs
    callableParams: CallableLociParams = CallableLociParams(),
    onProgress: (String, Long, Long) => Unit
  ): Either[String, CoverageCallableResult]
}

/**
 * Parameters for callable loci determination.
 * Defaults match GATK CallableLoci.
 */
case class CallableLociParams(
  minDepth: Int = 4,
  maxDepth: Option[Int] = None,
  minMappingQuality: Int = 10,
  minBaseQuality: Int = 20,
  maxLowMapQ: Int = 1,
  maxFractionLowMapQ: Double = 0.1
)

/**
 * Combined result from coverage and callable analysis.
 */
case class CoverageCallableResult(
  // Global coverage metrics (replaces WgsMetrics)
  genomeTerritory: Long,
  meanCoverage: Double,
  medianCoverage: Double,
  sdCoverage: Double,
  coverageHistogram: Array[Long],  // depth 0-255+
  pct1x: Double,
  pct5x: Double,
  pct10x: Double,
  pct15x: Double,
  pct20x: Double,
  pct25x: Double,
  pct30x: Double,
  pct40x: Double,
  pct50x: Double,

  // Callable loci summary (replaces CallableLociResult)
  callableBases: Long,
  contigSummaries: List[ContigSummary],

  // Per-contig coverage for visualizations
  contigCoverage: Map[String, ContigCoverageMetrics]
)
```

##### Memory Management for Position-Level Iteration

The pileup iterator naturally handles memory by not materializing all positions at once. Key considerations:

1. **Coverage histogram**: Global `Array[Long](256)` = 2KB, trivial
2. **Per-contig callable counts**: 6 `Long` counters per contig × 25 contigs = 1.2KB
3. **Per-contig coverage histogram** (optional): 2KB × 25 = 50KB
4. **Running statistics**: Welford's algorithm for mean/variance, O(1) space

Total memory: ~60KB + `SamLocusIterator` buffer (configurable)

##### Implementation Approach

Use HTSJDK's `SamLocusIterator` which provides pileup without full GATK overhead:

```scala
import htsjdk.samtools.util.SamLocusIterator
import htsjdk.samtools.util.IntervalList

val samReader = SamReaderFactory.makeDefault()
  .referenceSequence(new File(referencePath))
  .open(new File(bamPath))

// Create interval list for main assembly contigs only
val header = samReader.getFileHeader
val intervalList = new IntervalList(header)
for (seq <- header.getSequenceDictionary.getSequences.asScala) {
  if (isMainAssemblyContig(seq.getSequenceName)) {
    intervalList.add(new Interval(seq.getSequenceName, 1, seq.getSequenceLength))
  }
}

// Use indexed lookup for better performance with intervals
val locusIterator = new SamLocusIterator(samReader, intervalList, true /* useIndex */)
locusIterator.setEmitUncoveredLoci(true)  // Need zeros for coverage calculation
locusIterator.setIncludeIndels(false)     // We only need depth, not base-level detail
locusIterator.setMappingQualityScoreCutoff(0)  // Handle MAPQ filtering ourselves for callable logic

for (locus <- locusIterator.iterator().asScala) {
  val contig = locus.getSequenceName
  val position = locus.getPosition
  val pileup = locus.getRecordAndOffsets

  // 1. Update coverage histogram
  val depth = pileup.size()
  coverageHistogram(math.min(depth, 255)) += 1

  // 2. Determine callable state (requires reference base)
  val refBase = getRefBase(contig, position)
  val state = determineCallableState(refBase, pileup, callableParams)
  updateContigCallableCounts(contig, state)

  // 3. Optional: emit BED intervals on state transitions
  if (state != previousState) {
    emitBedInterval(contig, intervalStart, position - 1, previousState)
    intervalStart = position
    previousState = state
  }
}
```

**Important SamLocusIterator configuration notes:**

1. **`setEmitUncoveredLoci(true)`**: Required for accurate coverage calculation. Without this, positions with zero coverage are skipped, making it impossible to calculate mean coverage or PCT_0X.

2. **`setMappingQualityScoreCutoff(0)`**: We handle MAPQ filtering ourselves in the callable state logic. GATK's CallableLoci uses a more nuanced approach (fraction of low-MAPQ reads), not a simple cutoff.

3. **`setIncludeIndels(false)`**: Indel tracking adds overhead and we don't need it for coverage/callable analysis.

4. **Use `IntervalList` constructor**: Restricts iteration to main assembly contigs (chr1-22, X, Y, M). This avoids processing alt contigs, decoys, and HLA which would slow down iteration significantly.

5. **Use index (`useIndex=true`)**: When specifying intervals, indexed lookup is recommended for better performance.

##### Reference Base Access

**Important**: `SamLocusIterator` does NOT provide reference bases - it only provides pileup information. For callable loci, we need to know if the reference is 'N', so we must access the reference separately.

Two approaches:

1. **Pre-load contigs**: Load each contig's reference as we enter it. ~250MB for largest chromosome.
2. **On-demand lookup**: Use `ReferenceSequenceFile.getSubsequenceAt()`. Slower but constant memory.

Recommendation: Pre-load approach, loading one contig at a time as the iterator moves through contigs.

```scala
import htsjdk.samtools.reference.ReferenceSequenceFileFactory

val referenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(new File(referencePath))
var currentContigBases: Array[Byte] = null
var currentContigName: String = ""

def getRefBase(contig: String, position: Long): Byte = {
  if (contig != currentContigName) {
    // Load new contig reference
    val seq = referenceFile.getSequence(contig)
    currentContigBases = seq.getBases
    currentContigName = contig
  }
  currentContigBases((position - 1).toInt)  // 1-based position to 0-based index
}
```

##### BED File Output

Unlike the current `CallableLociProcessor` which gets BED from GATK, we need to generate it ourselves:

```scala
class CallableRegionWriter(outputPath: Path) {
  private var currentContig: String = ""
  private var intervalStart: Long = 0
  private var currentState: CallableState = CallableState.NoCoverage

  def update(contig: String, position: Long, state: CallableState): Unit = {
    if (contig != currentContig) {
      flush()
      currentContig = contig
      intervalStart = position
      currentState = state
    } else if (state != currentState) {
      emitInterval(currentContig, intervalStart, position - 1, currentState)
      intervalStart = position
      currentState = state
    }
  }

  def flush(): Unit = { /* emit final interval */ }
}
```

##### SVG Visualization Compatibility

The current `CallableLociProcessor.binIntervals()` reads BED files to generate SVG visualizations. The new walker should either:

1. **Generate BED files** (as above) so existing SVG code works unchanged
2. **Generate binned data directly** during traversal (more efficient, but requires refactoring viz code)

Recommendation: Generate BED files for compatibility, add direct binning later as optimization.

##### Phase 2 Implementation Tasks ✅ COMPLETE

**2.1 Core Data Structures**
- [x] Reuse existing `CallableState` enum from `CallableLociQueryService`
- [x] Create `CallableLociParams` case class with GATK defaults
- [x] Create `CoverageCallableResult` unified result type
- [x] Create `ContigCoverageMetrics` for per-chromosome data

**2.2 CoverageCallableWalker**
- [x] Implement `SamLocusIterator`-based traversal
- [x] Implement callable state determination logic
- [x] Implement Welford's algorithm for online mean/variance
- [x] Implement coverage histogram accumulation
- [x] Implement per-contig callable counts
- [x] Implement progress reporting with estimated completion

**2.3 Reference Handling**
- [x] Implement per-contig reference loading
- [x] Handle N-base detection for REF_N state

**2.4 Output Generation**
- [x] Implement BED file writer with interval coalescing
- [x] Implement per-contig summary file output (GATK format)
- [x] Implement coverage percentile calculations from histogram
- [x] SVG visualization generation (reuses existing pattern)

**2.5 Testing & Validation**
- [ ] Create test suite comparing output to GATK tools
- [ ] Validate coverage metrics against CollectWgsMetrics
- [ ] Validate callable counts against CallableLoci
- [ ] Performance benchmarking vs 2-pass GATK approach

### Phase 3: Integration
- [x] Create `UnifiedMetricsProcessor` following existing processor pattern (Phase 1)
- [x] Create `CoverageCallableProcessor` wrapping the walker (Phase 2)
- [ ] Replace `MultipleMetricsProcessor` usage (remove R dependency)
- [ ] Update `WorkbenchViewModel` to use new unified processors
- [ ] Ensure compatibility with existing `ContigSummary` visualization
- [ ] Update `AnalysisCache` to handle new result types
- [ ] Deprecate `WgsMetricsProcessor` and `CallableLociProcessor`

#### Phase 3 Implementation Strategy

The integration should be **gradual and backward-compatible**:

1. **Add new processors alongside existing ones** - Don't remove anything yet
2. **Feature toggle** - Use `feature_toggles.conf` to switch between old/new
3. **Validation period** - Run both and compare results
4. **Migration** - Once validated, make new processors the default
5. **Cleanup** - Remove old processors after migration complete

```scala
// Example feature toggle usage
if (FeatureToggles.isEnabled("unified_metrics")) {
  // Use CoverageCallableProcessor
} else {
  // Use WgsMetricsProcessor + CallableLociProcessor
}
```

### Phase 4: Optimization
- [ ] Multi-threaded contig processing (process contigs in parallel)
- [ ] Memory profiling and tuning window sizes
- [ ] Progress reporting with accurate ETA based on genome position
- [ ] Direct SVG binning during traversal (avoid BED parsing)
- [ ] Investigate CRAM-specific optimizations (reference caching)

#### Phase 4 Notes

**Parallel contig processing** is the most impactful optimization. Since each contig can be processed independently:

```scala
val contigResults = contigs.par.map { contig =>
  processContig(bamPath, referencePath, contig)
}.toList
```

However, this requires careful management of:
- Reference sequence file handles (thread-safe or per-thread)
- BAM index access (typically thread-safe in HTSJDK)
- Memory for concurrent contig reference sequences

## Open Questions & Decisions

### Q1: Should Phase 1 and Phase 2 share a single traversal?

**Current decision: No (separate passes)**

While tempting to collect read-level and position-level metrics in a single pass, this would require:
- Tracking read pileup state manually (complex)
- Significantly more memory for reads spanning multiple positions
- Slower iteration due to pileup construction overhead

The two-pass approach is simpler, more maintainable, and still achieves the core goal (3 passes → 2 passes).

### Q2: How to handle CRAM files with remote references?

The current `ReferenceGateway` caches reference genomes locally. For `SamLocusIterator`, we need to:
- Ensure reference is fully downloaded before starting
- Pass the local reference path to both `SamReader` and `ReferenceSequenceFile`

### Q3: Should we generate per-contig BED files or a single combined BED?

**Decision: Per-contig BED files**

While a single combined BED would be simpler, per-contig files are needed for:
- **chrY callable regions** are used for downstream Y-chromosome analysis workflows
- Compatibility with existing SVG visualization code
- Parallel processing in future optimizations

### Q4: What about HiFi/long-read specific settings?

The current processors have special handling for PacBio HiFi:
- `minDepth = 2` instead of 4 (lower coverage is still callable with high-accuracy reads)
- `countUnpaired = true` for single-molecule reads

The unified walker should accept these as parameters via `CallableLociParams`.

### Q5: Progress reporting granularity?

Options:
1. **Per-contig**: Report after each contig completes (~25 updates for WGS)
2. **Per-million-bases**: Report every 1M positions (~3000 updates)
3. **Time-based**: Report every N seconds regardless of position

**Recommendation: Per-million-bases** for smooth progress bar movement with reasonable overhead.

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
