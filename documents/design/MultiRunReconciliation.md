# Multi-Run Reconciliation Design

## Problem Statement

A single subject (person) may have multiple sequencing runs from different sources:
- Different vendors (Nebula, Dante, Ancestry, 23andMe)
- Different sequencing technologies (WGS, WES, SNP arrays)
- Different coverage depths
- Different reference builds
- Re-sequencing over time

These runs may produce slightly different results due to:
1. **Technical variation** - Coverage gaps, sequencing errors, different quality thresholds
2. **Legitimate biological variation** - Somatic mutations, heteroplasmy in mtDNA
3. **Sample mix-ups** - Wrong person entirely

We need to:
- Reconcile minor differences to produce a consensus
- Detect major incompatibilities that indicate sample problems
- Flag runs that don't belong to the same individual

## Haplogroup Reconciliation

### Compatible Differences (Expected)

Runs may legitimately differ in haplogroup depth due to coverage:

```
Run A (30x WGS):     R-FGC11134 > R-CTS4466 > R-BY18291
Run B (15x WGS):     R-FGC11134 > R-CTS4466
Run C (SNP array):   R-M269
```

These are **compatible** - each is a valid ancestor of the most specific call.

**Reconciliation Strategy:**
- Accept the deepest haplogroup that has supporting evidence from at least one high-confidence run
- Track "confirmed depth" vs "best estimate depth"
- Show path with confidence indicators per branch

### Incompatible Differences (Problem)

Runs that diverge at a branch point indicate a serious problem:

```
Run A: R-DF13 > R-FGC11134 > R-CTS4466
Run B: R-DF13 > R-L1065 > R-S5668
              ^--- divergence point
```

These haplogroups share R-DF13 ancestry but split into incompatible branches.

**Detection Strategy:**
- Find the Last Common Ancestor (LCA) of reported haplogroups
- If LCA is significantly shallower than individual calls, flag incompatibility
- Calculate "divergence depth" = (average depth of calls) - (LCA depth)
- High divergence depth = likely sample mix-up

### Severity Levels

| Level | Description | Example | Action |
|-------|-------------|---------|--------|
| **None** | Calls are on same branch | R-BY18291 vs R-CTS4466 | Accept deepest |
| **Minor** | Calls differ at tips only | Sibling terminal branches | Review manually |
| **Major** | Calls diverge at ancient branch | R-DF13 children split | Flag for review |
| **Critical** | Calls are entirely incompatible | R1b vs I2a | Reject as different person |

### Quantitative Metrics

**Branch Compatibility Score:**
```
compatibility = LCA_depth / max(depth_A, depth_B)
```
- 1.0 = One is ancestor of other (fully compatible)
- 0.8+ = Minor tip differences
- 0.5-0.8 = Significant divergence, needs review
- <0.5 = Likely different individuals

**SNP-level Concordance:**
```
concordance = matching_calls / (matching_calls + conflicting_calls)
```
- Ignore no-calls (missing data)
- 0.99+ = Same individual, technical variation
- 0.95-0.99 = Same individual, possible somatic/heteroplasmy
- <0.95 = Likely different individuals

## MT-DNA Specific Considerations

### Heteroplasmy
mtDNA can legitimately have mixed populations within one individual:
- Age-related heteroplasmy accumulation
- Tissue-specific differences
- Runs from blood vs saliva may differ

**Handling:**
- Track variant allele frequency (VAF) at each position
- Flag positions with VAF 10-90% as heteroplasmic
- Don't treat heteroplasmic differences as incompatible
- Report heteroplasmy separately

### Haplogroup Assignment with Heteroplasmy
- If a defining SNP is heteroplasmic, weight by VAF
- Example: Position 16519 at 30% C / 70% T
  - 70% confidence in derived state
  - Propagate uncertainty to haplogroup confidence

## Y-DNA Specific Considerations

### Coverage Variation in Difficult Regions
chrY has many repetitive and difficult-to-sequence regions:
- Palindromic regions
- Ampliconic sequences
- Pseudoautosomal regions (PAR)

**Handling:**
- Weight SNPs by mappability/complexity score
- Don't penalize no-calls in known difficult regions
- Track "high-confidence callable bases" per run

### Big Y vs WGS Comparison
FTDNA Big Y targets specific regions; WGS covers everything but may miss depth.

**Handling:**
- Normalize by expected coverage per technology
- Big Y: expect calls at targeted sites
- WGS: expect broader but potentially shallower coverage

## Implementation Architecture

### Data Model Extensions

```scala
case class HaplogroupCall(
  haplogroup: String,
  confidence: Double,
  supportingSnps: Int,
  conflictingSnps: Int,
  noCalls: Int,
  runId: String,
  technology: SequencingTechnology,
  meanCoverage: Double
)

case class ReconciliationResult(
  consensusHaplogroup: String,
  confidence: Double,
  compatibilityLevel: CompatibilityLevel,
  allCalls: List[HaplogroupCall],
  divergencePoint: Option[String],  // LCA if divergent
  snpConflicts: List[SnpConflict],
  warnings: List[String]
)

case class SnpConflict(
  position: Long,
  snpName: String,
  calls: Map[String, AlleleCall],  // runId -> call
  resolution: ConflictResolution
)

enum ConflictResolution:
  case AcceptMajority
  case AcceptHigherQuality
  case AcceptHigherCoverage
  case Unresolved
  case Heteroplasmy  // mtDNA only

enum CompatibilityLevel:
  case Compatible      // Same branch, different depths
  case MinorDivergence // Tip-level differences
  case MajorDivergence // Branch-level split
  case Incompatible    // Different individuals
```

### Reconciliation Algorithm

```
function reconcileHaplogroups(calls: List[HaplogroupCall]): ReconciliationResult
  1. Build path-to-root for each haplogroup call
  2. Find LCA of all calls
  3. Calculate compatibility score
  4. If compatible:
     - Return deepest call as consensus
     - Merge confidence from supporting runs
  5. If divergent:
     - Flag divergence point
     - Calculate SNP-level concordance
     - If concordance < threshold:
       - Mark as likely different individuals
     - Else:
       - Mark for manual review
       - Suggest possible causes (coverage, heteroplasmy)
```

### UI/UX Considerations

**Subject Dashboard:**
- Show consensus haplogroup prominently
- Traffic light indicator for cross-run compatibility
  - Green: All runs agree
  - Yellow: Minor differences, reconciled
  - Red: Major incompatibility detected

**Run Comparison View:**
- Side-by-side haplogroup paths
- Highlight divergence points
- SNP-level conflict table

**Alerts:**
- "Sample Verification Needed" for major divergences
- Suggest re-extraction or re-sequencing if incompatible
- Link to specific conflicting positions

## Identity Verification Metrics

Beyond haplogroups, use autosomal data for identity verification:

### Kinship Coefficient
If autosomal data available:
- Self-vs-self should have kinship ~0.5
- Different individuals typically <0.05

### STR Profiles
Y-STR profiles should match exactly (or within mutation rate):
- Calculate genetic distance between runs
- Flag if distance > 2 steps on any marker

### Fingerprint SNPs
Curated set of highly polymorphic, easy-to-call SNPs:
- 50-100 autosomal SNPs
- Should match 100% between runs of same individual
- Mismatch rate indicates sample problem

## Future Considerations

### Machine Learning Approach
Train a classifier on known same/different individual pairs:
- Features: SNP concordance, haplogroup compatibility, coverage correlation
- Output: Probability of same individual

### Cross-Reference with Relatives
If relatives are in system:
- Verify expected relationships still hold
- Parent-child pairs should share haplogroup branch
- Siblings should have consistent Y-DNA (if male)

### Audit Trail
Maintain history of:
- All runs associated with subject
- Reconciliation decisions made
- Manual overrides by user
- Confidence changes over time
