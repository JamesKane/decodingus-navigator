# Y Chromosome Profile & Region-Aware Haplogroup Reporting

## Executive Summary

Decoding-Us Navigator now provides two powerful capabilities for Y chromosome analysis:

1. **Unified Y Chromosome Profile**: A single, authoritative view of a person's Y chromosome variants that intelligently combines data from multiple sources—whole genome sequencing (WGS), targeted Y-DNA tests (Big Y), and SNP chips. The system automatically resolves conflicts when different tests report different results for the same position.

2. **Region-Aware Quality Reporting**: Haplogroup reports now show read depth at each variant position and annotate variants with their genomic region (palindromes, pseudoautosomal regions, etc.). Quality scores are adjusted based on how reliable each region is for variant calling, helping users understand which calls are rock-solid and which should be interpreted with caution.

**Why this matters**: These features help researchers and genetic genealogists make better-informed decisions about their Y-DNA data by providing context that was previously hidden in technical details.

---

## For Genetic Genealogists: A Gentle Introduction

### What is a Y Chromosome Profile?

Think of your Y chromosome profile as your complete Y-DNA fingerprint. It contains thousands of genetic markers (SNPs) that define your paternal lineage going back thousands of years.

The challenge is that different DNA tests capture different pieces of this fingerprint:
- **WGS (Whole Genome Sequencing)**: Sees everything but at varying quality
- **Big Y / Targeted Tests**: Focuses on Y-DNA specifically, very accurate
- **SNP Chips**: Tests specific known markers, fast but limited

Previously, if you had multiple tests, you had to manually reconcile them. Now, the Y Profile system does this automatically, keeping track of where each variant call came from and how confident we are in each result.

### Understanding the New Report Columns

Your haplogroup reports now include three new pieces of information:

| Column | What It Tells You |
|--------|-------------------|
| **Depth** | How many times this position was read (e.g., "42x" means 42 reads). Higher is better. |
| **Region** | The type of genomic region (e.g., "P6 Palindrome", "Yq11.223"). Some regions are trickier to sequence accurately. |
| **Quality** | Star rating with "(adj)" if adjusted for region or depth concerns. |

### A Note on Ancestral vs Derived

When working with Y-DNA SNPs, we use specific terminology:

- **Ancestral**: The original allele—the state shared by most humans, inherited from distant ancestors
- **Derived**: The mutated allele—a change that occurred in a specific ancestor and defines your branch on the Y tree

For example, if SNP M269 shows "G→A", this means:
- **G** is the ancestral allele (the original state)
- **A** is the derived allele (the mutation that defines haplogroup R1b-M269)
- If you carry **A**, you are "derived" or "positive" for M269
- If you carry **G**, you are "ancestral" or "negative" for M269

This differs from general genomics terminology where "REF/ALT" simply refers to what's in a reference genome assembly.

### Why Do Regions Matter?

The Y chromosome isn't uniform—some parts are easy to sequence accurately, others are notoriously tricky:

| Region Type | What It Means | Reliability |
|-------------|---------------|-------------|
| **X-degenerate** | Stable, single-copy regions—the gold standard | Most reliable |
| **Normal regions** | Standard callable DNA | Very reliable |
| **Palindromes (P1-P8)** | Mirror-image sequences that can "gene convert" | Moderate concern |
| **STRs** | Short tandem repeats that can change between generations | Higher concern |
| **PAR** | Pseudoautosomal regions that recombine with X | Moderate concern |
| **XTR** | X-transposed region, 99% identical to X chromosome | Higher concern |
| **Ampliconic** | Multiple similar copies, mapping is difficult | Higher concern |

When you see a quality rating marked "(adj)", it means the raw quality was good, but the region warrants extra caution.

---

## System Overview

### Y Chromosome Profile Architecture

The Y Profile system consists of three main layers:

```
┌─────────────────────────────────────────────────────────────────┐
│                        Y Profile                                 │
│  ┌──────────────┐ ┌──────────────┐ ┌──────────────┐             │
│  │   Variant    │ │   Source     │ │   Source     │             │
│  │   (M269)     │◄┤   Call #1    │ │   Call #2    │             │
│  │              │ │   (WGS)      │ │   (Big Y)    │             │
│  │  Consensus:  │ └──────────────┘ └──────────────┘             │
│  │  DERIVED     │                                                │
│  └──────────────┘                                                │
└─────────────────────────────────────────────────────────────────┘
```

**Key Concepts:**
- **Profile**: One per biosample, contains all Y chromosome variants
- **Variant**: A specific SNP position with ancestral and derived alleles, has a consensus state
- **Source**: Where we got variant data (WGS alignment, chip, etc.)
- **Source Call**: What a specific source reported for a variant (Ancestral, Derived, or No Call)

### Region Annotation Pipeline

```
┌──────────────┐     ┌──────────────┐     ┌──────────────┐
│   ybrowse    │────▶│  Download &  │────▶│   Annotate   │
│   GIAB       │     │  Liftover    │     │   Variants   │
│   UCSC       │     │              │     │              │
└──────────────┘     └──────────────┘     └──────────────┘
     GFF3/BED         YRegionGateway       YRegionAnnotator
```

**Data Sources:**
- **ybrowse.org**: Cytobands, palindromes, STRs (GFF3 format)
- **GIAB genome-stratifications**: PAR, XTR, ampliconic regions (BED format)
- **Hardcoded**: Heterochromatin boundaries (Yq12)

---

## Y Profile: Detailed Design

### Data Model

> **Note**: The entity definitions below are simplified for documentation purposes. The actual implementations in `yprofile/model/YChromosomeProfile.scala` contain additional fields for multi-reference support, STR metadata, and migration compatibility. See the Code Reference section for file paths.

#### YChromosomeProfileEntity

The root entity for a biosample's Y chromosome data:

```scala
case class YChromosomeProfileEntity(
  id: UUID,
  biosampleId: UUID,
  consensusHaplogroup: Option[String],
  consensusHaplogroupScore: Option[Double],
  status: YProfileStatus,  // DRAFT, ACTIVE, ARCHIVED
  createdAt: Instant,
  updatedAt: Instant
)
```

#### YProfileVariantEntity

Individual variant positions:

```scala
case class YProfileVariantEntity(
  id: UUID,
  yProfileId: UUID,
  position: Long,
  ancestralAllele: String,          // The ancestral (original) allele
  derivedAllele: String,            // The derived (mutated) allele
  variantName: Option[String],      // e.g., "M269"
  consensusState: YConsensusState,  // ANCESTRAL, DERIVED, HETEROZYGOUS, NO_CALL, CONFLICT
  consensusConfidence: Double,      // 0.0 - 1.0
  concordanceWeight: Double,        // Accumulated evidence weight
  sourceCount: Int,                 // Number of sources reporting
  conflictCount: Int                // Sources in disagreement
)
```

> **Terminology Note**: In Y-DNA phylogenetics, we use "ancestral" and "derived" rather than "ref/alt". The ancestral allele is the original state inherited from ancient ancestors, while the derived allele is the mutation that defines a branch. This differs from VCF terminology where REF/ALT simply refer to the reference genome assembly.

#### YProfileSourceEntity

Data sources contributing to the profile:

```scala
case class YProfileSourceEntity(
  id: UUID,
  yProfileId: UUID,
  sourceType: YProfileSourceType,   // SANGER, WGS_SHORT_READ, WGS_LONG_READ, TARGETED_NGS, CHIP, etc.
  sourceRef: Option[String],        // External reference ID
  vendor: Option[String],           // Test vendor (e.g., "FTDNA", "Nebula")
  testName: Option[String],         // Specific test name (e.g., "Big Y-700")
  testDate: Option[LocalDateTime],  // When test was performed
  methodTier: Int,                  // Quality tier (0-5 from sourceType.snpTier)
  meanReadDepth: Option[Double],    // Average coverage
  meanMappingQuality: Option[Double],
  coveragePct: Option[Double],      // Percentage of Y covered
  variantCount: Int,                // SNP count from this source
  strMarkerCount: Int,              // STR marker count
  novelVariantCount: Int,           // Private/novel variants
  alignmentId: Option[UUID],        // Link to alignment if from BAM/CRAM
  referenceBuild: Option[String],   // GRCh38, GRCh37, CHM13v2
  importedAt: LocalDateTime
)
```

Source type weights (from `YProfileSourceType` enum, separate SNP/STR weights):
| Source Type | SNP Weight | STR Weight | Rationale |
|-------------|------------|------------|-----------|
| Sanger Sequencing | 1.0 | 0.9 | Gold standard for SNPs, good for STRs |
| Capillary Electrophoresis | 0.5 | 1.0 | Not for SNPs, gold standard for STRs |
| WGS Long Read | 0.95 | 0.7 | Excellent for SNPs, good for repeats |
| WGS Short Read | 0.85 | 0.5 | Good for SNPs, repeat estimation error-prone |
| Targeted NGS (Big Y) | 0.75 | 0.4 | Good but limited regions |
| Chip | 0.5 | 0.3 | Probe-based, limited |
| Manual | 0.3 | 0.2 | User-provided, lowest confidence |

#### YVariantSourceCallEntity

What each source reported for a variant:

```scala
case class YVariantSourceCallEntity(
  id: UUID,
  variantId: UUID,
  sourceId: UUID,
  calledAllele: String,
  callState: YConsensusState,
  quality: Option[Double],
  readDepth: Option[Int],
  evidenceWeight: Double,           // Source weight * quality modifier
  isConflicting: Boolean
)
```

### Consensus Algorithm

The concordance system resolves conflicts using weighted voting:

```
For each variant position:
  1. Collect all source calls
  2. Weight each call: source_weight × quality_modifier × region_modifier
  3. Sum weights for each call state (ANCESTRAL, DERIVED, etc.)
  4. Consensus = state with highest total weight
  5. Confidence = winning_weight / total_weight
  6. Conflict = any source disagrees with consensus
```

**Weight Calculation Formula** (implemented in `YVariantConcordance.calculateWeight()`):
```
weight = methodWeight × (1 + depthBonus) × mapQFactor × callableFactor × regionFactor
```

Where:
- `methodWeight`: From `YProfileSourceType` enum (SNP vs STR weights)
- `depthBonus`: min(sqrt(depth)/10, 1.0) - rewards higher coverage
- `mapQFactor`: min(MQ/60, 1.0) - mapping quality normalized to 60
- `callableFactor`: From `YCallableState` enum (callable=1.0, non-callable=0.5)
- `regionFactor`: From `RegionType.modifier` (X-degenerate=1.0, palindrome=0.4, etc.)

**Mapping Quality Factor** (`mapQFactor`):
- Calculated as: min(MQ / 60.0, 1.0)
- MQ60: 1.0 (highest confidence mapping)
- MQ30: 0.5
- MQ0 or missing: 1.0 (default for non-sequencing sources)

> **Note**: The concordance algorithm uses mapping quality (MQ), not variant quality (PHRED). Mapping quality indicates how confidently a read aligns to the reference, while variant quality indicates confidence in the called allele. For haplogroup *display*, star ratings are based on variant quality—see `HaplogroupReportWriter`.

**Region Modifier** (from `RegionType.modifier`):
- X-degenerate: 1.0 (gold standard single-copy)
- Normal: 1.0
- PAR: 0.5 (recombines with X)
- Palindrome: 0.4 (gene conversion risk)
- XTR: 0.3 (99% X-identical)
- Ampliconic: 0.3 (high copy number)
- STR: 0.25 (recLOH risk)
- Centromere/Heterochromatin: 0.1 (unmappable)

### Database Schema

```sql
CREATE TABLE y_chromosome_profiles (
    id UUID PRIMARY KEY,
    biosample_id UUID NOT NULL REFERENCES biosamples(id),
    consensus_haplogroup VARCHAR(100),
    consensus_haplogroup_score DOUBLE,
    status VARCHAR(20) NOT NULL DEFAULT 'DRAFT',
    created_at TIMESTAMP NOT NULL,
    updated_at TIMESTAMP NOT NULL
);

CREATE TABLE y_profile_variants (
    id UUID PRIMARY KEY,
    y_profile_id UUID NOT NULL REFERENCES y_chromosome_profiles(id),
    position BIGINT NOT NULL,
    ancestral_allele VARCHAR(100) NOT NULL,  -- Original allele state
    derived_allele VARCHAR(100) NOT NULL,    -- Mutated allele defining the branch
    variant_name VARCHAR(100),
    consensus_state VARCHAR(20) NOT NULL,    -- ANCESTRAL, DERIVED, etc.
    consensus_confidence DOUBLE NOT NULL,
    concordance_weight DOUBLE NOT NULL DEFAULT 0,
    source_count INT NOT NULL DEFAULT 0,
    conflict_count INT NOT NULL DEFAULT 0,
    UNIQUE(y_profile_id, position, ancestral_allele, derived_allele)
);

CREATE TABLE y_profile_sources (
    id UUID PRIMARY KEY,
    y_profile_id UUID NOT NULL REFERENCES y_chromosome_profiles(id),
    source_type VARCHAR(50) NOT NULL,
    alignment_id UUID REFERENCES alignments(id),
    reference_build VARCHAR(20),
    evidence_weight DOUBLE NOT NULL DEFAULT 1.0,
    variant_count INT NOT NULL DEFAULT 0,
    processed_at TIMESTAMP NOT NULL
);

CREATE TABLE y_variant_source_calls (
    id UUID PRIMARY KEY,
    variant_id UUID NOT NULL REFERENCES y_profile_variants(id),
    source_id UUID NOT NULL REFERENCES y_profile_sources(id),
    called_allele VARCHAR(100) NOT NULL,
    call_state VARCHAR(20) NOT NULL,
    quality DOUBLE,
    read_depth INT,
    evidence_weight DOUBLE NOT NULL DEFAULT 1.0,
    is_conflicting BOOLEAN NOT NULL DEFAULT FALSE,
    UNIQUE(variant_id, source_id)
);
```

---

## Region Annotations: Detailed Design

### Region Types and Quality Modifiers

```scala
enum RegionType(val modifier: Double, val displayName: String, val description: String):
  case Cytoband       extends RegionType(1.0, "Cytoband", "Chromosome band (display only)")
  case XDegenerate    extends RegionType(1.0, "X-degenerate", "Stable, single-copy regions - gold standard")
  case Normal         extends RegionType(1.0, "Normal", "Normal callable region")
  case PAR            extends RegionType(0.5, "PAR", "Pseudoautosomal region - recombines with X")
  case Palindrome     extends RegionType(0.4, "Palindrome", "Palindromic region - gene conversion risk")
  case XTR            extends RegionType(0.3, "XTR", "X-transposed region - 99% X-identical, contamination risk")
  case Ampliconic     extends RegionType(0.3, "Ampliconic", "Ampliconic region - high copy number, mapping artifacts")
  case STR            extends RegionType(0.25, "STR", "Short tandem repeat - recLOH risk")
  case Centromere     extends RegionType(0.1, "Centromere", "Centromeric region - nearly unmappable")
  case Heterochromatin extends RegionType(0.1, "Heterochromatin", "Heterochromatic region (Yq12) - unmappable")
  case NonCallable    extends RegionType(0.5, "Non-callable", "Failed callable loci criteria")
  case LowDepth       extends RegionType(0.7, "Low depth", "Read depth below threshold (<10x)")
```

Modifiers combine **multiplicatively**. A variant in a palindrome (0.4) with low depth (0.7) has combined modifier: 0.4 × 0.7 = 0.28

### GRCh38 Y Chromosome Structure

```
Position (Mbp)    Region
0 ─────────────── PAR1 (10kb - 2.78 Mbp)
  │
3 ─────────────── X-degenerate regions begin
  │               (interspersed with ampliconic)
  │
10 ────────────── XTR (X-transposed region, ~3.4 Mbp)
  │               99% identical to Xq21
  │
13 ────────────── Ampliconic regions (with palindromes P1-P8)
  │               Gene conversion risk
  │
27 ────────────── Heterochromatin begins (Yq12)
  │               ~30 Mbp of unmappable sequence
  │
57 ────────────── PAR2 (56.89 - 57.22 Mbp)
```

### Data Sources

The system uses **native files for each reference build** when available, avoiding liftover errors and providing the best quality annotations.

#### GRCh38 (hg38)

| Type | Source | Format |
|------|--------|--------|
| Cytobands | ybrowse.org `cytobands_hg38.gff3` | GFF3 |
| Palindromes | ybrowse.org `palindromes_hg38.gff3` | GFF3 |
| STRs | ybrowse.org `str_hg38.gff3` | GFF3 |
| PAR | GIAB genome-stratifications `GRCh38_chrY_PAR.bed` | BED |
| XTR | GIAB genome-stratifications `GRCh38_chrY_XTR.bed` | BED |
| Ampliconic | GIAB genome-stratifications `GRCh38_chrY_ampliconic.bed` | BED |

#### GRCh37 (hg19)

| Type | Source | Format |
|------|--------|--------|
| Cytobands | ybrowse.org `cytobands_hg19.gff3` | GFF3 |
| Palindromes | ybrowse.org `palindromes_hg19.gff3` | GFF3 |
| STRs | ybrowse.org `str_hg19.gff3` | GFF3 |
| PAR | GIAB genome-stratifications `GRCh37_chrY_PAR.bed` | BED |
| XTR | GIAB genome-stratifications `GRCh37_chrY_XTR.bed` | BED |
| Ampliconic | GIAB genome-stratifications `GRCh37_chrY_ampliconic.bed` | BED |

#### CHM13v2.0 (hs1) - T2T Reference

CHM13v2.0 uses native annotations from the [T2T Consortium](https://github.com/marbl/CHM13), providing the highest quality Y chromosome annotations with 30+ Mbp more sequence than GRCh38.

| Type | Source | Format | Notes |
|------|--------|--------|-------|
| Cytobands | T2T S3 `chm13v2.0_cytobands_allchrs.bed` | BED | All chromosomes |
| Palindromes | T2T S3 `chm13v2.0Y_inverted_repeats_v1.bed` | BED | Native T2T annotation |
| Ampliconic | T2T S3 `chm13v2.0Y_amplicons_v1.bed` | BED | Native T2T annotation |
| Sequence Class | T2T S3 `chm13v2.0_chrXY_sequence_class_v1.bed` | BED | X-DEG, AMPL, HET, etc. |
| AZF/DYZ | T2T S3 `chm13v2.0Y_AZF_DYZ_v1.bed` | BED | Clinical regions |
| PAR | GIAB `CHM13v2.0_chrY_PAR.bed` | BED | |
| XTR | GIAB `CHM13v2.0_chrY_XTR.bed` | BED | |
| STRs | Liftover from GRCh38 | GFF3 | No native T2T STR file |

**Why CHM13v2.0 matters:** The T2T-Y assembly is the first complete human Y chromosome sequence (62,460,029 bp), revealing the full structure of ampliconic gene families (TSPY, DAZ, RBMY) and resolving the heterochromatic Yq12 region that was missing from GRCh38.

Files are downloaded once and cached at `~/.decodingus/cache/yregions/`.

### File Parsing

**GFF3 Format:**
```
##gff-version 3
chrY    ybrowse palindrome  14969754    15077740    .   .   .   Name=P6;Note=Gene conversion hotspot
```

**BED Format (0-based, half-open):**
```
chrY    10000   2781479 PAR1
chrY    56887902    57217415    PAR2
```

The `RegionFileParser` handles both formats and converts BED coordinates to 1-based inclusive (matching VCF/GFF3 conventions).

### Binary Search Lookup

Regions are stored sorted by start position. Lookup uses binary search for O(log n) performance:

```scala
private def findOverlapping(regions: IndexedSeq[GenomicRegion], position: Long): Option[GenomicRegion] = {
  var lo = 0
  var hi = regions.length - 1
  var result: Option[GenomicRegion] = None

  while (lo <= hi) {
    val mid = lo + (hi - lo) / 2
    val region = regions(mid)

    if (region.start <= position) {
      if (position <= region.end) result = Some(region)
      lo = mid + 1
    } else {
      hi = mid - 1
    }
  }
  result
}
```

### Integration with Reports

The `HaplogroupReportWriter` accepts optional parameters:
- `snpCallInfo: Option[Map[Long, SnpCallInfo]]` - depth and quality per position
- `yRegionAnnotator: Option[YRegionAnnotator]` - region lookup

When provided, reports include:
1. **Depth column** showing read coverage (e.g., "42x")
2. **Region column** showing the region name (e.g., "P6 Palindrome")
3. **Adjusted quality** stars with "(adj)" suffix when modifiers apply
4. **Legend** explaining the quality modifiers

---

## Code Reference

### Key Files

| Component | File | Purpose |
|-----------|------|---------|
| Y Profile Models | `yprofile/model/YChromosomeProfile.scala` | Entity definitions, enums |
| Y Profile Entity | `repository/YChromosomeProfileRepository.scala` | Profile CRUD |
| Y Variant Entity | `repository/YProfileVariantRepository.scala` | Variant CRUD |
| Y Source Entity | `repository/YProfileSourceRepository.scala` | Source tracking |
| Source Calls | `repository/YVariantSourceCallRepository.scala` | Per-source calls |
| Concordance | `yprofile/concordance/YVariantConcordance.scala` | Consensus calculation |
| Profile Service | `yprofile/service/YProfileService.scala` | High-level operations |
| Region Parser | `refgenome/RegionFileParser.scala` | GFF3/BED parsing |
| Region Cache | `refgenome/YRegionCache.scala` | File caching |
| Region Gateway | `refgenome/YRegionGateway.scala` | Download & liftover |
| Region Annotator | `refgenome/YRegionAnnotator.scala` | Position lookup |
| Enriched Call | `haplogroup/model/EnrichedVariantCall.scala` | Call + annotations |
| Report Writer | `haplogroup/report/HaplogroupReportWriter.scala` | Text report |
| Report Dialog | `ui/components/HaplogroupReportDialog.scala` | UI display |

### Tests

| Test File | Coverage |
|-----------|----------|
| `refgenome/RegionFileParserSpec.scala` | GFF3/BED parsing, coordinate conversion |
| `refgenome/YRegionAnnotatorSpec.scala` | Region lookup, modifier calculation |
| `repository/YChromosomeProfileRepositorySpec.scala` | Profile persistence |
| `repository/YProfileVariantRepositorySpec.scala` | Variant CRUD |
| `repository/YProfileSourceRepositorySpec.scala` | Source tracking |
| `repository/YVariantSourceCallRepositorySpec.scala` | Call persistence |
| `yprofile/concordance/YVariantConcordanceSpec.scala` | Consensus algorithm |
| `yprofile/service/YProfileServiceSpec.scala` | Integration tests |

---

## Future Enhancements

| Enhancement | Status | Notes |
|-------------|--------|-------|
| **Region-Aware Concordance** | ✅ Implemented | `regionModifier` parameter added to `YVariantConcordance.calculateWeight()` and `SourceCallInput` |
| **Callable Loci Integration** | ✅ Implemented | Uses `callable_loci.bed` from GATK via `YCallableState` enum |
| **X-Degenerate Annotation** | ✅ Implemented | Parsed from T2T `sequence_class` file for CHM13v2; marked as `RegionType.XDegenerate` |
| **Interactive Region Visualization** | ❌ Planned | Show regions on a chromosome ideogram in the UI |
| **Export with Annotations** | ✅ Implemented | `HaplogroupReportWriter.writeCsvReport()` exports CSVs with region annotations |

### CSV Export Details

The `writeCsvReport()` method creates three CSV files with full region annotation:

| File | Contents |
|------|----------|
| `ydna_snp_details.csv` | SNPs along the predicted haplogroup path |
| `ydna_novel_snps.csv` | Novel/unplaced SNPs |
| `ydna_novel_indels.csv` | Novel indels with STR marker info |

Each CSV includes columns: `region_type`, `region_name`, `quality_modifier`, `adjusted_quality`
