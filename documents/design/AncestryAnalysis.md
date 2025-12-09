# Ancestry Analysis Design

## Overview

Population percentage estimation using autosomal DNA (atDNA) with reference panels from 1000 Genomes and HGDP/SGDP. Provides ADMIXTURE-style ancestry breakdowns at sub-continental granularity.

## Algorithm

### Approach: PCA Projection + Gaussian Mixture Model

Instead of running ADMIXTURE directly (GPL-licensed, computationally intensive), we use a projection-based supervised ancestry estimation:

1. **Pre-compute reference structure** (offline, one-time):
   - Run PCA on merged 1000G + HGDP reference genotypes
   - Calculate population centroids in PCA space
   - Calculate population covariance matrices (diagonal approximation)
   - Store as binary files for efficient loading

2. **Per-sample analysis**:
   - Call genotypes at panel SNP positions using GATK HaplotypeCaller
   - Center genotypes using reference means
   - Project onto PCA space using pre-computed loadings
   - Calculate Mahalanobis distance to each population centroid
   - Convert distances to probabilities via multivariate Gaussian PDF
   - Normalize to percentage assignments

### Two-Tier Panel Strategy

| Panel | SNP Count | Use Case | Runtime |
|-------|-----------|----------|---------|
| **AIMs** | ~5,000 | Quick screening, low coverage | 2-5 min |
| **Genome-Wide** | ~500,000 | Detailed analysis, high coverage WGS | 20-30 min |

## Reference Populations

33 populations organized into 9 super-populations:

### European (5)
- CEU: Northwestern European
- FIN: Finnish
- GBR: British
- IBS: Iberian
- TSI: Tuscan

### African (5)
- YRI: Yoruba (West Africa)
- LWK: Luhya (East Africa)
- ESN: Esan (West Africa)
- MSL: Mende (West Africa)
- GWD: Gambian (West Africa)

### East Asian (5)
- CHB: Han Chinese
- JPT: Japanese
- KHV: Kinh Vietnamese
- CHS: Southern Han Chinese
- CDX: Dai Chinese

### South Asian (5)
- GIH: Gujarati
- PJL: Punjabi
- BEB: Bengali
- STU: Sri Lankan Tamil
- ITU: Indian Telugu

### Americas (4)
- MXL: Mexican
- PUR: Puerto Rican
- PEL: Peruvian
- CLM: Colombian

### West Asian (3) - from HGDP
- Druze
- Palestinian
- Bedouin

### Oceanian (2) - from HGDP
- Papuan
- Melanesian

### Central Asian (1) - from HGDP
- Yakut

### Native American (3) - from HGDP
- Maya
- Pima
- Karitiana

## Package Structure

```
src/main/scala/com/decodingus/ancestry/
├── model/
│   ├── Population.scala           # Population definitions with colors
│   ├── AncestryResult.scala       # Result model with percentages + CIs
│   └── ReferenceData.scala        # AlleleFrequencyMatrix, PCALoadings
├── reference/
│   ├── AncestryReferenceCache.scala   # Cache management
│   └── AncestryReferenceGateway.scala # Download/resolve reference data
├── processor/
│   ├── AncestryProcessor.scala    # Main pipeline orchestration
│   ├── GenotypeExtractor.scala    # GATK HaplotypeCaller wrapper
│   └── AncestryEstimator.scala    # PCA projection + GMM
└── report/
    └── AncestryReportWriter.scala # Text/JSON/HTML reports
```

## Data Flow

```
┌─────────────┐     ┌──────────────────┐     ┌─────────────────┐
│  BAM/CRAM   │────▶│ GenotypeExtractor │────▶│ Genotype Map    │
│   File      │     │ (HaplotypeCaller) │     │ chr:pos -> 0/1/2│
└─────────────┘     └──────────────────┘     └────────┬────────┘
                                                      │
┌─────────────┐     ┌──────────────────┐              │
│  Reference  │────▶│   PCA Loadings   │──────────────┤
│   Panel     │     │   + Pop Stats    │              │
└─────────────┘     └──────────────────┘              │
                                                      ▼
                                             ┌────────────────┐
                                             │ AncestryEstimator│
                                             │  - Project PCA   │
                                             │  - Calc distances│
                                             │  - Normalize %   │
                                             └────────┬─────────┘
                                                      │
                                                      ▼
                                             ┌────────────────┐
                                             │ AncestryResult │
                                             │  - Percentages │
                                             │  - Conf intervals│
                                             │  - Super-pop sum│
                                             └────────────────┘
```

## Cache Structure

```
~/.decodingus/cache/ancestry/
└── v1/
    ├── populations.json
    ├── aims/
    │   ├── GRCh38_sites.vcf.gz[.tbi]
    │   ├── GRCh37_sites.vcf.gz[.tbi]
    │   ├── allele_freqs.bin          # ~500KB
    │   └── pca_loadings.bin
    └── genome-wide/
        ├── GRCh38_sites.vcf.gz[.tbi]
        ├── GRCh37_sites.vcf.gz[.tbi]
        ├── allele_freqs.bin          # ~50MB
        └── pca_loadings.bin
```

## Binary File Formats

### AlleleFrequencyMatrix (`allele_freqs.bin`)

```
Header:
  magic:    4 bytes  (0x41464D58 = "AFMX")
  version:  2 bytes  (1)
  numPops:  2 bytes
  numSnps:  4 bytes

Population codes: numPops × 32 bytes (null-padded strings)
SNP IDs:          numSnps × 32 bytes (null-padded strings)
Frequencies:      numPops × numSnps × 4 bytes (floats, SNP-major)
```

### PCALoadings (`pca_loadings.bin`)

```
Header:
  magic:         4 bytes  (0x50434C44 = "PCLD")
  version:       2 bytes  (1)
  numSnps:       4 bytes
  numComponents: 2 bytes
  numPops:       2 bytes

SNP IDs:      numSnps × 32 bytes
SNP means:    numSnps × 4 bytes (floats)
Loadings:     numSnps × numComponents × 4 bytes (floats)
Pop codes:    numPops × 32 bytes
Centroids:    numPops × numComponents × 4 bytes
Variances:    numPops × numComponents × 4 bytes
```

## Configuration

```hocon
# feature_toggles.conf
ancestry-analysis {
  enabled = true
  default-panel = "aims"
  min-snps-aims = 3000
  min-snps-genome-wide = 100000
  display-threshold = 0.5
  reference-version = "v1"
}
```

## Result Model

```scala
case class AncestryResult(
  panelType: String,                     // "aims" or "genome-wide"
  snpsAnalyzed: Int,
  snpsWithGenotype: Int,
  snpsMissing: Int,
  percentages: List[PopulationPercentage],
  superPopulationSummary: List[SuperPopulationPercentage],
  confidenceLevel: Double,               // 0-1 based on data quality
  analysisVersion: String,
  referenceVersion: String,
  pcaCoordinates: Option[List[Double]]   // First 3 PCs for visualization
)

case class PopulationPercentage(
  populationCode: String,
  populationName: String,
  percentage: Double,                    // 0-100
  confidenceLow: Double,
  confidenceHigh: Double,
  rank: Int
)

case class SuperPopulationPercentage(
  superPopulation: String,               // "European", "African", etc.
  percentage: Double,
  populations: List[String]
)
```

## Usage Example

```scala
val processor = new AncestryProcessor()

processor.analyze(
  bamPath = "/path/to/sample.bam",
  libraryStats = libraryStats,
  panelType = AncestryPanelType.Aims,
  onProgress = (msg, curr, total) => println(s"$msg: ${(curr/total*100).toInt}%"),
  artifactContext = Some(artifactCtx)
) match {
  case Right(result) =>
    println(s"Primary ancestry: ${result.superPopulationSummary.head.superPopulation}")
    result.percentages.filter(_.percentage > 1.0).foreach { p =>
      println(f"  ${p.populationName}: ${p.percentage}%.1f%%")
    }
  case Left(error) =>
    println(s"Analysis failed: $error")
}
```

## Reference Data Preparation

The reference panel preparation is an offline process (not part of the application):

1. Download 1000 Genomes phase 3 VCFs
2. Download HGDP/SGDP genotypes
3. Merge and filter to high-quality biallelic SNPs (MAF > 1%)
4. Select AIMs panel using Fst-based ranking
5. Run PCA on reference samples (first 20 components)
6. Calculate per-population allele frequencies
7. Calculate population centroids and variances in PCA space
8. Serialize to binary format
9. Create GRCh37 and GRCh38 coordinate versions via liftover

## Future Enhancements

- **Phased haplotype analysis**: Use phased data for improved resolution
- **Chromosome painting**: Segment-by-segment ancestry assignment
- **Ancient DNA support**: Reference populations for aDNA samples
- **Custom reference panels**: User-defined population groups
- **IBD-based refinement**: Use IBD segments for recent ancestry
