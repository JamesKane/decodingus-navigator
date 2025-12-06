# Decoding-Us Navigator

Decoding-Us Navigator is an edge-computing companion application to [decoding.us](https://decoding.us). It leverages the Genome Analysis Toolkit (GATK) to analyze BAM/CRAM files directly on your local machine, empowering citizen scientists with advanced bioinformatics capabilities while preserving privacy.

## Privacy-Preserving Analysis

The application ensures user privacy by performing all analysis locally. Only anonymized summary information is optionally shared, including:
- Haplogroup assignments
- General coverage statistics for quality control
- Autosomal DNA matches with other researchers in the Federation (coming soon)

Data sharing uses the AT Protocol Personal Data Store (PDS) for user-controlled data ownership.

## Goal

Decoding-Us Navigator simplifies complex bioinformatics command-line tools by wrapping them in an intuitive interface. It is designed for hobbyists and citizen scientists, making advanced genetic analysis accessible without requiring programming expertise.

## Cross-Platform Compatibility

Built on the Java Virtual Machine (JVM) with ScalaFX, Decoding-Us Navigator runs on macOS, Windows, and Linux with a consistent user experience.

## Features

### Workspace Management
- Create and manage multiple projects and biosamples
- Drag-and-drop project membership management
- Persistent workspace saved locally
- Search and filter across projects and subjects

### Sequencing Data Management
- Import BAM/CRAM files via file picker or drag-and-drop
- Support for local files and cloud URLs (HTTP/S3)
- Automatic SHA-256 checksum calculation
- Platform detection (Illumina, PacBio, Oxford Nanopore, MGI, Ion Torrent, Complete Genomics)
- Test type classification (WGS, WES, HiFi, CLR, Nanopore, Targeted Panel)

### Analysis Capabilities
- **Library Statistics**: Rapid BAM/CRAM scanning for sample metadata, reference build detection, read length distribution, and insert size metrics
- **WGS Metrics**: Comprehensive coverage analysis including mean coverage, coverage distribution, and depth thresholds (1x-100x)
- **Callable Loci**: Per-contig analysis identifying callable bases, coverage gaps, and mapping quality issues with SVG visualizations
- **Haplogroup Determination**: Y-DNA and MT-DNA haplogroup analysis with multiple tree providers (FTDNA, DecodingUs)
- **Private SNP Detection**: Identify novel SNPs unique to an individual after haplogroup determination
- **Liftover**: Automatic coordinate conversion between reference builds (GRCh38, GRCh37, CHM13v2)

### Reference Genome Management
- Automatic reference genome download and caching
- Support for GRCh38, GRCh37, and CHM13v2
- Configurable local paths and cache directory
- Download prompts with size estimates

### Analysis Caching
- SHA-256 based result caching to prevent redundant analysis
- Cached results for coverage, WGS metrics, library stats, and contig summaries

### Optional Cloud Integration
- AT Protocol authentication and PDS integration
- Workspace sync from personal data store
- Optional upload of summary data with user consent

## Requirements

- Java 17 or later
- 4GB RAM minimum (8GB recommended for large BAM files)

## Building

```bash
# Compile
sbt compile

# Run
sbt run

# Create fat JAR
sbt assembly

# Run tests
sbt test
```