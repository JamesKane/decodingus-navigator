# Decoding-Us Navigator

'Decoding-Us Navigator' is an Edge-computing companion application to [https://decoding.us.com](https://decoding.us.com). It leverages the Genome Analysis Toolkit (GATK) to analyze BAM/CRAM files directly on a user's local machine, empowering citizen scientists with advanced bioinformatics capabilities while preserving privacy.

## Privacy-Preserving Analysis

The application ensures user privacy by performing local analysis. Only anonymized summary information is shared, including:
- Haplogroup assignments
- Autosomal DNA matches with other researchers in the Federation (feature coming soon)
- General coverage statistics for shared quality control using the AT Protocol PDS

## Goal

The primary goal of Decoding-Us Navigator is to simplify complex bioinformatics command-line tools by wrapping them in an intuitive and user-friendly interface. It is specifically designed for hobbyists, making advanced genetic analysis accessible to a broader audience.

## Cross-Platform Compatibility

Developed on the Java Virtual Machine (JVM), Decoding-Us Navigator offers easy migration to various operating systems, ensuring a consistent user experience across different platforms.

## Current Features

- **Haplogroup Analysis:** Determine your haplogroup from BAM/CRAM files.
- **WGS Metrics:** Generate whole-genome sequencing metrics for quality control.
- **Library Statistics:** Obtain detailed statistics about your sequencing libraries.
- **Private SNP Processing:** Analyze private SNPs within your genome.
- **Callable Loci Processing:** Identify callable loci in your sequencing data.
- **Liftover:** Convert genomic coordinates between different reference genome assemblies.