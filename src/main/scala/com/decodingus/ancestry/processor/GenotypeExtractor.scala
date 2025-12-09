package com.decodingus.ancestry.processor

import com.decodingus.analysis.GatkRunner
import htsjdk.variant.vcf.VCFFileReader

import java.io.File
import java.nio.file.{Files, Path}
import scala.jdk.CollectionConverters.*

/**
 * Result of genotype extraction from a BAM/CRAM file.
 *
 * @param genotypes Map of SNP ID (chr:pos format) to genotype value:
 *                  0 = homozygous reference
 *                  1 = heterozygous
 *                  2 = homozygous alternate
 *                  -1 = no call / missing
 * @param snpsWithGenotype Number of SNPs with valid genotype calls
 * @param snpsMissing Number of SNPs with no call
 * @param outputVcf Path to the output VCF with called genotypes
 */
case class GenotypeExtractionResult(
  genotypes: Map[String, Int],
  snpsWithGenotype: Int,
  snpsMissing: Int,
  outputVcf: File
)

/**
 * Extracts genotypes from BAM/CRAM at specified SNP positions.
 * Uses GATK HaplotypeCaller in force-call mode with diploid ploidy.
 */
class GenotypeExtractor {

  /**
   * Extract genotypes at specified SNP sites from a BAM/CRAM file.
   *
   * @param bamPath Path to the BAM/CRAM file
   * @param referencePath Path to the reference genome FASTA
   * @param sitesVcf VCF file containing SNP positions to genotype
   * @param onProgress Progress callback (message, current, total)
   * @param outputDir Optional directory for output VCF (uses temp if None)
   * @return Either error message or extraction result
   */
  def extractGenotypes(
    bamPath: String,
    referencePath: String,
    sitesVcf: File,
    onProgress: (String, Double, Double) => Unit,
    outputDir: Option[Path]
  ): Either[String, GenotypeExtractionResult] = {

    onProgress("Preparing genotype extraction...", 0.0, 1.0)

    // Ensure BAM index exists
    GatkRunner.ensureIndex(bamPath) match {
      case Left(error) => return Left(s"Failed to index BAM: $error")
      case Right(_) => // continue
    }

    // Create output VCF path
    val outputVcf = outputDir match {
      case Some(dir) =>
        Files.createDirectories(dir)
        dir.resolve("ancestry_genotypes.vcf").toFile
      case None =>
        val temp = File.createTempFile("ancestry_genotypes", ".vcf")
        temp.deleteOnExit()
        temp
    }

    onProgress("Indexing sites VCF...", 0.05, 1.0)

    // Index the sites VCF if needed
    val sitesVcfIndex = new File(sitesVcf.getAbsolutePath + ".tbi")
    if (!sitesVcfIndex.exists()) {
      GatkRunner.run(Array(
        "IndexFeatureFile",
        "-I", sitesVcf.getAbsolutePath
      )) match {
        case Left(error) => return Left(s"Failed to index sites VCF: $error")
        case Right(_) => // continue
      }
    }

    onProgress("Calling genotypes at ancestry-informative sites...", 0.1, 1.0)

    // Call genotypes with GATK HaplotypeCaller
    // Using diploid ploidy for autosomal SNPs
    val args = Array(
      "HaplotypeCaller",
      "-I", bamPath,
      "-R", referencePath,
      "-O", outputVcf.getAbsolutePath,
      "-L", sitesVcf.getAbsolutePath,
      "--alleles", sitesVcf.getAbsolutePath,
      "--disable-sequence-dictionary-validation", "true",
      "--sample-ploidy", "2",
      "--standard-min-confidence-threshold-for-calling", "10.0",
      "--force-call-filtered-alleles", "true",
      "--output-mode", "EMIT_ALL_SITES"  // Include ref calls and no-calls
    )

    val progressCallback: Option[(String, Double) => Unit] = Some { (msg, done) =>
      onProgress(msg, 0.1 + done * 0.8, 1.0)
    }

    GatkRunner.runWithProgress(
      args,
      progressCallback,
      None,
      None
    ) match {
      case Left(error) => Left(s"Genotype calling failed: $error")
      case Right(_) =>
        onProgress("Parsing genotypes...", 0.9, 1.0)
        parseGenotypes(outputVcf)
    }
  }

  /**
   * Parse called genotypes from VCF into a map.
   */
  private def parseGenotypes(vcfFile: File): Either[String, GenotypeExtractionResult] = {
    try {
      val reader = new VCFFileReader(vcfFile, false)
      var withGenotype = 0
      var missing = 0

      val genotypes = reader.iterator().asScala.map { vc =>
        val snpId = s"${vc.getContig}:${vc.getStart}"
        val genotype = vc.getGenotypes.get(0)

        val genotypeValue = if (genotype.isNoCall || genotype.isFiltered) {
          missing += 1
          -1
        } else if (genotype.isHomRef) {
          withGenotype += 1
          0
        } else if (genotype.isHet) {
          withGenotype += 1
          1
        } else if (genotype.isHomVar) {
          withGenotype += 1
          2
        } else {
          // Partial call or other edge case
          missing += 1
          -1
        }

        snpId -> genotypeValue
      }.toMap

      reader.close()

      Right(GenotypeExtractionResult(
        genotypes = genotypes,
        snpsWithGenotype = withGenotype,
        snpsMissing = missing,
        outputVcf = vcfFile
      ))
    } catch {
      case e: Exception =>
        Left(s"Failed to parse genotypes: ${e.getMessage}")
    }
  }

  /**
   * Extract genotypes for a subset of chromosomes (for parallelization).
   */
  def extractGenotypesForChromosomes(
    bamPath: String,
    referencePath: String,
    sitesVcf: File,
    chromosomes: List[String],
    onProgress: (String, Double, Double) => Unit,
    outputDir: Option[Path]
  ): Either[String, GenotypeExtractionResult] = {

    // Create interval list for specified chromosomes
    val intervalArg = chromosomes.mkString(" -L ", " -L ", "")

    val outputVcf = outputDir match {
      case Some(dir) =>
        Files.createDirectories(dir)
        dir.resolve(s"ancestry_genotypes_${chromosomes.head}.vcf").toFile
      case None =>
        val temp = File.createTempFile("ancestry_genotypes", ".vcf")
        temp.deleteOnExit()
        temp
    }

    val args = Array(
      "HaplotypeCaller",
      "-I", bamPath,
      "-R", referencePath,
      "-O", outputVcf.getAbsolutePath,
      "--alleles", sitesVcf.getAbsolutePath,
      "--disable-sequence-dictionary-validation", "true",
      "--sample-ploidy", "2",
      "--standard-min-confidence-threshold-for-calling", "10.0",
      "--force-call-filtered-alleles", "true",
      "--output-mode", "EMIT_ALL_SITES"
    ) ++ chromosomes.flatMap(chr => Array("-L", chr))

    GatkRunner.run(args) match {
      case Left(error) => Left(s"Genotype calling failed: $error")
      case Right(_) => parseGenotypes(outputVcf)
    }
  }
}
