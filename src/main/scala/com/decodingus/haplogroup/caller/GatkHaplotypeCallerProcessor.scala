package com.decodingus.haplogroup.caller

import com.decodingus.analysis.GatkRunner
import java.io.{BufferedReader, File, FileReader, PrintWriter}
import java.nio.file.{Files, Path}
import scala.util.Using

/**
 * Result of SNP calling including output file and GATK logs.
 */
case class CallerResult(vcfFile: File, logFile: Option[File])

/**
 * Result of two-pass calling for haplogroup assignment and private SNP discovery.
 */
case class TwoPassCallerResult(
  treeSitesVcf: File,
  privateVariantsVcf: File,
  treeSitesLog: Option[File],
  privateVariantsLog: Option[File]
)

/**
 * GATK-based variant caller for haplogroup analysis.
 *
 * Uses different calling strategies based on chromosome:
 * - Y-DNA (chrY): HaplotypeCaller - positions are spread out, local assembly works well
 * - mtDNA (chrM): Mutect2 --mitochondria - optimized for dense positions and high depth
 */
class GatkHaplotypeCallerProcessor {

  /**
   * Detect the primary contig from a VCF file by reading the first data line.
   */
  private def detectContigFromVcf(vcfFile: File): String = {
    Using.resource(new BufferedReader(new FileReader(vcfFile))) { reader =>
      var line = reader.readLine()
      while (line != null && line.startsWith("#")) {
        line = reader.readLine()
      }
      if (line != null) {
        line.split("\t").headOption.getOrElse("chrM")
      } else {
        "chrM" // Default fallback
      }
    }
  }

  /**
   * Call SNPs at specified allele sites.
   *
   * @param bamPath Path to the BAM/CRAM file
   * @param referencePath Path to the reference genome
   * @param allelesVcf Sites VCF file specifying positions to call
   * @param onProgress Progress callback
   * @param outputDir Optional directory to save the called VCF and logs (if None, uses temp files)
   * @param outputPrefix Optional prefix for output files (e.g., "mtdna" or "ydna")
   * @return Either error message or CallerResult with VCF and optional log file
   */
  def callSnps(
    bamPath: String,
    referencePath: String,
    allelesVcf: File,
    onProgress: (String, Double, Double) => Unit,
    outputDir: Option[Path] = None,
    outputPrefix: Option[String] = None
  ): Either[String, CallerResult] = {
    // Ensure BAM index exists
    onProgress("Checking BAM index...", 0.0, 1.0)
    GatkRunner.ensureIndex(bamPath) match {
      case Left(error) => return Left(error)
      case Right(_) => // index exists or was created
    }

    onProgress("Calling SNPs with GATK HaplotypeCaller...", 0.1, 1.0)

    // Determine output file location
    val (vcfFile, logFile) = outputDir match {
      case Some(dir) =>
        Files.createDirectories(dir)
        val prefix = outputPrefix.getOrElse("called")
        (dir.resolve(s"${prefix}_calls.vcf").toFile, Some(dir.resolve(s"${prefix}_haplotypecaller.log").toFile))
      case None =>
        val tempFile = File.createTempFile("haplotypes", ".vcf")
        tempFile.deleteOnExit()
        (tempFile, None)
    }

    // Index the allelesVcf file (assuming it's sorted by the caller)
    val indexArgs = Array(
      "IndexFeatureFile",
      "-I", allelesVcf.getAbsolutePath
    )

    GatkRunner.run(indexArgs) match {
      case Left(error) => return Left(s"Failed to index alleles VCF: $error")
      case Right(_) => // continue
    }

    val args = Array(
      "HaplotypeCaller",
      "-I", bamPath,
      "-R", referencePath,
      "-O", vcfFile.getAbsolutePath,
      // Only call at tree sites for haplogroup assignment
      "-L", allelesVcf.getAbsolutePath,
      // Force genotyping at known tree sites
      "--alleles", allelesVcf.getAbsolutePath,
      // Relax reference validation - allows GRCh38 with/without alts, etc.
      "--disable-sequence-dictionary-validation", "true",
      // Haploid calling for mtDNA and Y-DNA
      "--sample-ploidy", "1",
      // Lower confidence threshold to capture more calls
      "--standard-min-confidence-threshold-for-calling", "10.0",
      // Include filtered alleles in output
      "--force-call-filtered-alleles", "true"
    )

    GatkRunner.run(args) match {
      case Left(error) =>
        // Save error log even on failure
        logFile.foreach { lf =>
          Using(new PrintWriter(lf)) { writer =>
            writer.println(s"GATK HaplotypeCaller failed")
            writer.println(s"Arguments: ${args.mkString(" ")}")
            writer.println(s"Error: $error")
          }
        }
        Left(s"HaplotypeCaller failed: $error")
      case Right(result) =>
        // Save logs on success
        logFile.foreach { lf =>
          Using(new PrintWriter(lf)) { writer =>
            writer.println(s"GATK HaplotypeCaller completed successfully")
            writer.println(s"Arguments: ${args.mkString(" ")}")
            writer.println(s"Exit code: ${result.exitCode}")
            writer.println("\n=== STDOUT ===")
            writer.println(result.stdout)
            writer.println("\n=== STDERR ===")
            writer.println(result.stderr)
          }
        }
        onProgress("SNP calling complete.", 1.0, 1.0)
        Right(CallerResult(vcfFile, logFile))
    }
  }

  /**
   * Pass 2: Call all variants in a contig to discover private/novel SNPs.
   * Only emits variants (sites that differ from reference).
   */
  def callPrivateVariants(
    bamPath: String,
    referencePath: String,
    contig: String,
    onProgress: (String, Double, Double) => Unit,
    outputDir: Option[Path] = None,
    outputPrefix: Option[String] = None
  ): Either[String, CallerResult] = {
    // Ensure BAM index exists
    onProgress("Checking BAM index...", 0.0, 1.0)
    GatkRunner.ensureIndex(bamPath) match {
      case Left(error) => return Left(error)
      case Right(_) => // index exists or was created
    }

    onProgress(s"Calling variants in $contig...", 0.1, 1.0)

    // Determine output file location
    val (vcfFile, logFile) = outputDir match {
      case Some(dir) =>
        Files.createDirectories(dir)
        val prefix = outputPrefix.getOrElse(contig)
        (dir.resolve(s"${prefix}_private_variants.vcf").toFile, Some(dir.resolve(s"${prefix}_private_variants.log").toFile))
      case None =>
        val tempFile = File.createTempFile(s"private-variants-$contig", ".vcf")
        tempFile.deleteOnExit()
        (tempFile, None)
    }

    val args = Array(
      "HaplotypeCaller",
      "-I", bamPath,
      "-R", referencePath,
      "-O", vcfFile.getAbsolutePath,
      "-L", contig,
      // Relax reference validation - allows GRCh38 with/without alts, etc.
      "--disable-sequence-dictionary-validation", "true",
      // Haploid calling for mtDNA and Y-DNA
      "--sample-ploidy", "1",
      // Standard confidence for variant discovery
      "--standard-min-confidence-threshold-for-calling", "20.0"
      // Default output-mode is EMIT_VARIANTS_ONLY which is what we want
    )

    GatkRunner.run(args) match {
      case Left(error) =>
        logFile.foreach { lf =>
          Using(new PrintWriter(lf)) { writer =>
            writer.println(s"GATK HaplotypeCaller (private variants) failed for $contig")
            writer.println(s"Arguments: ${args.mkString(" ")}")
            writer.println(s"Error: $error")
          }
        }
        Left(s"Private variant calling failed for $contig: $error")
      case Right(result) =>
        logFile.foreach { lf =>
          Using(new PrintWriter(lf)) { writer =>
            writer.println(s"GATK HaplotypeCaller (private variants) completed for $contig")
            writer.println(s"Arguments: ${args.mkString(" ")}")
            writer.println(s"Exit code: ${result.exitCode}")
            writer.println("\n=== STDOUT ===")
            writer.println(result.stdout)
            writer.println("\n=== STDERR ===")
            writer.println(result.stderr)
          }
        }
        onProgress(s"Variant calling for $contig complete.", 1.0, 1.0)
        Right(CallerResult(vcfFile, logFile))
    }
  }

  /**
   * Two-pass calling: first at tree sites for haplogroup assignment,
   * then full contig for private variant discovery.
   *
   * Automatically selects the appropriate caller:
   * - chrM/MT: Mutect2 --mitochondria (optimized for dense mtDNA positions)
   * - chrY: HaplotypeCaller (better for spread-out Y-DNA positions)
   */
  def callTwoPass(
    bamPath: String,
    referencePath: String,
    allelesVcf: File,
    onProgress: (String, Double, Double) => Unit,
    outputDir: Option[Path] = None,
    outputPrefix: Option[String] = None
  ): Either[String, TwoPassCallerResult] = {
    val contig = detectContigFromVcf(allelesVcf)
    val isMtDna = contig.equalsIgnoreCase("chrM") || contig.equalsIgnoreCase("MT")

    if (isMtDna) {
      callTwoPassMutect2(bamPath, referencePath, allelesVcf, contig, onProgress, outputDir, outputPrefix)
    } else {
      callTwoPassHaplotypeCaller(bamPath, referencePath, allelesVcf, contig, onProgress, outputDir, outputPrefix)
    }
  }

  /**
   * Two-pass calling using HaplotypeCaller (for Y-DNA).
   * Checks for cached VCF files and skips GATK if they exist.
   */
  private def callTwoPassHaplotypeCaller(
    bamPath: String,
    referencePath: String,
    allelesVcf: File,
    contig: String,
    onProgress: (String, Double, Double) => Unit,
    outputDir: Option[Path],
    outputPrefix: Option[String]
  ): Either[String, TwoPassCallerResult] = {
    // Check for cached results
    outputDir match {
      case Some(dir) =>
        val prefix = outputPrefix.getOrElse("called")
        val cachedTreeSites = dir.resolve(s"${prefix}_calls.vcf").toFile
        val cachedPrivateVariants = dir.resolve(s"${prefix}_private_variants.vcf").toFile

        if (cachedTreeSites.exists() && cachedPrivateVariants.exists() &&
            cachedTreeSites.length() > 0 && cachedPrivateVariants.length() > 0) {
          println(s"[GatkHaplotypeCallerProcessor] Using cached VCFs (Y-DNA): ${cachedTreeSites.getName}, ${cachedPrivateVariants.getName}")
          onProgress("Using cached VCF files from previous analysis...", 1.0, 1.0)
          return Right(TwoPassCallerResult(
            treeSitesVcf = cachedTreeSites,
            privateVariantsVcf = cachedPrivateVariants,
            treeSitesLog = Some(dir.resolve(s"${prefix}_haplotypecaller.log").toFile).filter(_.exists()),
            privateVariantsLog = Some(dir.resolve(s"${prefix}_private_variants.log").toFile).filter(_.exists())
          ))
        }
      case None => // No output dir, can't cache
    }

    // Phase 1: Resolve overlapping reference reversed SNPs (tree sites for haplogroup assignment)
    onProgress(s"Phase 1: Resolving reference reversed SNPs...", 0.0, 1.0)
    callSnps(
      bamPath,
      referencePath,
      allelesVcf,
      (msg, done, total) => onProgress(s"Phase 1: $msg", done * 0.4, 1.0),
      outputDir,
      outputPrefix
    ) match {
      case Left(error) => Left(s"Phase 1 failed: $error")
      case Right(treeSitesResult) =>
        // Phase 2: Resolve remaining callable SNPs (private variant discovery)
        onProgress(s"Phase 2: Resolving remaining callable SNPs...", 0.4, 1.0)
        callPrivateVariants(
          bamPath,
          referencePath,
          contig,
          (msg, done, total) => onProgress(s"Phase 2: $msg", 0.4 + done * 0.6, 1.0),
          outputDir,
          outputPrefix
        ) match {
          case Left(error) => Left(s"Phase 2 failed: $error")
          case Right(privateResult) =>
            onProgress("SNP resolution complete.", 1.0, 1.0)
            Right(TwoPassCallerResult(
              treeSitesVcf = treeSitesResult.vcfFile,
              privateVariantsVcf = privateResult.vcfFile,
              treeSitesLog = treeSitesResult.logFile,
              privateVariantsLog = privateResult.logFile
            ))
        }
    }
  }

  // ============================================================================
  // Mutect2 Mitochondria Mode - for mtDNA haplogroup calling
  // ============================================================================

  /**
   * Call SNPs at specified allele sites using Mutect2 mitochondria mode.
   * Much faster than HaplotypeCaller for dense mtDNA positions.
   */
  private def callSnpsMutect2(
    bamPath: String,
    referencePath: String,
    allelesVcf: File,
    onProgress: (String, Double, Double) => Unit,
    outputDir: Option[Path],
    outputPrefix: Option[String]
  ): Either[String, CallerResult] = {
    onProgress("Calling SNPs with Mutect2 (mitochondria mode)...", 0.1, 1.0)

    // Determine output file location
    val (vcfFile, logFile) = outputDir match {
      case Some(dir) =>
        Files.createDirectories(dir)
        val prefix = outputPrefix.getOrElse("mtdna")
        (dir.resolve(s"${prefix}_calls.vcf").toFile, Some(dir.resolve(s"${prefix}_mutect2.log").toFile))
      case None =>
        val tempFile = File.createTempFile("mtdna_calls", ".vcf")
        tempFile.deleteOnExit()
        (tempFile, None)
    }

    // Index the allelesVcf file
    val indexArgs = Array(
      "IndexFeatureFile",
      "-I", allelesVcf.getAbsolutePath
    )

    GatkRunner.run(indexArgs) match {
      case Left(error) => return Left(s"Failed to index alleles VCF: $error")
      case Right(_) => // continue
    }

    // Mutect2 with mitochondria mode - optimized for mtDNA
    // Key differences from HaplotypeCaller:
    // - --mitochondria-mode flag auto-tunes: initial-tumor-lod=0, tumor-lod-to-emit=0,
    //   af-of-alleles-not-in-resource=4e-3, pruning-lod-threshold=-4*ln(10)
    // - Also sets --recover-all-dangling-branches true
    // - Much faster on dense positions due to optimized active region handling
    val args = Array(
      "Mutect2",
      "-I", bamPath,
      "-R", referencePath,
      "-O", vcfFile.getAbsolutePath,
      "-L", allelesVcf.getAbsolutePath,
      "--alleles", allelesVcf.getAbsolutePath,
      "--mitochondria-mode",
      "--disable-sequence-dictionary-validation", "true",
      // Force output at all allele sites
      "--genotype-germline-sites", "true",
      "--force-call-filtered-alleles", "true"
    )

    GatkRunner.run(args) match {
      case Left(error) =>
        logFile.foreach { lf =>
          Using(new PrintWriter(lf)) { writer =>
            writer.println(s"Mutect2 (mitochondria mode) failed")
            writer.println(s"Arguments: ${args.mkString(" ")}")
            writer.println(s"Error: $error")
          }
        }
        Left(s"Mutect2 failed: $error")
      case Right(result) =>
        logFile.foreach { lf =>
          Using(new PrintWriter(lf)) { writer =>
            writer.println(s"Mutect2 (mitochondria mode) completed successfully")
            writer.println(s"Arguments: ${args.mkString(" ")}")
            writer.println(s"Exit code: ${result.exitCode}")
            writer.println("\n=== STDOUT ===")
            writer.println(result.stdout)
            writer.println("\n=== STDERR ===")
            writer.println(result.stderr)
          }
        }
        onProgress("mtDNA SNP calling complete.", 1.0, 1.0)
        Right(CallerResult(vcfFile, logFile))
    }
  }

  /**
   * Call all variants in mtDNA using Mutect2 mitochondria mode for private variant discovery.
   */
  private def callPrivateVariantsMutect2(
    bamPath: String,
    referencePath: String,
    contig: String,
    onProgress: (String, Double, Double) => Unit,
    outputDir: Option[Path],
    outputPrefix: Option[String]
  ): Either[String, CallerResult] = {
    onProgress(s"Calling variants in $contig...", 0.1, 1.0)

    val (vcfFile, logFile) = outputDir match {
      case Some(dir) =>
        Files.createDirectories(dir)
        val prefix = outputPrefix.getOrElse("mtdna")
        (dir.resolve(s"${prefix}_private_variants.vcf").toFile, Some(dir.resolve(s"${prefix}_private_variants.log").toFile))
      case None =>
        val tempFile = File.createTempFile(s"mtdna-private-variants", ".vcf")
        tempFile.deleteOnExit()
        (tempFile, None)
    }

    val args = Array(
      "Mutect2",
      "-I", bamPath,
      "-R", referencePath,
      "-O", vcfFile.getAbsolutePath,
      "-L", contig,
      "--mitochondria-mode",
      "--disable-sequence-dictionary-validation", "true"
    )

    GatkRunner.run(args) match {
      case Left(error) =>
        logFile.foreach { lf =>
          Using(new PrintWriter(lf)) { writer =>
            writer.println(s"Mutect2 (private variants) failed for $contig")
            writer.println(s"Arguments: ${args.mkString(" ")}")
            writer.println(s"Error: $error")
          }
        }
        Left(s"Mutect2 private variant calling failed for $contig: $error")
      case Right(result) =>
        logFile.foreach { lf =>
          Using(new PrintWriter(lf)) { writer =>
            writer.println(s"Mutect2 (private variants) completed for $contig")
            writer.println(s"Arguments: ${args.mkString(" ")}")
            writer.println(s"Exit code: ${result.exitCode}")
            writer.println("\n=== STDOUT ===")
            writer.println(result.stdout)
            writer.println("\n=== STDERR ===")
            writer.println(result.stderr)
          }
        }
        onProgress(s"Variant calling for $contig complete.", 1.0, 1.0)
        Right(CallerResult(vcfFile, logFile))
    }
  }

  /**
   * Two-pass calling using Mutect2 mitochondria mode (for mtDNA).
   * Significantly faster than HaplotypeCaller for dense mtDNA positions.
   * Checks for cached VCF files and skips GATK if they exist.
   */
  private def callTwoPassMutect2(
    bamPath: String,
    referencePath: String,
    allelesVcf: File,
    contig: String,
    onProgress: (String, Double, Double) => Unit,
    outputDir: Option[Path],
    outputPrefix: Option[String]
  ): Either[String, TwoPassCallerResult] = {
    // Check for cached results
    outputDir match {
      case Some(dir) =>
        val prefix = outputPrefix.getOrElse("mtdna")
        val cachedTreeSites = dir.resolve(s"${prefix}_calls.vcf").toFile
        val cachedPrivateVariants = dir.resolve(s"${prefix}_private_variants.vcf").toFile

        if (cachedTreeSites.exists() && cachedPrivateVariants.exists() &&
            cachedTreeSites.length() > 0 && cachedPrivateVariants.length() > 0) {
          println(s"[GatkHaplotypeCallerProcessor] Using cached VCFs (mtDNA): ${cachedTreeSites.getName}, ${cachedPrivateVariants.getName}")
          onProgress("Using cached VCF files from previous analysis...", 1.0, 1.0)
          return Right(TwoPassCallerResult(
            treeSitesVcf = cachedTreeSites,
            privateVariantsVcf = cachedPrivateVariants,
            treeSitesLog = Some(dir.resolve(s"${prefix}_mutect2.log").toFile).filter(_.exists()),
            privateVariantsLog = Some(dir.resolve(s"${prefix}_private_variants.log").toFile).filter(_.exists())
          ))
        }
      case None => // No output dir, can't cache
    }

    // Ensure BAM index exists
    onProgress("Checking BAM index...", 0.0, 1.0)
    GatkRunner.ensureIndex(bamPath) match {
      case Left(error) => return Left(error)
      case Right(_) => // continue
    }

    // Phase 1: Resolve overlapping reference reversed SNPs (mtDNA tree sites)
    onProgress(s"Phase 1: Resolving overlapping reference reversed SNPs...", 0.0, 1.0)
    callSnpsMutect2(
      bamPath,
      referencePath,
      allelesVcf,
      (msg, done, total) => onProgress(s"Phase 1: $msg", done * 0.4, 1.0),
      outputDir,
      outputPrefix
    ) match {
      case Left(error) => Left(s"Phase 1 failed: $error")
      case Right(treeSitesResult) =>
        // Phase 2: Resolve remaining callable SNPs (private variant discovery)
        onProgress(s"Phase 2: Resolving remaining callable SNPs...", 0.4, 1.0)
        callPrivateVariantsMutect2(
          bamPath,
          referencePath,
          contig,
          (msg, done, total) => onProgress(s"Phase 2: $msg", 0.4 + done * 0.6, 1.0),
          outputDir,
          outputPrefix
        ) match {
          case Left(error) => Left(s"Phase 2 failed: $error")
          case Right(privateResult) =>
            onProgress("SNP resolution complete.", 1.0, 1.0)
            Right(TwoPassCallerResult(
              treeSitesVcf = treeSitesResult.vcfFile,
              privateVariantsVcf = privateResult.vcfFile,
              treeSitesLog = treeSitesResult.logFile,
              privateVariantsLog = privateResult.logFile
            ))
        }
    }
  }
}
