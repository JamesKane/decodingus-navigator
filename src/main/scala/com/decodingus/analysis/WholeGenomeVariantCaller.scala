package com.decodingus.analysis

import htsjdk.samtools.reference.ReferenceSequenceFileFactory

import java.io.{File, PrintWriter}
import java.nio.file.{Files, Path}
import scala.jdk.CollectionConverters._
import scala.util.Using

/**
 * Orchestrates whole-genome variant calling for a single alignment.
 *
 * Calling strategy by contig type:
 * - Autosomes (chr1-22): HaplotypeCaller with diploid ploidy (2)
 * - chrX: HaplotypeCaller with ploidy based on inferred sex (1 for male, 2 for female)
 * - chrY: HaplotypeCaller with haploid ploidy (1), skipped for females
 * - chrM: Mutect2 with --mitochondria flag (preserves existing implementation)
 */
class WholeGenomeVariantCaller {

  // Contig patterns
  private val autosomePattern = "^(chr)?([1-9]|1[0-9]|2[0-2])$".r
  private val chrXPattern = "^(chr)?X$".r
  private val chrYPattern = "^(chr)?Y$".r
  private val chrMPattern = "^(chr)?(M|MT)$".r

  // Main assembly contigs only
  private val mainAssemblyPattern = "^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$".r

  private def isMainAssemblyContig(name: String): Boolean =
    mainAssemblyPattern.findFirstIn(name).isDefined

  /**
   * Generate a whole-genome VCF for an alignment.
   *
   * @param bamPath Path to BAM/CRAM file
   * @param referencePath Path to reference genome
   * @param outputDir Directory to write output files
   * @param referenceBuild Reference build name (e.g., "GRCh38")
   * @param onProgress Progress callback (message, current, total)
   * @param sexInferenceResult Optional pre-computed sex inference result (avoids re-running)
   * @return Either error or the generated VCF metadata
   */
  def generateWholeGenomeVcf(
    bamPath: String,
    referencePath: String,
    outputDir: Path,
    referenceBuild: String,
    onProgress: (String, Int, Int) => Unit,
    sexInferenceResult: Option[SexInference.SexInferenceResult] = None
  ): Either[String, CachedVcfInfo] = {
    // Ensure BAM index exists
    onProgress("Checking BAM index...", 0, 1)
    GatkRunner.ensureIndex(bamPath) match {
      case Left(error) => return Left(error)
      case Right(_) =>
    }

    // Get contigs from reference
    onProgress("Reading reference dictionary...", 0, 1)
    val refFile = new File(referencePath)
    val dictionary = ReferenceSequenceFileFactory.getReferenceSequenceFile(refFile).getSequenceDictionary
    val allContigs = dictionary.getSequences.asScala.toList
    val mainContigs = allContigs.filter(c => isMainAssemblyContig(c.getSequenceName))

    if (mainContigs.isEmpty) {
      return Left("No main assembly contigs found in reference")
    }

    // Use provided sex result or infer if not provided
    val sexResult = sexInferenceResult.getOrElse {
      onProgress("Inferring sex from coverage ratios...", 0, 1)
      SexInference.inferFromBam(bamPath, (msg, _) => onProgress(msg, 0, 1)) match {
        case Left(error) =>
          // Continue with unknown sex (will use conservative defaults)
          println(s"[WholeGenomeVariantCaller] Sex inference failed: $error, using default ploidy")
          SexInference.SexInferenceResult(
            SexInference.InferredSex.Unknown,
            xAutosomeRatio = 0.0,
            autosomeMeanCoverage = 0.0,
            xCoverage = 0.0,
            confidence = "low"
          )
        case Right(result) =>
          println(s"[WholeGenomeVariantCaller] Inferred sex: ${result.inferredSex}, confidence: ${result.confidence}")
          result
      }
    }
    println(s"[WholeGenomeVariantCaller] Using sex: ${sexResult.inferredSex}, confidence: ${sexResult.confidence}")

    // Create output directory
    Files.createDirectories(outputDir)

    // Call variants per contig
    val perContigVcfs = scala.collection.mutable.ListBuffer[File]()
    val processedContigs = scala.collection.mutable.ListBuffer[String]()
    val totalContigs = mainContigs.size

    for ((contig, idx) <- mainContigs.zipWithIndex) {
      val contigName = contig.getSequenceName

      // Determine ploidy for this contig
      val ploidyOpt = SexInference.ploidyForContig(contigName, sexResult)

      ploidyOpt match {
        case None =>
          // Skip this contig (e.g., chrY for females)
          onProgress(s"Skipping $contigName (not applicable for ${sexResult.inferredSex})", idx + 1, totalContigs)
        case Some(ploidy) =>
          onProgress(s"Calling variants on $contigName (ploidy=$ploidy)...", idx + 1, totalContigs)

          val vcfResult = if (chrMPattern.findFirstIn(contigName).isDefined) {
            // Use Mutect2 for mitochondria
            callMitochondrialVariants(bamPath, referencePath, contigName, outputDir)
          } else {
            // Use HaplotypeCaller for autosomes/X/Y
            callContigVariants(bamPath, referencePath, contigName, ploidy, outputDir)
          }

          vcfResult match {
            case Left(error) =>
              println(s"[WholeGenomeVariantCaller] Warning: Failed to call $contigName: $error")
              // Continue with other contigs
            case Right(vcfFile) =>
              perContigVcfs += vcfFile
              processedContigs += contigName
          }
      }
    }

    if (perContigVcfs.isEmpty) {
      return Left("No variants called on any contig")
    }

    // Merge per-contig VCFs
    onProgress("Merging per-contig VCFs...", totalContigs, totalContigs)
    val finalVcfPath = outputDir.resolve("whole_genome.vcf.gz")
    val finalIndexPath = outputDir.resolve("whole_genome.vcf.gz.tbi")

    mergeVcfs(perContigVcfs.toList, finalVcfPath, referencePath) match {
      case Left(error) => return Left(s"Failed to merge VCFs: $error")
      case Right(_) =>
    }

    // Count variants
    val variantCount = countVariants(finalVcfPath)

    // Create metadata
    val metadata = VcfCache.createMetadata(
      vcfPath = finalVcfPath,
      indexPath = finalIndexPath,
      referenceBuild = referenceBuild,
      callerVersion = "WholeGenomeVariantCaller/1.0",
      gatkVersion = getGatkVersion,
      contigs = processedContigs.toList,
      variantCount = variantCount,
      inferredSex = Some(sexResult.inferredSex.toString)
    )

    // Save metadata
    val metadataPath = outputDir.resolve("vcf_metadata.json")
    import io.circe.syntax._
    Files.writeString(metadataPath, metadata.asJson.spaces2)

    // Clean up per-contig VCFs
    perContigVcfs.foreach { vcf =>
      try {
        vcf.delete()
        new File(vcf.getAbsolutePath + ".idx").delete()
        new File(vcf.getAbsolutePath + ".tbi").delete()
      } catch {
        case _: Exception => // Ignore cleanup errors
      }
    }

    onProgress("Whole-genome VCF generation complete", totalContigs, totalContigs)
    Right(metadata)
  }

  /**
   * Call variants on a single contig using HaplotypeCaller.
   */
  private def callContigVariants(
    bamPath: String,
    referencePath: String,
    contigName: String,
    ploidy: Int,
    outputDir: Path
  ): Either[String, File] = {
    val vcfFile = outputDir.resolve(s"$contigName.vcf").toFile

    val args = Array(
      "HaplotypeCaller",
      "-I", bamPath,
      "-R", referencePath,
      "-O", vcfFile.getAbsolutePath,
      "-L", contigName,
      "--sample-ploidy", ploidy.toString,
      "--disable-sequence-dictionary-validation", "true",
      // Emit all confident calls
      "--standard-min-confidence-threshold-for-calling", "10.0"
    )

    GatkRunner.run(args) match {
      case Left(error) => Left(error)
      case Right(_) =>
        if (vcfFile.exists() && vcfFile.length() > 0) {
          Right(vcfFile)
        } else {
          Left(s"HaplotypeCaller produced empty or missing VCF for $contigName")
        }
    }
  }

  /**
   * Call variants on mitochondrial DNA using Mutect2.
   * Preserves the existing Mutect2 --mitochondria implementation.
   */
  private def callMitochondrialVariants(
    bamPath: String,
    referencePath: String,
    contigName: String,
    outputDir: Path
  ): Either[String, File] = {
    val vcfFile = outputDir.resolve(s"$contigName.vcf").toFile

    // Mutect2 with mitochondria mode
    val args = Array(
      "Mutect2",
      "-I", bamPath,
      "-R", referencePath,
      "-O", vcfFile.getAbsolutePath,
      "-L", contigName,
      "--mitochondria-mode",
      "--disable-sequence-dictionary-validation", "true"
    )

    GatkRunner.run(args) match {
      case Left(error) => Left(error)
      case Right(_) =>
        if (vcfFile.exists() && vcfFile.length() > 0) {
          Right(vcfFile)
        } else {
          Left(s"Mutect2 produced empty or missing VCF for $contigName")
        }
    }
  }

  /**
   * Merge multiple VCF files into a single gzipped VCF with tabix index.
   */
  private def mergeVcfs(
    vcfFiles: List[File],
    outputPath: Path,
    referencePath: String
  ): Either[String, Unit] = {
    if (vcfFiles.isEmpty) {
      return Left("No VCF files to merge")
    }

    if (vcfFiles.size == 1) {
      // Single VCF - just compress and index
      return compressAndIndex(vcfFiles.head.toPath, outputPath)
    }

    // Use GATK MergeVcfs
    val inputArgs = vcfFiles.flatMap(f => Array("-I", f.getAbsolutePath))
    val args = Array("MergeVcfs") ++ inputArgs ++ Array(
      "-O", outputPath.toString,
      "--CREATE_INDEX", "true"
    )

    GatkRunner.run(args) match {
      case Left(error) => Left(error)
      case Right(_) => Right(())
    }
  }

  /**
   * Compress a VCF file with bgzip and create tabix index.
   */
  private def compressAndIndex(inputVcf: Path, outputPath: Path): Either[String, Unit] = {
    // Use GATK to sort and compress
    val args = Array(
      "SortVcf",
      "-I", inputVcf.toString,
      "-O", outputPath.toString,
      "--CREATE_INDEX", "true"
    )

    GatkRunner.run(args) match {
      case Left(error) => Left(error)
      case Right(_) => Right(())
    }
  }

  /**
   * Count variants in a VCF file.
   */
  private def countVariants(vcfPath: Path): Long = {
    try {
      import htsjdk.variant.vcf.VCFFileReader
      val reader = new VCFFileReader(vcfPath, false)
      try {
        var count = 0L
        val iter = reader.iterator()
        while (iter.hasNext) {
          iter.next()
          count += 1
        }
        count
      } finally {
        reader.close()
      }
    } catch {
      case _: Exception => 0L
    }
  }

  /**
   * Get GATK version string.
   */
  private def getGatkVersion: String = {
    try {
      val pkg = classOf[org.broadinstitute.hellbender.Main].getPackage
      Option(pkg.getImplementationVersion).getOrElse("4.x")
    } catch {
      case _: Exception => "unknown"
    }
  }
}

object WholeGenomeVariantCaller {
  private lazy val instance = new WholeGenomeVariantCaller()
  def apply(): WholeGenomeVariantCaller = instance
}
