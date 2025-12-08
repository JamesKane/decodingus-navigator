package com.decodingus.haplogroup.caller

import com.decodingus.analysis.GatkRunner
import java.io.{File, PrintWriter}
import java.nio.file.{Files, Path}
import scala.util.Using

/**
 * Result of SNP calling including output file and GATK logs.
 */
case class CallerResult(vcfFile: File, logFile: Option[File])

class GatkHaplotypeCallerProcessor {

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
      "-L", allelesVcf.getAbsolutePath,
      "--alleles", allelesVcf.getAbsolutePath,
      // Relax reference validation - allows GRCh38 with/without alts, etc.
      "--disable-sequence-dictionary-validation", "true"
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

  def callAllVariantsInContig(
    bamPath: String,
    referencePath: String,
    contig: String,
    onProgress: (String, Double, Double) => Unit,
    outputDir: Option[Path] = None
  ): Either[String, CallerResult] = {
    // Ensure BAM index exists
    onProgress("Checking BAM index...", 0.0, 1.0)
    GatkRunner.ensureIndex(bamPath) match {
      case Left(error) => return Left(error)
      case Right(_) => // index exists or was created
    }

    onProgress(s"Calling all variants in $contig with GATK HaplotypeCaller...", 0.1, 1.0)

    // Determine output file location
    val (vcfFile, logFile) = outputDir match {
      case Some(dir) =>
        Files.createDirectories(dir)
        (dir.resolve(s"${contig}_variants.vcf").toFile, Some(dir.resolve(s"${contig}_haplotypecaller.log").toFile))
      case None =>
        val tempFile = File.createTempFile(s"variants-$contig", ".vcf")
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
      "--disable-sequence-dictionary-validation", "true"
    )

    GatkRunner.run(args) match {
      case Left(error) =>
        logFile.foreach { lf =>
          Using(new PrintWriter(lf)) { writer =>
            writer.println(s"GATK HaplotypeCaller failed for $contig")
            writer.println(s"Arguments: ${args.mkString(" ")}")
            writer.println(s"Error: $error")
          }
        }
        Left(s"HaplotypeCaller failed for $contig: $error")
      case Right(result) =>
        logFile.foreach { lf =>
          Using(new PrintWriter(lf)) { writer =>
            writer.println(s"GATK HaplotypeCaller completed successfully for $contig")
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
}
