package com.decodingus.haplogroup.caller

import com.decodingus.analysis.GatkRunner
import java.io.File

class GatkHaplotypeCallerProcessor {

  def callSnps(
    bamPath: String,
    referencePath: String,
    allelesVcf: File,
    onProgress: (String, Double, Double) => Unit
  ): Either[String, File] = {
    onProgress("Calling SNPs with GATK HaplotypeCaller...", 0.0, 1.0)

    val vcfFile = File.createTempFile("haplotypes", ".vcf")
    vcfFile.deleteOnExit()

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
      case Left(error) => Left(s"HaplotypeCaller failed: $error")
      case Right(_) =>
        onProgress("SNP calling complete.", 1.0, 1.0)
        Right(vcfFile)
    }
  }

  def callAllVariantsInContig(
    bamPath: String,
    referencePath: String,
    contig: String,
    onProgress: (String, Double, Double) => Unit
  ): Either[String, File] = {
    onProgress(s"Calling all variants in $contig with GATK HaplotypeCaller...", 0.0, 1.0)

    val vcfFile = File.createTempFile(s"variants-$contig", ".vcf")
    vcfFile.deleteOnExit()

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
      case Left(error) => Left(s"HaplotypeCaller failed for $contig: $error")
      case Right(_) =>
        onProgress(s"Variant calling for $contig complete.", 1.0, 1.0)
        Right(vcfFile)
    }
  }
}
