package com.decodingus.liftover

import com.decodingus.analysis.GatkRunner

import java.io.File
import java.nio.file.Path

class LiftoverProcessor {

  def liftoverVcf(
                   vcfFile: File,
                   chainFile: Path,
                   targetReference: Path,
                   onProgress: (String, Double, Double) => Unit
                 ): Either[String, File] = {
    onProgress("Performing VCF liftover...", 0.0, 1.0)

    val liftedVcfFile = File.createTempFile("lifted_alleles", ".vcf")
    liftedVcfFile.deleteOnExit()
    val rejectFile = File.createTempFile("rejected_liftover", ".vcf")
    rejectFile.deleteOnExit()

    val args = Array(
      "LiftoverVcf",
      "-I", vcfFile.getAbsolutePath,
      "-O", liftedVcfFile.getAbsolutePath,
      "-C", chainFile.toString,
      "-R", targetReference.toString,
      "--REJECT", rejectFile.getAbsolutePath,
      // Relax validation - allows minor reference mismatches
      "--VALIDATION_STRINGENCY", "SILENT",
      "--WARN_ON_MISSING_CONTIG", "true"
    )

    GatkRunner.run(args) match {
      case Right(_) =>
        onProgress("VCF liftover complete.", 1.0, 1.0)
        Right(liftedVcfFile)
      case Left(error) =>
        Left(s"LiftoverVcf failed: $error")
    }
  }
}
