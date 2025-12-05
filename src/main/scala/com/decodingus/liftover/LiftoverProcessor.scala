package com.decodingus.liftover

import org.broadinstitute.hellbender.Main

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
      // Relax reference validation - allows GRCh38 with/without alts, etc.
      "--VALIDATION_STRINGENCY", "LENIENT",
      "--disable-sequence-dictionary-validation", "true"
    )

    try {
      Main.main(args)
      onProgress("VCF liftover complete.", 1.0, 1.0)
      Right(liftedVcfFile)
    } catch {
      case e: Exception =>
        e.printStackTrace()
        Left(s"GATK LiftoverVcf failed: ${e.getMessage}")
    }
  }
}
