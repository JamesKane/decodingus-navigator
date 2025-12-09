package com.decodingus.liftover

import com.decodingus.analysis.GatkRunner

import java.io.{File, PrintWriter}
import java.nio.file.Path
import scala.io.Source
import scala.util.Using

class LiftoverProcessor {

  /**
   * Liftover a VCF file to a new reference build.
   *
   * @param vcfFile Input VCF file
   * @param chainFile Chain file for coordinate conversion
   * @param targetReference Target reference genome
   * @param onProgress Progress callback
   * @param filterToContig Optional contig to filter results to. Useful for reverse liftover
   *                       where some positions may map to unexpected contigs (e.g., chrY -> chrX in PAR regions).
   */
  def liftoverVcf(
                   vcfFile: File,
                   chainFile: Path,
                   targetReference: Path,
                   onProgress: (String, Double, Double) => Unit,
                   filterToContig: Option[String] = None
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
      "--WARN_ON_MISSING_CONTIG", "true",
      // Recover SNPs where REF/ALT are swapped due to reverse-complement mapping
      // This is critical for chrY where GRCh38 and CHM13v2 have large inverted regions
      "--RECOVER_SWAPPED_REF_ALT", "true"
    )

    GatkRunner.run(args) match {
      case Right(_) =>
        // If filtering requested, filter the output VCF to only include the expected contig
        val finalVcf = filterToContig match {
          case Some(contig) =>
            onProgress(s"Filtering to $contig...", 0.9, 1.0)
            filterVcfToContig(liftedVcfFile, contig)
          case None =>
            liftedVcfFile
        }
        onProgress("VCF liftover complete.", 1.0, 1.0)
        Right(finalVcf)
      case Left(error) =>
        Left(s"LiftoverVcf failed: $error")
    }
  }

  /**
   * Filter a VCF file to only include variants on a specific contig.
   * This is used after reverse liftover to remove variants that mapped to unexpected contigs
   * (e.g., Y-DNA positions that lifted to chrX due to PAR regions or assembly differences).
   */
  private def filterVcfToContig(vcfFile: File, contig: String): File = {
    val filteredVcf = File.createTempFile("filtered_liftover", ".vcf")
    filteredVcf.deleteOnExit()

    var keptCount = 0
    var filteredCount = 0

    Using.resources(
      Source.fromFile(vcfFile),
      new PrintWriter(filteredVcf)
    ) { (source, writer) =>
      for (line <- source.getLines()) {
        if (line.startsWith("#")) {
          // Keep all header lines
          writer.println(line)
        } else {
          // Data line - check if contig matches
          val lineContig = line.split("\t").headOption.getOrElse("")
          if (lineContig.equalsIgnoreCase(contig)) {
            writer.println(line)
            keptCount += 1
          } else {
            filteredCount += 1
          }
        }
      }
    }

    if (filteredCount > 0) {
      println(s"[LiftoverProcessor] Filtered $filteredCount variants not on $contig (kept $keptCount)")
    }

    filteredVcf
  }
}
