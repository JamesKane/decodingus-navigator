package com.decodingus.liftover

import com.decodingus.analysis.GatkRunner
import htsjdk.samtools.reference.ReferenceSequenceFileFactory

import java.io.{File, PrintWriter}
import java.nio.file.Path
import scala.io.Source
import scala.jdk.CollectionConverters.*
import scala.util.Using

class LiftoverProcessor {

  // Mapping from chr-prefixed to non-prefixed contig names
  private val chrToNcbi: Map[String, String] = Map(
    "chrY" -> "Y",
    "chrX" -> "X",
    "chrM" -> "MT"
  ) ++ (1 to 22).map(i => s"chr$i" -> i.toString).toMap

  // Reverse mapping
  private val ncbiToChr: Map[String, String] = chrToNcbi.map(_.swap)

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
        // Normalize contig names to match target reference (handles UCSC vs NCBI naming)
        onProgress("Normalizing contig names...", 0.8, 1.0)
        val normalizedVcf = normalizeContigNames(liftedVcfFile, targetReference)

        // If filtering requested, filter the output VCF to only include the expected contig
        val finalVcf = filterToContig match {
          case Some(contig) =>
            onProgress(s"Filtering to $contig...", 0.9, 1.0)
            filterVcfToContig(normalizedVcf, contig)
          case None =>
            normalizedVcf
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
   *
   * Handles both UCSC (chrY) and NCBI (Y) naming - will match either style.
   */
  private def filterVcfToContig(vcfFile: File, contig: String): File = {
    val filteredVcf = File.createTempFile("filtered_liftover", ".vcf")
    filteredVcf.deleteOnExit()

    // Build set of acceptable contig names (both UCSC and NCBI forms)
    val acceptableContigs: Set[String] = {
      val normalized = contig.toLowerCase
      val base = Set(contig)
      // Add both chr-prefixed and non-prefixed variants
      if (normalized.startsWith("chr")) {
        base + chrToNcbi.getOrElse(contig, contig)
      } else {
        base + ncbiToChr.getOrElse(contig, contig)
      }
    }.map(_.toLowerCase)

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
          // Data line - check if contig matches (either naming convention)
          val lineContig = line.split("\t").headOption.getOrElse("").toLowerCase
          if (acceptableContigs.contains(lineContig)) {
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

  /**
   * Normalize contig names in a VCF to match the target reference.
   * Handles UCSC (chr1, chrY) vs NCBI (1, Y) naming conventions.
   *
   * This is needed because UCSC chain files use chr-prefixed names for both
   * hg38 and hg19, but some references (like hs37d5) use NCBI naming.
   */
  private def normalizeContigNames(vcfFile: File, targetReference: Path): File = {
    // Detect target reference naming style by checking for common contigs
    val targetContigs = getTargetContigNames(targetReference)
    if (targetContigs.isEmpty) {
      println(s"[LiftoverProcessor] Could not read target reference contigs, skipping normalization")
      return vcfFile
    }

    // Determine if target uses chr prefix
    val targetUsesChrPrefix = targetContigs.exists(_.startsWith("chr"))

    // Build the renaming map based on target style
    val renameMap: Map[String, String] = if (targetUsesChrPrefix) {
      // Target uses chr prefix (UCSC style) - rename NCBI to UCSC
      ncbiToChr.filter { case (ncbi, chr) => targetContigs.contains(chr) }
    } else {
      // Target uses no prefix (NCBI style) - rename UCSC to NCBI
      chrToNcbi.filter { case (chr, ncbi) => targetContigs.contains(ncbi) }
    }

    if (renameMap.isEmpty) {
      // No renaming needed
      return vcfFile
    }

    println(s"[LiftoverProcessor] Normalizing contig names (target uses ${if (targetUsesChrPrefix) "chr prefix" else "NCBI style"})")

    val normalizedVcf = File.createTempFile("normalized_liftover", ".vcf")
    normalizedVcf.deleteOnExit()

    var renamedCount = 0

    Using.resources(
      Source.fromFile(vcfFile),
      new PrintWriter(normalizedVcf)
    ) { (source, writer) =>
      for (line <- source.getLines()) {
        if (line.startsWith("##contig=")) {
          // Update contig header lines
          val updatedLine = renameMap.foldLeft(line) { case (l, (from, to)) =>
            l.replace(s"ID=$from,", s"ID=$to,").replace(s"ID=$from>", s"ID=$to>")
          }
          writer.println(updatedLine)
        } else if (line.startsWith("#")) {
          // Keep other header lines unchanged
          writer.println(line)
        } else {
          // Data line - rename contig in first column
          val fields = line.split("\t", 2)
          if (fields.length >= 2) {
            val contig = fields(0)
            val rest = fields(1)
            renameMap.get(contig) match {
              case Some(newContig) =>
                writer.println(s"$newContig\t$rest")
                renamedCount += 1
              case None =>
                writer.println(line)
            }
          } else {
            writer.println(line)
          }
        }
      }
    }

    if (renamedCount > 0) {
      println(s"[LiftoverProcessor] Renamed contigs in $renamedCount variants")
    }

    normalizedVcf
  }

  /**
   * Get the contig names from a reference sequence dictionary.
   */
  private def getTargetContigNames(referencePath: Path): Set[String] = {
    try {
      val refFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(referencePath)
      try {
        val dict = refFile.getSequenceDictionary
        if (dict != null) {
          dict.getSequences.asScala.map(_.getSequenceName).toSet
        } else {
          Set.empty
        }
      } finally {
        refFile.close()
      }
    } catch {
      case e: Exception =>
        println(s"[LiftoverProcessor] Failed to read reference dictionary: ${e.getMessage}")
        Set.empty
    }
  }
}
