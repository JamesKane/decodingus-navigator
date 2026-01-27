package com.decodingus.analysis

import com.decodingus.analysis.util.BioVisualizationUtil
import com.decodingus.model.ContigSummary
import htsjdk.samtools.reference.ReferenceSequenceFileFactory

import java.io.File
import java.nio.file.{Files, Path}
import scala.collection.mutable.ListBuffer
import scala.io.Source
import scala.util.{Either, Left, Right, Using, boundary}

case class CallableLociResult(
                               callableBases: Long,
                               contigAnalysis: List[ContigSummary]
                             )

class CallableLociProcessor {

  private val ARTIFACT_SUBDIR_NAME = "callable_loci"

  // Main assembly contigs only - excludes alts, decoys, HLA, etc.
  private val mainAssemblyPattern = "^(chr)?([1-9]|1[0-9]|2[0-2]|X|Y|M|MT)$".r

  private def isMainAssemblyContig(name: String): Boolean = {
    mainAssemblyPattern.findFirstIn(name).isDefined
  }

  /**
   * Process a BAM/CRAM file to compute callable loci.
   *
   * @param bamPath         Path to the BAM/CRAM file
   * @param referencePath   Path to the reference genome
   * @param onProgress      Progress callback
   * @param artifactContext Optional context for organizing output artifacts by subject/run/alignment
   * @param minDepth        Minimum depth to consider a position callable (default 4, use 2 for HiFi)
   */
  def process(
               bamPath: String,
               referencePath: String,
               onProgress: (String, Int, Int) => Unit,
               artifactContext: Option[ArtifactContext] = None,
               minDepth: Int = 4
             ): Either[Throwable, (CallableLociResult, List[String])] = {
    // Ensure BAM index exists
    onProgress("Checking BAM index...", 0, 1)
    GatkRunner.ensureIndex(bamPath) match {
      case Left(error) => return Left(new RuntimeException(error))
      case Right(_) => // index exists or was created
    }

    val referenceFile = new File(referencePath)
    val dictionary = ReferenceSequenceFileFactory.getReferenceSequenceFile(referenceFile).getSequenceDictionary
    val allContigs = dictionary.getSequences.toArray.map(_.asInstanceOf[htsjdk.samtools.SAMSequenceRecord])
    // Filter to main assembly contigs only (chr1-22, X, Y, M/MT)
    val contigs = allContigs.filter(c => isMainAssemblyContig(c.getSequenceName))
    val totalContigs = contigs.length
    val contigLengths = contigs.map(s => s.getSequenceName -> s.getSequenceLength).toMap
    val maxGenomeLength = if (contigLengths.values.isEmpty) 0 else contigLengths.values.max

    // Use artifact cache directory if context provided, otherwise use local directory
    val outputDir: File = artifactContext match {
      case Some(ctx) => ctx.getSubdir(ARTIFACT_SUBDIR_NAME).toFile
      case None =>
        val dir = new File(ARTIFACT_SUBDIR_NAME)
        if (!dir.exists()) Files.createDirectories(dir.toPath)
        dir
    }

    val allSvgStrings = ListBuffer[String]()
    val allContigSummaries = ListBuffer[ContigSummary]()

    boundary[Either[Throwable, (CallableLociResult, List[String])]] {
      for ((contig, index) <- contigs.zipWithIndex) {
        val contigName = contig.getSequenceName
        val contigLength = contig.getSequenceLength

        onProgress(s"Analyzing contig: $contigName (${index + 1} of $totalContigs)", index + 1, totalContigs)

        val bedFile = new File(outputDir, s"$contigName.callable.bed")
        val summaryFile = new File(outputDir, s"$contigName.table.txt")

        val args = Array(
          "CallableLoci",
          "-I", bamPath,
          "-R", referencePath,
          "-O", bedFile.getAbsolutePath,
          "--summary", summaryFile.getAbsolutePath,
          "-L", contigName,
          // Minimum depth to consider callable (default 4, lower for high-accuracy long reads)
          "--min-depth", minDepth.toString,
          // Relax reference validation - allows GRCh38 with/without alts, etc.
          "--disable-sequence-dictionary-validation", "true"
        )

        GatkRunner.run(args) match {
          case Right(_) =>
            val binData = BioVisualizationUtil.binIntervalsFromBed(bedFile.toPath, contigName, contigLength)
            val svgString = BioVisualizationUtil.generateSvgForContig(contigName, contigLength, maxGenomeLength, binData)
            allSvgStrings += svgString

            // Write SVG to file
            val svgFile = new File(outputDir, s"$contigName.callable.svg")
            Files.writeString(svgFile.toPath, svgString)

            val contigSummary = parseSummary(summaryFile.getAbsolutePath, contigName)
            allContigSummaries += contigSummary
          case Left(error) =>
            boundary.break(Left(new RuntimeException(s"GATK CallableLoci failed for contig $contigName: $error")))
        }
      }

      val callableBases = allContigSummaries.map(_.callable).sum

      val result = CallableLociResult(
        callableBases = callableBases,
        contigAnalysis = allContigSummaries.toList
      )

      Right((result, allSvgStrings.toList))
    }
  }

  private def parseSummary(summaryPath: String, contigName: String): ContigSummary = {
    val summaryMap = scala.collection.mutable.Map[String, Long]()
    Using(Source.fromFile(summaryPath)) { source =>
      for (line <- source.getLines()) {
        if (!line.strip.startsWith("state nBases") && line.strip.nonEmpty) {
          val fields = line.strip.split("\\s+")
          if (fields.length == 2) {
            summaryMap(fields(0)) = fields(1).toLong
          }
        }
      }
    }

    ContigSummary(
      contigName = contigName,
      refN = summaryMap.getOrElse("REF_N", 0L),
      callable = summaryMap.getOrElse("CALLABLE", 0L),
      noCoverage = summaryMap.getOrElse("NO_COVERAGE", 0L),
      lowCoverage = summaryMap.getOrElse("LOW_COVERAGE", 0L),
      excessiveCoverage = summaryMap.getOrElse("EXCESSIVE_COVERAGE", 0L),
      poorMappingQuality = summaryMap.getOrElse("POOR_MAPPING_QUALITY", 0L)
    )
  }
}

object CallableLociProcessor {

  /**
   * Load CallableLociResult from cached artifacts.
   * Reads the .table.txt summary files from the callable_loci directory.
   *
   * @param callableLociDir Path to the callable_loci artifact directory
   * @return CallableLociResult if successful, None if not found or invalid
   */
  def loadFromCache(callableLociDir: Path): Option[CallableLociResult] = {
    if (!Files.exists(callableLociDir)) return None

    import scala.jdk.CollectionConverters.*

    val tableFiles = Files.list(callableLociDir).iterator().asScala
      .filter(_.toString.endsWith(".table.txt"))
      .toList

    if (tableFiles.isEmpty) return None

    val contigSummaries = ListBuffer[ContigSummary]()

    for (tableFile <- tableFiles) {
      val fileName = tableFile.getFileName.toString
      val contigName = fileName.stripSuffix(".table.txt")

      Using(Source.fromFile(tableFile.toFile)) { source =>
        val summaryMap = scala.collection.mutable.Map[String, Long]()
        for (line <- source.getLines()) {
          if (!line.strip.startsWith("state nBases") && line.strip.nonEmpty) {
            val fields = line.strip.split("\\s+")
            if (fields.length == 2) {
              summaryMap(fields(0)) = fields(1).toLong
            }
          }
        }

        contigSummaries += ContigSummary(
          contigName = contigName,
          refN = summaryMap.getOrElse("REF_N", 0L),
          callable = summaryMap.getOrElse("CALLABLE", 0L),
          noCoverage = summaryMap.getOrElse("NO_COVERAGE", 0L),
          lowCoverage = summaryMap.getOrElse("LOW_COVERAGE", 0L),
          excessiveCoverage = summaryMap.getOrElse("EXCESSIVE_COVERAGE", 0L),
          poorMappingQuality = summaryMap.getOrElse("POOR_MAPPING_QUALITY", 0L)
        )
      }
    }

    if (contigSummaries.isEmpty) return None

    // Sort contigs by standard order (chr1, chr2, ..., chrX, chrY, chrM)
    val sortedSummaries = contigSummaries.toList.sortBy { cs =>
      val name = cs.contigName.replaceFirst("^chr", "")
      name match {
        case "X" => 23
        case "Y" => 24
        case "M" | "MT" => 25
        case n if n.forall(_.isDigit) => n.toInt
        case _ => 100
      }
    }

    val callableBases = sortedSummaries.map(_.callable).sum

    Some(CallableLociResult(
      callableBases = callableBases,
      contigAnalysis = sortedSummaries
    ))
  }

  /**
   * Load from cache using artifact context IDs.
   */
  def loadFromCache(
                     sampleAccession: String,
                     runId: String,
                     alignmentId: String
                   ): Option[CallableLociResult] = {
    val callableLociDir = SubjectArtifactCache.getArtifactSubdir(
      sampleAccession, runId, alignmentId, "callable_loci"
    )
    loadFromCache(callableLociDir)
  }

  /**
   * Check if callable loci data exists in cache.
   */
  def existsInCache(callableLociDir: Path): Boolean = {
    if (!Files.exists(callableLociDir)) return false
    import scala.jdk.CollectionConverters.*
    Files.list(callableLociDir).iterator().asScala.exists(_.toString.endsWith(".table.txt"))
  }
}