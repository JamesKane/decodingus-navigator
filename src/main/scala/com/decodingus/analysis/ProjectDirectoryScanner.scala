package com.decodingus.analysis

import com.decodingus.util.Logger

import java.io.File

/**
 * A file discovered within a sample directory, classified by role.
 */
enum DiscoveredFileType:
  case Alignment   // .bam, .cram
  case Index       // .bai, .crai, .tbi, .csi
  case Variant     // .vcf, .vcf.gz, .g.vcf.gz, .gvcf.gz
  case Flagstat    // samtools flagstat output
  case WgsMetrics  // GATK/Picard CollectWgsMetrics output
  case InsertSize  // GATK/Picard CollectInsertSizeMetrics output
  case AlignmentSummary // GATK/Picard CollectAlignmentSummaryMetrics output
  case Other       // anything else

/**
 * A single discovered file with its classified type.
 */
case class DiscoveredFile(
  file: File,
  fileType: DiscoveredFileType
)

/**
 * A sample directory containing discovered data and analysis files.
 *
 * @param sampleId       Directory name, typically a sample alias (e.g., HG02759)
 * @param directory      The sample directory
 * @param alignmentFiles BAM/CRAM files
 * @param variantFiles   VCF files
 * @param metricsFiles   Pre-existing analysis output (flagstat, wgs_metrics, etc.)
 * @param indexFiles     Index files (.bai, .crai, .tbi)
 * @param allFiles       Complete list of discovered files with types
 */
case class DiscoveredSample(
  sampleId: String,
  directory: File,
  alignmentFiles: List[File],
  variantFiles: List[File],
  metricsFiles: List[DiscoveredFile],
  indexFiles: List[File],
  allFiles: List[DiscoveredFile]
) {
  def hasAlignments: Boolean = alignmentFiles.nonEmpty
  def hasVariants: Boolean = variantFiles.nonEmpty
  def hasPrecomputedMetrics: Boolean = metricsFiles.nonEmpty

  def flagstatFiles: List[File] = metricsFiles.collect {
    case DiscoveredFile(f, DiscoveredFileType.Flagstat) => f
  }

  def wgsMetricsFiles: List[File] = metricsFiles.collect {
    case DiscoveredFile(f, DiscoveredFileType.WgsMetrics) => f
  }
}

/**
 * A project directory containing discovered samples.
 *
 * @param projectId ENA project accession or directory name (e.g., PRJEB31736)
 * @param directory The project root directory
 * @param samples   Discovered sample subdirectories
 */
case class DiscoveredProject(
  projectId: String,
  directory: File,
  samples: List[DiscoveredSample]
) {
  def sampleCount: Int = samples.size
  def totalAlignmentFiles: Int = samples.map(_.alignmentFiles.size).sum
  def totalVariantFiles: Int = samples.map(_.variantFiles.size).sum
  def samplesWithMetrics: Int = samples.count(_.hasPrecomputedMetrics)
}

/**
 * Scans NAS directory trees following the convention:
 *   {projectRoot}/{sampleId}/files...
 *
 * The scanner expects the user to select a project directory (e.g., /Volumes/nas/Genomics/PRJEB31736).
 * Each immediate subdirectory is treated as a sample. Files within each sample directory are classified
 * by type — alignments, variants, pre-existing analysis metrics, and indices.
 */
object ProjectDirectoryScanner {

  private val log = Logger[ProjectDirectoryScanner.type]

  // Project accession pattern (ENA/SRA/dbGaP)
  private val ProjectAccessionPattern = """^(PRJ[A-Z]{2}\d+|SRP\d+|phs\d+)$""".r

  // File classification patterns
  private val AlignmentExtensions = Set(".bam", ".cram")
  private val IndexExtensions = Set(".bai", ".crai", ".tbi", ".csi")
  private val VariantPatterns = List(".g.vcf.gz", ".gvcf.gz", ".vcf.gz", ".vcf")
  private val FlagstatPatterns = List(".flagstat", "_flagstat", ".flagstat.txt")
  private val WgsMetricsPatterns = List("wgs_metrics", "collect_wgs_metrics", "wgs.metrics",
    "collectwgsmetrics", ".wgs_metrics.txt")
  private val InsertSizePatterns = List("insert_size_metrics", "insert_size", "insertsizemetrics")
  private val AlignSummaryPatterns = List("alignment_summary_metrics", "alignmentsummary",
    "alignment_summary")

  /**
   * Scan a project directory for sample subdirectories and their files.
   *
   * @param projectDir The project root (e.g., /Volumes/nas/Genomics/PRJEB31736)
   * @return Either an error message or the discovered project structure
   */
  def scan(projectDir: File): Either[String, DiscoveredProject] = {
    if (!projectDir.exists())
      return Left(s"Directory does not exist: ${projectDir.getAbsolutePath}")
    if (!projectDir.isDirectory)
      return Left(s"Not a directory: ${projectDir.getAbsolutePath}")

    val projectId = projectDir.getName
    log.info(s"Scanning project directory: $projectId at ${projectDir.getAbsolutePath}")

    val subdirs = projectDir.listFiles()
      .filter(f => f.isDirectory && !f.getName.startsWith("."))
      .sortBy(_.getName)
      .toList

    if (subdirs.isEmpty)
      return Left(s"No sample subdirectories found in: ${projectDir.getAbsolutePath}")

    val samples = subdirs.flatMap(scanSampleDirectory)

    if (samples.isEmpty)
      return Left(s"No samples with data files found in: ${projectDir.getAbsolutePath}")

    log.info(s"Discovered ${samples.size} samples in project $projectId")
    Right(DiscoveredProject(projectId, projectDir, samples))
  }

  /**
   * Check if a directory name looks like an ENA project accession.
   */
  def isProjectAccession(name: String): Boolean =
    ProjectAccessionPattern.findFirstIn(name).isDefined

  /**
   * Scan a single sample subdirectory for data files.
   * Returns None if the directory contains no alignment or variant files.
   */
  private def scanSampleDirectory(dir: File): Option[DiscoveredSample] = {
    val files = listFilesRecursive(dir, maxDepth = 2)
    val classified = files.map(f => DiscoveredFile(f, classifyFile(f)))

    val alignments = classified.collect { case DiscoveredFile(f, DiscoveredFileType.Alignment) => f }
    val variants = classified.collect { case DiscoveredFile(f, DiscoveredFileType.Variant) => f }
    val indices = classified.collect { case DiscoveredFile(f, DiscoveredFileType.Index) => f }
    val metrics = classified.filter(df => df.fileType match {
      case DiscoveredFileType.Flagstat | DiscoveredFileType.WgsMetrics |
           DiscoveredFileType.InsertSize | DiscoveredFileType.AlignmentSummary => true
      case _ => false
    })

    // Only include directories that have at least one alignment or variant file
    if (alignments.isEmpty && variants.isEmpty) {
      log.debug(s"Skipping ${dir.getName}: no alignment or variant files")
      None
    } else {
      Some(DiscoveredSample(
        sampleId = dir.getName,
        directory = dir,
        alignmentFiles = alignments,
        variantFiles = variants,
        metricsFiles = metrics,
        indexFiles = indices,
        allFiles = classified
      ))
    }
  }

  /**
   * Classify a file based on its name and extension.
   */
  private def classifyFile(file: File): DiscoveredFileType = {
    val name = file.getName.toLowerCase

    // Check multi-part extensions first (order matters — .g.vcf.gz before .vcf.gz)
    if (VariantPatterns.exists(name.endsWith)) return DiscoveredFileType.Variant
    if (FlagstatPatterns.exists(p => name.contains(p) || name.endsWith(p))) return DiscoveredFileType.Flagstat
    if (WgsMetricsPatterns.exists(p => name.contains(p))) return DiscoveredFileType.WgsMetrics
    if (InsertSizePatterns.exists(p => name.contains(p))) return DiscoveredFileType.InsertSize
    if (AlignSummaryPatterns.exists(p => name.contains(p))) return DiscoveredFileType.AlignmentSummary

    // Simple extensions
    val ext = "." + name.split("\\.").lastOption.getOrElse("")
    if (AlignmentExtensions.contains(ext)) return DiscoveredFileType.Alignment
    if (IndexExtensions.contains(ext)) return DiscoveredFileType.Index

    DiscoveredFileType.Other
  }

  /**
   * List files recursively up to a max depth.
   * Avoids symlink loops and hidden directories.
   */
  private def listFilesRecursive(dir: File, maxDepth: Int, currentDepth: Int = 0): List[File] = {
    if (currentDepth > maxDepth || !dir.isDirectory) return List.empty

    val entries = Option(dir.listFiles()).getOrElse(Array.empty[File])
    val files = entries.filter(_.isFile).toList
    val subdirFiles = if (currentDepth < maxDepth) {
      entries
        .filter(f => f.isDirectory && !f.getName.startsWith("."))
        .flatMap(d => listFilesRecursive(d, maxDepth, currentDepth + 1))
        .toList
    } else List.empty

    files ++ subdirFiles
  }
}
