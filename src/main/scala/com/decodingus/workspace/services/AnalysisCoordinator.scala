package com.decodingus.workspace.services

import com.decodingus.analysis.*
import com.decodingus.config.{FeatureToggles, ReferenceConfigService, UserPreferencesService}
import com.decodingus.haplogroup.model.HaplogroupResult as AnalysisHaplogroupResult
import com.decodingus.haplogroup.tree.{TreeProviderType, TreeType}
import com.decodingus.model.{LibraryStats, WgsMetrics}
import com.decodingus.refgenome.{ReferenceGateway, ReferenceResolveResult}
import com.decodingus.workspace.WorkspaceState
import com.decodingus.workspace.model.*
import htsjdk.samtools.SamReaderFactory

import java.io.File
import scala.concurrent.{ExecutionContext, Future}

/**
 * Progress state for analysis operations.
 */
case class AnalysisProgress(
  message: String,
  percent: Double,
  isComplete: Boolean = false,
  error: Option[String] = None
)

/**
 * Coordinates all analysis operations (library stats, WGS metrics, haplogroups, callable loci).
 *
 * This service runs analyses and returns results along with workspace state updates.
 * The caller is responsible for persisting changes and updating UI state.
 */
class AnalysisCoordinator(implicit ec: ExecutionContext) {

  private val workspaceOps = new WorkspaceOperations()

  // --- Initial Analysis (Library Stats) ---

  /**
   * Runs library stats analysis on a BAM/CRAM file.
   * Returns the stats and updates needed to the workspace.
   */
  def runLibraryStatsAnalysis(
    state: WorkspaceState,
    sampleAccession: String,
    sequenceRunIndex: Int,
    onProgress: AnalysisProgress => Unit
  ): Future[Either[String, (WorkspaceState, LibraryStats, Alignment)]] = Future {
    workspaceOps.findSubject(state, sampleAccession) match {
      case None =>
        Left(s"Subject not found: $sampleAccession")

      case Some(subject) =>
        val sequenceRuns = state.workspace.main.getSequenceRunsForBiosample(subject)
        sequenceRuns.lift(sequenceRunIndex) match {
          case None =>
            Left(s"Sequence run not found at index $sequenceRunIndex")

          case Some(seqRun) =>
            seqRun.files.headOption match {
              case None =>
                Left("No alignment file associated with this sequence run")

              case Some(fileInfo) =>
                val bamPath = fileInfo.location.getOrElse("")
                runLibraryStatsInternal(state, subject, seqRun, bamPath, onProgress)
            }
        }
    }
  }

  private def runLibraryStatsInternal(
    state: WorkspaceState,
    subject: Biosample,
    seqRun: SequenceRun,
    bamPath: String,
    onProgress: AnalysisProgress => Unit
  ): Either[String, (WorkspaceState, LibraryStats, Alignment)] = {
    try {
      // Step 1: Detect reference build
      onProgress(AnalysisProgress("Reading BAM/CRAM header...", 0.1))
      val header = SamReaderFactory.makeDefault().open(new File(bamPath)).getFileHeader
      val libraryStatsProcessor = new LibraryStatsProcessor()
      val referenceBuild = libraryStatsProcessor.detectReferenceBuild(header)

      if (referenceBuild == "Unknown") {
        return Left("Could not determine reference build from BAM/CRAM header.")
      }

      // Step 2: Resolve reference genome
      onProgress(AnalysisProgress(s"Resolving reference: $referenceBuild", 0.2))
      val referenceGateway = new ReferenceGateway((done, total) => {
        val pct = if (total > 0) 0.2 + (done.toDouble / total) * 0.3 else 0.2
        onProgress(AnalysisProgress(s"Downloading reference: ${done / 1024 / 1024}MB", pct))
      })

      val referencePath = referenceGateway.resolve(referenceBuild) match {
        case Right(path) => path.toString
        case Left(error) => return Left(s"Failed to resolve reference: $error")
      }

      // Step 3: Collect library stats
      onProgress(AnalysisProgress("Analyzing library statistics...", 0.5))
      val libraryStats = libraryStatsProcessor.process(bamPath, referencePath, (message, current, total) => {
        val pct = 0.5 + (current.toDouble / total) * 0.4
        onProgress(AnalysisProgress(s"Library Stats: $message", pct))
      })

      // Step 4: Create or update alignment
      onProgress(AnalysisProgress("Saving results...", 0.95))

      val alignUri = seqRun.alignmentRefs.headOption.getOrElse(
        s"local:alignment:${subject.sampleAccession}:${java.util.UUID.randomUUID().toString.take(8)}"
      )
      val existingAlignment = state.workspace.main.alignments.find(_.atUri.contains(alignUri))

      val newAlignment = Alignment(
        atUri = Some(alignUri),
        meta = existingAlignment.map(_.meta.updated("analysis")).getOrElse(RecordMeta.initial),
        sequenceRunRef = seqRun.atUri.getOrElse(""),
        biosampleRef = Some(subject.atUri.getOrElse(s"local:biosample:${subject.sampleAccession}")),
        referenceBuild = libraryStats.referenceBuild,
        aligner = libraryStats.aligner,
        files = seqRun.files,
        metrics = existingAlignment.flatMap(_.metrics)
      )

      // Update sequence run with inferred metadata
      val updatedSeqRun = seqRun.copy(
        meta = seqRun.meta.updated("analysis"),
        platformName = if (seqRun.platformName == "Unknown" || seqRun.platformName == "Other") libraryStats.inferredPlatform else seqRun.platformName,
        instrumentModel = seqRun.instrumentModel.orElse(Some(libraryStats.mostFrequentInstrument)),
        testType = inferTestType(libraryStats),
        libraryLayout = Some(if (libraryStats.pairedReads > libraryStats.readCount / 2) "Paired-End" else "Single-End"),
        totalReads = Some(libraryStats.readCount.toLong),
        readLength = calculateMeanReadLength(libraryStats.lengthDistribution).orElse(seqRun.readLength),
        maxReadLength = libraryStats.lengthDistribution.keys.maxOption.orElse(seqRun.maxReadLength),
        meanInsertSize = calculateMeanInsertSize(libraryStats.insertSizeDistribution).orElse(seqRun.meanInsertSize),
        alignmentRefs = if (seqRun.alignmentRefs.contains(alignUri)) seqRun.alignmentRefs else seqRun.alignmentRefs :+ alignUri
      )

      // Update workspace state
      val updatedSequenceRuns = state.workspace.main.sequenceRuns.map { sr =>
        if (sr.atUri == seqRun.atUri) updatedSeqRun else sr
      }
      val updatedAlignments = if (existingAlignment.isDefined) {
        state.workspace.main.alignments.map { a =>
          if (a.atUri.contains(alignUri)) newAlignment else a
        }
      } else {
        state.workspace.main.alignments :+ newAlignment
      }
      val updatedContent = state.workspace.main.copy(
        sequenceRuns = updatedSequenceRuns,
        alignments = updatedAlignments
      )
      val newState = state.copy(workspace = state.workspace.copy(main = updatedContent))

      onProgress(AnalysisProgress("Analysis complete", 1.0, isComplete = true))
      Right((newState, libraryStats, newAlignment))

    } catch {
      case e: Exception =>
        Left(e.getMessage)
    }
  }

  // --- WGS Metrics Analysis ---

  def runWgsMetricsAnalysis(
    state: WorkspaceState,
    sampleAccession: String,
    sequenceRunIndex: Int,
    onProgress: AnalysisProgress => Unit
  ): Future[Either[String, (WorkspaceState, WgsMetrics)]] = Future {
    workspaceOps.findSubject(state, sampleAccession) match {
      case None =>
        Left(s"Subject not found: $sampleAccession")

      case Some(subject) =>
        val sequenceRuns = state.workspace.main.getSequenceRunsForBiosample(subject)
        sequenceRuns.lift(sequenceRunIndex) match {
          case None =>
            Left(s"Sequence run not found at index $sequenceRunIndex")

          case Some(seqRun) =>
            val alignments = state.workspace.main.getAlignmentsForSequenceRun(seqRun)
            alignments.headOption match {
              case None =>
                Left("Please run initial analysis first to detect reference build")

              case Some(alignment) =>
                seqRun.files.headOption match {
                  case None =>
                    Left("No alignment file associated with this sequence run")

                  case Some(fileInfo) =>
                    val bamPath = fileInfo.location.getOrElse("")
                    runWgsMetricsInternal(state, subject, seqRun, alignment, bamPath, onProgress)
                }
            }
        }
    }
  }

  private def runWgsMetricsInternal(
    state: WorkspaceState,
    subject: Biosample,
    seqRun: SequenceRun,
    alignment: Alignment,
    bamPath: String,
    onProgress: AnalysisProgress => Unit
  ): Either[String, (WorkspaceState, WgsMetrics)] = {
    try {
      onProgress(AnalysisProgress("Resolving reference genome...", 0.1))

      val referenceGateway = new ReferenceGateway((done, total) => {
        val pct = if (total > 0) 0.1 + (done.toDouble / total) * 0.2 else 0.1
        onProgress(AnalysisProgress(s"Downloading reference: ${done / 1024 / 1024}MB", pct))
      })

      val referencePath = referenceGateway.resolve(alignment.referenceBuild) match {
        case Right(path) => path.toString
        case Left(error) => return Left(s"Failed to resolve reference: $error")
      }

      onProgress(AnalysisProgress("Running WGS metrics analysis...", 0.3))

      val artifactCtx = ArtifactContext(
        sampleAccession = subject.sampleAccession,
        sequenceRunUri = seqRun.atUri,
        alignmentUri = alignment.atUri
      )

      // Use WgsMetricsProcessor for standard analysis
      // Pass maxReadLength for long-read data (PacBio HiFi, Nanopore) - GATK defaults to 150bp
      // Enable COUNT_UNPAIRED for single-end libraries (both short and long read)
      val processor = new WgsMetricsProcessor()
      val isSingleEnd = seqRun.libraryLayout.exists(_.equalsIgnoreCase("Single-End"))
      val wgsMetricsResult = processor.process(
        bamPath = bamPath,
        referencePath = referencePath,
        onProgress = (message, current, total) => {
          val pct = 0.3 + (current.toDouble / total) * 0.6
          onProgress(AnalysisProgress(s"Coverage: $message", pct))
        },
        readLength = seqRun.maxReadLength,
        artifactContext = Some(artifactCtx),
        totalReads = seqRun.totalReads,
        countUnpaired = isSingleEnd
      )

      wgsMetricsResult match {
        case Left(error) =>
          Left(s"WGS metrics failed: ${error.getMessage}")

        case Right(wgsMetrics) =>
          onProgress(AnalysisProgress("Saving results...", 0.95))

          // Update alignment with metrics
          val alignmentMetrics = AlignmentMetrics(
            genomeTerritory = Some(wgsMetrics.genomeTerritory),
            meanCoverage = Some(wgsMetrics.meanCoverage),
            sdCoverage = Some(wgsMetrics.sdCoverage),
            medianCoverage = Some(wgsMetrics.medianCoverage),
            pct10x = Some(wgsMetrics.pct10x),
            pct20x = Some(wgsMetrics.pct20x),
            pct30x = Some(wgsMetrics.pct30x),
            callableBases = None,
            contigs = List.empty // Per-contig coverage not available in WgsMetrics
          )

          val updatedAlignment = alignment.copy(
            metrics = Some(alignmentMetrics),
            meta = alignment.meta.updated("metrics")
          )
          val newState = workspaceOps.updateAlignment(state, updatedAlignment)

          onProgress(AnalysisProgress("WGS metrics complete", 1.0, isComplete = true))
          Right((newState, wgsMetrics))
      }

    } catch {
      case e: Exception =>
        Left(e.getMessage)
    }
  }

  // --- Haplogroup Analysis ---

  def runHaplogroupAnalysis(
    state: WorkspaceState,
    sampleAccession: String,
    sequenceRunIndex: Int,
    treeType: TreeType,
    onProgress: AnalysisProgress => Unit
  ): Future[Either[String, (WorkspaceState, AnalysisHaplogroupResult)]] = Future {
    workspaceOps.findSubject(state, sampleAccession) match {
      case None =>
        Left(s"Subject not found: $sampleAccession")

      case Some(subject) =>
        val sequenceRuns = state.workspace.main.getSequenceRunsForBiosample(subject)
        sequenceRuns.lift(sequenceRunIndex) match {
          case None =>
            Left(s"Sequence run not found at index $sequenceRunIndex")

          case Some(seqRun) =>
            val alignments = state.workspace.main.getAlignmentsForSequenceRun(seqRun)
            alignments.headOption match {
              case None =>
                Left("Please run initial analysis first to detect reference build")

              case Some(alignment) =>
                seqRun.files.headOption match {
                  case None =>
                    Left("No alignment file associated with this sequence run")

                  case Some(fileInfo) =>
                    val bamPath = fileInfo.location.getOrElse("")
                    runHaplogroupInternal(state, subject, seqRun, alignment, bamPath, treeType, onProgress)
                }
            }
        }
    }
  }

  private def runHaplogroupInternal(
    state: WorkspaceState,
    subject: Biosample,
    seqRun: SequenceRun,
    alignment: Alignment,
    bamPath: String,
    treeType: TreeType,
    onProgress: AnalysisProgress => Unit
  ): Either[String, (WorkspaceState, AnalysisHaplogroupResult)] = {
    try {
      // Build LibraryStats from existing data
      val libraryStats = LibraryStats(
        readCount = seqRun.totalReads.map(_.toInt).getOrElse(0),
        pairedReads = 0,
        lengthDistribution = Map.empty,
        insertSizeDistribution = Map.empty,
        aligner = alignment.aligner,
        referenceBuild = alignment.referenceBuild,
        sampleName = subject.donorIdentifier,
        flowCells = Map.empty,
        instruments = Map.empty,
        mostFrequentInstrument = seqRun.instrumentModel.getOrElse("Unknown"),
        inferredPlatform = seqRun.platformName,
        platformCounts = Map.empty
      )

      onProgress(AnalysisProgress("Loading haplogroup tree...", 0.1))

      val processor = new HaplogroupProcessor()
      val artifactCtx = ArtifactContext(
        sampleAccession = subject.sampleAccession,
        sequenceRunUri = seqRun.atUri,
        alignmentUri = alignment.atUri
      )

      // Select tree provider based on user preferences
      val treeProviderType = treeType match {
        case TreeType.YDNA =>
          if (UserPreferencesService.getYdnaTreeProvider.equalsIgnoreCase("decodingus"))
            TreeProviderType.DECODINGUS
          else TreeProviderType.FTDNA
        case TreeType.MTDNA =>
          if (UserPreferencesService.getMtdnaTreeProvider.equalsIgnoreCase("decodingus"))
            TreeProviderType.DECODINGUS
          else TreeProviderType.FTDNA
      }

      val result = processor.analyze(
        bamPath,
        libraryStats,
        treeType,
        treeProviderType,
        (message, current, total) => {
          val pct = if (total > 0) current / total else 0.0
          onProgress(AnalysisProgress(message, pct))
        },
        Some(artifactCtx)
      )

      result match {
        case Right(results) if results.nonEmpty =>
          val topResult = results.head

          // Convert to workspace model and update subject
          val workspaceResult = HaplogroupResult(
            haplogroupName = topResult.name,
            score = topResult.score,
            matchingSnps = Some(topResult.matchingSnps),
            mismatchingSnps = Some(topResult.mismatchingSnps),
            ancestralMatches = Some(topResult.ancestralMatches),
            treeDepth = Some(topResult.depth),
            lineagePath = None
          )

          val currentAssignments = subject.haplogroups.getOrElse(HaplogroupAssignments(None, None))
          val updatedAssignments = treeType match {
            case TreeType.YDNA => currentAssignments.copy(yDna = Some(workspaceResult))
            case TreeType.MTDNA => currentAssignments.copy(mtDna = Some(workspaceResult))
          }
          val updatedSubject = subject.copy(
            haplogroups = Some(updatedAssignments),
            meta = subject.meta.updated("haplogroups")
          )
          val newState = workspaceOps.updateSubjectDirect(state, updatedSubject)

          onProgress(AnalysisProgress("Haplogroup analysis complete", 1.0, isComplete = true))
          Right((newState, topResult))

        case Right(_) =>
          Left("No haplogroup matches found")

        case Left(error) =>
          Left(error)
      }

    } catch {
      case e: Exception =>
        Left(e.getMessage)
    }
  }

  // --- Callable Loci Analysis ---

  def runCallableLociAnalysis(
    state: WorkspaceState,
    sampleAccession: String,
    sequenceRunIndex: Int,
    onProgress: AnalysisProgress => Unit
  ): Future[Either[String, (WorkspaceState, CallableLociResult)]] = Future {
    workspaceOps.findSubject(state, sampleAccession) match {
      case None =>
        Left(s"Subject not found: $sampleAccession")

      case Some(subject) =>
        val sequenceRuns = state.workspace.main.getSequenceRunsForBiosample(subject)
        sequenceRuns.lift(sequenceRunIndex) match {
          case None =>
            Left(s"Sequence run not found at index $sequenceRunIndex")

          case Some(seqRun) =>
            val alignments = state.workspace.main.getAlignmentsForSequenceRun(seqRun)
            alignments.headOption match {
              case None =>
                Left("Please run initial analysis first to detect reference build")

              case Some(alignment) =>
                seqRun.files.headOption match {
                  case None =>
                    Left("No alignment file associated with this sequence run")

                  case Some(fileInfo) =>
                    val bamPath = fileInfo.location.getOrElse("")
                    runCallableLociInternal(state, subject, seqRun, alignment, bamPath, onProgress)
                }
            }
        }
    }
  }

  private def runCallableLociInternal(
    state: WorkspaceState,
    subject: Biosample,
    seqRun: SequenceRun,
    alignment: Alignment,
    bamPath: String,
    onProgress: AnalysisProgress => Unit
  ): Either[String, (WorkspaceState, CallableLociResult)] = {
    try {
      onProgress(AnalysisProgress("Resolving reference genome...", 0.1))

      val referenceGateway = new ReferenceGateway((done, total) => {
        val pct = if (total > 0) 0.1 + (done.toDouble / total) * 0.2 else 0.1
        onProgress(AnalysisProgress(s"Downloading reference: ${done / 1024 / 1024}MB", pct))
      })

      val referencePath = referenceGateway.resolve(alignment.referenceBuild) match {
        case Right(path) => path.toString
        case Left(error) => return Left(s"Failed to resolve reference: $error")
      }

      onProgress(AnalysisProgress("Running callable loci analysis...", 0.3))

      val artifactCtx = ArtifactContext(
        sampleAccession = subject.sampleAccession,
        sequenceRunUri = seqRun.atUri,
        alignmentUri = alignment.atUri
      )

      val processor = new CallableLociProcessor()
      val resultEither = processor.process(bamPath, referencePath, (message, current, total) => {
        val pct = 0.3 + (current.toDouble / total) * 0.6
        onProgress(AnalysisProgress(s"Callable Loci: $message", pct))
      }, Some(artifactCtx))

      resultEither match {
        case Left(error) =>
          Left(s"Callable loci failed: ${error.getMessage}")

        case Right((result, _)) =>
          onProgress(AnalysisProgress("Saving results...", 0.95))

          // Update alignment metrics with callable bases
          val currentMetrics = alignment.metrics.getOrElse(AlignmentMetrics())
          val updatedMetrics = currentMetrics.copy(callableBases = Some(result.callableBases))
          val updatedAlignment = alignment.copy(
            metrics = Some(updatedMetrics),
            meta = alignment.meta.updated("callableLoci")
          )
          val newState = workspaceOps.updateAlignment(state, updatedAlignment)

          onProgress(AnalysisProgress("Callable loci analysis complete", 1.0, isComplete = true))
          Right((newState, result))
      }

    } catch {
      case e: Exception =>
        Left(e.getMessage)
    }
  }

  // --- Helper Methods ---

  private def inferTestType(stats: LibraryStats): String = {
    val avgReadLength = if (stats.lengthDistribution.nonEmpty) {
      val total = stats.lengthDistribution.map { case (len, count) => len.toLong * count }.sum
      val count = stats.lengthDistribution.values.sum
      if (count > 0) total / count else 0
    } else 0

    if (stats.inferredPlatform == "PacBio" && avgReadLength > 10000) "WGS_HIFI"
    else if (stats.inferredPlatform == "PacBio") "WGS_CLR"
    else if (stats.inferredPlatform == "Oxford Nanopore") "WGS_NANOPORE"
    else "WGS"
  }

  private def calculateMeanReadLength(distribution: Map[Int, Int]): Option[Int] = {
    if (distribution.isEmpty) None
    else {
      val totalReads = distribution.values.sum.toDouble
      val weightedSum = distribution.map { case (len, count) => len.toLong * count }.sum
      if (totalReads > 0) Some((weightedSum / totalReads).round.toInt) else None
    }
  }

  private def calculateMeanInsertSize(distribution: Map[Long, Int]): Option[Double] = {
    if (distribution.isEmpty) None
    else {
      val total = distribution.map { case (size, count) => size * count }.sum
      val count = distribution.values.sum
      if (count > 0) Some(total.toDouble / count) else None
    }
  }

}
