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
      // Fallback: if libraryLayout not yet set, check if readsPaired < 50% of totalReads
      val processor = new WgsMetricsProcessor()
      val isSingleEnd = seqRun.libraryLayout.exists(_.equalsIgnoreCase("Single-End")) ||
        (seqRun.libraryLayout.isEmpty && seqRun.totalReads.exists(total =>
          seqRun.readsPaired.forall(_ < total / 2)))
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

  // --- Whole-Genome Variant Calling ---

  /**
   * Runs whole-genome variant calling on an alignment.
   * Uses HaplotypeCaller for autosomes/X/Y and Mutect2 for mtDNA.
   * Results are cached in the alignment's artifact directory.
   */
  def runWholeGenomeVariantCalling(
    state: WorkspaceState,
    sampleAccession: String,
    sequenceRunIndex: Int,
    alignmentIndex: Int,
    onProgress: AnalysisProgress => Unit
  ): Future[Either[String, (WorkspaceState, CachedVcfInfo)]] = Future {
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
            alignments.lift(alignmentIndex) match {
              case None =>
                Left(s"Alignment not found at index $alignmentIndex")

              case Some(alignment) =>
                // Get BAM path from alignment's files first, then fall back to seqRun's files
                val bamPathOpt = alignment.files.headOption.orElse(seqRun.files.headOption).flatMap(_.location)
                bamPathOpt match {
                  case None =>
                    Left("No alignment file associated with this alignment")

                  case Some(bamPath) =>
                    runWholeGenomeVariantCallingInternal(state, subject, seqRun, alignment, bamPath, onProgress)
                }
            }
        }
    }
  }

  private def runWholeGenomeVariantCallingInternal(
    state: WorkspaceState,
    subject: Biosample,
    seqRun: SequenceRun,
    alignment: Alignment,
    bamPath: String,
    onProgress: AnalysisProgress => Unit
  ): Either[String, (WorkspaceState, CachedVcfInfo)] = {
    try {
      onProgress(AnalysisProgress("Resolving reference genome...", 0.05))

      val referenceGateway = new ReferenceGateway((done, total) => {
        val pct = if (total > 0) 0.05 + (done.toDouble / total) * 0.1 else 0.05
        onProgress(AnalysisProgress(s"Downloading reference: ${done / 1024 / 1024}MB", pct))
      })

      val referencePath = referenceGateway.resolve(alignment.referenceBuild) match {
        case Right(path) => path.toString
        case Left(error) => return Left(s"Failed to resolve reference: $error")
      }

      onProgress(AnalysisProgress("Starting whole-genome variant calling...", 0.15))

      // Get output directory
      val runId = seqRun.atUri.map(SubjectArtifactCache.extractIdFromUri).getOrElse("unknown-run")
      val alignId = alignment.atUri.map(SubjectArtifactCache.extractIdFromUri).getOrElse("unknown-alignment")
      val outputDir = VcfCache.getVcfDir(subject.sampleAccession, runId, alignId)

      // Run variant calling
      val caller = WholeGenomeVariantCaller()
      val result = caller.generateWholeGenomeVcf(
        bamPath = bamPath,
        referencePath = referencePath,
        outputDir = outputDir,
        referenceBuild = alignment.referenceBuild,
        onProgress = (msg, current, total) => {
          val pct = 0.15 + (current.toDouble / total) * 0.80
          onProgress(AnalysisProgress(msg, pct))
        }
      )

      result match {
        case Left(error) =>
          Left(s"Variant calling failed: $error")

        case Right(vcfInfo) =>
          onProgress(AnalysisProgress("Saving results...", 0.98))

          // Update alignment metrics with VCF info
          val existingMetrics = alignment.metrics.getOrElse(AlignmentMetrics())
          val updatedMetrics = existingMetrics.copy(
            vcfPath = Some(vcfInfo.vcfPath),
            vcfCreatedAt = Some(vcfInfo.createdAt),
            vcfVariantCount = Some(vcfInfo.variantCount),
            vcfReferenceBuild = Some(vcfInfo.referenceBuild),
            inferredSex = vcfInfo.inferredSex
          )

          val updatedAlignment = alignment.copy(
            metrics = Some(updatedMetrics),
            meta = alignment.meta.updated("vcf")
          )
          val newState = workspaceOps.updateAlignment(state, updatedAlignment)

          onProgress(AnalysisProgress("Whole-genome variant calling complete", 1.0, isComplete = true))
          Right((newState, vcfInfo))
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

          // Determine technology based on test type
          val technology = seqRun.testType match {
            case t if t.startsWith("BIGY") || t.contains("Y_ELITE") || t.contains("Y_PRIME") =>
              HaplogroupTechnology.BIG_Y
            case _ => HaplogroupTechnology.WGS
          }

          // Create a RunHaplogroupCall for the reconciliation system
          val runCall = RunHaplogroupCall(
            sourceRef = seqRun.atUri.getOrElse(s"local:sequencerun:unknown"),
            haplogroup = topResult.name,
            confidence = topResult.score,
            callMethod = CallMethod.SNP_PHYLOGENETIC,
            score = Some(topResult.score),
            supportingSnps = Some(topResult.matchingSnps),
            conflictingSnps = Some(topResult.mismatchingSnps),
            noCalls = None,
            technology = Some(technology),
            meanCoverage = None,
            treeProvider = Some(treeProviderType.toString.toLowerCase),
            treeVersion = None
          )

          // Convert TreeType to DnaType
          val dnaType = treeType match {
            case TreeType.YDNA => DnaType.Y_DNA
            case TreeType.MTDNA => DnaType.MT_DNA
          }

          // Add to reconciliation - this automatically updates biosample haplogroups with consensus
          workspaceOps.addHaplogroupCall(state, subject.sampleAccession, dnaType, runCall) match {
            case Right((newState, _)) =>
              onProgress(AnalysisProgress("Haplogroup analysis complete", 1.0, isComplete = true))
              Right((newState, topResult))
            case Left(err) =>
              Left(s"Failed to update reconciliation: $err")
          }

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

  // --- Comprehensive Batch Analysis ---

  /**
   * Runs a comprehensive batch analysis on a new sample:
   * 1. Read/WGS Metrics - Coverage depth analysis
   * 2. Callable Loci Metrics - Base-level coverage assessment
   * 3. Sex Inference - Determine biological sex from X:autosome ratio
   * 4. mtDNA Haplogroup - Maternal lineage determination
   * 5. Y-DNA Haplogroup - Paternal lineage (if male)
   * 6. Ancestral Composition - Stub for future implementation
   *
   * @param state Current workspace state
   * @param sampleAccession Sample accession identifier
   * @param sequenceRunIndex Index of sequence run
   * @param alignmentIndex Index of alignment
   * @param onProgress Progress callback with step information
   * @return Updated state and batch analysis results
   */
  def runComprehensiveAnalysis(
    state: WorkspaceState,
    sampleAccession: String,
    sequenceRunIndex: Int,
    alignmentIndex: Int,
    onProgress: AnalysisProgress => Unit
  ): Future[Either[String, (WorkspaceState, BatchAnalysisResult)]] = Future {
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
            alignments.lift(alignmentIndex) match {
              case None =>
                Left(s"Alignment not found at index $alignmentIndex")

              case Some(alignment) =>
                val bamPathOpt = alignment.files.headOption.orElse(seqRun.files.headOption).flatMap(_.location)
                bamPathOpt match {
                  case None =>
                    Left("No alignment file associated with this alignment")

                  case Some(bamPath) =>
                    runComprehensiveAnalysisInternal(state, subject, seqRun, alignment, bamPath, onProgress)
                }
            }
        }
    }
  }

  private def runComprehensiveAnalysisInternal(
    state: WorkspaceState,
    subject: Biosample,
    seqRun: SequenceRun,
    alignment: Alignment,
    bamPath: String,
    onProgress: AnalysisProgress => Unit
  ): Either[String, (WorkspaceState, BatchAnalysisResult)] = {
    var currentState = state
    var result = BatchAnalysisResult()

    try {
      // Step 0: Resolve reference genome (shared across steps)
      onProgress(AnalysisProgress("Resolving reference genome...", 0.02))
      val referenceGateway = new ReferenceGateway((done, total) => {
        val pct = if (total > 0) 0.02 + (done.toDouble / total) * 0.03 else 0.02
        onProgress(AnalysisProgress(s"Downloading reference: ${done / 1024 / 1024}MB", pct))
      })

      val referencePath = referenceGateway.resolve(alignment.referenceBuild) match {
        case Right(path) => path.toString
        case Left(error) => return Left(s"Failed to resolve reference: $error")
      }

      val artifactCtx = ArtifactContext(
        sampleAccession = subject.sampleAccession,
        sequenceRunUri = seqRun.atUri,
        alignmentUri = alignment.atUri
      )

      // Step 1: WGS Metrics (0.05 - 0.20)
      onProgress(AnalysisProgress("Step 1/6: Running coverage depth analysis...", 0.05))
      val wgsMetricsResult = runWgsMetricsStep(bamPath, referencePath, seqRun, artifactCtx, { pct =>
        onProgress(AnalysisProgress(s"Step 1/6: Coverage analysis (${(pct * 100).toInt}%)", 0.05 + pct * 0.15))
      })
      wgsMetricsResult match {
        case Right(wgsMetrics) =>
          result = result.copy(wgsMetrics = Some(wgsMetrics))
          // Update alignment metrics
          val alignmentMetrics = alignment.metrics.getOrElse(AlignmentMetrics()).copy(
            genomeTerritory = Some(wgsMetrics.genomeTerritory),
            meanCoverage = Some(wgsMetrics.meanCoverage),
            sdCoverage = Some(wgsMetrics.sdCoverage),
            medianCoverage = Some(wgsMetrics.medianCoverage),
            pct10x = Some(wgsMetrics.pct10x),
            pct20x = Some(wgsMetrics.pct20x),
            pct30x = Some(wgsMetrics.pct30x)
          )
          val updatedAlignment = alignment.copy(metrics = Some(alignmentMetrics))
          currentState = workspaceOps.updateAlignment(currentState, updatedAlignment)
        case Left(error) =>
          println(s"[BatchAnalysis] WGS metrics warning: $error")
          // Continue even if this step fails
      }

      // Step 2: Callable Loci (0.20 - 0.45)
      onProgress(AnalysisProgress("Step 2/6: Running callable loci analysis...", 0.20))
      val callableLociResult = runCallableLociStep(bamPath, referencePath, seqRun, artifactCtx, { pct =>
        onProgress(AnalysisProgress(s"Step 2/6: Callable loci (${(pct * 100).toInt}%)", 0.20 + pct * 0.25))
      })
      callableLociResult match {
        case Right((clResult, _)) =>
          result = result.copy(callableLociResult = Some(clResult))
          // Update alignment with callable bases
          val alignmentMetrics = currentState.workspace.main.alignments
            .find(_.atUri == alignment.atUri)
            .flatMap(_.metrics)
            .getOrElse(AlignmentMetrics())
            .copy(
              callableBases = Some(clResult.callableBases),
              callableLociComplete = Some(true)
            )
          val updatedAlignment = alignment.copy(metrics = Some(alignmentMetrics))
          currentState = workspaceOps.updateAlignment(currentState, updatedAlignment)
        case Left(error) =>
          println(s"[BatchAnalysis] Callable loci warning: $error")
      }

      // Step 3: Sex Inference (0.45 - 0.50)
      onProgress(AnalysisProgress("Step 3/6: Inferring biological sex...", 0.45))
      val sexResult = SexInference.inferFromBam(bamPath, (msg, pct) => {
        onProgress(AnalysisProgress(s"Step 3/6: Sex inference - $msg", 0.45 + pct * 0.05))
      })
      sexResult match {
        case Right(sr) =>
          result = result.copy(sexInferenceResult = Some(sr))
          // Update alignment metrics with sex inference
          val alignmentMetrics = currentState.workspace.main.alignments
            .find(_.atUri == alignment.atUri)
            .flatMap(_.metrics)
            .getOrElse(AlignmentMetrics())
            .copy(
              inferredSex = Some(sr.inferredSex.toString),
              sexInferenceConfidence = Some(sr.confidence),
              xAutosomeRatio = Some(sr.xAutosomeRatio)
            )
          val updatedAlignment = alignment.copy(metrics = Some(alignmentMetrics))
          currentState = workspaceOps.updateAlignment(currentState, updatedAlignment)
        case Left(error) =>
          println(s"[BatchAnalysis] Sex inference warning: $error")
      }

      // Step 4: mtDNA Haplogroup (0.50 - 0.70)
      onProgress(AnalysisProgress("Step 4/6: Determining mtDNA haplogroup...", 0.50))
      val mtDnaResult = runHaplogroupStep(bamPath, subject, seqRun, alignment, TreeType.MTDNA, artifactCtx, { pct =>
        onProgress(AnalysisProgress(s"Step 4/6: mtDNA haplogroup (${(pct * 100).toInt}%)", 0.50 + pct * 0.20))
      })
      mtDnaResult match {
        case Right(haplogroupResult) =>
          result = result.copy(mtDnaHaplogroup = Some(haplogroupResult))
        case Left(error) =>
          println(s"[BatchAnalysis] mtDNA haplogroup warning: $error")
      }

      // Step 5: Y-DNA Haplogroup (0.70 - 0.90) - Only if male
      val isMale = result.sexInferenceResult.exists(_.isMale)
      val sexUnknown = result.sexInferenceResult.exists(_.isUnknown)

      if (isMale || sexUnknown) {
        onProgress(AnalysisProgress("Step 5/6: Determining Y-DNA haplogroup...", 0.70))
        val yDnaResult = runHaplogroupStep(bamPath, subject, seqRun, alignment, TreeType.YDNA, artifactCtx, { pct =>
          onProgress(AnalysisProgress(s"Step 5/6: Y-DNA haplogroup (${(pct * 100).toInt}%)", 0.70 + pct * 0.20))
        })
        yDnaResult match {
          case Right(haplogroupResult) =>
            result = result.copy(yDnaHaplogroup = Some(haplogroupResult))
          case Left(error) =>
            println(s"[BatchAnalysis] Y-DNA haplogroup warning: $error")
            result = result.copy(skippedYDna = true, skippedYDnaReason = Some(error))
        }
      } else {
        onProgress(AnalysisProgress("Step 5/6: Skipping Y-DNA (female sample)...", 0.70))
        result = result.copy(skippedYDna = true, skippedYDnaReason = Some("Sample inferred as female"))
      }

      // Step 6: Ancestral Composition Stub (0.90 - 1.0)
      onProgress(AnalysisProgress("Step 6/6: Preparing for ancestry analysis...", 0.90))
      // This is a stub for future implementation
      // Will use autosomal SNPs from the VCF for ancestry composition
      result = result.copy(ancestryStub = true)
      onProgress(AnalysisProgress("Ancestry analysis will be available in a future update", 0.95))

      onProgress(AnalysisProgress("Comprehensive analysis complete!", 1.0, isComplete = true))
      Right((currentState, result))

    } catch {
      case e: Exception =>
        Left(s"Batch analysis failed: ${e.getMessage}")
    }
  }

  /**
   * Run WGS metrics step for batch analysis.
   */
  private def runWgsMetricsStep(
    bamPath: String,
    referencePath: String,
    seqRun: SequenceRun,
    artifactCtx: ArtifactContext,
    onProgress: Double => Unit
  ): Either[String, WgsMetrics] = {
    val processor = new WgsMetricsProcessor()
    val isSingleEnd = seqRun.libraryLayout.exists(_.equalsIgnoreCase("Single-End")) ||
      (seqRun.libraryLayout.isEmpty && seqRun.totalReads.exists(total =>
        seqRun.readsPaired.forall(_ < total / 2)))

    processor.process(
      bamPath = bamPath,
      referencePath = referencePath,
      onProgress = (_, current, total) => {
        if (total > 0) onProgress(current.toDouble / total)
      },
      readLength = seqRun.maxReadLength,
      artifactContext = Some(artifactCtx),
      totalReads = seqRun.totalReads,
      countUnpaired = isSingleEnd
    ).left.map(_.getMessage)
  }

  /**
   * Run callable loci step for batch analysis.
   */
  private def runCallableLociStep(
    bamPath: String,
    referencePath: String,
    seqRun: SequenceRun,
    artifactCtx: ArtifactContext,
    onProgress: Double => Unit
  ): Either[String, (CallableLociResult, List[String])] = {
    val processor = new CallableLociProcessor()
    val minDepth = seqRun.testType.toUpperCase match {
      case t if t.contains("HIFI") => 2
      case _ => 4
    }

    processor.process(
      bamPath = bamPath,
      referencePath = referencePath,
      onProgress = (_, current, total) => {
        if (total > 0) onProgress(current.toDouble / total)
      },
      artifactContext = Some(artifactCtx),
      minDepth = minDepth
    ).left.map(_.getMessage)
  }

  /**
   * Run haplogroup analysis step for batch analysis.
   * Returns just the haplogroup result - state updates are handled by the caller.
   */
  private def runHaplogroupStep(
    bamPath: String,
    subject: Biosample,
    seqRun: SequenceRun,
    alignment: Alignment,
    treeType: TreeType,
    artifactCtx: ArtifactContext,
    onProgress: Double => Unit
  ): Either[String, AnalysisHaplogroupResult] = {
    // Build LibraryStats from existing data for the processor
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

    val processor = new HaplogroupProcessor()
    processor.analyze(
      bamPath,
      libraryStats,
      treeType,
      treeProviderType,
      (_, current, total) => {
        val pct = if (total > 0) current / total else 0.0
        onProgress(pct)
      },
      Some(artifactCtx)
    ).map(_.headOption.getOrElse(
      // Return a default result if no haplogroups found
      AnalysisHaplogroupResult(
        name = "Unknown",
        score = 0.0,
        matchingSnps = 0,
        mismatchingSnps = 0,
        ancestralMatches = 0,
        noCalls = 0,
        totalSnps = 0,
        cumulativeSnps = 0,
        depth = 0
      )
    ))
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

/**
 * Result of comprehensive batch analysis.
 */
case class BatchAnalysisResult(
  wgsMetrics: Option[WgsMetrics] = None,
  callableLociResult: Option[CallableLociResult] = None,
  sexInferenceResult: Option[SexInference.SexInferenceResult] = None,
  mtDnaHaplogroup: Option[AnalysisHaplogroupResult] = None,
  yDnaHaplogroup: Option[AnalysisHaplogroupResult] = None,
  skippedYDna: Boolean = false,
  skippedYDnaReason: Option[String] = None,
  ancestryStub: Boolean = false  // Stub for future ancestry composition
)
