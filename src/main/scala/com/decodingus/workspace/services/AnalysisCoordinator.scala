package com.decodingus.workspace.services

import com.decodingus.analysis.*
import com.decodingus.haplogroup.processor.HaplogroupProcessor
import com.decodingus.analysis.SexInference.{InferredSex, SexInferenceResult}
import com.decodingus.config.{FeatureToggles, ReferenceConfigService, UserPreferencesService}
import com.decodingus.haplogroup.model.HaplogroupResult as AnalysisHaplogroupResult
import com.decodingus.haplogroup.tree.{TreeProviderType, TreeType}
import com.decodingus.model.{LibraryStats, WgsMetrics}
import com.decodingus.refgenome.{ReferenceGateway, ReferenceResolveResult}
import com.decodingus.util.Logger
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

  private val log = Logger[AnalysisCoordinator]
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

      // Find existing alignment for THIS reference build (not just the first one)
      val existingAlignment = seqRun.alignmentRefs.flatMap { ref =>
        state.workspace.main.alignments.find(a => a.atUri.contains(ref) && a.referenceBuild == libraryStats.referenceBuild)
      }.headOption

      val alignUri = existingAlignment.flatMap(_.atUri).getOrElse(
        s"local:alignment:${subject.sampleAccession}:${libraryStats.referenceBuild}:${java.util.UUID.randomUUID().toString.take(8)}"
      )

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

      // Check if we have cached artifacts from previous analysis
      val runId = seqRun.atUri.map(SubjectArtifactCache.extractIdFromUri).getOrElse("unknown-run")
      val alignId = alignment.atUri.map(SubjectArtifactCache.extractIdFromUri).getOrElse("unknown-alignment")
      val vcfDir = VcfCache.getVcfDir(subject.sampleAccession, runId, alignId)
      val cachedVcfPath = vcfDir.resolve("whole_genome.vcf.gz")

      // Also check for contig-specific VCF (generated by previous haplogroup analysis)
      val contigName = if (treeType == TreeType.YDNA) "chrY" else "chrM"
      val contigVcfPath = vcfDir.resolve(s"$contigName.vcf.gz")

      // Check for callable loci BED file (indicates callable loci analysis was completed)
      val callableLociDir = artifactCtx.getSubdir("callable_loci")
      val callableBedPath = callableLociDir.resolve(s"$contigName.callable.bed")

      val treeTypeStr = if (treeType == TreeType.YDNA) "Y-DNA" else "mtDNA"

      // Check for vendor-provided VCF first (e.g., FTDNA Big Y)
      // Check both alignment-level (for BAM-based imports) and run-level (for VCF-only imports)
      val vendorVcf: Option[VendorVcfInfo] = if (treeType == TreeType.YDNA) {
        VcfCache.findYDnaVendorVcf(subject.sampleAccession, runId, alignId)
          .orElse(VcfCache.findYDnaRunVendorVcf(subject.sampleAccession, runId))
      } else {
        VcfCache.findMtDnaVendorVcf(subject.sampleAccession, runId, alignId)
          .orElse(VcfCache.findMtDnaRunVendorVcf(subject.sampleAccession, runId))
      }

      // Check for vendor-provided FASTA (mtDNA only - e.g., FTDNA mtFull Sequence, YSEQ mtDNA)
      val vendorFasta: Option[VendorFastaInfo] = if (treeType == TreeType.MTDNA) {
        VcfCache.findMtDnaRunFasta(subject.sampleAccession, runId)
      } else {
        None
      }

      val result: Either[String, List[AnalysisHaplogroupResult]] =
        if (vendorFasta.isDefined && treeType == TreeType.MTDNA) {
          // Use vendor-provided FASTA for mtDNA - highest priority for mtDNA
          val vfasta = vendorFasta.get
          log.info(s"Using ${vfasta.vendor.displayName} FASTA for mtDNA haplogroup analysis: ${vfasta.fastaPath}")
          onProgress(AnalysisProgress(s"Using ${vfasta.vendor.displayName} FASTA for mtDNA analysis...", 0.2))
          val haplogroupOutputDir = artifactCtx.getSubdir("haplogroup")
          processor.analyzeFromFasta(
            fastaPath = vfasta.fastaPath,
            treeProviderType = treeProviderType,
            onProgress = (message, current, total) => {
              val pct = if (total > 0) 0.2 + (current / total) * 0.7 else 0.2
              onProgress(AnalysisProgress(message, pct))
            },
            outputDir = Some(haplogroupOutputDir)
          )
        } else if (vendorVcf.isDefined) {
          // Use vendor-provided VCF - highest priority
          val vvcf = vendorVcf.get
          log.info(s"Using ${vvcf.vendor.displayName} VCF for $treeTypeStr haplogroup analysis: ${vvcf.vcfPath}")
          onProgress(AnalysisProgress(s"Using ${vvcf.vendor.displayName} VCF for $treeTypeStr analysis...", 0.2))
          // Analyze directly from the vendor VCF path
          val haplogroupOutputDir = artifactCtx.getSubdir("haplogroup")
          processor.analyzeFromVcfFile(
            vcfPath = vvcf.vcfPath,
            referenceBuild = vvcf.referenceBuild,
            treeType = treeType,
            treeProviderType = treeProviderType,
            onProgress = (message, current, total) => {
              val pct = if (total > 0) 0.2 + (current / total) * 0.7 else 0.2
              onProgress(AnalysisProgress(message, pct))
            },
            outputDir = Some(haplogroupOutputDir)
          )
        } else if (java.nio.file.Files.exists(cachedVcfPath)) {
          // Use cached whole-genome VCF - second priority
          log.info(s"Using cached whole-genome VCF for $treeTypeStr haplogroup analysis: $cachedVcfPath")
          onProgress(AnalysisProgress(s"Using cached VCF for $treeTypeStr analysis...", 0.2))
          processor.analyzeFromCachedVcf(
            sampleAccession = subject.sampleAccession,
            runId = runId,
            alignmentId = alignId,
            referenceBuild = alignment.referenceBuild,
            treeType = treeType,
            treeProviderType = treeProviderType,
            onProgress = (message, current, total) => {
              val pct = if (total > 0) 0.2 + (current / total) * 0.7 else 0.2
              onProgress(AnalysisProgress(message, pct))
            }
          )
        } else if (java.nio.file.Files.exists(contigVcfPath)) {
          // Use contig-specific VCF from previous haplogroup analysis
          log.info(s"Using cached contig VCF for $treeTypeStr haplogroup analysis: $contigVcfPath")
          onProgress(AnalysisProgress(s"Using cached $contigName VCF for analysis...", 0.2))
          processor.analyzeFromCachedVcf(
            sampleAccession = subject.sampleAccession,
            runId = runId,
            alignmentId = alignId,
            referenceBuild = alignment.referenceBuild,
            treeType = treeType,
            treeProviderType = treeProviderType,
            onProgress = (message, current, total) => {
              val pct = if (total > 0) 0.2 + (current / total) * 0.7 else 0.2
              onProgress(AnalysisProgress(message, pct))
            }
          )
        } else {
          // No cached VCF - fall back to BAM-based variant calling
          if (java.nio.file.Files.exists(callableBedPath)) {
            log.info(s"Callable loci BED found at $callableBedPath - will use for quality filtering")
          }
          log.info(s"No cached VCF found, using BAM-based calling for $treeTypeStr haplogroup analysis")
          onProgress(AnalysisProgress(s"Calling variants from BAM for $treeTypeStr...", 0.2))

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

          val analysisResult = processor.analyze(
            bamPath,
            libraryStats,
            treeType,
            treeProviderType,
            (message, current, total) => {
              val pct = if (total > 0) 0.2 + (current / total) * 0.7 else 0.2
              onProgress(AnalysisProgress(message, pct))
            },
            Some(artifactCtx)
          )

          // Save the generated VCF to cache for future use
          analysisResult.foreach { _ =>
            val haplogroupDir = artifactCtx.getSubdir("haplogroup")
            val prefix = if (treeType == TreeType.YDNA) "ydna" else "mtdna"
            val sourceVcf = haplogroupDir.resolve(s"${prefix}_calls.vcf")

            if (java.nio.file.Files.exists(sourceVcf)) {
              try {
                java.nio.file.Files.createDirectories(vcfDir)
                val destVcf = vcfDir.resolve(s"$contigName.vcf.gz")
                GatkRunner.run(Array(
                  "SortVcf",
                  "-I", sourceVcf.toString,
                  "-O", destVcf.toString,
                  "--CREATE_INDEX", "true"
                )) match {
                  case Right(_) =>
                    log.info(s"Saved $treeTypeStr VCF to cache: $destVcf")
                  case Left(err) =>
                    log.warn(s"Failed to cache VCF: $err")
                }
              } catch {
                case e: Exception =>
                  log.warn(s"Failed to cache VCF: ${e.getMessage}")
              }
            }
          }

          analysisResult
        }

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

      // Determine minDepth based on test type AND coverage (same logic as batch analysis)
      val meanCoverage = alignment.metrics.flatMap(_.meanCoverage).getOrElse(30.0)
      val isLowPass = meanCoverage <= 5.0
      val isHiFi = seqRun.testType.toUpperCase.contains("HIFI")
      val isLongRead = seqRun.testType.toUpperCase.contains("NANOPORE") ||
                       seqRun.testType.toUpperCase.contains("CLR") ||
                       seqRun.maxReadLength.exists(_ > 10000)

      val minDepth = if (isHiFi) {
        2  // HiFi: high accuracy, minDepth=2 is fine
      } else if (isLowPass) {
        // Low-pass data: use minDepth proportional to coverage
        math.max(1, (meanCoverage / 2).toInt)
      } else if (isLongRead) {
        3  // ONT/CLR: moderate accuracy, minDepth=3
      } else {
        4  // Illumina WGS at normal depth: standard minDepth=4
      }

      log.info(s"[CallableLoci] Using minDepth=$minDepth (testType=${seqRun.testType}, meanCov=${f"$meanCoverage%.1f"}x, isHiFi=$isHiFi, isLowPass=$isLowPass)")

      val resultEither = processor.process(bamPath, referencePath, (message, current, total) => {
        val pct = 0.3 + (current.toDouble / total) * 0.6
        onProgress(AnalysisProgress(s"Callable Loci: $message", pct))
      }, Some(artifactCtx), minDepth)

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
   * 1. Read Metrics - Read length and alignment metrics
   * 2. WGS Metrics - Coverage depth analysis
   * 3. Callable Loci - Base-level coverage assessment
   * 4. Sex Inference - Determine biological sex from X:autosome ratio
   * 5. Variant Calling - Generate whole-genome VCF (uses sex for ploidy)
   * 6. mtDNA Haplogroup - Maternal lineage determination
   * 7. Y-DNA Haplogroup - Paternal lineage (if male)
   * 8. Ancestral Composition - Stub for future implementation
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

    // IMPORTANT: Track the current sequence run separately from the parameter.
    // The parameter `seqRun` is immutable, but we update currentState throughout.
    // We need to track the latest version to avoid using stale data in later steps.
    var currentSeqRun = seqRun

    // Helper to get the latest alignment from state (since we update it in multiple steps)
    def getCurrentAlignment: Alignment =
      currentState.workspace.main.alignments
        .find(_.atUri == alignment.atUri)
        .getOrElse(alignment)

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

      // Load checkpoint - resume from last successful step if BAM hasn't changed
      val artifactDir = artifactCtx.getArtifactDir
      var checkpoint = AnalysisCheckpoint.loadAndValidate(artifactDir, bamPath)

      // Track the read length - from checkpoint or will be set by step 1
      var effectiveReadLength: Option[Int] = checkpoint.maxReadLength.orElse(seqRun.maxReadLength)

      // Step 1: Read Metrics (0.02 - 0.07) - MUST run before WGS Metrics to get read length
      if (!checkpoint.readMetricsCompleted) {
        onProgress(AnalysisProgress("Step 1/8: Collecting read metrics...", 0.02))
        val readMetricsProcessor = new UnifiedMetricsProcessor()
        val readMetricsResult = readMetricsProcessor.process(
          bamPath = bamPath,
          referencePath = referencePath,
          onProgress = (msg, current, total) => {
            val pct = if (total > 0) current / total else 0.0
            onProgress(AnalysisProgress(s"Step 1/8: Read metrics - $msg", 0.02 + pct * 0.05))
          },
          artifactContext = Some(artifactCtx)
        )
        readMetricsResult match {
          case Right(readMetrics) =>
            result = result.copy(readMetrics = Some(readMetrics))
            effectiveReadLength = Some(readMetrics.maxReadLength)

            // Update sequence run with read length and other metrics
            // IMPORTANT: Use currentSeqRun (not the stale seqRun parameter) to preserve any prior updates

            // Infer testType based on collected metrics if not already set or is generic WGS
            // Platform detection uses read length as proxy: maxReadLength > 10000 implies long-read (PacBio/ONT)
            val inferredTestType = currentSeqRun.testType match {
              case t if t.contains("HIFI") || t.contains("CLR") || t.contains("NANOPORE") =>
                // Already a specific long-read type, keep it
                t
              case "WGS" | "Unknown" | "" =>
                // Generic type - infer from read length
                if (readMetrics.maxReadLength > 10000) {
                  // Long reads - use platform hint if available, otherwise assume HiFi (most common now)
                  if (currentSeqRun.platformName.toUpperCase.contains("NANOPORE") ||
                      currentSeqRun.platformName.toUpperCase.contains("ONT")) {
                    "WGS_NANOPORE"
                  } else {
                    "WGS_HIFI"  // Default to HiFi for PacBio or unknown long-read
                  }
                } else {
                  currentSeqRun.testType  // Keep existing for short reads
                }
              case other => other  // Keep any other specific type
            }

            // Infer library layout from read pairing
            val inferredLayout = if (readMetrics.readsAlignedInPairs > readMetrics.totalReads / 2) "Paired-End" else "Single-End"

            val updatedSeqRun = currentSeqRun.copy(
              maxReadLength = Some(readMetrics.maxReadLength),
              readLength = Some(readMetrics.meanReadLength.toInt),
              totalReads = Some(readMetrics.totalReads),
              pfReads = Some(readMetrics.pfReads),
              pfReadsAligned = Some(readMetrics.pfReadsAligned),
              meanInsertSize = Some(readMetrics.meanInsertSize),
              medianInsertSize = Some(readMetrics.medianInsertSize.toInt),
              stdInsertSize = Some(readMetrics.stdInsertSize),
              testType = inferredTestType,
              libraryLayout = Some(inferredLayout)
            )
            currentState = workspaceOps.updateSequenceRunByUri(currentState, updatedSeqRun)
            currentSeqRun = updatedSeqRun  // Keep local var in sync with state

            checkpoint = AnalysisCheckpoint.markReadMetricsComplete(artifactDir, checkpoint, readMetrics.maxReadLength)
            log.info(s"Read metrics complete: maxReadLength=${readMetrics.maxReadLength}, testType=$inferredTestType, layout=$inferredLayout")
          case Left(error) =>
            log.warn(s"Read metrics warning: ${error.getMessage}")
            // Mark complete to continue, but WGS metrics will use seqRun.maxReadLength fallback
            checkpoint = AnalysisCheckpoint.markStepComplete(artifactDir, checkpoint, 1)
        }
      } else {
        onProgress(AnalysisProgress("Step 1/8: Read metrics (cached)", 0.07))
        log.debug(s"Using cached read length: ${checkpoint.maxReadLength}")
      }

      // Step 2: WGS Metrics (0.07 - 0.17) - uses read length from step 1
      if (!checkpoint.wgsMetricsCompleted) {
        onProgress(AnalysisProgress("Step 2/8: Running coverage depth analysis...", 0.07))

        // Use currentSeqRun which has updated read length from Step 1
        // Also patch effectiveReadLength in case Step 1 was cached but checkpoint has it
        val seqRunForWgs = currentSeqRun.copy(maxReadLength = effectiveReadLength)
        log.debug(s"WGS Metrics using maxReadLength: ${seqRunForWgs.maxReadLength}")

        val wgsMetricsResult = runWgsMetricsStep(bamPath, referencePath, seqRunForWgs, artifactCtx, { pct =>
          onProgress(AnalysisProgress(s"Step 2/8: Coverage analysis (${(pct * 100).toInt}%)", 0.07 + pct * 0.10))
        })
        wgsMetricsResult match {
          case Right(wgsMetrics) =>
            result = result.copy(wgsMetrics = Some(wgsMetrics))
            // Update alignment metrics - use getCurrentAlignment to get latest state
            val currentAlign = getCurrentAlignment
            val alignmentMetrics = currentAlign.metrics.getOrElse(AlignmentMetrics()).copy(
              genomeTerritory = Some(wgsMetrics.genomeTerritory),
              meanCoverage = Some(wgsMetrics.meanCoverage),
              sdCoverage = Some(wgsMetrics.sdCoverage),
              medianCoverage = Some(wgsMetrics.medianCoverage),
              pct10x = Some(wgsMetrics.pct10x),
              pct20x = Some(wgsMetrics.pct20x),
              pct30x = Some(wgsMetrics.pct30x)
            )
            val updatedAlignment = currentAlign.copy(metrics = Some(alignmentMetrics))
            currentState = workspaceOps.updateAlignment(currentState, updatedAlignment)
            checkpoint = AnalysisCheckpoint.markStepComplete(artifactDir, checkpoint, 2)
          case Left(error) =>
            log.warn(s"WGS metrics warning: $error")
            // Mark complete anyway to allow continuing (metrics are optional)
            checkpoint = AnalysisCheckpoint.markStepComplete(artifactDir, checkpoint, 2)
        }
      } else {
        onProgress(AnalysisProgress("Step 2/8: Coverage analysis (cached)", 0.17))
      }

      // Step 3: Callable Loci (0.17 - 0.30)
      if (!checkpoint.callableLociCompleted) {
        onProgress(AnalysisProgress("Step 3/8: Running callable loci analysis...", 0.17))
        // Use currentSeqRun for testType and getCurrentAlignment for coverage metrics in minDepth calculation
        val callableLociResult = runCallableLociStep(bamPath, referencePath, currentSeqRun, getCurrentAlignment, artifactCtx, { pct =>
          onProgress(AnalysisProgress(s"Step 3/8: Callable loci (${(pct * 100).toInt}%)", 0.17 + pct * 0.13))
        })
        callableLociResult match {
          case Right((clResult, _)) =>
            result = result.copy(callableLociResult = Some(clResult))
            // Update alignment with callable bases - use getCurrentAlignment for latest metrics
            val currentAlign = getCurrentAlignment
            val alignmentMetrics = currentAlign.metrics
              .getOrElse(AlignmentMetrics())
              .copy(
                callableBases = Some(clResult.callableBases),
                callableLociComplete = Some(true)
              )
            val updatedAlignment = currentAlign.copy(metrics = Some(alignmentMetrics))
            currentState = workspaceOps.updateAlignment(currentState, updatedAlignment)
            checkpoint = AnalysisCheckpoint.markStepComplete(artifactDir, checkpoint, 3)
          case Left(error) =>
            log.warn(s"Callable loci warning: $error")
            checkpoint = AnalysisCheckpoint.markStepComplete(artifactDir, checkpoint, 3)
        }
      } else {
        onProgress(AnalysisProgress("Step 3/8: Callable loci (cached)", 0.30))
      }

      // Step 4: Sex Inference (0.30 - 0.35)
      // Check if user already provided sex on the biosample - skip inference if so
      val userProvidedSex: Option[SexInferenceResult] = subject.sex.flatMap { sexStr =>
        sexStr.toLowerCase match {
          case "male" => Some(SexInferenceResult(
            inferredSex = InferredSex.Male,
            xAutosomeRatio = 0.5, // Expected male ratio
            autosomeMeanCoverage = 0.0,
            xCoverage = 0.0,
            confidence = "user-provided"
          ))
          case "female" => Some(SexInferenceResult(
            inferredSex = InferredSex.Female,
            xAutosomeRatio = 1.0, // Expected female ratio
            autosomeMeanCoverage = 0.0,
            xCoverage = 0.0,
            confidence = "user-provided"
          ))
          case _ => None // "Other", "Unknown", or unrecognized - need to infer
        }
      }

      if (!checkpoint.sexInferenceCompleted) {
        userProvidedSex match {
          case Some(userSex) =>
            // Use user-provided sex - skip BAM scanning
            onProgress(AnalysisProgress(s"Step 4/8: Using user-provided sex (${userSex.inferredSex})...", 0.30))
            result = result.copy(sexInferenceResult = Some(userSex))
            log.info(s"Using user-provided sex: ${userSex.inferredSex}")
            // Update alignment metrics
            val alignmentMetrics = currentState.workspace.main.alignments
              .find(_.atUri == alignment.atUri)
              .flatMap(_.metrics)
              .getOrElse(AlignmentMetrics())
              .copy(
                inferredSex = Some(userSex.inferredSex.toString),
                sexInferenceConfidence = Some(userSex.confidence),
                xAutosomeRatio = Some(userSex.xAutosomeRatio)
              )
            val updatedAlignment = alignment.copy(metrics = Some(alignmentMetrics))
            currentState = workspaceOps.updateAlignment(currentState, updatedAlignment)
            checkpoint = AnalysisCheckpoint.markStepComplete(artifactDir, checkpoint, 4)

          case None =>
            // Infer sex from BAM coverage
            onProgress(AnalysisProgress("Step 4/8: Inferring biological sex...", 0.30))

            // Ensure BAM index exists before sex inference (it requires indexed BAM)
            GatkRunner.ensureIndex(bamPath) match {
              case Left(indexError) =>
                log.warn(s"Failed to create BAM index for sex inference: $indexError")
              case Right(_) => // Index exists or was created
            }

            val sexResult = SexInference.inferFromBam(bamPath, (msg, pct) => {
              onProgress(AnalysisProgress(s"Step 4/8: Sex inference - $msg", 0.30 + pct * 0.05))
            })
            sexResult match {
              case Right(sr) =>
                result = result.copy(sexInferenceResult = Some(sr))
                log.info(s"Sex inferred: ${sr.inferredSex} (confidence: ${sr.confidence}, X:autosome ratio: ${sr.xAutosomeRatio})")
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
                checkpoint = AnalysisCheckpoint.markStepComplete(artifactDir, checkpoint, 4)
              case Left(error) =>
                log.warn(s"Sex inference failed: $error - using Unknown sex (will attempt Y-DNA analysis)")
                // Set Unknown sex instead of leaving as None - this allows Y-DNA to proceed
                val unknownSex = SexInferenceResult(
                  inferredSex = InferredSex.Unknown,
                  xAutosomeRatio = 0.0,
                  autosomeMeanCoverage = 0.0,
                  xCoverage = 0.0,
                  confidence = "failed"
                )
                result = result.copy(sexInferenceResult = Some(unknownSex))
                // Update alignment metrics to indicate inference failed
                val alignmentMetrics = currentState.workspace.main.alignments
                  .find(_.atUri == alignment.atUri)
                  .flatMap(_.metrics)
                  .getOrElse(AlignmentMetrics())
                  .copy(
                    inferredSex = Some(InferredSex.Unknown.toString),
                    sexInferenceConfidence = Some("failed"),
                    xAutosomeRatio = None
                  )
                val updatedAlignment = alignment.copy(metrics = Some(alignmentMetrics))
                currentState = workspaceOps.updateAlignment(currentState, updatedAlignment)
                checkpoint = AnalysisCheckpoint.markStepComplete(artifactDir, checkpoint, 4)
            }
        }
      } else {
        onProgress(AnalysisProgress("Step 4/8: Sex inference (cached)", 0.35))
        // Load cached sex inference result from alignment metrics
        val cachedSex = currentState.workspace.main.alignments
          .find(_.atUri == alignment.atUri)
          .flatMap(_.metrics)
          .flatMap(m => m.inferredSex.map(s => SexInferenceResult(
            inferredSex = InferredSex.valueOf(s),
            xAutosomeRatio = m.xAutosomeRatio.getOrElse(0.0),
            autosomeMeanCoverage = 0.0, // Not stored in metrics, use default
            xCoverage = 0.0, // Not stored in metrics, use default
            confidence = m.sexInferenceConfidence.getOrElse("unknown")
          )))
        // If cached lookup fails, use Unknown sex (allows Y-DNA to proceed)
        val sexResult = cachedSex.getOrElse {
          log.warn("Could not load cached sex inference - using Unknown (will attempt Y-DNA analysis)")
          SexInferenceResult(
            inferredSex = InferredSex.Unknown,
            xAutosomeRatio = 0.0,
            autosomeMeanCoverage = 0.0,
            xCoverage = 0.0,
            confidence = "cache-miss"
          )
        }
        result = result.copy(sexInferenceResult = Some(sexResult))
      }

      // Step 5: Variant Calling (0.35 - 0.55) - uses sex for ploidy
      if (!checkpoint.variantCallingCompleted) {
        onProgress(AnalysisProgress("Step 5/8: Running whole-genome variant calling...", 0.35))

        // Get output directory for VCF - use currentSeqRun for consistency
        val runId = currentSeqRun.atUri.map(SubjectArtifactCache.extractIdFromUri).getOrElse("unknown-run")
        val currentAlign = getCurrentAlignment
        val alignId = currentAlign.atUri.map(SubjectArtifactCache.extractIdFromUri).getOrElse("unknown-alignment")
        val outputDir = VcfCache.getVcfDir(subject.sampleAccession, runId, alignId)

        val caller = WholeGenomeVariantCaller()
        val vcfResult = caller.generateWholeGenomeVcf(
          bamPath = bamPath,
          referencePath = referencePath,
          outputDir = outputDir,
          referenceBuild = currentAlign.referenceBuild,
          onProgress = (msg, current, total) => {
            val pct = if (total > 0) current.toDouble / total else 0.0
            onProgress(AnalysisProgress(s"Step 5/8: Variant calling - $msg", 0.35 + pct * 0.20))
          },
          sexInferenceResult = result.sexInferenceResult  // Pass sex result from Step 4 to avoid re-computing
        )
        vcfResult match {
          case Right(vcfInfo) =>
            result = result.copy(vcfInfo = Some(vcfInfo))
            // Update alignment metrics with VCF info - re-fetch to get latest
            val latestAlign = getCurrentAlignment
            val alignmentMetrics = latestAlign.metrics
              .getOrElse(AlignmentMetrics())
              .copy(
                vcfPath = Some(vcfInfo.vcfPath),
                vcfCreatedAt = Some(vcfInfo.createdAt),
                vcfVariantCount = Some(vcfInfo.variantCount),
                vcfReferenceBuild = Some(vcfInfo.referenceBuild)
              )
            val updatedAlignment = latestAlign.copy(metrics = Some(alignmentMetrics))
            currentState = workspaceOps.updateAlignment(currentState, updatedAlignment)
            checkpoint = AnalysisCheckpoint.markStepComplete(artifactDir, checkpoint, 5)
            log.info(s"Variant calling complete: ${vcfInfo.variantCount} variants")
          case Left(error) =>
            log.warn(s"Variant calling warning: $error")
            checkpoint = AnalysisCheckpoint.markStepComplete(artifactDir, checkpoint, 5)
        }
      } else {
        onProgress(AnalysisProgress("Step 5/8: Variant calling (cached)", 0.55))
      }

      // Step 6: mtDNA Haplogroup (0.55 - 0.70)
      if (!checkpoint.mtDnaHaplogroupCompleted) {
        onProgress(AnalysisProgress("Step 6/8: Determining mtDNA haplogroup...", 0.55))
        // Use currentSeqRun and getCurrentAlignment for latest state
        val mtDnaResult = runHaplogroupStep(bamPath, subject, currentSeqRun, getCurrentAlignment, TreeType.MTDNA, artifactCtx, { pct =>
          onProgress(AnalysisProgress(s"Step 6/8: mtDNA haplogroup (${(pct * 100).toInt}%)", 0.55 + pct * 0.15))
        })
        mtDnaResult match {
          case Right(haplogroupResult) =>
            result = result.copy(mtDnaHaplogroup = Some(haplogroupResult))
            checkpoint = AnalysisCheckpoint.markStepComplete(artifactDir, checkpoint, 6)
          case Left(error) =>
            log.warn(s"mtDNA haplogroup warning: $error")
            checkpoint = AnalysisCheckpoint.markStepComplete(artifactDir, checkpoint, 6)
        }
      } else {
        onProgress(AnalysisProgress("Step 6/8: mtDNA haplogroup (cached)", 0.70))
      }

      // Step 7: Y-DNA Haplogroup (0.70 - 0.92) - Only if male or unknown
      // Note: We run Y-DNA for male, unknown, OR if sex inference is completely missing (None)
      // This ensures we don't incorrectly skip Y-DNA due to inference failures
      if (!checkpoint.yDnaHaplogroupCompleted && !checkpoint.yDnaSkipped) {
        val isMale = result.sexInferenceResult.exists(_.isMale)
        val sexUnknown = result.sexInferenceResult.exists(_.isUnknown)
        val sexMissing = result.sexInferenceResult.isEmpty  // Defensive: treat missing as "try Y-DNA"
        val isFemale = result.sexInferenceResult.exists(_.isFemale)

        if (isMale || sexUnknown || sexMissing) {
          val reason = if (isMale) "male sample" else if (sexMissing) "sex unknown (missing)" else "sex unknown"
          onProgress(AnalysisProgress(s"Step 7/8: Determining Y-DNA haplogroup ($reason)...", 0.70))
          // Use currentSeqRun and getCurrentAlignment for latest state
          val yDnaResult = runHaplogroupStep(bamPath, subject, currentSeqRun, getCurrentAlignment, TreeType.YDNA, artifactCtx, { pct =>
            onProgress(AnalysisProgress(s"Step 7/8: Y-DNA haplogroup (${(pct * 100).toInt}%)", 0.70 + pct * 0.22))
          })
          yDnaResult match {
            case Right(haplogroupResult) =>
              result = result.copy(yDnaHaplogroup = Some(haplogroupResult))
              checkpoint = AnalysisCheckpoint.markStepComplete(artifactDir, checkpoint, 7)
            case Left(error) =>
              log.warn(s"Y-DNA haplogroup warning: $error")
              result = result.copy(skippedYDna = true, skippedYDnaReason = Some(error))
              checkpoint = AnalysisCheckpoint.markStepComplete(artifactDir, checkpoint, 7, skipped = true)
          }
        } else if (isFemale) {
          onProgress(AnalysisProgress("Step 7/8: Skipping Y-DNA (female sample)...", 0.70))
          result = result.copy(skippedYDna = true, skippedYDnaReason = Some("Sample inferred as female"))
          checkpoint = AnalysisCheckpoint.markStepComplete(artifactDir, checkpoint, 7, skipped = true)
        } else {
          // Shouldn't happen with current logic, but handle defensively
          log.warn("Unexpected sex inference state - attempting Y-DNA analysis")
          onProgress(AnalysisProgress("Step 7/8: Determining Y-DNA haplogroup (unknown state)...", 0.70))
          val yDnaResult = runHaplogroupStep(bamPath, subject, currentSeqRun, getCurrentAlignment, TreeType.YDNA, artifactCtx, { pct =>
            onProgress(AnalysisProgress(s"Step 7/8: Y-DNA haplogroup (${(pct * 100).toInt}%)", 0.70 + pct * 0.22))
          })
          yDnaResult match {
            case Right(haplogroupResult) =>
              result = result.copy(yDnaHaplogroup = Some(haplogroupResult))
              checkpoint = AnalysisCheckpoint.markStepComplete(artifactDir, checkpoint, 7)
            case Left(error) =>
              result = result.copy(skippedYDna = true, skippedYDnaReason = Some(error))
              checkpoint = AnalysisCheckpoint.markStepComplete(artifactDir, checkpoint, 7, skipped = true)
          }
        }
      } else {
        onProgress(AnalysisProgress("Step 7/8: Y-DNA haplogroup (cached/skipped)", 0.92))
        if (checkpoint.yDnaSkipped) {
          result = result.copy(skippedYDna = true)
        }
      }

      // Step 8: Ancestral Composition Stub (0.92 - 1.0)
      if (!checkpoint.ancestryCompleted) {
        onProgress(AnalysisProgress("Step 8/8: Preparing for ancestry analysis...", 0.92))
        // This is a stub for future implementation
        // Will use autosomal SNPs from the VCF for ancestry composition
        result = result.copy(ancestryStub = true)
        onProgress(AnalysisProgress("Ancestry analysis will be available in a future update", 0.96))
        checkpoint = AnalysisCheckpoint.markStepComplete(artifactDir, checkpoint, 8)
      } else {
        onProgress(AnalysisProgress("Step 8/8: Ancestry analysis (cached)", 1.0))
        result = result.copy(ancestryStub = true)
      }

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
   * Note: This now takes alignment to access coverage metrics for minDepth calculation.
   */
  private def runCallableLociStep(
    bamPath: String,
    referencePath: String,
    seqRun: SequenceRun,
    alignment: Alignment,
    artifactCtx: ArtifactContext,
    onProgress: Double => Unit
  ): Either[String, (CallableLociResult, List[String])] = {
    val processor = new CallableLociProcessor()

    // Determine minDepth based on test type AND coverage
    // HiFi reads are highly accurate, so minDepth=2 is appropriate
    // Low-pass WGS (<=5x) should also use minDepth=2 to avoid excessive no-calls
    val meanCoverage = alignment.metrics.flatMap(_.meanCoverage).getOrElse(30.0)
    val isLowPass = meanCoverage <= 5.0
    val isHiFi = seqRun.testType.toUpperCase.contains("HIFI")
    val isLongRead = seqRun.testType.toUpperCase.contains("NANOPORE") ||
                     seqRun.testType.toUpperCase.contains("CLR") ||
                     seqRun.maxReadLength.exists(_ > 10000)

    val minDepth = if (isHiFi) {
      2  // HiFi: high accuracy, minDepth=2 is fine
    } else if (isLowPass) {
      // Low-pass data: use minDepth proportional to coverage
      // At 4x, use minDepth=2; at 2x, use minDepth=1
      math.max(1, (meanCoverage / 2).toInt)
    } else if (isLongRead) {
      3  // ONT/CLR: moderate accuracy, minDepth=3
    } else {
      4  // Illumina WGS at normal depth: standard minDepth=4
    }

    log.info(s"[CallableLoci] Using minDepth=$minDepth (testType=${seqRun.testType}, meanCov=${f"$meanCoverage%.1f"}x, isHiFi=$isHiFi, isLowPass=$isLowPass)")

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
   * Prefers using the cached whole-genome VCF from Step 5 (variant calling) if available.
   * Falls back to BAM-based calling if no cached VCF exists.
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

    // Check if we have a cached VCF from Step 5 (variant calling)
    val runId = seqRun.atUri.map(SubjectArtifactCache.extractIdFromUri).getOrElse("unknown-run")
    val alignId = alignment.atUri.map(SubjectArtifactCache.extractIdFromUri).getOrElse("unknown-alignment")
    val vcfDir = VcfCache.getVcfDir(subject.sampleAccession, runId, alignId)
    val cachedVcfPath = vcfDir.resolve("whole_genome.vcf.gz")

    val treeTypeStr = if (treeType == TreeType.YDNA) "Y-DNA" else "mtDNA"

    // Check for vendor-provided FASTA (mtDNA only)
    val vendorFasta = if (treeType == TreeType.MTDNA) {
      VcfCache.findMtDnaRunFasta(subject.sampleAccession, runId)
    } else {
      None
    }

    // Check for vendor-provided VCF
    val vendorVcf = if (treeType == TreeType.YDNA) {
      VcfCache.findYDnaRunVendorVcf(subject.sampleAccession, runId)
    } else {
      VcfCache.findMtDnaRunVendorVcf(subject.sampleAccession, runId)
    }

    if (vendorFasta.isDefined && treeType == TreeType.MTDNA) {
      // Use vendor-provided FASTA for mtDNA
      val vfasta = vendorFasta.get
      log.info(s"Using ${vfasta.vendor.displayName} FASTA for mtDNA haplogroup analysis")
      processor.analyzeFromFasta(
        fastaPath = vfasta.fastaPath,
        treeProviderType = treeProviderType,
        onProgress = (_, current, total) => {
          val pct = if (total > 0) current / total else 0.0
          onProgress(pct)
        },
        outputDir = None
      ).map(_.headOption.getOrElse(defaultHaplogroupResult))
    } else if (vendorVcf.isDefined) {
      // Use vendor-provided VCF
      val vvcf = vendorVcf.get
      log.info(s"Using ${vvcf.vendor.displayName} VCF for $treeTypeStr haplogroup analysis")
      processor.analyzeFromVcfFile(
        vcfPath = vvcf.vcfPath,
        referenceBuild = vvcf.referenceBuild,
        treeType = treeType,
        treeProviderType = treeProviderType,
        onProgress = (_, current, total) => {
          val pct = if (total > 0) current / total else 0.0
          onProgress(pct)
        },
        outputDir = None
      ).map(_.headOption.getOrElse(defaultHaplogroupResult))
    } else if (java.nio.file.Files.exists(cachedVcfPath)) {
      // Use cached VCF from variant calling step - much faster!
      log.info(s"Using cached whole-genome VCF for $treeTypeStr haplogroup analysis")
      processor.analyzeFromCachedVcf(
        sampleAccession = subject.sampleAccession,
        runId = runId,
        alignmentId = alignId,
        referenceBuild = alignment.referenceBuild,
        treeType = treeType,
        treeProviderType = treeProviderType,
        onProgress = (_, current, total) => {
          val pct = if (total > 0) current / total else 0.0
          onProgress(pct)
        }
      ).map(_.headOption.getOrElse(defaultHaplogroupResult))
    } else {
      // Fall back to BAM-based calling (slower, but works without Step 5)
      // This generates contig-specific VCFs which we'll save to the common VCF location
      log.info(s"No cached VCF found, using BAM-based calling for $treeTypeStr haplogroup analysis")

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

      val result = processor.analyze(
        bamPath,
        libraryStats,
        treeType,
        treeProviderType,
        (_, current, total) => {
          val pct = if (total > 0) current / total else 0.0
          onProgress(pct)
        },
        Some(artifactCtx)
      )

      // Copy the generated contig VCF to the common VCF location for future reuse
      // This ensures the VCF is available for subsequent analysis without re-calling
      result.foreach { _ =>
        val contigName = if (treeType == TreeType.YDNA) "chrY" else "chrM"
        val haplogroupDir = artifactCtx.getSubdir("haplogroup")
        val prefix = if (treeType == TreeType.YDNA) "ydna" else "mtdna"
        val sourceVcf = haplogroupDir.resolve(s"${prefix}_calls.vcf")

        if (java.nio.file.Files.exists(sourceVcf)) {
          // Save to common VCF location as contig-specific VCF
          val destVcf = vcfDir.resolve(s"$contigName.vcf.gz")
          try {
            java.nio.file.Files.createDirectories(vcfDir)
            // Compress and copy to standard location
            GatkRunner.run(Array(
              "SortVcf",
              "-I", sourceVcf.toString,
              "-O", destVcf.toString,
              "--CREATE_INDEX", "true"
            )) match {
              case Right(_) =>
                log.info(s"Saved $treeTypeStr VCF to common location: $destVcf")
              case Left(err) =>
                log.warn(s"Failed to copy VCF to common location: $err")
            }
          } catch {
            case e: Exception =>
              log.warn(s"Failed to save VCF to common location: ${e.getMessage}")
          }
        }
      }

      result.map(_.headOption.getOrElse(defaultHaplogroupResult))
    }
  }

  /** Default result when no haplogroup matches found */
  private val defaultHaplogroupResult = AnalysisHaplogroupResult(
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
  readMetrics: Option[ReadMetrics] = None,
  wgsMetrics: Option[WgsMetrics] = None,
  callableLociResult: Option[CallableLociResult] = None,
  sexInferenceResult: Option[SexInference.SexInferenceResult] = None,
  vcfInfo: Option[CachedVcfInfo] = None,
  mtDnaHaplogroup: Option[AnalysisHaplogroupResult] = None,
  yDnaHaplogroup: Option[AnalysisHaplogroupResult] = None,
  skippedYDna: Boolean = false,
  skippedYDnaReason: Option[String] = None,
  ancestryStub: Boolean = false  // Stub for future ancestry composition
)
