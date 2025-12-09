package com.decodingus.ancestry.processor

import com.decodingus.analysis.ArtifactContext
import com.decodingus.ancestry.model.*
import com.decodingus.ancestry.reference.{AncestryReferenceGateway, AncestryReferenceResult}
import com.decodingus.ancestry.report.AncestryReportWriter
import com.decodingus.config.FeatureToggles
import com.decodingus.model.LibraryStats
import com.decodingus.refgenome.ReferenceGateway

import java.nio.file.Path

/**
 * Orchestrates ancestry analysis pipeline.
 *
 * Pipeline steps:
 * 1. Check/download ancestry reference data
 * 2. Call genotypes at ancestry-informative SNP positions
 * 3. Project onto PCA space and estimate population proportions
 * 4. Generate reports
 */
class AncestryProcessor {

  private val ARTIFACT_SUBDIR = "ancestry"

  /**
   * Run ancestry analysis on a BAM/CRAM file.
   *
   * @param bamPath Path to aligned reads (BAM/CRAM)
   * @param libraryStats Library statistics from initial analysis
   * @param panelType AIMs (quick) or GenomeWide (detailed)
   * @param onProgress Progress callback (message, current, total)
   * @param artifactContext Optional context for organizing output artifacts
   * @return Either error message or ancestry result
   */
  def analyze(
    bamPath: String,
    libraryStats: LibraryStats,
    panelType: AncestryPanelType,
    onProgress: (String, Double, Double) => Unit,
    artifactContext: Option[ArtifactContext] = None
  ): Either[String, AncestryResult] = {

    val panelName = panelType match {
      case AncestryPanelType.Aims => "aims"
      case AncestryPanelType.GenomeWide => "genome-wide"
    }

    onProgress(s"Starting $panelName ancestry analysis...", 0.0, 1.0)

    // Step 1: Check ancestry reference data availability
    val ancestryGateway = new AncestryReferenceGateway((_, _) => {})
    val referenceBuild = libraryStats.referenceBuild

    ancestryGateway.checkAvailability(panelType, referenceBuild) match {
      case AncestryReferenceResult.Available(sitesVcf, alleleFreqs, pcaLoadings) =>
        analyzeWithData(
          bamPath,
          libraryStats,
          panelType,
          sitesVcf,
          alleleFreqs,
          pcaLoadings,
          onProgress,
          artifactContext
        )

      case AncestryReferenceResult.DownloadRequired(_, _, sizeMB, _) =>
        Left(s"Ancestry reference panel not available. Download required (~${sizeMB}MB). " +
          "Please run the download first.")

      case AncestryReferenceResult.Error(msg) =>
        Left(s"Failed to access ancestry reference data: $msg")
    }
  }

  /**
   * Run analysis with pre-loaded reference data.
   */
  private def analyzeWithData(
    bamPath: String,
    libraryStats: LibraryStats,
    panelType: AncestryPanelType,
    sitesVcf: Path,
    alleleFreqs: AlleleFrequencyMatrix,
    pcaLoadings: PCALoadings,
    onProgress: (String, Double, Double) => Unit,
    artifactContext: Option[ArtifactContext]
  ): Either[String, AncestryResult] = {

    val panelName = panelType match {
      case AncestryPanelType.Aims => "aims"
      case AncestryPanelType.GenomeWide => "genome-wide"
    }

    val artifactDir = artifactContext.map(_.getSubdir(ARTIFACT_SUBDIR))

    // Step 2: Resolve genome reference for genotype calling
    onProgress("Resolving genome reference...", 0.05, 1.0)
    val refGateway = new ReferenceGateway((_, _) => {})

    refGateway.resolve(libraryStats.referenceBuild).flatMap { referencePath =>

      // Step 3: Extract genotypes at ancestry-informative sites
      onProgress("Calling genotypes at ancestry-informative sites...", 0.1, 1.0)
      val genotypeExtractor = new GenotypeExtractor()

      genotypeExtractor.extractGenotypes(
        bamPath,
        referencePath.toString,
        sitesVcf.toFile,
        (msg, done, total) => onProgress(msg, 0.1 + done * 0.6, 1.0),
        artifactDir
      ).flatMap { genotypeResult =>

        // Step 4: Validate data quality
        val minSnps = panelType match {
          case AncestryPanelType.Aims => FeatureToggles.ancestryAnalysis.minSnpsAims
          case AncestryPanelType.GenomeWide => FeatureToggles.ancestryAnalysis.minSnpsGenomeWide
        }

        if (genotypeResult.snpsWithGenotype < minSnps) {
          Left(s"Insufficient data: ${genotypeResult.snpsWithGenotype} SNPs genotyped, " +
            s"minimum required is $minSnps. Consider using higher coverage data.")
        } else {
          // Step 5: Estimate ancestry proportions
          onProgress("Estimating ancestry proportions...", 0.75, 1.0)
          val estimator = new AncestryEstimator()

          val result = estimator.estimate(
            genotypeResult.genotypes,
            alleleFreqs,
            pcaLoadings,
            panelName
          )

          // Step 6: Write reports
          artifactDir.foreach { dir =>
            onProgress("Writing ancestry report...", 0.9, 1.0)
            AncestryReportWriter.writeReport(dir.toFile, result, panelType)
            AncestryReportWriter.writeJsonReport(dir.toFile, result)
          }

          onProgress("Ancestry analysis complete.", 1.0, 1.0)
          Right(result)
        }
      }
    }
  }

  /**
   * Run both AIMs and genome-wide analysis.
   * Returns results for both if successful.
   */
  def analyzeAll(
    bamPath: String,
    libraryStats: LibraryStats,
    onProgress: (String, Double, Double) => Unit,
    artifactContext: Option[ArtifactContext] = None
  ): Either[String, (AncestryResult, Option[AncestryResult])] = {

    // Run AIMs first (quick)
    onProgress("Running quick AIMs analysis...", 0.0, 1.0)
    analyze(bamPath, libraryStats, AncestryPanelType.Aims,
      (msg, done, total) => onProgress(msg, done * 0.3, 1.0),
      artifactContext
    ).flatMap { aimsResult =>

      // Try genome-wide (may not have reference data)
      onProgress("Attempting genome-wide analysis...", 0.3, 1.0)
      val genomeWideResult = analyze(bamPath, libraryStats, AncestryPanelType.GenomeWide,
        (msg, done, total) => onProgress(msg, 0.3 + done * 0.7, 1.0),
        artifactContext
      ).toOption

      Right((aimsResult, genomeWideResult))
    }
  }

  /**
   * Check if ancestry analysis can be performed for a sample.
   */
  def canAnalyze(referenceBuild: String): (Boolean, Boolean) = {
    val gateway = new AncestryReferenceGateway((_, _) => {})
    val aimsAvailable = gateway.checkAvailability(AncestryPanelType.Aims, referenceBuild) match {
      case _: AncestryReferenceResult.Available => true
      case _ => false
    }
    val genomeWideAvailable = gateway.checkAvailability(AncestryPanelType.GenomeWide, referenceBuild) match {
      case _: AncestryReferenceResult.Available => true
      case _ => false
    }
    (aimsAvailable, genomeWideAvailable)
  }
}
