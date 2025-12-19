package com.decodingus.analysis.sv

import com.decodingus.analysis.ArtifactContext

import java.io.{File, PrintWriter}
import java.nio.file.{Files, Path}
import java.time.Instant
import scala.util.Using

/**
 * Main orchestrator for structural variant calling.
 *
 * Coordinates the SV calling pipeline:
 * 1. Evidence collection (SvEvidenceWalker)
 * 2. Depth segmentation (DepthSegmenter)
 * 3. Evidence clustering (SvEvidenceClusterer)
 * 4. VCF output (SvVcfWriter)
 *
 * References:
 * - GATK-SV pipeline: Collins et al. "A structural variation reference for medical
 *   and population genetics." Nature 581.7809 (2020): 444-451.
 *   https://doi.org/10.1038/s41586-020-2287-8
 *
 * - Integrated SV calling: Kosugi et al. "Comprehensive evaluation of structural
 *   variation detection algorithms for whole genome sequencing."
 *   Genome Biology 20.1 (2019): 117.
 *   https://doi.org/10.1186/s13059-019-1720-5
 */
class SvCaller(config: SvCallerConfig = SvCallerConfig.default) {

  private val walker = new SvEvidenceWalker(config)
  private val segmenter = new DepthSegmenter(config)
  private val clusterer = new SvEvidenceClusterer(config)

  /**
   * Call structural variants from a BAM/CRAM file.
   *
   * @param bamPath           Path to BAM/CRAM file
   * @param referencePath     Path to reference genome
   * @param contigLengths     Map of contig names to lengths
   * @param referenceBuild    Reference build name (e.g., "GRCh38")
   * @param meanCoverage      Mean genome coverage
   * @param meanInsertSize    Mean insert size from library
   * @param insertSizeSd      Insert size standard deviation
   * @param meanReadLength    Mean read length
   * @param onProgress        Progress callback (message, fraction 0-1)
   * @param artifactContext   Optional context for artifact storage
   * @return Either error message or SV analysis result
   */
  def callStructuralVariants(
    bamPath: String,
    referencePath: String,
    contigLengths: Map[String, Long],
    referenceBuild: String,
    meanCoverage: Double,
    meanInsertSize: Double,
    insertSizeSd: Double,
    meanReadLength: Double = 150.0,
    onProgress: (String, Double) => Unit,
    artifactContext: Option[ArtifactContext] = None
  ): Either[String, SvAnalysisResult] = {

    // Check minimum coverage threshold
    if (meanCoverage < 10.0) {
      return Left(s"Coverage too low for SV calling (${meanCoverage}x, minimum 10x required)")
    }

    onProgress("Phase 1: Collecting SV evidence from BAM...", 0.0)

    // Phase 1: Evidence collection (0% - 60%)
    val progressAdapter: (String, Long, Long) => Unit = (msg, current, total) => {
      val fraction = if (total > 0) current.toDouble / total else 0.0
      onProgress(s"Collecting evidence: $msg", fraction * 0.6)
    }

    val evidence = walker.collectEvidence(
      bamPath = bamPath,
      referencePath = referencePath,
      contigLengths = contigLengths,
      expectedInsertSize = meanInsertSize,
      insertSizeSd = insertSizeSd,
      onProgress = progressAdapter
    ) match {
      case Right(ev) => ev
      case Left(error) => return Left(error)
    }

    onProgress("Phase 2: Segmenting read depth for CNV detection...", 0.6)

    // Phase 2: Depth segmentation (60% - 75%)
    val rawSegments = segmenter.segment(
      depthBins = evidence.depthBins,
      contigLengths = contigLengths,
      meanCoverage = meanCoverage,
      readLength = meanReadLength
    )

    val mergedSegments = segmenter.mergeNearbySegments(rawSegments)

    onProgress("Phase 3: Clustering PE/SR evidence...", 0.75)

    // Phase 3: Evidence clustering (75% - 90%)
    val svCalls = clusterer.cluster(evidence, mergedSegments)

    onProgress("Phase 4: Writing results...", 0.9)

    // Phase 4: Write outputs (90% - 100%)
    artifactContext.foreach { ctx =>
      writeArtifacts(ctx, evidence, mergedSegments, svCalls, referenceBuild)
    }

    onProgress("SV calling complete.", 1.0)

    Right(SvAnalysisResult(
      svCalls = svCalls.filter(_.filter == "PASS"),
      totalDiscordantPairs = evidence.totalDiscordantPairs,
      totalSplitReads = evidence.totalSplitReads,
      cnvSegments = mergedSegments.size,
      analysisTimestamp = Instant.now(),
      referenceBuild = referenceBuild,
      meanCoverage = meanCoverage
    ))
  }

  /**
   * Write analysis artifacts to the cache directory.
   */
  private def writeArtifacts(
    ctx: ArtifactContext,
    evidence: SvEvidenceCollection,
    segments: List[DepthSegment],
    calls: List[SvCall],
    referenceBuild: String
  ): Unit = {
    val outputDir = ctx.getSubdir("sv")
    outputDir.toFile.mkdirs()

    // Write depth segments
    writeDepthSegments(outputDir.resolve("depth_segments.tsv"), segments)

    // Write evidence summary
    writeEvidenceSummary(outputDir.resolve("evidence_summary.txt"), evidence)

    // Write VCF
    val vcfPath = outputDir.resolve("structural_variants.vcf.gz")
    val vcfWriter = new SvVcfWriter(config)
    vcfWriter.write(calls, vcfPath, evidence.sampleName, referenceBuild)

    // Write metadata
    writeMetadata(outputDir.resolve("sv_metadata.json"), calls, referenceBuild, vcfPath.toString)
  }

  /**
   * Write depth segments to TSV file.
   */
  private def writeDepthSegments(path: Path, segments: List[DepthSegment]): Unit = {
    Using.resource(new PrintWriter(path.toFile)) { writer =>
      writer.println("chrom\tstart\tend\tmean_depth\tlog2_ratio\tz_score\tnum_bins\tsv_type")
      segments.foreach { seg =>
        writer.println(f"${seg.chrom}\t${seg.start}\t${seg.end}\t${seg.meanDepth}%.2f\t${seg.log2Ratio}%.3f\t${seg.zScore}%.2f\t${seg.numBins}\t${seg.svType}")
      }
    }
  }

  /**
   * Write evidence collection summary.
   */
  private def writeEvidenceSummary(path: Path, evidence: SvEvidenceCollection): Unit = {
    Using.resource(new PrintWriter(path.toFile)) { writer =>
      writer.println("## SV EVIDENCE SUMMARY")
      writer.println(s"Sample: ${evidence.sampleName}")
      writer.println(s"Expected insert size: ${evidence.expectedInsertSize}")
      writer.println(s"Insert size SD: ${evidence.insertSizeSd}")
      writer.println()
      writer.println("## DISCORDANT PAIRS")
      writer.println(s"Total: ${evidence.totalDiscordantPairs}")
      writer.println(s"Inter-chromosomal: ${evidence.interChromosomalPairs.size}")
      writer.println(s"Insert size outliers: ${evidence.discordantPairs.count(_.reason == DiscordantReason.InsertSizeOutlier)}")
      writer.println(s"Wrong orientation: ${evidence.discordantPairs.count(_.reason == DiscordantReason.WrongOrientation)}")
      writer.println()
      writer.println("## SPLIT READS")
      writer.println(s"Total: ${evidence.totalSplitReads}")
      writer.println()
      writer.println("## DEPTH BINS")
      evidence.depthBins.foreach { case (chrom, bins) =>
        val nonZero = bins.count(_ > 0)
        val total = bins.length
        writer.println(s"$chrom: $nonZero / $total bins with reads")
      }
    }
  }

  /**
   * Write SV metadata JSON.
   */
  private def writeMetadata(path: Path, calls: List[SvCall], referenceBuild: String, vcfPath: String): Unit = {
    val passCalls = calls.filter(_.filter == "PASS")
    val metadata = CachedSvInfo(
      vcfPath = vcfPath,
      indexPath = vcfPath + ".tbi",
      referenceBuild = referenceBuild,
      createdAt = Instant.now().toString,
      svCallCount = passCalls.size,
      deletionCount = passCalls.count(_.svType == SvType.DEL),
      duplicationCount = passCalls.count(_.svType == SvType.DUP),
      inversionCount = passCalls.count(_.svType == SvType.INV),
      translocationCount = passCalls.count(_.svType == SvType.BND)
    )

    import io.circe.syntax.*
    val json = metadata.asJson.spaces2
    Files.writeString(path, json)
  }
}
