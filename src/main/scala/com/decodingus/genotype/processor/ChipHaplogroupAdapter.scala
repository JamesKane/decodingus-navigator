package com.decodingus.genotype.processor

import com.decodingus.config.{FeatureToggles, UserPreferencesService}
import com.decodingus.genotype.model.GenotypeCall
import com.decodingus.haplogroup.model.{Haplogroup, HaplogroupResult}
import com.decodingus.haplogroup.scoring.{ConfidenceCalculator, HaplogroupScorer}
import com.decodingus.haplogroup.tree.{TreeProvider, TreeProviderType, TreeType}
import com.decodingus.haplogroup.vendor.{DecodingUsTreeProvider, FtdnaTreeProvider}
import com.decodingus.liftover.GenotypeLiftover

/**
 * Result of chip-based haplogroup analysis.
 *
 * @param treeType      Y-DNA or MT-DNA
 * @param results       Scored haplogroup results, sorted by score descending
 * @param snpsMatched   Number of chip positions that matched tree positions
 * @param snpsTotal     Total tree positions checked
 * @param topHaplogroup Best matching haplogroup name
 * @param confidence    Confidence level based on coverage
 */
case class ChipHaplogroupResult(
                                 treeType: TreeType,
                                 results: List[HaplogroupResult],
                                 snpsMatched: Int,
                                 snpsTotal: Int,
                                 topHaplogroup: String,
                                 confidence: Double
                               )

/**
 * Adapter to run haplogroup analysis on chip/array genotype data.
 *
 * Unlike WGS-based haplogroup analysis which uses GATK HaplotypeCaller to call
 * variants at tree positions, chip data already has genotypes at known positions.
 * The challenge is that chip positions are typically identified by rsID while
 * tree positions use genomic coordinates.
 *
 * Approach:
 * 1. Load the haplogroup tree (Y-DNA or MT-DNA)
 * 2. Build a map of tree positions to expected ref/alt alleles
 * 3. Match chip calls by position (most chips use GRCh37 coordinates)
 * 4. Convert to position->allele map for the scorer
 * 5. Run the standard HaplogroupScorer
 *
 * Limitations:
 * - Chip coverage of tree positions is typically 10-30% (most chips focus on common SNPs)
 * - Haplogroup resolution is limited compared to WGS (terminal branch may be upstream)
 * - No private variant detection (chips only call predefined positions)
 */
class ChipHaplogroupAdapter {

  private val scorer = new HaplogroupScorer()

  /**
   * Analyze chip data for haplogroup assignment.
   *
   * @param chipResult The processed chip data
   * @param treeType   Y-DNA or MT-DNA
   * @param onProgress Progress callback
   * @return Either error or haplogroup result
   */
  def analyze(
               chipResult: ChipProcessingResult,
               treeType: TreeType,
               onProgress: (String, Double, Double) => Unit
             ): Either[String, ChipHaplogroupResult] = {

    val typeName = if (treeType == TreeType.YDNA) "Y-DNA" else "mtDNA"
    onProgress(s"Loading $typeName haplogroup tree...", 0.0, 1.0)

    // Check if we have the appropriate calls
    val relevantCalls = treeType match {
      case TreeType.YDNA => chipResult.yDnaCalls
      case TreeType.MTDNA => chipResult.mtDnaCalls
    }

    if (relevantCalls.isEmpty) {
      return Left(s"No $typeName markers found in chip data")
    }

    // Most chip data uses GRCh37 coordinates
    val chipBuild = "GRCh37"

    // Get tree provider from user preferences
    val providerType = treeType match {
      case TreeType.YDNA =>
        if (UserPreferencesService.getYdnaTreeProvider.equalsIgnoreCase("decodingus")) TreeProviderType.DECODINGUS
        else TreeProviderType.FTDNA
      case TreeType.MTDNA =>
        if (UserPreferencesService.getMtdnaTreeProvider.equalsIgnoreCase("decodingus")) TreeProviderType.DECODINGUS
        else TreeProviderType.FTDNA
    }

    val treeProvider: TreeProvider = providerType match {
      case TreeProviderType.FTDNA => new FtdnaTreeProvider(treeType)
      case TreeProviderType.DECODINGUS => new DecodingUsTreeProvider(treeType)
    }

    // Determine which coordinate system to use for the tree
    // mtDNA: All chip vendors and trees use rCRS positions (1-16569), so use rCRS if available
    // Y-DNA: Need to match chip build to tree build, with liftover if necessary
    val treeBuild = treeType match {
      case TreeType.MTDNA =>
        // Prefer rCRS for mtDNA since that's what chip data uses
        if (treeProvider.supportedBuilds.contains("rCRS")) "rCRS"
        else treeProvider.sourceBuild // Fall back to source build (positions are same anyway)
      case TreeType.YDNA =>
        if (treeProvider.supportedBuilds.contains(chipBuild)) chipBuild
        else treeProvider.sourceBuild
    }

    // mtDNA positions use rCRS numbering (1-16569) which is consistent across builds
    // No liftover needed for mtDNA - chip companies and trees both use rCRS positions
    val needsLiftover = treeType match {
      case TreeType.MTDNA => false // rCRS coordinates are universal for mtDNA
      case TreeType.YDNA => GenotypeLiftover.needsLiftover(chipBuild, treeBuild)
    }

    if (needsLiftover) {
      onProgress(s"Chip data is $chipBuild, tree is $treeBuild. Preparing coordinate liftover...", 0.05, 1.0)
    }

    // Prepare lifted genotypes if needed
    val liftedCallsResult: Either[String, List[GenotypeCall]] = if (needsLiftover) {
      liftChipGenotypes(relevantCalls, chipBuild, treeBuild, onProgress)
    } else {
      Right(relevantCalls)
    }

    liftedCallsResult.flatMap { callsInTreeCoords =>
      treeProvider.loadTree(treeBuild).flatMap { tree =>
        onProgress("Building position map...", 0.3, 1.0)

        // Collect all tree positions with their ref/alt info
        val treeLoci = collectAllLoci(tree)
        val treePositions = treeLoci.map(l => l.position -> (l.ref, l.alt)).toMap

        onProgress(s"Matching ${callsInTreeCoords.size} chip markers to ${treePositions.size} tree positions...", 0.5, 1.0)

        // Build SNP calls map from chip data (now in tree coordinates)
        val snpCalls: Map[Long, String] = callsInTreeCoords.flatMap { call =>
          if (call.isNoCall) {
            None
          } else {
            // For haploid chromosomes (Y, MT), use allele1
            // The tree expects the actual called allele (ref or alt)
            Some(call.position.toLong -> call.allele1.toString)
          }
        }.toMap

        // Count how many tree positions we have coverage for
        val matchedPositions = treePositions.keys.count(snpCalls.contains)

        if (matchedPositions < getMinSnps(treeType)) {
          Left(s"Insufficient coverage: only $matchedPositions tree positions " +
            s"covered by chip data (minimum ${getMinSnps(treeType)} required). " +
            s"Chip-based $typeName haplogroup estimation may not be reliable.")
        } else {
          onProgress(s"Scoring haplogroups ($matchedPositions/${treePositions.size} positions covered)...", 0.7, 1.0)

          // Run the scorer
          val results = scorer.score(tree, snpCalls)

          if (results.isEmpty) {
            Left(s"No haplogroup matches found for $typeName")
          } else {
            onProgress("Analysis complete.", 1.0, 1.0)

            val topResult = results.head
            // Cap at 85% for chip data since we can't detect private variants or resolve terminal branches
            val confidence = ConfidenceCalculator.calculate(topResult, results, maxCap = 0.85)

            Right(ChipHaplogroupResult(
              treeType = treeType,
              results = results,
              snpsMatched = matchedPositions,
              snpsTotal = treePositions.size,
              topHaplogroup = topResult.name,
              confidence = confidence
            ))
          }
        }
      }
    }
  }

  /**
   * Lift chip genotypes from one reference build to another.
   * Handles reverse-complement when mapping to negative strand.
   */
  private def liftChipGenotypes(
                                 calls: List[GenotypeCall],
                                 fromBuild: String,
                                 toBuild: String,
                                 onProgress: (String, Double, Double) => Unit
                               ): Either[String, List[GenotypeCall]] = {
    onProgress(s"Lifting ${calls.size} genotypes from $fromBuild to $toBuild...", 0.1, 1.0)

    val liftover = new GenotypeLiftover(fromBuild, toBuild)

    liftover.initialize().map { liftOver =>
      val genotypeTuples = calls.map(c => (c.chromosome, c.position, c.allele1, c.allele2))

      val result = liftover.liftGenotypes(liftOver, genotypeTuples, (done, total) => {
        val pct = 0.1 + (done.toDouble / total) * 0.15
        onProgress(s"Lifting coordinates: $done/$total", pct, 1.0)
      })

      println(s"[ChipHaplogroupAdapter] Lifted ${result.lifted.size} genotypes, ${result.failedCount} failed")

      // Convert lifted results back to GenotypeCall objects
      val originalByPos = calls.map(c => c.position.toLong -> c).toMap

      result.lifted.flatMap { case (origPos, lifted) =>
        originalByPos.get(origPos).map { origCall =>
          origCall.copy(
            chromosome = lifted.chromosome,
            position = lifted.position,
            allele1 = lifted.allele1,
            allele2 = lifted.allele2
          )
        }
      }
    }
  }

  /**
   * Get minimum required SNPs for analysis.
   */
  private def getMinSnps(treeType: TreeType): Int = {
    treeType match {
      case TreeType.YDNA => FeatureToggles.chipData.minYMarkers
      case TreeType.MTDNA => FeatureToggles.chipData.minMtMarkers
    }
  }

  /**
   * Recursively collect all loci from the haplogroup tree.
   */
  private def collectAllLoci(tree: List[Haplogroup]): List[com.decodingus.haplogroup.model.Locus] = {
    tree.flatMap(collectLociFromHaplogroup)
  }

  private def collectLociFromHaplogroup(haplogroup: Haplogroup): List[com.decodingus.haplogroup.model.Locus] = {
    haplogroup.loci ++ haplogroup.children.flatMap(collectLociFromHaplogroup)
  }
}

object ChipHaplogroupAdapter {

  /**
   * Check if chip data is suitable for Y-DNA haplogroup analysis.
   */
  def canAnalyzeYDna(result: ChipProcessingResult): Boolean = {
    result.summary.yMarkersCalled.exists(_ >= FeatureToggles.chipData.minYMarkers)
  }

  /**
   * Check if chip data is suitable for mtDNA haplogroup analysis.
   */
  def canAnalyzeMtDna(result: ChipProcessingResult): Boolean = {
    result.summary.mtMarkersCalled.exists(_ >= FeatureToggles.chipData.minMtMarkers)
  }

  /**
   * Get a confidence description string.
   */
  def confidenceDescription(confidence: Double): String = {
    if (confidence >= 0.8) "High"
    else if (confidence >= 0.5) "Medium"
    else if (confidence >= 0.3) "Low"
    else "Very Low"
  }

  /**
   * Format haplogroup result for display.
   */
  def formatResult(result: ChipHaplogroupResult): String = {
    val treeLabel = if (result.treeType == TreeType.YDNA) "Y-DNA" else "mtDNA"
    f"$treeLabel: ${result.topHaplogroup} (${result.snpsMatched}/${result.snpsTotal} SNPs, " +
      f"${confidenceDescription(result.confidence)} confidence)"
  }
}
