package com.decodingus.workspace.model

import java.time.Instant

/**
 * DNA type for haplogroup reconciliation.
 */
enum DnaType:
  case Y_DNA, MT_DNA

/**
 * Compatibility level between runs.
 */
enum CompatibilityLevel:
  case COMPATIBLE // Same branch, different depths
  case MINOR_DIVERGENCE // Tip-level differences
  case MAJOR_DIVERGENCE // Branch-level split
  case INCOMPATIBLE // Different individuals

/**
 * Technology used for sequencing/testing.
 */
enum HaplogroupTechnology:
  case WGS, WES, BIG_Y, SNP_ARRAY, AMPLICON, STR_PANEL

/**
 * Method used to determine haplogroup.
 */
enum CallMethod:
  case SNP_PHYLOGENETIC, STR_PREDICTION, VENDOR_REPORTED

/**
 * How a conflict was resolved.
 */
enum ConflictResolution:
  case ACCEPT_MAJORITY, ACCEPT_HIGHER_QUALITY, ACCEPT_HIGHER_COVERAGE, UNRESOLVED, HETEROPLASMY

/**
 * A haplogroup call from a single sequencing run or chip profile.
 * Matches global Atmosphere schema: com.decodingus.atmosphere.defs#runHaplogroupCall
 */
case class RunHaplogroupCall(
                              sourceRef: String, // AT URI of SequenceRun, ChipProfile, or StrProfile
                              haplogroup: String,
                              confidence: Double,
                              callMethod: CallMethod,
                              score: Option[Double] = None,
                              supportingSnps: Option[Int] = None,
                              conflictingSnps: Option[Int] = None,
                              noCalls: Option[Int] = None,
                              technology: Option[HaplogroupTechnology] = None,
                              meanCoverage: Option[Double] = None,
                              treeProvider: Option[String] = None, // "ftdna", "decodingus"
                              treeVersion: Option[String] = None,
                              lineagePath: Option[List[String]] = None
                            )

/**
 * A single SNP call from one run.
 */
case class SnpCallFromRun(
                           runRef: String,
                           allele: String,
                           quality: Option[Double] = None,
                           depth: Option[Int] = None,
                           variantAlleleFrequency: Option[Double] = None // For heteroplasmy
                         )

/**
 * A conflict at a specific SNP position between runs.
 */
case class SnpConflict(
                        position: Long,
                        snpName: Option[String] = None,
                        contigAccession: Option[String] = None,
                        calls: List[SnpCallFromRun],
                        resolution: Option[ConflictResolution] = None,
                        resolvedValue: Option[String] = None
                      )

/**
 * Summary reconciliation status.
 * Matches global Atmosphere schema: com.decodingus.atmosphere.defs#reconciliationStatus
 */
case class ReconciliationStatus(
                                 compatibilityLevel: CompatibilityLevel,
                                 consensusHaplogroup: String,
                                 confidence: Double,
                                 divergencePoint: Option[String] = None,
                                 branchCompatibilityScore: Option[Double] = None,
                                 snpConcordance: Option[Double] = None,
                                 runCount: Int,
                                 warnings: List[String] = List.empty
                               )

/**
 * Reconciliation of haplogroup calls across multiple runs for a biosample.
 * Matches global Atmosphere schema: com.decodingus.atmosphere.haplogroupReconciliation
 *
 * Note: Global schema has this at specimen donor level, but for Edge App
 * we track at biosample level since users typically have one biosample per donor.
 */
case class HaplogroupReconciliation(
                                     atUri: Option[String] = None,
                                     meta: RecordMeta,
                                     biosampleRef: String, // AT URI of the biosample
                                     dnaType: DnaType,
                                     status: ReconciliationStatus,
                                     runCalls: List[RunHaplogroupCall],
                                     snpConflicts: List[SnpConflict] = List.empty,
                                     lastReconciliationAt: Option[Instant] = None
                                   ) {

  /**
   * Add or update a run call. Replaces existing call from same source.
   */
  def withRunCall(call: RunHaplogroupCall): HaplogroupReconciliation = {
    val filtered = runCalls.filterNot(_.sourceRef == call.sourceRef)
    copy(runCalls = call :: filtered)
  }

  /**
   * Remove a run call by source reference.
   */
  def removeRunCall(sourceRef: String): HaplogroupReconciliation = {
    copy(runCalls = runCalls.filterNot(_.sourceRef == sourceRef))
  }

  /**
   * Recalculate the consensus/status from current run calls.
   *
   * Uses lineage paths (when available) to compute branch compatibility via LCA.
   * Falls back to simple tier-based selection when lineage data is missing.
   */
  def recalculate(): HaplogroupReconciliation = {
    if (runCalls.isEmpty) {
      copy(
        status = status.copy(
          compatibilityLevel = CompatibilityLevel.COMPATIBLE,
          consensusHaplogroup = "",
          confidence = 0.0,
          runCount = 0,
          divergencePoint = None,
          branchCompatibilityScore = None,
          snpConcordance = None,
          warnings = List.empty
        ),
        snpConflicts = List.empty,
        lastReconciliationAt = Some(Instant.now())
      )
    } else {
      // Pick best call by quality tier then confidence
      val best = HaplogroupReconciliation.bestCall(runCalls)

      // Compute branch compatibility from lineage paths
      val (compatLevel, compatScore, divergence, warnings) =
        HaplogroupReconciliation.assessBranchCompatibility(runCalls, dnaType)

      // Compute SNP-level concordance and conflicts
      val (concordance, conflicts) = HaplogroupReconciliation.detectSnpConflicts(runCalls, dnaType)

      // Combine warnings from branch and SNP analysis
      val allWarnings = warnings ++ concordanceWarnings(concordance)

      copy(
        status = status.copy(
          consensusHaplogroup = best.haplogroup,
          confidence = best.confidence,
          runCount = runCalls.size,
          compatibilityLevel = compatLevel,
          branchCompatibilityScore = compatScore,
          divergencePoint = divergence,
          snpConcordance = concordance,
          warnings = allWarnings
        ),
        snpConflicts = conflicts,
        lastReconciliationAt = Some(Instant.now())
      )
    }
  }

  private def concordanceWarnings(concordance: Option[Double]): List[String] =
    concordance match {
      case Some(c) if c < 0.95 =>
        List(s"Low SNP concordance (${f"${c * 100}%.1f"}%%) — possible sample mix-up")
      case Some(c) if c < 0.99 =>
        List(s"SNP concordance ${f"${c * 100}%.1f"}%% — possible somatic variation or heteroplasmy")
      case _ => List.empty
    }
}

object HaplogroupReconciliation {

  /**
   * Select the best call by technology quality tier then confidence.
   */
  def bestCall(calls: List[RunHaplogroupCall]): RunHaplogroupCall =
    calls.maxBy { call =>
      val qualityTier = call.technology match {
        case Some(HaplogroupTechnology.WGS) => 3
        case Some(HaplogroupTechnology.BIG_Y) => 2
        case Some(HaplogroupTechnology.SNP_ARRAY) => 1
        case _ => 0
      }
      (qualityTier, call.confidence)
    }

  // ============================================
  // Branch Compatibility (LCA-based)
  // ============================================

  /**
   * Assess branch compatibility across all run calls using lineage paths.
   *
   * Computes the LCA (longest common prefix) of all lineage paths, then
   * calculates a compatibility score: LCA_depth / max(call depths).
   *
   * @return (compatibilityLevel, score, divergencePoint, warnings)
   */
  def assessBranchCompatibility(
                                  calls: List[RunHaplogroupCall],
                                  dnaType: DnaType
                                ): (CompatibilityLevel, Option[Double], Option[String], List[String]) = {
    if (calls.size <= 1) {
      return (CompatibilityLevel.COMPATIBLE, Some(1.0), None, List.empty)
    }

    val paths = calls.flatMap(_.lineagePath).filter(_.nonEmpty)
    if (paths.size < 2) {
      // Not enough lineage data for tree comparison — assume compatible
      return (CompatibilityLevel.COMPATIBLE, None, None,
        List("Insufficient lineage data for tree comparison"))
    }

    val lca = longestCommonPrefix(paths)
    val lcaDepth = lca.size
    val maxDepth = paths.map(_.size).max

    if (maxDepth == 0) {
      return (CompatibilityLevel.COMPATIBLE, None, None, List.empty)
    }

    val warnings = scala.collection.mutable.ListBuffer[String]()

    // Check for actual divergence: paths that go in different directions beyond the LCA
    val branchesBeyondLca = paths.map(_.drop(lcaDepth)).filter(_.nonEmpty).map(_.head).distinct
    val hasDivergence = branchesBeyondLca.size > 1

    // If no divergence, all paths are on the same branch (ancestor/descendant)
    // Score 1.0 = fully compatible, even if depths differ
    val score = if (hasDivergence) {
      lcaDepth.toDouble / maxDepth
    } else {
      1.0
    }

    val divergencePoint = if (hasDivergence) lca.lastOption else None

    if (hasDivergence) {
      warnings += s"Runs diverge after ${divergencePoint.getOrElse("root")}: ${branchesBeyondLca.mkString(", ")}"
    }

    val level = scoreToCompatibilityLevel(score)
    if (level == CompatibilityLevel.INCOMPATIBLE) {
      warnings += "Sample verification recommended — haplogroup calls are incompatible"
    }

    (level, Some(score), divergencePoint, warnings.toList)
  }

  /**
   * Find the longest common prefix of multiple lineage paths.
   * This gives the Last Common Ancestor (LCA) path.
   */
  def longestCommonPrefix(paths: List[List[String]]): List[String] = {
    if (paths.isEmpty) return List.empty
    val minLen = paths.map(_.size).min
    val reference = paths.head
    val prefixLen = (0 until minLen).takeWhile { i =>
      paths.forall(_(i) == reference(i))
    }.size
    reference.take(prefixLen)
  }

  /**
   * Map a branch compatibility score to a compatibility level.
   */
  def scoreToCompatibilityLevel(score: Double): CompatibilityLevel =
    if (score >= 0.8) CompatibilityLevel.COMPATIBLE
    else if (score >= 0.5) CompatibilityLevel.MINOR_DIVERGENCE
    else if (score >= 0.3) CompatibilityLevel.MAJOR_DIVERGENCE
    else CompatibilityLevel.INCOMPATIBLE

  // ============================================
  // SNP-Level Conflict Detection
  // ============================================

  /**
   * Detect SNP-level conflicts across runs.
   *
   * Compares supporting/conflicting SNP counts to estimate concordance.
   * When detailed SNP data isn't available, returns None.
   *
   * @return (concordance, conflicts)
   */
  def detectSnpConflicts(
                          calls: List[RunHaplogroupCall],
                          dnaType: DnaType
                        ): (Option[Double], List[SnpConflict]) = {
    if (calls.size < 2) return (None, List.empty)

    // Estimate concordance from supporting/conflicting SNP counts
    val callsWithSnpData = calls.filter(c => c.supportingSnps.isDefined || c.conflictingSnps.isDefined)
    if (callsWithSnpData.size < 2) return (None, List.empty)

    val totalSupporting = callsWithSnpData.flatMap(_.supportingSnps).sum
    val totalConflicting = callsWithSnpData.flatMap(_.conflictingSnps).sum
    val totalCalls = totalSupporting + totalConflicting

    val concordance = if (totalCalls > 0) {
      Some(totalSupporting.toDouble / totalCalls)
    } else {
      None
    }

    // For now we report aggregate concordance. Detailed per-SNP conflicts
    // require the full variant call data which isn't stored on RunHaplogroupCall.
    // The snpConflicts list will be populated when detailed SNP comparison
    // data is available from the analysis pipeline.
    (concordance, List.empty)
  }

  // ============================================
  // Factory Methods
  // ============================================

  /**
   * Create a new reconciliation record from a single run call.
   */
  def fromSingleRun(
                     biosampleRef: String,
                     dnaType: DnaType,
                     call: RunHaplogroupCall
                   ): HaplogroupReconciliation = {
    HaplogroupReconciliation(
      meta = RecordMeta.initial,
      biosampleRef = biosampleRef,
      dnaType = dnaType,
      status = ReconciliationStatus(
        compatibilityLevel = CompatibilityLevel.COMPATIBLE,
        consensusHaplogroup = call.haplogroup,
        confidence = call.confidence,
        runCount = 1,
        branchCompatibilityScore = Some(1.0)
      ),
      runCalls = List(call),
      lastReconciliationAt = Some(Instant.now())
    )
  }

  /**
   * Convert a RunHaplogroupCall to a HaplogroupResult for storage in Biosample.haplogroups.
   */
  def toHaplogroupResult(call: RunHaplogroupCall): HaplogroupResult = {
    HaplogroupResult(
      haplogroupName = call.haplogroup,
      score = call.score.getOrElse(call.confidence),
      matchingSnps = call.supportingSnps,
      mismatchingSnps = call.conflictingSnps,
      lineagePath = call.lineagePath,
      source = call.technology.map {
        case HaplogroupTechnology.WGS => "wgs"
        case HaplogroupTechnology.BIG_Y => "bigy"
        case HaplogroupTechnology.SNP_ARRAY => "chip"
        case _ => "other"
      },
      sourceRef = Some(call.sourceRef),
      treeProvider = call.treeProvider,
      treeVersion = call.treeVersion,
      analyzedAt = Some(Instant.now())
    )
  }
}
