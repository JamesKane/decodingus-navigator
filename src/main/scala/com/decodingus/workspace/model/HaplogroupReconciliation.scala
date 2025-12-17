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
                              treeVersion: Option[String] = None
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
   * Returns updated reconciliation with new status.
   */
  def recalculate(): HaplogroupReconciliation = {
    if (runCalls.isEmpty) {
      copy(status = status.copy(
        compatibilityLevel = CompatibilityLevel.COMPATIBLE,
        consensusHaplogroup = "",
        confidence = 0.0,
        runCount = 0
      ))
    } else {
      // Simple reconciliation: pick best by quality tier then confidence
      val sorted = runCalls.sortBy { call =>
        val qualityTier = call.technology match {
          case Some(HaplogroupTechnology.WGS) => 3
          case Some(HaplogroupTechnology.BIG_Y) => 2
          case Some(HaplogroupTechnology.SNP_ARRAY) => 1
          case _ => 0
        }
        (-qualityTier, -call.confidence)
      }
      val best = sorted.head

      copy(
        status = status.copy(
          consensusHaplogroup = best.haplogroup,
          confidence = best.confidence,
          runCount = runCalls.size,
          compatibilityLevel = CompatibilityLevel.COMPATIBLE // TODO: actual tree comparison
        ),
        lastReconciliationAt = Some(Instant.now())
      )
    }
  }
}

object HaplogroupReconciliation {

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
        runCount = 1
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
